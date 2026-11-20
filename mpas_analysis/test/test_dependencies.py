# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/main/LICENSE
"""
Check that the dependencies declared in the four places we list them stay in
sync:

* ``dev-spec.txt``      -- the conda development spec
* ``pixi.toml``         -- the pixi workspace
* ``ci/recipe/recipe.yaml`` -- the rattler-build recipe
* ``pyproject.toml``    -- the PyPI package metadata

``pyproject.toml`` is deliberately a subset: several dependencies are only
available from conda-forge and have no PyPI equivalent.  Those are listed in
``CONDA_ONLY`` below, and the tests check that list is exactly right so that a
new conda-only dependency has to be added consciously.
"""

import re
import tomllib
from pathlib import Path

import pytest

# Packages that appear in the conda dependency lists but deliberately not in
# `pyproject.toml`, because they are not distributed on PyPI (or `python`
# itself, which `pyproject.toml` expresses as `requires-python`).
CONDA_ONLY = {
    'cartopy-offlinedata',
    'esmf',
    'geometric-features',
    'mache',
    'mpas-tools',
    'nco',
    'pyremap',
    'python',
}

# conda calls it `matplotlib-base`, PyPI calls it `matplotlib`
NAME_ALIASES = {'matplotlib-base': 'matplotlib'}


def _repo_root():
    root = Path(__file__).resolve().parents[2]
    if not (root / 'pixi.toml').exists():
        pytest.skip('dependency files are not available from an installed '
                    'package; these tests only run from a source checkout')
    return root


def _canonical(name):
    """Normalize a package name for comparison across the four files."""
    name = name.strip().lower().replace('_', '-')
    return NAME_ALIASES.get(name, name)


def _norm_constraint(constraint):
    """Normalize a version constraint.  pixi spells "any version" as ``*``
    while the other three files simply leave the constraint off."""
    constraint = constraint.replace(' ', '')
    return '' if constraint == '*' else constraint


def _split_spec(spec):
    """Split e.g. ``cartopy >=0.18.0`` into ``('cartopy', '>=0.18.0')``."""
    spec = spec.strip()
    match = re.match(r'^([A-Za-z0-9._-]+)\s*(.*)$', spec)
    name, constraint = match.group(1), match.group(2)
    return _canonical(name), _norm_constraint(constraint)


def _parse_dev_spec(root):
    """Parse ``dev-spec.txt`` into ``{section: {name: constraint}}`` plus the
    conda build pins (lines of the form ``name=version=build``)."""
    sections = {}
    builds = {}
    current = None
    for line in (root / 'dev-spec.txt').read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('#'):
            heading = line.lstrip('#').strip().lower()
            if heading in ('base', 'development', 'documentation'):
                current = heading
                sections[current] = {}
            continue
        if current is None:
            continue
        # a conda build pin, e.g. `esmf=*=mpi_mpich_*`
        if line.count('=') >= 2 and ' ' not in line:
            name, _, build = line.split('=', 2)
            builds[_canonical(name)] = build
            continue
        name, constraint = _split_spec(line)
        sections[current][name] = constraint
    return sections, builds


def _parse_pixi(root):
    """Parse ``pixi.toml`` runtime/dev/docs dependencies and build pins."""
    data = tomllib.loads((root / 'pixi.toml').read_text())

    def convert(table):
        deps = {}
        builds = {}
        for name, value in table.items():
            key = _canonical(name)
            if isinstance(value, dict):
                deps[key] = _norm_constraint(value.get('version', '*'))
                if 'build' in value:
                    builds[key] = value['build']
            else:
                deps[key] = _norm_constraint(value)
        return deps, builds

    runtime, builds = convert(data['dependencies'])
    features = data.get('feature', {})
    dev, _ = convert(features.get('dev', {}).get('dependencies', {}))
    docs, _ = convert(features.get('docs', {}).get('dependencies', {}))
    return {'base': runtime, 'development': dev, 'documentation': docs}, builds


def _parse_recipe(root):
    """Parse the ``requirements: run:`` list out of the rattler-build recipe,
    expanding ``${{ ... }}`` references against the recipe ``context``."""
    text = (root / 'ci' / 'recipe' / 'recipe.yaml').read_text()

    context = {}
    match = re.search(r'^context:\n((?:[ \t]+\S.*\n)+)', text, re.MULTILINE)
    assert match is not None, 'could not find the `context` block in recipe.yaml'
    for line in match.group(1).splitlines():
        key, _, value = line.strip().partition(':')
        context[key.strip()] = value.strip().strip('"\'')

    match = re.search(
        r'^requirements:\n(?:.*\n)*?^  run:\n((?:^    - .*\n)+)',
        text, re.MULTILINE)
    assert match is not None, \
        'could not find the `requirements: run:` block in recipe.yaml'

    def expand(value):
        return re.sub(r'\$\{\{\s*(\w+)\s*\}\}',
                      lambda m: context[m.group(1)], value)

    deps = {}
    for line in match.group(1).splitlines():
        name, constraint = _split_spec(expand(line.strip()[2:]))
        deps[name] = constraint
    return deps


def _parse_pyproject(root):
    data = tomllib.loads((root / 'pyproject.toml').read_text())
    project = data['project']
    runtime = dict(_split_spec(dep) for dep in project['dependencies'])
    optional = {
        key: dict(_split_spec(dep) for dep in value)
        for key, value in project.get('optional-dependencies', {}).items()
    }
    build = dict(_split_spec(dep)
                 for dep in data['build-system']['requires'])
    return runtime, optional, build, project['requires-python'].replace(' ', '')


# Tests #######################################################################

def test_conda_runtime_deps_match():
    """dev-spec.txt, pixi.toml and recipe.yaml must list the same runtime
    dependencies with the same version constraints."""
    root = _repo_root()
    dev_spec, _ = _parse_dev_spec(root)
    pixi, _ = _parse_pixi(root)
    recipe = _parse_recipe(root)

    assert dev_spec['base'] == pixi['base'], \
        'dev-spec.txt and pixi.toml runtime dependencies differ'
    assert dev_spec['base'] == recipe, \
        'dev-spec.txt and ci/recipe/recipe.yaml runtime dependencies differ'


def test_conda_build_pins_match():
    """Build strings (e.g. the MPI flavor of ESMF) must agree between
    dev-spec.txt and pixi.toml."""
    root = _repo_root()
    _, dev_spec_builds = _parse_dev_spec(root)
    _, pixi_builds = _parse_pixi(root)
    assert dev_spec_builds == pixi_builds, \
        'conda build pins differ between dev-spec.txt and pixi.toml'


def test_pyproject_runtime_deps_match_conda():
    """Every dependency shared with the conda lists must carry the same
    version constraint, and the conda-only set must be exactly CONDA_ONLY."""
    root = _repo_root()
    dev_spec, _ = _parse_dev_spec(root)
    conda = dev_spec['base']
    pyproject, _, _, _ = _parse_pyproject(root)

    assert set(pyproject) - set(conda) == set(), \
        'pyproject.toml lists dependencies missing from the conda lists'
    assert set(conda) - set(pyproject) == CONDA_ONLY, \
        ('the set of conda-only dependencies changed; update CONDA_ONLY in '
         'this test if that is intended')

    mismatched = {name: (constraint, conda[name])
                  for name, constraint in pyproject.items()
                  if constraint != conda[name]}
    assert not mismatched, \
        f'version constraints differ (pyproject, conda): {mismatched}'


def test_requires_python_consistent():
    """`requires-python` must match the `python` constraint everywhere."""
    root = _repo_root()
    dev_spec, _ = _parse_dev_spec(root)
    pixi, _ = _parse_pixi(root)
    recipe = _parse_recipe(root)
    _, _, _, requires_python = _parse_pyproject(root)

    assert dev_spec['base']['python'] == requires_python
    assert pixi['base']['python'] == requires_python
    assert recipe['python'] == requires_python


@pytest.mark.parametrize('section,extra', [('development', 'dev'),
                                           ('documentation', 'docs')])
def test_extra_deps_match(section, extra):
    """The dev and docs dependency groups must agree across dev-spec.txt,
    pixi.toml and pyproject.toml.

    `pyproject.toml` may omit build-time requirements that it declares under
    `[build-system] requires` instead, so those are allowed to be absent.
    """
    root = _repo_root()
    dev_spec, _ = _parse_dev_spec(root)
    pixi, _ = _parse_pixi(root)
    pyproject_runtime, optional, build_requires, _ = _parse_pyproject(root)

    assert dev_spec[section] == pixi[section], \
        f'{section} dependencies differ between dev-spec.txt and pixi.toml'

    declared = optional.get(extra, {})
    missing = set(dev_spec[section]) - set(declared)
    assert missing <= set(build_requires), \
        (f'{extra} dependencies missing from pyproject.toml and not declared '
         f'under [build-system] requires: {sorted(missing - set(build_requires))}')

    mismatched = {name: (constraint, dev_spec[section][name])
                  for name, constraint in declared.items()
                  if constraint != dev_spec[section].get(name)}
    assert not mismatched, \
        f'{extra} version constraints differ (pyproject, dev-spec): {mismatched}'
