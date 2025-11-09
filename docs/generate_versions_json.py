#!/usr/bin/env python
import argparse
import json
import os
import re


def version_key(name):
    """Key function for sorting versions."""
    match = re.match(r'^(\d+)\.(\d+)\.(\d+)$', name)
    if match:
        # Sort by major, minor, patch
        return tuple(map(int, match.groups()))
    return ()


# Mode: local or production
parser = argparse.ArgumentParser(
    description='Generate versions.json for MPAS Analysis documentation.')
parser.add_argument(
    '--local',
    action='store_true',
    help='Generate versions.json for local build.'
)
args = parser.parse_args()
local = args.local
base_dir = '_build/html' if local else 'gh-pages'
shared_dir = os.path.join(base_dir, 'shared')

entries = []

if not os.path.exists(base_dir) or not os.listdir(base_dir):
    raise FileNotFoundError(
        f"Base directory '{base_dir}' does not exist or is empty.")

versions = os.listdir(base_dir)
numeric_versions = []
non_numeric_versions = []

for version in versions:
    # Check if it matches version pattern
    if re.match(r'^\d+\.\d+\.\d+$', version):
        numeric_versions.append(version)
    else:
        non_numeric_versions.append(version)

# Sort numeric versions by major, minor, patch in descending order
numeric_versions.sort(key=version_key, reverse=True)
# Sort non-numeric versions alphabetically
non_numeric_versions.sort()

# Combine the sorted lists
versions = non_numeric_versions + numeric_versions

if 'main' in versions:
    versions.insert(0, versions.pop(versions.index('main')))

for name in versions:
    path = os.path.join(base_dir, name)
    if os.path.isdir(path) and name not in ('shared', '.git'):
        entries.append({
            'version': name,
            'url': f'../{name}/' if local else f'/MPAS-Analysis/{name}/'
        })

os.makedirs(shared_dir, exist_ok=True)
with open(os.path.join(shared_dir, 'versions.json'), 'w') as f:
    json.dump(entries, f, indent=2)

