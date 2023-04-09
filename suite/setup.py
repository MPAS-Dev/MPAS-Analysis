#!/usr/bin/env python
import argparse
import os
from mache import discover_machine, MachineInfo
from jinja2 import Template
import shutil


def main():
    parser = argparse.ArgumentParser(
        description='Set up a run of the test suite')
    parser.add_argument('-p', dest='python', required=True,
                        help='the python version')
    parser.add_argument('-r', dest='run', required=True,
                        help='the run name')
    parser.add_argument('-b', dest='branch', required=True,
                        help='the branch name')
    parser.add_argument('-e', dest='conda_env', help='the conda environment')
    parser.add_argument('--no_polar_regions', dest='polar_regions',
                        action='store_false',
                        help='whether to run mpas_analysis with '
                             '--polar_regions flag')
    parser.add_argument('--copy_docs', dest='copy_docs', action='store_true',
                        help='whether to copy the docs to the html dir')
    parser.add_argument('--clean', dest='clean', action='store_true',
                        help='whether to delete existing test-suite dirs for a '
                             'fresh start')

    args = parser.parse_args()

    machine = discover_machine()

    machine_info = MachineInfo(machine=machine)
    account, partition, configuration, qos = \
        machine_info.get_account_defaults()

    use_e3sm_unified = 'E3SMU_SCRIPT' in os.environ
    if use_e3sm_unified:
        e3sm_unified_script = os.environ['E3SMU_SCRIPT']
        args.branch = \
            os.path.splitext(os.path.basename(e3sm_unified_script))[0]
    else:
        e3sm_unified_script = None

    if machine == 'chrysalis':
        # we don't want the default, which is 'debug'
        partition = 'compute'

    username = machine_info.username
    web_section = machine_info.config['web_portal']
    web_base = os.path.join(web_section['base_path'], web_section['username'])
    html_base = f'{web_base}/analysis_testing'
    if args.run == 'QU480':
        simulation = '20200305.A_WCYCL1850.ne4_oQU480.anvil'
        mesh = 'QU480'
    else:
        simulation = '20230406.GMPAS-IAF-ISMF.T62_oQU240wLI.chrysalis'
        mesh = 'oQU240wLI'
    if machine in ['anvil', 'chrysalis']:
        input_base = '/lcrc/group/e3sm/public_html/diagnostics/mpas_analysis/example_simulations'
        output_base = f'/lcrc/group/e3sm/{username}/analysis_testing'
    elif machine == 'pm-cpu':
        input_base = '/global/cfs/cdirs/e3sm/diagnostics/mpas_analysis/example_simulations'
        scratch = os.environ['SCRATCH']
        output_base = f'{scratch}/analysis_testing'
    elif machine == 'compy':
        input_base = '/compyfs/diagnostics/mpas_analysis/example_simulations'
        output_base = f'/compyfs/{username}/analysis_testing'
    else:
        raise ValueError(f'Machine {machine} is not set up for the test suite '
                         f'yet.')

    suite_path = f'{machine}_test_suite'
    if args.clean:

        paths = [suite_path, os.path.join(output_base, machine, args.branch),
                 os.path.join(html_base, machine, args.branch)]
        for path in paths:
            try:
                shutil.rmtree(path)
            except FileNotFoundError:
                pass

    try:
        os.makedirs(suite_path)
    except FileExistsError:
        pass

    if args.copy_docs:
        docs_path = os.path.join(html_base, machine, args.branch, 'docs')
        try:
            shutil.rmtree(docs_path)
        except FileNotFoundError:
            pass
        shutil.copytree(os.path.join('docs', '_build', 'html'), docs_path)

    if mesh == 'QU480':
        generate = "['all', 'no_BGC', 'no_icebergs', 'no_index', 'no_eke',\n" \
                   "            'no_landIceCavities', 'no_waves']"
        end_year = '5'
    elif mesh == 'oQU240wLI':
        generate = "['all', 'no_BGC', 'no_icebergs', 'no_index', 'no_eke', " \
                   "'no_waves']"
        end_year = '10'
    else:
        raise ValueError(f'Unexpected mesh: {mesh}')

    if args.run == 'mesh_rename':
        mesh = f'new_{mesh}'

    sbatch = list()
    if account is not None:
        sbatch.append(f'#SBATCH -A {account}')
    if configuration is not None:
        sbatch.append(f'#SBATCH -C {configuration}')
    if partition is not None:
        sbatch.append(f'#SBATCH -p {partition}')
    if qos is not None:
        sbatch.append(f'#SBATCH --qos {qos}')

    sbatch = '\n'.join(sbatch)

    conda_base = os.path.abspath(
        os.path.join(os.environ['CONDA_EXE'], '..', '..'))

    if args.conda_env is not None:
        conda_env = args.conda_env
    else:
        conda_env = f'test_mpas_analysis_py{args.python}'

    config = os.path.join(suite_path, f'{args.run}.cfg')
    config_from_job = os.path.join('..', f'{args.run}.cfg')

    if args.run in ['main', 'ctrl']:
        out_subdir = os.path.join(machine, args.branch, f'main_py{args.python}')
    else:
        out_subdir = os.path.join(machine, args.branch, args.run)
    out_common_dir = os.path.join(machine, args.branch)

    with open(os.path.join('suite', 'template.cfg')) as template_file:
        template_data = template_file.read()
    template = Template(template_data)
    config_text = template.render(
        use_e3sm_unified=use_e3sm_unified, run_name=args.run,
        input_base=input_base, simulation=simulation, mesh=mesh,
        output_base=output_base, html_base=html_base, out_subdir=out_subdir,
        generate=generate, end_year=end_year)
    with open(config, 'w') as config_file:
        config_file.write(config_text)

    if args.run in ['main_vs_ctrl', 'moc_am', 'no_ncclimo', 'wc_defaults']:
        # add the run-specific config second
        config_from_job = ' '.join(
            [config_from_job,
             os.path.join('..', '..', 'suite', f'{args.run}.cfg')])

    if args.run not in ['main', 'ctrl']:
        try:
            os.makedirs(os.path.join(suite_path, args.run))
        except FileExistsError:
            pass
        job = os.path.join(suite_path, args.run, 'job_script.bash')

        if args.polar_regions:
            flags = '--polar_regions'
        else:
            flags = ''

        with open(os.path.join('suite', 'job_script.bash')) as template_file:
            template_data = template_file.read()
        template = Template(template_data)
        job_text = template.render(
            sbatch=sbatch, conda_base=conda_base,
            use_e3sm_unified=use_e3sm_unified,
            e3sm_unified_script=e3sm_unified_script, conda_env=conda_env,
            machine=machine, flags=flags, config=config_from_job,
            html_base=html_base, out_subdir=out_subdir,
            out_common_dir=out_common_dir)
        with open(job, 'w') as job_file:
            job_file.write(job_text)


if __name__ == '__main__':
    main()
