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

    if machine == 'chrysalis':
        # we don't want the default, which is 'debug'
        partition = 'compute'

    if machine in ['anvil', 'chrysalis']:
        input_base = '/lcrc/group/e3sm/ac.xylar/acme_scratch/anvil'
        output_base = '/lcrc/group/e3sm/ac.xylar/analysis_testing'
        html_base = '/lcrc/group/e3sm/public_html/diagnostic_output/ac.xylar/analysis_testing'
        if args.run == 'QU480':
            simulation = '20200305.A_WCYCL1850.ne4_oQU480.anvil'
            mesh = 'QU480'
        else:
            simulation = '20201025.GMPAS-IAF.T62_oQU240wLI.anvil'
            mesh = 'oQU240wLI'
    elif machine == 'cori-haswell':
        input_base = '/global/cfs/cdirs/e3sm/xylar'
        output_base = '/global/cscratch1/sd/xylar/analysis_testing'
        html_base = '/global/cfs/cdirs/e3sm/www/xylar/analysis_testing'
        simulation = '20200305.A_WCYCL1850.ne4_oQU480.anvil'
        mesh = 'QU480'
    elif machine == 'compy':
        input_base = '/compyfs/asay932/analysis_testing/test_output'
        output_base = '/compyfs/asay932/analysis_testing'
        html_base = '/compyfs/www/asay932/analysis_testing'
        simulation = '20200305.A_WCYCL1850.ne4_oQU480.anvil'
        mesh = 'QU480'
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
                   "            'no_landIceCavities']"
        end_year = '5'
    elif mesh == 'oQU240wLI':
        generate = "['all', 'no_BGC', 'no_icebergs', 'no_index', 'no_eke']"
        end_year = '8'
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

    if args.run == 'ctrl':
        out_subdir = os.path.join(machine, args.branch, f'main_py{args.python}')
    else:
        out_subdir = os.path.join(machine, args.branch, args.run)

    if machine == 'cori-haswell':
        execute_options = \
            '# the number of MPI tasks to use in creating mapping files (1 means tasks run in\n' \
            '# serial, the default)\n' \
            'mapMpiTasks = 1\n' \
            '\n' \
            '# "None" if ESMF should perform mapping file generation in serial without a\n' \
            '# command, or one of "srun" or "mpirun" if it should be run in parallel (or ins\n' \
            '# serial but with a command)\n' \
            'mapParallelExec = None'
    else:
        execute_options = ''

    with open(os.path.join('suite', 'template.cfg')) as template_file:
        template_data = template_file.read()
    template = Template(template_data)
    config_text = template.render(
        run_name=args.run, input_base=input_base, simulation=simulation,
        mesh=mesh, output_base=output_base, html_base=html_base,
        out_subdir=out_subdir, generate=generate, end_year=end_year,
        execute_options=execute_options)
    with open(config, 'w') as config_file:
        config_file.write(config_text)

    if args.run in ['main_vs_ctrl', 'no_ncclimo', 'wc_defaults']:
        # add the run-specific config second
        config_from_job = ' '.join(
            [config_from_job,
             os.path.join('..', '..', 'suite', f'{args.run}.cfg')])

    if args.run != 'ctrl':
        try:
            os.makedirs(os.path.join(suite_path, args.run))
        except FileExistsError:
            pass
        job = os.path.join(suite_path, args.run, 'job_script.bash')

        if args.polar_regions:
            flags = '--polar_regions'
        else:
            flags = ''

        if machine == 'cori-haswell':
            parallel_exec = ''
        else:
            prallel_exec = 'srun -N 1 -n 1'

        with open(os.path.join('suite', 'job_script.bash')) as template_file:
            template_data = template_file.read()
        template = Template(template_data)
        job_text = template.render(
            sbatch=sbatch, conda_base=conda_base, conda_env=conda_env,
            machine=machine, flags=flags, config=config_from_job,
            parallel_exec=parallel_exec, html_base=html_base)
        with open(job, 'w') as job_file:
            job_file.write(job_text)


if __name__ == '__main__':
    main()
