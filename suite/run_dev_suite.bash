#!/usr/bin/env bash

set -e

env_name=mpas_dev

conda_base=$(dirname $(dirname $CONDA_EXE))
source $conda_base/etc/profile.d/conda.sh

export HDF5_USE_FILE_LOCKING=FALSE

branch=$(git symbolic-ref --short HEAD)

# test building the docs
conda activate ${env_name}
cd docs
make clean
make html
cd ..

machine=$(python -c "from mache import discover_machine; print(discover_machine())")

py=3.11
./suite/setup.py -p ${py} -r main_py${py} -b ${branch} --copy_docs --clean -e ${env_name}
./suite/setup.py -p ${py} -r wc_defaults -b ${branch} --no_polar_regions -e ${env_name}
./suite/setup.py -p ${py} -r moc_am -b ${branch} -e ${env_name}
./suite/setup.py -p ${py} -r no_ncclimo -b ${branch} -e ${env_name}
./suite/setup.py -p ${py} -r main -b ${branch} -e ${env_name}
./suite/setup.py -p ${py} -r ctrl -b ${branch} -e ${env_name}
./suite/setup.py -p ${py} -r main_vs_ctrl -b ${branch} -e ${env_name}
./suite/setup.py -p ${py} -r no_polar_regions -b ${branch} --no_polar_regions -e ${env_name}
./suite/setup.py -p ${py} -r mesh_rename -b ${branch} -e ${env_name}

# submit the jobs
cd ${machine}_test_suite

main_py=3.11
cd main_py${main_py}
echo main_py${main_py}
RES=$(sbatch job_script.bash)
cd ..

cd main_vs_ctrl
echo main_vs_ctrl
sbatch --dependency=afterok:${RES##* } --kill-on-invalid-dep=yes job_script.bash
cd ..

for run in wc_defaults moc_am no_ncclimo no_polar_regions \
    mesh_rename
do
    cd ${run}
    echo ${run}
    sbatch job_script.bash
    cd ..
done

cd ..

