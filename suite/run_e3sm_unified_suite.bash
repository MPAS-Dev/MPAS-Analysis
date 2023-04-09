#!/usr/bin/env bash

set -e

# placeholder that gets replaced
branch=test_e3sm_unified

# test building the docs
py=3.10
machine=${E3SMU_MACHINE}

./suite/setup.py -p ${py} -r main_py${py} -b ${branch} --clean
./suite/setup.py -p ${py} -r wc_defaults -b ${branch} --no_polar_regions
./suite/setup.py -p ${py} -r moc_am -b ${branch}
./suite/setup.py -p ${py} -r no_ncclimo -b ${branch}
./suite/setup.py -p ${py} -r ctrl -b ${branch}
./suite/setup.py -p ${py} -r main_vs_ctrl -b ${branch}
./suite/setup.py -p ${py} -r no_polar_regions -b ${branch} --no_polar_regions
./suite/setup.py -p ${py} -r mesh_rename -b ${branch}

# submit the jobs
cd ${machine}_test_suite

cd main_py${py}
echo main_py${py}
RES=$(sbatch job_script.bash)
cd ..

cd main_vs_ctrl
echo main_vs_ctrl
sbatch --dependency=afterok:${RES##* } --kill-on-invalid-dep=yes job_script.bash
cd ..

for run in wc_defaults moc_am no_ncclimo no_polar_regions mesh_rename
do
    cd ${run}
    echo ${run}
    sbatch job_script.bash
    cd ..
done

cd ..

# only LCRC machines have a separate QU480 run
if [[ "$machine" == "anvil" || "$machine" == "chrysalis" ]]
then
  ./suite/setup.py -p ${py} -r QU480 -b ${branch}
  cd ${machine}_test_suite/QU480
  echo QU480
  sbatch job_script.bash
  cd ../..
fi
