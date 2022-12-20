#!/usr/bin/env bash

set -e

conda_base=$(dirname $(dirname $CONDA_EXE))
source $conda_base/etc/profile.d/conda.sh

main_py=3.10
alt_py=3.9

export HDF5_USE_FILE_LOCKING=FALSE

branch=$(git symbolic-ref --short HEAD)

conda update -y conda conda-build mamba boa
conda mambabuild ci/recipe

# create the test conda envs
for py in ${main_py} ${alt_py}
do
    env=test_mpas_analysis_py${py}
    mamba create -y -n ${env} --use-local python=${py} mpas-analysis sphinx \
        mock sphinx_rtd_theme "tabulate>=0.8.2" m2r2 "mistune<2" pytest \
	      "mache>=1.1.2" jinja2
    conda activate ${env}
    pytest
    conda deactivate
done

# create another env for testing xarray main branch
py=${main_py}
env=test_mpas_analysis_xarray_main
mamba create --yes --quiet --name ${env} --use-local python=${py} \
    mpas-analysis pytest
conda activate ${env}
pip install git+https://github.com/pydata/xarray.git
pytest
conda deactivate

# test building the docs
py=${main_py}
conda activate test_mpas_analysis_py${py}
cd docs
make clean
make html
cd ..

machine=$(python -c "from mache import discover_machine; print(discover_machine())")

./suite/setup.py -p ${py} -r main_py${py} -b ${branch} --copy_docs --clean
./suite/setup.py -p ${py} -r wc_defaults -b ${branch} --no_polar_regions
./suite/setup.py -p ${py} -r no_ncclimo -b ${branch}
./suite/setup.py -p ${py} -r ctrl -b ${branch}
./suite/setup.py -p ${py} -r main_vs_ctrl -b ${branch}
./suite/setup.py -p ${py} -r no_polar_regions -b ${branch} --no_polar_regions
./suite/setup.py -p ${py} -r mesh_rename -b ${branch}
./suite/setup.py -p ${py} -r xarray_main -b ${branch} -e test_mpas_analysis_xarray_main
conda deactivate

py=${alt_py}
conda activate test_mpas_analysis_py${py}
./suite/setup.py -p ${py} -r main_py${py} -b ${branch}
conda deactivate

# submit the jobs
cd ${machine}_test_suite

cd main_py${main_py}
echo main_py${main_py}
RES=$(sbatch job_script.bash)
cd ..

cd main_vs_ctrl
echo main_vs_ctrl
sbatch --dependency=afterok:${RES##* } job_script.bash --kill-on-invalid-dep=yes
cd ..

for run in main_py${alt_py} wc_defaults no_ncclimo no_polar_regions \
    mesh_rename xarray_main
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
   py=${main_py}
   conda activate test_mpas_analysis_py${py}
  ./suite/setup.py -p ${py} -r QU480 -b ${branch}
  cd ${machine}_test_suite/QU480
  echo QU480
  sbatch job_script.bash
  cd ../..
fi
