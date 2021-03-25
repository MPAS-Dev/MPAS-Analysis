#!/usr/bin/env bash

set -e

branch=$(git symbolic-ref --short HEAD)

export HDF5_USE_FILE_LOCKING=FALSE

source ${HOME}/miniconda3/etc/profile.d/conda.sh

conda activate base
conda update -y conda conda-build
rm -rf ${HOME}/miniconda3/conda-bld

# create the test conda envs
for py in 3.7 3.8
do
    env=test_mpas_analysis_py${py}
    conda build -m ci/python${py}.yaml ci/recipe
    conda remove -y --all -n ${env}
    conda create -y -n ${env} --use-local python=${py} mpas-analysis sphinx \
        mock sphinx_rtd_theme "tabulate>=0.8.2" m2r pytest
    conda activate ${env}
    pytest
    conda deactivate
done

# create another env for testing xarray master branch
env=test_mpas_analysis_xarray_master
conda create --yes --quiet --name ${env} --use-local python=${py} \
    mpas-analysis pytest
conda activate ${env}
pip install git+https://github.com/pydata/xarray.git
pytest
conda deactivate

# test building the docs
py=3.8
conda activate test_mpas_analysis_py${py}
cd docs
make clean
make html
rm -rf /lcrc/group/acme/public_html/diagnostic_output/ac.xylar/analysis_testing/${branch}/docs
mkdir -p /lcrc/group/acme/public_html/diagnostic_output/ac.xylar/analysis_testing/${branch}/
cp -r _build/html /lcrc/group/acme/public_html/diagnostic_output/ac.xylar/analysis_testing/${branch}/docs
cd ..
conda deactivate

# move to a subdirectory so we use the conda package, not the local package
rm -rf anvil_test_suite
mkdir anvil_test_suite

cd anvil_test_suite

template_path=../configs/anvil/test_suite

for py in 3.7 3.8
do
    env=test_mpas_analysis_py${py}
    run=main_py${py}
    config=${run}.cfg
    mkdir ${run}
    job=${run}/job_script.bash
    sed "s/baseline/${branch}\/py${py}/g" ${template_path}/main.cfg > ${config}
    sed -e "s/main.cfg/${config}/g" -e "s/test_env/${env}/g" \
         ${template_path}/job_script.bash > ${job}
done


py=3.8
env=test_mpas_analysis_py${py}

run=wc_defaults
config=${run}.cfg
mkdir ${run}
job=${run}/job_script.bash
sed "s/baseline/${branch}\/${run}/g" ${template_path}/${config} > ${config}
sed -e "s/no_polar_regions.cfg/${config}/g" -e "s/test_env/${env}/g" \
     ${template_path}/job_script_no_polar_regions.bash > ${job}

run=no_ncclimo
config=${run}.cfg
mkdir ${run}
job=${run}/job_script.bash
sed "s/baseline/${branch}\/${run}/g" ${template_path}/${config} > ${config}
sed -e "s/main.cfg/${config}/g" -e "s/test_env/${env}/g" \
     ${template_path}/job_script.bash > ${job}

run=ctrl
config=${run}.cfg
sed "s/baseline/${branch}\/py${py}/g" ${template_path}/${config} > ${config}

run=main_vs_ctrl
config=${run}.cfg
mkdir ${run}
job=${run}/job_script.bash
sed "s/baseline/${branch}\/${run}/g" ${template_path}/${config} > ${config}
sed -e "s/main.cfg/${config}/g" -e "s/test_env/${env}/g" \
     ${template_path}/job_script.bash > ${job}

run=no_polar_regions
config=${run}.cfg
mkdir ${run}
job=${run}/job_script.bash
sed "s/baseline/${branch}\/${run}/g" ${template_path}/main.cfg > ${config}
sed -e "s/test_env/${env}/g" \
     ${template_path}/job_script_no_polar_regions.bash > ${job}

run=QU480
config=${run}.cfg
mkdir ${run}
job=${run}/job_script.bash
sed "s/baseline/${branch}\/${run}/g" ${template_path}/${config} > ${config}
sed -e "s/main.cfg/${config}/g" -e "s/test_env/${env}/g" \
     ${template_path}/job_script.bash > ${job}

run=mesh_rename
config=${run}.cfg
mkdir ${run}
job=${run}/job_script.bash
sed "s/baseline/${branch}\/${run}/g" ${template_path}/${config} > ${config}
sed -e "s/main.cfg/${config}/g" -e "s/test_env/${env}/g" \
     ${template_path}/job_script.bash > ${job}

env=test_mpas_analysis_xarray_master
run=xarray_master
config=${run}.cfg
mkdir ${run}
job=${run}/job_script.bash
sed "s/baseline/${branch}\/${run}/g" ${template_path}/main.cfg > ${config}
sed -e "s/main.cfg/${config}/g" -e "s/test_env/${env}/g" \
     ${template_path}/job_script.bash > ${job}


# submit the jobs
for run in main_py3.7 wc_defaults no_ncclimo no_polar_regions QU480 \
    mesh_rename xarray_master
do
    cd ${run}
    sbatch job_script.bash
    cd ..
done

cd main_py3.8
RES=$(sbatch job_script.bash)
cd ..

cd main_vs_ctrl
sbatch --dependency=afterok:${RES##* } job_script.bash
cd ..

cd ..

