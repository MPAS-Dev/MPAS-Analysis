#!/usr/bin/env bash

set -e

machine=compy

export HDF5_USE_FILE_LOCKING=FALSE

source ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate base

branch=$(git symbolic-ref --short HEAD)

conda update -y conda conda-build
rm -rf ${HOME}/miniconda3/conda-bld
conda build ci/recipe

# create the test conda envs
for py in 3.7 3.8
do
    env=test_mpas_analysis_py${py}
    conda remove -y --all -n ${env}
    conda create -y -n ${env} --use-local python=${py} mpas-analysis sphinx \
        mock sphinx_rtd_theme "tabulate>=0.8.2" m2r pytest
    conda activate ${env}
    pytest
    conda deactivate
done

# create another env for testing xarray master branch
py=3.8
env=test_mpas_analysis_xarray_master
conda create --yes --quiet --name ${env} --use-local python=${py} \
    mpas-analysis pytest
conda activate ${env}
pip install git+https://github.com/pydata/xarray.git
pytest
conda deactivate

# move to a subdirectory so we use the conda package, not the local package
rm -rf ${machine}_test_suite
mkdir ${machine}_test_suite

cd ${machine}_test_suite

template_path=../configs/${machine}/test_suite
job_template_path=${template_path}

for py in 3.7 3.8
do
    env=test_mpas_analysis_py${py}
    run=main_py${py}
    config=${run}.cfg
    mkdir ${run}
    job=${run}/job_script.bash
    sed "s/baseline/${machine}\/${branch}\/py${py}/g" ${template_path}/main.cfg > ${config}
    sed -e "s/main.cfg/${config}/g" -e "s/test_env/${env}/g" \
         ${job_template_path}/job_script.bash > ${job}
done


py=3.8
env=test_mpas_analysis_py${py}

run=wc_defaults
config=${run}.cfg
mkdir ${run}
job=${run}/job_script.bash
sed "s/baseline/${machine}\/${branch}\/${run}/g" ${template_path}/${config} > ${config}
sed -e "s/no_polar_regions.cfg/${config}/g" -e "s/test_env/${env}/g" \
     ${job_template_path}/job_script_no_polar_regions.bash > ${job}

run=no_ncclimo
config=${run}.cfg
mkdir ${run}
job=${run}/job_script.bash
sed "s/baseline/${machine}\/${branch}\/${run}/g" ${template_path}/${config} > ${config}
sed -e "s/main.cfg/${config}/g" -e "s/test_env/${env}/g" \
     ${job_template_path}/job_script.bash > ${job}

run=ctrl
config=${run}.cfg
sed "s/baseline/${machine}\/${branch}\/py${py}/g" ${template_path}/${config} > ${config}

run=main_vs_ctrl
config=${run}.cfg
mkdir ${run}
job=${run}/job_script.bash
sed "s/baseline/${machine}\/${branch}\/${run}/g" ${template_path}/${config} > ${config}
sed -e "s/main.cfg/${config}/g" -e "s/test_env/${env}/g" \
     ${job_template_path}/job_script.bash > ${job}

run=no_polar_regions
config=${run}.cfg
mkdir ${run}
job=${run}/job_script.bash
sed "s/baseline/${machine}\/${branch}\/${run}/g" ${template_path}/main.cfg > ${config}
sed -e "s/test_env/${env}/g" \
     ${job_template_path}/job_script_no_polar_regions.bash > ${job}

run=mesh_rename
config=${run}.cfg
mkdir ${run}
job=${run}/job_script.bash
sed "s/baseline/${machine}\/${branch}\/${run}/g" ${template_path}/${config} > ${config}
sed -e "s/main.cfg/${config}/g" -e "s/test_env/${env}/g" \
     ${job_template_path}/job_script.bash > ${job}

env=test_mpas_analysis_xarray_master
run=xarray_master
config=${run}.cfg
mkdir ${run}
job=${run}/job_script.bash
sed "s/baseline/${machine}\/${branch}\/${run}/g" ${template_path}/main.cfg > ${config}
sed -e "s/main.cfg/${config}/g" -e "s/test_env/${env}/g" \
     ${job_template_path}/job_script.bash > ${job}


# submit the jobs
for run in main_py3.7 wc_defaults no_ncclimo no_polar_regions \
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

