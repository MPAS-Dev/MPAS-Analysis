#!/usr/bin/env bash

set -euo pipefail

main_py=3.13
alt_py=3.12
mode=package
pixi_env=${PIXI_ENVIRONMENT_NAME:-default}

usage() {
    cat <<EOF
Usage: ./suite/run_suite.bash [--dev] [--e3sm-unified] [--pixi-env ENV]

Without flags, run the package-build suite in fresh conda environments.
Use --dev for the Pixi-based developer workflow.
Use --e3sm-unified when running inside an E3SM-Unified deployment.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dev)
            mode=dev
            ;;
        --e3sm-unified)
            mode=e3sm-unified
            ;;
        --pixi-env)
            if [[ $# -lt 2 ]]; then
                echo "--pixi-env requires an environment name" >&2
                exit 1
            fi
            shift
            pixi_env=$1
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
    shift
done

export HDF5_USE_FILE_LOCKING=FALSE

branch=$(git symbolic-ref --short HEAD)

setup_run() {
    local py="$1"
    local run="$2"
    shift 2
    "${setup_cmd[@]}" -p "${py}" -r "${run}" -b "${branch}" "$@"
}

submit_jobs() {
    local machine="$1"
    local primary_py="$2"
    shift 2

    cd "${machine}_test_suite"

    cd "main_py${primary_py}"
    echo "main_py${primary_py}"
    RES=$(sbatch job_script.bash)
    cd ..

    cd main_vs_ctrl
    echo main_vs_ctrl
    sbatch --dependency=afterok:${RES##* } --kill-on-invalid-dep=yes \
        job_script.bash
    cd ..

    for run in "$@"; do
        cd "${run}"
        echo "${run}"
        sbatch job_script.bash
        cd ..
    done

    cd ..
}

if [[ "${mode}" == "dev" ]]; then
    if ! command -v pixi >/dev/null 2>&1; then
        echo "pixi is required for --dev" >&2
        exit 1
    fi

    docs_cmd=(pixi run -e "${pixi_env}" bash -lc \
        "cd docs && DOCS_VERSION=test make clean versioned-html")
    setup_cmd=(pixi run -e "${pixi_env}" python ./suite/setup.py \
        --pixi-env "${pixi_env}")

    "${docs_cmd[@]}"

    machine=$(pixi run -e "${pixi_env}" python -c \
        "from mache import discover_machine; print(discover_machine())")
    py=$(pixi run -e "${pixi_env}" python -c \
        'import sys; print(f"{sys.version_info[0]}.{sys.version_info[1]}")')

    setup_run "${py}" "main_py${py}" --copy_docs --clean
    setup_run "${py}" wc_defaults --no_polar_regions
    setup_run "${py}" moc_am
    setup_run "${py}" no_ncclimo
    setup_run "${py}" ctrl
    setup_run "${py}" main_vs_ctrl
    setup_run "${py}" no_polar_regions --no_polar_regions
    setup_run "${py}" mesh_rename

    submit_jobs "${machine}" "${py}" \
        wc_defaults moc_am no_ncclimo no_polar_regions mesh_rename
    exit 0
fi

if [[ "${mode}" == "e3sm-unified" ]]; then
    setup_cmd=(python ./suite/setup.py)
    py=$(python -c 'import sys; print(f"{sys.version_info[0]}.{sys.version_info[1]}")')
    machine=${E3SMU_MACHINE}
    branch=test_e3sm_unified

    setup_run "${py}" "main_py${py}" --clean
    setup_run "${py}" wc_defaults --no_polar_regions
    setup_run "${py}" moc_am
    setup_run "${py}" no_ncclimo
    setup_run "${py}" ctrl
    setup_run "${py}" main_vs_ctrl
    setup_run "${py}" no_polar_regions --no_polar_regions
    setup_run "${py}" mesh_rename

    submit_jobs "${machine}" "${py}" \
        wc_defaults moc_am no_ncclimo no_polar_regions mesh_rename
    exit 0
fi

conda_base=$(dirname "$(dirname "${CONDA_EXE}")")
source "${conda_base}/etc/profile.d/conda.sh"

conda update -y conda conda-build
conda build ci/recipe

for py in "${main_py}" "${alt_py}"; do
    env="test_mpas_analysis_py${py}"
    conda create -y -n "${env}" --use-local python="${py}" mpas-analysis \
        sphinx mock sphinx_rtd_theme "tabulate>=0.8.2" \
        "sphinx-mdinclude>=0.6.2" pytest "mache>=1.11.0" \
        "esmf=*=mpi_mpich_*" jinja2
    conda activate "${env}"
    pytest
    conda deactivate
done

py=${main_py}
env=test_mpas_analysis_xarray_main
conda create --yes --quiet --name "${env}" --use-local python="${py}" \
    mpas-analysis pytest
conda activate "${env}"
pip install git+https://github.com/pydata/xarray.git
pytest
conda deactivate

conda activate "test_mpas_analysis_py${py}"
(
    cd docs
    DOCS_VERSION=test make clean versioned-html
)

machine=$(python -c "from mache import discover_machine; print(discover_machine())")
setup_cmd=(./suite/setup.py)

setup_run "${py}" "main_py${py}" --copy_docs --clean
setup_run "${py}" wc_defaults --no_polar_regions
setup_run "${py}" moc_am
setup_run "${py}" no_ncclimo
setup_run "${py}" ctrl
setup_run "${py}" main_vs_ctrl
setup_run "${py}" no_polar_regions --no_polar_regions
setup_run "${py}" mesh_rename
setup_run "${py}" xarray_main -e test_mpas_analysis_xarray_main
conda deactivate

py=${alt_py}
conda activate "test_mpas_analysis_py${py}"
setup_run "${py}" "main_py${py}"
conda deactivate

submit_jobs "${machine}" "${main_py}" \
    "main_py${alt_py}" wc_defaults moc_am no_ncclimo no_polar_regions \
    mesh_rename xarray_main
