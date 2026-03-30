#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=2:00:00
{{ sbatch }}
#SBATCH --job-name=mpas_analysis
#SBATCH --output=mpas_analysis.o%j
#SBATCH --error=mpas_analysis.e%j

set -e

{% if use_e3sm_unified %}
source {{ e3sm_unified_script }}

echo E3SM-Unified: {{ e3sm_unified_script }}
{% elif pixi_env %}
export HDF5_USE_FILE_LOCKING=FALSE
export E3SMU_MACHINE={{ machine }}

run_mpas_analysis() {
    pixi run --manifest-path ../../pixi.toml -e {{ pixi_env }} mpas_analysis "$@"
}

echo pixi env: {{ pixi_env }}
{% else %}
source {{ conda_base }}/etc/profile.d/conda.sh
conda activate {{ conda_env }}
export HDF5_USE_FILE_LOCKING=FALSE
export E3SMU_MACHINE={{ machine }}

run_mpas_analysis() {
    mpas_analysis "$@"
}

echo env: {{ conda_env }}
{% endif %}
echo configs: {{ flags }} {{ config }}

run_mpas_analysis --list
run_mpas_analysis --plot_colormaps
run_mpas_analysis --setup_only {{ flags }} {{ config }}
run_mpas_analysis --purge {{ flags }} {{ config }} --verbose
run_mpas_analysis --html_only {{ flags }} {{ config }}

chmod ugo+rx {{ html_base }}/{{ out_common_dir }}
chmod -R ugo+rX {{ html_base }}/{{ out_subdir }}
