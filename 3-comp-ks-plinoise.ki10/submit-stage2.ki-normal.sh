#!/bin/bash
#SBATCH --job-name=ki-ks-norm
#SBATCH -p wzhcnormal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --constraint=32core
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --output=ks_ki_norm_%j.out
#SBATCH --error=ks_ki_norm_%j.err

RESULTS_DIR="${SLURM_SUBMIT_DIR}/job_${SLURM_JOB_ID}"
LOCAL_TMP="/dev/shm/${USER}_job_${SLURM_JOB_ID}"
mkdir -p "${RESULTS_DIR}"
cleanup() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] triggering the cleanup mechanism..."
    if [ -d "${LOCAL_TMP}" ]; then
        rm -rf "${LOCAL_TMP}"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] temporary directory successfully cleared: ${LOCAL_TMP}"
    fi
}
trap cleanup EXIT SIGINT SIGTERM SIGKILL

mkdir -p "${LOCAL_TMP}"
cd "${LOCAL_TMP}" || exit 1
echo "[$(date '+%Y-%m-%d %H:%M:%S')] current working directory has been switched to the local RAM disk: $(pwd)"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] start executing the computation task..."

module purge
module load apps/anaconda3/2023.09

source activate r_env

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cp "${SLURM_SUBMIT_DIR}/res_ki_ki10.rda" .
cp "${SLURM_SUBMIT_DIR}/ki-stage2-ks.normal.R" .
Rscript ki-stage2-ks.normal.R

echo "Running computation on $HOSTNAME..." > output.txt
sleep 10  
echo "[$(date '+%Y-%m-%d %H:%M:%S')] calculation complete; writing results back to network storage..."
cp -r "${LOCAL_TMP}"/* "${RESULTS_DIR}/"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] data synchronization complete; results saved to: ${RESULTS_DIR}"
