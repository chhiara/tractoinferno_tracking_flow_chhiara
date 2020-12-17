#!/bin/bash
#SBATCH --job-name=LocalTracking
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err

#SBATCH --mail-user=philippe.poulin2@usherbrooke.ca
#SBATCH --mail-type=ALL


#SBATCH --account=def-descotea
#SBATCH --nodes=4
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=120:00:00

module load httpproxy/1.0
export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)

# Nextflow modules
module load java/1.8
module load singularity/3.6
module load nextflow

WORK_DIR=$HOME/projects/rrg-descotea/TractoInferno/derivatives/LocalTracking
SINGULARITY_IMG=$HOME/projects/rrg-descotea/containers/tractoflow-scilpy-1.0.0.sif
TRACTOFLOW_DIR=$HOME/projects/rrg-descotea/TractoInferno/derivatives/TractoFlow/output/results
FLOW_DIR=$HOME/git/tractoinferno_tracking_flow

cd "${WORK_DIR}"
srun nextflow -c "${FLOW_DIR}/nextflow.config" run "${FLOW_DIR}/main.nf" --input "${TRACTOFLOW_DIR}" -with-singularity "${SINGULARITY_IMG}" -resume -with-mpi
