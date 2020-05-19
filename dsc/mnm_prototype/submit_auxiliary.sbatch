#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=dsc-submitter
#SBATCH --output=dsc-log-%J.out
#SBATCH --error=dsc-log-%J.err

# 1. extract data from `simulate_only` pipeline in DSC 
#for m in gtex_mixture_identity; do 
#    sos run analysis/20200502_Prepare_ED_prior.ipynb extract --model $m -c midway2.yml -q midway2 -s build
#done
# 2. Run FLASH for data extracted
#for m in artificial_mixture_identity gtex_mixture_identity; do
#sos run ~/GIT/gtexresults/workflows/mashr_flashr_workflow.ipynb flash \
#    --cwd /project2/mstephens/gaow/mvarbvs/dsc/mnm_prototype/mnm_sumstats/ \
#    --data /project2/mstephens/gaow/mvarbvs/dsc/mnm_prototype/mnm_sumstats/$m.rds \
#    --effect-model EE -c midway2.yml -q midway2
#done
# 3. Run ED for artificial mixture
# sos run analysis/20200502_Prepare_ED_prior.ipynb ed --model artificial_mixture_identity -c midway2.yml -q midway2
# 4. Run ED for GTEx mixture
sos run analysis/20200502_Prepare_ED_prior.ipynb ed --model gtex_mixture_identity -c midway2.yml -q midway2
