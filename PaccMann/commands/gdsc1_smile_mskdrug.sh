#!/bin/bash 

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem 51200
#SBATCH -n 1
#SBATCH -t 200:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=qhz@email.unc.edu
#SBATCH --array=1-10

source paccmann_predict/bin/activate

python /nas/longleaf/home/qhz/paccmann_predictor/examples/IC50/train_paccmann.py \
gdsc_old/mask_drug/maskdrug_train_${SLURM_ARRAY_TASK_ID}.csv \
gdsc_old/mask_drug/maskdrug_test_${SLURM_ARRAY_TASK_ID}.csv \
gdsc_old/gdsc_gene_exp.csv \
gdsc_old/gdsc_smile.smi \
data/2128_genes.pkl \
single_pytorch_model/smiles_language \
gdsc1_maskdrug \
single_pytorch_model/model_params_50.json \
"gdsc1_maskdrug_${SLURM_ARRAY_TASK_ID}"
