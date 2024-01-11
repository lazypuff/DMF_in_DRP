#!/bin/bash 

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem 12800
#SBATCH -n 1
#SBATCH --mail-type=end
#SBATCH --mail-user=qhz@email.unc.edu
#SBATCH --array=1-10

source paccmann_predict/bin/activate

python /nas/longleaf/home/qhz/paccmann_predictor/examples/IC50/test_paccmann.py \
gdsc_old/mask_comb/maskcomb_valid_${SLURM_ARRAY_TASK_ID}.csv \
gdsc_old/gdsc_gene_exp.csv \
gdsc_old/gdsc_smile_random.smi \
data/2128_genes.pkl \
single_pytorch_model/smiles_language \
gdsc1_maskcomb_random/gdsc1_maskcomb_random_${SLURM_ARRAY_TASK_ID}/weights/best_mse_paccmann_v2.pt \
gdsc1_mskcomb_random_${SLURM_ARRAY_TASK_ID} \
gdsc1_maskcomb_random/gdsc1_maskcomb_random_${SLURM_ARRAY_TASK_ID}/model_params.json

