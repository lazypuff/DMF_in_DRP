# Drug sensitivity prediction with PaccMann

- This folder contains pretrained weights of our PaccMann model as detailed in the paper [_Toward explainable anticancer compound sensitivity prediction via multimodal attention-based convolutional encoders_](https://doi.org/10.1021/acs.molpharmaceut.9b00520) (*Molecular Pharmaceutics*, 2019).

- Please note that the original models (used to compile the results in the paper) were developed and evaluated in `tensorflow`. Please see [here](https://github.com/drugilsberg/paccmann) for the original code.

- This model is the `pytorch` implementation of the best PaccMann architecture (multiscale convolutional encoder). Note that there is one architectural advancement compared to the model detailed in the original paper. Instead of using self-attention on the gene expression profiles we now use context-attention, just like for the SMILES. The rest is the same.

- 


- PaccMann acronyms "_Prediction of AntiCancer Compound sensitivity with Multi-modal Attention-based Neural Networks_".

# Evaluating PaccMann on your own data

## Installation

To re-evaluate this 

```sh
git clone https://github.com/PaccMann/paccmann_predictor
cd paccmann_predictor
git checkout dev
conda env create -f examples/IC50/conda.yml
conda activate paccmann_predictor
pip install -e .
pip install torch==1.7.0 # since the model in this example was trained with that version
```

## Example usage


```console
(paccmann_predictor) $ python examples/IC50/test_paccmann.py -h
usage: test_paccmann.py [-h]
                        test_sensitivity_filepath gep_filepath smi_filepath
                        gene_filepath smiles_language_filepath model_filepath
                        predictions_filepath params_filepath

positional arguments:
  test_sensitivity_filepath
                        Path to the drug sensitivity (IC50) data.
  gep_filepath          Path to the gene expression profile data.
  smi_filepath          Path to the SMILES data.
  gene_filepath         Path to a pickle object containing list of genes.
  smiles_language_filepath
                        Path to a folder with SMILES language .json files.
  model_filepath        Path to the stored model.
  predictions_filepath  Path to the predictions.
  params_filepath       Path to the parameter file.

optional arguments:
  -h, --help            show this help message and exit
```

For example, assuming that you downloaded this model in a directory called `single_pytorch_model`, the data from https://ibm.box.com/v/paccmann-pytoda-data in folders `data` and `splitted_data` the following command should work:
```console
(paccmann_predictor) $ python examples/IC50/test_paccmann.py \
splitted_data/gdsc_cell_line_ic50_test_fraction_0.1_id_997_seed_42.csv \
data/gene_expression/gdsc-rnaseq_gene-expression.csv \
data/smiles/gdsc.smi \
data/2128_genes.pkl \
single_pytorch_model/smiles_language \
single_pytorch_model/weights/best_mse_paccmann_v2.pt \
results \
single_pytorch_model/model_params.json
```


## References

If you use this model in your projects, please cite the following:

```bib
@article{manica2019paccmann,
  title={Toward explainable anticancer compound sensitivity prediction via multimodal attention-based convolutional encoders},
  author={Manica, Matteo and Oskooei, Ali and Born, Jannis and Subramanian, Vigneshwari and S{\'a}ez-Rodr{\'\i}guez, Julio and Rodríguez Martínez, María},
  journal={Molecular pharmaceutics},
  volume={16},
  number={12},
  pages={4797--4806},
  year={2019},
  publisher={ACS Publications},
  doi = {10.1021/acs.molpharmaceut.9b00520},
  note = {PMID: 31618586}
}
@article{cadow2020paccmann,
  title={PaccMann: a web service for interpretable anticancer compound sensitivity prediction},
  author={Cadow, Joris and Born, Jannis and Manica, Matteo and Oskooei, Ali and Rodr{\'\i}guez Mart{\'\i}nez, Mar{\'\i}a},
  journal={Nucleic acids research},
  volume={48},
  number={W1},
  pages={W502--W508},
  year={2020},
  publisher={Oxford University Press}
}
```
