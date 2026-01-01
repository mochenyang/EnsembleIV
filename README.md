# EnsembleIV and ForestIV
An R package that implements [EnsembleIV](https://arxiv.org/abs/2303.02820) and [ForestIV](https://pubsonline.informs.org/doi/abs/10.1287/ijds.2022.0019).

To install, run 
```r
devtools::install_github("mochenyang/EnsembleIV")
```

Please consider citing
```
@article{burtch2023ensembleiv,
  title={EnsembleIV: Creating Instrumental Variables from Ensemble Learners for Robust Statistical Inference},
  author={Burtch, Gordon and McFowland III, Edward and Yang, Mochen and Adomavicius, Gediminas},
  journal={Management Science, forthcoming},
  year={2025}
}
```

## Replication Materials for EnsembleIV paper

The ```replication``` folder contains necessary replication code and data for the EnsembleIV paper. Specifically,

- The ```simulations``` folder contains code and data to replicate the simulation results on both Bike Sharing data and Bank Marketing data, as well as benchmarking and sensitivity analyses on the same datasets.
- The ```FongTyler Replication``` folder contains (author-provided) implementation of Fong and Tyler's GMM method (used in benchmarking analysis).
- The ```process results.R``` file contains R code to post-process results and generate tables / figures in the paper.
- To access data / code for the ```facebook``` study, please contact Mochen Yang (yang3653@umn.edu).