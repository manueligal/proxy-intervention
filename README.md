# Transferring Causal Effects using Proxies
This repository contains code for "Transferring Causal Effects using Proxies" <add link to OpenReview or arxiv> by Iglesias-Alonso, Schur, von Kügelgen, and Peters (NeurIPS 2025).

### Main files
* `BaselineComparison.R`: Comparison of the point estimates of the two proposed estimators with different baselines.
* `Bootstrap.R`: Comparison of the asymptotic and bootstrap confidence intervals in terms of coverage and median interval length.
* `ConditionNumber.R`: Comparison of the point estimates from the causal and reduced parametrisations.
* `HotelsExample.R`: Application of the model and the estimators to a real example.
* `TimeComparison.R`: Runtime analysis of the proposed estimators and baselines.

### Auxiliary files
The auxiliary files can be found in the `utils` folder.
* `misc.R`: General helper functions.
* `Lfunction.R`: Definition of the likelihood function.
* `FunctionsSampling.R`: Functions to generate the matrices that describe the model and samples from the entailed distribution.
* `FunctionsEstimation.R`: Functions to apply the estimators from the causal and reduced parametrisations.
* `FunctionsCI.R`: Function to calculate the asymptotic standard deviation in the reduced parametrisation.

If you use this code in your research, please cite:

@inproceedings{iglesias2025transferring, \
title={Transferring Causal Effects using Proxies}, \
author={Iglesias-Alonso, Manuel and Schur, Felix and {von K{"u}gelgen}, Julius and Peters, Jonas}, \
booktitle={Advances in Neural Information Processing Systems}, \
year={2025} \
}
