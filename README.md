# Transferring Causal Effects using Proxies
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
