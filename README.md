# schengen
Code and data for the paper "Retrospective causal inference via matrix completion, with an evaluation of the effect of European integration on labour market outcomes"

Prerequsites
------

* **R** >= 3.5.0 (tested on 3.6.0 and 3.6.1)

Install the forked MCPanel repo and NNLM:
```R
install.packages("devtools")
library(devtools) 
install_github("jvpoulos/MCPanel")
install_github('linxihui/NNLM')
```

Create folders to store outputs:

```bash
mkdir data
mkdir results
mkdir plots
```

Run order
------

0. package-list.R # required **R** packages

1. schengen_MCM_data.R # prepare data for analyses
	* MCEst.R

2. Schengen_MCM_covars.R # model with covariates 
	* MCEst.R
	* PolitisWhite.R
	* MCEstBoot.R
	* MCEstBootTraj.R
	
3. Schengen_MCM.R # model without covariates

4. mc-plot.R # plot estimates
	* TsPlot.R

5. schengen-rmse-placebo.R # placebo tests: compare different estimators in terms of RMSE (model without covariates)
6. schengen-rmse-placebo-plot.R  # plot placebo test results

7. mc-placebo.R # placebo test of null hypothesis (model without covariates)

8. Schengen_DID.R # DID, SCM, or SCM-L1 estimates for comparison 
	* DIDEstBoot.R