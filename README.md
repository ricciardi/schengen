# schengen
Code and data for the paper "Retrospective causal inference via matrix completion, with an evaluation of the effect of European integration on labour market outcomes"

Prerequsites
------

* **R** >= 3.5.0 (tested on 3.6.0 and 3.6.1)

Install the forked MCPanel repo:
```R
install.packages("devtools")
library(devtools) 
install_github("jvpoulos/MCPanel")
```

Create folders to store outputs:

```bash
mkdir data
mkdir results
mkdir plots
```

Run order
------

0. package-list.R # required packages
1. schengen_MCM_data.R # prepare data for analyses
2. Schengen_MCM_covars.R # estimates with covariates 
  (Schengen_MCM.R # estimates without covariates)
	* MCEst.R
	* PolitisWhite.R
	* MCEstBoot.R
	* MCEstBootTraj.R

3. mc-plot.R # plot estimates
	* TsPlot.R

4. schengen-rmse-placebo.R # placebo tests: compare different estimators in terms of RMSE (without covariates)
5. schengen-rmse-placebo-plot.R  # plot placebo test results

6. mc-placebo.R # placebo test of null hypothesis
7. mc-placebo-plot.R # placebo test p-values