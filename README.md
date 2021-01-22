# schengen
Code for the paper "Retrospective causal inference via matrix completion, with an evaluation of the effect of European integration on labour market outcomes."

Microdata access
------

The data used for the paper is constructed from European Labour Force Survey (ELFS) microdata, which requires access from Eurostat. 

Instructions on how to gain access to the microdata is found [here](https://ec.europa.eu/eurostat/web/microdata/european-union-labour-force-survey).

After obtaining the microdata, follow the commands in BUILDING DATA.do for constructing the dataset FINAL_21.dta.

Prerequsites
------

* **R** >= 3.5.0 (tested on 3.6.3)

* package-list.R # required **R** packages

Install the forked MCPanel repo:
```R
install.packages("devtools")
library(devtools) 
install_github("jvpoulos/MCPanel")
```
**Note:** fitting the matrix completion model with covariates (*mcnnm_wc_cv*) is computationally expensive and will likely make a laptop crash. The code below is run on a high-performance compute cluster (384 GB RAM and 40 CPU-cores) with an allocated memory of 30G. 

Create folders to store outputs:

```bash
mkdir data
mkdir results
mkdir plots
```

Run order
------

1. schengen_MCM_data.R # prepare data for analyses
	* MCEst.R

2. Schengen_MCM_covars.R # model with covariates 
	* MCEst.R
	* PolitisWhite.R
	* MCEstBoot.R
	* MCEstBootTraj.R

3. mc-plot.R # plot estimates
	* TsPlot.R
	* TsPlotTrends.R
	* getTrends.R

4. schengen-rmse-placebo.R # placebo tests: compare different estimators in terms of RMSE (model without covariates)
5. schengen-rmse-placebo-plot.R  # plot placebo test results

6. mc-placebo.R # placebo test of null hypothesis (model without covariates)

7. Schengen_DID.R # DID and SCM estimates for comparison 
	* DIDEstBoot.R