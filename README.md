# schengen
Code for the paper "Retrospective causal inference via matrix completion, with an evaluation of the effect of European integration on labour market outcomes."

Microdata access
------

The data used for the paper is constructed from European Labour Force Survey (ELFS) microdata, which requires access from Eurostat. 

Instructions on how to gain access to the microdata is found [here](https://ec.europa.eu/eurostat/web/microdata/european-union-labour-force-survey).

After obtaining the microdata, follow the commands in `BUILDING DATA.do` for constructing the dataset `FINAL_21.`.

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
**Note:** fitting the matrix completion model with covariates (*mcnnm_wc_cv*) is computationally expensive and will likely make a laptop crash. The code below is run on a single node with 30G RAM and 6 CPU-cores on a high-performance compute cluster.  

Create folders to store outputs:

```bash
mkdir data
mkdir results
mkdir plots
```

The file FINAL_21.dta is needed for steps 1-7 and should be in the `data/` directory. 

Run order
------

0. mc-simulation.R # simulated data experiments
	* run with command line argument  `Rscript mc-simulation.R [arg1]`, where `[arg1]` is a number specifying the simulation setting
	* mc-simulation-plot.R # plot matrix completion simulation study results

1. schengen_MCM_data.R # prepare data for analyses

2. mc-simulation-placebo.R # placebo tests experiments on cross-border worker data
	* run with command line argument  `Rscript mc-simulation-placebo.R [arg1]`, where `[arg1]` is a number specifying the simulation setting
	* mc-simulation-placebo-plot.R  # plot placebo test results

3. mc-placebo.R # placebo test of null hypothesis

4. Schengen_MCM.R # fit model with covariates 
	* mc-plot.R # plot estimates

5. Schengen_Compare.R # DID, SCM, and IFE estimates for comparison 