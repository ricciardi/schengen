# schengen
Code and data for the paper "Retrospective causal inference via elapsed time-weighted
matrix completion, with an evaluation of the effect of the
Schengen Area on the labour market of border regions"

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

0. package-list.R
1. schengen_MCM_data.R
2. Schengen_MCM_covars.R 
	* MCEst.R
	* PolitisWhite.R
	* MCEstBoot.R

3. mc-plot.R
	* TsPlot.R

4. mc-placebo_covars.R
5. mc-placebo-plot.R 