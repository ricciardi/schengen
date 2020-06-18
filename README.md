# schengen
 

Prerequsites
------

* **R** >= 3.5.0 (tested on 3.6.1)
 * dplyr, readstata13, glmnet, boot, doParallel, caret

Install the forked MCPanel repo:
```R
install.packages("devtools")
install.packages("latex2exp")
library(devtools) 
install_github("jvpoulos/MCPanel")
```

Run order
------

1. schengen_MCM_data.R
2. Schengen_MCM.R
	* MCEst.R
	* PolitisWhite.R
	* MCEstBoot.R
	* ChernoTest.R
3. mc-plot.R

4. mc-placebo.R
5. mc-placebo-plot.R