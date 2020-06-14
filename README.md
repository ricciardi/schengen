# schengen
 

Prerequsites
------

* **R** >= 3.5.0 (tested on 3.6.1)
 * dplyr, readstata13, glmnet, boot, doParallel

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
3. mc-plot.R

4. mc-placebo.R
5. mc-placebo-plot.R


Example usage:
------

```R
library(MCPanel)
estimated_obj <- mcnnm_cv(M, mask, W, to_estimate_u = 0, to_estimate_v = 0, num_lam_L = 40)
                  
best_lam_L <- estimated_obj$best_lambda
estimated_mat <- estimated_obj$L

```