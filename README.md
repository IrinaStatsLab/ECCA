# ECCA (Exponential canonical correlation analysis)
A repository to apply ECCA method and reproduce the results from the following manuscript:

 * Yuan D, Zhang Y, Guo S, Wang W, and Gaynanova I 'Exponential canonical correlation analysis with orthogonal variation'.

To reproduce the results, clone this project locally and open **ECCA_code.Rproj** in RStudio to ensure all relative paths in scripts work as expected.

## 1. Getting started
### Main functions
*Function* folder contains function scripts needed for implementation of ECCA. 

**exponentialCCA.R** - *ecca_v2* is the main function that performs ECCA, which estimates every joint and individual parameter in ECCA model given the ranks.

**AnalyticalUpdate.R** - Includes functions that performs closed-form updates for Gaussian cases.

**exp_family.R** - This file includes the preliminary processing functions for exponential families. Part of this file is a modification of [generalizedPCA](https://github.com/andland/generalizedPCA/blob/master/R/exponential_family_functions.R). 

**newton_method.R** - Contains functions of performing damped Newton's method for updating loading matrices.

**optU.R** - Contains functions of performing SOC algorithm for updating joint score matrices.

**SOC.R** - Contains functions of performing SOC algorithm for updating individual score matrices.

### Example
To illustrate application of ECCA, we will simulate data based on ECCA model using scripts from **data_generator.R** (located in **Simulation** folder)

```{r}
rm(list = ls())
source("Simulation/data_generator.R")
source("Function/exp_family.R")
set.seed(37)

# Set parameters
n = 10
p1 = 8
p2 = 7
p = c(p1, p2)
family1 = 'gaussian'
family2 = 'binomial'
lambdas = c(1, 0.7) # covariance between joint structures
trials = 100

# Set ranks (joint and individuals)
r0 = 2
r1 = 1
r2 = 2
SNR = 5

# Generate data
sampledata = dataGenerator(r0 = r0, r1 = r1, r2 = r2, n = n, p = p, lambda = lambdas, family1, family2, SNR = SNR, size = trials, main_effect = TRUE)

# Get true signal
theta1_true = sampledata$theta1
theta2_true = sampledata$theta2
Z1_true = sampledata$Z1
Z2_true = sampledata$Z2
U1_true = sampledata$U1
U2_true = sampledata$U2

# Sanity check
print(crossprod(U1_true, U2_true))
```

We then apply ECCA and measure estimation accuracy.
```{r}
# Source all necessary functions to apply ECCA
source('Function/exponentialCCA.R')
source('Function/optU.R')
source('Function/SOC.R')
source('Function/newton_method.R')
source('Function/AnalyticalUpdate.R')

# Apply ECCA
ECCA_result = ecca_v2(x1 = sampledata$X1, x2 = sampledata$X2, r0 = r0, r1 = r1, r2 = r2, family1 = family1, family2 = family2, trials = trials)

# Record ECCA results
U1_ecca = ECCA_result$U1
U2_ecca = ECCA_result$U2
Z1_ecca = ECCA_result$Z1
Z2_ecca = ECCA_result$Z2 

# Get estimated natural parameters
one_vec = rep(1, n)
theta1_ecca = tcrossprod(one_vec, ECCA_result$mu1) + tcrossprod(U1_ecca, ECCA_result$V1) + tcrossprod(Z1_ecca, ECCA_result$A1)
theta2_ecca = tcrossprod(one_vec, ECCA_result$mu2) + tcrossprod(U2_ecca, ECCA_result$V2) + tcrossprod(Z2_ecca, ECCA_result$A2)

# Check the relative error between estimated signal of DMMD with true signal 
Fnorm(theta1_ecca - theta1_true)^2/Fnorm(theta1_true)^2
Fnorm(theta2_ecca - theta2_true)^2/Fnorm(theta2_true)^2

# Check the chordal distance between estimated joint space
1/sqrt(2)*Fnorm(projection(U1_ecca) - projection(U1_true))
1/sqrt(2)*Fnorm(projection(U2_ecca) - projection(U2_true))

# Check the estimated canonical correlations
angles = angle_cal(U1_ecca, U2_ecca)$angle
cos(angles)
```

## 2. Simulations
*Simulation* folder contains the necessary files to reproduce the results. 

**data_generator.R** - Generates simulated data, including Gaussian and Binomial proportion data.

**GG_Simulation_SNR5**, **GB_Simulation_SNR5**, **BB_Simulation_Trials100** folders within *Simulation*:
Each folder represents one setting in the paper. To be specific: **GG_Simulation_SNR5** represents Gaussian-Gaussian setting 1; **GB_Simulation_SNR5** represents Gaussian-Binomial setting 2; **BB_Simulation_Trials** represents Binomial-Binomial setting 3.

- *data_generation.R* within each respective folder generates simulation data for corresponding setting. The data is stored as *Data.RData*.
- *ECCA* folder within each respective folder performs ECCA algorithm on the simulated data and stores the result at *ECCA_result.RData*. Notice for setting 2 and 3, we use parallel computing. Readers might need to adjust their working directory.
- *EPCA-DCCA* folder within each respective folder performs heuristic DCCA and EPCA-DCCA methods on the simulated data and stores the result at *DCCA_result.RData*. For EPCA, please refer to [generalizedPCA Github repo](https://github.com/andland/generalizedPCA).
- *GAS* folder within each respective folder performs GAS algorithm on the simulated data and stores the result at a RData file starting with 'GAS'. For GAS, please refer to [Li, G. and Gaynanova, I. (2018). A general framework for association analysis of heterogeneous data. The Annals of Applied Statistics 12, 1700–1726.](https://www.jstor.org/stable/26542591?seq=1) and [GAS Github repo](https://github.com/reagan0323/GAS). 
- *DraftPlot.R* within each respective folder plots the results within each simulation setting and save the figures at *Figures* folder within each respective folder. These figures are not used in the paper.

**dcca.py** The python file for performing DCCA algorithm. See [Shu, Hai, Xiao Wang, and Hongtu Zhu. "D-CCA: A decomposition-based canonical correlation analysis for high-dimensional datasets." Journal of the American Statistical Association 115.529 (2020): 292-306](https://amstat.tandfonline.com/doi/full/10.1080/01621459.2018.1543599?casa_token=HA13MS9KztkAAAAA%3A1Q_j0Z1DWQ-32p83DDooAf1SxI318fE5HglIgRj1YyNpZY_Kv6BJ-0RTkIajA3t6vIA_QHmhuw).

**DraftPlotAll.R** - R file that combines all the results and produce graphs. The figures are saved in *Figures* folder.

## 3. Applications
*Application* folder contains the necessary files to reproduce the results shown in application section. 

**Mouse** folder - code for analyses of data from nutrigenomic study of mice

- *Rank* folder includes the files for determining the total rank and joint rank based on cross-validation method and principal angles. The main file is *RankEst.R*.
- *Apply_ECCA.R* runs ECCA on the nutrimouse data and stores results at *ECCA_result.RData*.
- *PlotCluster.R* file uses ECCA joint and individual results and illustrates the low-dimensional representation of 40 mice on the ecca joint/individual basis. The graphs are stored at *Figures* folder.
- *Compare_ECCA_vs_GAS.* files provides comparison of two methods with respect to genotype and diet separation

**DeMixT_Timer_Open** folder - scripts for analyses of data from cellular heterogeneity study. 

- *RunECCA.R* call to ECCA on selected ranks
- *DeMixT_Timer_Reproduce_Downstream.pdf* rendered Rmd output of downstream analyses.
To replicate, the original sequencing data for prostate cancer patients can be obtained from [the Genomic Data Commons Data Portal](https://portal.gdc.cancer.gov/), see [Wang et al. 2018](https://doi.org/10.1016/j.isci.2018.10.028) for DeMixT application, and [Li et al., 2017](https://doi.org/10.1158/0008-5472.CAN-17-0307) for TIMER application

<!--**ItsTimer** folder includes code for analyses of Tumor heterogeneity in prostate cancer.

- *RawData* folder includes raw data and the pre-processing R file. The processed data is stored at *prostate_timer_pit3.RData*.
- *Rank* folder includes the files for determining the total rank and joint rank based on bi-cross-validation method and principal angles. The main file is *RankEst.R*.
- *RunECCA.R* runs ECCA on processed data and stores ECCA results at *ecca_result_without_main_effect.RData*.
- *ItsTimer_Analysis.R* file uses ECCA result and fit Cox Proportional Hazards model.
- *Heatmap_loading.R* file uses ECCA result and plot the loadings with heatmaps. The heatmaps are stored at *Figures* folder.
- *Survival_etc* folder includes the survival analysis. These results are not used in the paper.

-->
