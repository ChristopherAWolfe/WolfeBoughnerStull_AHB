# Overview
This github repository provides the code necessary to reproduce the results for the following article:
> Wolfe CA, Boughner JB, Stull KE. (2024). What use are ontogenetic data anyway? Challenges in multivariate modelling of primate tooth formation. *Annals of Human Biology*. In Review.

Reproduction of all results assumes three things:
1.  An updated version of Git to clone materials to your local machine,
2.  An updated version of R. Preferably, RStudio is necessary to completely restore the working environment utilized in the present results. If not, the directory can be easily cloned to a local folder and teh scripts runs individually. The instructions here will demonstrate how to complete the full workflow in RStudio.
3.  Access to the data sources used in the paper. All data are available open-access at the sources cited below. The code and instructions are written in such away to ensure any dataset could be used.

 &nbsp;

 **It will be easy for a user to click through this `README.md` and try to reproduce the results. We recommend against this naive approach**

&nbsp;

We strongly encourage one to open the code files, manipulate them, adjust to your own coding language, and adjust syntax as need be. As long as the data are from the sources described below and import according to "01_data_prep.R" and the same general steps are followed, the results will be identical. We have purposely written the code in a very general manner in order to allow for flexibility in later approaches.

## Workflow

### Prepare Environment
**Note: We assume knowledge of the command line or the terminal. One can use the terminal through RStudio.**  

&nbsp;

We can download all material in two ways:
1. For those comfortable with git, set your directory to the place you'd like the files to be stored locally, and clone the repository.

```console
cd "file/path/to/desired/repository/location"
git clone https://github.com/ChristopherAWolfe/WolfeBoughnerStull_AHB.git
```
2. For those less comfortable with git, download the repository as a zipped file by clicking the green "CODE" button on the top-right side and select 'Download ZIP'. Place the extracted files in a folder you'd like to work from on your machine.

&nbsp;

Afterwards, find the folder where you cloned / extracted the files and double click on the `.RProj` file. This will open an R Project in RStudio that will allow you to work through all analyses. Before we begin, a user must run the script `00_set_workspace.R`. This will ensure the user has the R environment and packages necessary to complete the analyses. However, there is a caveat here, see below...

```r
source("code/00_set_workspace.R")
```

&nbsp;

### Install RTools and Cmdstan
The analyses require that one compute the joint log probability density in `CmdstanR` - a lightweight interface to Stan on the command line `Cmdstan`. A user will require updated versions of `Rtools` and ` Cmdstan`.  The instructions to complete this step across several operating systems can be found [here](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) and [here](https://mc-stan.org/docs/cmdstan-guide/installation.html). 

&nbsp;

### Import Data
We have made the concious choice to exclude the data from this repository. This is not done to hide the data from others users, but instead to encourage others to seek out the initial sources. We do this to ensure those sources get the requisite credit and recognition for collecting, housing, and disseminating the data. We stress that all data are in fact available open-access. Should there be issues, please contact wolfec23@ecu.edu. The data sources are:

&nbsp;

1.  __*H. sapiens*__ - Available at Zenodo. 
> Stull KE, Corron LK. 2021. SVAD_US (1.0.0) [Data Set]. Zenodo. https://doi.org/10.5281/zenodo.5193208
2.  __*Pan*__ - Supporting Information, Table S2 and Table S4.
> Boughner JC, Dean MC, Wilgenbusch CS. 2012. Permanent tooth mineralization in bonobos (Pan paniscus) and chimpanzees (P. troglodytes). American Journal of Physical Anthropology. 149(4):560–571. https://doi.org/10.1002/ajpa.22166
3.  __*Papio*__ - Supporting Information, Table S4
> Boughner JC, Der J, Kuykendall KL. 2015. A multivariate approach to assess variation in tooth mineralization using free-lived and captive-raised chimpanzees (P. troglodytes). American Journal of Physical Anthropology. 158(3):452–462. https://doi.org/10.1002/ajpa.22800.
4. __*Hylobates*__ - Initial Publication in Review, available upon completion of the review process
> Cofran Z, Boughner JC. In Review. Permanent dental development in the white-handed gibbon (Hylobates lar carpenteri).

&nbsp;
This next code chunk assumes one has accessed the requisite data, saved it in `data` folder, and named it accordingly: *Homo* = "SVAD_US.csv", *Pan* = "pan.csv", *Papio* = "papio.csv", and *Hylobates* = "hylobates.csv". Assuming this is completed, one can run:
```r
source("code/01_data_prep.R")
```
This will import all data, clean it up for analyses, and output 4 objects in your R environment. Note, it is possible to **not** `source` the script and run it line by line.

&nbsp;

### Fit a Multivariate Probit Model
This is the most crucial step to recreating the analyses in the main text. This script takes in the data described above and fits it to a multivariate cumulative probit model in Stan. This can be done with continuous predictor information ("02a_fit_model.R") or ordinal predictor information ("02b_fit_model.R"). This step will take awhile, anywhere from a few hours to a few days depending on the sample size. These scripts are created assuming the *age* model is best (see "03_..._CV.R" scripts). Each model will output a `CmdStanMCMC ` object into the R environemnt and will save an RDS file to the "fitted_models" folder. We recommend going to this file and renaming it per your own use. We also recommend running each species seperately, ensuring each fit is saved, and moving on to the remainder of the scripts. One example for how to complete the analyses:
```r
# Set your analysis type
# analysis <- "human", analysis <- "pan", analysis <- "papio", analysis <- "hylobates"

if(analysis == "human"){

dat <- human_dat
source("code/02a_fit_model.R")

} else if (analysis == "pan"){

dat <- pan_dat
source("code/02b_fit_model.R")

} else if(analysis == "papio"){

dat <- papio_dat
source("code/02b_fit_model.R")

} else{

dat <- hylobates_dat
source("code/02b_fit_model.R")

}
```
&nbsp;

We recommend running this 4 times, ensuring the model is saved in the "fitted_models" folder, renaming the RDS, and moving on to the next set of analyses. 

&nbsp;

### Secondary Analyses
After the model is fit, a user can complete the remainder of the analyses in any order. However, the scripts are designed for files to be names a certain way. As such, again we recommend a user at least become familiar with the code before running `source(...)` and receiving an error because there is a naming mismtach. See general instructions below:

&nbsp;

#### Model Selection
This series of scripts performs model selection using 5-fold cross-validation and the expected log pointwise predictive density (EPLD). Note, this whole process may take a few days to a week because we are refitting each model 5 times. For example, in *Pan* we fit a model with and without age - this is 10 total models as we perform cross-validation seperately and than compare the ELPD. There are 4 total scripts. To do this for a continuous predictor we need "03a_continuous_age_CV.R", "03c_intercept_age_CV.R", and "04_model_selection.R". To do this for an ordinal predictor we need "03b_ordinal_age_CV.R", "03c_intercept_age_CV.R", and "04_model_selection.R". In each case, the "..._CV.R" files import data, fit a model, test the model, and export log-likelihood values. The "04_model_selection.R" takes these values and performs the comparison. One example of how to complete this could be:

```r
# Set your analysis type
# analysis <- "human", analysis <- "pan", analysis <- "papio", analysis <- "hylobates"

if(analysis == "human"){

dat <- human_dat
source("code/03a_continuous_age_CV.R")
source("code/03c_intercept_age_CV.R")
source("code/04_model_selection.R")

} else if (analysis == "pan"){

dat <- pan_dat
source("code/03b_ordinal_age_CV.R")
source("code/03c_intercept_age_CV.R")
source("code/04_model_selection.R")

} else if(analysis == "papio"){

dat <- papio_dat
source("code/03b_ordinal_age_CV.R")
source("code/03c_intercept_age_CV.R")
source("code/04_model_selection.R")

} else{

dat <- hylobates_dat
source("code/03b_ordinal_age_CV.R")
source("code/03c_intercept_age_CV.R")
source("code/04_model_selection.R")

}

```

&nbsp;

#### Posterior Correlations
Here we visualize the posterior correlation matrices from the main text. Note, the script assumes the initial models are saved as "fit.rds". This can be adjusted as needed. 

```r
# Assuming the model from "02_model_fit.R is saved as 'fit$save_object("fitted_models/fit.rds")'
source("05_posterior_corrs.R")
```

&nbsp;

#### Eigendecomposition
Here we perform eigendecomposition of the posterior correlation matrices in the main text. Note, the script assumes the initial models are saved as "fit.rds". This can be adjusted as needed. 

```r
# Assuming the model from "02_model_fit.R is saved as 'fit$save_object("fitted_models/fit.rds")'
source("06_eigendecomposition.R")
```

&nbsp;

#### Matrix Comparisons
This final script compares the posterior correlations between **ALL** models. To keep it general, we assume the *Homo* model has been saved as `fit$save_object("fitted_models/fit1.rds"`, the *Pan* model has been saved as `fit$save_object("fitted_models/fit2.rds"`, the *Papio* model has been saved as `fit$save_object("fitted_models/fit3.rds"`, and the *Hylobates* model has been saved as `fit$save_object("fitted_models/fit4.rds"`. If this is the case, the below script will reproduce exactly the figure from the main text. Here, we compute the Frobenius norm between all 4000 posterior correlation matrices and then find the posterior means of this difference. We do this is parallel - to avoid computational bottlenecks, we recommend a user adjust the # of cores in this script. If not run in parallel, this script will take days to complete. 

```r
source("07_matrix_comparisons.R")
```

#### Other Analyses
The main text does present boxplots for comparative purposes and posterior predictive checking on a single variable using the `bayesplot` package. We have left this code out of this repository because we felt it was not germane to the general reproduction of results. A user can easily perform their own posterior checks using a combination of `bayesplot` and `generated_quantities` in Stan. The boxplots can be recreated in both base R and `ggplot2`. 

## Additional Information
For additional information regarding the main text or the code published herein, please contact Dr. Christopher Wolfe at wolfec23@ecu.edu. To cite this repository please use:
>
