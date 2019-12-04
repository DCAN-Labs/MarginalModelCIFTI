# title: "MarginalModelReliability_example"
# author: "Eric Feczko"
# date: "10/10/2019"

# To install MarginalModelCifti, one should do the following:
# 1) make a directory for the MarginalModelCifti package `mkdir ~/MarginalModelCifti`
# 2) enter the directory `cd ~/MarginalModelCifti`
# 3) clone the MarginalModelCifti repository `git clone https://gitlab.com/Fair_lab/marginalmodelcifti.git ./`
# 4) return to your initial home directory `cd ..`
# 5) Type `R`
# 6) After a prompt appears, make sure devtools is installed by typing `install.packages("devtools")`
# 7) Load devtools `library(devtools)`
# 8) install the MarginalModelCifti package `install("MarginalModelCifti/")`
# 
# *NOTE: You may also want to clone the SurfConnectivity package, in case you do not have access to it.*
# a) make a directory for SurfConnectivity `mkdir ~/SurfConnectivity`
# b) go into SurfConnectivity folder `cd ~/SurfConnectivity`
# c) clone the SurfConnectivity repository here `git clone https://gitlab.com/Fair_lab/surfconnectivity.git ./`
# 
# *NOTE: You may also want to clone the CommunityChisquared package, in case you do not have access to it.*
# i) make a directory for CommunityChisquared `mkdir ~/CommunityChiSquaredAnalysis`
# ii) go into CommunityChisquared folder `cd ~/CommunityChiSquaredAnalysis`
# iii) clone the CommunityChiSquared repostory here `git clone https://github.com/DCAN-Labs/CommunityChiSquaredAnalysis.git ./`

# To run this package on a cluster, e.g. slurm, one can use Rscript to call this script:
# srun -t 36:00:00 -c 12 --mem-per-cpu=8GB -e /path/to/log.err -o /path/to/log.out Rscript /path/to/parameter_script.R

### Call the MarginalModelCifti library -- if this errors you will need to install it using devtools

library(MarginalModelCifti)


### Set your project folder, which is where you plan to run the analysis, then go to the folder

projectsfolder="/mnt/rose/shared/projects/ABCD/avg_pconn_maker/"
setwd(projectsfolder)
getwd()

### Now declare the needed variables to run a marginal model

### Set the below variable to your subset file. The subset file is a two-column file with the name of the subject contained within. This name should be unique and contained within the conc files as well.

subset_file="/mnt/rose/shared/projects/ABCD/avg_pconn_maker/subset_analyses/raw/group1_needed_subsets/group1_subset_10_with_50_subjects.csv"

### Set the below variable to your external (i.e. non-imaging) dataset. The dataset should be a csv with headers representing the variables.

external_df <- "/mnt/rose/shared/projects/ABCD/avg_pconn_maker/cordova_analysis_margmod_pcs/gp1_10min_pconn_reordered_nonan.csv"

### Set the below variable to a single column textfile, where each row contains a path to each participant's metric file.
### The metric file should be ordered in the same order as the external_df file

concfile <- "/mnt/rose/shared/projects/ABCD/avg_pconn_maker/cordova_analysis_margmod_pcs/group1_10min_nonan.conc"

### Set the below variable structtype to the type of brain structure (i.e. "surface", "volume", or "pconn") for the metric file.
### This is needed for the cluster detection to work properly.


structtype="pconn"

### If the `structtype` is set to "surface", set the below variable structfile to the corresponding surface file.
### This is needed for surface-based cluster detection. The value can be set to NULL for volumes.


structfile=NULL


### If the structtype is set to "surface" or "pconn", set the below variable matlab_path to the matlab2016b compiler.

matlab_path="/mnt/max/shared/code/external/utilities/Matlab2016bRuntime/v91"


### If the structtype is set to "surface", set the below variable to the SurfConnectivity script

surf_command="/mnt/max/shared/projects/FAIR_users/Feczko/code_in_dev/SurfConnectivity/"


### Specify the model you want to run in the below variable notation.
### Use the function `formula` to make a formal notation within R. 
### The predicted variable will always be y representing the values in the metric file. 
### The predictor variables should use the column names within the `external_df` csv header.


notation = formula(y~pc1_new)

### Set the below variable `corstr` to the correlation structure of the cases. Usually this should just be "independence".


corstr="independence"

### Set the below variable `family_dist` to the appropriate distribution of your data, "gaussian" is the default


family_dist="gaussian"

### Set the below variable `dist_type` to the distribution used for wild bootstrapping.
### Acceptable values are "radenbacher", "webb4", "webb6", and "mammen".


dist_type="radenbacher"

### set the below variable `z_thresh` to the z statistic threshold used for determining observed and permuted cluster sizes


z_thresh = 2.3

### Set the below variable `nboot` to the number of wild bootstraps to perform. 
### The precision of the p value is equal to 1/`nboot`.
### For example, if 1000 bootstraps are selected, the smallest p value can be 0.001. 
### WARNING, this part can be slow.


nboot=4

### Set the below variable `p_thresh` to the p value threshold for assessing significant clusters.
### Currently this has no functionality


p_thresh=0.05

### Set the below variable `sigtype` to determine how to perform multiple comparison correction. 
### Acceptable values are "point" (FWE for voxels), "enrichment", and "cluster".


sigtype="enrichment"

### Set the below variable `id_subjects` to the column header containing the subject id in the `external_df` file.


id_subjects="subjectkey"

### Set the below variable `output_directory` to where you want to save your outputs


output_directory="/mnt/rose/shared/projects/ABCD/avg_pconn_maker/subset_analyses/analyses/scratch/group1_observed//perm5"

### Set the below variable `zcor` to a custom covariance matrix to denote participant similarity (e.g. a kinship or site matrix)

zcor=NULL

### The below variable `fastSwE` enables the fast sandwich estimator. 
### Generally, it is recommended to set this to true to reduce overhead computations.


fastSwE=TRUE

### The below variable handles sample size adjustments when dealing with small sample sizes. 
### Acceptable transforms are currently "HC2" and "HC3". "HC2" and "HC3". 
### Correction applied as the homogenous version, per Bryan Guillaume and Tom Nichols's work.


adjustment=NULL

### Set the below variable `wave` to a csv file that denotes how subjects should be grouped and nested

wave <- "/mnt/rose/shared/projects/ABCD/avg_pconn_maker/cordova_analysis_margmod_pcs/gp1_marg_nested_reordered_nonan.csv"

### The below variable will normalize the external_df data per variable if set to true


norm_external_data=TRUE

### The below variable will normalize the imaging data per datapoint if set to true


norm_internal_data=TRUE

### The below variable will output marginal values as statistical maps if set to true


marginal_outputs = FALSE

### The below variable represents a numeric matrix for drawing the map, only useful for marginal outputs


marginal_matrix = NULL

### The below variable is a string representing the path to the enrichment repository and compiled code


enrichment_path = "/mnt/max/shared/projects/FAIR_users/Feczko/code_in_dev/CommunityChisquaredAnalysis/"

### The below variable is a string that represents the path to a csv containing modules


modules = "/mnt/max/shared/projects/FAIR_users/Feczko/code_in_dev/CommunityChisquaredAnalysis/gordon_modules.csv"

### The below variable is a string that represents the path to the workbench command


wb_command = "/usr/local/bin/wb_command"

### The below variable represents a path to a subset file for sub-selecting cases for analysis. If set to NULL, this will be ignored.

subsetfile="/mnt/rose/shared/projects/ABCD/avg_pconn_maker/subset_analyses/raw/group1_needed_subsets/group1_subset_10_with_50_subjects.csv"

### The below variable represents a path to the output permutation file. 

output_permfile = '/mnt/rose/shared/projects/ABCD/avg_pconn_maker/subset_analyses/analyses/scratch/perms_group1/permfile.txt1'

### With all the variables determined, you can now run the MarginalModel package using the `ConstructMarginalModel` command 

PermuteMarginalModel(external_df=external_df,
                     concfile=concfile,
                     structtype=structtype,
                     structfile=structfile,
                     matlab_path=matlab_path,
                     surf_command=surf_command,
                     wave=wave,
                     notation=notation,
                     zcor=zcor,
                     corstr=corstr,
                     family_dist=family_dist,
                     dist_type=dist_type,
                     z_thresh=z_thresh,
                     nboot=nboot,
                     p_thresh=p_thresh,
                     sigtype=sigtype,
                     id_subjects=id_subjects,
                     output_directory=output_directory,
                     fastSwE=fastSwE,
                     adjustment=adjustment,
                     norm_external_data=norm_external_data,
                     norm_internal_data=norm_internal_data,
                     marginal_outputs=marginal_outputs,
                     marginal_matrix=marginal_matrix,
                     enrichment_path=enrichment_path,
                     modules=modules,
                     wb_command=wb_command,
                     subsetfile=subsetfile,
                     output_permfile = output_permfile)

