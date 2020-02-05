<!-- README.md is generated from README.Rmd. Please edit that file -->
MarginalModelCifti
==================

The goal of MarginalModelCifti is to perform marginal models on CIFTI processed datasets. The package contains a single main function that runs multiple subfunctions. Advanced users can use the subfunctions to construct their own analytic pipeline. However, this is not recommend for beginning users.

Dependencies
------------

1) MarginalModelCifti requires R version 3.5.1 or better. 
2) In addition, two DCAN github packages are also required:

	a)https://github.com/DCAN-Labs/SurfConnectivity
	b)https://github.com/DCAN-Labs/CommunityChiSquaredAnalysis
	These packages can be cloned to any folder. 
3) Using these packages requires the matlab runtime compiler engine, version 91:

https://www.mathworks.com/products/compiler/matlab-runtime.html

Installation
------------

Follow these steps to install the MarginalModelCifti package on your system, ignore the backtics if you copy and paste the commands:

1) make a directory for the MarginalModelCifti package `mkdir ~/MarginalModelCifti`
2) enter the directory `cd ~/MarginalModelCifti`
3) clone the MarginalModelCifti repository `git clone https://gitlab.com/Fair_lab/marginalmodelcifti.git ./`
4) return to your initial home directory `cd ..`
5) Type `R`
6) After a prompt appears, make sure devtools is installed by typing `install.packages("devtools")`
7) Load devtools `library(devtools)`
8) install the MarginalModelCifti package `install("MarginalModelCifti/")`
9) install RFast by typing `install.packages("RFast")`

Example
-------

To use the MarginalModelCifti package, please refer to the example scripts found in the examples folder, such as `examples/MarginalModelCifti_PCONN_analysis.R` the example Rscripts contain the documented parameters that control the analysis.

In order to setup your analysis you will need to create the following files:

1) A concfile. A single column text file where each line contains the path to an imaging metric file, such as a CIFTI connectivity matrix (pconn/dconn), CIFTI scalar (pscalar/dscalar), statistical volume(.nii), or surface metric (func.gii) 
2) A structural file. A file representing the underlying structural surface or binary volume mask -- used for cluster detection.
3) An external data file. A csv file where each row corresponds to each imaging metric file, and each column contains an external non-imaging measure. The first row represents the variable names. Please note that the program will only use the columns specified by the user, so the csv file can contain additional columns without error.
4) A wave file. A csv file specifying the nested structure for any nested variables. Each row corresponds to an imaging metric file and each column contains a grouping variable. As many grouping variables can be specified as needed. All variables contained in the file will be used.


The main wrapper of the package can be called as follows:

``` r
map_objects <- ConstructMarginalModel(external_df=external_df,concfile=concfile,structtype=structtype,structfile=structfile,matlab_path=matlab_path,surf_command=surf_command,wave=wave,notation=notation,zcor=zcor,
corstr=corstr,family_dist=family_dist,dist_type=dist_type,z_thresh=z_thresh,nboot=nboot,p_thresh=p_thresh,sigtype=sigtype,id_subjects=id_subjects,output_directory=output_directory,
ncores=ncores,fastSwE=fastSwE,adjustment=adjustment,norm_external_data=norm_external_data,norm_internal_data=norm_internal_data,marginal_outputs=marginal_outputs,
marginal_matrix=marginal_matrix,enrichment_path=enrichment_path,modules=modules,wb_command=wb_command)

```
