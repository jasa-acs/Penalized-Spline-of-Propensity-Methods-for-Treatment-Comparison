# Penalized Spline of Propensity Methods for Treatment Comparison

# Author Contributions Checklist Form

## Data

### Abstract

We used data from the Multicenter AIDS Cohort study (MACS) (Kaslow et al, 1987). The
MACS was started in 1984, and a total of 4,954 gay and bisexual men were enrolled in the study
and followed up semi-annually. At each visit, data from physical examination, questionnaires
about medical and behavioral history, and blood test results were collected. For this analysis, we
restricted the analysis to data from visit 7 and 21 and focused on blood test results. The CD4
count was the intermediate and final outcome of interest. The CD4 count is a time-dependent
confounder because it is both an intermediate outcome of past antiretroviral treatments and a
confounder of future antiretroviral treatments. The goal of the analysis to estimate the short term
effect of antiretroviral treatments, i.e. every moving three-visit window from visit 7. See the
study website for more information https://statepi.jhsph.edu/macs/pdt.html

### Availability

We requested the dataset from Johns Hopkins University Bloomberg School of Public Health
Department of Epidemiology. The version we used is P25, released October 2016. Some request
forms were completed in order to obtain it. The entire data set can be obtained by contacting the
current research program manager of the Multicenter AIDS Cohort Study, Jeremy LaMaster at
jlamaster@jhu.edu. However, the original data that we analyzed for the paper are included in the
supplementary materials.

### Description

The version we used is P25, released October 2016. The data files are in .dat format but we
converted them to csv files before analysis. The datasets that we analyzed for the paper are
included in the supplementary materials. Inside the subfolder dataset in the application folder,
there is a zip file named rawDatasets.7z, which contains the original datasets (macsid.csv,
section2.csv, outcome.csv, hivstats.csv, section4.csv, lab_rslt.csv). We used the questions in
these files to determine HIV status and blood counts and demographics for each subject. The
dataProcessing.R script produces the final dataset we used-data.csv.

### Optional Information

This is a public data set, so there are no unique identifiers.

## Code

### Abstract

All the codes are inside three folders: Application, oneTimePointSimulation, and
twoTimePointSimulation. Inside each folder, there is a ReproduceTablesandFigures.pdf detailing
the codes and how to reproduce the tables and figures.

### Description

The supplementary materials we submitted will be deposited on GitHub and made publicly
available.

### Optional Information 

All the simulations and applications were run on a high performance computing cluster.
R packages nlme, specifically lme function (most recent version 3.1-131.1), and R package rootS
olve (version1.7) are needed to run the codes. All the work was done in R version 3.4.4

## Instructions for use

### Reproducibility 

All the codes to reproduce all the tables and figures from paper are inside the subfolder
FiguresandPlots in each of the main folders: Application, oneTimePointSimulation, and
twoTimePointSimulation. The code file name includes the table and figure numbers that it
reproduces. Also the ReproduceTablesandFigures.pdf inside each of the main folders describes
how to reproduce the tables and figures.

### Notes

If you have questions regarding implementation of the method and codes, contact the authors.
