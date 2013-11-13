# Prerequisite
You need to install the following packages :

- R
- knitr (from CRAN : install.packages('knitr', dependencies = TRUE) )
- Any latex package to generate PDF

# Generate reports
In order to generate chimCT-report as PDF run :

  knit chimCT-report.Rnw

If you only want to generate a the *tex* file use (--no-convert option) :

  knit -n chimCT-report.Rnw

# Others
Figures will be saved in the /figure folder. 
