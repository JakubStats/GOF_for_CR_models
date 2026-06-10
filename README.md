-------------------------------------------------------------------------------
[![DOI](https://zenodo.org/badge/1068249756.svg)](https://doi.org/10.5281/zenodo.20621527)

README: Reproducible R-code for "Residual diagnostics for assessing closed population 
capture-recapture models.” 
-------------------------------------------------------------------------------

Jakub Stoklosa, Bernice Laitly, and Pierre Lafaye de Micheaux. 

Accepted in Methods in Ecology and Evolution.

All real examples and simulation studies are presented in the accompanying supplementary 
materials (Web Appendices C and D).

-------------------------------------------------------------------------------
FILES INCLUDED
-------------------------------------------------------------------------------

Main analysis file:
  - CR_GOF_Examples.R

Simulation files:
  - CR_GOF_Simulation_Scenario1.R
  - CR_GOF_Simulation_Scenario2.R
  - CR_GOF_Simulation_Scenario3.R
  - CR_GOF_Simulation_Scenario4.R
  - CR_GOF_Simulation_Scenario5.R
  - CR_GOF_Simulation_Scenario6.R
  - CR_GOF_Simulation_Scenario7.R

Custom functions:
  - CR_GOF_Functions.R

Data files:
  - poss2003.RData (Obtained from Stoklosa et al. (2012)). It contains capture histories as well as
    body weight measurements for each capture possum. The number of unique possums captured is D = 43.
    The number of capture occasions is 5. 
  - TaxiM.inp (obtained from the publicly avaiable software MARK: White, G. C. & Burnham, K. P. (1999). 
    Program MARK: survival estimation from populations of marked animals. Bird study, 
    46, S120–S139). It only contains capture histories across 10 capture occasions. The number of unique
    taxi cabs identified (or captured) is D = 283.

Additional C++ files:
  - TMB_eps_NEW2.cpp
  - TMB_Mh_RE.cpp

-------------------------------------------------------------------------------
REQUIRED R PACKAGES
-------------------------------------------------------------------------------

The following R-packages are required:

  library(DHARMa)
  library(mgcv)
  library(tidyverse)
  library(secr)
  library(VGAM)
  library(TMB)
  library(RMark)
  library(marked)
  library(glmmTMB)

If these packages are not installed, they can be installed using:

  install.packages(c(
    "DHARMa", "mgcv", "tidyverse", "secr",
    "VGAM", "TMB", "RMark", "marked", "glmmTMB"
  ))

-------------------------------------------------------------------------------
TMB / C++ REQUIREMENTS
-------------------------------------------------------------------------------

Example 3 and Simulation Scenario 6 additionally require the Template Model
Builder (TMB) package and the accompanying C++ files:

  - TMB_eps_NEW2.cpp
  - TMB_Mh_RE.cpp

These files must be located in the working directory when running the
corresponding analyses.

Users must also have a working C++ compiler configured for R to compile the 
TMB models. See:

  https://cran.r-project.org/package=TMB

for installation and setup instructions.

-------------------------------------------------------------------------------
WORKING DIRECTORY SETUP
-------------------------------------------------------------------------------

To ensure the code runs correctly, place all R scripts, function files,
data files, and C++ files in the same working directory before running
the analyses.

The scripts assume that:

  - CR_GOF_Functions.R
  - poss2003.RData
  - TaxiM.inp
  - TMB_eps_NEW2.cpp
  - TMB_Mh_RE.cpp

are located in the current working directory.

-------------------------------------------------------------------------------
RUNNING THE REAL DATA EXAMPLES
-------------------------------------------------------------------------------

The file:

  CR_GOF_Examples.R

contains the code used to reproduce the real-data analyses presented in the
manuscript and Web Appendix C.

-------------------------------------------------------------------------------
RUNNING THE SIMULATION STUDIES
-------------------------------------------------------------------------------

Simulation scenarios are provided as separate R scripts:

  - CR_GOF_Simulation_Scenario1.R
  ...
  - CR_GOF_Simulation_Scenario7.R

Each file reproduces the corresponding simulation scenario described in
Web Appendix D.

Random seeds (set.seed()) are included in the scripts to facilitate
reproducibility of the reported results.

-------------------------------------------------------------------------------
NOTES
-------------------------------------------------------------------------------

The code was developed and tested in R (version 4.6.0).

The custom functions contained in:

  CR_GOF_Functions.R

which prepares capture-recapture data under the partial likelihood framework. 
For additional details, see Web Appendix B.

-------------------------------------------------------------------------------
CONTACT
-------------------------------------------------------------------------------

For questions regarding the code or methodology, please contact:

  Dr. Jakub Stoklosa (j.stoklosa@unsw.edu.au)
