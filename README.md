-This repository contaiancs codes and data used in Uddin et al. 2023
(arXiv: 2308.01875)

 titled as:  Carnegie Supernova Project-I & -II: Measurements of H0 using Cepheid, TRGB, and SBF Distance Calibration to Type Ia Supernovae

-The primary code /scripts/H0CSP.py can be executed to sample parameteres using
the MCMC sampler EMCEE. The input argument is a file name that
contains both the distant SNe Ia and calibrators.

-Files are located at data/working/. We have blended calibrators and
distanct SNe Ia into individual files for easy use. For example:
data/working/B_ceph_update3.csv contains distant SNe Ia and Cepheid calirators
in B band maxmimum magnitudes. Column name 'dist' or 'caltype' are
used to separate distant SNe Ia from calibrattors.

-Resulting corner plots with marginalized distributions are located at
plots/ and parameters values at results/.

-Other codes used to generate figure and analysis are located at
scripts/misc/.

-Please contact at saushuvo@gmail.com for any issue.


 

