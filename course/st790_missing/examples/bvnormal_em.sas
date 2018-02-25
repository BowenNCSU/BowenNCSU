/*********************************************************************

    Use PROC MI to carry out the EM algorithm on bivariate normal
    data given in EXAMPLE 2 of Chapter 3 of the notes

    We read in a simulated data set created in the R program
    bvnormal_em.R with the "NAs" replaced by "." (the SAS missing value
    indicator)

***********************************************************************/

options ls=80 nodate; run;

/*********************************************************************

  Read in the data set

***********************************************************************/

data bivariate;  infile "bvnormal.dat";
    input y1 y2;
run;

/**********************************************************************

  Call proc mi to carry out the EM algorithm.  Here, we use the
  complete case usual sample mean and sample covariance matrix as
  the starting values (initial=cc), set the maximum number of
  iterations to 100, and specify the convergence criterion as 1e-4.
  The nimpute=0 options suppresses multiple imputation, so that the
  procedure does only the EM.  The seed for the random number generator
  is for multiple imputation and thus is superfluous here. The output
  data set "outem" contains the final estimates.  See the documentation
  for more options.
    
***********************************************************************/

proc mi data=bivariate seed=1518971 simple nimpute=0;
    em initial=cc maxiter=100 converge=1e-4 itprint outem=outem;
    var y1 y2;
run;




