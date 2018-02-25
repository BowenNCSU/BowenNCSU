/*********************************************************************

    Use PROC MI to carry out the EM algorithm on bivariate normal
    data given in EXAMPLE 2 of Chapter 3 of the notes

    Also use PROC MIXED to maximize the observed data likelihood
    directly

***********************************************************************/

options ls=80 nodate; run;

/*********************************************************************

  Read in the data set

***********************************************************************/

data bivariate;  infile "bvnormal.dat";
    input y1 y2;
run;

/*********************************************************************

  Read in the data set

***********************************************************************/

data bivariate; set bivariate; id = _n_; run;

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

/**********************************************************************

  Now reconfigure the data set from the "wide" format of one record
  per individual to the "long" format of one record per observation
  as required by proc mixed
    
***********************************************************************/

data bivariate2; set bivariate;
  array wt(2) y1 y2;
  do time = 1 to 2;
     y = wt(time);
     output;
  end;
  drop y1 y2;
run;

/**********************************************************************

  Call proc mixed to maximize the observed data likelihood.  Be
  sure to specify "method=ML" in the proc mixed statement; otherwise,
  it will use the default of restricted maximum likelihood.  
    
***********************************************************************/

proc mixed data=bivariate2 method=ml;  
  class time id;
  model y = time / noint solution;
  repeated / subject=id type=un r=1 rcorr=1;
run;  



