/******************************************************************

  CHAPTER 9, Argatroban pharmacokinetic study

  Fit a subject-specific nonlinear mixed effects model to the
  argatroban PK data in Section 9.5 of Davidian and Giltinan (1995)
  by fitting a "linear mixed effects model" using transformed
  individual-specific estimates as "data."  These "data" are 
  read in from the first stage pooled GLS algorithm in the 
  program arg_mixsetup.R, which creates the needed file "argmixsig.dat"

******************************************************************/

options ps=59 ls=80 nodate;

/*  read in the "data"  */
    
data arg; 
infile "argmixsig.dat";  
  input id y x1 x2 z1 z2 @@;
run;

/*  now run PROC MIXED on the "data".  The PARMS statement
    specifies initial values for the covariance parameters, 
    which here are the distinct elements of D (in the order
    PROC MIXED will print them in the Covariance Parameters
    table of the output) followed by the within-individual
    variance sigma^2 -- here, we are fitting the default within-
    individual covariance model sigma^2 I.  The HOLD option
    allows the user to hold some or all of these values fixed at the
    given values - here, we hold the 4th, corresponding to 
    sigma^2, = 1.  The values for the elements of D are thus 
    starting values and were chosen based on other fits of
    the model using other methods.
  
*/
    
proc mixed data=arg method=ml;
  class id;
  model y = x1 x2 / noint solution chisq;
  random z1 z2 / subject=id type=un  g gcorr gc;
  parms (0.14) (0.006) (0.006) (1) / hold=4;
run;

