/******************************************************************

  CHAPTER 8, EXAMPLE 6, Six Cities Study
  
  Fit a logistic regression model to the "wheezing" data.
  These are binary data, thus, we use the Bernoulli (bin) 
  mean/variance assumptions.  The model is fitted with different
  working correlation matrices using PROC GENMOD, PROC GEE, 
  PROC GLIMMIX
    
******************************************************************/

options ls=80 ps=59 nodate; run;

/******************************************************************

  The data look like (first 4 records):

    1 portage   9 0 1  10 0 1  11 0 1  12 0 0
    2 kingston  9 1 1  10 2 1  11 2 0  12 2 0
    3 kingston  9 0 1  10 0 0  11 1 0  12 1 0
    4 portage   9 0 0  10 0 1  11 0 1  12 1 0

          .
          .
          .

  column 1      child
  column 2      city
  columns 3-5   age=9, smoking indiciator, wheezing response
  columns 6-8   age=10, smoking indiciator, wheezing response
  columns 9-11  age=11, smoking indiciator, wheezing response
  columns 12-14 age=12, smoking indiciator, wheezing response

  Some of the children have missing values for smoking and wheezing,
  as shown in Chapter 1.  There are 32 children all together.  See the
  output for the full data printed out one observation per line.

  We read in the data using the "@@" symbol so that SAS will continue
  to read for data on the same line and the OUTPUT statement to 
  write each block of three observations for each age in as a separate
  data record.  The resulting data set is one with a separate line for
  each observation.  City is a character variable.

******************************************************************/

data wheeze; infile 'wheeze.dat';
  input child city$ @@;
  do i=1 to 4;
    input age smoke wheeze @@;
    output;
  end;
run;

proc print data=wheeze(obs=10); run;

/*****************************************************************

  Fit the logistic regression model using PROC GENMOD and PROC GEE
  three different working correlation matrix assumptions:

  -  unstructured
  -  compound symmetry
  -  AR(1)

  We fit a model with linear predictor allowing effects of 
  city and maternal smoking status subject to the caveats noted
  in Sections 8.2, 8.6, and 8.8. This is a marginal model, so
  given concerns over the observational nature of this study and the
  likelihood that the covariate smoking status is endogenous, it
  might be saftest to report the fit with the independence working assumption.
  
  The DIST=BIN option in the MODEL statement specifies that the
  Bernoulli mean-variance relationship be assumed.  The LINK=LOGIT
  option asks for the logistic mean model.

  The REPEATED statement specifies the working correlation
  structure to be assumed.    The CORRW option in the REPEATED 
  statement prints out the estimated working correlation matrix
  under the assumption given in the TYPE= option.  The COVB
  option prints out the estimated covariance matrix of the estimate
  of beta -- both the usual estimate and the robust version
  are printed.  The MODELSE option specifies that the standard
  error estimates printed for the elements of betahat are based
  on the usual theory.  By default, the ones based on the "robust"
  version of the sampling covariance matrix are printed, too.  

  The dispersion parameter phi is held fixed at 1 by default.

  The missing values are coded in the usual SAS way by periods (.).
  We delete these from the full data set, so that the data set input
  to PROC GENMOD contains only the observed data.  As noted in Section 
  8.7, the missingness mechanism must be MCAR for the analysis to be
  valid, which is very likely not true.  Thus, this analysis is for
  illustrative purposes only to demonstrate how to deal with unbalanced
  data.  We must use the WITHIN option of the REPEATED
  statement to give SAS the time variable AGE as a
  classification variable so that it can figure out where the missing
  values are and use this information in estimating the correlation 
  matrix.

  PROC GENMOD models by default the probability that Y=0 rather
  than Y=1.  To make PROC GENMOD model probability that 
  Y=1, as is standard, one must include the DESCENDING option in
  the PROC GENMOD statement.  An explicit statement about 
  what is being modeled will appear in the log file and output.

 ******************************************************************/

data wheeze; set wheeze; 
  if wheeze=. then delete;
  time=age;
run;
 
title "PROC GENMOD, UNSTRUCTURED CORRELATION";
proc genmod data=wheeze descending;
  class child city smoke time;
  model wheeze = city smoke / dist=bin link=logit;
  repeated subject=child / type=un corrw covb modelse within=time;
run;

title "PROC GENMOD, INDEPENDENCE";
proc genmod data=wheeze descending;
  class child city smoke time;
  model wheeze = city smoke / dist=bin link=logit;
  repeated subject=child / type=ind covb modelse within=time;
run;

title "PROC GENMOD, COMPOUND SYMMETRY (EXCHANGEABLE) CORRELATION";
proc genmod data=wheeze descending;
  class child city smoke time;
  model wheeze = city smoke / dist=bin link=logit;
  repeated subject=child / type=cs corrw covb modelse within=time;
run;

title "PROC GENMOD, AR(1) CORRELATION";
proc genmod data=wheeze descending;
  class child city smoke time;
  model wheeze = city smoke / dist=bin link=logit;
  repeated subject=child / type=ar(1) corrw covb modelse within=time;
run;

*  Using PROC GEE;

title "PROC GEE, INDEPENDENCE";
proc gee data=wheeze descending;
  class child city smoke time;
  model wheeze = city smoke / dist=bin link=logit;
  repeated subject=child / type=ind covb modelse within=time;
run;

title "PROC GEE, COMPOUND SYMMETRY (EXCHANGEABLE) CORRELATION";
proc gee data=wheeze descending;
  class child city smoke time;
  model wheeze = city smoke / dist=bin link=logit;
  repeated subject=child / type=cs corrw covb modelse within=time;
run;

/*****************************************************************

  We demonstrate the use of PROC GLIMMIX to instead estimate
  the correlation parameters using the quadratic estimating 
  equation.  This procedure is really meant for fitting    
  SS generalized linear mixed effects models as in Chapter 9 of the 
  course.  However, just as PROC MIXED can be used to fit PA linear 
  models, GLIMMIX can be used to fit PA models of the "generalized
  linear models type." Some of the syntax is indeed similar to that
  for PROC MIXED as a result.

  The DIST, LINK, and OFFSET options in the model statement are the
  same as for GENMOD and GEE.  As with MIXED, one must use the SOLUTION
  option in the MODEL statement to get the parameter estimates. The 
  EMPIRICAL option in the PROC statement requests robust empiriacal
  standard errors.  The METHOD=MMPL option in the PROC statement requests
  that the quadratic estimating equation for estimating the correlation
  parameter and scale parameter be used.  Another option, RMPL, uses a
  REML version of this.

  Specification of the overall covariance matrix, and particularly the 
  working correlation structure, is a bit different from the above
  procedures.  Because this procedure really isn't directly designed
  to fit PA models, the syntax is a bit arcane.  Instead of having a 
  REPEATED statement that can be tricked into fitting an overall
  covariance matrix as we did with PROC MIXED, GLIMMIX uses the "residual"
  option to do this in a RANDOM statement.  The way it is used below
  requests that what SAS calls a "R-side" covariance matrix (what we 
  would call the overall covariance matrix in this context) be 
  constructed.  You should consult the documentation for details. 

******************************************************************/

title "PROC GLIMMIX, COMPOUND SYMMETRY (EXCHANGEABLE) CORRELATION";
proc glimmix data=wheeze method=mmpl empirical;
  class child city smoke time;
  model wheeze = city smoke / dist=bin link=logit covb solution;
  random time / subject=child type=cs vcorr residual;
run;

