 /*******************************************************************

  CHAPTER 5, EXAMPLE 1, Dental Study

  Population-averaged model
     
  We use the REPEATED statement of PROC MIXED with the
  TYPE= options to fit the model assuming several different 
  covariance structures.

*******************************************************************/

options ls=80 ps=59 nodate; run;

/******************************************************************

  Read in the data set
    
*******************************************************************/

data dent1; infile 'dental.dat';
  input obsno child age distance gender;
  ag = age*gender;
run;

/*******************************************************************

  Fit the basic model with separate intercepts and slopes for 
  both genders with several different assumptions on the overall
  covariance structure.  Overall covariance structure is specified
  in the REPEATED statement.  We take this to be the same for both genders
  as well as different for each gender, which is specified using the
  GROUP = option in the REPEATED statement.  We also include fitting
  of a heterogeneous compound symmetric structure, which allows
  different variances at each time point.  See the SAS documentation
  for the REPEATED statement for a list of all the possible covariance
  model specifications, of which there are many.

  We use ML to fit all these models, but in the first call to MIXED
  we use REML (default) and show the option for getting the "robust" sandwich 
  covariance standard errors (EMPIRICAL option)
    
  We can compare the various fits of the same mean model using AIC and BIC.
  When ML is used, PROC MIXED defines BIC using what we have called N, and
  AND gls() uses what we have called m.  When REML is used, the programs
  use different conventions on both nunber of observations and number 
  of parameters.  Thus, these will differ across implementations
  but can be compared within implementation.
    
  The R and RCORR options in the REPEATED statement have an option to
  print out the estimated covariance and correlation matrices for
  specific -- we print out these for the first girl (1) and first boy (12).
  
  *******************************************************************/

title "COMMON UNSTRUCTURED USING REML (default) AND ROBUST SES";
proc mixed empirical data=dent1;
  class gender child;
  model distance = gender gender*age / noint solution ;
  repeated / type = un subject = child r rcorr;
run;

title "(a) COMMON UNSTRUCTURED";
proc mixed data=dent1 method=ml;
  class gender child;
  model distance = gender gender*age / noint solution ;
  repeated / type = un subject = child r rcorr;
run;

title "(b) SEPARATE UNSTRUCTURED BY GENDER";
proc mixed data=dent1 method=ml;
  class gender child;
  model distance = gender gender*age / noint solution ;
  repeated / type = un subject = child r=1,12 rcorr=1,12 group=gender;
run;

title "(c) COMMON COMPOUND SYMMETRY STRUCTURE";
proc mixed data=dent1 method=ml;
  class gender child;
  model distance = gender gender*age / noint solution ;
  repeated / type = cs subject = child r rcorr;
run;

title "(d) COMMON AR(1) STRUCTURE";
proc mixed data=dent1 method=ml;
  class gender child ;
  model distance = gender age*gender / noint solution chisq;
  repeated / type = ar(1)  subject=child r rcorr;
run;

title "(e) COMMON ONE-DEPENDENT STRUCTURE";
proc mixed  data=dent1 method=ml;
  class gender child ;
  model distance = gender age*gender / noint solution chisq;
  repeated / type = toep(2)  subject=child r rcorr;
run;

title "(f) COMMON HETEROGENEOUS COMPOUND SYMMETRY STRUCTURE";
proc mixed data=dent1 method=ml;
  class gender child;
  model distance = gender gender*age / noint solution chisq;
  repeated / type = csh subject = child r rcorr;
run;

title "(g) SEPARATE COMPOUND SYMMETRY FOR EACH GENDER";
proc mixed  data=dent1 method=ml;
  class gender child ;
  model distance = gender age*gender / noint solution chisq;
  repeated / type = cs subject=child r=1,12 rcorr=1,12 group=gender;
run;

title "(h) SEPARATE AR(1) FOR EACH GENDER";
proc mixed  data=dent1 method=ml;
  class gender child ;
  model distance = gender age*gender / noint solution chisq;
  repeated / type = ar(1)  subject=child r=1,12 rcorr=1,12 group=gender;
run;

title "(i) SEPARATE ONE-DEPENDENT FOR EACH GENDER";
proc mixed data=dent1 method=ml;
  class gender child;
  model distance = gender age*gender / noint solution chisq;
  repeated / type = toep(2) subject=child r=1,12 rcorr=1,12 group=gender;
run;

title "(j) SEPARATE HETEROGENEOUS COMPOUND SYMMETRY STRUCTURE";
proc mixed  data=dent1 method=ml;
  class gender child;
  model distance = gender gender*age / noint solution chisq;
  repeated / type = csh subject = child r=1,12 rcorr=1,12 group=gender;
run;

/*******************************************************************

    AIC and BIC for these fits are as follows:

              AIC      BIC    
    (a)      447.5    465.6
    (b)      443.0    474.1
    (c)      440.6    448.4
    (d)      452.7    460.5
    (e)      469.4    477.2
    (f)      444.7    456.4
    (g)      424.8*   435.2*
    (h)      431.4    441.8
    (i)      460.6    471.0
    (j)      432.3    450.5
    
    Examination of the AIC, BIC suggests that a compound symmetric 
    that is different for each gender with different constant variance 
    across time for each gender is preferred.  

    This model is adopted in further analyes shown below.  We use ML
    to fit a full model with different slopes and reduced model with
    same slope, so that the likelihood ratio test statistic can be calculatd.
    The COVB option prints out the covariance matrix of the fixed
    effects estimates.  Here, the TESTS OF FIXED EFFECTS test the null hypotheses
    of both intercepts = 0 and both slopes = 0.
    
    We then reparameterize the model to get the Wald and F tests for
    differences in intercept and slope.  Here, the TESTS OF FIXED EFFECTS test
    the null hypotheses of no difference in intercepts and slopes between
    genders.

    
*******************************************************************/

*  full model again with covariance matrix of betahat printed;
 
title "FULL MODEL WITH COMPOUND SYMMETRY FOR EACH GENDER";
proc mixed method=ml data=dent1;
  class gender child;
  model distance = gender gender*age  / noint solution covb;
  repeated / type=cs subject=child r=1,12 rcorr=1,12 group=gender;
run; 

*  reduced model;

title "REDUCED MODEL WITH COMPOUND SYMMETRY FOR EACH GENDER";
proc mixed method=ml data=dent1;
  class gender child;
  model distance = gender age  / noint solution covb;
  repeated / type=cs subject=child r=1,12 rcorr=1,12 group=gender;
run; 

*  full model using REML;
*  use ESTIMATE statement to estimate the mean for a boy of age 11;
*  use CONTRAST statement to test whether slopes coincide;
*  use CONTRAST statement to test whether lines coincide;
*  CHISQ option gets Wald test in addition to F test;

title "FULL MODEL WITH COMPOUND SYMMETRY FOR EACH GENDER, REML";
proc mixed data=dent1;
  class gender child;
  model distance = gender gender*age  / noint solution covb;
  repeated / type=cs subject=child r=1,12 rcorr=1,12 group=gender;
  estimate 'boy at 11' gender 0 1 gender*age 0 11;
  contrast 'diff in slp' gender 0 0 gender*age 1 -1 /chisq;
  contrast 'both diff' gender 1 -1 gender*age 0  0,
                       gender 0  0 gender*age 1 -1 / chisq;
run; 

*  also fit full model in alt parameterization and use CHISQ option;
*  to get Wald test in addition to F test of differences in intercept;
*  (gender) and slope (age*gender);

title "FULL MODEL, ALTERNATIVE PARAMETERIZATION";
proc mixed data=dent1;
  class gender child;
  model distance = gender age gender*age  /  solution chisq covb;
  repeated / type=cs subject=child r=1,12 rcorr=1,12 group=gender;
run; 
