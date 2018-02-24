/*******************************************************************

   Cdystonia data - three weightloss programs

   0 = Placebo
   5000 = 5000 units botulinum toxin type B (BotB)
   10000 = 10000 units BotB

  Population averaged model analysis using proc mixed
  using a spline model 

*******************************************************************/

options ls=80 ps=59 nodate; run;

/******************************************************************

  The data set is in "long" format, so can be read in directly
   
*******************************************************************/

data cdystonia; infile "cdystonia.dat";
    input id week dose age gender twstrs;
run;

data cdystonia; set cdystonia;
    weekplus=week-4;
    if weekplus<0 then weekplus=0;
    time=week;
run;


/******************************************************************

  Assume a spline model - we allow the intercepts to vary here, but
  you could have taken them to be the same given that this is a
  randomized trial.  Because there are missing data here, we are
  careful to use maximum likelihood.  
    
  Fit various covariance models.  We first fit a model that allows a
  different unstructured covariance matrix for each dose group to
  look at the pattern of correlation and variance.  Remember that 
  in proc mixed, we can fit unstructured covariance that is either same or
  different in each group, but we cannot force the variances to be the same
  over time.  We fit unstructured models that are different for each group and
  then the same for each group. 

  The suggestion is that if there is a structure, it appears approximately compound
  symmetric, with possible different correlation for each group and
  variances that possibly are diffrent over time.  So we fit several
  different compound symmetric models with and without variance
  changing over time and common or not by dose.  

    The AIC values from the fits below, in order of appearance, are

                     AIC (Smaller is Better)        4107.8
                     AIC (Smaller is Better)        4109.0
                     AIC (Smaller is Better)        4138.0
                     AIC (Smaller is Better)        4140.7
                     AIC (Smaller is Better)        4160.4
                     AIC (Smaller is Better)        4161.0
   
    and the BIC values are

                     BIC (Smaller is Better)        4301.5
                     BIC (Smaller is Better)        4189.8
                     BIC (Smaller is Better)        4218.7
                     BIC (Smaller is Better)        4183.7
                     BIC (Smaller is Better)        4200.7
                     BIC (Smaller is Better)        4190.6

    AIC and BIC do not agree - AIC favors larger models and prefers  
    separate unstructured matrices for each dose group, while BIC
    prefers the common compound symmetric specification with variances
    that are different at different time points.  For parsimmony, we
    choose to go with BIC's choice of the common compound symmetric 
    specification with  variance changing over time. This model can 
    be fitted in gls() in R, too.  Because R and SAS define BIC
    differently, SAS prefers this model, while R prefers the common
    compound symmetric model with the same variance at all time points.
    
    If one prefers to switch to R and gls() at this point, one can, using 
    this model (the separate unstructured model is not possible in gls()).
    
*******************************************************************/

title "UNSTRUCTURED SEPARATE BY DOSE"; 
proc mixed data=cdystonia method=ML; 
  class dose time; 
  model twstrs = dose dose*week dose*weekplus / noint solution ; 
  repeated time / type = un subject = id r=1,6,8 rcorr=1,6,8 group=dose; 
run;

title "COMMON UNSTRUCTURED";
proc mixed data=cdystonia method=ML;
  class dose time;
  model twstrs = dose dose*week dose*weekplus / noint solution ;
  repeated time / type = un subject = id r rcorr;
run;

title "HETEROGENEOUS COMPOUND SYMMETRIC SEPARATE BY DOSE";
proc mixed data=cdystonia method=ML;
  class dose time;
  model twstrs =dose dose*week dose*weekplus / noint solution ;
  repeated time / type = csh subject = id r=1,6,8 rcorr=1,6,8 group=dose;
run;

title "COMMON HETEROGENEOUS COMPOUND SYMMETRIC";
proc mixed data=cdystonia method=ML;
  class dose time;
  model twstrs = dose dose*week dose*weekplus / noint solution ;
  repeated time / type = csh subject = id r rcorr;
run;

title "COMPOUND SYMMETRIC SEPARATE BY DOSE";
proc mixed data=cdystonia method=ML;
  class dose time;
  model twstrs = dose dose*week dose*weekplus / noint solution ;
  repeated time / type = cs subject = id r=1,6,8 rcorr=1,6,8 group=dose;
run;

title "COMMON COMPOUND SYMMETRIC";
proc mixed data=cdystonia method=ML;
  class dose time;
  model twstrs =dose dose*week dose*weekplus / noint solution ;
  repeated time / type = cs subject = id r rcorr;
run;

/******************************************************************

  Adopting the common compound symmetric structure with variances
  that differ over time, we fit again. We don't request robust
  sandwich standard errors given the missing data. Again, we leave
  in separate intercepts for each treatment; you may have
  simplified the model further to have a common intercept.
    
*******************************************************************/

title "ANALYSES WITH FINAL MODEL";
proc mixed data=cdystonia method=ML;
  class dose time;
  model twstrs = dose dose*week dose*weekplus / noint solution ;
  repeated time / type = csh subject = id;
  
*  Tests of common intercept and 1st and second phase slopes;

  contrast 'common intercept' dose 1 -1 0,
                              dose 1 0 -1  / chisq;

  contrast 'common 1st stage slp'  week*dose 1 -1 0,
                                   week*dose 1 0 -1 /chisq;

  contrast 'common 2nd stage slp'  week*dose 1 -1 0 weekplus*dose 1 -1 0,
                              week*dose 1 0 -1 weekplus*dose 1 0 -1 / chisq;

*  Test of different means at week 16;
  
  contrast 'diff means @ week 16'
      dose 1 -1 0 week*dose 16 -16 0 weekplus*dose 12 -12 0,
      dose 1 0 -1 week*dose 16 0 -16 weekplus*dose 12 0 -12 /chisq;

*  Estimates of slope in second stage;

  estimate '2nd stg slp, placebo' week*dose 1 0 0 weekplus*dose 1 0 0;
  estimate '2nd stg slp, 5000 U'  week*dose 0 1 0 weekplus*dose 0 1 0;
  estimate '2nd stg slp, 10000 U' week*dose 0 0 1 weekplus*dose 0 0 1;

*  Estimates of mean at end of study (week 16);
  
  estimate 'mean @ wk 16, placebo' dose 1 0 0 week*dose 16 0 0 weekplus*dose
      12 0 0 ;
  estimate 'mean @ wk 16, 5000 U'  dose 0 1 0 week*dose 0 16 0 weekplus*dose
      0 12 0;
  estimate 'mean @ wk 16, 10000 U' dose 0 0 1 week*dose 0 0 16  weekplus*dose
      0 0 12;
run;

title "ADD AGE AND GENDER, DOSE-SPECIFIC";
proc mixed data=cdystonia method=ML;
  class dose time;
  model twstrs = dose dose*gender dose*age dose*week dose*week*gender
      dose*week*age dose*weekplus / noint solution ;
  repeated time / type = csh subject = id;
run;

title "ADD AGE AND GENDER, NOT DOSE-SPECIFIC";
proc mixed data=cdystonia method=ML;
  class dose time;
  model twstrs = dose gender age dose*week week*gender
      week*age dose*weekplus / noint solution ;
  repeated time / type = csh subject = id;
run;
