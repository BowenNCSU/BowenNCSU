/*******************************************************************

   ACTG 193 CD4 data

   1 = zidovudine alternating monthly with 400mg didanosine
   2 = zidovudine plus 2.25mg of zalcitabine
   3 = zidovudine plus 400mg of didanosine
   4 = zidovudine plus 400mg of didanosine plus 400mg of nevirapine

  These analyses take 1,2,3 to be "dual" therapy and combine the
  three groups, and compare to "triple" therapy 4
    
  Subject-specific model analysis using proc mixed

*******************************************************************/

options ls=80 ps=59 nodate; run;

/******************************************************************

  The data set is in "long" format, which is what proc mixed wants
   
*******************************************************************/

data thedat; infile "cd4.dat";
    input id trt age gender week logcd4;
run;

data thedat; set thedat;
    weekplus=week-16;
    if weekplus<0 then weekplus=0;
    dual=0; if trt=4 then dual=1;
run;

/******************************************************************

  We assume a spline model - preliminary investigation shows no
  evidence that the intercept is different across the two groups,
  consistent with the fact this is a randomized study, so we adopt
  a common intercept.  We allow a random effect for each of
  intercept, phase 1 slope, and phase 2 change from phase 1 slope.

  Tbe five models that converged below have, in order of appearance

                     AIC (Smaller is Better)       11912.1
                     AIC (Smaller is Better)       11913.6
                     AIC (Smaller is Better)       11907.3
                     AIC (Smaller is Better)       11912.1
                     AIC (Smaller is Better)       11914.1

                     BIC (Smaller is Better)       11974.2
                     BIC (Smaller is Better)       11980.9
                     BIC (Smaller is Better)       12005.6
                     BIC (Smaller is Better)       11974.2
                     BIC (Smaller is Better)       11981.4

  AIC perfers the different D/different within variance model, 
  while BIC prefers the common D/common within variance model.
  On the basis of parsimony, recognizing that AIC often chooses
  larger models, we will go with the simpler common D/variance
  model for further analyses.
    
*******************************************************************/

title "COMMON D, COMMON DIAGONAL WITHIN"; 
proc mixed data=thedat method=ML; 
  class dual id; 
  model logcd4 = dual*week dual*weekplus /  solution ; 
  random intercept week weekplus / type=un subject=id g gcorr v=1 vcorr=1;
run;

title "COMMON D, DIFFERENT DIAGONAL WITHIN"; 
proc mixed data=thedat method=ML; 
  class dual id ; 
  model logcd4 = dual*week dual*weekplus /  solution ; 
  random intercept week weekplus / type=un subject=id g gcorr v=1 vcorr=1;the
  repeated / group=dual subject=id;
run;

/*

This fit would not converge

title "DIFFERENT D, SAME DIAGONAL WITHIN"; 
proc mixed data=thedat method=ML; 
  class id; 
  model logcd4 = dual dual*week dual*weekplus /  solution ; 
  random intercept week weekplus / type=un group=dual subject=id g gcorr v=1 vcorr=1;
run;
*/

title "DIFFERENT D, DIFFERENT DIAGONAL WITHIN"; 
proc mixed data=thedat method=ML; 
  class dual id ; 
  model logcd4 =  dual*week dual*weekplus / solution ; 
  random intercept week weekplus / type=un group=dual subject=id g gcorr v=1 vcorr=1;
  repeated / group=dual subject=id;
run;

title "COMMON D, COMMON EXPONENTIAL WITHIN"; 
proc mixed data=thedat method=ML; 
  class dual id ; 
  model logcd4 =  dual*week dual*weekplus /  solution ; 
  random intercept week weekplus / type=un subject=id g gcorr v=1 vcorr=1;
  repeated / type=sp(exp)(week) subject=id;
run;

title "COMMON D, COMMON EXPONENTIAL WITHIN + MEAS ERROR"; 
proc mixed data=thedat method=ML; 
  class dual id ; 
  model logcd4 =  dual*week dual*weekplus /  solution ; 
  random intercept week weekplus / type=un subject=id g gcorr v=1 vcorr=1;
  repeated / type=sp(exp)(week) local subject=id;
run;

/*

Could also try these, but given the lack of evidence of
within individual correlation, don't bother.

title "DIFFERENT D, COMMON EXPONENTIAL WITHIN"; 
proc mixed data=thedat method=ML; 
  class dual id ; 
  model logcd4 =  dual*week dual*weekplus /  solution ; 
  random intercept week weekplus / group=dual type=un subject=id g gcorr v=1 vcorr=1;
  repeated / type=sp(exp)(week) subject=id;
run;

title "DIFFERENT D, COMMON EXPONENTIAL WITHIN + MEAS ERROR"; 
proc mixed data=thedat method=ML; 
  class dual id ; 
  model logcd4 =  dual*week dual*weekplus / solution ; 
  random intercept week weekplus /  group=dual type=un subject=id g gcorr v=1 vcorr=1;
  repeated  / type=sp(exp)(week) local subject=id;
run;

*/

/******************************************************************

  We adopt the common D, common within-individual variance, no
  within-individual correlation. Given the missing data, we need
  to recognize that standard errors and tests based on the usual
  asymptotic theory might be unreliable; likelihood ratio tests
  of nested models might be a better approach.  Here, we fit three
  different population models, allowing intercept, phase 1 slope
  and phase 2 difference from phase 1 slope to depend on age and
  gender.  From this series of models, it seems gender doesn't
  matter, but age does.  
    
*******************************************************************/

title "INTERCEPT AGE AND GENDER"; 
proc mixed data=thedat method=ML; 
  class dual id dual; 
  model logcd4 = age gender dual*week dual*weekplus / solution ; 
  random intercept week weekplus / type=un subject=id g gcorr;
run;

title "INTERCEPT AND PHASE 1 AGE AND GENDER"; 
proc mixed data=thedat method=ML; 
  class dual id dual; 
  model logcd4 =  age gender dual*week dual*week*age dual*week*gender dual*weekplus / solution ; 
  random intercept week weekplus / type=un subject=id g gcorr;
run;

title "AGE AND GENDER EVERYWHERE"; 
proc mixed data=thedat method=ML; 
  class dual id dual; 
  model logcd4 = age gender dual*week dual*week*age dual*week*gender dual*weekplus
                 dual*week*age dual*weekplus*gender / solution ; 
  random intercept week weekplus / type=un subject=id g gcorr;
run;

title "AGE EVERYWHERE"; 
proc mixed data=thedat method=ML; 
  class dual id dual; 
  model logcd4 = age  dual*week dual*week*age dual*weekplus dual*weekplus*age / solution ; 
  random intercept week weekplus / type=un subject=id g gcorr;
run;

title "INTERCEPT AGE"; 
proc mixed data=thedat method=ML; 
  class dual id dual; 
  model logcd4 = age dual*week dual*weekplus / solution ; 
  random intercept week weekplus / type=un subject=id g gcorr;
run;
