/******************************************************************

  CHAPTER 9, EXAMPLE 5, Epileptic Seizure Study

  Fit a subject-specific (generalized linear mixed effects) loglinear model 
  to the epileptic seizure data of Thall and Vail (1990) using 
  the GLIMMIX macro (not to be confused with PROC GLIMMIX).  This 
  macro implements a method that expands about current empirical
  Bayes estimates of random effects and is roughly similar to
  PQL.  It also implements a method that expands about zero and is
  similar to MQL. 
    
  See the macro file glmm800.sas and the SAS website (search on "glimmix
  macro" in customer resources) for details on options and usage.  The
  macro relies on linearizing the model at each iteration so that
  the linearized model can be fitted using PROC MIXED.
    
******************************************************************/

options ls=80 ps=59 nodate; run;

%inc 'glmm800.sas'  / nosource;

/******************************************************************

    The data look like (first 8 records on first 2 subjects)

       104  11   0   0  11  31
       104   5   1   0  11  31
       104   3   2   0  11  31
       104   3   3   0  11  31
       104   3   4   0  11  31
       106  11   0   0  11  30
       106   3   1   0  11  30
       106   5   2   0  11  30
       106   3   3   0  11  30
       106   3   4   0  11  30

  column 1      subject
  column 2      number of seizures
  column 3      visit (baseline 0, 1--4 biweekly visits)
  column 4      = 0 if placebo, = 1 if progabide
  column 5      baseline number of seizures in 8 weeks prior to study
  column 6      age

  Note that the baseline seizures are included in the dataset as both
  a separate variable and as the response at time 0.

******************************************************************/

data seizure; infile 'seize.dat';
    input subject seize visit trt base age;
*    if subject=207 then delete;
run;

*  create additional indicator variables;

data seizure; set seizure;
   logage=log(age);
   o=2; v=1;
   if visit=0 then o=8;
   if visit=0 then v=0;
   logo=log(o);
   visit4=1;
   if visit<4 then visit4=0;   
run;

/*****************************************************************

  Here, we fit two models: 

  - population model contains an intercept with random effect and 
    no dependence on treatment; the linear (visit) term depends on 
    treatment but has not random effect
    
  - a fancier model as above but with with random effects both
    intercept and linear (visit) term 
    
  See the macro for information on syntax, required statements, 
  and options.  See the PROC MIXED documentation for information
  on specifying the model using the "stmts=" statement.  The default 
  method is a version of PQL.  By adding

  options = MQL 

  you can get MQL, which linearizes the model about the random effects  
  = 0.

  The macro creates several data sets, including the random effects  
  estimates.  See the macro.

  The user should check the .log file for information on convergence
  (or not) of the algorithm.  
    
******************************************************************/

title "RANDOM INTERCEPT AND TIME";
%glimmix(data=seizure,
  procopt=method=ml,
  stmts=%str(
     class subject;
     model seize = v trt trt*v / solution;
     random intercept v / subject=subject type=un g gcorr;
  ),
  error=poisson,
  link=log,
  offset=logo     
);
run;


title "RANDOM INTERCEPT AND TIME, USING MQL";
%glimmix(data=seizure,
  procopt=method=ml,
  stmts=%str(
     class subject;
     model seize = v trt trt*v / solution;
     random intercept v / subject=subject type=un g gcorr;
  ),
  error=poisson,
  link=log,
  offset=logo,
  options=MQL      
  );
run;

title "RANDOM INTERCEPT ONLY";
%glimmix(data=seizure,
  procopt=method=ml,
  stmts=%str(
     class subject;
     model seize = v trt trt*v / solution;
     random intercept / subject=subject type=un;
  ),
  error=poisson,
  link=log,
  offset=logo     
);
run;

title "RANDOM INTERCEPT ONLY USING MQL";
%glimmix(data=seizure,
  procopt=method=ml,
  stmts=%str(
     class subject;
     model seize = v trt trt*v / solution;
     random intercept / subject=subject type=un;
  ),
  error=poisson,
  link=log,
  offset=logo,
  options=MQL     
);
run;






