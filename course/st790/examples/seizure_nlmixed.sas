/******************************************************************

  CHAPTER 9, EXAMPLE 5, Epileptic Seizure Study

  Fit a subject-specific (generalized linear mixed effects) loglinear model 
  to the epileptic seizure data of Thall and Vail (1990) using PROC NLMIXED.

 ******************************************************************/

options ls=80 ps=59 nodate; run;

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
   trtv=trt*v; 
run;

/*****************************************************************

   Fit two models: 

   Here, we fit two models: 

  - population model contains an intercept with random effect and 
    no dependence on treatment; the linear (visit) term depends on 
    treatment but has not random effect
    
  - a fancier model as above but with with random effects both
    intercept and linear (visit) term 

   PROC NLMIXED is very general, so requires a lot of intervention
   by the user to set up the model.  Unlike PROC GLIMMIX and the
   GLIMMIX macro, which "automatically" calculate starting values 
   for beta by fitting a generalized linear mixed model assuming
   independence and starting values for D in various ways, it is 
   necessary to specify starting values for all parameters.  Starting
   values could be obtained by fitting the desired model using one of 
   these methods first.
    
   The default integation method is adaptive Gaussian quadrature.  The
   number of quadrature points (abscissae) is calculated adaptively.  If
   you want to specify a fixed number of points L, specify

   METHOD=GAUSS QPOINTS=L

   in the PROC NLMIXED statement.  The Laplace method can be obtained
   with QPOINTS=1.  The quadraure can be made nonadaptive with the
   NOAD option. 
     
   As shown below, the user has to specify starting values for all 
   fixed effects and covariance parameters in the PARMS statement.
   Programming statements allow the user to specify whatever model
   s/he wants, with any sort of nonlinearity.  The MODEL statement
   specifies the within-individual conditional distribution.  See 
   the documentation for the built-in conditional distributions that
   are possible - these assume within-individual independence of the
   Y_ij.  The user can specify his/her own within-individual conditional 
   likelihood with some programming.  The RANDOM statement is where one
   specifies the random effects - the syntax is pretty obvious.
    
   See the proc nlmixed documentation for information on syntax,
   required statements, and options.

   Compare the results to those from PROC GLIMMIX.
    
******************************************************************/

title "RANDOM INTERCEPT AND TIME, ADAPTIVE QUADRATURE";
proc nlmixed data=seizure method=gauss qpoints=10;
  parms b0=1 b1=0 b2=0.05 b3=-0.2 sb12=0.5 cb12=-0.05 sb22=0.2;
  eta = logo + b0 + b1*v + b2*trt + b3*trtv + b1i + b2i*v;
  f = exp(eta);
  model seize ~ poisson(f);
  random b1i b2i ~ normal([0,0],[sb12,cb12,sb22]) subject=subject;
run;

title "RANDOM INTERCEPT AND TIME, LAPLACE APPROXIMATION";
proc nlmixed data=seizure method=gauss qpoints=1;
  parms b0=1 b1=0 b2=0.05 b3=-0.2 sb12=0.5 cb12=-0.05 sb22=0.2;
  eta = logo + b0 + b1*v + b2*trt + b3*trtv + b1i + b2i*v;
  f = exp(eta);
  model seize ~ poisson(f);
  random b1i b2i ~ normal([0,0],[sb12,cb12,sb22]) subject=subject;
run;

title "RANDOM INTERCEPT ONLY, ADAPTIVE QUADRATURE";
proc nlmixed data=seizure method=gauss qpoints=25;
  parms b0=1 b1=0 b2=0.05 b3=-0.2 sb2=0.6;
  eta = logo + b0 + b1*v + b2*trt + b3*trtv + bi;
  f = exp(eta);
  model seize ~ poisson(f);
 random bi ~ normal(0,sb2) subject=subject;
run;

title "RANDOM INTERCEPT ONLY, LAPLACE APPROXIMATION";
proc nlmixed data=seizure method=gauss qpoints=1;
  parms b0=1 b1=0 b2=0.05 b3=-0.2 sb2=0.6;
  eta = logo + b0 + b1*v + b2*trt + b3*trtv + bi;
  f = exp(eta);
  model seize ~ poisson(f);
 random bi ~ normal(0,sb2) subject=subject;
run;

