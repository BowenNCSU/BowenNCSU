/******************************************************************

  CHAPTER 9, EXAMPLE 5, Epileptic Seizure Study

  Fit a subject-specific (generalized linear mixed effects) loglinear model 
  to the epileptic seizure data of Thall and Vail (1990) using 
  PROC GLIMMIX (not to be confused with the GLIMMIX macro).

  This procedure implements several methods, which can be specified
  in the PROC GLIMMIX statement using METHOD=

  -  MMPL, RMPL -- these are versions of "MQL," where the model is linearized
     about the random effects = 0.  MMPL implements the linear GEE for beta
     and the quadratic equation for covariance parameters with the Gaussian
     working assumption; RMPL replaces the quadratic equation by a REML
     version

  -  MSPL, RSPL -- these are as MMPL and RMPL but where the model is
     linearized about the current empirical Bayes estimates of the random
     effets.  These are similar to what has been called "PQL."

  -  Laplace -- This is a full Laplace aproximation to the integrals
     in the likelihood without a further approximation to a linear mixed
     effects model as in PQL

  -  quad -- Adaptive Gaussian quadrature to "do" the integrals.  The number
     of abscissae (quadrature points) is calculated in a data-adaptive way
     (see the documentation) -- if you want to specify a fixed number L,
     use quad(qpoints=L)
    
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
run;

/*****************************************************************

  Here, we fit two models: 

  - population model contains an intercept with random effect and 
    no dependence on treatment; the linear (visit) term depends on 
    treatment but has not random effect
    
  - a fancier model as above but with with random effects both
    intercept and linear (visit) term 

  The RANDOM statement functions just as in PROC MIXED, so there are
  options G and GCORR to print out the estimated among-individual
  covariance matrix D and GROUP= to allow it to differ by levels of
  a grouping variable.
    
  See the proc glimmix documentation for more information on syntax,
  required statements, and options.
 
  Note that with a random intercept only several of the analyses give
  virtually identical estimates; see if you can figure out why.

  The empirical Bayes estimates of the random effects are printed
  if the SOLUTION option is included in the RANDOM statement as with 
  PROC MIXED.  These can be extracted using ODS as in PROC MIXED for
  further analysis.
    
******************************************************************/

title "RANDOM INTERCEPT AND TIME, ADAPTIVE QUADATURE";
proc glimmix data=seizure method=quad(qpoints=10);
  class subject;
  model seize = v trt trt*v / solution link=log 
        dist=poisson offset=logo;
  random int v / subject=subject type=un g gcorr;
run;

title "RANDOM INTERCEPT AND TIME, LAPLACE APPROXIMATION";
proc glimmix data=seizure method=laplace;
  class subject;
  model seize = v trt trt*v / solution link=log 
        dist=poisson offset=logo;
  random int v / subject=subject type=un g gcorr;
run;

title "RANDOM INTERCEPT AND TIME, LINEARIZATION ABOUT EB ESTIMATES";
proc glimmix data=seizure method=mspl;
  class subject;
  model seize = v trt trt*v / solution link=log 
        dist=poisson offset=logo;
  random int v / subject=subject type=un g gcorr;
run;

title "RANDOM INTERCEPT AND TIME, LINEARIZATION ABOUT 0";
proc glimmix data=seizure method=mmpl;
  class subject;
  model seize = v trt trt*v / solution link=log 
        dist=poisson offset=logo;
  random int v / subject=subject type=un g gcorr;
run;

title "RANDOM INTERCEPT ONLY, ADAPTIVE QUADATURE";
proc glimmix data=seizure method=quad(qpoints=25);
  class subject;
  model seize = v trt trt*v / solution link=log 
        dist=poisson offset=logo;
  random int / subject=subject type=un;
run;

title "RANDOM INTERCEPT ONLY, LAPLACE APPROXIMATION";
proc glimmix data=seizure method=laplace;
  class subject;
  model seize = v trt trt*v / solution link=log 
        dist=poisson offset=logo;
  random int / subject=subject type=un;
run;

title "RANDOM INTERCEPT ONLY, LINEARIZATION ABOUT EB ESTIMATES";
proc glimmix data=seizure method=mspl;
  class subject;
  model seize = v trt trt*v / solution link=log 
        dist=poisson offset=logo;
  random int / subject=subject type=un;
run;

title "RANDOM INTERCEPT ONLY, LINEARIZATION ABOUT 0";
proc glimmix data=seizure method=mmpl;
  class subject;
  model seize = v trt trt*v / solution link=log 
        dist=poisson offset=logo;
  random int / subject=subject type=un;
run;
