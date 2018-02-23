/******************************************************************

  CHAPTER 8, EXAMPLE 5, Epileptic Seizure Study

  Fit a population-averaged loglinear model to the epileptic
  seizure data of Thall and Vail (1990) using several SAS procedures:

  PROC GENMOD and PROC GEE, for which the working correlation matrix
  parameter alpha is estimated using moment methods.  

  PROC GLIMMIX can also be used without random effects to fit PA
  models with working correlation parameter alpha estimated by
  solving the quadratic estimating equation with Gaussian working
  assumption.  
    
  These are count data, so we use the Poisson mean/variance
  assumptions and allow for overdispersion.  The model is fitted with
  several different working assumptions on the overall correlation
  matrix:  completely unstructured, compound symmetric, and AR(1).

  Subject 207 has what appear to be very unusual data -- for
  this subject, both baseline and study-period numbers of seizures
  are huge, much larger than for any other subject.  We do not delete 
  207 for the analyses here, but many authors have done so.  See
  Diggle, Heagerty, Liang, and Zeger (2002) and Thall and Vail (1990)
  for more on this subject.

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

  We fit the basic and modified models presented in Section 8.8 of
  the course notes. The definitions of o, v, and v4 correspond to those 
  in that section.

  We first use PROC GENMOD and PROC GEE.  The syntax for each is the same.  
  
  The DIST=POISSON option in the model statement specifies
  that the Poisson requirement that mean = variance, be used.  
  The LINK=LOG option asks for the loglinear model.  Other 
  LINK= choices are available.

  The REPEATED statement specifies the "working" correlation
  structure to be assumed.    The CORRW option in the REPEATED 
  statement prints out the estimated working correlation matrix
  under the assumption given in the TYPE= option.  By default,  
  robust/empirical standard errors are used; the MODELSE 
  option specifies that the standard errors based on assuming
  that the working correlation nodel is correctly specified are also
  printed.  The REPEATED statement has an option WITHIN= to align the
  time points correctly in the event of unbalanced/missing data; we
  illustrate this here even though we do not need to use it, as the data
  are balanced.
    
  The dispersion parameter is estimated rather then being held 
  fixed at 1 -- this allows for the possibility of overdispersion.

  The OFFSET option implements the division by the length of the time 
  period, as discussed in the course notes.  In principle, this could 
  also be accomplished by dividing the raw counts by 8 or 2 as appropriate;
  however, many programs for solving GEEs or fitting generalized linear
  models restrict the response to be an integer if the distribution is
  specified as POISSON, so this doesn't necessarily always work.  
      
******************************************************************/

* proc genmod, fits with three different working correlation models;

title "PROC GENMOD, UNSTRUCTURED CORRELATION";
proc genmod data=seizure;
  class subject visit;
  model seize = v trt trt*v /  dist = poisson link = log offset=logo;
  repeated subject=subject / within=visit type=un corrw covb modelse;
run;

title "PROC GENMOD, COMPOUND SYMMETRIC CORRELATION";
proc genmod data=seizure;
  class subject visit;
  model seize = v trt trt*v / dist = poisson link = log offset=logo;
  repeated subject=subject / within=visit type=cs corrw covb modelse;
run;

title "PROC GENMOD, AR(1) CORRELATION";
proc genmod data=seizure;
  class subject visit;
  model seize = v trt trt*v / dist = poisson link = log offset=logo;
  repeated subject=subject / within=visit  type=ar(1) corrw covb modelse;
run;

*  fancier model, compound symmetric only;

title "PROC GENMOD, FANCIER MODEL";
proc genmod data=seizure;
  class subject visit;
  model seize = logage v visit4 trt trt*v trt*visit4 /
                 dist = poisson link = log offset=logo;
  repeated subject=subject / within=visit type=cs corrw covb modelse;
run;

*  proc gee;

title "PROC GEE, COMPOUND SYMMETRIC CORRELATION";
proc gee data=seizure;
  class subject visit;
  model seize = v trt trt*v /  dist = poisson link = log offset=logo;
  repeated subject=subject / within=visit  type=cs corrw covb modelse;
run;

title "PROC GEE, FANCIER MODEL";
proc gee data=seizure;
  class subject visit;
  model seize = logage v visit4 trt trt*v trt*visit4 /
                 dist = poisson link = log offset=logo;
  repeated subject=subject / within=visit type=cs corrw covb modelse;
run;

/*****************************************************************

  PROC GENMOD and PROC GEE use moment methods to estimate the 
  correlation parameters.  We demonstrate the use of PROC GLIMMIX to
  instead estimate the correlation parameters using the quadratic
  estimating equation.  This procedure is really meant for fitting    
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

title "PROC GLIMMIX, COMPOUND SYMMETRIC CORRELATION";
proc glimmix data=seizure method=mmpl empirical;
  class subject visit;
  model seize = v trt trt*v / dist = poisson link = log offset=logo
      covb solution;
  random visit / subject=subject type=cs vcorr residual;

*  Alternatively, the following statement is equivalent;
*  random _residual_ / subject=subject type=cs vcorr;
run;

title "PROC GLIMMIX, FANCIER MODEL";
proc glimmix data=seizure method=mmpl empirical;
  class subject visit;
  model seize = logage v visit4 trt trt*v trt*visit4 / dist = poisson link = log offset=logo
      covb solution;
  random visit / subject=subject type=cs vcorr residual;
*  random _residual_ / subject=subject type=cs vcorr;
run;

