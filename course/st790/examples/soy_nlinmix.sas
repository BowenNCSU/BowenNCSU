/******************************************************************

  CHAPTER 9, Soybean study

  Fit a subject-specific nonlinear mixed effects model to the
  soybean data using the nlinmix macro.  This macro implements approximate
  maximum likelihood based on models that are linearized about 
  the empirical Bayes estimates of the random effects and about
  zero.
    
  See the macro file nlmm801.sas and the SAS website (search on "nlinmix
  macro" in customer resources) for details on options and usage.  The
  macro relies on linearizing the model at each iteration so that
  the linearized model can be fitted using PROC MIXED.
    
******************************************************************/
    
options ps=55 ls=80 nodate;

%inc 'nlmm801.sas' / nosource;

/******************************************************************

    The data on the first 2 plots are 

    1 1988F1 F 1988 14 0.106
    1 1988F1 F 1988 21 0.261
    1 1988F1 F 1988 28 0.666
    1 1988F1 F 1988 35 2.11
    1 1988F1 F 1988 42 3.56
    1 1988F1 F 1988 49 6.23
    1 1988F1 F 1988 56 8.71
    1 1988F1 F 1988 63 13.35
    1 1988F1 F 1988 70 16.3417
    1 1988F1 F 1988 77 17.75083
    2 1988F2 F 1988 14 0.104
    2 1988F2 F 1988 21 0.269
    2 1988F2 F 1988 28 0.778
    2 1988F2 F 1988 35 2.12
    2 1988F2 F 1988 42 2.93
    2 1988F2 F 1988 49 5.29
    2 1988F2 F 1988 56 9.5
    2 1988F2 F 1988 70 16.9667
    2 1988F2 F 1988 77 17.7467

    column 1      plot id (numeric, 48 plots total)
    column 2      alternate plot id
    column 4      year
    column 5      day
    column 6      average leaf weight/plant

    NOTE:  These are the data from R, so rounded to 5 decimal places.
    
******************************************************************/

data soy; infile 'soyR.dat';
  input plot id $ geno $ year day lwt;
run;

data soy; set soy;
  d1=0; if geno="F" and year=1988 then d1=1;
  d2=0; if geno="P" and year=1988 then d2=1;
  d3=0; if geno="F" and year=1989 then d3=1;
  d4=0; if geno="P" and year=1989 then d4=1;
  d5=0; if geno="F" and year=1990 then d5=1;
  d6=0; if geno="P" and year=1990 then d6=1;
run;

/*****************************************************************

  First fit the logistic growth model with no covariates.

  Specification of the model is self-evident and is done in the MODEL
  section - the final expression for f(z_ij,beta_i) should be in the
  variable PREDV.

  The user can specify analytical derivatives in the DERIVS section 
  if desired; see the documentation.  We use numerical derivatives  
  here.

  The WEIGHT section allows specification of a within-individual 
  variance model, but it does not support estimation of unknown 
  parameters delta in the variance function.  We use the estimated
  value obtained from nlme(), which do allow this.
    
  The PARMS section is for specifying starting values for beta

  The STMTS section calls PROC MIXED to fit the linearized model
  at each interation.  The response variable must be the response
  variable in the data set with "PSEUDO_" in front.  The MODEL
  statement includes the fixed effects, appended by "D_" and the 
  RANDOM statement includes the random effects similarly. 

  EXPAND=ZERO does the first order linearization about zero;
  EXPAND=EBLUP expands about current empirical Bayes estimates of
  random effects.

  The PROCOPT section allows various options to be specified; see
  the documentation in the macro.

   It would be possible to extract the estimated random effects from 
   this fit and plot them against the covariates genotype and year,
   as in the nlme() fits; see the macro documentation.  We do not do 
   this here. 
      
   Note from the .log file that the estimated D matrix is not   
   positive definite.  As in the R fits, it appears that the among-plot
   variation in some of the parameters may be negligible relative to 
   others.   The fit of the model with covariates with random effects 
   for all parameters, not shown here, ends up with similar problems.
   Thus, below,  n fits with covariates, we take the random effect
   associated with the "half life" and growth rate constant parameters
   beta2i and beta3i to be negligible.
      
******************************************************************/

title "NO COVARIATES, LINEARIZATION ABOUT EMPIRICAL BAYES ESTIMATES";
%nlinmix(data=soy,
   model=%str(
     beta1i=beta1+b1;
     beta2i=beta2+b2;
     beta3i=beta3+b3;
     denom=1+exp(-beta3i*(day-beta2i));
     predv=beta1i/denom;
     ),
   weight=%str(
     _weight_= 1/predv**(2*0.95);
   ),
   parms=%str(
       beta1=20 beta2=50 beta3=0.125),
   stmts=%str(
       class plot;
       model pseudo_lwt = d_beta1 d_beta2 d_beta3 /
                          noint notest solution;
      random d_b1 d_b2 d_b3 / subject=plot type=un solution;
      weight _weight_;
   ),
   expand=eblup,
       procopt=%str(maxiter=500 method=ml)
)
run;

title "COVARIATES, NO RANDOM EFFECT FOR HALF LIFE OR GROWTH RATE";   
%nlinmix(data=soy,
   model=%str(
     beta1i=beta11*d1+beta12*d2+beta13*d3+beta14*d4+beta15*d5+beta16*d6+b1;
     beta2i=beta21*d1+beta22*d2+beta23*d3+beta24*d4+beta25*d5+beta26*d6;
     beta3i=beta31*d1+beta32*d2+beta33*d3+beta34*d4+beta35*d5+beta36*d6;
     denom=1+exp(-beta3i*(day-beta2i));
     predv=beta1i/denom;
     ),
   weight=%str(
     _weight_= 1/predv**(2*0.95);
   ),
   parms=%str(beta11=20 beta12=20 beta13=20 beta14=20 beta15=20 beta16=20
              beta21=50 beta22=50 beta23=50 beta24=50 beta25=50 beta26=50
              beta31=0.125 beta32=0.125 beta33=0.125 beta34=0.125 beta35=0.125 beta36=0.125),
   stmts=%str(
       class plot;
       model pseudo_lwt = d_beta11 d_beta12 d_beta13 d_beta14 d_beta15 d_beta16
                          d_beta21 d_beta22 d_beta23 d_beta24 d_beta25 d_beta26
                          d_beta31 d_beta32 d_beta33 d_beta34 d_beta35 d_beta36 /
                          noint notest solution;
      random d_b1  / subject=plot type=un solution;
      weight _weight_;
   ),
   expand=eblup,
   procopt=%str(maxiter=500 method=ml)
)
run;

 
