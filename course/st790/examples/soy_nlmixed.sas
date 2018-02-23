/******************************************************************

  CHAPTER 9, Soybean study

  Fit a subject-specific nonlinear mixed effects model to the
  soybean data using PROC NLIMIXED.  he procedure implments the model using
  the first order linearization about zero and solving GEE-2 
  equations (METHOD=FIRO) and adaptive Gaussian quadature (the default).
      
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

 PROC NLMIXED is very general, so requires a lot of intervention
   by the user to set up the model.  It is necessary to specify 
   starting values for all parameters.  Starting values for the fixed
   effects might be found by fitting the model by OLS taking all 
   observations to be independent, and then using the NLINMIX macro
   or R nlme() to get starting values for the covariance parameters.
   
   The default integation method is adaptive Gaussian quadrature.  The
   number of quadrature points (abscissae) is calculated adaptively.  If
   you want to specify a fixed number of points L, specify

   METHOD=GAUSS QPOINTS=L

   in the PROC NLMIXED statement.  The full Laplace method with no 
   with no further linearization with QPOINTS=1.  The quadrature
   can be made nonadaptive with the NOAD option. 
     
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

   Here, we fit only the model with covariates with no random effects
   for the half-life and growth rate parameters; see the results using 
   R nlme() and the NLINMIX macro.
    
******************************************************************/

title "COVARIATES, NO RANDOM EFFECT FOR HALF LIFE OR GROWTH RATE";   
proc nlmixed data=soy;
    parms beta11=20 beta12=20 beta13=20 beta14=20 beta15=20 beta16=20
          beta21=50 beta22=50 beta23=50 beta24=50 beta25=50 beta26=50
          beta31=0.125 beta32=0.125 beta33=0.125 beta34=0.125 beta35=0.125
          beta36=0.125 s12=1 s2=0.04 delta=1;
    beta1i=beta11*d1+beta12*d2+beta13*d3+beta14*d4+beta15*d5+beta16*d6+b1;
    beta2i=beta21*d1+beta22*d2+beta23*d3+beta24*d4+beta25*d5+beta26*d6;
    beta3i=beta31*d1+beta32*d2+beta33*d3+beta34*d4+beta35*d5+beta36*d6;
    denom=1+exp(-beta3i*(day-beta2i));
    pred=beta1i/denom;
    model lwt ~ normal(pred,(pred**(2*delta))*s2);
    random b1 ~ normal(0,s12) subject=plot;
run;

 
