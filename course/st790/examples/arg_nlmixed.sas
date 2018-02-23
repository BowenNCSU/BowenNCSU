/******************************************************************

  CHAPTER 9, Argatroban pharmacokinetic study

  Fit a subject-specific nonlinear mixed effects model to the
  argatroban PK data in Section 9.5 of Davidian and Giltinan (1995)
  using SAS PROC NLMIXED.  The procedure implments the model using
  the first order linearization about zero and solving GEE-2 
  equations (METHOD=FIRO) and adaptive Gaussian quadature (the default).
    
******************************************************************/
 
options ps=55 ls=80 nodate;

%inc 'nlmm801.sas'
     / nosource;

/******************************************************************

    The data look like (data on first subject)

    1 1 1 30 95.7
    2 1 1 60 122
    3 1 1 90 133
    4 1 1 160 162
    5 1 1 200 200
    6 1 1 240 172
    7 1 1 245 122
    8 1 1 250 120
    9 1 1 260 60.6
   10 1 1 275 70
   11 1 1 295 47.3

  column 1      observation no
  column 2      subject id
  column 3      infusion rate
  column 4      time (min)
  column 5      argatroban concentration

******************************************************************/

data arg; infile 'argconc.dat';
  input obsno indiv dose time conc;
  tinf=240;
  t1=1; if time>tinf then t1=0;
  t2=tinf*(1-t1)+t1*time;
run;
data arg; infile 'argconc.dat';
  input obsno id dose time conc;
run;

data arg; set arg;
  tinf=240;
  t1=1;
  if time>tinf then t1=0;
  t2=tinf*(1-t1)+t1*time;
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

   Compare the results to those from NLINMIX and nlme() - they are 
   somewhat different, particularly for estimation of the covariance
   parameters.
    
******************************************************************/

title "ARGATROBAN PK USING ADAPTIVE GAUSSIAN QUADRATURE";
proc nlmixed data=arg method=gauss qpoints=10;
  parms beta1=-6.0 beta2=-2.0 s2b1=0.14 cb12=0.006 s2b2=0.005 s2=20.0,delta=0.5;
  cl=exp(beta1+b1);
  v=exp(beta2+b2);
  pred=(dose/cl)*(1-exp(-cl*t2/v))*exp(-cl*(1-t1)*(time-tinf)/v);
  model conc ~ normal(pred,(pred**(2*delta))*s2);
  random b1 b2 ~ normal([0,0],[s2b1,cb12,s2b2])
     subject=id;
run;

*  Check the .log file -- this fit is not reliable;

title "ARGATROBAN PK USING FO METHOD WITH GEE-2";
proc nlmixed data=arg method=firo;
  parms beta1=-6.0 beta2=-2.0 s2b1=0.14 cb12=0.006 s2b2=0.005 s2=20.0,delta=0.5;
  cl=exp(beta1+b1);
  v=exp(beta2+b2);
  pred=(dose/cl)*(1-exp(-cl*t2/v))*exp(-cl*(1-t1)*(time-tinf)/v);
  model conc ~ normal(pred,(pred**(2*delta))*s2);
  random b1 b2 ~ normal([0,0],[s2b1,cb12,s2b2])
     subject=id;
run;
