 /******************************************************************

  CHAPTER 9, Argatroban pharmacokinetic study

  Fit a subject-specific nonlinear mixed effects model to the
  argatroban PK data in Section 9.5 of Davidian and Giltinan (1995)
  using the SAS NLINMIX macro.  This macro implements approximate
  maximum likelihood based on models that are linearized about 
  the empirical Bayes estimates of the random effects and about
  zero.
    
  See the macro file nlmm801.sas and the SAS website (search on "nlinmix
  macro" in customer resources) for details on options and usage.  The
  macro relies on linearizing the model at each iteration so that
  the linearized model can be fitted using PROC MIXED.
    
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

/*****************************************************************

  We fit the mode in the course notes.  Specification of the model
  is self-evident and is done in the MODEL section - the final
  expression for f(z_ij,beta_i) should be in the variable PREDV.

  The user can specify analytical derivatives in the DERIVS section 
  if desired; see the documentation.  We use numerical derivatives  
  here.

  The WEIGHT section allows specification of a within-individual 
  variance model, but it does not support estimation of unknown 
  parameters delta in the variance function.  We use the estimated
  value obtained from nlme() and proc nlmixed, which do allow this.
    
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
      
******************************************************************/

title "FIRST ORDER LINEARIZATION ABOUT ZERO";
%nlinmix(data=arg,
   model=%str( 
     cl=exp(beta1+b1);
     v=exp(beta2+b2);
     predv=(dose/cl)*(1-exp(-cl*t2/v))*exp(-cl*(1-t1)*(time-tinf)/v);
   ),
   weight=%str(
     _weight_= 1/predv**(2*0.24);
   ),
   parms=%str(beta1=-6.0 beta2=-2.0),
   stmts=%str(
      class indiv;
      model pseudo_conc = d_beta1 d_beta2 / noint notest solution;
      random d_b1 d_b2 / subject=indiv type=un solution;
      weight _weight_;
   ),
      expand=zero,
   procopt=%str(maxiter=500 method=ml)
)
run;

title "FIRST ORDER CONDITIONAL LINEARIZATION ABOUT EMPIRICAL BAYES ESTIMATES";   
%nlinmix(data=arg,
   model=%str( 
     cl=exp(beta1+b1);
     v=exp(beta2+b2);
     predv=(dose/cl)*(1-exp(-cl*t2/v))*exp(-cl*(1-t1)*(time-tinf)/v);
   ),
   weight=%str(
     _weight_=1/predv**(2*0.24);
   ),
   parms=%str(beta1=-6.0 beta2=-2.0),
   stmts=%str(
      class indiv;
      model pseudo_conc = d_beta1 d_beta2 / noint notest solution ;
      random d_b1 d_b2 / subject=indiv type=un solution;
      weight _weight_;
   ),
      expand=eblup,
   procopt=%str(maxiter=500 method=ml)
)
run;
