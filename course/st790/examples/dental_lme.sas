 /*******************************************************************

  CHAPTER 6, EXAMPLE 1, Dental Study

  Subject-specific linear mixed effects model

  We use the RANDOM statement of PROC MIXED with TYPE=UN
  and possibly the GROUP= option to specify the among-individual
  covariance matrix D (possibly different by levels of GROUP)  

  We now use the REPEATED statement with various TYPE= specifications,
  the LOCAL option, and possibly the GROUP = option with SUBJECT = to 
  specify the within-individual covariance matrix R_i:

     - TYPE = specifies the form of the component of R_i due to the
       realization process; typically, one would choose a correlation
       structure that "damps out" over time.  Leaving out TYPE =
       fits R_i = sigma^2 I_{n_i} where sigma^2 can be made different
       by GROUP with the GROUP = option
     
     - LOCAL specifies that the user wishes to include a measurement
       error process component with constant variance over time.
       If LOCAL is not present, the overall form of R_i is determined 
       by the TYPE = specification

     - Leaving out the REPEATED statement produces the default
       R_i = sigma^2 I_{n_i}
     
*******************************************************************/

options ls=80 ps=59 nodate; run;

/******************************************************************

  Read in the data set

*******************************************************************/

data dent1; infile 'dental.dat';
  input obsno child age distance gender;
run;

data dent1; set dent1;
  time=age;
run;

/*******************************************************************

  Use PROC MIXED to fit the random coefficient model with different
  assumptions about among- and within-individual covariance
  struture.  We use ML -- leave off method=ML to get the default
  REML fit.  In all cases, we parametrize the random intercepts and
  slopes in terms of separate mean intercept and slope for each
  gender.

  The G and GCORR options in the RANDOM statement ask that the 
  D matrix and its corresponding correlation matrix it implies
  be printed.  The V and VCORR options ask that the overall
  V matrix be printed (for the first subject or particular
  subjects as specified).

  To fit a random coefficient model, we must specify that both
  intercept and slope are random in the RANDOM statement.

*******************************************************************/

/******************************************************************

  MODEL (a)
    
  R_i = diagonal with constant variance sigma^2 same in both genders
  No REPEATED statement necessary to fit this R_i (default)
      
  D = (2x2) unstructured matrix same for both genders, specified
  in the RANDOM statement

*******************************************************************/
    
title '(a) DIAGONAL WITHIN-CHILD COVARIANCE MATRIX R_i';
title2 'WITH CONSTANT VARIANCE SAME FOR EACH GENDER';
title3 'SAME D MATRIX FOR BOTH GENDERS';
proc mixed method=ml data=dent1;
  class gender child;
  model distance = gender gender*age / noint solution;
  random intercept age / type=un subject=child g gcorr v vcorr;

* comparison of mean slopes and whether mean intercepts and slopes;
* coincide for each gender;

  estimate 'diff in mean slope' gender 0 0 gender*age 1 -1;
  contrast 'overall gender diff' gender 1 -1, gender*age 1 -1 /chisq;
run;

/******************************************************************

   MODEL (b)

    Fit the same model as (a) but with a separate diagonal Ri matrix for
   each gender.  Thus, there are 2 separate variances sigma^2_(G and B)
   specified using GROUP=GENDER in the REPEATED statement

   D still = (2x2) unstructured matrix same for both genders as in (a)
   (specified in the RANDOM statement)

*******************************************************************/
    
title '(b) DIAGONAL WITHIN-CHILD COVARIANCE MATRIX R_i';
title2 'WITH SEPARATE CONSTANT VARIANCE FOR EACH GENDER';
title3 'SAME D MATRIX FOR BOTH GENDERS';
proc mixed  method=ml data=dent1; 
  class  child gender;
  model distance = gender gender*age / noint solution;
  repeated / group=gender subject=child;
  random intercept age / type=un subject=child g gcorr v=1,12 vcorr=1,12;
  estimate 'diff in mean slope' gender 0 0 gender*age 1 -1;
  contrast 'overall gender diff' gender 1 -1, gender*age 1 -1 /chisq;
run;

/******************************************************************

   MODEL (c)

   Fit the same model as (b) but now with D a (2x2) unstructured
   matrix different across genders, so there are two matrices D_G and D_B
   This is specified in the RANDOM statement by the GROUP=GENDER option

*******************************************************************/
    
title '(c) DIAGONAL WITHIN-CHILD COVARIANCE MATRIX R_i';
title2 'WITH SEPARATE CONSTANT VARIANCE FOR EACH GENDER';
title3 'DIFFERENT D MATRIX FOR BOTH GENDERS';
proc mixed  method=ml data=dent1; 
  class  child gender;
  model distance = gender gender*age / noint solution;
  repeated / group=gender subject=child;
  random intercept age / type=un group=gender subject=child g gcorr v=1,12 vcorr=1,12;
  estimate 'diff in mean slope' gender 0 0 gender*age 1 -1;
  contrast 'overall gender diff' gender 1 -1, gender*age 1 -1 /chisq;
run;
/******************************************************************

   MODEL (d)

   R_ii is the sum of two components, an AR(1) component for realizations
   with common variance sigma^2_P for both genders and a diagonal
   component with variance sigma^2_M common to both genders for 
   measurerment error    

   In the REPEATED statement, TYPE = specifies the realization component

   Here, we use the AR(1) specification, as we have equally-spaced
   time points.  SAS allows the correlation parameter alpha to satisfy 
   | alpha | <= 1, so that negative values are possible.  In the "Covariance
   Parameter Estimates" table, "Variance" is sigma^2_P and "AR(1)" is alpha.
    
   The LOCAL option specifies the measurement error component; in the table,
   "Residual" is sigma^2_M.

   Thus, the estimated R_i matrix has diagonal elements sigma^2_P+sigma^_M
   = "Variance + Residual"  and off-diagonal elements sigma^2_P alpha^|j-j'|,
   where (j,j') are time indices of time.  Thus, the implied correlation
   between Y_ij and Y_ij' is alpha^|j-j'| sigma^2_P/(sigma^2_P+sigma^2_M).
    
   D = (2x2) unstructured matrix same for both genders is specified as 
   above in the RANDOM statement

   Even though from the .log file the optimization seems to have
   "converged,"  note that things are suspect.  The estimated lag-1 correlation 
   parameter alpha = -0.9935, or virtually -1.  Clearly, the optimization has
   been driven to put the estimate on the boundary of the parameter
   space, reflecting that there is some practical difficulty in identifying
   this model.  This is most likely because there really isn't any
   need for the AR(1) structure.
    
*******************************************************************/
    
title '(d) COMMON AR(1) WITHIN-CHILD REALIZATION COMPONENT AND';
title2 'COMMON WITHIN-CHILD MEASUREMENT ERROR COMPONENT FOR R_i';
title3 'SAME D MATRIX FOR BOTH GENDERS';
proc mixed method=ml data=dent1;
  class gender child time;
  model distance = gender gender*age / noint solution ;
  random intercept age / type=un subject=child g gcorr v vcorr;
  repeated time / type=ar(1) local subject=child rcorr r;
  estimate 'diff in mean slope' gender 0 0 gender*age 1 -1;
  contrast 'overall gender diff' gender 1 -1, gender*age 1 -1 /chisq;
run;


/******************************************************************

   MODEL (e)
    
   Re-fit Model (d) using the SP(EXP)(time) correlation model, which
   is just a reparameterization of the AR(1) model in our case, where
   where alpha in the AR(1) parameterization is replaced by exp(-1/alpha^*).
   Thus, this specification enforces the correlation to be > 0 -- if
   alpha^* starts heading off to zero, the estimated implied correlation
   is zero.  
    
   This fit results in a message in the .log file:

   NOTE: Convergence criteria met but final hessian is not positive definite.
   NOTE: A linear combination of covariance parameters is confounded with the 
         residual
    
   This is because this parameterization constrains the within-child
   correlation parameter alpha to be positive as noted above, and it is
   being driven to zero.  Thus, there is no way to identify both the 
   realization variance sigma^2_P ("Variance") and the measurement error
   variance sigma^2_M ("Residual").  Again, previous diagnostics suggested
   nonnegiligible within-child correlation, so it is moving the estimate
   of correlation parameter, and thus correlation, to 0.

 *******************************************************************/
    
title '(e) COMMON AR(1) WITHIN-CHILD REALIZATION COMPONENT AND';
title2 'COMMON WITHIN-CHILD MEASUREMENT ERROR COMPONENT FOR R_i';
title3 'USING EXPONENTIAL MODEL; SAME D MATRIX FOR BOTH GENDERS';
proc mixed method=ml data=dent1;
  class gender child time;
  model distance = gender gender*age / noint solution ;
  random intercept age / type=un subject=child g gcorr v vcorr;
  repeated  time / type=sp(exp)(age) local subject=child rcorr r;
  estimate 'diff in mean slope' gender 0 0 gender*age 1 -1;
  contrast 'overall gender diff' gender 1 -1, gender*age 1 -1 /chisq;
run;

/******************************************************************

   MODEL (f)
    
    This is the same as (d) except that we allow alpha and sigma^2_P
    to differ by gender -- apparently the LOCAL statement creates a 
    COMMON measurement error variance sigma^2_M ("Residual") despite
    the GROUP option.  See above for the description of the 
    parameterization and implied within-child correlation.
    
    This appears to converge, note that the implied overall R_i 
    correlation matrix for girls has estimated off-diagonals virtually  
    equal to zero.  The implied within-boy correlation, on the other
    hand, is -0.45 at lag 1.  Recall from the within-child examination
    of sample autocorrelation in Chapter 2 that this appeared strongly 
    negative, and we surmised that it is likely driven by an outlying
    observation for one boy.  
    
*******************************************************************/
    
title '(f) DIFFERENT AR(1) WITHIN-CHILD REALIZATION COMPONENT AND';
title2 'SAME WITHIN-CHILD MEASUREMENT ERROR COMPONENT BY GENDER';
title3 'DIFFERENT D MATRIX FOR EACH GENDER';
proc mixed method=ml data=dent1;
  class gender child time;
  model distance = gender gender*age / noint solution ;
  random intercept age / type=un subject=child g gcorr v=1,12 vcorr=1,12 group=gender;
  repeated time / type=ar(1) local subject=child r = 1,12 rcorr=1,12 group=gender;
  estimate 'diff in mean slope' gender 0 0 gender*age 1 -1;
  contrast 'overall gender diff' gender 1 -1, gender*age 1 -1 /chisq;
run;

/******************************************************************

   Here are the AIC and BIC values for each fit:

                 (a)  AIC (Smaller is Better)         443.8
                 (b)  AIC (Smaller is Better)         424.0*
                 (c)  AIC (Smaller is Better)         429.1
                 (d)  AIC (Smaller is Better)         444.0
                 (e)  AIC (Smaller is Better)         445.8
                 (f)  AIC (Smaller is Better)         432.6

                 (a)  BIC (Smaller is Better)         454.2
                 (b)  BIC (Smaller is Better)         435.7*
                 (c)  BIC (Smaller is Better)         444.7
                 (d)  BIC (Smaller is Better)         457.0
                 (e)  BIC (Smaller is Better)         457.5
                 (f)  BIC (Smaller is Better)         452.0

    From above, the AIC and BIC values for model (d) are meaningless
    because the optimization did not really "converge."

    Both criteria suggest, among the models that can be identified
    from the data, model (b), with common D and diagonal R_i with
    separate within-child variances for each gender fits best.  (These
    variances thus represent the sum of realization process and 
    measurement error variance for each gender.)  
    
    We can compare the results to the PA fit we did earlier.  The
    AIC and BIC using ML for the preferred PA model with overall
    covariance structure compound symmetric and different by gender  
    yielded AIC = 424.8 and BIC = 435.2.  These values are entirely
    comparable to those for model (b).  In fact, the fitted overall
    compound symmetric covariance and correlation matrix  from the PA
    fit has overall variance = 4.4704 and correlation = 0.8680  for
    girls and variance = 5.2041 and correlation = 0.4701 for boys.
    Here, the estimated V matrices are not exactly compound symmetric,
    but they are close: For girls
  
                       Estimated V Matrix for child 1
 
               Row        Col1        Col2        Col3        Col4

                 1      3.1426      2.7933      2.8889      2.9845
                 2      2.7933      3.4128      3.1426      3.3172
                 3      2.8889      3.1426      3.8411      3.6499
                 4      2.9845      3.3172      3.6499      4.4275

                      Estimated V Correlation Matrix for child 1
 
               Row        Col1        Col2        Col3        Col4

                 1      1.0000      0.8529      0.8315      0.8001
                 2      0.8529      1.0000      0.8680      0.8534
                 3      0.8315      0.8680      1.0000      0.8851
                 4      0.8001      0.8534      0.8851      1.0000

    and for boys

                      Estimated V Matrix for child 12
 
               Row        Col1        Col2        Col3        Col4

                 1      5.3271      2.7933      2.8889      2.9845
                 2      2.7933      5.5973      3.1426      3.3172
                 3      2.8889      3.1426      6.0256      3.6499
                 4      2.9845      3.3172      3.6499      6.6120

                      Estimated V Correlation Matrix for child 12
 
               Row        Col1        Col2        Col3        Col4

                 1      1.0000      0.5115      0.5099      0.5029
                 2      0.5115      1.0000      0.5411      0.5453
                 3      0.5099      0.5411      1.0000      0.5782
                 4      0.5029      0.5453      0.5782      1.0000

*******************************************************************/

/******************************************************************

   MODEL (b), revisited

   Make some diagnostic plots of the "residuals" -- there are two
   types of residuals:

  - marginal (PA)   Y_i - X_i betahat  for individual i

  - conditional (SS)   Y_i - X_i betahat - Z_i bhat_i  for individual i

  where bhat_i are the empirical Bayes estimates of the random effect b_i.
  The RESIDUAL option in the MODEL statement requests that the residuals
  be computed.  The PLOTS option in the PROC MIXED statement requests
  that plots of the residuals be constructed:  an ordinary residuals
  vs. predicted values, a Q-Q plot of the residuals to assess normality,
  and a histogram.  The OUTPM option along with the RESIDUAL option requests
  that the PA predicted values X_i betahat and the marginal residuals
  be output to a data set.  The OUTP option does the same for the SS
  predicted values X_i betahat + Z_i bhat_i and conditional residuals.

  Note that the conditional predicted values and residuals will be 
  subject to shrinkage, as discussed in Section 6.4.

  The SOLUTION option in the RANDOM statement produces the empirical  
  Bayes estimates bhat_i. 
    
  We use REML for this call.  The STUDENTPANEL graph name in the PLOTS
  option requests that the residuals be standardized by dividing by 
  estimates of their variances; see the PROC MIXED documentation.

  Can get estimated semivariogram of the residuals using PROC VARIOGRAM.  

  We request that the output, which includes the plots, be put in  
  the file dentresid.pdf.
    
*******************************************************************/
    
title '(b) DIAGONAL WITHIN-CHILD COVARIANCE MATRIX R_i';
title2 'WITH SEPARATE CONSTANT VARIANCE FOR EACH GENDER';
title3 'SAME D MATRIX FOR BOTH GENDERS';

ods graphics on;
ods pdf file="dentalresid.sas.pdf";
proc mixed  data=dent1 plots=studentpanel(marginal conditional) ; 
  class  child gender;
  model distance = gender gender*age / noint solution residual
      outp=dentss outpm=dentpa;
  repeated / group=gender subject=child;
  random intercept age / type=un subject=child g gcorr v=1,12
      vcorr=1,12 solution;
run;
ods pdf close;
ods graphics off;

*  Marginal predicted values and residuals for the first 5 girls;

proc print data=dentpa(obs=20);
run;

*  Conditional predicted values and residuals for the first 5 girls;

proc print data=dentss(obs=20);
run;

