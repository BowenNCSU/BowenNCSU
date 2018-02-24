/*******************************************************************

  Weightloss data

  Population averaged model analysis using proc mixed

*******************************************************************/

options ls=80 ps=59 nodate; run;

/******************************************************************

  The data set is in "wide" format, so must convert to "long" format,
  which is what proc mixed wants
   
*******************************************************************/

data weight1; infile "weightloss.dat";
    input id month0 month3 month6 month9 month12 diet;
run;

data weight2; set weight1;
    array wt(5) month0 month3 month6 month9 month12;
    do m = 1 to 5;
        weight = wt(m);
        output;
        end;
    drop month0 month3 month6 month9 month12;
run;

data weight2; set weight2;
    if m=1 then month=0;
    if m=2 then month=3;
    if m=3 then month=6;
    if m=4 then month=9;
    if m=5 then month=12;
    drop m;
run;

data weight2; set weight2;
    month2=month*month;
    time=month;
run;
    
/******************************************************************

  We assume a basic mean model with separate intercepts, linear,
  and quadratic terms for each agent.  You may have decided to
  collapse this model to have a common intercept after investigation
  of whether the intercepts are the same across treatments.  We don't
  do that in the analyses here, but you might have.
    
  Fit various covariance models.  From the crude empirical evidence, 
  compound symmetry seems plausible except possibly for the control group.
  It also seems like variances don't change very much over time.  In
  proc mixed, we can fit unstructured covariance that is either same or
  different in each group, but we cannot force the variances to be the same
  over time.  We fit unstructured models that are different for each group and
  then the same for each group.  We can also fit compound symmetric models that
  are different or the same for each group, with variances that change over
  time (type=csh) or don't (type=cs). We fit all of these using ML.
  The AIC values from the fits below, in order of appearance, are

                     AIC (Smaller is Better)        4154.9
                     AIC (Smaller is Better)        4124.4
                     AIC (Smaller is Better)        4146.5
                     AIC (Smaller is Better)        4136.1
                     AIC (Smaller is Better)        4140.1
                     AIC (Smaller is Better)        4135.7       

    and the BIC values are

                     BIC (Smaller is Better)        4295.6
                     BIC (Smaller is Better)        4186.9
                     BIC (Smaller is Better)        4216.9
                     BIC (Smaller is Better)        4175.1
                     BIC (Smaller is Better)        4179.2
                     BIC (Smaller is Better)        4164.3

    Based on both AIC and BIC, the model with common unstructured overall 
    covariance structure is preferred. The sample correlation matrices for
    groups 2 and 3 seemed approximately compound symmetric, while that for
    group 1 did not, so this fit is probably a compromose.  Using gls() in R,
    a common unstructured model that further restricts the variances to be
    the same over time, which seemed supported by the sample variances in each 
    in each group, is preferred.  Given we can't fit that model in proc mixed,
    we settle on the common unstructured model.  T
    
    If one prefers to switch to R and gls() at this point, one can, using either 
    this model or the one with constant variance over time, which can't be fit here.
    
*******************************************************************/

title "UNSTRUCTURED SEPARATE BY PROGRAM"; 
proc mixed data=weight2 method=ML; 
  class diet id time; 
  model weight = diet diet*month diet*month2 / noint solution ; 
  repeated time / type = un subject = id r=1,35,63 rcorr=1,35,63 group=diet; 
run;

title "COMMON UNSTRUCTURED";
proc mixed data=weight2 method=ML;
  class diet id time;
  model weight = diet diet*month diet*month2 / noint solution ;
  repeated time / type = un subject = id r rcorr;
run;

title "HETEROGENEOUS COMPOUND SYMMETRIC SEPARATE BY PROGRAM";
proc mixed data=weight2 method=ML;
  class diet id time;
  model weight = diet diet*month diet*month2 / noint solution ;
  repeated time / type = csh subject = id r=1,35,63 rcorr=1,35,63 group=diet;
run;

title "COMMON HETEROGENEOUS COMPOUND SYMMETRIC";
proc mixed data=weight2 method=ML;
  class diet id time;
  model weight = diet diet*month diet*month2 / noint solution ;
  repeated time / type = csh subject = id r rcorr;
run;

title "COMPOUND SYMMETRIC SEPARATE BY PROGRAM";
proc mixed data=weight2 method=ML;
  class diet id time;
  model weight = diet diet*month diet*month2 / noint solution ;
  repeated time / type = cs subject = id r rcorr group=diet;
run;

title "COMMON COMPOUND SYMMETRIC";
proc mixed data=weight2 method=ML;
  class diet id time;
  model weight = diet diet*month diet*month2 / noint solution ;
  repeated time / type = cs subject = id r rcorr;
run;

/******************************************************************

  Adopting the common UN structure as above, we fit again. We don't request 
  robust sandwich standard errors given the structure is unstructured.
  Again, we leave in separate intercepts for each treatment; you may have
  simplified the model further to have a common intercept.
    
*******************************************************************/

title "ANALYSES WITH COMMON UNSTRUCTURED COVARIANCE";
proc mixed data=weight2 method=ML;
  class diet id time;
  model weight = diet diet*month diet*month2 / noint solution ;
  repeated time / type = un subject = id;

*  get the test of quadratic effects manually;

  contrast 'quadratic effect' diet*month2 1 0 0,
                              diet*month2 0 1 0,
                              diet*month2 0 0 1 / chisq;

  contrast 'diff quad effect' diet*month2 1 -1 0,
      diet*month2 1 0 -1 /chisq;

  contrast 'diff means @ month 12'
      diet 1 -1 0 diet*month 12 -12 0 diet*month2 144 -144 0,
      diet 1 0 -1 diet*month 12 0 -12 diet*month2 144 0 -144 /chisq;
           
  estimate 'mean @ month 12, diet 1' diet 1 0 0  diet*month 12 0 0 diet*month2
      144 0 0 ;
  estimate 'mean @ month 12, diet 2' diet 0 1 0  diet*month 0 12 0 diet*month2
      0 144 0;
  estimate 'mean @ month 12, diet 3' diet 0 0 1 diet*month 0 0 12 diet*month2
      0 0 144;

  estimate 'rate of chng @ month 8 diet 1' diet*month 1 0 0 diet*month2 16 0 0 ;
  estimate 'rate of chng @ month 8 diet 2' diet*month 0 1 0 diet*month2 0 16 0 ;
  estimate 'rate of chng @ month 8 diet 3' diet*month 0 0 1 diet*month2 0 0 16 ;
run;

*  with alternative parametrization to get tests if quadratic;
*  and linear terms differ across groups;

proc mixed data=weight2 method=ML;
  class diet id time;
  model weight = diet month diet*month month2 diet*month2 / solution chisq covb;
  repeated time / type = un subject = id;
run;
 
