/***********************************************************

   Homework 4, Problem 4

   Population-averaged analysis of the arthritis clinical 
   trial using SAS PROC GENMOD and PROC GEE
    
***********************************************************/

options ls=80 ps=59 nodate; run;

data arth; infile 'arthritis.dat';
  input id trt age month scale dscale;
  run;

/*****************************************************************

  Sample means at each week-drug combination

******************************************************************/

proc means mean std noprint data=arth;
   class trt month;
   var dscale;
   output out=arthsumm mean=mean std=sd;
run;

data arthsumm; set arthsumm;
  if trt=. then delete;
  if month=. then delete;
  drop _TYPE_ _FREQ_;
run;

title "sample proportions and SDs by week";
proc print data=arthsumm;
run;

/*****************************************************************

  Fit the logistic model using proc genmod and proc glimmix.
  The raw proportions suggest that the probability of severe
  arthritic symptoms trends downward slightly on placebo but
  appears to show a steeper, although more erratic, decline on
  auranofin.  So adopt a basic model that approximates
  this behavior, with a linear effect of time on the logit scale.
 
******************************************************************/

*  proc genmod with unstructured to look at likely structure;

title "Unstructured correlation, PROC GENMOD";
proc genmod data=arth descending;
  class id;
  model dscale = month trt*month / dist=binomial link=logit;
  repeated subject=id / type=un modelse corrw;
run;

*  proc genmod with compound symmetric working correlation matrix;

title "Compound symmetric correlation, PROC GENMOD";
proc genmod data=arth descending;
  class id;
  model dscale = month trt*month / dist=binomial link=logit;
  repeated subject=id / type=cs modelse corrw;
run;

*  proc glimmix with compound symmetric working correlation matrix;

title "Compound symmetric correlation, PROC GLIMMIX";
proc glimmix data=arth empirical method=mmpl;
  class id;
  model dscale = month trt*month / dist=binomial link=logit solution;
  random _residual_ / subject=id type=cs vcorr;
run;

* refit the compound symmetric models in alternative parameterization;

title "Compound symmetric correlation, PROC GENMOD";
proc genmod data=arth descending;
  class id trt;
  model dscale =  trt*month / dist=binomial link=logit;
  repeated subject=id / type=cs modelse corrw;
run;

*  proc glimmix with compound symmetric working correlation matrix;

title "Compound symmetric correlation, PROC GLIMMIX";
proc glimmix data=arth empirical method=mmpl;
  class id trt;
  model dscale = trt*month / dist=binomial link=logit solution;
  random _residual_ / subject=id type=cs vcorr;
run;

/*****************************************************************

  Now incorporate age in the model to examine dependence
  of probability at baseline and of decline of probability
  post-baseline.  Based on above, adopt compound symmetric as
  a reasonable model
    
******************************************************************/

title "Age at baseline, PROC GENMOD";
proc genmod data=arth descending;
  class id;
  model dscale = age month trt*month / dist=binomial link=logit;
  repeated subject=id / type=cs modelse corrw;
run;

title "Age at baseline, PROC GLIMMIX";
proc glimmix data=arth empirical method=mmpl;
  class id;
  model dscale = age month trt*month / dist=binomial
      link=logit solution;
  random _residual_ / subject=id type=cs vcorr;
run;

title "Age at baseline and in slope, PROC GENMOD";
proc genmod data=arth descending;
  class id;
  model dscale = age month month*age trt*month trt*month*age / dist=binomial link=logit;
  repeated subject=id / type=cs modelse corrw;
run;

title "Age at baseline and in slope, PROC GLIMMIX";
proc glimmix data=arth empirical method=mmpl;
  class id;
  model dscale = age month month*age trt*month trt*month*age / dist=binomial
      link=logit solution;
  random _residual_ / subject=id type=cs vcorr;
run;


