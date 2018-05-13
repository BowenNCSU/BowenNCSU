/******************************************************************

  Hmework 5, Problem 4
    
  Fit a subject specific model to the arthritis data
  with "intercept" random effect only
  
******************************************************************/

options ls=80 ps=59 nodate; run;

data arth; infile 'arthritis.dat';
  input id trt age month scale dscale;
  run;

/******************************************************************

 PA model with compound symmetric covariance structure
  
******************************************************************/

title "PA MODEL USING GLIMMIX";
proc glimmix data=arth empirical method=mmpl;
  class id;
  model dscale = month trt*month / dist=binomial link=logit solution;
  random _residual_ / subject=id type=cs vcorr;
run;

title "PA MODEL USING GENMOD";
proc genmod data=arth descending;
  class id;
  model dscale = month trt*month / dist=binomial link=logit;
  repeated subject=id / type=cs modelse corrw;
run;

/*****************************************************************

  Fit the logistic model using proc glimmix, proc nlmixed, and the
  glimmix macro
 
******************************************************************/

title "LINEARIZATION ABOUT 0, GLIMMIX";
proc glimmix data=arth method=mmpl;
  class id;
  model dscale = month trt*month / dist=binomial link=logit solution;
  random intercept / subject=id;
run;

title "LINEARIZATION ABOUT EMP BAYES ESTS, GLIMMIX";
proc glimmix data=arth method=mspl;
  class id;
  model dscale = month trt*month / dist=binomial link=logit solution;
  random intercept / subject=id;
run;

title "FULL LAPLACE APPROX, GLIMMIX";
proc glimmix data=arth method=laplace;
  class id;
  model dscale = month trt*month / dist=binomial link=logit solution;
  random intercept / subject=id;
run;

title "QUADRATURE 10, GLIMMIX";
proc glimmix data=arth method=quad(qpoints=10);
  class id;
  model dscale = month trt*month / dist=binomial link=logit solution;
  random intercept / subject=id;
run;

title "QUADRATURE 20, GLIMMIX";
proc glimmix data=arth method=quad(qpoints=20);
  class id;
  model dscale = month trt*month / dist=binomial link=logit solution;
  random intercept / subject=id;
run;

title "QUADRATURE 50, GLIMMIX";
proc glimmix data=arth method=quad(qpoints=50);
  class id;
  model dscale = month trt*month / dist=binomial link=logit solution;
  random intercept / subject=id;
run;

title "FULL LAPLACE APPROX, NLMIXED";
proc nlmixed data=arth method=gauss qpoints=1;
  parms b1=-0.8 b2=-0.06 b3=-0.1 su2=3.5;
  num = exp(b1+b2*month+b3*month*trt+u);
  prob = num/(1+num);
  model dscale ~ binary(prob);
  random u ~ normal(0,su2) subject=id;
run;

title "QUADRATURE 10, NLMIXED";
proc nlmixed data=arth method=gauss qpoints=10;
  parms b1=-0.8 b2=-0.06 b3=-0.1 su2=3.5;
  num = exp(b1+b2*month+b3*month*trt+u);
  prob = num/(1+num);
  model dscale ~ binary(prob);
  random u ~ normal(0,su2) subject=id;
run;

title "QUADRATURE 20, NLMIXED";
proc nlmixed data=arth method=gauss qpoints=20;
  parms b1=-0.8 b2=-0.06 b3=-0.1 su2=3.5; 
  num = exp(b1+b2*month+b3*month*trt+u);
  prob = num/(1+num);
  model dscale ~ binary(prob);
  random u ~ normal(0,su2) subject=id;
run;

title "QUADRATURE 50, NLMIXED";
proc nlmixed data=arth method=gauss qpoints=50;
  parms b1=-0.8 b2=-0.06 b3=-0.1 su2=3.5; 
  num = exp(b1+b2*month+b3*month*trt+u);
  prob = num/(1+num);
  model dscale ~ binary(prob);
  random u ~ normal(0,su2) subject=id;
run;

%inc 'glmm800.sas'  / nosource;

title "LINEARIZATION USING GLIMMIX MACRO";
%glimmix(data=arth,
    procopt=method=ml,
    stmts=%str(
     class id;
     model dscale = month trt*month  / solution;
     random intercept / subject=id type=un;
  ), 
  error=binomial,
  link=logit
)

