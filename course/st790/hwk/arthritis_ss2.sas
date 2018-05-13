/******************************************************************

  Hmework 5, Problem 4
    
  Fit a subject specific model to the arthritis data
  with "intercept" and "slope" random effects
  
******************************************************************/

options ls=80 ps=59 nodate; run;

data arth; infile 'arthritis.dat';
  input id trt age month scale dscale;
  run;

/*****************************************************************

  Fit the logistic model using proc glimmix, proc nlmixed, and the
  glimmix macro
 
******************************************************************/

title "LINEARIZATION ABOUT 0, GLIMMIX";
proc glimmix data=arth method=mmpl;
  class id;
  model dscale = month trt*month / dist=binomial link=logit solution;
  random intercept month / subject=id type=un g gcorr;
run;

title "LINEARIZATION ABOUT EMP BAYES ESTS, GLIMMIX";
proc glimmix data=arth method=mspl;
  class id;
  model dscale = month trt*month / dist=binomial link=logit solution;
  random intercept month / subject=id type=un g gcorr;
run;

title "FULL LAPLACE APPROX, GLIMMIX";
proc glimmix data=arth method=laplace;
  nloptions gtol=1E-6;    *  relax the convergence criterion;
  class id;
  model dscale = month trt*month / dist=binomial link=logit solution;
  random intercept month / subject=id type=un g gcorr;
run;

title "QUADRATURE 10, GLIMMIX";
proc glimmix data=arth method=quad(qpoints=10);
  class id;
  model dscale = month trt*month / dist=binomial link=logit solution;
  random intercept month / subject=id type=un g gcorr;
run;

title "QUADRATURE 20, GLIMMIX";
proc glimmix data=arth method=quad(qpoints=20);
  class id;
  model dscale = month trt*month / dist=binomial link=logit solution;
  random intercept month / subject=id type=un g gcorr;
run;

title "QUADRATURE 50, GLIMMIX";
proc glimmix data=arth method=quad(qpoints=50);
  class id;
  model dscale = month trt*month / dist=binomial link=logit solution;
  random intercept month / subject=id type=un g gcorr;
run;

title "FULL LAPLACE APPROX, NLMIXED";
proc nlmixed data=arth method=gauss qpoints=1;
  parms b1=-0.8 b2=-0.06 b3=-0.1 su12=3 c12=0 su22=0.1;
  num = exp(b1+b2*month+b3*month*trt+u1+month*u2);
  prob = num/(1+num);
  model dscale ~ binary(prob);
  random u1 u2 ~ normal([0,0],[su12,c12,su22]) subject=id;
run;

title "QUADRATURE 10, NLMIXED";
proc nlmixed data=arth method=gauss qpoints=10;
  parms b1=-0.8 b2=-0.06 b3=-0.1 su12=3 c12=0 su22=0.1;
  num = exp(b1+b2*month+b3*month*trt+u1+month*u2);
  prob = num/(1+num);
  model dscale ~ binary(prob);
  random u1 u2 ~ normal([0,0],[su12,c12,su22]) subject=id;
run;

title "QUADRATURE 20, NLMIXED";
proc nlmixed data=arth method=gauss qpoints=20;
   parms b1=-0.8 b2=-0.06 b3=-0.1 su12=3 c12=0 su22=0.1;
  num = exp(b1+b2*month+b3*month*trt+u1+month*u2);
  prob = num/(1+num);
  model dscale ~ binary(prob);
  random u1 u2 ~ normal([0,0],[su12,c12,su22]) subject=id;
run;

title "QUADRATURE 50, NLMIXED";
proc nlmixed data=arth method=gauss qpoints=50;
   parms b1=-0.8 b2=-0.06 b3=-0.1 su12=3 c12=0 su22=0.1;
  num = exp(b1+b2*month+b3*month*trt+u1+month*u2);
  prob = num/(1+num);
  model dscale ~ binary(prob);
  random u1 u2 ~ normal([0,0],[su12,c12,su22]) subject=id;
run;

%inc 'glmm800.sas'  / nosource;

title "LINEARIZATION USING GLIMMIX MACRO";
%glimmix(data=arth,
    procopt=method=ml,
    stmts=%str(
     class id;
     model dscale = month trt*month  / solution;
     random intercept month / subject=id type=un;
  ), 
  error=binomial,
  link=logit
)


