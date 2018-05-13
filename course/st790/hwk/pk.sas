* Fit the PK data using the nlinmix macro and proc nlmixed;

options ps=55 ls=80 nodate;

/*  include the macro from file nlinmix.sas */

%inc 'nlmm801.sas'  / nosource;

data pkdata; infile 'pk.dat';
    input id time conc weight age metab;
run;

data pkdata; set pkdata;
    meta1=0; meta2=0; meta3=0;
    if metab=1 then meta1=1;
    if metab=2 then meta2=1;
    if metab=3 then meta3=1;
run;

*  Model with no covariates;

title "Model with No Covariates, nlinmix";
%nlinmix(data=pkdata,
    model=%str(
     ka=exp(beta1+b1);
     cl=exp(beta2+b2);
     v=exp(beta3+b3);
     ke=cl/v; kk=ka-ke;
     dose=1000;
     predv=dose*ka*(exp(-ke*time)-exp(-ka*time))/(v*kk);
     ),
     derivs=%str(
     tt=dose*time*ka/(v*kk);
     d_beta1=-predv*(ke/kk)+tt*ka*exp(-ka*time);
     d_beta2=predv*(ke/kk)-tt*ke*exp(-ke*time);
     d_beta3=-predv*(ka/kk)+tt*ke*exp(-ke*time);
     d_b1=-predv*(ke/kk)+tt*ka*exp(-ka*time);
     d_b2=predv*(ke/kk)-tt*ke*exp(-ke*time);
     d_b3=-predv*(ka/kk)+tt*ke*exp(-ke*time);
    ), 
    weight=%str(
     _weight_= 1/predv**(2*1.06);
   ),
   parms=%str(beta1=-0.3 beta2=-0.2 beta3=1.6),
   stmts=%str(
      class id;
      model pseudo_conc = d_beta1 d_beta2 d_beta3 / noint notest solution ;
      random d_b1 d_b2 d_b3/ subject=id  type=un solution g gcorr;
      weight _weight_;
   ),
   expand=eblup,
   procopt=%str(maxiter=500 method=ml)
)
run;

*  model with all covariates;

title "Model with All Covariates, nlinmix"; 
   %nlinmix(data=pkdata,
    model=%str(
     ka=exp(beta1+beta4*meta2+beta5*meta3+beta6*age+beta7*weight+b1);
     cl=exp(beta2+beta8*meta2+beta9*meta3+beta10*age+beta11*weight+b2);
     v=exp(beta3+beta12*meta2+beta13*meta3+beta14*age+beta15*weight+b3);
     ke=cl/v; kk=ka-ke;
     dose=1000;
     predv=dose*ka*(exp(-ke*time)-exp(-ka*time))/(v*kk);
     ),
     derivs=%str(
     tt=dose*time*ka/(v*kk);
     d_beta1=-predv*(ke/kk)+tt*ka*exp(-ka*time);
     d_beta4=(-predv*(ke/kk)+tt*ka*exp(-ka*time))*meta2;
     d_beta5=(-predv*(ke/kk)+tt*ka*exp(-ka*time))*meta3;
     d_beta6=(-predv*(ke/kk)+tt*ka*exp(-ka*time))*age;
     d_beta7=(-predv*(ke/kk)+tt*ka*exp(-ka*time))*weight;
     
     d_beta2=predv*(ke/kk)-tt*ke*exp(-ke*time);
     d_beta8=(predv*(ke/kk)-tt*ke*exp(-ke*time))*meta2;
     d_beta9=(predv*(ke/kk)-tt*ke*exp(-ke*time))*meta3;
     d_beta10=(predv*(ke/kk)-tt*ke*exp(-ke*time))*age;
     d_beta11=(predv*(ke/kk)-tt*ke*exp(-ke*time))*weight;

     d_beta3=-predv*(ka/kk)+tt*ke*exp(-ke*time);
     d_beta12=(-predv*(ka/kk)+tt*ke*exp(-ke*time))*meta2;
     d_beta13=(-predv*(ka/kk)+tt*ke*exp(-ke*time))*meta3;
     d_beta14=(-predv*(ka/kk)+tt*ke*exp(-ke*time))*age;
     d_beta15=(-predv*(ka/kk)+tt*ke*exp(-ke*time))*weight;

     d_b1=-predv*(ke/kk)+tt*ka*exp(-ka*time);
     d_b2=predv*(ke/kk)-tt*ke*exp(-ke*time);
     d_b3=-predv*(ka/kk)+tt*ke*exp(-ke*time);
    ), 
    weight=%str(
     _weight_= 1/predv**(2*1.06);
   ),
   parms=%str(beta1=-0.3 beta2=-0.2 beta3=1.6 beta4=0 beta5=0 beta6=0
              beta7=0 beta8=0 beta9=0 beta10=0 beta11=0 beta12=0
              beta13=0 beta14=0 beta15=0),
   stmts=%str(
      class id;
   model pseudo_conc = d_beta1 d_beta2 d_beta3 d_beta4 d_beta5 d_beta6
       d_beta7 d_beta8 d_beta9 d_beta10 d_beta11 d_beta12 d_beta13
       d_beta14 d_beta15 / noint notest solution ;
      random d_b1 d_b2 d_b3/ subject=id  type=un solution g gcorr;
      weight _weight_;
   ),
   expand=eblup,
   procopt=%str(maxiter=500 method=ml)
)
   
*  model with some covariates;

title "Model with Apparent Important Covariates, nlinmix";
   %nlinmix(data=pkdata,
    model=%str(
     ka=exp(beta1+b1);
     cl=exp(beta2+beta4*meta2+beta5*meta3+beta6*age+b2);
     v=exp(beta3+beta7*weight+b3);
     ke=cl/v; kk=ka-ke;
     dose=1000;
     predv=dose*ka*(exp(-ke*time)-exp(-ka*time))/(v*kk);
     ),
     derivs=%str(
     tt=dose*time*ka/(v*kk);
     d_beta1=-predv*(ke/kk)+tt*ka*exp(-ka*time);
          
     d_beta2=predv*(ke/kk)-tt*ke*exp(-ke*time);
     d_beta4=(predv*(ke/kk)-tt*ke*exp(-ke*time))*meta2;
     d_beta5=(predv*(ke/kk)-tt*ke*exp(-ke*time))*meta3;
     d_beta6=(predv*(ke/kk)-tt*ke*exp(-ke*time))*age;

     d_beta3=-predv*(ka/kk)+tt*ke*exp(-ke*time);
     d_beta7=(-predv*(ka/kk)+tt*ke*exp(-ke*time))*weight;

     d_b1=-predv*(ke/kk)+tt*ka*exp(-ka*time);
     d_b2=predv*(ke/kk)-tt*ke*exp(-ke*time);
     d_b3=-predv*(ka/kk)+tt*ke*exp(-ke*time);
    ), 
    weight=%str(
     _weight_= 1/predv**(2*1.06);
   ),
   parms=%str(beta1=-0.3 beta2=-0.2 beta3=1.6 beta4=0 beta5=0 beta6=0 beta7=0),
   stmts=%str(
      class id;
   model pseudo_conc = d_beta1 d_beta2 d_beta3 d_beta4 d_beta5 d_beta6 d_beta7 / noint notest solution ;
      random d_b1 d_b2 d_b3/ subject=id  type=un solution g gcorr;
      weight _weight_;
    ),
   expand=eblup,
   procopt=%str(maxiter=500 method=ml)
)
       
*  Model with no covariates, full Laplace;

title "Model with No Covariates, proc nlmixed, Laplace";   
proc nlmixed data=pkdata method=gauss qpoints=1;
    parms beta1=0.2 beta2=-0.1 beta3=1.6 s2b1=0.01 s2b2=0.3 s2b3=0.02
        cb12=0.02 cb13=0.009 cb23=0.026 s2=0.01 delta=1;
     ka=exp(beta1+b1);
     cl=exp(beta2+b2);
     v=exp(beta3+b3);
     ke=cl/v; kk=ka-ke;
     dose=1000;
     pred=dose*ka*(exp(-ke*time)-exp(-ka*time))/(v*kk);
     model conc ~ normal(pred,(pred**(2*delta))*s2);
     random b1 b2 b3 ~ normal([0,0,0],[s2b1,cb12,s2b2,cb13,cb23,s2b3])
         subject=id;
run;

*  Model with some covariates, full Laplace;

title "Model with Apparent Important Covariates, proc nlmixed, Laplace";
proc nlmixed data=pkdata method=gauss qpoints=1;
    parms beta1=0.2 beta2=-0.66 beta3=1.97 beta4=0.9 beta5=1.4 beta6=-0.01
        beta7=-0.005 s2b1=0.011 s2b2=0.012 s2b3=0.013
        cb12=0.004 cb13=0.007 cb23=0.009 s2=0.01 delta=1;
     ka=exp(beta1+b1);
     cl=exp(beta2+beta4*meta2+beta5*meta3+beta6*age+b2);
     v=exp(beta3+beta7*weight+b3);
     ke=cl/v; kk=ka-ke;
     dose=1000;
     pred=dose*ka*(exp(-ke*time)-exp(-ka*time))/(v*kk);
     model conc ~ normal(pred,(pred**(2*delta))*s2);
     random b1 b2 b3 ~ normal([0,0,0],[s2b1,cb12,s2b2,cb13,cb23,s2b3])
         subject=id;
run;

/*
*  Model with some covariates, quadrature;

*  I couldn't get this to converge; probably need better starting values;

proc nlmixed data=pkdata method=gauss qpoints=10;
    parms beta1=0.2 beta2=-0.66 beta3=1.97 beta4=0.9 beta5=1.4 beta6=-0.01
        beta7=-0.005 s2b1=0.011 s2b2=0.012 s2b3=0.013
        cb12=0.004 cb13=0.007 cb23=0.009 s2=0.01 delta=1;
     ka=exp(beta1+b1);
     cl=exp(beta2+beta4*meta2+beta5*meta3+beta6*age+b2);
     v=exp(beta3+beta7*weight+b3);
     ke=cl/v; kk=ka-ke;
     dose=1000;
     pred=dose*ka*(exp(-ke*time)-exp(-ka*time))/(v*kk);
     model conc ~ normal(pred,(pred**(2*delta))*s2);
     random b1 b2 b3 ~ normal([0,0,0],[s2b1,cb12,s2b2,cb13,cb23,s2b3])
         subject=id;
run;
*/
