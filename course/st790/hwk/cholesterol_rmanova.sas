/*******************************************************************

  Weightloss Data 

  Univariate and multivariate repeated measures analysis of variance   

*******************************************************************/

options ls=80 ps=59 nodate; run;

/******************************************************************

  The data set is in "wide" format.  Read in the data and create an
  alternative data set in "long" form with one record per observation.

*******************************************************************/

data chol1; infile "cholesterol.dat";
    input trt id month0 month6 month12 month20 month24;
run;

data chol2; set chol1;
    array wt(5) month0 month6 month12 month20 month24;
    do m = 1 to 5;
        chol = wt(m);
        output;
        end;
    drop month0 month6 month12 month20 month24;
run;

data chol2; set chol2;
    if m=1 then month=0;
    if m=2 then month=6;
    if m=3 then month=12;
    if m=4 then month=20;
    if m=5 then month=24;
    drop m;
run;

proc print data=chol2(obs=10); run;

/******************************************************************

  Get the univariate ANOVA different ways

*******************************************************************/

title "UNIVARIATE ANALYSIS VIA SPLIT PLOT SPECIFICATION USING PROC GLM";
proc glm data=chol2;
  class month trt id;
  model chol = trt id(trt) month month*trt;
  random id(trt) / test;
run;

title "UNIVARIATE ANALYSIS USING PROC MIXED";
proc mixed data=chol2 method=type3;
  class month trt id;
  model chol = trt month month*trt;
  random id(trt);
run;

title "UNIVARIATE ANALYSIS USING REPEATED STATEMENT IN PROC GLM";
proc glm data=chol1;
  class trt;
  model month0 month6 month12 month20 month24 = trt / nouni;
  repeated month / printe nom;
run;

title "ORTHOGONAL POLYNOMIAL TRANSFORMATION";
proc glm data=chol1;
  class trt;
  model month0 month6 month12 month20 month24 = trt / nouni;
  repeated month 5 (0 6 12 20 24) polynomial /summary printm nom;
run;

title "UNNORMALIZED PROFILE TRANSFORMATION";
proc glm data=chol1;
  class trt; 
  model month0 month6 month12 month20 month24 = trt / nouni;
  repeated month 5 (0 6 12 20 24) profile /summary printm nom;
run;

title "UNNORMALIZED HELMERT TRANSFORMATION";
proc glm data=chol1;
  class trt; 
  model month0 month6 month12 month20 month24 = trt / nouni;
  repeated month 5 (0 6 12 20 24) helmert /summary printm nom;
run;

/******************************************************************

  Get the multiariate ANOVA 

*******************************************************************/

title "MULTIVARIATE ANALYSIS USING PROC GLM MANOVA STATEMENT";
proc glm data=chol1;
  class trt;
  model month0 month6 month12 month20 month24 = trt / nouni;
  means trt;
  manova h=trt / printh printe;
run;

title "MULTIVARIATE PROFILE ANALYSIS";
proc glm data=chol1;
  class trt;
  model month0 month6 month12 month20 month24 = trt / nouni;
  repeated month / printe nou;
run;




