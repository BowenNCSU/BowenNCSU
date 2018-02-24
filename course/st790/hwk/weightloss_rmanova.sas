/*******************************************************************

  Weightloss Data 

  Univariate and multivariate repeated measures analysis of variance   

*******************************************************************/

options ls=80 ps=59 nodate; run;

/******************************************************************

  The data set is in "wide" format.  Read in the data and create an
  alternative data set in "long" form with one record per observation.

*******************************************************************/

data weight1; infile "weightloss.dat";
    input id month0 month3 month6 month9 month12 program;
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

proc print data=weight2(obs=10); run;

/******************************************************************

  Get the univariate ANOVA different ways

*******************************************************************/

title "UNIVARIATE ANALYSIS VIA SPLIT PLOT SPECIFICATION USING PROC GLM";
proc glm data=weight2;
  class month program id;
  model weight = program id(program) month month*program;
  random id(program) / test;
run;

title "UNIVARIATE ANALYSIS USING PROC MIXED";
proc mixed data=weight2 method=type3;
  class month program id;
  model weight = program month month*program;
  random id(program);
run;

title "UNIVARIATE ANALYSIS USING REPEATED STATEMENT IN PROC GLM";
proc glm data=weight1;
  class program;
  model  month0 month3 month6 month9 month12 = program / nouni;
  repeated month / printe nom;
run;

title "ORTHOGONAL POLYNOMIAL TRANSFORMATION";
proc glm data=weight1;
  class program;
  model month0 month3 month6 month9 month12 = program / nouni;
  repeated month 5 (0 3 6 9 12) polynomial /summary printm nom;
run;

title "UNNORMALIZED PROFILE TRANSFORMATION";
proc glm data=weight1;
  class program; 
  model month0 month3 month6 month9 month12 = program / nouni;
  repeated month 5 (0 3 6 9 12) profile /summary printm nom;
run;

title "UNNORMALIZED HELMERT TRANSFORMATION";
proc glm data=weight1;
  class program; 
  model month0 month3 month6 month9 month12 = program / nouni;
  repeated month 5 (0 3 6 9 12) helmert /summary printm nom;
run;

/******************************************************************

  Get the multiariate ANOVA 

*******************************************************************/

title "MULTIVARIATE ANALYSIS USING PROC GLM MANOVA STATEMENT";
proc glm data=weight1;
  class program;
  model month0 month3 month6 month9 month12 = program / nouni;
  means program;
  manova h=program / printh printe;
run;

title "MULTIVARIATE PROFILE ANALYSIS";
proc glm data=weight1;
  class program;
  model month0 month3 month6 month9 month12 = program / nouni;
  repeated month / printe nou;
run;




