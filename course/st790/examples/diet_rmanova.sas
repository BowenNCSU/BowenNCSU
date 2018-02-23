/*******************************************************************

  CHAPTER 3, EXAMPLE 2, Guinea Pig Vitamin E Diet Study

  Univariate and multivariate repeated measures analysis of variance   

*******************************************************************/

options ls=80 ps=59 nodate; run;

/******************************************************************

  The data set is in "wide" format with first 5 lines 

1	455	460	510	504	436	466  1
2	467	565	610	596	542	587  1
3	445	530	580	597	582	619  1
4	485	542	594	583	611	612  1
5	480	500	550	528	562	576  1

  column 1      pig number
  columns 2-7   body weights at weeks 1, 3, 4, 5, 6, 7
  column 8      dose group  (1=zero, 2 = low, 3 = high dose

  Read in the data and create an alternative data set in "long"
  form with one record per observation.

*******************************************************************/

data pigs1; infile 'diet.dat';
  input pig week1 week3 week4 week5 week6 week7 dose;

data pigs2; set pigs1;
  array wt(6) week1 week3 week4 week5 week6 week7;
  do week = 1 to 6;
     weight = wt(week);
     output;
  end;
  drop week1 week3-week7;
run;

data pigs2; set pigs2; 
  if week>1 then week=week+1;
run;

proc print data=pigs2(obs=6); run;

/*******************************************************************

  Construct the univariate analysis of variance using PROC GLM  
  via a "split plot" specification, which requires the data in
  "long" form (PIGS2).

  The F ratio that PROC GLM prints out automatically for the 
  dose effect (averaged across week) uses the MSE in the 
  denominator, so is incorrect.

  The RANDOM statement produces the expected mean squares

  The TEST option request the test for the dose effect
  treating the pig(dose) efffect as random, yielding the correct
  F ratio.  Other F-ratios are correct.
    
*******************************************************************/  

title "UNIVARIATE ANALYSIS VIA SPLIT PLOT SPECIFICATION USING PROC GLM";
proc glm data=pigs2;
  class week dose pig;
  model weight = dose pig(dose) week week*dose;
  random pig(dose) / test;
run;

/*******************************************************************

  This analysis can also be done using PROC MIXED.  The 
  MODEL and RANDOM statements have DIFFERENT interpretations in
  PROC GLM and PROC MIXED.  The METHOD=TYPE3 option produces the
  correct tests.

*******************************************************************/

title "UNIVARIATE ANALYSIS USING PROC MIXED";
proc mixed data=pigs2 method=type3;
  class week dose pig;
  model weight = dose week week*dose;
  random pig(dose);
run;

/*******************************************************************

  Now carry out the same analysis using the REPEATED statement in
  PROC GLM.  This requires that the data be represented in the 
  form of data set PIGS1.

  The option NOUNI suppresses individual analyses of variance
  at each week value from being printed.

  The PRINTE option asks for the test of sphericity to be performed.

  The NOM option means "no multivariate," which means univariate
  tests under the assumption that the compound symmetry model 
  is correct.
  
*******************************************************************/  

title "UNIVARIATE ANALYSIS USING REPEATED STATEMENT IN PROC GLM";
proc glm data=pigs1;
  class dose;
  model week1 week3 week4 week5 week6 week7 = dose / nouni;
  repeated week / printe nom;
run;

/*******************************************************************

  Using the REPEATED statement, one can obtain specialized within-
  individual tests; the REPEATED statement allows one to specify
  diffrent contrast matrices U.  We demonstrate with two calls to
  PROC GLM like the one above specifying the Helmert and polynomial
  transformations.  

  In both, the SUMMARY option asks that PROC GLM print out the results 
  of tests corresponding to the contrasts in each column of the U
  matrix.

  The PRINTM option prints the U matrix; SAS calls the transpose
  of this matrix M, so the M printed is our U'.

  The NOU option suppresses the results of the univariate analysis,
  which we already did above.

  The following analyses request the orthogonal polynomial
  transformation (which is normalized) and the profile and Helmert
  transformations (which are not).

*******************************************************************/  

title "ORTHOGONAL POLYNOMIAL TRANSFORMATION";
proc glm data=pigs1;
  class dose;
  model week1 week3 week4 week5 week6 week7 = dose / nouni;
  repeated week 6 (1 3 4 5 6 7) polynomial /summary printm nom;
run;

title "UNNORMALIZED PROFILE TRANSFORATION";
proc glm data=pigs1;
  class dose; 
  model week1 week3 week4 week5 week6 week7 = dose / nouni;
  repeated week 6 (1 3 4 5 6 7) profile /summary printm nom ;
run;

title "UNNORMALIZED HELMERT TRANSFORATION";
proc glm data=pigs1;
  class dose; 
  model week1 week3 week4 week5 week6 week7 = dose / nouni;
  repeated week 6  helmert /summary printm nom;
run;

/*******************************************************************

  Use PROC GLM to carry out the multivariate analysis and profile
  analysis, respectively.  The description is exactly the same as
  for Example 1 (the dental study).  In the first call, we also show
  use of the MEANS statement to calculate the means for each dose
  group at each time.
  
*******************************************************************/

title "MULTIVARIATE ANALYSIS USING PROC GLM MANOVA STATEMENT";
proc glm data=pigs1;
  class dose;
  model week1 week3 week4 week5 week6 week7 = dose / nouni;
  means dose;
  manova h=dose / printh printe;
run;

title "MULTIVARIATE PROFILE ANALYSIS";
proc glm data=pigs1;
  class dose;
  model week1 week3 week4 week5 week6 week7 = dose / nouni;
  repeated week / printe nou;
run;

