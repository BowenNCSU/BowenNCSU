/*******************************************************************

  CHAPTER 3, EXAMPLE 1, Dental Study

  Univariate and multivariate repeated measures analysis of variance

*******************************************************************/

options ls=80 ps=59 nodate; run;

/******************************************************************

  Read in the data (in "long" form of one record per observation)
  and create an alternative data set in "wide" form with one record
  per individual

*******************************************************************/

data dent1; infile 'dental.dat';
  input obsno child age distance gender;
run;

proc transpose data=dent1 out=dent2 prefix=age;
   by gender child notsorted;
   var distance; 
run;

/*******************************************************************

  Construct the univariate analysis of variance using PROC GLM  
  via a "split plot" specification, which requires the data in
  "long" form (DENT1).

  The F ratio that PROC GLM prints out automatically for the 
  gender effect (averaged across age) uses the MSE in the 
  denominator, so is incorrect.

  The RANDOM statement produces the expected mean squares

  The TEST option request the test for the gender effect
  treating the child(gender) efffect as random, yielding the correct
  F ratio.  Other F-ratios are correct.
  
*******************************************************************/  

title "UNIVARIATE ANALYSIS VIA SPLIT PLOT SPECIFICATION USING PROC GLM";
proc glm data=dent1;
  class age gender child;
  model distance = gender child(gender) age age*gender;
  random child(gender) / test;
run;

/*******************************************************************

  This analysis can also be done using PROC MIXED.  The 
  MODEL and RANDOM statements have DIFFERENT interpretations in
  PROC GLM and PROC MIXED.  The METHOD=TYPE3 option produces the
  correct tests.

*******************************************************************/  

title "UNIVARIATE ANALYSIS USING PROC MIXED";
proc mixed data=dent1 method=type3;
  class age gender child;
  model distance = gender age age*gender;
  random child(gender);
run;

/*******************************************************************

  Now carry out the same analysis using the REPEATED statement in
  PROC GLM.  This requires that the data be represented in the "wide" 
  form as in DENT2.

  The option NOUNI suppresses individual analyses of variance
  for the data at each age value from being printed.

  The PRINTE option asks for the test of sphericity to be performed.

  The NOM option means "no multivariate," which means just do 
  the univariate repeated measures analysis under the assumption 
  that the exchangable (compound symmetry) model is correct.
  
*******************************************************************/  

title "UNIVARIATE ANALYSIS USING REPEATED STATEMENT IN PROC GLM";
proc glm data=dent2;
  class gender;
  model age1 age2 age3 age4 = gender / nouni;
  repeated age / printe nom;
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

  First request the orthogonal polynomial transformation.  Here, SAS
  DOES use the normalized version of the U matrix.  Thus, the SSs 
  from the individual ANOVAs for each column will add up to the Gender 
  by Age interaction SS (and similarly for the within-unit error SS).
    
*******************************************************************/  

title "ORTHOGONAL POLYNOMIAL TRANSFORMATION";
proc glm data=dent2;
  class gender;
  model age1 age2 age3 age4 = gender / nouni;
  repeated age 4 (8 10 12 14) polynomial /summary nou nom printm;
run;

/*******************************************************************

  Now request the Helmert transformation.  SAS does NOT use the
  normalized version of the Helmert transformation matrix.  Thus,
  the SSs from the individual ANOVAs for each column will NOT add
  up to the Gender by Age interaction SS (similarly for within-unit
  error).  However, the F ratios are correct.  

*******************************************************************/  

title "UNNORMALIZED HELMERT TRANSFORATION";
proc glm data=dent2;
  class gender;
  model age1 age2 age3 age4 = gender / nouni;
  repeated age 4 (8 10 12 14) helmert /summary nou nom printm;
run;

/*******************************************************************

  Here, we manually perform the same analysis using the
  NORMALIZED version of the Helmert transformation matrix.
  We get each individual test separately using the PROC GLM
  MANOVA statement.  

********************************************************************/

title "NORMALIZED HELMERT TRANSFORATION";
proc glm data=dent2; 
  model age1 age2 age3 age4 = gender /nouni;
  manova h=gender m=0.866025404*age1 - 0.288675135*age2- 0.288675135*age3 - 0.288675135*age4;
  manova h=gender m=  0.816496581*age2-0.40824829*age3-0.40824829*age4;
  manova h=gender m=  0.707106781*age3-  0.707106781*age4;
run;

/*******************************************************************

  To compare, we apply the contrasts (normalized version) to each
  child's data.  We thus get a single value for each child corresponding
  to each contrast.  These are in the variables AGE1P -- AGE3P.
  We then use PROC GLM to perform each separate ANOVA.  It can be
  verified that the separate gender sums of squares add up to
  the interaction SS in the analysis above.  

********************************************************************/

data dent3; set dent2;
  age1p = sqrt(3/4)*(age1-age2/3-age3/3-age4/3);
  age2p = sqrt(2/3)*(age2-age3/2-age4/2);
  age3p = sqrt(1/2)*(age3-age4);
run;

proc glm; class gender; model age1p age2p age3p = gender;
run;

/*******************************************************************

  Now carry out the multivariate analysis.

  The MANOVA statement yields the test of equality of gender mean
  vectors, which is equivalent to Hotelling's T^2 test in this case
  because there are 2 groups.

  The PRINTH and PRINTE options print the SS&CP matrices
  Q_H and Q_E corresponding to the null hypothesis of equal means.
 
  Without the NOUNI option, the individual analyses of variance
  for the data at each age value are printed.   
  
********************************************************************/

title "MULTIVARIATE ANALYSIS USING PROC GLM MANOVA STATEMENT";
proc glm data=dent2;
  class gender;
  model age1 age2 age3 age4 = gender;
  manova h=gender / printh printe;
run;

/*******************************************************************
 
  Use the REPEATED option to do profile analysis.  The
  "between subjects" (units) test is that for coincidence assuming
  profiles are parallel, based on averaging across times, so
  is th same as the univariate test.
  
  The tests for age and age*gender are the multivariate tests for
  profile constancy and parallelism, respectively.  The test for
  constancy (age effect here) is the multivariate test for constancy
  assuming that the profiles are parallel.  Both are different from
  the univariate tests, which are prediated on the assumption of
  compound symmetry.  

  The within-unit analyses using different contrast matrices are
  be the same as in the univariate case above so are not repeated.

*******************************************************************/  

title "MULTIVARIATE PROFILE ANALYSIS";
proc glm data=dent2;
  class gender;
  model age1 age2 age3 age4 = gender / nouni;
  repeated age / nou;
run;
