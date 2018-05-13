/***********************************************************

   Homework 4, Problem 3

   Population-averaged analysis of the leprosy data using
   SAS PROC GENMOD, PROC GEE, and PROC GLIMMIX
    
***********************************************************/

options ps=59 ls=80 nodate; run;

data leprosy2; 
  infile "leprosy.dat";
  input id trt pre post;
run;

*  reconfigure in long format;

data leprosy; set leprosy2;
    array wt(2) pre post;
    do time = 1 to 2;
        count = wt(time);
        output;
        end;
    drop pre post;

data leprosy; set leprosy;
    time = time -1;
run;

proc print data = leprosy(obs=10); run;

*   sample means and variances for each drug at each time;

proc means mean var noprint data=leprosy;
   class trt time;
   var count;
   output out=lepsumm mean=mean var=var;
run;

data lepsumm; set lepsumm;
  if trt=. then delete;
  if time=. then delete;
  drop _TYPE_ _FREQ_;
run;

title "Pre/PostSample Means and Variances";
proc print data=lepsumm;
run;
    
/***********************************************************

  We fit the model with each procedures; the syntax is
  identical or very similar.  We use the log link to get
  the loglinear model and the Poisson distribution to get
  variance proportional to the mean (by a dispersion or 
  or scale parameter).

  Because there are only 2 observations per individual, we
  can fit the correlation structure using compound symmetric
  or unstructured; we try both.

  By default, SAS will make the treatment with the "largest"
  factor level (2 in this case) the reference treatment,
  but the model we want to fit takes the placebo (0) as the
  reference.  We specify the reference treatment as 0 in
  the CLASS statement.
    
***********************************************************/

title "Unstructured correlation, PROC GENMOD";
proc genmod data=leprosy;
  class id trt(ref="0");
  model count = time time*trt / link=log dist=poisson;
  repeated subject=id / type=un corrw modelse;
run;

title "Compound symmetric correlation, PROC GENMOD";
proc genmod data=leprosy;
  class id trt(ref="0");
  model count = time time*trt / link=log dist=poisson;
  repeated subject=id / type=cs corrw modelse;

*  2 df Wald test of antibiotic effect;
  
  contrast "No antibiotic effect" time*trt 1 0 -1,
                                  time*trt 0 1 -1 / wald;
run;

title "Compound symmetric correlation, PROC GEE";    
proc gee data=leprosy;
  class id trt(ref="0");
  model count = time time*trt / link=log dist=poisson;
  repeated subject=id / type=cs corrw modelse;
run;

*  Using PROC GLIMMIX with quadratic estimating equation;

title "Compound symmetric correlation, PROC GLIMMIX";
proc glimmix data=leprosy method=mmpl;
  class id trt(ref="0");
  model count = time trt*time / solution link=log dist=poisson;
  random _residual_ / subject=id type=cs vcorr;
run;

title"Using REML type equation instead";
proc glimmix data=leprosy method=rmpl;
  class id trt(ref="0");
  model count = time trt*time / solution link=log dist=poisson;
  random _residual_ / subject=id type=cs vcorr;
run;




