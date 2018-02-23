#########################################################################
#
#   CHAPTER 1, EXAMPLE 1, Dental Study 
# 
#   How to make a spaghetti plot.  We demonstrate using the R function
#   interaction.plot() and the ggplot2 package.  There are also
#   numerous examples on the web.  Just google "spaghetti plot in R"
#
#########################################################################

#  Read in the data -- they are in the form of one record per observation

thedat <- read.table("dental.dat")
thedat <- thedat[,2:5]      #  remove the first column
colnames(thedat) <- c("id","age","distance","gender")

#  Total number of individuals

m <- max(thedat$id)

#########################################################################
# 
#  Use the interaction.plot function to plot all the data on one plot
#  and use the gender indicator (0 for girls, 1 for boys) as the plotting
#  symbol.
#
#  By default, only lines are plotted, with a different line
#  type for each individual.  Here, we make all the lines solid with lty=1
#  
#  By default, a legend is displayed; here, we suppress the legend.
#
#  We specify plotting of both lines and symbols with type="b"
#
#  Below, gensym is a vector of length m with the gender indicator for each
#  child -- these will be the plotting symbols
#
#########################################################################

gensym <- thedat[thedat$age==8,4]

pdf("dental_spaghetti.pdf",width=8)
interaction.plot(thedat$age,thedat$id,thedat$distance,xlab="Age (years)",
                 ylab="Distance (mm)",lty=1,legend=FALSE,type="b",
                 pch=as.character(gensym))
title("Dental Study Data")
graphics.off()

#########################################################################
#
#   Now make separate plots by gender with sample means at each age
#   superimposed for each gender.  We illustrate how the ggplot2 package
#   can be used to build this plot.  
#
#########################################################################

library(ggplot2)

#   Make the gender and id variables factors for use with ggplot 

thedat <- within(thedat, {  id <- factor(id)
    gender <- factor(gender,levels=0:1,labels=c("Girls","Boys"))
})

#   First create the basic plot object

pp <- ggplot(thedat,aes(x=age,y=distance,group=id))

#########################################################################
#
#   Now add features to the basic plot object
#
#   facet_grid makes separate plots for each level of the factor gender
#   and labels them as specified for the factor gender above
#
#   geom_line() creates the "spaghetti connecting the points for each 
#   child -- by default solid lines are used
#
#   geom_point() adds plotting symbols - shape=18 specifies the diamond
#   shape; the default is solid circles
#
#   stat_summary() adds the means for each gender and connects them by 
#   lines -- the size option specifies thicker lines than the default
#
#   the axes are labeled in the obvious way
#
#######################################################################

pdf("dental_bygender.pdf",width=8)
pp + geom_line() + geom_point(shape=18) + stat_summary(aes(group = 1),
   geom = "line", fun.y = mean, size = 1.5) + facet_grid(. ~ gender) +
        xlab("Age (years)") + ylab("Distance (mm)")
dev.off()


