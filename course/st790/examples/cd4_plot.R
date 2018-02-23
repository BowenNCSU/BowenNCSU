#########################################################################
#
#   CHAPTER 5, ACTG 193A Data
# 
#   Data from ACTG 193A with 4 treatment groups   
#
#   1 = zidovudine alternating monthly with 400mg didanosine
#   2 = zidovudine plus 2.25mg of zalcitabine
#   3 = zidovudine plus 400mg of didanosine
#   4 = zidovudine plus 400mg of didanosine plus 400mg of nevirapine
#
#########################################################################

#  Read in the data -- they are in the "long" format of one record per observation

thedat <- read.table("cd4.dat")
colnames(thedat) <- c("id","trt","age","gender","week","logcd4")

#  Total number of individuals

m <- max(thedat$id)

#########################################################################
#
#   Make separate plots by treatment group
#
#########################################################################

library(ggplot2)

#   Make the treatment and id variables factors for use with ggplot 

thedat <- within(thedat, {  id <- factor(id)
    trt <- factor(trt,levels=1:4,labels=c("ZDV+ alt ddI","ZDV+ZAL","ZDV+ddI","ZDV+ddI+NVP"))
})

#   First create the basic plot object

pp <- ggplot(thedat,aes(x=week,y=logcd4,group=id)) +
    theme(panel.background = element_rect(fill = "white"),
  panel.grid.major = element_line(colour = "white"),
  panel.grid.minor = element_line(colour = "white"))

#########################################################################
#
#   Now add features to the basic plot object
#
#   facet_wrap makes separate plots for each level of the factor trt
#   and labels them as specified for the factor trt above
#
#   geom_line() creates the spaghetti connecting the points for each 
#   subject -- by default solid lines are used
#
#   geom_point() adds plotting symbols - shape=18 specifies the diamond
#   shape; the default is solid circles
#
#   geom_smooth with the "group=1" aesthetic creates a loess smooth
#   through all data for each treatment (specified by facet_wrap)
#
#   the axes are labeled in the obvious way
#
#######################################################################

pdf("cd4_bytrt.pdf",width=8)
pp +  geom_line(color='dark gray') + geom_point(shape=18,color='dark gray') +
    geom_smooth(aes(group = 1),method="loess",color="black",se=FALSE) + facet_wrap(~ trt, ncol=2) + xlab("Week") + ylab("log CD4")
dev.off()


