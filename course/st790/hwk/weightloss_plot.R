#########################################################################
#
#   Weightloss data - three weightloss programs
#
#   1 = Control
#   2 = Diet + Exercise
#   3 = Diet Alone
#
#   Make spaghetti plots
#
#########################################################################

#  Read in the data -- they are in the form of one record per individual
#  aka the "wide" format

thedat.wide <- read.table("weightloss.dat",row.names=NULL)
colnames(thedat.wide) <- c("id","month0","month1","month2","month3","month4","program")

#  Total number of individuals

m <- nrow(thedat.wide)

#  Reconfigure as one observation per record ("long" format) and sort
#  in id order

thedat.long <- reshape(thedat.wide,varying=c("month0","month1","month2","month3","month4"),
            v.names="weight",idvar="id",times=c(0,3,6,9,12),timevar="month",direction="long")

thedat <- thedat.long[order(thedat.long$id),]

#########################################################################
#
#   Make separate plots by treatment group
#
#########################################################################

library(ggplot2)

#   Make the treatment and id variables factors for use with ggplot 

thedat <- within(thedat, {  id <- factor(id)
    program <- factor(program,levels=1:3,labels=c("Control","Diet+Exercise","Diet Alone"))
})

#   First create the basic plot object

pp <- ggplot(thedat,aes(x=month,y=weight,group=id)) +
    theme(panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "gray94"),
    panel.grid.minor = element_line(colour = "gray94"),
    panel.spacing=unit(1, "lines")) +scale_x_continuous(breaks=seq(0,12,3))

#########################################################################
#
#   Now add features to the basic plot object
#
#   facet_wrap makes separate plots for each level of the factor program
#   and labels them as specified for the factor program above
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

pdf("weightloss_byprogram.pdf",width=8)
qq <- pp +  geom_line(color='dark gray') + geom_point(shape=18,color='dark gray') +
#    geom_smooth(aes(group = 1),method="loess",color="black",se=FALSE)
     stat_summary(aes(group = 1),geom = "line", fun.y = mean, size = 1.5) +                       
     facet_wrap(~ program, ncol=2) + xlab("Month") + ylab("Weight (lb)")
    print(qq)
dev.off()


