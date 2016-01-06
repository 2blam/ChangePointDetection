#Sys.setlocale("LC_ALL","English") #ensure the plot output is in English
library(R.matlab) #read mat data file
library(cpm) #change point model library
library(ggplot2)
setwd("/home/user/Data/Research/Experiments")
hosts = c('BU','CityU', 'CU', 'HKU', 'IED', 'LN', 'PolyU', 'UST');

#multiplot function - http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# generate four different changepoints in time series
# change in mean
# change in variance
# change in regression
# change in dependence

set.seed(0)
t = seq(1, 500)

#change in mean
x_chg_mean = c(rnorm(250, 0, 1), rnorm(250, 4, 1))
data_chg_mean = data.frame(time=t, val=x_chg_mean)
#plot(x_chg_mean, type="l", xlab="t", ylab="x", bty="l",  col="blue")
#abline(v = 250, lty=2,  col="red")

#change in variance
x_chg_var = c(rnorm(150, 0, 1), rnorm(200, 0, 5), rnorm(150, 0, 1))
data_chg_var = data.frame(time=t, val=x_chg_var)
#plot(x_chg_var, type="l", xlab="t", ylab="x", bty="l",  col="blue")
#abline(v = c(150, 350), lty=2, col="red")

#change in regression
x_chg_reg = c(rnorm(150, 0, 1), rnorm(200, 0, 1) + seq(0, 6, length.out = 200) , rnorm(150, 6, 1))
data_chg_reg = data.frame(time=t, val=x_chg_reg)
#plot(x_chg_reg, type="l", xlab="t", ylab="x", b`ty="l",  col="blue")
#abline(v = c(150, 350), lty=2, col="red")

#change in dependence
x1 = rnorm(125, 0, 1)
x2 = x1 + rnorm(125, 0, 0.00005)
x_dep =unlist(Map(c, x1, x2))
x_chg_dep = c(rnorm(250, 0, 1), x_dep)
data_chg_dep = data.frame(time=t, val=x_chg_dep)
#plot(x_chg_dep, type="l", xlab="t", ylab="x", bty="l",  col="blue")
#abline(v = 250, lty=2,  col="red")


#multiplot
p1 = ggplot(data_chg_mean, aes(x=time, y=val)) + 
  geom_line(color="blue") + 
  geom_vline(xintercept=250, color="red", linetype="longdash") + 
  ggtitle("change in mean")

p2 = ggplot(data_chg_var, aes(x=time, y=val)) + 
  geom_line(color="blue") + 
  geom_vline(xintercept=c(150, 350), color="red", linetype="longdash") + 
  ggtitle("change in variance")

p3 = ggplot(data_chg_reg, aes(x=time, y=val)) + 
  geom_line(color="blue") + 
  geom_vline(xintercept=c(150, 350), color="red", linetype="longdash") + 
  ggtitle("change in regression")

p4 = ggplot(data_chg_dep, aes(x=time, y=val)) + 
  geom_line(color="blue") + 
  geom_vline(xintercept=250, color="red", linetype="longdash") + 
  ggtitle("change in dependence")

multiplot(p1, p2, p3, p4, cols=2)