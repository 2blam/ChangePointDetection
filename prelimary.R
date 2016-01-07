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

#change in variance
x_chg_var = c(rnorm(150, 0, 1), rnorm(200, 0, 5), rnorm(150, 0, 1))
data_chg_var = data.frame(time=t, val=x_chg_var)

#change in regression
x_chg_reg = c(rnorm(150, 0, 1), rnorm(200, 0, 1) + seq(0, 6, length.out = 200) , rnorm(150, 6, 1))
data_chg_reg = data.frame(time=t, val=x_chg_reg)

#change in dependence
x1 = rnorm(125, 0, 1)
x2 = x1 + rnorm(125, 0, 0.00005)
x_dep =unlist(Map(c, x1, x2))
x_chg_dep = c(rnorm(250, 0, 1), x_dep)
data_chg_dep = data.frame(time=t, val=x_chg_dep)



#multiplot
p1 = ggplot(data_chg_mean, aes(x=time, y=val)) + 
  labs(x="Time", y="Value") +
  geom_line(color="blue") + 
  geom_vline(xintercept=250, color="red", linetype="longdash") + 
  ggtitle("Change in mean")

p2 = ggplot(data_chg_var, aes(x=time, y=val)) + 
  labs(x="Time", y="Value") +
  geom_line(color="blue") + 
  geom_vline(xintercept=c(150, 350), color="red", linetype="longdash") + 
  ggtitle("Change in variance")

p3 = ggplot(data_chg_reg, aes(x=time, y=val)) + 
  labs(x="Time", y="Value") +
  geom_line(color="blue") + 
  geom_vline(xintercept=c(150, 350), color="red", linetype="longdash") + 
  ggtitle("Change in regression")

p4 = ggplot(data_chg_dep, aes(x=time, y=val)) + 
  labs(x="Time", y="Value") +
  geom_line(color="blue") + 
  geom_vline(xintercept=250, color="red", linetype="longdash") + 
  ggtitle("Change in dependence")

multiplot(p1, p2, p3, p4, cols=2)


#################################################################
# get the sample data
filelist = dir(paste0(getwd(), '/RouteDataPreprocessed/', sep=""), pattern="*.mat");

ifn = filelist[1];#get the filename of the first file

#read the mat data
data = readMat(paste0(getwd(), '/RouteDataPreprocessed/', ifn , sep=''));
RTT = as.data.frame(data$data[1]);
t = unlist(data$data[2]);
colnames(RTT) = hosts;
start = t[1];
finish = tail(t, n=1); #get the last element

t = seq(start, finish, 300);#tick - 5 minutes

#we use BU data at the moment
data_BU = data.frame(time=t, val=RTT$BU)


startup_list = c(20, 40, 60, 80)
ARL_list = c(100, 500, 1000, 5000)
algo_list = c("Mann-Whitney", "Mood", "Lepage", "Kolmogorov-Smirnov")
algo_abbr_list = c("MW", "M", "L", "KS")

#for each startup
for (startup in startup_list){
  #for each ARL 
  for (ARL in ARL_list){
    #for each algo 
    for (algo in algo_list){
      algo_abbr = algo_abbr_list[match(algo, algo_list)]
      cat(sprintf(paste(startup, ARL, algo_abbr, " ")))
      #initalization
      detectiontimes = numeric()
      changepoints = numeric()
      
      cpm = makeChangePointModel(cpmType=algo, ARL0=ARL, startup=startup)
      
      i = 0
      #get BU RTT
      diffs = diff(RTT$BU) #first order difference
      
      while(i < length(diffs)){
        i = i + 1
        cpm = processObservation(cpm, diffs[i])
        if (changeDetected(cpm)){
          #cat(sprintf("change detected at %d\n", i))
          detectiontimes = c(detectiontimes, i)
          Ds = getStatistics(cpm)
          tau = which.max(Ds)
          if (length(changepoints) > 0){
            tau = tau + changepoints[length(changepoints)]
          }
          changepoints = c(changepoints, tau)
          cpm = cpmReset(cpm)
          i = tau
        }
      }
      
      #print number of changepoints
      cat(sprintf("%d ", length(changepoints)))
      
      #calculate the mean detection delay
      cat(sprintf("%f \n", mean(detectiontimes - changepoints) * 5))
      
      #for each change points, plot the related time series
      #for (idx in 1:length(changepoints)){
      # plot the first 10 changepoints
      for (idx in 1:10){
        out_fn = paste0("BU", "_", gsub(".mat", "", ifn), "_ARL0_", ARL, "_s_", startup, "_", algo_abbr, "_")
        
        out_fn = paste0(out_fn, "_", idx)
        #get the changepoints location & detection time
        cp = changepoints[idx]
        dt = detectiontimes[idx]
        delta = dt - cp
        
        #extract time series for plotting purpose
        ts = seq((cp-50), (cp+50))
        
        #ensure ts is positive
        ts = ts[(ts > 0)]
        
        #check which cp with the time series
        cps = intersect(changepoints, ts)
        dts = detectiontimes[match(cps, changepoints)]
        
        #if the last element is not within the ts, remove it
        if (is.na(match(tail(dts, 1), ts))){
          ridx = which(is.na(match(dts, ts)))
          dts = dts[-ridx]
        }
        
        data = data_BU[ts, ]
        x_axis = seq(1, length(ts))
        data$x_axis = x_axis
        #plot the graph
        png(out_fn)
        p = ggplot(data, aes(x=x_axis, y=val)) + 
          labs(x="Time", y="Value") +
          geom_line(color="black") + 
          geom_vline(xintercept=match(cps, ts), color="red", linetype="longdash", size=1.5) +
          geom_vline(xintercept=match(dts, ts), color="blue", linetype="longdash") +
          ggtitle(paste0("Mann-Whitney Algorithm", " ARL0:", ARL, " startup:", startup,
                         "\nchangepoint at ", as.POSIXct(t[cp], origin="1970-01-01"), "(red dashline)", 
                         "\ndetection time at ", as.POSIXct(t[dt], origin="1970-01-01"), "(blue dashline)",
                         "\ntime diff:", delta * 5, " mins"
          ))
        
        print(p)
        dev.off()
      }
    }
  }
}