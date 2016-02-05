#eid,nid, tid
#54, 302, 391
#66, 301, 391
#78, 303, 391
setwd("~/Data/CSC_Backend/experimental")
library(zoo)
library(ggplot2)
library(reshape)

#hpf - high pass filter
hpf = function (X){
  p = length(X)
  if (p %% 2 == 0){
    #prepare resultant list 
    d_i = rep(0, (p/2))
    
    for (i in 1:(p/2)){
      d_i[i]  = (X[(2*i)-1] - X[(2*i)]) / sqrt(2)
    }
    
    unlist(d_i)
  }else{
    sprintf("%s", "The number of element must be even!")
    NA
  }
}
#lpf - low pass filter
lpf = function (X){
  p = length(X)
  if (p %% 2 == 0){
    #prepare resultant list 
    a_i = rep(0, (p/2))
    
    for (i in 1:(p/2)){
      a_i[i]  = (X[(2*i)-1] + X[(2*i)]) / sqrt(2)
    }
    
    unlist(a_i)
  }else{
    sprintf("%s", "The number of element must be even!")
    NA
  }
}
#ihpf - inverse high pass filter
ihpf = function(X){
  p = length(X) * 2
  D = rep(0, p)
  for (i in 1:(p/2)){
    D[2*i-1] = D[2*i-1] + X[i] * 1/sqrt(2)
    D[2*i] = D[2*i] - X[i] * 1/sqrt(2)
  }
  D
}
#ilpf - inverse low pass filter
ilpf = function(X){
  p = length(X) * 2
  A = rep(0, p)
  for (i in 1:(p/2)){
    A[2*i-1] = A[2*i-1] + X[i] * 1/sqrt(2)
    A[2*i] = A[2*i] + X[i] * 1/sqrt(2)
  }
  A
}

energyRatio = function(D_t, A_t){
  sum(D_t *D_t) / (sum(A_t*A_t) + sum(D_t*D_t))
}

#read data
fn1 = "e51_n302_t391_0900_15Jan_0900_21Jan.csv"
fn2 = "e66_n301_t391_0900_15Jan_0900_21Jan.csv"
fn3 = "e78_n303_t391_0900_15Jan_0900_21Jan.csv"

fn = c(fn1, fn2, fn3)

start_record_time = as.numeric(as.POSIXct("2016-01-15 09:00:00"))
end_record_time   = as.numeric(as.POSIXct("2016-01-21 09:00:00"))
ticks = seq(start_record_time, end_record_time, by=600)

for (i in 1:3){
  #read file
  csv = read.csv(fn[i])
  colnames(csv) = c("time", "RTT")
  
  #preprocessing
  csv = rbind(csv, data.frame(time=ticks, RTT=NA))
  csv = csv[order(csv$time), ]
  
  #interpolate missing data
  csv = as.data.frame(na.fill(na.approx(csv), "extend"))
  #remove duplicate
  csv = csv[!duplicated(csv), ]
  #extract aligned ticks only
  idx =which(csv$time %in% ticks)
  if (i == 1){
    df = csv[idx, ]
  }else{
     df = cbind(df, csv[idx, "RTT"])
  }
}

colnames(df) = c("time", "RTT", "RTT2", "RTT3")
meltDf = melt(df, id=c("time"))

#prepare x axis label
idx = seq(from=1, to=length(ticks), by=144)
bks = ticks[idx]
xlbls = as.POSIXct(ticks, origin="1970-01-01")[idx]

meltDf$time = as.POSIXct(meltDf$time, origin="1970-01-01")

#plot
ggplot(meltDf, aes(x =time, y=value)) + 
  geom_line(aes(color=variable, group=variable)) +
  scale_shape_discrete(name  ="variable",
                        breaks=c("RTT", "RTT2","RTT3"),
                        labels=c("RTT", "RTT2","RTT3")) +
  theme(axis.text.x = element_text(angle=90, hjust=1))
  

#ZERO padding - ensure to have even number of elements
df$RTT4 = 0 #because we only have 3 paths to the destination, add 1 more with ZERO value
alpha = 0.2
beta = 0.2
w = 10
k = 1.96
E_hat = rep(0, w)
R_i = rep(0, w)
numOfPath = 3
r = 2
alert_time = c()
alert_idx = c()
#harr wavelet transform
#prepare LPF
for (t in 1:nrow(df)){
  X_t = df[t, 2:ncol(df)]
  current_time = df[t, 1]
  #decompose
  a_i = lpf(X_t)
  d_i = hpf(X_t)
  
  #reconstruction
  A_t = ilpf(a_i)
  D_t = ihpf(d_i)
  
  
  #monitor the energy of fluctutation signal (D) 
  # with respect to the total energy of both signals 
  #E = sum(D *D) / (sum(A*A) + sum(D*D))
  E_i = energyRatio(D_t, A_t)
  
  #Holt - Winters  (non-seasonal series)
  #reference: http://www.slideshare.net/vishalkukreja376/chap19-time-seriesanalysisandforecasting
  if (t == 1){
    E_1 = E_i
    df_history = X_t
  }else if (t == 2){
    E_hat_prev = E_i
    T_prev = E_i - E_1
    df_history = rbind(df_history, X_t)
  }else{ #t >=3
    E_hat_i = alpha * (E_hat_prev + T_prev) + (1 - alpha) * E_i
    T_i = beta * T_prev + (1 - beta) * (E_hat_i - E_hat_prev)
    #calculate the residual
    residual = E_hat_i - E_i
    
    #store in list
    if (t <= w){
      R_i[t] = residual
      df_history = rbind(df_history, X_t)
    }else{
      #enough element for checking
      phi_w = median(R_i)
      std_w = 1.4826 * median(abs(R_i - phi_w))
      
      LHS = abs(residual - phi_w)
      RHS = (k * std_w)     
      #print(LHS)
      #print(RHS)
      #print("--")
      if (LHS >= RHS){
        #spatial anomaly, further checking temporal 
        #for each of the attribute / channel
        #print("spatial")
        alarm = 0
        for (j in 1:numOfPath){
          myboxplot <- boxplot(df_history[, j], plot = FALSE) 
          
          if (length(myboxplot$out) > 1){
            alarm = alarm + 1
          }
        }
        
        #print(alarm)
        if (alarm >= r){
          print(t)
          print("alert")
          alert_time = append(alert_time, current_time)
          alert_idx = append(alert_idx, t)
        }
      }
      
      #append the data and remove the 1st element
      R_i = append(R_i, residual)
      R_i = R_i[-1]
      
      df_history = rbind(df_history, X_t)
      df_history = df_history[-c(1), ]
    }
    
    E_hat_prev = E_hat_i
    T_prev = T_i
  }
}




#plot

alert_time_p = as.POSIXct(alert_time, origin="1970-01-01")
p = ggplot(meltDf, aes(x =time, y=value)) + 
  geom_line(aes(color=variable, group=variable)) +
  scale_shape_discrete(name  ="variable",
                       breaks=c("RTT", "RTT2","RTT3"),
                       labels=c("RTT", "RTT2","RTT3")) +
  theme(axis.text.x = element_text(angle=90, hjust=1)) 
  p + geom_vline(aes(xintercept = alert_time))
  