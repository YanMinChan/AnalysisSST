# F70DB Statistical Dissertation B
# Due date: 30 March 2023
#--------------------
# Necessary library 
#--------------------
library(raster)
library(ncdf4)
library(forecast)
library(zoo)
library(tseries)
library(ggplot2)
library(astsa)
#--------------------------------------------------------
# Simulating SARIMA model and their ACF and PACF plots
#--------------------------------------------------------
s1 = sarima.sim(ar = c(1.5,-0.75), sar = -0.8, n = 500, S = 12)
par(mfrow=c(1, 2))
acf(s1, xaxt = "n"); axis(1, at=0:500/12, labels=0:500)
pacf(s1, xaxt = "n"); axis(1, at=0:500/12, labels=0:500)

#-------------------------------
# Whittle's estimate functions
#-------------------------------

# General spectral density of SGAR with seasonality s
spec = function(omega, alpha, delta, sigsq, s){
  val = ((1-2*alpha*cos(s*omega)+alpha^2)^-delta)*(sigsq/(2*pi))
  return(val)
}

# Periodogram
per = function(omega){
  sI = colSums((sst.ts[1:n]-mean(sst.ts))*exp(1i*omega%*%t(c(0:(n-1)))))
  per = (1/(2*pi*n))*(abs(sI)^2)
  return(per)
}

# Whittle's estimation of parameter alpha, delta and sigma^2
whit = function(init){
  s1 = sum(log(spec(2*pi*c(0:(n-1))/n, init[1], init[2], init[3], s)))
  s2 = sum(per(2*pi*c(0:(n-1))/n)/spec(2*pi*c(0:(n-1))/n, init[1], init[2], init[3], s))
  Est = 2*n*log(2*pi)+s1+s2
  return(Est)
}

neg_whit = function(init){
  s1 = sum(log(spec(2*pi*c(0:(n-1))/n, init[1], init[2], init[3], s)))
  s2 = sum(per(2*pi*c(0:(n-1))/n)/spec(2*pi*c(0:(n-1))/n, init[1], init[2], init[3], s))
  Est = 2*n*log(2*pi)+s1+s2
  return(-Est)
}

#------------------------------------------
# Functions for simulation of SGAR process
#------------------------------------------

# Simulate SGAR process
ft = function(t, k, alpha, delta, s, wn){
  val = numeric(k+1)
  for (j in 0:k){
    if ((t - s*j) > 0){
      val[j+1] = ((gamma(j + delta)*(alpha^j))/(gamma(delta)*gamma(j + 1)))*wn[t-s*j]
    } else{
      val[j+1] = 0
    }
  }
  return(val)
}

gen_SGAR = function(n, alpha, delta, sig, s){
  Xt = numeric(n)
  wn = rnorm(10000, 0, sqrt(sig))
  for (i in 1:n){
    Xt[i] = sum(ft(i, 100, alpha, delta, s, wn))
  }
  return(Xt)
}

# Estimate white noise
ft2 = function(t, k, alpha, delta, s, wn){
  val = numeric(k)
  for (j in 1:k){
    if ((t - s*j) > 0){
      val[j] = ((gamma(j + delta)*(alpha^j))/(gamma(delta)*gamma(j + 1)))*wn[t-s*j]
    } else{
      val[j] = 0
    }
  }
  return(val)
}

est_wn = function(n, alpha, delta, s, Xt){
  wn = numeric(n)
  for (i in 1:n){
    wn[i] = Xt[i] - sum(ft2(i, 100, alpha, delta, s, wn))
  }
  return(wn)
}

# Forecast of future data points
forecast_SGAR = function(n, Yt, alpha, delta, sig, s, et){
  Xt = numeric(n)
  for (i in 1:n){
    et = append(et, rnorm(1, sd = sqrt(sig)))
    Xt[i] = sum(ft(length(Yt)+i, 100, alpha, delta, s, et))
  }
  return(Xt + mean(Yt))
}

#---------------------------------------
# Simulation to check Whittle function
#---------------------------------------
# SGAR check
init.gen = c(0.9, 0.8, 1)
nsim = 200
param = matrix(0, ncol = length(init.gen), nrow = nsim)
for (i in 1:nsim){
  # Simulate a SGAR process
  ngen = 1000; s = 12
  sst = gen_SGAR(ngen, init.gen[1], init.gen[2], init.gen[3], s)
  sst.ts = ts(sst[500:1000]); n = length(sst.ts)
  
  # Parameter estimate
  init.test = c(0, 0.5, var(sst.ts))
  opti1 = nlminb(init.test, whit, lower = c(-0.9999, 0.0001, 0.0001), upper = c(0.9999, 1, Inf)); opti1
  param[i,] = c(opti1$par[1], opti1$par[2], opti1$par[3])
}

# Numerical summaries
num_summary = function(param, true){
  n = ncol(param)
  q = matrix(0, ncol = 7, nrow = n)
  mean_theta = numeric(n)
  bias_theta = numeric(n)
  est_se = numeric(n)
  est_rmse = numeric(n)
  for (i in 1:n){
    q[i,] = quantile(param[,i], probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.90, 0.95))
    mean_theta[i] = mean(param[,i])
    bias_theta[i] = mean_theta[i]-true[i]
    est_se[i] = sqrt((1/(n-1))*sum((param[,i]-mean_theta[i])^2))
    est_rmse[i] = sqrt((1/n)*sum((param[,i]-true[i])^2))
  }
  
  ret = list("5th, 10th, 90th & 95th quantile and median" = q, "Mean" = mean_theta, "Estimated bias" = bias_theta, 
             "Estimated standard error" = est_se, "Estimated root mean square error" = est_rmse)
  return(ret)
}

num_summary(param, init.gen)

# Plot of density
plot(density(param[,1]), xlim = c(0, 2), ylim = c(0, 10), col = "blue", main = expression(paste("Density plot of estimated" ~ alpha, ","  ~ delta, "," ~ sigma^2)))
lines(density(param[,2]), col = "red")
lines(density(param[,3]), col = "green")

#--------------------
# Necessary data set
#--------------------
hadisst = brick(x="HadISST_sst.nc", values=T, varname='sst'); hadisst

# Labuan coordinates
llon = 115
llat = 5
extract.coor = matrix(c(llon, llat), nrow = 1, ncol = 2)

# Extract train and test data set
dt = hadisst[[961:1824]]
ext = extract(dt, extract.coor, method = "bilinear")
sst.train = ext[1:804]
sst.test = ext[805:length(ext)]
sst.ts = ts(sst.train, frequency = 12, start = c(1950, 1))
sst.ts.test = ts(sst.test, frequency = 12, start = c(2017, 1))

n = length(sst.ts); s = 12

#-----------------------------------------
# Application of SARIMA model on SST data
#-----------------------------------------
# If running this section after running seasonal differenced GAR section,
# necessary dataset needed to be run again before proceeding 

sst.ts.df1 = diff(sst.ts, lag = 12) # Deseasonalized and detrended data

# Observe increasing trend
plot(sst.ts, main = "Time series plot for SST", ylab = "Temperature", xlab = "Year")
lt = lm(sst.ts~index(sst.ts))
abline(lt$coefficients[[1]], lt$coefficients[[2]], col = "red")

plot(sst.ts.df1, main = "Deseasonalized and detrended time series plot for SST", ylab = expression(X[t]), xlab = "Year")
lt2 = lm(sst.ts.df1~index(sst.ts.df1))
abline(lt2$coefficients[[1]], lt2$coefficients[[2]], col = "red")

# Checking stationality with unit root test
adf.test(sst.ts.df1) #stationary
nsdiffs(sst.ts)
ndiffs(nsdiffs(sst.ts))

# Observe acf and pacf
par(mfrow= c(1, 2))
acf(sst.ts.df1, xaxt = "n", xlab = "Lag(months)")
axis(1, at=0:804/12, labels=0:804)
pacf(sst.ts.df1, xaxt = "n", xlab = "Lag(months)")
axis(1, at=0:804/12, labels=0:804)

# Fit sarima model
#------------------
# AR(p) and SAR(P) look at pacf, MA(q) and SMA(Q) look at acf
t1 = Arima(sst.ts, order=c(4, 1, 1), seasonal = list(order=c(2, 1, 1), period = 12), optim.control = list(maxit = 1000));t1
t2 = Arima(sst.ts, order=c(4, 1, 2), seasonal = list(order=c(2, 1, 1), period = 12), optim.control = list(maxit = 1000));t2
t3 = Arima(sst.ts, order=c(3, 1, 1), seasonal = list(order=c(2, 1, 1), period = 12), optim.control = list(maxit = 1000));t3
t4 = Arima(sst.ts, order=c(3, 1, 2), seasonal = list(order=c(2, 1, 1), period = 12), optim.control = list(maxit = 1000));t4
t5 = Arima(sst.ts, order=c(2, 1, 1), seasonal = list(order=c(2, 1, 1), period = 12), optim.control = list(maxit = 1000));t5
t6 = Arima(sst.ts, order=c(2, 1, 2), seasonal = list(order=c(2, 1, 1), period = 12), optim.control = list(maxit = 1000));t6
t7 = Arima(sst.ts, order=c(4, 1, 1), seasonal = list(order=c(2, 1, 2), period = 12), optim.control = list(maxit = 1000));t7
t8 = Arima(sst.ts, order=c(4, 1, 2), seasonal = list(order=c(2, 1, 2), period = 12), optim.control = list(maxit = 1000));t8
t9 = Arima(sst.ts, order=c(3, 1, 1), seasonal = list(order=c(2, 1, 2), period = 12), optim.control = list(maxit = 1000));t9
t10 = Arima(sst.ts, order=c(3, 1, 2), seasonal = list(order=c(2, 1, 2), period = 12), optim.control = list(maxit = 1000));t10
t11 = Arima(sst.ts, order=c(2, 1, 1), seasonal = list(order=c(2, 1, 2), period = 12), optim.control = list(maxit = 1000));t11
t12 = Arima(sst.ts, order=c(2, 1, 2), seasonal = list(order=c(2, 1, 2), period = 12), optim.control = list(maxit = 1000));t12
t13 = Arima(sst.ts, order=c(4, 0, 1), seasonal = list(order=c(2, 1, 1), period = 12), optim.control = list(maxit = 1000));t13
t14 = Arima(sst.ts, order=c(4, 0, 2), seasonal = list(order=c(2, 1, 1), period = 12), optim.control = list(maxit = 1000));t14
t15 = Arima(sst.ts, order=c(3, 0, 1), seasonal = list(order=c(2, 1, 1), period = 12), optim.control = list(maxit = 1000));t15
t16 = Arima(sst.ts, order=c(3, 0, 2), seasonal = list(order=c(2, 1, 1), period = 12), optim.control = list(maxit = 1000));t16
t17 = Arima(sst.ts, order=c(2, 0, 1), seasonal = list(order=c(2, 1, 1), period = 12), optim.control = list(maxit = 1000));t17
t18 = Arima(sst.ts, order=c(2, 0, 2), seasonal = list(order=c(2, 1, 1), period = 12), optim.control = list(maxit = 1000));t18
t19 = Arima(sst.ts, order=c(4, 0, 1), seasonal = list(order=c(2, 1, 2), period = 12), optim.control = list(maxit = 1000));t19
t20 = Arima(sst.ts, order=c(4, 0, 2), seasonal = list(order=c(2, 1, 2), period = 12), optim.control = list(maxit = 1000));t20
t21 = Arima(sst.ts, order=c(3, 0, 1), seasonal = list(order=c(2, 1, 2), period = 12), optim.control = list(maxit = 1000));t21
t22 = Arima(sst.ts, order=c(3, 0, 2), seasonal = list(order=c(2, 1, 2), period = 12), optim.control = list(maxit = 1000));t22
t23 = Arima(sst.ts, order=c(2, 0, 1), seasonal = list(order=c(2, 1, 2), period = 12), optim.control = list(maxit = 1000));t23
t24 = Arima(sst.ts, order=c(2, 0, 2), seasonal = list(order=c(2, 1, 2), period = 12), optim.control = list(maxit = 1000));t24
#------------------

# Comparing best fit with error measure
# ARIMA(3, 1, 1)(2, 1, 1)[12]
forecast1 = forecast(t3, h = 60)
accuracy(forecast1, sst.ts.test)

# ARIMA(3, 1, 2)(2, 1, 1)[12]
forecast2 = forecast(t4, h = 60)
accuracy(forecast2, sst.ts.test)

# ARIMA(2, 1, 1)(2, 1, 1)[12]
forecast3 = forecast(t5, h = 60)
accuracy(forecast3, sst.ts.test)

# ARIMA(2, 1, 1)(2, 1, 2)[12]
forecast4 = forecast(t11, h = 60)
accuracy(forecast4, sst.ts.test)

# ARIMA(2, 0, 2)(2, 1, 1)[12]
forecast5 = forecast(t18, h = 60)
accuracy(forecast5, sst.ts.test)

# Best model
sst.fit = Arima(sst.ts, order=c(2, 0, 2), seasonal = list(order=c(2, 1, 1), period = 12), optim.control = list(maxit = 1000)); sst.fit
checkresiduals(sst.fit)

# Model diagnostic
par(mfrow = c(1, 1))
plot(sst.fit$residuals, main = "Residuals from ARIMA(2,0,2)(2,1,1)[12]")
abline(h = 0, col = "red", lty = "dashed")

par(mfrow = c(1, 2))
acf(sst.fit$residuals, xaxt = "n", xlab = "Lag(months)", main = "ACF of residuals from ARIMA(2,0,2)(2,1,1)[12]")
axis(1, at=0:804/12, labels=0:804)
pacf(sst.fit$residuals, xaxt = "n", xlab = "Lag(months)", main = "PACF of residuals from ARIMA(2,0,2)(2,1,1)[12]")
axis(1, at=0:804/12, labels=0:804)

par(mfrow = c(1, 1))
hist(sst.fit$residuals, main = "Histogram of residuals from ARIMA(2,0,2)(2,1,1)[12]", xlab = "Residuals")

# Ljung-box test check for residual independence
Box.test(sst.fit$residuals, type = "Ljung-Box")

# Forecast data
par(mfrow = c(1, 1))
sst.forecast = forecast(sst.fit, h = 60) 
plot(sst.forecast, PI = FALSE)

# Calculate accuracy with test data
accuracy(sst.forecast, sst.ts.test)

# Plot of forecasted data and test data from 1980 to 2022
autoplot(sst.forecast, PI = FALSE) + 
  xlim(1980, 2022) +
  autolayer(sst.forecast$mean, series = "Forecasted value") +
  autolayer(sst.ts.test, series = "True value") +
  ggtitle("Comparison of SARIMA forecast and true value of SST") +
  labs(y = "Temperature", x = "Year")

# Comparing data before 1950 and after 1950
c = hadisst[[500:1831]] #sst data of most recent year
extract.pts = matrix(c(llon, llat), nrow = 1, ncol = 2)
ext2 = extract(c,extract.pts, method="bilinear") #extract sst data of labuan
sst2 = ext2[1:length(ext2)]
sst.ts2 = ts(sst2, frequency = 12, start = c(1911, 8))
plot(sst.ts2, main = "Time series plot for SST")
abline(v = 1950, col = "blue")

#---------------------------------------
# Application of NNAR model on SST data
#---------------------------------------
# If running this section after running seasonal differenced GAR section,
# necessary dataset needed to be run again before proceeding 

# Fit into neural network
sst.nn = nnetar(sst.ts);sst.nn

# Experimenting with different model
m = 1
RMSE.nn = numeric(48)
MAE.nn = numeric(48)
MAPE.nn = numeric(48)
for (i in 0:5){
  for (j in 0:7){
    set.seed(1508)
    nn = nnetar(sst.ts, p = 23 + i, P = 1, size = 4 + j, period = 12)
    f = forecast(nn, h = 60)
    a = accuracy(f,sst.ts.test)
    RMSE.nn[m] = a[4]
    MAE.nn[m] = a[6]
    MAPE.nn[m] = a[10]
    m = m + 1
  }
}

# Best fit
set.seed(1508)
sst.nn.best = nnetar(sst.ts, p = 28, P = 1, size = 8, period = 12)

# Forecast
sst.nnforecast = forecast(sst.nn.best, h = 60)

# Calculate forecast accuracy
accuracy(sst.nnforecast, sst.ts.test)
checkresiduals(sst.nn.best)

# Plot of forecast data and test data from 1980 to 2022
autoplot(sst.nnforecast, PI = FALSE) + 
  xlim(1980, 2022) +
  autolayer(sst.nnforecast$mean, series = "Forecasted value") +
  autolayer(sst.ts.test, series = "True value") +
  ggtitle("Comparison of NNAR forecast and true value of SST") +
  labs(y = "Temperature", x = "Year")

#---------------------------------------
# Application of SGAR model on SST data
#---------------------------------------
# If running this section after running seasonal differenced GAR section,
# necessary dataset needed to be run again before proceeding 

# Preliminary estimation
library("hypergeo")
k = 200

ahat = function(k, s, ts){
  x = acf(ts, type = "covariance", plot = FALSE, lag.max = 450)
  ret = x$acf[k+s+1]/x$acf[k+1]
  return(ret)
}

lag0 = acf(sst.ts, type = "covariance", plot = FALSE)$acf[1] #ACVF at lag0

alpha_hat = ahat(k, s, sst.ts);alpha_hat
delta_hat = 1;delta_hat
sig2_hat = lag0/hypergeo(delta_hat, delta_hat, 1, alpha_hat^2);sig2_hat

# Minimize Whittle's estimate for sst
init = c(alpha_hat, delta_hat, sig2_hat)
whit(init)
whit_est = nlminb(init, whit, lower = c(-0.9999, 0.0001, 0.0001), upper = c(0.9999, 1, 10))
opti = c(whit_est$par[1], whit_est$par[2], whit_est$par[3]); opti

# Forecasting sst dataset
# White noise
et = ts(est_wn(n, opti[1], opti[2], 12, sst.ts-mean(sst.ts)), frequency = 12, start = c(1950,1))
plot(et, type = "l", main = "Plot of extracted white noise process for SST", ylab = expression(paste(epsilon[t])), xlab = "Year")
v = 1950+(300/12)
abline(v = v, col = "red")
et[0:299] = 0 # remove earlier estimate of inaccurate white noise

l = lm(et[300:n]~index(et)[300:n])
clip(v,2022, -50, 50)
abline(l$coefficients[[1]], l$coefficients[[2]], xlim = c(v, 2022), col = "blue")

# Forecast
set.seed(1234)
sst.forecast = ts(forecast_SGAR(60, sst.ts, opti[1], opti[2], opti[3], 12, et), frequency = 12, start = c(2017, 1))
plot(sst.forecast, main = "Seasonal GAR forecast of SST from year 2017 to year 2022", ylab = "Temperature", xlab = "Year")


# Comparison plot
autoplot(ts(sst.train, frequency = 12, start = c(1950,1))) + 
  xlim(1980, 2022) +
  autolayer(sst.forecast, series = "Forecasted value") +
  autolayer(sst.ts.test, series = "True value") +
  ggtitle("Comparison of SGAR forecast and true value of SST") +
  labs(y = "Temperature", x = "Year")


# Compare accuracy
library(forecast)
acry = accuracy(sst.forecast, sst.ts.test); acry

# BIC
BIC = 3*log(n)+whit_est$objective; BIC

#-------------------------------------------------------
# Application of SGAR model on seasonal differenced SST 
#-------------------------------------------------------
sst.ts = diff(ts(sst.train, frequency = 12, start = c(1950, 1)), lag = 12) 
n = length(sst.ts); s = 12
plot(sst.ts, main = "Plot of seasonal differenced SST", ylab = expression(X[t]), xlab = "Year")
l2 = lm(sst.ts~index(sst.ts))
abline(l2$coefficients[[1]], l2$coefficients[[2]], col = "blue")

adf = adf.test(sst.ts, alternative = "stationary", k = s); adf

# Preliminary estimation
ahat = function(k, s, ts){
  x = acf(ts, type = "covariance", plot = FALSE, lag.max = 450)
  ret = x$acf[k+s+1]/x$acf[k+1]
  if (abs(ret) > 1){
    ret = sign(ret)*0.999
  }
  return(ret)
}

lag0 = acf(sst.ts, type = "covariance", plot = FALSE)$acf[1] #ACVF at lag0

alpha_hat = ahat(k, s, sst.ts);alpha_hat
delta_hat = 1;delta_hat
sig2_hat = lag0/hypergeo(delta_hat, delta_hat, 1, alpha_hat^2);sig2_hat

# Minimize Whittle's estimate for sst
init = c(alpha_hat, delta_hat, sig2_hat)
whit(init)
whit_est = nlminb(init, whit, lower = c(-0.9999, 0.0001, 0.0001), upper = c(0.9999, 1, 10))
opti = c(whit_est$par[1], whit_est$par[2], whit_est$par[3]); opti

# Forecasting sst dataset
# White noise
et = ts(est_wn(n, opti[1], opti[2], 12, sst.ts-mean(sst.ts)), frequency = 12, start = c(1950, 1))
plot(et, type = "l", main = "Plot of extracted white noise process for seasonal differenced SST", ylab = expression(paste(epsilon[t])), xlab = "Year")
v = 1950+(300/12)
abline(v = v, col = "red")
et[0:299] = 0 # remove earlier estimate of inaccurate white noise

l = lm(et[300:n]~index(et)[300:n])
clip(v,2022, -50, 50)
abline(l$coefficients[[1]], l$coefficients[[2]], xlim = c(v, 2022), col = "blue")

# Forecast
set.seed(1234)
fct.diff = forecast_SGAR(60, sst.ts, opti[1], opti[2], opti[3], 12, et)
fct = diffinv(fct.diff, lag = 12, xi = sst.train[793:804]) #undifference forecast
sst.forecast = ts(fct[13:72], frequency = 12, start = c(2017, 1))
plot(sst.forecast, main = "Undifferenced Seasonal GAR forecast for SST data from year 2017 to year 2021", ylab = "Temperature", xlab = "Year")

# Comparison plot of forecasted value and true value
autoplot(ts(sst.train, frequency = 12, start = c(1950,1))) + 
  xlim(1980, 2022) +
  autolayer(sst.forecast, series = "Forecasted value") +
  autolayer(sst.ts.test, series = "True value") +
  ggtitle("Comparison of SGAR forecast and true value of SST") +
  labs(y = "Temperature", x = "Year")


# Compare accuracy
library(forecast)
acry2 = accuracy(sst.forecast, sst.ts.test); acry2

# BIC
BIC2 = 3*log(n)+whit_est$objective; BIC2

#----------------------------------------------------------------------------
# Comparison of Box-Jenkins, Neural Network and Seasonal GAR model estimates
#----------------------------------------------------------------------------

BIC(sst.fit)
BIC(sst.nn.best)
