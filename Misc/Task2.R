library(dplyr)
library(survival)
library(ggplot2)

#---------Loading Data--------------------#
data_patients = read.csv("../project5-jandj-project5_group19/Data/msk_chord_2024/data_clinical_patient.txt", 
                 header = TRUE, sep = "\t", skip = 4)
data_samples = read.csv("../project5-jandj-project5_group19/Data/msk_chord_2024/data_clinical_sample.txt", 
                        header = TRUE, sep = "\t", skip = 4)

#Getting patient ids of only lung cancer patients
required_patients = unique(data_samples$PATIENT_ID[data_samples$CANCER_TYPE == "Non-Small Cell Lung Cancer"])

#getting all the final required data
data = data_patients %>% filter(PATIENT_ID %in% required_patients) %>% select(PATIENT_ID, OS_MONTHS, OS_STATUS)
data

#Doing Some Pre-Processing
data$OS_STATUS = gsub(":.*", "", data$OS_STATUS)
data$OS_MONTHS = as.numeric(data$OS_MONTHS)
data$OS_STATUS = as.numeric(data$OS_STATUS)
data

#-------Estimating Survival Function-----------#
km_fit = survfit(Surv(OS_MONTHS,OS_STATUS)~1, data = data)

# Plot Kaplan-Meier survival curve
plot(km_fit, 
     main = "Kaplan-Meier Survival Curve", 
     xlab = "Time", 
     ylab = "Survival Probability")

#Part (a): Median Survival Time
median_survival_time = as.numeric(summary(km_fit)$table['median'])
median_survival_time

#Part (b): Hazard Rate at Median Survival Time
HZ = log(2)/median_survival_time
HZ

#Part (c): Simulate Control Group Survival
n = 1e4
control_times = rexp(n, HZ)
control_times
#Here we assume they died at this time

#Part (D): Simulate from Treatment Group
treatment_times = rexp(n, 0.7 * HZ)
treatment_times

#Part (E):log-rank Test for H0: HZ < 1
data_combined <- data.frame(time=c(control_times, treatment_times),
                          status=rep(1, 2*n),
                          group=rep(c("Control", "Treatment"), each=n))
#Here we have assume the trial duration to be longer then maximum survival time.
#This basically means that the trial duration is long enough that by the end 
#Everyone is deceased. Hence we have 1 eveywhere. 

log_rank_test <- survdiff(Surv(time, status) ~ group, data = data_combined)
log_rank_test

#P-Value is very small -> This implies that there is a difference between the distributions
log_rank_test$pvalue

#Checking with plot
km_fit_comb = survfit(Surv(time, status) ~ group, data = data_combined)
km_fit_comb

plot(km_fit_comb, 
     main = "Kaplan-Meier Survival Curve", 
     xlab = "Time", 
     ylab = "Survival Probability",
     col = c("blue", "red"),
     fun = "surv")

legend("topright", 
       legend = names(km_fit_comb$strata), 
       lty = 1, 
       col = c("blue", "red"), 
       bty = "n")

# We can see that the treatment group has a better survival probability. Using this + the fact that 
# the p-value calcaulated above is less thes zero i.e. both groups are different, we can say that the
# hazard ratio is less then 1

#Defining Functions for the next part. 
simu = function(N, hr, HZ = log(2)/median_survival_time)
{
  Ct = rexp(N, HZ)
  Tt = rexp(N, hr * HZ)
  #print(length(Tt), length(Ct))
  df <- data.frame(time=c(Ct, Tt), status=rep(1, 2*N), group=rep(c("C", "T"), each=N))
  
  log_rt <- survdiff(Surv(time, status) ~ group, data = df)
  P_Val = log_rt$pvalue
  return(P_Val)
}

power = function(N, hr)
{
  P = replicate(1000, simu(N, hr))
  pow = sum(P<0.05)/length(P)
  return(pow)
}


#Part(F)
power(100, 0.7)

#Part(G)
X = as.integer(10^(seq(1, 3, length.out = 20)))
Y = numeric(length(X))
for(i in 1:length(X))
{
  Y[i] = power(X[i], 0.7)
}
plot(X, Y)

#Get Value of n where power >= 90%
index = which.max(Y >= 0.9)
X[index]

#Part(H)
simu_HR = function()
{
  HR = seq(0.6, 1, length.out = 10)
  Num = numeric(length(HR))
  for(j in 1:length(HR))
  {
    print(HR[j])
    X = as.integer(10^(seq(0, 3.5, length.out = 10)))
    Y = numeric(length(X))
    for(i in 1:length(X))
    {
      print(i)
      Y[i] = power(X[i], HR[j])
    }
    index = which.max(Y >= 0.9)
    Num[j] = X[index]
    print(X[index])
  }
  return(list(HR, Num))
}
out = simu_HR() 
out
