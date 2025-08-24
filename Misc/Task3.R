library(dplyr)
library(survival)

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

#---------Uniform Recruitment-----------------#

#Part (a): 

### Estimating Survival Function
km_fit = survfit(Surv(OS_MONTHS,OS_STATUS)~1, data = data)

### Median Survival Time (Same as in Task 2)
median_survival_time = as.numeric(summary(km_fit)$table['median'])
median_survival_time 
#This is medium survival time after people have come into the study!
#Actual survival time == Time Until Death(if they did die) will
#have recruitment time added to this

#Hazard Rate at Median Survival Time
HZ = log(2)/median_survival_time
HZ

### Simulate Control Group Survival
n = 1e4
control_times_survival = rexp(n, HZ)
control_times_survival
#This is survival after they have been recruited

recruitment_rate = 20 #(given)
recruitment_months <- n / recruitment_rate -1
recruitment_times_control <- unlist(lapply(0:recruitment_months, function(m) runif(recruitment_rate, m, m+1)))
#These are the recruitment times of the patients

control_TOD = recruitment_times_control + control_times_survival
#This is the survival time of the controls

### Simulate Treatment Group Survival
treatment_times_survival = rexp(n, HZ * 0.7)
treatment_times_survival
#This is survival after they have been recruited

recruitment_times_treatment <- unlist(lapply(0:recruitment_months, function(m) runif(recruitment_rate, m, m+1)))
#These are the recruitment times of the patients

treatment_TOD = recruitment_times_treatment + treatment_times_survival
#This is the survival time of the treatment group i,e time till death

# Part(b): Total Follow Up Time and Power 

### Outlining the process
Follow_Up_Time = 100
control_status = ifelse(control_TOD <= Follow_Up_Time, 1, 0)
treatment_status = ifelse(treatment_TOD <= Follow_Up_Time, 1, 0)

data_combined <- data.frame(time=c(control_TOD, treatment_TOD),
                            status=c(control_status, treatment_status),
                            group=rep(c("Control", "Treatment"), each=n))
log_rank_test <- survdiff(Surv(time, status) ~ group, data = data_combined)
log_rank_test

#P-value
log_rank_test$pvalue


### Defining Functions
simu = function(FT, hr, N = 1000, HZ = log(2)/median_survival_time)
{
  Cs = rexp(N, HZ)
  Ts = rexp(N, hr * HZ)
  
  RC <- unlist(lapply(0:(N/20 -1), function(m) runif(20, m, m+1)))
  RT <- unlist(lapply(0:(N/20 -1), function(m) runif(20, m, m+1)))
  
  Ct = Cs + RC
  Tt = Ts + RT
  
  control_status = ifelse(Ct <= FT, 1, 0)
  treatment_status = ifelse(Tt <= FT, 1, 0)
  
  #print(length(Tt), length(Ct))
  df <- data.frame(time=c(Ct, Tt), status=c(control_status, treatment_status), 
                   group=rep(c("C", "T"), each=N))
  
  log_rt <- survdiff(Surv(time, status) ~ group, data = df)
  P_Val = log_rt$pvalue
  return(P_Val)
}

power = function(FT, hr)
{
  P = replicate(1000, simu(FT, hr))
  pow = sum(P<0.05)/length(P)
  return(pow)
}

X = as.integer(seq(1, 50, length.out = 20))
Y = numeric(length(X))
for(i in 1:length(X))
{
  print(i)
  Y[i] = power(X[i], 0.7)
}
#plot(X, Y)
#abline(y = 0.9)
#abline(x = 40)

# Plot Power vs Sample Size
plot(X, Y, type = "b", xlab = "Follow Up Time", ylab = "Power", 
     main = "Power vs Follow Up Time", col = "blue", pch = 19, lwd = 2)

# Add a horizontal line at 90% power
abline(h = 0.9, col ="red",lty=2)
abline(v = 39, col ="red",lty=2)

#Get Value of n where power >= 90%
index = which.max(Y >= 0.9)
X[index]

power(40, 0.7)
# A follow up time of around 40months (around 3.5 years) is enough. 

# Part(c): Hazard Ratio and Follow Up Time

###Defining Function
simu_HR = function()
{
  HR = seq(0.6, 1, length.out = 10)
  Num = numeric(length(HR))
  for(j in 1:length(HR))
  {
    print(HR[j])
    X = as.integer(seq(1, 50, length.out = 10))
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


#---------Uniform Recruitment and Dropout-----------------#
# Part(d): 

#Up until status's is same as above. Now based on follow up time we decide how many 
#people drop out. 

#Index of people that have dropped out
control_drop_index = sample(1:n, size = Follow_Up_Time, replace = FALSE)
treatment_drop_index = sample(1:n, size = Follow_Up_Time, replace = FALSE)

#Drop those indexes
control_status = control_status[-control_drop_index]
treatment_status = treatment_status[-treatment_drop_index]

control_TOD = control_TOD[-control_drop_index]
treatment_TOD = treatment_TOD[-treatment_drop_index]

#Get P-Value
data_combined2 <- data.frame(time=c(control_TOD, treatment_TOD),
                            status=c(control_status, treatment_status),
                            group=rep(c("Control", "Treatment"), each=n-Follow_Up_Time))
log_rank_test2 <- survdiff(Surv(time, status) ~ group, data = data_combined2)
log_rank_test2

#P-value
log_rank_test2$pvalue

### Defining Function
simu2 = function(FT, hr, N = 1000, HZ = log(2)/median_survival_time)
{
  Cs = rexp(N, HZ)
  Ts = rexp(N, hr * HZ)
  
  RC <- unlist(lapply(0:(N/20 -1), function(m) runif(20, m, m+1)))
  RT <- unlist(lapply(0:(N/20 -1), function(m) runif(20, m, m+1)))
  
  Ct = Cs + RC
  Tt = Ts + RT
  
  cont_stat = ifelse(Ct <= FT, 1, 0)
  treat_stat = ifelse(Tt <= FT, 1, 0)
  
  CI= sample(1:N, size = FT, replace = FALSE)
  TI = sample(1:N, size = FT, replace = FALSE)
  
  cont_stat = cont_stat[-CI]
  treat_stat = treat_stat[-TI]
  
  Ct = Ct[-CI]
  Tt = Tt[-TI]
  
  #print(length(Tt), length(Ct))
  df <- data.frame(time=c(Ct, Tt), status=c(cont_stat, treat_stat), 
                   group=rep(c("C", "T"), each=N-FT))
  
  log_rt <- survdiff(Surv(time, status) ~ group, data = df)
  P_Val = log_rt$pvalue
  return(P_Val)
}

power2 = function(FT, hr)
{
  P = replicate(1000, simu2(FT, hr))
  pow = sum(P<0.05)/length(P)
  return(pow)
}

simu_HR2 = function()
{
  HR = seq(0.6, 1, length.out = 10)
  Num = numeric(length(HR))
  for(j in 1:length(HR))
  {
    print(HR[j])
    X = as.integer(seq(1, 50, length.out = 10))
    Y = numeric(length(X))
    for(i in 1:length(X))
    {
      print(i)
      Y[i] = power2(X[i], HR[j])
    }
    index = which.max(Y >= 0.9)
    Num[j] = X[index]
    print(X[index])
  }
  return(list(HR, Num))
}
out2 = simu_HR2()
out2
