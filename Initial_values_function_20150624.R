# Function to produce good initial values for Ben_Tsunami_20150612 jags code

# The idea is to first sample a reasonable age for T12 (theta[1]), then T11, etc

init= function() {
  theta = vector(length=14)
  
  # Choose one of the dates for T12
  choose_T12 = sample(length(x_T12),1)
  T12_cal = BchronCalibrate(x_T12[choose_T12],sig_T12[choose_T12],calCurves='intcal13')
  theta[1] = sample(T12_cal$date1$ageGrid,1,prob=T12_cal$date1$densities)
  
  # Same again for T11
  choose_T11 = sample(length(x_T11),1)
  T11_cal = BchronCalibrate(x_T11[choose_T11],sig_T11[choose_T11],calCurves='intcal13')
  theta[2] = sample(T11_cal$date1$ageGrid,1,prob=T11_cal$date1$densities)
  
  # T10
  choose_T10 = sample(length(x_T10),1)
  T10_cal = BchronCalibrate(x_T10[choose_T10],sig_T10[choose_T10],calCurves='intcal13')
  theta[3] = sample(T10_cal$date1$ageGrid,1,prob=T10_cal$date1$densities)
  
  # T19
  choose_T9 = sample(length(x_T9),1)
  T9_cal = BchronCalibrate(x_T9[choose_T9],sig_T9[choose_T9],calCurves='intcal13')
  theta[4] = sample(T9_cal$date1$ageGrid,1,prob=T9_cal$date1$densities)
  
  # T8
  choose_T8 = sample(length(x_T8),1)
  T8_cal = BchronCalibrate(x_T8[choose_T8],sig_T8[choose_T8],calCurves='intcal13')
  theta[5] = sample(T8_cal$date1$ageGrid,1,prob=T8_cal$date1$densities)
    
  # No direct dates for T6 (theta[6]), T6 (theta[7]), T5 (theta[8]), or T4 ((theta[9]) so simulate theta[10] (between T3 and T4)
  # and then uniformly simulate between them
  
  # T4_min_max
  choose_T4_min_max = sample(length(x_T4_min_max),1)
  T4_min_max_cal = BchronCalibrate(x_T4_min_max[choose_T4_min_max],sig_T4_min_max[choose_T4_min_max],calCurves='intcal13')
  theta[10] = sample(T4_min_max_cal$date1$ageGrid,1,prob=T4_min_max_cal$date1$densities)
  
  # Now simulate T7, T6, T5, T4
  theta[6:9] = round(runif(4,theta[5],theta[10]),0)
  
  # No dates for T3 (theta[11]) or T2 (theta[12]) so simulate theta[13] (between T1 and T2) and then fill in
  choose_T2_min_max = sample(length(x_T2_min_max),1)
  T2_min_max_cal = BchronCalibrate(x_T2_min_max[choose_T2_min_max],sig_T2_min_max[choose_T2_min_max],calCurves='intcal13')
  theta[13] = sample(T2_min_max_cal$date1$ageGrid,1,prob=T2_min_max_cal$date1$densities)
  
  # Now simulate T3 and T2
  theta[11:12] = round(runif(2,theta[10],theta[13]),0)
  
  # No direct ages for T1 but there are two max ages so use the lowest of these
  T1_max_cal = BchronCalibrate(x_T1_max_1,sig_T1_max_1,calCurves='intcal13')
  theta[14] = round(runif(1,theta[13],sample(T1_max_cal$date1$ageGrid,1,prob=T1_max_cal$date1$densities)),0)
  
  return(list(theta_raw=sort(theta)))
}
  
  
  
  
  
