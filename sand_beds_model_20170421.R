# A new model which concentrates on estimating only those sand bed ages which are unknown based on only the surrounding radiocarbon dates
# Updated 20170421 to change the calibration curve and use 


# The details of the model are:
# There are 8 sand beds with different dates surrounding them
# The first 5 sand beds have known ages/locations:
# Sand bed 1 - AD 2010 - northern rupture
# Sand bed 2 - AD 1960 - southern rupture
# Sand bed 3 - AD 1835 - northern rupture
# Sand bed 4- AD 1751 - northern rupture
# Sand bed 5- AD 1575 - southern rupture
# The bottom three sand beds are unknown ages but known location
# Sand bed 6 - prehistoric (ages in attached table) - northern rupture
# Sand bed 7 - prehistoric (ages in attached table) - southern rupture
# Sand bed 8 - prehistoric (ages in attached table) - northern rupture
# Basal unit - base of section 

# The questions are:
# 1 - what are the ages of the 8 events? (first 5 known)
# 2 - what is the reccurrence interval between events?
# 3 - what is the recurrence interval between northern and southern earthquakes?

# Set up ------------------------------------------------------------------

# Clear the workspace if necessary
rm(list=ls())

# Load in packages
pkgs = c('tidyverse','Bchron')
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

# Data carpentry ----------------------------------------------------------

# Load in the calibration curves
cal = read.table('shcal13.14c',sep=',')
colnames(cal) = c('calbp','c14age','errorBP','Delta14C','Sigma_per_mil') 

# SB6 ---------------------------------------------------------------------

# Sand bed 6 is the first one that isn't known. It has two restrictions, an upper and a lower date. They are:
x_SB6_max = 420  
sig_SB6_max = 20
x_SB6_min = 245  
sig_SB6_min = 55

# If we calibrate them and plot we get
SB6_min_cal = BchronCalibrate(x_SB6_min,
                              sig_SB6_min,
                              calCurves='shcal13') # about 0 to 450
SB6_max_cal = BchronCalibrate(x_SB6_max,
                              sig_SB6_max,
                              calCurves='shcal13') # about 500

par(mfrow=c(2,1))
plot(SB6_min_cal)
plot(SB6_max_cal)
par(mfrow=c(1,1))

# Must be older than SB6 so theta_SB6 must be greater than 376
xlim_range = c(0, 600)
par(mfrow=c(2,1))
plot(SB6_min_cal, xlim = xlim_range)
plot(SB6_max_cal, xlim = xlim_range)
par(mfrow=c(1,1))

# Consider a model for these data you have:
# x_1 ~ N(r(theta_1), sigma_1^2)
# x_2 ~ N(r(theta_2), sigma_2^2)
# theta_1 ~ U(0, 600)
# theta_2 ~ U(0, 600)
# theta_SB6 ~ U(min(theta_1,376), theta_2)

# Work on a grid (say 10 years?) would give ~216000 possibilities 
grid_size = 5
theta_1_grid = theta_2_grid = theta_SB6_grid = seq(xlim_range[1], xlim_range[2], by = grid_size)
theta_all = expand.grid(theta_1_grid, theta_2_grid, theta_SB6_grid)
colnames(theta_all) = c('theta_1', 'theta_2', 'theta_SB6')
good_rows1 = with(theta_all, theta_SB6 < theta_2)
good_rows2 = with(theta_all, theta_SB6 > theta_1)
good_rows3 = with(theta_all, theta_SB6 > 376)
theta_use = theta_all[good_rows1 & good_rows2 & good_rows3,]

# Now loop through each row and score on the log likelihood
post = rep(NA, nrow(theta_use))
print("SB6")
for(i in 1:nrow(theta_use)) {
  #print(nrow(theta_use) - i)
  curr_mean_min = approx(cal$calbp,cal$c14,theta_use[i,'theta_1'])$y
  curr_sd_min = approx(cal$calbp,cal$errorBP,theta_use[i,'theta_1'])$y
  curr_mean_max = approx(cal$calbp,cal$c14,theta_use[i,'theta_2'])$y
  curr_sd_max = approx(cal$calbp,cal$errorBP,theta_use[i,'theta_2'])$y
  post[i] = dnorm(x_SB6_min, 
                     mean = curr_mean_min, 
                     sd = sqrt(curr_sd_min^2 + sig_SB6_min^2)) *
                     #log = TRUE) +
    dnorm(x_SB6_max, 
          mean = curr_mean_max, 
          sd = sqrt(curr_sd_max^2 + sig_SB6_max^2))
          #log = TRUE)
}

# Add up the posts for each unique value of theta_SB6
theta_SB6_final = aggregate(post, by = list(theta_use[,'theta_SB6']), sum)
theta_SB6_dens = theta_SB6_final[,2]/sum(theta_SB6_final[,2])

pdf(file='SB6_20170421.pdf',width=8,height=8)
par(mfrow=c(3,1))
par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01,las=1)
plot(SB6_min_cal, xlim = xlim_range)
plot(theta_SB6_final[,1], theta_SB6_dens, type = 'l', xlab = 'SB6', xlim = xlim_range)
plot(SB6_max_cal, xlim = xlim_range)
par(mfrow=c(1,1))
dev.off()

# SB7 ---------------------------------------------------------------------

# Sand bed 7 is the next one. It has 8(!) restrictions, an upper and 4 lower date. They are:
x_SB7_max = 460 # Quidico - OS-106356  
sig_SB7_max = 30
x_SB7_min = 460 # OS-102372
sig_SB7_min = 20
x_SB7_min_2 = 1630 # OS-111756    
sig_SB7_min_2 = 20
x_SB7_min_3 = 570 # Quidico - OS-103409
sig_SB7_min_3 = 45
x_SB7_min_4 = 495 # Quidico - OS-115318
sig_SB7_min_4 = 15
x_SB7_max_2 = 425 #425 +/- 15 (above basal sand) - OS-106268???
sig_SB7_max_2 = 15
x_SB7_max_3 = 475 #475 +/- 20 (above basal sand) - OS-103174???
sig_SB7_max_3 = 20
x_SB7_min_5 = 760 #760 +/- 20 (within basal sand) - OS-103173???
sig_SB7_min_5 = 20

# If we calibrate them and plot we get
SB7_min_cal = BchronCalibrate(x_SB7_min,
                              sig_SB7_min,
                              calCurves='shcal13')
SB7_min_2_cal = BchronCalibrate(x_SB7_min_2,
                                sig_SB7_min_2,
                                calCurves='shcal13')
SB7_min_3_cal = BchronCalibrate(x_SB7_min_3,
                                sig_SB7_min_3,
                                calCurves='shcal13')
SB7_min_4_cal = BchronCalibrate(x_SB7_min_4,
                                sig_SB7_min_4,
                                calCurves='shcal13')
SB7_min_5_cal = BchronCalibrate(x_SB7_min_5,
                                sig_SB7_min_5,
                                calCurves='shcal13')
SB7_max_cal = BchronCalibrate(x_SB7_max,
                              sig_SB7_max,
                              calCurves='shcal13')
SB7_max_2_cal = BchronCalibrate(x_SB7_max_2,
                                sig_SB7_max_2,
                                calCurves='shcal13')
SB7_max_3_cal = BchronCalibrate(x_SB7_max_3,
                                sig_SB7_max_3,
                                calCurves='shcal13')

pdf(file = 'SB7_plots.pdf', width = 12, height = 8)
par(mfrow=c(8,1))
par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01,las=1)
plot(SB7_min_cal, xlim = c(200, 1600), main = 'Tirua - OS-102372 - maximum', xlab = 'Age (years BP)')
plot(SB7_min_2_cal, xlim = c(200, 1600), main = 'Tirua - OS-111756 - maximum', xlab = 'Age (years BP)')
plot(SB7_min_3_cal, xlim = c(200, 1600), main = 'Quidico - OS-103409 - maximum', xlab = 'Age (years BP)')
plot(SB7_min_4_cal, xlim = c(200, 1600), main = 'Quidico - OS-115318 - maximum', xlab = 'Age (years BP)')
plot(SB7_min_5_cal, xlim = c(200, 1600), main = 'Quidico - OS-103173 - maximum/within? (new)', xlab = 'Age (years BP)')
plot(SB7_max_cal, xlim = c(200, 1600), main = 'Quidico - OS-106356 - minimum', xlab = 'Age (years BP)')
plot(SB7_max_2_cal, xlim = c(200, 1600), main = 'Quidico - OS-106268 - minimum (new)', xlab = 'Age (years BP)')
plot(SB7_max_3_cal, xlim = c(200, 1600), main = 'Quidico - OS-103174 - minimum (new)', xlab = 'Age (years BP)')
par(mfrow=c(1,1))
dev.off()

# You might as well drop min_2 as it provides no useful bound

# Consider a model for these data you have:
# x_min_1 ~ N(r(theta_1), sigma_1^2)
# x_min_3 ~ N(r(theta_1), sigma_3^2)
# x_min_4 ~ N(r(theta_1), sigma_4^2)
# x_2 ~ N(r(theta_2), sigma_2^2)
# theta_1 ~ U(400, 700)
# theta_2 ~ U(400, 700)
# theta_SB7 ~ U(theta_1, theta_2)

# Work on a grid (say 10 years?) 
grid_size = 5
theta_1_grid = theta_2_grid = theta_SB7_grid = seq(400, 700, by = grid_size)
theta_all = expand.grid(theta_1_grid, theta_2_grid, theta_SB7_grid)
colnames(theta_all) = c('theta_1', 'theta_2', 'theta_SB7')
good_rows1 = with(theta_all, theta_SB7 < theta_2)
good_rows2 = with(theta_all, theta_SB7 > theta_1)
theta_use = theta_all[good_rows1 & good_rows2,]

# Now loop through each row and score on the log likelihood
post = rep(NA, nrow(theta_use))
print("SB7")
for(i in 1:nrow(theta_use)) {
  #print(nrow(theta_use) - i)
  curr_mean_min = approx(cal$calbp,cal$c14,theta_use[i,'theta_1'])$y
  curr_sd_min = approx(cal$calbp,cal$errorBP,theta_use[i,'theta_1'])$y
  curr_mean_max = approx(cal$calbp,cal$c14,theta_use[i,'theta_2'])$y
  curr_sd_max = approx(cal$calbp,cal$errorBP,theta_use[i,'theta_2'])$y
  post[i] = dnorm(x_SB7_min, 
                     mean = curr_mean_min, 
                     sd = sqrt(curr_sd_min^2 + sig_SB7_min^2)) *
    dnorm(x_SB7_min_3, 
          mean = curr_mean_min, 
          sd = sqrt(curr_sd_min^2 + sig_SB7_min_3^2)) *
    dnorm(x_SB7_min_4, 
          mean = curr_mean_min, 
          sd = sqrt(curr_sd_min^2 + sig_SB7_min_4^2)) *
    dnorm(x_SB7_max, 
          mean = curr_mean_max, 
          sd = sqrt(curr_sd_max^2 + sig_SB7_max^2))
  #log = TRUE)
}

# Add up the posts for each unique value of theta_3
theta_SB7_final = aggregate(post, by = list(theta_use[,'theta_SB7']), sum)
theta_SB7_dens = theta_SB7_final[,2]/sum(theta_SB7_final[,2])

xlim_range = c(400,700)
pdf(file='SB7_20170421.pdf',width=8,height=8)
par(mfrow=c(6,1))
par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01,las=1)
plot(SB7_min_cal, xlim = xlim_range)
plot(SB7_min_2_cal, xlim = xlim_range)
plot(SB7_min_3_cal, xlim = xlim_range)
plot(SB7_min_4_cal, xlim = xlim_range)
plot(theta_SB7_final[,1], theta_SB7_dens, type = 'l', xlab = 'SB7', xlim = xlim_range)
plot(SB7_max_cal, xlim = xlim_range)
par(mfrow=c(1,1))

# SB8 ---------------------------------------------------------------------

# Sand bed 8 is next. It has two restrictions, an upper and a lower date. They are:
x_SB8_max = 1630  
sig_SB8_max = 20
x_SB8_min = 1790  
sig_SB8_min = 20

# If we calibrate them and plot we get
SB8_min_cal = BchronCalibrate(x_SB8_min,
                              sig_SB8_min,
                              calCurves='shcal13') # about 0 to 450
SB8_max_cal = BchronCalibrate(x_SB8_max,
                              sig_SB8_max,
                              calCurves='shcal13') # about 500

xlim_range = c(1400, 1850)
par(mfrow=c(2,1))
plot(SB8_min_cal, xlim = xlim_range)
plot(SB8_max_cal, xlim = xlim_range)
par(mfrow=c(1,1))

# Consider a model for these data you have:
# x_1 ~ N(r(theta_1), sigma_1^2)
# x_2 ~ N(r(theta_2), sigma_2^2)
# theta_1 ~ U(1400, 1850)
# theta_2 ~ U(1400, 1850)
# theta_SB8 ~ U(theta_1, theta_2)

# Work on a grid (say 10 years?)
grid_size = 5
theta_1_grid = theta_2_grid = theta_SB8_grid = seq(xlim_range[1], xlim_range[2], by = grid_size)
theta_all = expand.grid(theta_1_grid, theta_2_grid, theta_SB8_grid)
colnames(theta_all) = c('theta_1', 'theta_2', 'theta_SB8')
good_rows1 = with(theta_all, theta_SB8 < theta_2)
good_rows2 = with(theta_all, theta_SB8 > theta_1)
theta_use = theta_all[good_rows1 & good_rows2,]

# Now loop through each row and score on the log likelihood
post = rep(NA, nrow(theta_use))
print("SB8")
for(i in 1:nrow(theta_use)) {
  #print(nrow(theta_use) - i)
  curr_mean_min = approx(cal$calbp,cal$c14,theta_use[i,'theta_1'])$y
  curr_sd_min = approx(cal$calbp,cal$errorBP,theta_use[i,'theta_1'])$y
  curr_mean_max = approx(cal$calbp,cal$c14,theta_use[i,'theta_2'])$y
  curr_sd_max = approx(cal$calbp,cal$errorBP,theta_use[i,'theta_2'])$y
  post[i] = dnorm(x_SB8_min, 
                  mean = curr_mean_min, 
                  sd = sqrt(curr_sd_min^2 + sig_SB8_min^2)) *
    dnorm(x_SB8_max, 
          mean = curr_mean_max, 
          sd = sqrt(curr_sd_max^2 + sig_SB8_max^2))
}

# Add up the posts for each unique value of theta_SB8
theta_SB8_final = aggregate(post, by = list(theta_use[,'theta_SB8']), sum)
theta_SB8_dens = theta_SB8_final[,2]/sum(theta_SB8_final[,2])

pdf(file='SB8_20170421.pdf',width=8,height=8)
par(mfrow=c(3,1))
par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01,las=1)
plot(SB8_min_cal, xlim = xlim_range)
plot(theta_SB8_final[,1], theta_SB8_dens, type = 'l', xlab = 'SB8', xlim = xlim_range)
plot(SB8_max_cal, xlim = xlim_range)
par(mfrow=c(1,1))
dev.off()

# basal ---------------------------------------------------------------------

# The basal slice ('SB9' but not really) is next. It has two restrictions, an upper and a lower date. They are:
x_basal_max = 1760  
sig_basal_max = 25
x_basal_min = 1850  
sig_basal_min = 20

# If we calibrate them and plot we get
basal_min_cal = BchronCalibrate(x_basal_min,
                              sig_basal_min,
                              calCurves='shcal13') # about 0 to 450
basal_max_cal = BchronCalibrate(x_basal_max,
                              sig_basal_max,
                              calCurves='shcal13') # about 500

xlim_range = c(1550, 1900)
par(mfrow=c(2,1))
plot(basal_min_cal, xlim = xlim_range)
plot(basal_max_cal, xlim = xlim_range)
par(mfrow=c(1,1))

# Consider a model for these data you have:
# x_1 ~ N(r(theta_1), sigma_1^2)
# x_2 ~ N(r(theta_2), sigma_2^2)
# theta_1 ~ U(1550, 1900)
# theta_2 ~ U(1550, 1900)
# theta_basal ~ U(theta_1, theta_2)

# Work on a grid (say 10 years?)
grid_size = 5
theta_1_grid = theta_2_grid = theta_basal_grid = seq(xlim_range[1], xlim_range[2], by = grid_size)
theta_all = expand.grid(theta_1_grid, theta_2_grid, theta_basal_grid)
colnames(theta_all) = c('theta_1', 'theta_2', 'theta_basal')
good_rows1 = with(theta_all, theta_basal < theta_2)
good_rows2 = with(theta_all, theta_basal > theta_1)
theta_use = theta_all[good_rows1 & good_rows2,]

# Now loop through each row and score on the log likelihood
post = rep(NA, nrow(theta_use))
print("Basal")
for(i in 1:nrow(theta_use)) {
  #print(nrow(theta_use) - i)
  curr_mean_min = approx(cal$calbp,cal$c14,theta_use[i,'theta_1'])$y
  curr_sd_min = approx(cal$calbp,cal$errorBP,theta_use[i,'theta_1'])$y
  curr_mean_max = approx(cal$calbp,cal$c14,theta_use[i,'theta_2'])$y
  curr_sd_max = approx(cal$calbp,cal$errorBP,theta_use[i,'theta_2'])$y
  post[i] = dnorm(x_basal_min, 
                  mean = curr_mean_min, 
                  sd = sqrt(curr_sd_min^2 + sig_basal_min^2)) *
    dnorm(x_basal_max, 
          mean = curr_mean_max, 
          sd = sqrt(curr_sd_max^2 + sig_basal_max^2))
}

# Add up the posts for each unique value of theta_basal
theta_basal_final = aggregate(post, by = list(theta_use[,'theta_basal']), sum)
theta_basal_dens = theta_basal_final[,2]/sum(theta_basal_final[,2])

pdf(file='basal_20170421.pdf',width=8,height=8)
par(mfrow=c(3,1))
par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01,las=1)
plot(basal_min_cal, xlim = xlim_range)
plot(theta_basal_final[,1], theta_basal_dens, type = 'l', xlab = 'basal', xlim = xlim_range)
plot(basal_max_cal, xlim = xlim_range)
par(mfrow=c(1,1))
dev.off()

# Create nice plot --------------------------------------------------------

# So we now have dates for all layers - need to create a nice plot for them

# First set of sand beds are known
# Sand bed 1 - AD 2010 - northern rupture
# Sand bed 2 - AD 1960 - southern rupture
# Sand bed 3 - AD 1835 - northern rupture
# Sand bed 4- AD 1751 - northern rupture
# Sand bed 5- AD 1575 - southern rupture

# Create a list of all of them together in BC/AD format
SB_ages = c(2010, 1960, 1835, 1751, 1575)
SB_ages_shift = sort(c(SB_ages-2, SB_ages-1, SB_ages, SB_ages + 1, SB_ages+2), decreasing = TRUE)
df = data.frame(ages = c(SB_ages_shift,
                         c(1950-theta_SB6_final[1,1] - 1, 1950-theta_SB6_final[,1], 1950-theta_SB6_final[nrow(theta_SB6_final),1] + 1), 
                         c(1950-theta_SB7_final[1,1] - 1, 1950-theta_SB7_final[,1], 1950-theta_SB7_final[nrow(theta_SB7_final),1] + 1), 
                         c(1950-theta_SB8_final[1,1] - 1, 1950-theta_SB8_final[,1], 1950-theta_SB8_final[nrow(theta_SB8_final),1] + 1), 
                         c(1950-theta_basal_final[1,1] - 1, 1950-theta_basal_final[,1], 1950-theta_basal_final[nrow(theta_basal_final),1] + 1)), 
                Densities = c(rep(c(0, 1, 1, 1, 0), 5),
                              c(0, theta_SB6_dens, 0),
                              c(0, theta_SB7_dens, 0),
                              c(0, theta_SB8_dens, 0),
                              c(0, theta_basal_dens, 0)),
                labels = c(rep(c('SB1','SB2','SB3','SB4','SB5'), each = 5),
                           rep('SB6', length(theta_SB6_dens)+2),
                           rep('SB7', length(theta_SB7_dens)+2),
                           rep('SB8', length(theta_SB8_dens)+2),
                           rep('Basal', length(theta_basal_dens)+2)))

df$labels = factor(df$labels, 
                   levels = c('SB1','SB2','SB3','SB4','SB5','SB6','SB7','SB8','Basal'),
                   ordered = TRUE)

# save the object
saveRDS(df, file = 'SB_posteriors.rds')
df = readRDS(file = 'SB_posteriors.rds')

#Create plot 
ggplot(df, aes(x = ages,y = Densities, fill=labels, group = labels)) + 
  geom_polygon() + 
  facet_grid(labels ~ .,scales='free') + 
  ylab('Probability Density') +
  theme_bw()+
  theme(legend.position='None', axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        #axis.title.y=element_text(angle=0,vjust=1,hjust=0),+
  ggtitle('Age of sand beds')+xlab('Age (years AD)')
ggsave('Sand_beds_fig1.pdf',height=8,width=8)


# Recurrence relations ----------------------------------------------------

# First difference between sand bed ages

# Simulate ages for each SB 
n_sim = 10000
all_sim = matrix(NA, ncol = 9, nrow = n_sim)
colnames(all_sim) = c(paste0('SB',1:8), 'Basal')
for(i in 1:9) {
  # Get the current densities
  curr_dens = df[df$labels == levels(df$labels)[i],]
  all_sim[,i] = with(curr_dens, sample(curr_dens$ages, size = n_sim, replace = TRUE, prob = curr_dens$Densities))
}

# Now create differences
all_sim_diff = matrix(NA, ncol = 8, nrow = n_sim)
for(i in 1:ncol(all_sim_diff)) {
  all_sim_diff[,i] = all_sim[,i] - all_sim[,i+1]
}
df2 = data.frame(all_sim_diff) 
colnames(df2) = c('SB2 - SB1',
                  'SB3 - SB2',
                  'SB4 - SB3',
                  'SB5 - SB4',
                  'SB6 - SB5',
                  'SB7 - SB6',
                  'SB8 - SB7',
                  'Basal - SB8')
df2a = gather(df2)
df2a$key = factor(df2a$key, 
                   levels = c('SB2 - SB1',
                              'SB3 - SB2',
                              'SB4 - SB3',
                              'SB5 - SB4',
                              'SB6 - SB5',
                              'SB7 - SB6',
                              'SB8 - SB7',
                              'Basal - SB8'),
                   ordered = TRUE)

# Now plot the differences
ggplot(df2a, aes(x = value , fill=key, group = key)) + 
  geom_density(colour = NA) + 
  facet_grid(key ~ .,scales='free') + 
  ylab('Probability Density') +
  theme_bw()+
  theme(legend.position='None', axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  #axis.title.y=element_text(angle=0,vjust=1,hjust=0),+
  ggtitle('Age difference of sand beds')+xlab('Age difference (years AD)')
ggsave('Sand_beds_fig2.pdf',height=8,width=8)

# Recurrence within segment -----------------------------------------------

# Next need reurrences within north/south
df$segment = 'Northern'
df$segment[df$labels == 'SB2'] = 'Southern'
df$segment[df$labels == 'SB5'] = 'Southern'
df$segment[df$labels == 'SB7'] = 'Southern'
df$segment[df$labels == 'Basal'] = 'Neither'
table(df$label, df$segment)

# How many of each
n_north = 5
n_south = 3

# Create separate data frames
df_north = df[df$segment == 'Northern', ]
df_south = df[df$segment == 'Southern', ]

# Now create differences, samples, and plot
n_sim = 10000
all_sim_north = matrix(NA, ncol = n_north, nrow = n_sim)
for(i in 1:n_north) {
  # Get the current densities
  curr_dens = df_north[df_north$labels == unique(df_north$labels)[i],]
  all_sim_north[,i] = with(curr_dens, sample(curr_dens$ages, size = n_sim, replace = TRUE, prob = curr_dens$Densities))
}

# Now create differences
all_sim_diff_north = matrix(NA, ncol = ncol(all_sim_north)-1, nrow = n_sim)
for(i in 1:ncol(all_sim_diff_north)) {
  all_sim_diff_north[,i] = all_sim_north[,i] - all_sim_north[,i+1]
}
df2_north = data.frame(all_sim_diff_north) 
colnames(df2_north) = c('SB3 - SB1',
                  'SB4 - SB3',
                  'SB6 - SB4',
                  'SB8 - SB6')
df2a_north = gather(df2_north)
df2a_north$key = factor(df2a_north$key, 
                  levels = c('SB3 - SB1',
                             'SB4 - SB3',
                             'SB6 - SB4',
                             'SB8 - SB6'),
                  ordered = TRUE)

# Now plot the differences
ggplot(df2a_north, aes(x = value , fill=key, group = key)) + 
  geom_density(colour = NA) + 
  facet_grid(key ~ .,scales='free') + 
  ylab('Probability Density') +
  theme_bw()+
  theme(legend.position='None', axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  #axis.title.y=element_text(angle=0,vjust=1,hjust=0),+
  ggtitle('Age difference of sand beds - northern segment')+xlab('Age difference (years AD)')
ggsave('Sand_beds_fig3.pdf',height=8,width=8)


# Same again for south ----------------------------------------------------

# Now create differences, samples, and plot
n_sim = 10000
all_sim_south = matrix(NA, ncol = n_south, nrow = n_sim)
for(i in 1:n_south) {
  # Get the current densities
  curr_dens = df_south[df_south$labels == unique(df_south$labels)[i],]
  all_sim_south[,i] = with(curr_dens, sample(curr_dens$ages, size = n_sim, replace = TRUE, prob = curr_dens$Densities))
}

# Now create differences
all_sim_diff_south = matrix(NA, ncol = ncol(all_sim_south)-1, nrow = n_sim)
for(i in 1:ncol(all_sim_diff_south)) {
  all_sim_diff_south[,i] = all_sim_south[,i] - all_sim_south[,i+1]
}
df2_south = data.frame(all_sim_diff_south) 
colnames(df2_south) = c('SB5 - SB2',
                        'SB7 - SB5')
df2a_south = gather(df2_south)
df2a_south$key = factor(df2a_south$key, 
                        levels = c('SB5 - SB2',
                                   'SB7 - SB5'),
                        ordered = TRUE)

# Now plot the differences
ggplot(df2a_south, aes(x = value , fill=key, group = key)) + 
  geom_density(colour = NA) + 
  facet_grid(key ~ .,scales='free') + 
  ylab('Probability Density') +
  theme_bw()+
  theme(legend.position='None', axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  #axis.title.y=element_text(angle=0,vjust=1,hjust=0),+
  ggtitle('Age difference of sand beds - southern segment')+xlab('Age difference (years AD)')
ggsave('Sand_beds_fig4.pdf',height=8,width=8)


# Create tabular output ---------------------------------------------------

# First based on all_sim
ans1 = apply(all_sim, 2, quantile, probs = c(0.025, 0.25, 0.50, 0.75, 0.975))

# Now on all diffs
ans2 = apply(df2, 2, quantile, probs = c(0.025, 0.25, 0.50, 0.75, 0.975))

# Now on north diffs
ans3 = apply(df2_north, 2, quantile, probs = c(0.025, 0.25, 0.50, 0.75, 0.975))

# Now on south diffs
ans4 = apply(df2_south, 2, quantile, probs = c(0.025, 0.25, 0.50, 0.75, 0.975))

write.csv(ans1, file = 'SB_ages.csv')
write.csv(ans2, file = 'SB_age_diffs.csv')
write.csv(ans3, file = 'SB_age_diffs_north.csv')
write.csv(ans4, file = 'SB_age_diffs_south.csv')
