# A JAGS model for a constrained radiocarbon chronology for some data 
# from two cores - Tirúa and Quidico

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
pkgs = c('rjags','ggplot2','Bchron')
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

# Data carpentry ----------------------------------------------------------

# Load in the calibration curves
cal = read.table('intcal13.14c',sep=',')
colnames(cal) = c('calbp','c14age','errorBP','Delta14C','Sigma_per_mil') 

# Now write out the dates - all for core Tirúa except otherwise stated
x_SB2_min = 65 # Quidico
sig_SB2_min = 30
x_SB3_min = 105 # Quidico  
sig_SB3_min = 25
x_SB4_min = 360  # Quidico
sig_SB4_min = 40
x_SB4_min_2 = 185 # Already one min date above  # Quidico
sig_SB4_min_2 = 25
x_SB4_min_3 = 245 # Quidico  
sig_SB4_min_3 = 15
x_SB4_max = 250  
sig_SB4_max = 25
x_SB5_max = 430  
sig_SB5_max = 40
x_SB5_min = 530 
sig_SB5_min = 25
x_SB5_min_2 = 365  
sig_SB5_min_2 = 15
x_SB6_max = 420  
sig_SB6_max = 20
x_SB6_min = 245  
sig_SB6_min = 55
x_SB7_max = 460 # Quidico  
sig_SB7_max = 30
x_SB7_min = 460 
sig_SB7_min = 20
x_SB7_min_2 = 1630  
sig_SB7_min_2 = 20
x_SB7_min_3 = 570 # Quidico  
sig_SB7_min_3 = 45
x_SB7_min_4 = 495 # Quidico  
sig_SB7_min_4 = 15
x_SB8_max = 1630  
sig_SB8_max = 20
x_SB8_min = 1790  
sig_SB8_min = 20
x_basal_max = 1760  
sig_basal_max = 25
x_basal_min = 1850  
sig_basal_min = 20


o = order(cal$calbp) # - JAGS needs the calibration curve to be in order
data = list(c14=cal$c14age[o],calbp=cal$calbp[o],err=cal$errorBP[o])
rm(cal)
for(i in 4:(length(ls())+1)) {
  data[[i]] = get(ls()[i])
}
names(data) = c('c14','calbp','err',ls()[4:length(ls())])

# To get goo starting values, calibrate some radiocarbon dates and sort them
# Give some initial values
#source('Initial_values_function_20150624.R')
  
# Model definition
modelstring ='
model{
  # Priors for the age of the sand beds that are known
  theta_raw[1] ~ dnorm(-60, 10)
  theta_raw[2] ~ dnorm(-10, 10)
  theta_raw[3] ~ dnorm(115, 10)
  theta_raw[4] ~ dnorm(199, 10)
  theta_raw[5] ~ dnorm(375, 10)
  # Now prior distributions for the others
  for(i in 6:9) {
    theta_raw[i] ~ dunif(376,10000)
  }
  theta[1:9] <- sort(theta_raw) #

  # First do sand bed 2 correponding to theta[2]
  theta_SB2_min <- theta[2] - extra_SB2
  extra_SB2 ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_SB2_min <- interp.lin(theta_SB2_min,calbp,c14)
  sig_cal_SB2_min <- interp.lin(theta_SB2_min,calbp,err)
  x_SB2_min ~ dnorm(mu_cal_SB2_min,tau_all_SB2_min)
  tau_all_SB2_min <- 1/sig_sq_all_SB2_min
  sig_sq_all_SB2_min <- pow(sig_SB2_min,2) + pow(sig_cal_SB2_min,2)

  # Now sand bed 3 correponding to theta[3]
  theta_SB3_min <- theta[3] - extra_SB3
  extra_SB3 ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_SB3_min <- interp.lin(theta_SB3_min,calbp,c14)
  sig_cal_SB3_min <- interp.lin(theta_SB3_min,calbp,err)
  x_SB3_min ~ dnorm(mu_cal_SB3_min,tau_all_SB3_min)
  tau_all_SB3_min <- 1/sig_sq_all_SB3_min
  sig_sq_all_SB3_min <- pow(sig_SB3_min,2) + pow(sig_cal_SB3_min,2)

  # Now sand bed 4 correponding to theta[4]
  theta_SB4_min <- theta[4] - extra_SB4
  extra_SB4 ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_SB4_min <- interp.lin(theta_SB4_min,calbp,c14)
  sig_cal_SB4_min <- interp.lin(theta_SB4_min,calbp,err)
  x_SB4_min ~ dnorm(mu_cal_SB4_min,tau_all_SB4_min)
  tau_all_SB4_min <- 1/sig_sq_all_SB4_min
  sig_sq_all_SB4_min <- pow(sig_SB4_min,2) + pow(sig_cal_SB4_min,2)

  # Another SB4 min
  theta_SB4_min_2 <- theta[4] - extra_SB4_2
  extra_SB4_2 ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_SB4_min_2 <- interp.lin(theta_SB4_min_2,calbp,c14)
  sig_cal_SB4_min_2 <- interp.lin(theta_SB4_min_2,calbp,err)
  x_SB4_min_2 ~ dnorm(mu_cal_SB4_min_2,tau_all_SB4_min_2)
  tau_all_SB4_min_2 <- 1/sig_sq_all_SB4_min_2
  sig_sq_all_SB4_min_2 <- pow(sig_SB4_min_2,2) + pow(sig_cal_SB4_min_2,2)

  # And another SB4 min
  theta_SB4_min_3 <- theta[4] - extra_SB4_3
  extra_SB4_3 ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_SB4_min_3 <- interp.lin(theta_SB4_min_3,calbp,c14)
  sig_cal_SB4_min_3 <- interp.lin(theta_SB4_min_3,calbp,err)
  x_SB4_min_3 ~ dnorm(mu_cal_SB4_min_3,tau_all_SB4_min_3)
  tau_all_SB4_min_3 <- 1/sig_sq_all_SB4_min_3
  sig_sq_all_SB4_min_3 <- pow(sig_SB4_min_3,2) + pow(sig_cal_SB4_min_3,2)

  # Now an SB4_max
  theta_SB4_max <- theta[4] + extra_SB4_max
  extra_SB4_max ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_SB4_max <- interp.lin(theta_SB4_max,calbp,c14)
  sig_cal_SB4_max <- interp.lin(theta_SB4_max,calbp,err)
  x_SB4_max ~ dnorm(mu_cal_SB4_max,tau_all_SB4_max)
  tau_all_SB4_max <- 1/sig_sq_all_SB4_max
  sig_sq_all_SB4_max <- pow(sig_SB4_max,2) + pow(sig_cal_SB4_max,2)

  # SB5_max
  theta_SB5_max <- theta[5] + extra_SB5_max
  extra_SB5_max ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_SB5_max <- interp.lin(theta_SB5_max,calbp,c14)
  sig_cal_SB5_max <- interp.lin(theta_SB5_max,calbp,err)
  x_SB5_max ~ dnorm(mu_cal_SB5_max,tau_all_SB5_max)
  tau_all_SB5_max <- 1/sig_sq_all_SB5_max
  sig_sq_all_SB5_max <- pow(sig_SB5_max,2) + pow(sig_cal_SB5_max,2)

  # SB5 min
  theta_SB5_min <- theta[5] - extra_SB5
  extra_SB5 ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_SB5_min <- interp.lin(theta_SB5_min,calbp,c14)
  sig_cal_SB5_min <- interp.lin(theta_SB5_min,calbp,err)
  x_SB5_min ~ dnorm(mu_cal_SB5_min,tau_all_SB5_min)
  tau_all_SB5_min <- 1/sig_sq_all_SB5_min
  sig_sq_all_SB5_min <- pow(sig_SB5_min,2) + pow(sig_cal_SB5_min,2)

  # Another SB5 min
  theta_SB5_min_2 <- theta[5] - extra_SB5_2
  extra_SB5_2 ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_SB5_min_2 <- interp.lin(theta_SB5_min_2,calbp,c14)
  sig_cal_SB5_min_2 <- interp.lin(theta_SB5_min_2,calbp,err)
  x_SB5_min_2 ~ dnorm(mu_cal_SB5_min_2,tau_all_SB5_min_2)
  tau_all_SB5_min_2 <- 1/sig_sq_all_SB5_min_2
  sig_sq_all_SB5_min_2 <- pow(sig_SB5_min_2,2) + pow(sig_cal_SB5_min_2,2)

  # SB6 max
  theta_SB6_max <- theta[6] + extra_SB6_max
  extra_SB6_max ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_SB6_max <- interp.lin(theta_SB6_max,calbp,c14)
  sig_cal_SB6_max <- interp.lin(theta_SB6_max,calbp,err)
  x_SB6_max ~ dnorm(mu_cal_SB6_max,tau_all_SB6_max)
  tau_all_SB6_max <- 1/sig_sq_all_SB6_max
  sig_sq_all_SB6_max <- pow(sig_SB6_max,2) + pow(sig_cal_SB6_max,2)

  # SB6 min
  theta_SB6_min <- theta[6] - extra_SB6
  extra_SB6 ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_SB6_min <- interp.lin(theta_SB6_min,calbp,c14)
  sig_cal_SB6_min <- interp.lin(theta_SB6_min,calbp,err)
  x_SB6_min ~ dnorm(mu_cal_SB6_min,tau_all_SB6_min)
  tau_all_SB6_min <- 1/sig_sq_all_SB6_min
  sig_sq_all_SB6_min <- pow(sig_SB6_min,2) + pow(sig_cal_SB6_min,2)

  # SB7 max
  theta_SB7_max <- theta[7] + extra_SB7_max
  extra_SB7_max ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_SB7_max <- interp.lin(theta_SB7_max,calbp,c14)
  sig_cal_SB7_max <- interp.lin(theta_SB7_max,calbp,err)
  x_SB7_max ~ dnorm(mu_cal_SB7_max,tau_all_SB7_max)
  tau_all_SB7_max <- 1/sig_sq_all_SB7_max
  sig_sq_all_SB7_max <- pow(sig_SB7_max,2) + pow(sig_cal_SB7_max,2)

  # SB7 min
  theta_SB7_min <- theta[7] - extra_SB7
  extra_SB7 ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_SB7_min <- interp.lin(theta_SB7_min,calbp,c14)
  sig_cal_SB7_min <- interp.lin(theta_SB7_min,calbp,err)
  x_SB7_min ~ dnorm(mu_cal_SB7_min,tau_all_SB7_min)
  tau_all_SB7_min <- 1/sig_sq_all_SB7_min
  sig_sq_all_SB7_min <- pow(sig_SB7_min,2) + pow(sig_cal_SB7_min,2)

  # Another SB7 min
  theta_SB7_min_2 <- theta[7] - extra_SB7_2
  extra_SB7_2 ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_SB7_min_2 <- interp.lin(theta_SB7_min_2,calbp,c14)
  sig_cal_SB7_min_2 <- interp.lin(theta_SB7_min_2,calbp,err)
  x_SB7_min_2 ~ dnorm(mu_cal_SB7_min_2,tau_all_SB7_min_2)
  tau_all_SB7_min_2 <- 1/sig_sq_all_SB7_min_2
  sig_sq_all_SB7_min_2 <- pow(sig_SB7_min_2,2) + pow(sig_cal_SB7_min_2,2)
  
  # And another SB7 min
  theta_SB7_min_3 <- theta[7] - extra_SB7_3
  extra_SB7_3 ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_SB7_min_3 <- interp.lin(theta_SB7_min_3,calbp,c14)
  sig_cal_SB7_min_3 <- interp.lin(theta_SB7_min_3,calbp,err)
  x_SB7_min_3 ~ dnorm(mu_cal_SB7_min_3,tau_all_SB7_min_3)
  tau_all_SB7_min_3 <- 1/sig_sq_all_SB7_min_3
  sig_sq_all_SB7_min_3 <- pow(sig_SB7_min_3,2) + pow(sig_cal_SB7_min_3,2)

  # And yet another SB7 min
  theta_SB7_min_4 <- theta[7] - extra_SB7_4
  extra_SB7_4 ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_SB7_min_4 <- interp.lin(theta_SB7_min_4,calbp,c14)
  sig_cal_SB7_min_4 <- interp.lin(theta_SB7_min_4,calbp,err)
  x_SB7_min_4 ~ dnorm(mu_cal_SB7_min_4,tau_all_SB7_min_4)
  tau_all_SB7_min_4 <- 1/sig_sq_all_SB7_min_4
  sig_sq_all_SB7_min_4 <- pow(sig_SB7_min_4,2) + pow(sig_cal_SB7_min_4,2)

  # SB8 max
  theta_SB8_max <- theta[8] + extra_SB8_max
  extra_SB8_max ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_SB8_max <- interp.lin(theta_SB8_max,calbp,c14)
  sig_cal_SB8_max <- interp.lin(theta_SB8_max,calbp,err)
  x_SB8_max ~ dnorm(mu_cal_SB8_max,tau_all_SB8_max)
  tau_all_SB8_max <- 1/sig_sq_all_SB8_max
  sig_sq_all_SB8_max <- pow(sig_SB8_max,2) + pow(sig_cal_SB8_max,2)
  
  # SB8 min
  theta_SB8_min <- theta[8] - extra_SB8
  extra_SB8 ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_SB8_min <- interp.lin(theta_SB8_min,calbp,c14)
  sig_cal_SB8_min <- interp.lin(theta_SB8_min,calbp,err)
  x_SB8_min ~ dnorm(mu_cal_SB8_min,tau_all_SB8_min)
  tau_all_SB8_min <- 1/sig_sq_all_SB8_min
  sig_sq_all_SB8_min <- pow(sig_SB8_min,2) + pow(sig_cal_SB8_min,2)

  # basal max
  theta_basal_max <- theta[9] + extra_basal_max
  extra_basal_max ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_basal_max <- interp.lin(theta_basal_max,calbp,c14)
  sig_cal_basal_max <- interp.lin(theta_basal_max,calbp,err)
  x_basal_max ~ dnorm(mu_cal_basal_max,tau_all_basal_max)
  tau_all_basal_max <- 1/sig_sq_all_basal_max
  sig_sq_all_basal_max <- pow(sig_basal_max,2) + pow(sig_cal_basal_max,2)
  
  # basal min
  theta_basal_min <- theta[9] - extra_basal
  extra_basal ~ dgamma(4,0.02) # Gamma distribution with mean 200 and sd 100
  mu_cal_basal_min <- interp.lin(theta_basal_min,calbp,c14)
  sig_cal_basal_min <- interp.lin(theta_basal_min,calbp,err)
  x_basal_min ~ dnorm(mu_cal_basal_min,tau_all_basal_min)
  tau_all_basal_min <- 1/sig_sq_all_basal_min
  sig_sq_all_basal_min <- pow(sig_basal_min,2) + pow(sig_cal_basal_min,2)
}'

# Run the model in jags
set.seed(111) # Set the seed for repeatable results
model=jags.model(textConnection(modelstring), data=data,n.chains=4)#, inits=init
stop()
update(model,n.iter=1e6) # Burnin period
output=coda.samples(model=model,
                    variable.names=c("theta"),
                    n.iter=1e6,
                    thin=5e2)
saveRDS(output,file='output_20170302.rds')
load('output_20150624.rda')
pdf(file='trace_plot_20150624.pdf',width=8,height=8)
par(mar=c(2,2,2,2))
plot(output)
dev.off()
summary(output)

# Congergence diagnostics - all good
gelman.diag(output)

##########################################################################################

## Extract the parameters
output_2 = do.call(rbind,output)

# Now create upper and lower CIs
CIs = t(round(apply(output_2,2,'quantile',c(0.025,0.5,0.975)),0))

##########################################################################################

# First plot - tsunamis on y axis with CIs on x axis
Tsunami = factor(paste('T',12:1,sep=''),levels=paste('T',1:12,sep=''))

df_1 = data.frame(Tsunami,CIs=CIs[4+c(1:9,11:12,14),])
rownames(df_1)=NULL
colnames(df_1)=c('Tsunami','pc2.5','pc50','pc97.5')
  
#ggplot(df_1,aes(x=pc50,y=Tsunami,colour=Tsunami))+ geom_point(size=2) + geom_errorbarh(aes(xmax = pc97.5, xmin =pc2.5),height=0.5)+theme_bw()+theme(legend.position='None',axis.title.y=element_SBext(angle=0,vjust=1,hjust=0))+ggtitle('Age of Tsunamis with 95% error bounds\n')+xlab('Age (thousands of years')+ scale_x_continuous(breaks=seq(2500,7600,by=500),limits=c(2500,7600))

##########################################################################################

# An alternative version - plot densities but facet by Tsunami
output_3 = c(output_2[,'theta[1]'],output_2[,'theta[2]'],output_2[,'theta[3]'],output_2[,'theta[4]'],output_2[,'theta[5]'],output_2[,'theta[6]'],output_2[,'theta[7]'],output_2[,'theta[8]'],output_2[,'theta[9]'],output_2[,'theta[11]'],output_2[,'theta[12]'],output_2[,'theta[14]'])
Tsunami = rep(factor(paste('T',12:1,sep=''),levels=paste('T',12:1,sep='')),each=4*2000) # num_chains*num_iterations
df_2 = data.frame(Tsunami=Tsunami,Age=output_3)

ggplot(df_2,aes(x=Age,fill=Tsunami))+ geom_density(colour=NA) + facet_grid(Tsunami ~ .,scales='free')+ theme_bw()+theme(legend.position='None',axis.title.y=element_SBext(angle=0,vjust=1,hjust=0),axis.text.y = element_blank(),axis.ticks.y = element_blank())+ggtitle('Age of tsunamis\n')+xlab('Age (years BP)')+ scale_x_continuous(breaks=seq(7500,2500,by=-500),limits=c(2500,7600))+ylab("")+coord_SBrans(x="reverse", y="reverse")+theme(strip.text.y = element_SBext(size = 8, angle = 0),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())
ggsave('Tsunami_fig1.pdf',height=8,width=8)

##########################################################################################

# Create a csv file of ages for Ben
CIs_2 = t(round(apply(output_2,2,'quantile',c(0.01,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.99)),0))[4+c(1:9,11:12,14),]
rownames(CIs_2) = paste('T',12:1,sep='')
write.csv(CIs_2,file='Tsunami_ages_20150624.csv',quote=FALSE)

##########################################################################################

# Finally calculate Age differences between Tsunamis - need to be careful here
output_4 = output_2[,4+c(1:9,11:12,14)]
colnames(output_4) = paste('T',12:1,sep='')
output_diff = t(apply(output_4,1,'diff'))
colnames(output_diff) = c('T11_12','T10_11','T9_10','T8_9','T7_8','T6_7','T5_6','T4_5','T3_4','T2_3','T1_2')
output_diff_2 = as.vector(output_diff)
Tsunami = rep(factor(c('T11_12','T10_11','T9_10','T8_9','T7_8','T6_7','T5_6','T4_5','T3_4','T2_3','T1_2'),c('T11_12','T10_11','T9_10','T8_9','T7_8','T6_7','T5_6','T4_5','T3_4','T2_3','T1_2')),each=4*2000)
df_3 = data.frame(Tsunami=Tsunami,Age=output_diff_2)
ggplot(df_3,aes(x=Age,fill=Tsunami))+ geom_density(colour=NA) + facet_grid(Tsunami ~ .,scales='free')+ theme_bw() +theme(legend.position='None',axis.title.y=element_SBext(angle=0,vjust=1,hjust=0),axis.text.y = element_blank(),axis.ticks.y = element_blank())+ggtitle('Age gaps between consecutive tsunamis\n')+xlab('Age gap') + scale_x_continuous(breaks=seq(0,2600,by=200),limits=c(0,2500))+ylab("")+theme(strip.text.y = element_SBext(size = 8, angle = 0),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())
ggsave('Tsunami_fig2_20150624.pdf',height=8,width=8)

# Create a csv file of age gaps
CIs_3 = t(round(apply(output_diff,2,'quantile',c(0.01,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.99)),0))
write.csv(CIs_3,file='Tsunami_age_gaps_20150624.csv',quote=FALSE)
