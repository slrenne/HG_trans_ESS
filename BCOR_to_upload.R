#libraries
library(rethinking)
library(tidyverse)
library(survival)
library(readxl)

########################################################################################
#data import
survHESS <- read.csv("~/DB_to_upload.csv", 
                        sheet = "DB", na = "NA")
survHESS$DOD <- as.integer(survHESS$DOD)
survHESS$mo_fup <- as.integer(survHESS$mo_fup)
survHESS$trans <- as.factor((survHESS$BCOR+1))
#######################################################################################
survHESS %>%
  ggplot(aes(Age, y = ..count.., fill= trans))+
  geom_density(alpha = 0.2, position = "stack")+
  labs(x = "Age", 
       y = "Number of Cases",
       fill = "Translocation",
        title ="HG Translocated ESS")+
  scale_fill_discrete(breaks = 1:4,
                      labels  = c("YWHAE","BCOR","PHF1","SUZ12"))

FIGO.lev <- c("Stage I","Stage II",  "Stage III","Stage IV")
survHESS %>% filter(!is.na(FIGO_simpl)) %>%
  ggplot(aes(Age, y = ..count.., fill= trans))+
  geom_density(alpha = 0.2, position = "stack")+
  labs(x = "Age", 
       y = "Number of Cases",
       fill = "Translocation",
       title ="Age of HG Translocated ESS")+
  scale_fill_discrete(breaks = 1:4,
                      labels  = c("YWHAE","BCOR","PHF1","SUZ12"))+
  facet_grid(as.factor(FIGO_simpl)~., 
             labeller = labeller(FIGO_simpl = FIGO.lev))

survHESS %>%
  ggplot(aes(Mitotic_count, y = ..count.., fill= trans))+
  geom_density(alpha = 0.2, position = "stack")+
  labs(x = "Mitosis", 
       y = "Number of Cases",
       fill = "Translocation",
       title ="HG Translocated ESS")+
  scale_fill_discrete(breaks = 1:4,
                      labels  = c("YWHAE","BCOR","PHF1","SUZ12"))

########################################################################################
#plotting Kaplan-Meier of raw data
##by translocation
datasurv <-  data.frame(list( time = survHESS$mo_fup, 
                 status = survHESS$DOD,
                 group = survHESS$trans))
formula <- Surv(time, status == 1) ~ group
surv_by_group = survfit(Surv(time, status == 1) ~ group, data = datasurv)

#saving the plot
jpeg("~/Dropbox/BCOR/km_trans_raw.jpeg", units = "in", width = 5, height = 5, res = 300)
plot(surv_by_group, xlab = "Time", ylab="Overall Survival", lty = c(1:4), 
     col = c("black", "red", "blue","grey"), main="Kaplan-Meier Survival by Translocation", bty = "L")
legend( 90, 0.95, legend = c("YWAHE", "BCOR","PHF1","SUZ12"), lty = c(1:4), 
        col = c("black", "red", "blue","grey") )
dev.off()

##by stage
survHESS$FIGO_simpl <- as.factor((survHESS$FIGO_simpl))
datasurv <-  data.frame(list( time = survHESS$mo_fup, 
                              status = survHESS$DOD,
                              group = survHESS$FIGO_simpl))
formula <- Surv(time, status == 1) ~ group
surv_by_group = survfit(Surv(time, status == 1) ~ group, data = datasurv)
#saving the plot
jpeg("~/Dropbox/BCOR/km_stage_raw.jpeg", units = "in", width = 5, height = 5, res = 300)
plot(surv_by_group, xlab = "Time", ylab="Overall Survival", lty = c(1:4), 
     col = c("black", "blue","dark green","red"), 
     main="Kaplan-Meier Survival by Stage", bty = "L")
legend( 95, 1.03, legend = c("Stage I", "Stage II", "Stage III", "Stage IV"),
        col = c("black", "blue","dark green","red"), lty = c(1:4) )
dev.off()

########################################################################################
#DAG for variables
dag1 <- dagitty( "dag {
Age -> Surv
Age -> Mitosis -> Stage -> Surv
Age -> Trans -> Stage -> Surv
Trans -> Mitosis -> Surv
Trans -> Surv
}")
#saving the plot
jpeg("~/Dropbox/BCOR/DAG.jpeg", units = "in", width = 5, height = 5, res = 300)
drawdag(dag1)
dev.off()
#total effect of translocation condition on AGE
#effect of translocation not mediated by proliferation and stage at diagnosis condition on
#AGE, stage, Mitosis
impliedConditionalIndependencies (dag1)
adjustmentSets( dag1 , exposure="Trans" , outcome="Surv")

########################################################################################
#survival model conditioning on age
dat <- list(
  months_to_event = as.numeric( survHESS$mo_fup ), 
  R = as.integer(survHESS$trans) , 
  DOD = as.integer(survHESS$DOD),
  age = standardize(survHESS$Age)
)

m_tr_age <- ulam( alist(
  months_to_event|DOD==1 ~ exponential( lambda ), 
  months_to_event|DOD==0 ~ custom(exponential_lccdf( !Y | lambda )), 
  lambda <- 1.0/mu,
  log(mu) <- a_bar+z[R]*sigma_a + b[R]*age,
  age ~ dnorm (nu_a, sigma_age),
  b[R] ~ dnorm (b_bar, sigma_b),
  z[R] ~ dnorm(0,1),
  c(a_bar, b_bar,nu_a) ~ dnorm(0,1),
  c(sigma_a, sigma_b, sigma_age) ~ dexp(1)
), data=dat , chains=4 , cores=4, iter= 2000, control=list(adapt_delta=0.97)) 

plot(precis(m_tr_age,3,pars="z"),lab=c("YWAHE", "BCOR","PHF1","SUZ12"))

#Prior predictive simulation
prior <- extract.prior(m_tr_age)
pz_BCOR <-  prior$z[,2]
pz_YWHAE <- prior$z[,1]
pa_bar <- prior$a_bar
psigma_a <- prior$sigma_a
plambda_BCOR_mean <- mean(1/exp(pz_BCOR+pa_bar*psigma_a)) # @ the mean age beta= 0
plambda_YWHAE_mean <- mean(1/exp(pz_YWHAE+pa_bar*psigma_a))
plambda_BCOR <- 1/exp(pz_BCOR+pa_bar*psigma_a)
plambda_YWHAE <- 1/exp(pz_YWHAE+pa_bar*psigma_a)
N <- 1:200
jpeg("~/Dropbox/BCOR/km_prPS_tr.jpeg", units = "in", width = 5, height = 5, res = 300)
plot(dexp(N,plambda_YWHAE_mean)/plambda_YWHAE_mean, type = "n", 
     ylim = c(0,1), xlim = c(0,40), 
     xlab = "Months", ylab="Simulated Overall survival",
     main="Kaplan-Meier Survival by Translocation, Prior")
lines(dexp(N,plambda_BCOR_mean)/plambda_BCOR_mean, col = "red", lty = 2)
lines(dexp(N,plambda_YWHAE_mean)/plambda_YWHAE_mean, col = "black")
for ( i in 1:30 )
  lines(dexp(N,plambda_BCOR[i])/plambda_BCOR[i], col=col.alpha("red",0.2), lty = 2)
for ( i in 1:30 )
  lines(dexp(N,plambda_YWHAE[i])/plambda_YWHAE[i], col=col.alpha("black",0.2))
legend( 25, 0.95, legend = c("YWHAE","BCOR"), lty = 1:2, col = c("black", "red"))
dev.off()

#Posterior predictive simulation
post <- extract.samples(m_tr_age)
z_BCOR <-  post$z[,2]
z_YWHAE <- post$z[,1]
z_PHF1 <- post$z[,3]
z_SUZ12 <- post$z[,4]
a_bar <- post$a_bar
sigma_a <- post$sigma_a
lambda_BCOR_mean <- mean(1/exp(a_bar+z_BCOR*sigma_a)) # @ the mean age beta= 0
lambda_YWHAE_mean <- mean(1/exp(a_bar+z_YWHAE*sigma_a))
lambda_PHF1_mean <- mean(1/exp(a_bar+z_PHF1*sigma_a))
lambda_SUZ12_mean <- mean(1/exp(a_bar+z_SUZ12*sigma_a))
lambda_BCOR <- 1/exp(a_bar+z_BCOR*sigma_a)
lambda_YWHAE <- 1/exp(a_bar+z_YWHAE*sigma_a)

#calculating the median difference among the two grups
p_median <- data.frame(matrix(ncol = 2, nrow = 4000))
median_b <-c(log(2)/lambda_YWHAE_median, log(2)/lambda_BCOR_median) #calculate the median for exponential
names(p_median)[1] <- "ywhae"
names(p_median)[2] <- "bcor"
p_median$bcor <-  log(2)/(lambda_BCOR)
p_median$ywhae <- log2(2)/(lambda_YWHAE)
p_median$diff <- p_median$ywhae - p_median$bcor

#plotting the survival curves
N <- 1:200
jpeg("~/Dropbox/BCOR/km_postPS_tr.jpeg", units = "in", width = 5, height = 5, res = 300)
plot(dexp(N,lambda_YWHAE_mean)/lambda_YWHAE_mean, type = "n", ylim = c(0,1), xlim = c(0,150), 
     xlab = "Months", ylab="Overall survival",
     main="Kaplan-Meier Survival by Translocation")
lines(dexp(N,lambda_YWHAE_mean)/lambda_YWHAE_mean, col = "black")
lines(dexp(N,lambda_BCOR_mean)/lambda_BCOR_mean, col = "red", lty = 2) #survival distribution of exp 
lines(dexp(N,lambda_PHF1_mean)/lambda_PHF1_mean, col = "blue", lty = 3)
lines(dexp(N,lambda_SUZ12_mean)/lambda_SUZ12_mean, col = "gray", lty = 4)
legend( 90, 0.95, legend = c("YWHAE", "BCOR","PHF1","SUZ12"), lty = c(1:4), 
        col = c("black", "red", "blue","grey") )
#for ( i in 1:30 )
# lines(dexp(N,lambda_BCOR[i])/lambda_BCOR[i], col=col.alpha("red",0.2), lty = 2)
#for ( i in 1:30 )
#  lines(dexp(N,lambda_YWHAE[i])/lambda_YWHAE[i], col=col.alpha("black",0.2))
#abline(h = 0.5, col = "gray", lty = 3)
#segments(-1, 0.5, median_b[1],0.5, col = "gray", lty = 3)
#for(i in 1:2) segments(median_b[i], -2, median_b[i],0.5, col = "gray", lty = 3)
dev.off()

#calculating the lambda HPDI
lambda_bcor_hpdi <- HPDI(lambda_BCOR)
lambda_ywhae_hpdi <-HPDI(lambda_YWHAE)
#Plotting the survival curves with 0.89 HPDI intervals
x <- seq(0,200, by=0.01)
jpeg("~/Dropbox/BCOR/km_postPS_tr_v2.jpeg", units = "in", width = 5, height = 5, res = 300)
plot(NULL, bty="L", ylim = c(0,1), xlim = c(0,200), 
     xlab = "Months", ylab="Overall survival",
     main="Kaplan-Meier Survival by Translocation")
ym <- dexp(x,lambda_BCOR_mean)/lambda_BCOR_mean
y1 <- dexp(x,lambda_bcor_hpdi[2])/lambda_bcor_hpdi[2]
y2 <- dexp(x,lambda_bcor_hpdi[1])/lambda_bcor_hpdi[1]
polygon(c(x,rev(x)),c(y2,rev(y1)),col=col.alpha("red"), border = NA)
ym <- dexp(x,lambda_YWHAE_mean)/lambda_YWHAE_mean
y1 <- dexp(x,lambda_ywhae_hpdi[2])/lambda_ywhae_hpdi[2]
y2 <- dexp(x,lambda_ywhae_hpdi[1])/lambda_ywhae_hpdi[1]
polygon(c(x,rev(x)),c(y2,rev(y1)),col=col.alpha("black"), border = NA)
lines(dexp(N,lambda_BCOR_mean)/lambda_BCOR_mean, col = "red", lty = 2) #survival distribution of exp 
lines(dexp(N,lambda_YWHAE_mean)/lambda_YWHAE_mean, col = "black")
legend( 0, 0.3, legend = c("YWHAE","BCOR"), lty = 1:2, col = c("black", "red"))
dev.off()

#plotting the median difference
jpeg("~/Dropbox/BCOR/Mdiff_tr.jpeg", units = "in", width = 5, height = 5, res = 300)
dens(p_median$diff, main = "Difference of Median Survival", 
     ylab = "Posterior Probability Density",
     xlab = "Months lost with BCOR rearrangement", 
     xlim = c(-50, 400), show.HPDI=.89, show.zero = TRUE)
#arrows(300,0.003,380,0.003)
#text(300,0.003, labels = "BCOR worse", pos = c(3))
legend(270, 0.004, pch = 15, legend = "89% HPDI", 
       col = col.alpha("black", 0.15), pt.cex = 2)
dev.off()
plot(precis(as.data.frame(p_median)))

########################################################################################
#multivariate model including mitosis and stage
i= !is.na(survHESS$FIGO_simpl)
dat <- list(
  months_to_event = as.numeric( survHESS$mo_fup )[i], 
  R = as.integer(survHESS$trans)[i] , 
  S = as.integer(survHESS$FIGO_simpl)[i],
  DOD = as.integer(survHESS$DOD)[i],
  M = standardize(survHESS$Mitotic_count)[i],
  age= standardize(survHESS$Age)[i]
)

m_tr_st_mit_age <- ulam( alist(
  months_to_event|DOD==1 ~ exponential( lambda ), 
  months_to_event|DOD==0 ~ custom(exponential_lccdf( !Y | lambda )), 
  lambda <- 1.0/mu,
  log(mu) <- alpha[R,S] + beta[R,S]*M + gamm[R,S]*age, 
  # adaptive priors - non-centered 
  transpars> matrix[R,4]:alpha <-
    compose_noncentered( sigma_a , L_Rho_a , z_a ), 
  transpars> matrix[R,4]:beta <-
    compose_noncentered( sigma_b , L_Rho_b , z_b ), 
  transpars> matrix[R,4]:gamm <-
    compose_noncentered( sigma_g , L_Rho_g , z_g ), 
  matrix[4,R]:z_a ~ normal( 0 , 1 ), 
  matrix[4,R]:z_b ~ normal( 0 , 1 ),
  matrix[4,R]:z_g ~ normal( 0 , 1 ),
  # fixed priors
  vector[4]:sigma_a ~ dexp(1), 
  cholesky_factor_corr[4]:L_Rho_a ~ lkj_corr_cholesky( 2 ),
  vector[4]:sigma_b ~ dexp(1), 
  cholesky_factor_corr[4]:L_Rho_b ~ lkj_corr_cholesky( 2 ),
  vector[4]:sigma_g ~ dexp(1), 
  cholesky_factor_corr[4]:L_Rho_g ~ lkj_corr_cholesky( 2 ),
  # compute ordinary correlation matrixes from Cholesky factors
  gq> matrix[4,4]:Rho_a <<- multiply_lower_tri_self_transpose(L_Rho_a),
  gq> matrix[4,4]:Rho_b <<- multiply_lower_tri_self_transpose(L_Rho_b), 
  gq> matrix[4,4]:Rho_g <<- multiply_lower_tri_self_transpose(L_Rho_g), 
  #impute the NA in age & mitosis
  M ~ dnorm (nu_m,sigma_m),
  nu_m ~ dnorm( 0,1),
  sigma_m ~ dexp(1),
  age ~ dnorm (nu_age,sigma_age),
  nu_age ~ dnorm( 0,1),
  sigma_age ~ dexp(1)
), data=dat , chains=4 , cores=4, iter = 2000 ) 


#Posterior predictive simulation
post <- extract.samples(m_tr_st_mit_age)

#plotting surv by stage
jpeg("~/Dropbox/BCOR/km_postPS_by_stage.jpeg", 
     units = "in", width = 7.5, height = 7.5, res = 300)
R_col <- c("black","red")
N <- 1:200
x <- seq(0,200, by=0.01)
par(mfrow=c(2,2))
for(i in 1:4) {
  plot(NULL, bty= "L", ylim = c(0,1), xlim = c(0,150), 
       xlab = "Months", ylab="Overall survival",
       main=paste("Stage",i))
  for(n in 1:2) {
  lambda <- 1/exp(post$alpha[,n,i])
  hdpi_lambda <- HPDI(lambda)
  lambda_m <- mean(lambda)
  y1 <- dexp(x,hdpi_lambda[2])/hdpi_lambda[2]
  y2 <- dexp(x,hdpi_lambda[1])/hdpi_lambda[1]
  polygon(c(x,rev(x)),c(y2,rev(y1)),col=col.alpha(R_col[n]), border = NA)
  lines(dexp(N,lambda_m)/lambda_m, col = R_col[n], lty = n)}}
dev.off()

#plotting median difference by stage
jpeg("~/Dropbox/BCOR/median_diff_by_stage.jpeg", 
     units = "in", width = 7.5, height = 7.5, res = 300)
par(mfrow=c(2,2))
for(i in 1:4) {
  R <- matrix(ncol = 2, nrow = 4000)
  for(n in 1:2) {
    lambda <- 1/exp(post$alpha[,n,i])
    R[,n] <-  log(2)/(lambda)}
  diff <- R[,1] - R[,2]
  dens(diff, main=paste("Stage",i), 
       ylab = "Posterior Probability Density",
       xlab = "Months lost with BCOR rearrangement", 
       xlim = c(-500, 500), 
       show.HPDI=.89, show.zero = TRUE)}
dev.off()

#precis model 2a

jpeg("~/Dropbox/BCOR/precis_tr_st_mit_age.jpg", 
     units = "in", width = 7.5, height = 7.5, res = 300)
lab_m2_sub <- c("YWHAE","BCOR","PHF1","SUZ12")
labs_m2 <- c("Stage I","Stage II","Stage III","Stage IV",
             "Stage I","Stage II","Stage III","Stage IV",
             "Stage I","Stage II","Stage III","Stage IV",
             "Stage I","Stage II","Stage III","Stage IV")
plot(precis(m_tr_st_mit_age,3,pars = "alpha"), 
     lab=labs_m2, main = "Coefficients once Age and Mitosis are known")
abline(h= c(4.5,8.5,12.5))
text(x=c(0,0,0,0), y=c(16,12,8,4)-0.4, labels = lab_m2_sub)
arrows(6, 0.3, 4, 0.3, length = 0.1, angle = 30,
       code = 2, col = "red", lty = par("lty"),
       lwd = par("lwd"))
text(x=c(3,7), y=0.3, labels = c("worse","better"), col = "red")
text(x=5,y=0.7,labels = "survival", col = "Red")
dev.off()

# Precis model 2b
jpeg("~/Dropbox/BCOR/precis_tr_st_mit_age_Mitosis.jpg", 
     units = "in", width = 7.5, height = 7.5, res = 300)
lab_m2_sub <- c("YWHAE","BCOR","PHF1","SUZ12")
labs_m2 <- c("Stage I","Stage II","Stage III","Stage IV",
             "Stage I","Stage II","Stage III","Stage IV",
             "Stage I","Stage II","Stage III","Stage IV",
             "Stage I","Stage II","Stage III","Stage IV")
plot(precis(m_tr_st_mit_age,3,pars = "beta"), 
     lab=labs_m2, main = "Coefficients of Mitosis once Age and Stage are known")
abline(h= c(4.5,8.5,12.5))
text(x=c(0,0,0,0)+2, y=c(16,12,8,4)-0.4, labels = lab_m2_sub)
arrows(6, 0.3, 4, 0.3, length = 0.1, angle = 30,
       code = 2, col = "red", lty = par("lty"),
       lwd = par("lwd"))
text(x=c(3,7), y=0.3, labels = c("worse","better"), col = "red")
text(x=5,y=0.7,labels = "survival", col = "Red")
dev.off()