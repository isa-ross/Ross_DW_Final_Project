# Code for running the full Full Bayesian GLMM relating the illegal killing of elephants 
#(represented by PIKE) to covariates (armed conflict, household wealth, etc). 
# Timothy Kuiper timothykuiper@gmail.com
# April 2022 

##NOTE - to run the supplementary models with HDI instead of IWI (household wealth), 
#and with the separate HDI indices (education, income, health), remove the "#"
#before these variable (e.g. "#beta.hdi") below and add the # to the variabels you want
#to exclude (see NOTES below)

#Set the folder where the data files are:
setwd("path to folder")
# Load the data frame with PIKE and covariate data
all.data.final<-read.csv("all.data.final.csv")


#Convert key variables factors
all.data.final$yr<-as.factor(all.data.final$yr)
all.data.final$ccode<-as.factor(all.data.final$ccode)
all.data.final$sitecode<-as.factor(all.data.final$sitecode)

#### load R packages for JAGS####
set.seed(12345)
library(R2jags)
library(BayesianTools)

# Read the PIKE and covariate data into JAGS format
jags.data.full <- list("illegal" = all.data.final$illegal, 
                       "total" = all.data.final$totcarc, 
                       "Nobs" = NROW(all.data.final),
                       "year" = all.data.final$year, 
                       "Y" = nlevels(all.data.final$yr), 
                       "country" = as.numeric(all.data.final$ccode),
                       "nc" = nlevels(all.data.final$ccode), 
                       "site"  = as.numeric(all.data.final$sitecode),
                       "ns" = nlevels(all.data.final$sitecode),
                       "wgi" = c(scale(all.data.final$wgi)),  #scale variables Z transform
                       "ivory.price" = c(scale(all.data.final$ivory.price.toan)),
                       "precip.anom" = c(scale(all.data.final$anom)),
                       "ndvi" = c(scale(all.data.final$ndvi)),
                       "hdi" = c(scale(all.data.final$hdi)),
                       "travel" = c(scale(all.data.final$travel.log)),
                       "iwi" = c(scale(all.data.final$iwi)),
                       "lec" = c(scale(all.data.final$lec)),
                       #NOTE - change the below to conflict.2yrs for the supplementary model
                       "conflict" = c(scale(all.data.final$conflict.1yr)),
                       "area" = c(scale(all.data.final$site.area.log)),
                       "ele.dens" = c(scale(all.data.final$ele.density.log)),
                       "ele.pop" = c(scale(all.data.final$ele.pop.log)),
                       "health" = c(scale(all.data.final$health)),
                       "edu" = c(scale(all.data.final$edu)),
                       "income" = c(scale(all.data.final$income)),
                       "e" = 0.0001 #error term to avoid dividing by zero in the goodness of fit test
                       
)

######################################################################################################
# PIKE model with Bayesian Lasso
######################################################################################################

pike.model <- function(){                                                                             
  
  for (i in 1:Nobs) {
    illegal[i] ~ dbin(p.illegal[i], total[i])
    p.logit[i] <-  beta0 + 
      
      #random effects
      yr.effect[year[i]-2001] +
      site.effect[site[i]] +
      c.effect[country[i]] +
      site.year.effect[site[i],year[i]-2001]+
      
      #fixed effects
      
      beta.wgi *wgi[i] +
      beta.ivory.price* ivory.price[i] +
      beta.precip.anom *precip.anom[i] +
      #beta.hdi* hdi[i] + #NOTE - remove the # to activate this variable
      beta.ndvi *ndvi[i] +
      beta.travel* travel[i]+
      beta.iwi* iwi[i]+ #NOTE - add a # to deactivate this variable
      beta.lec* lec[i]+
      beta.conflict*conflict[i]+
      beta.area*area[i]+
      beta.ele.d *ele.dens[i] +
      beta.ele.p* ele.pop[i]+
      beta.health* health[i]#+#NOTE - remove the # to activate the below variables
      #beta.edu* edu[i]+ #NOTE - remove the # to activate this variable
      #beta.income* income[i] #NOTE - remove the # to activate this variable
    
    p.illegal[i] <- 1/(1+exp(-p.logit[i]))
    
    #impute missing values
    
    ivory.price[i] ~ dnorm(0, 1) 
    precip.anom[i] ~ dnorm(0,1) #no rainfall anomaly data for 2020
    wgi[i] ~ dnorm(0,1) #no WGI data for 2020
    lec[i] ~ dnorm(0, 1) 
    
    #Model fit assess - Chi-square and posterior predictive checks
    #expected value of binomial = n x p = tot.carc * pike
    chi2[i] <- pow((illegal[i]-total[i]*p.illegal[i]),2) / (sqrt(total[i]*p.illegal[i])+e) # obs.
    illegal.new[i] ~ dbin(p.illegal[i], total[i]) # Replicate (new) data set
    chi2.new[i] <- pow((illegal.new[i]-total[i]*p.illegal[i]),2) / (sqrt(total[i]*p.illegal[i])+e) # obs
    
  } # end of i loop
  
  #These loops set up the random intercept for each year, country and site
  #These values are then fed into the p.logit[i] equation above
  
  for (yr in 1:Y) {  
    yr.effect[yr] ~ dnorm(0, year.prec)
  }
  
  for (c in 1:nc) {   
    c.effect[c] ~ dnorm(0, country.prec)
  }
  
  for (s in 1:ns) {   
    site.effect[s] ~ dnorm(0, site.prec)
  }  
  
  # Add hierarchical model for site.year effects
  for(s in 1:ns){
    for(y in 1:Y){
      site.year.effect[s,y] ~ dnorm(0, site.year.prec)
    }
  }
  
  #Priors
  year.sigma ~ dgamma(1,1)#dunif(0.001, 100)
  year.prec <- pow(year.sigma, -2)
  
  site.sigma ~ dgamma(1,1)#dunif(0.001, 100)
  site.prec <- pow(site.sigma, -2)
  
  country.sigma ~ dgamma(1,1)#dunif(0.001, 100)
  country.prec <- pow(country.sigma, -2)
  
  site.year.sigma ~ dgamma(1,1)
  site.year.prec<-pow(site.year.sigma,-2)
  
  beta0 ~ dnorm(0, 0.001)
  
  # L1 regularization == a Laplace (double exponential) prior 
  # Lasso
  
  beta.wgi  ~ ddexp(0, lambda)
  beta.ivory.price   ~ ddexp(0, lambda)
  beta.precip.anom  ~ ddexp(0, lambda)
  #beta.hdi       ~ ddexp(0, lambda) #NOTE - remove the # to activate this variable
  beta.ndvi      ~ ddexp(0, lambda)
  beta.travel     ~ ddexp(0, lambda)
  beta.iwi        ~ ddexp(0, lambda)
  beta.lec        ~ ddexp(0, lambda)
  beta.conflict  ~ ddexp(0, lambda)
  beta.area      ~ ddexp(0, lambda)
  beta.ele.d      ~ ddexp(0, lambda)
  beta.ele.p     ~ ddexp(0, lambda)
  beta.health        ~ ddexp(0, lambda)
  #beta.edu       ~ ddexp(0, lambda) #NOTE - remove the # to activate this variable
  #beta.income       ~ ddexp(0, lambda) #NOTE - remove the # to activate this variable
  
  lambda2 ~ dgamma(1,1) # shape, rate
  lambda <- sqrt(lambda2)
  
  # Add up discrepancy measures for entire data set (Goodness of fit)
  fit <- sum(chi2[]) # Omnibus test statistic actual data
  fit.new <- sum(chi2.new[]) # Omnibus test statistic replicate data
  # range of data as a second discrepancy measure
  obs.range <- max(illegal[]) - min(illegal[])
  exp.range <- max(illegal.new[]) - min(illegal.new[])
  
}

#Initialise
jags.inits.lasso <- function(){
  list(beta0 = rnorm(1, 0, 0.01),
       year.sigma = dgamma(1,1,1),#runif(1, 0.1, 1),
       site.sigma = dgamma(1,1,1),#runif(1, 0.1, 1),
       country.sigma = dgamma(1,1,1),#runif(1, 0.1, 1)
       site.year.sigma = dgamma(1,1,1)#runif(1, 0.1, 1)
  )}


# Run the final JAGS Bayesian GLMM:
pike.lasso<-jags.parallel(model.file = pike.model, working.directory = getwd(),
                          data = jags.data.full, inits = jags.inits.lasso,
                          DIC = F, n.thin = 50, n.cluster = 5,
                          n.chains = 5, n.iter = 100000, n.burnin = 50000,
                          parameters.to.save = c("beta0", "beta.wgi","beta.ivory.price", 
                                                 "beta.precip.anom",#"beta.hdi", #NOTE - remove the # to activate this variable
                                                 "beta.ndvi", "beta.travel","beta.iwi", 
                                                 "beta.lec",'beta.conflict',
                                                 "beta.area","beta.ele.p",'beta.ele.d',
                                                 "beta.health",
                                                 #"beta.edu","beta.income",#NOTE - remove the # to activate this variable
                                                 "p.illegal", "lambda",
                                                 "site.effect", "yr.effect", "c.effect",
                                                 "site.sigma", "year.sigma", "country.sigma",
                                                 "site.year.effect","site.year.sigma",
                                                 "illegal.new", "chi2.new", "chi2", "fit", 
                                                 "fit.new", "obs.range", "exp.range"))

# Optional - save the model output
save(pike.lasso, file = "name_file.RData")

##Now check fit and fit.new (Goodness of fit)
output<-pike.lasso$BUGSoutput$sims.list
par(mfrow = c(1, 1), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
plot(output$fit, output$fit.new, xlim = c(200, 1000), ylim =
       c(200, 1000), main = "", xlab = "Discrepancy observed data", ylab = "Discrepancy
expected data", frame.plot = F, cex = 1.5)
abline(0,1, lwd = 2)

#Compute Bayesian p-value
bpv <- mean(output$fit.new > output$fit)
bpv

#Check convergence

#install.packages('MCMCvis')
library(MCMCvis)
MCMCtrace(pike.lasso, 
          params = c("beta0", "beta.wgi","beta.ivory.price", 
                     "beta.precip.anom",#"beta.hdi",
                     "beta.ndvi", "beta.travel", "beta.iwi", 
                     "beta.lec",'beta.conflict',
                     "beta.area","beta.ele.p",'beta.ele.d',
                     "beta.health",#"beta.edu", "beta.income",#NOTE - remove the # to activate this variable
                     "site.sigma", "year.sigma", "country.sigma",
                     "site.year.sigma"),
          ISB = FALSE, 
          exact = TRUE,
          pdf = FALSE)

###Table of results (check Rhat convergence values for all parameters)

library(dplyr)
output<-pike.lasso$BUGSoutput$summary # a matrix
output; rownames(output);colnames(output)
covariate<-rownames(output);s<-output[,c(1,3,7,8,9)]
res<-data.frame(covariate,s)
!grepl("]",res$covariate) #rows with covariates that do not have "[" in name
res<-filter(res,!grepl("]",res$covariate));row.names(res)<-NULL
res<-rename(res, lower=X2.5.,upper=X97.5.)
#includes zero?
res$sig<-res$lower<0&res$upper>0
res

# Create the plot showing 90% CI for covariate effects

output<-pike.lasso$BUGSoutput
a<-output$sims.list
mean<-unlist(lapply(a,mean))
lower<-unlist(lapply(a,quantile,probs=0.05))
upper<-unlist(lapply(a,quantile,probs=0.95))
res<-data.frame(covariate=names(a),mean,lower,upper)
row.names(res)<-NULL;res
#includes zero?
res$sig<-res$lower<0&res$upper>0

#rename covariates

res$covariate[res$covariate=="beta.hdi"]<-"Subnational Human Development Index"
res$covariate[res$covariate=="beta.income"]<-"Household Income (subnational HDI)"
res$covariate[res$covariate=="beta.edu"]<-"Household Education (subnational HDI)"
res$covariate[res$covariate=="beta0"]<-"Intercept"
res$covariate[res$covariate=="beta.area"]<-"Area of Site "
res$covariate[res$covariate=="beta.travel"]<-"Travel Time to Site"
res$covariate[res$covariate=="beta.precip.anom"]<-"Precipitation Anomaly"
res$covariate[res$covariate=="beta.ndvi"]<-"Vegetation Density (NDVI)"
res$covariate[res$covariate=="beta.lec"]<-"Law Enforcement Capacity"
res$covariate[res$covariate=="beta.iwi"]<-"Household Wealth"
res$covariate[res$covariate=="beta.ivory.price"]<-"Global Ivory Price"
res$covariate[res$covariate=="beta.ele.p"]<-"Elephant Pop. Size"
res$covariate[res$covariate=="beta.ele.d"]<-"Elephant Pop. Density"
res$covariate[res$covariate=="beta.conflict"]<-"Armed Conflict Intensity"
res$covariate[res$covariate=="beta.health"]<-"Household Health (Life Expectancy Index)"
res$covariate[res$covariate=="beta.wgi"]<-"National Governance Quality (WGI)"
res$covariate[res$covariate=="country.sigma"]<-"Country random effect sigma"
res$covariate[res$covariate=="site.sigma"]<-"Site random effect sigma"
res$covariate[res$covariate=="year.sigma"]<-"Year random effect sigma"
res$covariate[res$covariate=="site.year.sigma"]<-"Site-year random effect sigma"
res$covariate[res$covariate=="lambda"]<-"Lambda param (LASSO-regularisation)"

#produce plot
library(ggplot2)
res2<-res[1:13,] #here you can choose which results to plot (just want the 13 main covariates here)
ggplot(res2,aes(y=covariate,x=mean,col=sig))+
  geom_point()+
  geom_errorbarh(aes(xmax = upper, xmin=lower,height=0.2),show.legend = FALSE)+
  scale_color_manual(values=c("coral3", "grey30"))+
  geom_vline(xintercept = 0, linetype="dotted",size=0.8)+
  ylab("Covariate")+
  xlab("Mean and 90% Credible Interval for Covariate Coefficient")+
  theme(legend.position = "none",
        axis.title.x=element_text(vjust=0.5),
        axis.title.y=element_text(hjust=0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 7),)




