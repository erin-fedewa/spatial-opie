# notes ----
#Erin Fedewa

#Goal: Mixed (GLMM) vrs fixed (GLS) effects models to test effects of cold pool/bottom temp on snow crab  
  #spatial extent (D95) and Lat center of dist (COD) across size/sex categories 
##See DFA script for DFA_Trend covariate development- trend is reduced from cold pool extent/center and avg temp 

#Initial Approach: 
  #Proceed with Zuur lme approach with sex/sex as random effect, DFA trend as fixed effect:
      #select random effects structure/best model with REML/AICc and assess goodness of fit/model output 
      #Model selection for correlation structure or using abundance as fixed effect for D95 model were not 
          #done for initial approach!
#Final Approach: 
  #1) Center (subtract by avg) response variables to account for different intercepts for size/sex category
  #2) Fit full fixed effects model and compare correlation structures with AICc using REML
  #3) Re-fit lme random slope and gls full models with selected correlation structure and compare w/ AICc/REML
  #4) D95 full fixed effects structure includes abundance and interaction as fixed effects
      # so re-fit with ML and compare nested models to evaluate support for fixed terms 
  #5) Refit selected model with REML


#packages ----

library(tidyverse)
library(nlme)
library(MuMIn)
library(ggeffects)
library(cowplot)
library(forecast)

#data ----

setwd('//Nmfs/akc-kod/Research/Spatial Opie MS/Datasets')
setwd("C:/Users/erin.fedewa/Work/NBS/Distribution MS/Datasets")

dat <- read.csv("Spatial_Env_bySizeSex.csv")
head(dat)
str(dat)

dat %>%
  filter(SIZESEX != "POP") %>%
  group_by(SIZESEX) %>%
  mutate(center.D95 = jtools::center(D95)) %>%
  mutate(center.COD = jtools::center(LAT_COD)) ->dat #could also scale, but that would standardize variance 


#############################################

# D95 Models ----

#Initial approach: Mixed Effects Models----

#Fit full Fixed Effects Model 
m0<-gls(D95~1+DFA_Trend, method="REML", data = dat) 
  resid<-resid(m0, type="normalized") #extract normalized residuals 
  acf(resid) #autocorrelation of residuals! Add AR1 correlation structure 

m1<-gls(D95~1+DFA_Trend, method="REML", correlation=corAR1(form = ~1|SIZESEX), data = dat)
#Correlation structure form: temporal order of data is specified by Size/Sex

#Select best random effects structure via REML:
  #Random intercept for SIZESEX
m2<-lme(D95~1+DFA_Trend, random= ~1|SIZESEX, method="REML", correlation=corAR1(form = ~1|SIZESEX), data = dat)

#Random intercept and slope for size/sex
m3<-lme(D95~1+DFA_Trend, random= ~1+DFA_Trend|SIZESEX, method="REML", correlation=corAR1(form = ~1|SIZESEX), data = dat)
  #no degrees of freedom left- model over-parameterized 

#Random slope for size/sex
m4<-lme(D95~1+DFA_Trend, random= ~-1+DFA_Trend|SIZESEX, method="REML", correlation=corAR1(form = ~1|SIZESEX), data = dat)        

AIC<- AICc(m1, m2, m4)
AIC <- AIC %>%
  mutate(model=rownames(AIC)) %>%
  mutate(dAICc=AICc-min(AICc)) %>%
  arrange(AICc)
AIC   #AICc favors random intercepts for each size/sex, model m2
  summary(m2)
  random.effects(m2)
  qqnorm(m2, ~ranef(., level=1))

#Run final model with REML estimation 
mfinal<- lme(D95~1+DFA_Trend, random= ~1|SIZESEX, method="REML", correlation=corAR1(form = ~1|SIZESEX), data = dat)
  summary(mfinal)
  intervals(mfinal) #Error
#Calculate conditional (fixed+random effects) and marginal rsq (fixed effects) for mixed models (Nakagawa and Schielzeth 2013)
  r.squaredGLMM(mfinal) #Proportion of variance explained by fixed effect alone is very low!!

#Diagnostics
E2 <- resid(mfinal, type = "normalized")
F2 <- fitted(mfinal)
plot(x = F2, y = E2, xlab = "Fitted values", ylab = "Residuals")

#Plots 
# Extract the prediction data frame
pred.mm <- ggpredict(mfinal, terms = c("DFA_Trend"))  # this gives overall predictions for the model
  pred.mm

# Plot the model predicted fit  
ggplot(pred.mm) + 
  geom_line(aes(x = x, y = predicted)) +    # slope
  geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = dat,                      # adding the raw data (scaled values)
             aes(x = DFA_Trend, y = D95, colour = SIZESEX)) + 
  labs(y="Snow Crab Areal Extent", x= "Env Trend") +
  scale_color_discrete(labels = c("Immature Females", "Mature Females", "31 to 60mm Males", "61 to 90mm Males","91 to 120mm Males")) +
  theme_bw() +theme(legend.title=element_blank())  

########################
#Final approach: centered D95 response, model selection for fixed vrs mixed effects models ----

  #Rationale: given that random effects are explaining most of the variation in lme, run stepwise model selection
    #with D95 mean-centered by size/sex (i.e. controlling for random intercept) b/c we aren't interested in 
    #variation in intercept across size/sex categories 

#1) Compare correlation structures for full Fixed Effects Model (Fixed: Trend and Abundance)
corr0<-gls(center.D95~1+DFA_Trend*TOTAL_ABUN_mil, data = dat) #OLS base model 
  acf(resid(corr0)) 
corr1<-gls(center.D95~1+DFA_Trend*TOTAL_ABUN_mil, correlation=corAR1(), data = dat) 

#Specify Autocorrelation nested by Size/Sex (vrs default order of data by year) and compare
corr2<-gls(center.D95~1+DFA_Trend*TOTAL_ABUN_mil, correlation=corAR1(form = ~1|SIZESEX), data = dat) 
  AICc(corr0, corr1, corr2) #model 1 favored 
  
#Check whether first order autocorrelation is sufficient 
  auto.arima(dat$center.D95, trace=1) #p=# of terms, q=# of lagged forecast errors 
corr3<-gls(center.D95 ~1+DFA_Trend*TOTAL_ABUN_mil, correlation = corARMA(p=2,q=1), data=dat)
  AICc(corr1, corr3) #model 3 favored, although only marginal improvement from AR1

#2)Re-fit lme random slope and gls full models with selected correlation structure and compare w/ AICc/REML
  
#Full Fixed Effects Model --Mixed effects model is overparamaterized with interaction- leave out? 
m1<-gls(center.D95~1+DFA_Trend+TOTAL_ABUN_mil, method="REML", correlation=corARMA(p=2,q=1), data = dat) 

#Full Mixed Effects model with random slope for sizesex
m2<-lme(center.D95~1+DFA_Trend+TOTAL_ABUN_mil, random= ~-1+DFA_Trend|SIZESEX, 
        method="REML", correlation=corARMA(p=2,q=1), data = dat)
  AICc(m1, m2)
  
#No need to look at random intercept  b/c response already scaled for different intercepts by size/sex
  #Model selection favors mixed effects model with random slope  

#3) Test for best-supported fixed structure using likelihood ratio test comparing nested models fit with ML
 
#Full model 
m1<-lme(center.D95~1+DFA_Trend+TOTAL_ABUN_mil, random= ~-1+DFA_Trend|SIZESEX, 
        method="ML", correlation=corARMA(p=2,q=1), data = dat) 
summary(m1)

#Drop abundance fixed effect
m2<-lme(center.D95~1+DFA_Trend, random= ~-1+DFA_Trend|SIZESEX, 
        method="ML", correlation=corARMA(p=2,q=1), data = dat)
  AICc(m1,m2) #model without abundance is supported, though very marginal difference in AIC scores 
  
#4) Re-fit final model with REML estimation 
  
mod_final <-lme(center.D95~1+DFA_Trend, random= ~-1+DFA_Trend|SIZESEX, 
                     method="REML", correlation=corARMA(p=2,q=1), data = dat) 
summary(mod_final)
  #Env trend fixed effect significant

#Diagnostics----
plot(mod_final)
  resid<-resid(mod_final, type="normalized") #extract normalized residuals 
  fit<-fitted(mod_final)
plot(fit, resid)
abline(h=0)
  hist(resid) #Looks good
  acf(resid) #Also looks good 

#Plots ----

#Color-blind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Raw data plot with smoothed trend line 
dat %>%
  ggplot(aes(DFA_Trend, D95, col=SIZESEX, shape=SIZESEX)) +
  geom_point() +
  geom_smooth(method="lm", se=F) + #raw data smoothed
  labs(y="D95 (nm²)", x= "Environmental Trend") +
  scale_color_manual(labels = c("Immature Females", "Mature Females", "31 to 60mm Males", "61 to 90mm Males","91 to 120mm Males"),
            values=cbPalette) +
   theme_bw() +theme(legend.title=element_blank())+
  theme(panel.grid = element_blank(),legend.position = "none")->g1


#Plots 
# Extract the prediction data frame
pred.mm <- ggpredict(mod_final, terms = c("DFA_Trend"))  # this gives overall predictions for the model
pred.mm

# Plot the model predicted fit  
ggplot(pred.mm) + 
  geom_line(aes(x = x, y = predicted)) +    # slope
  geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = dat,                      # adding the raw data (centered values)
             aes(x = DFA_Trend, y = center.D95, colour = SIZESEX)) + 
  labs(y="Mean-centered Snow Crab D95 (nm²)", x= "Environmental Trend") +
  scale_color_manual(labels = c("Immature Females", "Mature Females", "31 to 60mm Males", "61 to 90mm Males","91 to 120mm Males"),
      values=cbPalette) +
  theme_bw() +theme(legend.title=element_blank())-> g2

#########################################################

# Center of Distribution Models ----

#Initial approach: Mixed Effects Models----

#Fit full Fixed Effects Model 
m0<-gls(LAT_COD~1+DFA_Trend, method="REML", data = dat) 
  resid<-resid(m0, type="normalized") #extract normalized residuals 
  acf(resid) #autocorrelation of residuals! Add AR1 correlation structure 

m1<-gls(LAT_COD~1+DFA_Trend, method="REML", correlation=corAR1(form = ~1|SIZESEX), data = dat)
AICc(m0,m1)

#Random intercept for SIZESEX
m2<-lme(LAT_COD~1+DFA_Trend, random= ~1|SIZESEX, method="REML", correlation=corAR1(form = ~1|SIZESEX), data = dat)

#Random intercept and slope on bottom temp for size/sex
m3<-lme(LAT_COD~1+DFA_Trend, random= ~1+DFA_Trend|SIZESEX, method="REML", correlation=corAR1(form = ~1|SIZESEX), data = dat)

#Random slope on bottom temp for size/sex
m4<-lme(LAT_COD~1+DFA_Trend, random= ~-1+DFA_Trend|SIZESEX, method="REML", correlation=corAR1(form = ~1|SIZESEX), data = dat)        

AIC<- AICc(m1, m2, m3, m4)
AIC <- AIC %>%
  mutate(model=rownames(AIC)) %>%
  mutate(dAICc=AICc-min(AICc)) %>%
  arrange(AICc)
AIC   #AICc favors random intercepts for each size/sex, model m2
  summary(m2) 
  random.effects(m2)
  qqnorm(m2, ~ranef(., level=1))

#Run final model with REML estimation 
mfinal<- lme(LAT_COD~1+DFA_Trend, random= ~1|SIZESEX, method="REML", correlation=corAR1(form = ~1|SIZESEX), data = dat)
  summary(mfinal) #Env trend is not significant 
  intervals(mfinal) #Error
  
#Calculate conditional (fixed+random effects) and marginal rsq (fixed effects) for mixed models (Nakagawa and Schielzeth 2013)
r.squaredGLMM(mfinal) #Proportion of variance explained by fixed effect alone is very low!!

######################################

#Final approach: centered COD response, model selection for fixed vrs mixed effects models ----

#1) Compare correlation structures for full Fixed Effects Model (Fixed: Trend)
corr0<-gls(center.COD~1+DFA_Trend, data = dat) #OLS base model 
  acf(resid(corr0)) #Try AR1 
corr1<-gls(center.COD~1+DFA_Trend, correlation=corAR1(), data = dat) 

#Specify Autocorrelation nested by Size/Sex (vrs default order of data by year) and compare
corr2<-gls(center.COD~1+DFA_Trend, correlation=corAR1(form = ~1|SIZESEX), data = dat) 
AICc(corr0, corr1, corr2) #model 1 favored 

#Check whether first order autocorrelation is sufficient 
auto.arima(dat$center.COD, trace=1) #looks like ARMA(2,0,2) is best correlation structure 
corr3<-gls(center.COD ~1+DFA_Trend, correlation = corARMA(p=2,q=2), data=dat)
AICc(corr1, corr3) #model 3 favored

#2)Re-fit lme random slope and gls full models with selected correlation structure and compare w/ AICc/REML

#Full Fixed Effects Model
mcentered1<-gls(center.COD~1+DFA_Trend, method="REML", correlation = corARMA(p=2,q=2), data = dat)

#Random slope for SIZESEX
mcentered2<-lme(center.COD~1+DFA_Trend, random= ~-1+DFA_Trend|SIZESEX, method="REML", correlation = corARMA(p=2,q=2), data = dat)
  AICc(mcentered1, mcentered2) #Model selection favors fixed effects model 

#3) Final model summary 
mod_final <- gls(center.COD~1+DFA_Trend, method="REML", correlation = corARMA(p=2,q=2), data = dat)
  summary(mod_final)
  intervals(mod_final)
#DFA environmental trend fixed effect NOT significant

#Diagnostics----
plot(mod_final)
  resid<-resid(mod_final, type="normalized") #extract normalized residuals 
  fit<-fitted(mod_final)
plot(fit, resid)
  hist(resid) #Looks good
  acf(resid) #Also looks good 

#Plots----
  dat %>%
    ggplot(aes(DFA_Trend, LAT_COD, col=SIZESEX, shape=SIZESEX)) +
    geom_point() +
    geom_smooth(method="lm", se=F) + #raw data smoothed
    labs(y="Snow Crab Center of Distribution", x= "Environmental Trend") +
    scale_color_manual(labels = c("Immature Females", "Mature Females", "31 to 60mm Males", "61 to 90mm Males","91 to 120mm Males"),
                       values=cbPalette) +
    theme_bw() +theme(legend.title=element_blank())+
    theme(panel.grid = element_blank(),legend.position = "none")->g3
  
# Extract the prediction data frame
  pred.mm <- ggpredict(mod_final, terms = c("DFA_Trend"))  # this gives overall predictions for the model
  pred.mm
  
# Plot the model predicted fit  
  ggplot(pred.mm) + 
    geom_line(aes(x = x, y = predicted)) +    # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "lightgrey", alpha = 0.5) +  # error band
    geom_point(data = dat,                      # adding the raw data (centered values)
               aes(x = DFA_Trend, y = center.COD, colour = SIZESEX)) + 
    labs(y="Mean-centered Snow Crab \n Center of Distribution", x= "Environmental Trend") +
    scale_color_manual(labels = c("Immature Females", "Mature Females", "31 to 60mm Males", "61 to 90mm Males","91 to 120mm Males"),
                       values=cbPalette) +
    theme_bw() +theme(legend.title=element_blank())-> g4
  
##Plot for Fig 4 

## combine plots
  plot_grid(g2 + theme(legend.position="none"), 
            g4 + theme(legend.position="none"),
            align = 'vh',
            labels = c("a)", "b)")) -> ts_plot
  
  # extract a legend that is laid out horizontally
  legend_b <- get_legend(g2 + 
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  
  # add the legend underneath the plot
  plot_grid(ts_plot, legend_b, ncol = 1, rel_heights = c(4, .2))
  
## write plot
  ggsave(filename = "modelfits.png", device = "png", width = 10, height = 5, 
         dpi = 300)
