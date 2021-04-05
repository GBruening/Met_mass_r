library(ggplot2)
library(ggthemes)
library(tidyr)
library(multcomp)
library(lme4)
library(nlstools)
library(knitr)
library("RColorBrewer")
library('R.matlab')
library(boot)
library(viridis)
library(dplyr)

#====================== Load Fonts
# font_import()
# loadfonts(device='win')
windowsFonts(Times=windowsFont("TT Times New Roman"))
extrafont::loadfonts()

theme_set(theme_classic()+
            theme(axis.text         = element_text(color='black'),
                  axis.ticks        = element_line(color='black'),
                  legend.text       = element_text(size=15),
                  legend.text.align = 0,
                  legend.position   = 'right',
                  plot.title        = element_text(hjust = 0.5),
                  axis.line         = element_line(color='black')))

save_plots = 0

#====================== LOAD DATA ==============================
# Alaa Change the parent directory to your directory
parent_directory = 'D:/Google Drive/Preferred Mass/alaa code/Smoll_ut_mod'
setwd(parent_directory)

prefdata <-read.csv('2018Data-UnorderednoFML.csv')
metdata<-read.csv('Met_data_erik_test.csv')
mpdata <- read.csv('met_power_data.csv')
mpdata_rest <- readMat('baselinemetabolics.mat')$baseline
pilotdata <- read.csv('Pilot_Data-UnorderednoFML.csv')
smalltdata <- read.csv('DataArcT-UnorderednoFML.csv')

eff_masses_pref <- read.csv('eff_masses_pref.csv',header=0)
eff_masses_meta <- read.csv('eff_masses_meta.csv',header=0)

react_sim <-read.csv('react_sim.csv')
torque2_data <- read.csv('sum_torque2_only.csv')
torque2_mj_data <- read.csv('sum_torque2_scscript.csv')

# ===================== Add things =================================================
pilotdata$reacttoback <- (pilotdata$idxmoveback - pilotdata$idxonset)*.005
pilotdata$reaction_tanv <- .005*pilotdata$reaction_tanv
pilotdata$movedur <- pilotdata$reacttoback
pilotdata$peakvel_target <- pilotdata$peakvel_radv

prefdata$success = abs(prefdata$miss_dist)<.014
smalltdata$success = abs(smalltdata$maxex)<.11 & abs(smalltdata$missangle)<7

prefdata$reaction_tanv <- .005*prefdata$reaction_tanv
smalltdata$reaction_tanv <- .005*smalltdata$reaction_tanv

prefdata$reaction_erik <- 0.005*(prefdata$idxonsetErik-prefdata$idxtargetshow)
prefdata$reaction_vthresh <- 0.005*(prefdata$vthres_onset-prefdata$idxtargetshow)
prefdata$reaction_extrap <- 0.005*(prefdata$idxonset_extrap-prefdata$idxtargetshow)

smalltdata$reaction_erik <- 0.005*(smalltdata$idxonsetErik-smalltdata$idxtargetshow)
smalltdata$reaction_vthresh <- 0.005*(smalltdata$vthres_onset-smalltdata$idxtargetshow)
smalltdata$reaction_extrap <- 0.005*(smalltdata$idxonset_extrap-smalltdata$idxtargetshow)

pilotdata$reaction_erik <- 0.005*(pilotdata$idxonsetErik-pilotdata$idxtargetshow)
pilotdata$reaction_vthresh <- 0.005*(pilotdata$vthres_onset-pilotdata$idxtargetshow)
pilotdata$reaction_extrap <- 0.005*(pilotdata$idxonset_extrap-pilotdata$idxtargetshow)

prefdata$peakvel_target   <- prefdata$peakvel_target_vsign
smalltdata$peakvel_target <- smalltdata$peakvel_target_vsign
pilotdata$peakvel_target  <- pilotdata$peakvel_target_vsign

# ===================== Adding Eff_mass =======================
# Adding Column for estimated eff_mass to preferred experiment.
index <- c(0,3,5,8)
# values <- c(2.44,3.80,4.7,6.1)
values <- c(2.506,3.959,4.894,6.282)
subjects = c(1,2,3,4,5,6,7,8,9,10,11,12)
prefdata$eff_mass <- values[match(prefdata$condition,index)]
pilotdata$eff_mass <- values[match(pilotdata$condition,index)]
smalltdata$eff_mass <- values[match(smalltdata$condition,index)]

eff_mass = numeric(length(prefdata$movedur))
for (i in 1:length(prefdata$movedur)){
  eff_mass[i] = eff_masses_pref[match(prefdata$condition[i],index),match(prefdata$subj[i],subjects)]
}

prefdata$eff_mass2 <- eff_mass

#============================ Filtering Funciton =========================
filtering_func <- function(vars,data,dataname){
  for (cond in c(0,3,5,8)){
    filt_string = ""
    for (var in vars){
      eval(parse(text = paste(var,cond,' = boxplot.stats(filter(',dataname,',condition == ',cond,')$',var,')$out', sep = "")))
      eval(parse(text = paste(var,cond,'_high = min(',var,cond,'[',var,cond,'>median(filter(data,condition==',cond,')$',var,')])', sep = "")))
      
      lowstring <- paste('if (length(',var,cond,'[',var,cond,'<mean(filter(data,condition==',cond,')$',var,')])>0){
                                     ',var,cond,'_low = max(',var,cond,'[',var,cond,'<mean(filter(data,condition==',cond,')$',var,')])
                          } else{
                                     ',var,cond,'_low = 0
                          }',sep="")
      eval(parse(text = lowstring))
      filt_string <- paste(filt_string,',',var,'<',var,cond,'_high,',var,'>',var,cond,'_low', sep = "")
    }
    eval(parse(text = paste('a',cond,' <- filter(data,condition == ',cond,filt_string,')', sep = "")))
  }
  data <- rbind(a0,a3,a5,a8)
  return(data)
}

# ========================== Filtering ================================================
prefdata <- filtering_func(c('movedur','reaction_tanv','reaction_tanvel','miss_dist'),prefdata,'prefdata')
prefdata <- filter(prefdata,maxex<.14,odd_trial==1)

# pilotdata <- filter(pilotdata,reacttoback<3.0,reacttoback>.25,idxonset1<100)
pilotdata <- filter(pilotdata,trialtot>200)
pilotdata$miss_dist <- pilotdata$maxex
pilotdata <-filtering_func(c('movedur','reaction_tanv','reaction_tanvel','miss_dist'),pilotdata,'pilotdata')
pilotdata <- filter(pilotdata,maxex<.20)

smalltdata <- filtering_func(c('movedur','reaction_tanv','reaction_tanvel'),smalltdata,'smalltdata')
# smalltdata <- filter(pilotdata,max)

prefdata_factor <- prefdata
prefdata_factor$eff_mass <- as.factor(prefdata_factor$eff_mass)
prefdata_factor$condition <- as.factor(prefdata_factor$condition)
prefdata_factor$targetnum <- as.factor(prefdata_factor$targetnum)

smalltdata_factor <- smalltdata
smalltdata_factor$eff_mass <- as.factor(smalltdata_factor$eff_mass)
smalltdata_factor$condition <- as.factor(smalltdata_factor$condition)
smalltdata_factor$targetnum <- as.factor(smalltdata_factor$targetnum)

pilotdata_factor <- pilotdata
pilotdata_factor$eff_mass <- as.factor(pilotdata_factor$eff_mass)
pilotdata_factor$condition <- as.factor(pilotdata_factor$condition)
pilotdata_factor$targetnum <- as.factor(pilotdata_factor$targetnum)

met_orig_length = length(metdata$movedur)

metdata = filter(metdata,odd_trial==1)
met_length = length(metdata$movedur)
metdata <- filter(metdata,miss_dist<.1,movedur>.2,movedur<2,maxex<.2,reaction_tanv<.5,abs(missangle)<50)

index <- c(1,2,3,4)
values <- c(2.44,4.830,7.13,11.69)
metdata$eff_mass <- values[match(metdata$condition,index)]
metdata = cbind(metdata,metdata$trial%%2==0)
colnames(metdata) = append(colnames(metdata)[1:length(colnames(metdata))-1],'outward')
#metdata2 = metdata
metdata = filter(metdata,odd_trial==1)

# ========================== Adding normalized ================================================
add_norm <- function(data, variables){
  for (vari in variables){
    eval(parse(text = paste('zero_vals = aggregate(',vari,'~ subj + condition,data,mean)',sep = '')))
    # eval(parse(text = paste('zero_vals = aggregate(',vari,'~ condition,data,mean)',sep = '')))
    zero_vals = zero_vals[zero_vals$condition == '0',]
    for (k in 1:length(data$subj)){
      # data[k,paste(vari,'_norm',sep='')] = data[k,vari]/zero_vals[,vari]
      data[k,paste(vari,'_norm',sep='')] = data[k,vari]/zero_vals[zero_vals$subj == data[k,'subj'],vari]
    }
  }
  return(data)
}

prefdata <- add_norm(prefdata,c('movedur','peakvel_target','reaction_tanv','miss_dist'))
smalltdata <- add_norm(smalltdata,c('movedur','peakvel_target','reaction_tanv','miss_dist'))
pilotdata <- add_norm(pilotdata,c('movedur','peakvel_target','reaction_tanv','miss_dist'))

# ========================== Adding Metabolics things =========================================

index <- c(2.47,4.73,6.99,11.50)
values <- c(1,2,3,4)
mpdata$condition <- values[match(mpdata$effmass,index)]

index <- c(2.47,4.73,6.99,11.50)
values <- c(2.44,4.830,7.13,11.69)
mpdata$effmass <- values[match(mpdata$effmass,index)]

# Adding Column for estimated eff_mass to preferred experiment.
index <- c(1,2,3,4)
values <- c(2.44,4.830,7.13,11.69)
subjects = c(1,2,3,4,5,6,7,8)

eff_mass = numeric(length(metdata$subj))
for (i in 1:length(metdata$subj)){
  eff_mass[i] = eff_masses_meta[match(metdata$condition[i],index),match(metdata$subj[i],subjects)]
}

metdata$eff_mass2 <- eff_mass
metdata$effmass = metdata$eff_mass

index <- unique(mpdata$effmass)
eff_mass = numeric(length(mpdata$subj))
for (i in 1:length(mpdata$subject)){
  eff_mass[i] = eff_masses_meta[match(mpdata$effmass[i],index),match(mpdata$subj[i],subjects)]
  # mpdata$metpowergross[i] = mpdata$metpowernet[i]+mpdata_rest[mpdata$subject[i],mpdata$condition[i]]
  mpdata$metpowerrest[i] = mpdata_rest[mpdata$subject[i],mpdata$condition[i]]
}

mpdata$effmass2 <- eff_mass

# ========================== Adding Variance Metrics =========================================

stdev_trial_plot <- function(data, variable, add_or_plot){
  if (add_or_plot=='add'){
    max_trial = max(data$trial)
    
    for (k in 11:(length(data$subj)-10)){
      if (var(data[(k-10):(k+10),'eff_mass'])==0){
        if ('speed' %in% colnames(data) && var(data[(k-10):(k+10),'speed'])==0){
          data[k,paste(variable,'_sd',sep='')] = sd(data[(k-10):(k+10),variable])
        } else {
          data[k,paste(variable,'_sd',sep='')] = sd(data[(k-10):(k+10),variable])
        }
      } else {
        data[k,paste(variable,'_sd',sep='')] = NaN
      }
    }
    
    return(data)
  }
  if (add_or_plot=='plot'){
    eval(parse(text = paste('data$yvar = data$',variable,'_sd',sep='')))
    g<-ggplot(data=data)+
      geom_smooth(aes(x=trial,
                      y=yvar,
                      color = factor(subj)))+
      # theme_classic()
      labs(x = 'Trial Number', y = paste('StDev (',variable,')',sep=''),color = 'Subject')
    
    if ('speed' %in% colnames(data)){
      g<-g+facet_grid(rows = vars(eff_mass),cols = vars(speed))
    } else {
      g<-g+facet_grid(rows = vars(eff_mass))
    }
    return(g)
  }
}

average_trial_plot <- function(data, variable, add_or_plot){
  if (add_or_plot=='add'){
    max_trial = max(data$trial)
    
    for (k in 11:(length(data$subj)-10)){
      if (var(data[(k-10):(k+10),'eff_mass'])==0){
        if ('speed' %in% colnames(data) && var(data[(k-10):(k+10),'speed'])==0){
          data[k,paste(variable,'_avg',sep='')] = mean(data[(k-10):(k+10),variable])
        } else {
          data[k,paste(variable,'_avg',sep='')] = mean(data[(k-10):(k+10),variable])
        }
      } else {
        data[k,paste(variable,'_avg',sep='')] = NaN
      }
    }
    
    return(data)
  }
  if (add_or_plot=='plot'){
    eval(parse(text = paste('data$yvar = data$',variable,'_avg',sep='')))
    g<-ggplot(data=data)+
      geom_smooth(aes(x=trial,
                      y=yvar,
                      color = factor(subj)))+
      # theme_classic()
      labs(x = 'Trial Number', y = paste('Average (',variable,')',sep=''),color = 'Subject')
    
    if ('speed' %in% colnames(data)){
      g<-g+facet_grid(rows = vars(eff_mass),cols = vars(speed))
    } else {
      g<-g+facet_grid(rows = vars(eff_mass))
    }
    return(g)
  }
}


for (variable in c('miss_dist','miss_rad','missangle')){
  prefdata <- stdev_trial_plot(prefdata,variable,'add')
  smalltdata <- stdev_trial_plot(smalltdata,variable,'add')
  pilotdata <- stdev_trial_plot(pilotdata,variable,'add')
  metdata <- stdev_trial_plot(metdata,variable,'add')
  
  prefdata <- average_trial_plot(prefdata,variable,'add')
  smalltdata <- average_trial_plot(smalltdata,variable,'add')
  pilotdata <- average_trial_plot(pilotdata,variable,'add')
  metdata <- average_trial_plot(metdata,variable,'add')
}


# =========================== Generate the df for 3 experiments ===================
prefdata$exp <- 'pref'
pilotdata$exp <- 'pilot'
smalltdata$exp <- 'smallt'
smalltdata$subj <- smalltdata$subj+12
pilotdata$subj <- pilotdata$subj+12+12

allexp_plot <- function(vars,data,...){
  str = paste('a=as.data.frame(cbind(',data,'$subj,',data,'$eff_mass,',data,'$condition,',data,'$exp',sep='')
  for (var in vars){
    # print(var)
    str = paste(str,',',data,'$',var,sep='')
  }
  str = paste(str,'))',sep='')
  a= eval(parse(text=str))
  colnames(a) = c('subj','eff_mass','condition','exp',vars)
  a$eff_mass = as.numeric(as.character(a$eff_mass))
  a$condition = as.numeric(as.character(a$condition))
  for (var in vars){
    eval(parse(text = paste('a$',var,' = as.numeric(as.character(a$',var,'))',sep='')))
  }
  return(a)
}

vars = c('targetnum','movedur','movedur_norm',
         'peakvel_target','peakvel_target_norm','idxpeakv',
         'miss_dist','miss_dist_norm',
         'reaction_tanv','reaction_tanv_norm','reaction_tanvel',
         'pathltar','missangle','miss_rad','trial',
         'miss_dist_avg','miss_dist_sd',
         'missangle_avg','missangle_sd',
         'miss_rad_avg','miss_rad_sd')
a=allexp_plot(vars,'prefdata')
b=allexp_plot(vars,'pilotdata')
c=allexp_plot(vars,'smalltdata')

b$peakvel_target = pilotdata$peakvel_radv
# colnames(b)[colnames(b)=='peakvel_target'] <- 'peavel_radv'

prefpilot=data.frame(rbind(a,b,c))
prefpilot$exp <- factor(prefpilot$exp, levels = c('pref','smallt','pilot'))
prefpilot = filter(prefpilot,movedur<1.75)
prefpilot$subj <- as.numeric(prefpilot$subj)

# ===================== Setting plots sizes ===================
# Pref stuff
vplot_point_size = 0
vplot_avg_point_size = 3
v_plot_line_size = 1
v_plot_err_bar_size = 1

vel_plot_line_size = 1
vel_plot_point_size = 3

err_plot_point_size = 1
err_plot_line_size = 1

react_plot_point_size = 1

# Met stuff
met_point_size = 3
met_costmin_size = met_point_size-2
met_line_size = 1

met_errorbar_size = .5

# ===================== Fucntions ===================
#Another stat summary function
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


#Fischer Method Function
p_fisch_func <- function(var,data){
  p_fisch = c()
  chi_fisch = 0
  for (subjnum in unique(data$subj)){
    str = paste('p_fisch[subjnum] = anova(lm(',var,'~ targetnum + eff_mass2 + targetnum*eff_mass2,
             data=filter(data,subj == subjnum)))$`Pr(>F)`[2]',sep='')
    eval(parse(text = str))
    chi_fisch = chi_fisch + log(p_fisch[subjnum])
  }
  chi_fisch = -2*chi_fisch
  return(pchisq(chi_fisch, df=subjnum*2, lower.tail=FALSE))
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

number_ticks <- function(n) {function(limits) pretty(limits, n)}

fix_plot_theme <- function(plot_to_fix,font){
  string = paste(plot_to_fix,'<-',plot_to_fix,'+
                   # theme_classic()
                   theme(text              = element_text(family=',font,',color=\'black\'),
                   axis.text         = element_text(color=\'black\'),
                   axis.ticks        = element_line(color=\'black\'),
                   plot.title        = element_text(hjust = 0.5),
                   axis.line         = element_line(color=\'black\'),
                   legend.position   = \'none\')',sep = '')
  eval(parse(text = string))
  return(plot_to_fix)
}


color1 = c(255,255,0)/255
color2 = c(255,196,0)/255
color3 = c(216,117,1)/255
color4 = c(139,0,0)/255

mass_colors <- c(rgb(color1[1],color1[2],color1[3]),
                 rgb(color2[1],color2[2],color2[3]),
                 rgb(color3[1],color3[2],color3[3]),
                 rgb(color4[1],color4[2],color4[3]))

# ===================== Met Additions ===================


index = 1:8
values = c(125,156.199,194.7,129.8,173.8,159.7,136.4,127.6)
metdata$subjmass <- values[match(metdata$subj,index)]
mpdata$subjmass <- values[match(mpdata$subj,index)]

metdata_factor = metdata
metdata_factor$condition <- as.factor(metdata_factor$condition)

mpdata$netcost <- mpdata$movedur * mpdata$metpowernet
mpdata$grosscost <- mpdata$movedur * mpdata$metpowergross

mpdata <- filter(mpdata,metpowernet>0)
mpdata_factor <- mpdata
mpdata_factor$condition <- as.factor(mpdata_factor$condition)

tab = matrix(,nrow=4,ncol=1)
dcount = 0
for (data in c('pref')){#,'smallt','pilot')){
  dcount = dcount +1
  data_set = metdata
  ccount = 0
  for (c in c(1,2,3,4)){
    ccount = ccount + 1
    tab[ccount,dcount] = paste(mean(filter(data_set,condition == c)$eff_mass2),
                               '±',
                               sd(filter(data_set,condition == c)$eff_mass2)/sqrt(12))
  }
}


# Movement Duration =======
move_pref <- aggregate(movedur ~ eff_mass,prefdata,mean)$movedur
move_pilo <- aggregate(movedur ~ eff_mass,pilotdata,mean)$movedur
move_arct <- aggregate(movedur ~ eff_mass,smalltdata,mean)$movedur

move_pref_se <- aggregate(movedur ~ eff_mass,prefdata,sd)$movedur/sqrt(12)
move_pilo_se <- aggregate(movedur ~ eff_mass,pilotdata,sd)$movedur/sqrt(12)
move_arct_se <- aggregate(movedur ~ eff_mass,smalltdata,sd)$movedur/sqrt(18)

move_dur_avg <- round(cbind(move_pref,move_arct,move_pilo),digits=3)
move_dur_se <- round(cbind(move_pref_se,move_arct_se,move_pilo_se),digits=3)
move_dur_df = matrix(,nrow = 4, ncol = 3)
for (rows in c(1:4)){
  for (cols in c(1:3)){
    move_dur_df[rows,cols] = paste(move_dur_avg[rows,cols],'±',move_dur_se[rows,cols],sep='')
  }
}

move_dur <- cbind(c(2.5,3.8,4.7,6.1),move_dur_df)
colnames(move_dur) <- c('Effective Mass (kg)','2a','2b','2c')

# MP Gross ===============

MPGross_model=nls(metpowergross ~ a1 + 100*a2*(effmass2^a3)/(movedur^a4),
                  data=mpdata,start=list(a1=100,a2=.1,a3=1,a4=1))
modelsum_mpgross = summary(MPGross_model)
a_mpgross=coef(MPGross_model)[1]
b_mpgross=coef(MPGross_model)[2]*100
c_mpgross=coef(MPGross_model)[3]
d_mpgross=coef(MPGross_model)[4]

fun.1 <- function(t) a_mpgross+b_mpgross*(2.44^c_mpgross)/(t^d_mpgross)
fun.2 <- function(t) a_mpgross+b_mpgross*(4.83^c_mpgross)/(t^d_mpgross)
fun.3 <- function(t) a_mpgross+b_mpgross*(7.13^c_mpgross)/(t^d_mpgross)
fun.4 <- function(t) a_mpgross+b_mpgross*(11.69^c_mpgross)/(t^d_mpgross)

# MP Net ================
MPNet_model=nls(metpowernet ~ a1 + a2*(effmass2^a3)/(movedur^a4), data=mpdata,start=list(a1=1,a2=1,a3=1,a4=1))
modelsum_mpnet = summary(MPNet_model)

# plot(nlsResiduals(MPNet_model))

# a_mpnet=coef(MPNet_model)[1]
# a_mpnet=coef(MPNet_model)[1]
# b_mpnet=coef(MPNet_model)[2]
# c_mpnet=coef(MPNet_model)[3]
# d_mpnet=coef(MPNet_model)[4]

a_mpnet=coef(MPGross_model)[1]-mean(mpdata_rest)
b_mpnet=coef(MPNet_model)[2]
c_mpnet=coef(MPNet_model)[3]
d_mpnet=coef(MPNet_model)[4]

# Sum of Torque ================
values <- c(2.44,4.830,7.13,11.69)
index <- c(1,2,3,4)
torque2_data$effmass <- values[match(torque2_data$c,index)]
torque2_data$distance <- .1
torque2_data$sum_t2_rate <- torque2_data$sum_t2/torque2_data$movedur

a = 0

torque2_model=nls(sum_t2_rate ~ a + a1*(effmass^a2)/(movedur^a3),
                  data=torque2_data,
                  start=list(a1=1,a2=1,a3=1),
                  control = list(maxiter = 500))
modelsum_torque2 = summary(torque2_model)

b_torque=coef(torque2_model)[1]
c_torque=coef(torque2_model)[2]
d_torque=coef(torque2_model)[3]

fun.1 <- function(t) b_torque*(2.44^c_torque)/(t^(d_torque))
fun.2 <- function(t) b_torque*(4.830^c_torque)/(t^(d_torque))
fun.3 <- function(t) b_torque*(7.130^c_torque)/(t^(d_torque))
fun.4 <- function(t) b_torque*(11.69^c_torque)/(t^(d_torque))


# MCost Gross
MCostmin_dur = c(optimize(fun.1,interval=c(0,2))$minimum,
                 optimize(fun.2,interval=c(0,2))$minimum,
                 optimize(fun.3,interval=c(0,2))$minimum,
                 optimize(fun.4,interval=c(0,4))$minimum)
MCostmin_obj = c(optimize(fun.1,interval=c(0,2))$objective,
                 optimize(fun.2,interval=c(0,2))$objective,
                 optimize(fun.3,interval=c(0,2))$objective,
                 optimize(fun.4,interval=c(0,4))$objective)
temp = cbind(c(2.44,4.830,7.13,11.69),MCostmin_dur,MCostmin_obj)
colnames(temp) = c('effmass','MCostmin_dur','MCostmin_obj')


preffun.1 <- function(t) a_mpgross*t+b_mpgross*(2.44^c_mpgross)/(t^(d_mpgross-1))
preffun.2 <- function(t) a_mpgross*t+b_mpgross*(3.8^c_mpgross)/(t^(d_mpgross-1))
preffun.3 <- function(t) a_mpgross*t+b_mpgross*(4.7^c_mpgross)/(t^(d_mpgross-1))
preffun.4 <- function(t) a_mpgross*t+b_mpgross*(6.1^c_mpgross)/(t^(d_mpgross-1))

MCostmin_dur = c(optimize(preffun.1,interval=c(0,2))$minimum,
                 optimize(preffun.2,interval=c(0,2))$minimum,
                 optimize(preffun.3,interval=c(0,2))$minimum,
                 optimize(preffun.4,interval=c(0,4))$minimum)
MCostmin_obj = c(optimize(preffun.1,interval=c(0,2))$objective,
                 optimize(preffun.2,interval=c(0,2))$objective,
                 optimize(preffun.3,interval=c(0,2))$objective,
                 optimize(preffun.4,interval=c(0,4))$objective)
temppref = cbind(c(2.441,3.8,4.7,6.1),MCostmin_obj,MCostmin_dur,aggregate(movedur ~ eff_mass,prefdata,mean)$movedur)
colnames(temppref) = c('effmass','MCostmin_obj','MCostgross_min_dur','pref_avg_dur')

# Probability Functions ==========
# Alaa
pref_p_glm  = glm(success ~ movedur,
                  data=prefdata,
                  family=binomial(link="logit"))
prefdata$p_success = inv.logit(predict(pref_p_glm))

smallt_p_glm  = glm(success ~ movedur,
                    data=smalltdata,
                    family=binomial(link="logit"))
smalltdata$p_success = inv.logit(predict(smallt_p_glm))

# Metabolics Data
metdata$success_2a = abs(metdata$miss_dist)<.014
met_p_glm  = glm(success_2a ~ movedur+eff_mass,
                 data=metdata,
                 family=binomial(link="logit"))
metdata$p_success_2a = inv.logit(predict(met_p_glm))

metdata$success_2b = abs(metdata$maxex)<.110 & abs(metdata$missangle)<7
met_smallt_p_glm  = glm(success_2b ~ movedur+eff_mass,
                        data=metdata,
                        family=binomial(link="logit"))
metdata$p_success_2b = inv.logit(predict(met_smallt_p_glm))


# Utility ======================

# pref_glm_use = pref_p_glm
pref_glm_use = met_p_glm

# smallt_glm_use = smallt_p_glm
smallt_glm_use = met_smallt_p_glm

############ Modeling All the different Models

fun_optim_prx <- function(alpha,masses,mvttimes,p_alpha,rxtimes){
  err=0
  dur = rep(0,length(mvttimes))
  for (k in 1:length(mvttimes)){
    t = mvttimes[k]
    m = masses[k]
    p = p_alpha[k]
    rt = rxtimes[k]
    count = 1
    times = seq(0.1,5,0.01)
    ut = rep(0,length(times))
    for (t in times){
      ut[count] = (alpha*p-a0*rt-(a*t+b*(m^c)/(t^(d-1))))/(rt+t)
      count = count+1
    }
    dur[k] = times[which.max(ut)]
  }
  err = (mvttimes-dur)^2
  return(sum(err))
}

fun_dur_prx <- function(alpha,masses,mvttimes,p_alpha,rxtimes){
  err=0
  dur = rep(0,length(mvttimes))
  for (k in 1:length(mvttimes)){
    t = mvttimes[k]
    m = masses[k]
    p = p_alpha[k]
    rt = rxtimes[k]
    count = 1
    times = seq(0.2,5,0.01)
    ut = rep(0,length(times))
    for (t in times){
      ut[count] = (alpha*p-a0*rt-(a*t+b*(m^c)/(t^(d-1))))/(rt+t)
      count = count+1
    }
    dur[k] = times[which.max(ut)]
  }
  return(dur)
}

#================  Gross Metabolics ==========================
a=coef(MPGross_model)[1]
b=coef(MPGross_model)[2]*100
c=coef(MPGross_model)[3]
d=coef(MPGross_model)[4]

fun.1 <- function(t) a*t+b*(2.44^c)/(t^(d-1))
fun.2 <- function(t) a*t+b*(3.8^c)/(t^(d-1))
fun.3 <- function(t) a*t+b*(4.7^c)/(t^(d-1))
fun.4 <- function(t) a*t+b*(6.1^c)/(t^(d-1))

met_gross_dur = c(optimize(fun.1,interval=c(0,2))$minimum,
                  optimize(fun.2,interval=c(0,2))$minimum,
                  optimize(fun.3,interval=c(0,2))$minimum,
                  optimize(fun.4,interval=c(0,4))$minimum)

# ===================== Utility ===================
masses = unique(prefdata$eff_mass)
mvttimes_pref = aggregate(movedur~condition,prefdata,mean)$movedur
mvttimes_pref_se = aggregate(movedur~condition,prefdata,sd)$movedur/sqrt(8)
miss_dist_sd_pref = aggregate(miss_dist ~ eff_mass, prefdata,'sd')$miss_dist
rxtimes_pref = aggregate(reaction_tanv~condition,prefdata,mean)$reaction_tanv
p_alpha_pref = inv.logit(predict(pref_glm_use,data.frame(movedur=mvttimes_pref,eff_mass=masses)))

mvttimes_smallt = aggregate(movedur~condition,smalltdata,mean)$movedur
mvttimes_smallt_se = aggregate(movedur~condition,smalltdata,sd)$movedur/sqrt(8)
miss_dist_sd_smallt = aggregate(miss_dist ~ eff_mass, smalltdata,'sd')$miss_dist
rxtimes_smallt = aggregate(reaction_tanv~eff_mass,smalltdata,mean)$reaction_tanv
p_alpha_smallt_gross = inv.logit(predict(smallt_glm_use,data.frame(movedur=mvttimes_smallt,eff_mass=masses)))

mvttimes_pilot = aggregate(movedur~condition,pilotdata,mean)$movedur
mvttimes_pilot_se = aggregate(movedur~condition,pilotdata,sd)$movedur/sqrt(8)
miss_angle_sd_pilot = aggregate(absmissangle ~ eff_mass, pilotdata,'sd')$absmissangle
rxtimes_pilot = aggregate(reaction_tanv~eff_mass,smalltdata,mean)$reaction_tanv

fun_optim_prx <- function(alpha,masses,mvttimes,p_alpha,rxtimes,opt_or_dur){
  err=0
  dur = rep(0,length(mvttimes))
  times = seq(0.3,2,0.001)
  ut = matrix(,nrow = length(mvttimes),ncol = length(times))
  for (k in 1:length(mvttimes)){
    t = mvttimes[k]
    m = masses[k]
    p = p_alpha[k]
    rt = rxtimes[k]
    count = 1
    for (t in times){
      ut[k, count] = (alpha*p-a0*rt-(a*t+b*(m^c)/(t^(d-1))))/(rt+t)
      count = count+1
    }
    dur[k] = times[which.max(ut[k,])]
  }
  err = (mvttimes-dur)^2
  if (opt_or_dur == 'optimize'){
    return(sum(err))
  } else if (opt_or_dur == 'duration'){
    return(dur)
  }
}

fun_optim_prx_logit <- function(alpha,masses,mvttimes,p_alpha,rxtimes,opt_or_dur){
  err=0
  dur = rep(0,length(mvttimes))
  for (k in 1:length(mvttimes)){
    t = mvttimes[k]
    m = masses[k]
    p = inv.logit(predict(p_alpha,data.frame(movedur=t,eff_mass=m)))
    rt = rxtimes[k]
    times = seq(0.5,2,0.001)
    ut = rep(0,length(times))
    probs = inv.logit(predict(p_alpha,data.frame(movedur=times,eff_mass=m)))
    count = 1
    for (t in times){
      ut[count] = (alpha*probs[count]-a0*rt-(a*t+b*(m^c)/(t^(d-1))))/(rt+t)
      count = count+1
    }
    dur[k] = times[which.max(ut)]
    # print(max(ut))
  }
  err = (mvttimes-dur)^2
  if (opt_or_dur == 'optimize'){
    return(sum(err))
  } else if (opt_or_dur == 'duration'){
    return(dur)
  }
}

fun_optim_prx_comb <- function(alpha,masses,mvttimes,p_alpha,rxtimes,opt_or_dur){
  err=0
  dur = rep(0,length(mvttimes))
  for (k in 1:length(mvttimes)){
    t = mvttimes[k]
    m = masses[k]
    p = p_alpha[k]
    rt = rxtimes[k]
    count = 1
    times = seq(0.1,2,0.001)
    ut = rep(0,length(times))
    for (t in times){
      ut[count] = (alpha*p-a0*rt-(a*t+b*(m^c)/(t^(d-1))))/(rt+t)
      count = count+1
    }
    dur[k] = times[which.max(ut)]
  }
  err = (mvttimes-dur)^2
  if (opt_or_dur == 'optimize'){
    return(sum(err))
  } else if (opt_or_dur == 'duration'){
    return(dur)
  }
}

fun_optim_prx_logit_comb <- function(alpha,masses,mvttimes,glm1,glm2,rxtimes,opt_or_dur){
  err=0
  dur = rep(0,length(mvttimes))
  for (k in 1:length(mvttimes)){
    t = mvttimes[k]
    m = masses[k]
    rt = rxtimes[k]
    count = 1
    times = seq(0.1,2,0.001)
    ut = rep(0,length(times))
    if (k<5){
      probs = inv.logit(predict(glm1,data.frame(movedur=times,eff_mass=m)))
    } else {
      probs = inv.logit(predict(glm2,data.frame(movedur=times,eff_mass=m)))
    }
    t_count = 0
    for (t in times){
      t_count = t_count+1
      ut[count] = (alpha*probs[t_count]-a0*rt-(a*t+b*(m^c)/(t^(d-1))))/(rt+t)
      count = count+1
    }
    dur[k] = times[which.max(ut)]
  }
  err = (mvttimes-dur)^2
  if (opt_or_dur == 'optimize'){
    return(sum(err))
  } else if (opt_or_dur == 'duration'){
    return(dur)
  }
}

fun_p_alpha <- function(miss_dist, min_dist){
  p_alpha = c(0,0,0,0)
  count = 1
  for (k in miss_dist){
    p_alpha[count] = 2*(pnorm(min_dist, mean = 0, sd = k)-.5)
    count = count+1
  }
  return(p_alpha)
}
p_alpha_pilot_gross = fun_p_alpha(miss_angle_sd_pilot, 90)


# COMBINED UTILITY Torque Squared ===========================
# USING logit function
a0 = mean(mpdata$metpowerrest)
a0_se = sd(mpdata$metpowerrest)/sqrt(8)
a=coef(torque2_model)[1]
b=coef(torque2_model)[2]
c=coef(torque2_model)[3]
d=coef(torque2_model)[4]
alpha_torque_comb_gross=optimize(fun_optim_prx_logit_comb,
                                 masses     = c(masses,masses),
                                 mvttimes   = c(mvttimes_pref,mvttimes_smallt),
                                 glm1       = pref_glm_use,
                                 glm2       = smallt_glm_use,
                                 rxtimes    = c(rxtimes_pref,rxtimes_smallt),
                                 opt_or_dur = 'optimize',
                                 interval   = c(-100,200))

util_torque_pref_alpha_comb_gross = as.numeric(fun_optim_prx_logit_comb(
  alpha      = alpha_torque_comb_gross$minimum,
  masses     = c(masses,masses),
  mvttimes   = c(mvttimes_pref,mvttimes_smallt),
  glm1       = pref_glm_use,
  glm2       = smallt_glm_use,
  opt_or_dur = 'duration',
  rxtimes    = c(rxtimes_pref,rxtimes_smallt)))

ut1_torque_alpha_2a2b = alpha_torque_comb_gross
ut1_torque_dur_2a2b   = util_torque_pref_alpha_comb_gross

ut1_torque_dur_2a2b2c = as.numeric(fun_optim_prx_comb(
  alpha_torque_comb_gross$minimum,
  masses = masses,
  mvttimes = mvttimes_pilot,
  p_alpha  = c(1,1,1,1),
  opt_or_dur = 'duration',
  rxtimes = rxtimes_pilot))

# COMBINED UTILITY Gross ===========================
# USING logit function
a=coef(MPGross_model)[1]
b=coef(MPGross_model)[2]*100
c=coef(MPGross_model)[3]
d=coef(MPGross_model)[4]
alpha_pref_comb_gross=optimize(fun_optim_prx_logit_comb,
                               masses     = c(masses,masses),
                               mvttimes   = c(mvttimes_pref,mvttimes_smallt),
                               glm1       = pref_glm_use,
                               glm2       = smallt_glm_use,
                               rxtimes    = c(rxtimes_pref,rxtimes_smallt),
                               opt_or_dur = 'optimize',
                               interval   = c(-100,200))

# Alaa
util_dur_pref_alpha_comb_gross = as.numeric(fun_optim_prx_logit_comb(
  alpha      = alpha_pref_comb_gross$minimum,
  masses     = c(masses,masses),
  mvttimes   = c(mvttimes_pref,mvttimes_smallt),
  glm1       = pref_glm_use,
  glm2       = smallt_glm_use,
  opt_or_dur = 'duration',
  rxtimes    = c(rxtimes_pref,rxtimes_smallt)))

ut1_grs_alpha_2a2b = alpha_pref_comb_gross
ut1_grs_dur_2a2b   = util_dur_pref_alpha_comb_gross

util_dur_pilot_gross_alpha_comb_gross = as.numeric(fun_optim_prx_comb(
  alpha_pref_comb_gross$minimum,
  masses = masses,
  mvttimes = mvttimes_pilot,
  p_alpha  = c(1,1,1,1),
  opt_or_dur = 'duration',
  rxtimes = rxtimes_pilot))

ut1_grs_dur_2a2b2c = as.numeric(fun_optim_prx(
  alpha    = alpha_pref_comb_gross$minimum,
  masses   = masses,
  mvttimes = mvttimes_pilot,
  p_alpha  = p_alpha_pilot_gross,
  rxtimes  = rxtimes_pilot,
  opt_or_dur = 'duration'))

# Uitlity for just 2c
a0 = mean(mpdata$metpowerrest)
a0_se = sd(mpdata$metpowerrest)/sqrt(8)
a=coef(MPGross_model)[1]
b=coef(MPGross_model)[2]*100
c=coef(MPGross_model)[3]
d=coef(MPGross_model)[4]

masses = unique(prefdata$eff_mass)

alpha_pilot_gross=optimize(fun_optim_prx,
                           masses   = masses,
                           mvttimes = mvttimes_pilot,
                           p_alpha  = p_alpha_pilot_gross,
                           rxtimes  = rxtimes_pilot,
                           opt_or_dur = 'optimize',
                           interval = c(-100,1000))

util_dur_pilot_gross = as.numeric(fun_optim_prx(alpha_pilot_gross$minimum,
                                                masses   = masses,
                                                mvttimes = mvttimes_pilot,
                                                p_alpha  = p_alpha_pilot_gross,
                                                rxtimes  = rxtimes_pilot,
                                                opt_or_dur = 'duration'))

ut1_grs_alpha_2c = alpha_pilot_gross
ut1_grs_dur_2c   = util_dur_pilot_gross

# COMBINED UTILITY Net  ===========================
a=coef(MPNet_model)[1]
b=coef(MPNet_model)[2]
c=coef(MPNet_model)[3]
d=coef(MPNet_model)[4]


alpha_pref_comb_net=optimize(fun_optim_prx_logit_comb,
                             masses     = c(masses,masses),
                             mvttimes   = c(mvttimes_pref,mvttimes_smallt),
                             glm1       = pref_glm_use,
                             glm2       = smallt_glm_use,
                             rxtimes    = c(rxtimes_pref,rxtimes_smallt),
                             opt_or_dur = 'optimize',
                             interval   = c(-100,200))

util_dur_pref_alpha_comb_net = as.numeric(fun_optim_prx_logit_comb(alpha_pref_comb_net$minimum,
                                                                   masses     = c(masses,masses),
                                                                   mvttimes   = c(mvttimes_pref,mvttimes_smallt),
                                                                   glm1       = pref_glm_use,
                                                                   glm2       = smallt_glm_use,
                                                                   opt_or_dur = 'duration',
                                                                   rxtimes    = c(rxtimes_pref,rxtimes_smallt)))

ut1_net_alpha_2a2b = alpha_pref_comb_net
ut1_net_dur_2a2b   = util_dur_pref_alpha_comb_net

ut1_net_dur_2a2b2c = as.numeric(fun_optim_prx(
  alpha    = alpha_pref_comb_net$minimum,
  masses   = masses,
  mvttimes = mvttimes_pilot,
  p_alpha  = p_alpha_pilot_gross,
  rxtimes  = rxtimes_pilot,
  opt_or_dur = 'duration'))

# Plotting it ======
control_times = rbind(cbind(rep('smallt',4),rep('data',4),unique(prefdata$condition)/2.2,mvttimes_smallt,rep(3,4),mvttimes_smallt_se),
                      cbind(rep('pref',4),rep('data',4),unique(prefdata$condition)/2.2,mvttimes_pref,rep(1,4),mvttimes_pref_se),
                      cbind(rep('pilot',4),rep('data',4),unique(prefdata$condition)/2.2,mvttimes_pilot,rep(2,4),mvttimes_pilot_se),
                      
                      cbind(rep('smallt',4),rep('model',4),unique(prefdata$condition)/2.2,util_dur_pref_alpha_comb_gross[5:8],rep(3,4),rep(0,4)),
                      cbind(rep('pref',4),rep('model',4),unique(prefdata$condition)/2.2,util_dur_pref_alpha_comb_gross[1:4],rep(1,4),rep(0,4)),
                      cbind(rep('pilot',4),rep('model',4),unique(prefdata$condition)/2.2,util_dur_pilot_gross,rep(2,4),rep(0,4)))
# cbind(rep('pilot',4),rep('model',4),c(2.5,3.8,4.7,6.1),util_dur_pilot_gross_alpha_comb_gross,rep(2,4),rep(0,4)))
colnames(control_times) = c('exp','datatype','eff_mass','movedur','expnum','movedur_se')

control_times = as.data.frame(control_times)
control_times$movedur = as.numeric(as.character(control_times$movedur))
control_times$movedur_se = as.numeric(as.character(control_times$movedur_se))
control_times$eff_mass = as.numeric(as.character(control_times$eff_mass))

utilfits_by_exp_2alphas <- ggplot()+
  geom_errorbar(data=control_times,
                aes(x=eff_mass,
                    ymin=movedur-movedur_se,
                    ymax=movedur+movedur_se,
                    color=expnum,
                    alpha=datatype),
                size=3,
                width=.2)+
  geom_line(data=control_times,
            aes(x=eff_mass,
                y=movedur,
                linetype=datatype,
                color=expnum,
                alpha = datatype),
            size=3)+
  geom_line(data=data.frame(cbind(unique(control_times$eff_mass),
                                  util_dur_pilot_gross_alpha_comb_gross) %>% `colnames<-`(c('eff_mass','movedur'))),
            aes(x=eff_mass,
                y=movedur),
            linetype = 'dashed',
            color = 'black',
            size=3)+
  geom_line(data=data.frame(cbind(unique(control_times$eff_mass),
                                  met_gross_dur) %>% `colnames<-`(c('eff_mass','movedur'))),
            aes(x=eff_mass,
                y=movedur),
            linetype = 'dashed',
            color    = 'red',
            size=3)+
  scale_linetype_manual(values=c('solid','dashed'),
                        labels=c('Data (Solid)','Model (Dashed)'))+
  scale_alpha_discrete(range=c(.5,1),
                       labels=c('Data (Solid)','Model (Dashed)'))+
  labs(x='Added mass (kg)',
       y='Movement Duration (s)',
       title = paste('Util fits by Experiment\n2a,2b \U1D6FC = ',
                     round(alpha_pref_comb_gross$minimum,3),
                     ', 2c \U1D6FC = ',
                     round(alpha_pilot_gross$minimum,3),
                     '\nBlack = prob equal 1, red = metgross min',sep=''))+
  scale_color_manual(labels = c('Circle\n(Normal)',
                                'None\n(Large)',
                                'Arc\n(Small)',
                                paste('\U1D6FC = ',round(alpha_pref_comb_gross$minimum,3),sep=''),
                                'Gross Met Min'),
                     values = c(gg_color_hue(3),
                                'black',
                                'red'))+
  scale_x_continuous(breaks = c(unique(control_times$eff_mass)))+
  theme_classic()+
  theme(text              = element_text(size=20,color='black'),
        axis.text         = element_text(color='black'),
        axis.ticks        = element_line(color='black'),
        legend.text       = element_text(size=15),
        legend.text.align = 0,
        legend.position   = 'right',
        plot.title        = element_text(hjust = 0.5),
        axis.line         = element_line(color='black'))+
  guides(color=guide_legend(override.aes = list(size = 7),
                            title="Target type\n(size)",
                            keyheight=.5,
                            default.unit='inch'),
         linetype = guide_legend(title="Data Type",
                                 keyheight=.4,
                                 default.unit='inch'),
         alpha = FALSE)