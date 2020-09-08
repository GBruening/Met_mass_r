
# Some test code
# Using logit function
alpha_utility_gross_rx = optimize(fun_optim_prx_logit,
                                  masses     = masses,
                                  mvttimes   = mvttimes,
                                  p_alpha    = met_p_glm,
                                  rxtimes    = rxtimes,
                                  opt_or_dur = 'optimize',
                                  interval   = c(-100,300))

util_dur_gross_rx = as.numeric(fun_optim_prx_logit(
                                       alpha      = alpha_utility_gross_rx$minimum,
                                       masses     = masses,
                                       mvttimes   = mvttimes,
                                       p_alpha    = met_p_glm,
                                       opt_or_dur = 'duration',
                                       rxtimes    = rxtimes))

masses   = masses

masses = c(2.47,3.8,4.7,6.1)
mvttimes = aggregate(movedur~condition,prefdata,mean)$movedur
rxtimes = aggregate(reaction_tanv~condition,prefdata,mean)$reaction_tanv
p_alpha_pref = inv.logit(predict(met_p_glm,data.frame(movedur=mvttimes,eff_mass=masses)))

mvttimes_smallt = aggregate(movedur~condition,smalltdata,mean)$movedur
rxtimes_smallt = aggregate(reaction_tanv~eff_mass,smalltdata,mean)$reaction_tanv
p_alpha_smallt_gross = inv.logit(predict(met_smallt_p_glm,data.frame(movedur=mvttimes_smallt,eff_mass=masses)))

mvttimes = mvttimes_smallt
p_alpha  = smallt_p_glm
rxtimes  = rxtimes_smallt

p_alpha = met_smallt_p_glm

# mvttimes   = mvttimes
# p_alpha    = met_p_glm
# rxtimes    = rxtimes

dur = rep(0,length(mvttimes))
dur2 = rep(0,length(mvttimes))
dur100 = rep(0,length(mvttimes))
k=1

# pref_p_glm smallt_p_glm

m = masses[k]
rt = rxtimes[k]
count = 1
times = seq(0.25,2,0.001)
probs = inv.logit(predict(p_alpha,data.frame(movedur=times,eff_mass=m)))

alpha100  = 100
ut100 = rep(0,length(times))
count = 1
for (t in times){
  ut100[count] = (alpha100*probs[count]-a0*rt-(a*t+b*(m^c)/(t^(d-1))))/(rt+t)
  count = count+1
}
dur100 = times[which.max(ut100)]
p100 =  inv.logit(predict(p_alpha,data.frame(movedur=dur100,eff_mass=m)))

alpha1000    = 1000
ut1000 = rep(0,length(times))
count = 1
for (t in times){
  ut1000[count] = (alpha1000*probs[count]-a0*rt-(a*t+b*(m^c)/(t^(d-1))))/(rt+t)
  count = count+1
}
dur1000 = times[which.max(ut1000)]
p1000 =  inv.logit(predict(p_alpha,data.frame(movedur=dur1000,eff_mass=m)))

alpha5000    = 5000
ut5000 = rep(0,length(times))
count = 1
for (t in times){
  ut5000[count] = (alpha5000*probs[count]-a0*rt-(a*t+b*(m^c)/(t^(d-1))))/(rt+t)
  count = count+1
}
dur5000 = times[which.max(ut5000)]
p5000 =  inv.logit(predict(p_alpha,data.frame(movedur=dur5000,eff_mass=m)))

test_data = data.frame(cbind(c(rep(alpha100,length(times)),
                               rep(alpha1000,length(times)),
                               rep(alpha5000,length(times))),
                             c(times,times,times),
                             c((ut100+100)/max(ut100+100),
                               (ut1000+100)/(max(ut1000+100)),
                               (ut5000+100)/(max(ut5000+100)))))
  
colnames(test_data) = c('alpha','movedur','utility')
ggplot()+
  geom_line(data=filter(test_data),
            aes(x=movedur,
                y=utility,
                color=factor(alpha)),
            size = 2)+
  geom_segment(aes(x = times[which.max(ut100)],
                   y = 0,
                   xend = times[which.max(ut100)],
                   yend = 1))+
  annotate(geom = 'text',
           x = times[which.max(ut100)],
           y = 1.05,
           label = round(p100,4))+
  geom_segment(aes(x = times[which.max(ut1000)],
                   y = 0,
                   xend = times[which.max(ut1000)],
                   yend = 1))+
  annotate(geom = 'text',
           x = times[which.max(ut1000)],
           y = 1.1,
           label = round(p1000,4))+
  geom_segment(aes(x = times[which.max(ut5000)],
                   y = 0,
                   xend = times[which.max(ut5000)],
                   yend = 1))+
  annotate(geom = 'text',
           x = times[which.max(ut5000)],
           y = 1.15,
           label = round(p5000,4))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.20))+
  scale_x_continuous(expand = c(0, 0), limits = c(.3, 1.20))+
  labs(x='Movement Duration (s)',y='Utility',color='alpha',title='Experiment 2b')


inv.logit(predict(p_alpha,data.frame(movedur=util_dur_smallt_gross,eff_mass=masses)))

# More Test code
mvttimes = mvttimes_smallt
p_alpha  = met_smallt_p_glm
rxtimes  = rxtimes_smallt

alphas = exp(seq(1,20,.2))

dur = rep(0,length(alphas))
p = rep(0,length(alphas))
rew = rep(0,length(alphas))
eff = rep(0,length(alphas))

m = masses[1]
rt = rxtimes[1]
count = 1
times = seq(0.4,2,0.001)
probs = inv.logit(predict(p_alpha,data.frame(movedur=times,eff_mass=m)))

alph_c = 1
ut = matrix(,nrow = length(times), ncol = length(alphas))
for (alpha100 in alphas){
  count = 1
  for (t in times){
    ut[count,alph_c] = (alpha100*probs[count]-a0*rt-(a*t+b*(m^c)/(t^(d-1))))/(rt+t)
    count = count+1
  }
  dur[alph_c] = times[which.max(ut[,alph_c])]
  p[alph_c] =  inv.logit(predict(p_alpha,data.frame(movedur=dur[alph_c],eff_mass=m)))
  rew[alph_c] = (alpha100*p[alph_c])#/(rt+dur[alph_c])
  eff[alph_c] = a0*rt-(a*times[which.max(ut[,alph_c])]+b*(m^c)/(times[which.max(ut[,alph_c])]^(d-1)))
  alph_c = alph_c + 1
}
ut[ut<0] = 0

ut_col = c()
eff_col = c()
alph_col = c()
times_col = c()
dur_col = c()
for (k in c(1:length(alphas))){
  times_col = c(times_col,times)
  dur_col = c(dur_col,rep(dur[k],length(times)))
  ut_col = c(ut_col,ut[,k])
  alph_col = c(alph_col,rep(alphas[k],length(times)))
  eff_col = c(eff_col,rep(alphas[k],length(times)))
}
dataf = data.frame(times = times_col,
                  opt_dur = dur_col,
                  util = ut_col,
                  alph = alph_col)
dataf = filter(dataf,alph<5000)
ggplot(data = filter(aggregate(opt_dur ~ alph, dataf, mean),alph<5000))+geom_point(aes(x = opt_dur, y = alph))
ggplot(data = filter(dataf,alph<1000,alph>100))+geom_line(aes(x = times, y = util,color = factor(alph)))

