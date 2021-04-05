tempdf1 = aggregate(peakvel_radv ~ condition + targetnum,prefdata,mean)
tempdf1$exp = '2a'
tempdf2 = aggregate(peakvel_radv ~ condition + targetnum,smalltdata,mean)
tempdf2$exp = '2b'
tempdf3 = aggregate(peakvel_radv ~ condition + targetnum,pilotdata,mean)
tempdf3$exp = '2c'

tempdf = rbind(tempdf1,tempdf2,tempdf3)
tempdf4 = aggregate(peakvel_radv ~ condition + targetnum,tempdf,mean)
tempdf4$exp = 'Average'
tempdf = rbind(tempdf,tempdf4)
g <- ggplot()+
  geom_line(data=tempdf,
            aes(x = factor(targetnum),
                y = peakvel_radv, 
                group = factor(condition)))+
  geom_point(data=tempdf,
             aes(x = factor(targetnum),
                 y = peakvel_radv,
                 color = factor(condition),
                 size = 3))+
  labs(x = 'Target Number (1 = Top Right, 4 = Bottom right)',
       y = 'Radial Peak Velocity (m/s)',
       color = 'Added mass (lbs)')+
  facet_wrap(~exp)
setwd('D:/Google Drive/Preferred Mass/alaa code/Graphs_pref')
ggsave('Peakv_by_target.pdf',g)
ggsave('Peakv_by_target.png',g)



tempdf1 = aggregate(reaction_tanv ~ condition + targetnum,prefdata,mean)
tempdf1$exp = '2a'
tempdf2 = aggregate(reaction_tanv ~ condition + targetnum,smalltdata,mean)
tempdf2$exp = '2b'
tempdf3 = aggregate(reaction_tanv ~ condition + targetnum,pilotdata,mean)
tempdf3$exp = '2c'

tempdf = rbind(tempdf1,tempdf2,tempdf3)
tempdf4 = aggregate(reaction_tanv ~ condition + targetnum,tempdf,mean)
tempdf4$exp = 'Average'
tempdf = rbind(tempdf,tempdf4)
g <- ggplot()+
  geom_line(data=tempdf,
            aes(x = factor(targetnum),
                y = reaction_tanv, 
                group = factor(condition)))+
  geom_point(data=tempdf,
             aes(x = factor(targetnum),
                 y = reaction_tanv,
                 color = factor(condition),
                 size = 3))+
  labs(x = 'Target Number (1 = Top Right, 4 = Bottom right)',
       y = 'Reaction time (s)',
       color = 'Added mass (lbs)')+
  facet_wrap(~exp)
setwd('D:/Google Drive/Preferred Mass/alaa code/Graphs_pref')
ggsave('Reacttime_by_target.pdf',g)
ggsave('Reacttime_by_target.png',g)



