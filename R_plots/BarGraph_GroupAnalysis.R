library(haven)
library(ggplot2)
library(ggsignif)

v1NBR <- read_sav("~/Desktop/PhD/Negative_BOLD/Plot_data/V1_PBR_noncon.sav")

v1NBR$Group <- factor(v1NBR$Group, levels = c("Low","High"))
v1NBR$Hemifield <- factor(v1NBR$Hemifield, levels = c("Left","Right"))



b <- ggplot(v1NBR, aes((Group),Mean, fill=factor((Hemifield))))
b + geom_bar(aes(fill=Hemifield),stat='identity', position='dodge', color='black')+
  geom_errorbar(aes(ymin = (v1NBR$Mean-v1NBR$SE), ymax=(v1NBR$Mean+v1NBR$SE)), position=position_dodge(.9),width=0.3)+
  xlab("Group")+
  ylab("BOLD %")+
  theme_light()+
  theme(legend.position='right')+
  geom_hline(yintercept=0,linetype="dotted", colour='gray30')+
  scale_fill_manual("Stimulus Hemifield", values = c("gray45","gray70"))+
  theme(axis.title.x = element_text(margin = margin(b = 30),face = "bold", color = "black", size = 14, vjust=-5),axis.title.y = element_text(margin = margin(l = 30),face = "bold", color = "black", size = 14, vjust=20))+
  theme(axis.text.x = element_text(color = "black", size = 10,face='bold'),axis.text.y = element_text(color = "black", size = 10,face='bold'))+
  theme(legend.title = element_text(color = "black", size = 10), legend.key.size =unit(1.2, 'line'), legend.text = element_text(color = "black", size = 8))+
  theme(strip.text.x = element_text(face='bold',color = "black", size = 12))



