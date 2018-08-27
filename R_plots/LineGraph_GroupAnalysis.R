library(ggplot2)

ll <- matrix(nrow = 2, ncol = 12)
ll[1,] <-c(0.0533459001803819,-0.133059288365375,-0.242937069823166,-0.276521799379142,-0.226522930535688,-0.176117466816664,-0.143974158765330,0.146891461459652,0.368291785774637,0.333082207554500,	0.192918492609724,0.104602866106456)
ll[2,] <-c(0.0358639526853729,-0.0233894904925072,-0.152973529681280,-0.193053207012406,-0.206138106688492,-0.193819990385164,-0.0393970310389294,0.164119839928179,0.252726190818792,0.265765240362564,0.0887516721118536,0.00154445939190851)

llsem<- matrix(nrow = 2, ncol = 12)
llsem[1,] <- c(0.0465047669082459,0.0461771389883564,0.0430096296690938,0.0522376187230128,0.0522470689339353,0.0701101841473937,0.0712981917981542,0.0665854774247781,0.0952242303925798,0.0747186724990091,0.0728432769911141,0.0632209133511219)
llsem[2,] <- c(0.0607220576402135,0.0679870084373188,0.0527910674927370,0.0552950744385861,0.0871432890479457,0.0708429788891387,0.0614421373974704,0.0556880090102015,0.0451912620875006,0.0569066856220678,0.0611160730049797,0.0510586656173393)

x=seq(0,33,3)
y1 <- ll[1,]-ll[1,1]
y2 <- ll[2,]-ll[2,1]

lldata <-data.frame(X=x, y1=y1, y2=y2, sem1=llsem[1,], sem2=llsem[2,])
p <- ggplot(data=lldata,aes(x, y1))
p+  geom_rect(aes(xmin = 18,xmax = 36,
                  ymin = -Inf, ymax = Inf),
              fill = 'grey70',
              alpha = .07)+
  geom_line(data=lldata,aes(x, y2), colour='red')+
  geom_errorbar(aes(ymin = (lldata$y2-lldata$sem2), ymax=(lldata$y2+lldata$sem2)), position=position_dodge(.9),width=0.6, colour='red')+
  geom_line(data=lldata,aes(x, y1), colour='blue')+
  geom_errorbar(aes(ymin = (lldata$y1-lldata$sem1), ymax=(lldata$y1+lldata$sem1)), position=position_dodge(.9),width=0.6, colour='blue')+
  scale_x_continuous(name="Time (s)", breaks =seq(0,36,6))+
  ylim(-0.4,0.6)+
  ylab("BOLD %")+
  theme_light()+
  theme(axis.title.x = element_text(margin = margin(b = 30),face = "bold", color = "black", size = 14, vjust=-5),axis.title.y = element_text(margin = margin(l = 30),face = "bold", color = "black", size = 14, vjust=25))+
  theme(axis.text.x = element_text(color = "black", size = 12,face='bold'),axis.text.y = element_text(color = "black", size = 12,face='bold'))+
  geom_hline(yintercept=0,linetype="dotted", colour='gray30')

