library(ggplot2)
library(readxl)

# Read in the Excel file
Blinks <- read_excel("~/Desktop/PhD/Negative_BOLD/Plot_data/Blink_average.xls")

# Create a new time variables that shows the time each eye
# recording was taken. Frame rate ws 30Hz.
Blinks$time <- c(0)
Blinks$time <- (Blinks$Number)/30

# Create Group variable to re-code 1's and 2's
Blinks$Grouptxt <- c()
Blinks$Grouptxt[which(Blinks$Group=='1')] <- 'Low'
Blinks$Grouptxt[which(Blinks$Group=='2')] <- 'High'

# Turn into factor variable to use with facet grid plot
Blinks$Grouptxt <- factor(Blinks$Grouptxt, levels = c("Low","High"))

# If the blink value was above 0.027, classify it as a blink
Blinks$Binary <- c(0)
Blinks$Binary[which(Blinks$Corr>=0.027)] = 1

labels <- c(Low = "Low Group", High = "High Group")

p <- ggplot(data=Blinks,aes(time, Binary))
p+ 
  geom_rect(aes(xmin = 18,xmax = 36,
                 ymin = 0, ymax = Inf),
             fill = 'grey70',
             alpha = .07)+
  geom_rect(aes(xmin = 54,xmax = 72,
                ymin = 0, ymax = Inf),
            fill = 'grey70',
            alpha = .07)+
  geom_rect(aes(xmin = 90,xmax = 108,
                ymin = 0, ymax = Inf),
            fill = 'grey70',
            alpha = .07)+
  geom_rect(aes(xmin = 126,xmax = 144,
                ymin = 0, ymax = Inf),
            fill = 'grey70',
            alpha = .07)+
  geom_rect(aes(xmin = 162,xmax = 180,
                ymin = 0, ymax = Inf),
            fill = 'grey70',
            alpha = .07)+
  geom_rect(aes(xmin = 198,xmax = 216,
                ymin = 0, ymax = Inf),
            fill = 'grey70',
            alpha = .07)+
  geom_rect(aes(xmin = 234,xmax = 252,
                ymin = 0, ymax = Inf),
            fill = 'grey70',
            alpha = .07)+
  geom_rect(aes(xmin = 270,xmax = 288,
                ymin = 0, ymax = Inf),
            fill = 'grey70',
            alpha = .07)+
  facet_grid(Grouptxt~., labeller=labeller(Grouptxt = labels))+
  scale_x_continuous(name="Time (s)", breaks = c(0,288))+
  ylim(0,1.1)+
  ylab("Blink Events")+
  theme_classic()+
  theme(axis.title.x = element_text(margin = margin(b = 30),face = "bold", color = "black", size = 14, vjust=-5),axis.title.y = element_text(margin = margin(l = 30),face = "bold", color = "black", size = 14, hjust=0.5,vjust=0.5))+
  theme(axis.text.x = element_text(color = "black", size = 12,face='bold'),axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  geom_line(data=Blinks,aes(time, Binary), colour='black')

