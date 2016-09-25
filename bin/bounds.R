library(ggplot2)
mydata <- read.csv("data.csv")

ggplot(mydata, aes(iteration)) + 
  geom_line(aes(y = LB, colour = "LB"), size=1) + 
  geom_line(aes(y = UB, colour = "UB"), size=1) +
  geom_ribbon(aes(ymin=UBl, ymax=UBr), alpha=0.2) +
  xlab("iteration") +
  ylab("bounds") +
  ylim(2800, 8500) +
  labs(title="cut: L") +
  theme(legend.title=element_blank()) +
  theme(legend.position="none")

mydata <- read.csv("data2.csv")

ggplot(mydata, aes(iteration)) + 
  #geom_line(aes(y = B, colour = "B"), size=1) + 
  #geom_line(aes(y = SB, colour = "SB"), size=1) + 
  #geom_line(aes(y = L, colour = "L"), size=1) + 
  geom_line(aes(y = I, colour = "I"), size=1) + 
  geom_line(aes(y = BI, colour = "I + B"), size=1) +
  #geom_line(aes(y = BL, colour = "B + L"), size=1) +
  geom_line(aes(y = SBI, colour = "I + SB"), size=1) +
  #geom_line(aes(y = SBL, colour = "SB + L"), size=1) + 
  geom_line(aes(y = LI, colour = "I + L"), size=1) + 
  xlab("iteration") +
  ylab("lower bounds") +
  labs(title="lower bound improvement") +
  theme(legend.title=element_blank())
  
#theme(legend.position="none")

#guides(colour = guide_legend(override.aes = list(size=2))) + 
#theme(legend.text = element_text(size=15))
#theme(plot.title = element_text(size=20)) + 
#theme(axis.title = element_text(size=15),axis.text = element_text(size=12)) +