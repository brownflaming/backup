myData <- read.csv("comparison.csv")

head(myData)
rng <- 1:59
ggplot(data = myData, aes(x=rng)) +
  geom_line(aes(y=myData$B, color="B"),size=1.0) +
  geom_line(aes(y=myData$SB, color="SB"),size=1.0) +
  geom_line(aes(y=myData$L, color="L"),size=1.0) +
  geom_line(aes(y=myData$I, color="I"),size=1.0) +
  geom_line(aes(y=myData$BL, color="B + L"),size=1) +
  geom_line(aes(y=myData$BI, color="B + I"),size=1) +
  geom_line(aes(y=myData$SBL, color="SB + L"),size=1) +
  geom_line(aes(y=myData$SBI, color="SB + I"),size=1) +
  geom_line(aes(y=myData$SBLI, color="SB + L + I"),size=1) +
  geom_line(aes(y=myData$opt, color="Opt. Val."),size=1) +
  xlab("iteration") +
  ylab("objective function value") + 
  labs(title="Lower bound improvement comparison") +
  theme(plot.title = element_text(size=20)) + 
  theme(axis.title = element_text(size=15),axis.text = element_text(size=12)) +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  theme(legend.text = element_text(size=15))