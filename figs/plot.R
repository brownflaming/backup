myData <- read.csv("data_7.csv")

## plot for each sample size
df1 <- myData[myData$SampleSize==20, ]
rng <- 1:length(df1$lb)

p <- ggplot(data=df1,aes(x=rng)) +
  geom_line(aes(y = df1$ub), size=1, colour = "blue") +
  geom_point(aes(y = df1$ub), size = 3, color = "blue") +
  geom_errorbar(aes(ymin=df1$ub_l, ymax=df1$ub_r), width = 0.5, color = "blue") + 
  geom_line(aes(y = df1$lb), size=1, colour = "red") +
  geom_point(aes(y = df1$lb), size = 2, color = "red") +
  xlab("iteration") +
  ylab("objective function value") + 
  labs(title="95%-CI for UB and LB improvement (T=7; M=20)") +
  theme(plot.title = element_text(size=20)) + 
  theme(axis.title = element_text(size=15),axis.text = element_text(size=12))
p



# ## plot lower bound for different sample size
# df <- myData[,c("SampleSize","lb")]
# df1 <- df[df$SampleSize==1, ]
# df2 <- df[df$SampleSize==2, ]
# df3 <- df[df$SampleSize==3, ]
# df5 <- df[df$SampleSize==5, ]
# df10 <- df[df$SampleSize==10, ]
# df20 <- df[df$SampleSize==20, ]
# df50 <- df[df$SampleSize==50, ]
# 
# p <- ggplot() +
#   geom_line(data=df1, aes(x=1:length(df1$lb), y=df1$lb, color="M = 1")) +
#   geom_point(data=df1, aes(x=1:length(df1$lb), y=df1$lb, color="M = 1"), size=2) +
#   geom_line(data=df2, aes(x=1:length(df2$lb), y=df2$lb, color="M = 2")) +
#   geom_point(data=df2, aes(x=1:length(df2$lb), y=df2$lb, color="M = 2"), size=2) +
#   geom_line(data=df3, aes(x=1:length(df3$lb), y=df3$lb, color="M = 3")) +
#   geom_point(data=df3, aes(x=1:length(df3$lb), y=df3$lb, color="M = 3"), size=2) +
#   geom_line(data=df5, aes(x=1:length(df5$lb), y=df5$lb, color="M = 5")) +
#   geom_point(data=df5, aes(x=1:length(df5$lb), y=df5$lb, color="M = 5"), size=2) +
#   geom_line(data=df10, aes(x=1:length(df10$lb), y=df10$lb, color="M = 10")) +
#   geom_point(data=df10, aes(x=1:length(df10$lb), y=df10$lb, color="M = 10"), size=2) +
#   geom_line(data=df20, aes(x=1:length(df20$lb), y=df20$lb, color="M = 20")) +
#   geom_point(data=df20, aes(x=1:length(df20$lb), y=df20$lb, color="M = 20"), size=2) +
#   geom_line(data=df50, aes(x=1:length(df50$lb), y=df50$lb, color="M = 50")) +
#   geom_point(data=df50, aes(x=1:length(df50$lb), y=df50$lb, color="M = 50"), size=2)
# 
# p