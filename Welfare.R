
mydata = data.frame(months=c("-0.4 to -0.2", "-0.2 to 0", "0 to 0.5", "0.5 to 1"), values=c(0.7, 0.2, 0.0, 0.05))

mydata$months <-factor(mydata$months, 
                       levels = c("-0.4 to -0.2", "-0.2 to 0", "0 to 0.5", "0.5 to 1"))

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/welfare_dist.pdf",width=4,height=3)
ggplot(mydata, aes(months, values)) +
  geom_bar(stat = "identity") + geom_abline(slope=0, intercept=0.25,  col = "red",lty=2, size=1) + 
  ylab("% Consumption Equivalent") + xlab("% Equity in House")+
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom")) + 
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid.major = element_line(colour = "lightgray"),
        plot.title = element_text(hjust = 0.5)) + annotate("text", x = "0 to 0.5", y = 0.3, label = "Average", size=4.5, colour="red")
dev.off()





# Welfare equivalent gains by equity level
# Este ya esta mejor
# 
# mydata = data.frame(months=c("(0.4)-(0.2)", "(0.2)-0", "0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1", "Non-owners"), 
#                     values=c(2.5, 0.9, 0.05, 0.01, 0.0, -0.11, -0.11, 0.11))
# 
# mydata2 = data.frame(months=c("(0.4)-(0.2)", "(0.2)-0", "0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1", "Non-owners"), 
#                     values=c(5.1-2.5, 2.4-0.9, 0.06-0.05, 0.055-0.01, 0.05, 0.05+0.11, 0.02+0.11, -0.025-0.11))
# 
# df = data.frame(months=c("(0.4)-(0.2)", "(0.2)-0", "0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1", "Non-owners","(0.4)-(0.2)", "(0.2)-0", "0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1", "Non-owners"),
#                 num = c(rep("Better eligibility",8), rep("Subsidy-only",8)), values=c(mydata$values,mydata2$values))
# 
# df$months <-factor(mydata$months, 
#                    levels = c("(0.4)-(0.2)", "(0.2)-0", "0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1", "Non-owners"))
# 
# pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/welfare_gains_dist.pdf",width=6.5,height=3)
# # ggplot(df, aes(months, values, values2)) +
# #   geom_bar(stat = "dodge") + geom_abline(slope=0, intercept=0.21,  col = "red",lty=2, size=1) + 
# #   ylab("% Consumption Equivalent") + xlab("% Equity in House")+
# #   theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom")) + 
# #   theme(panel.background = element_rect(fill='white', colour='black'),
# #         panel.grid.major = element_line(colour = "lightgray"),
# #         plot.title = element_text(hjust = 0.5)) + annotate("text", x = "0.4-0.6", y = 0.4, label = "Average: 0.21%", size=4.5, colour="red")
# ggplot(df, aes(x=factor(months), y=values, fill=factor(num))) +
#   geom_bar(stat="identity") + geom_abline(slope=0, intercept=0.21,  col = "red",lty=2, size=1) + 
#   ylab("% Consumption Equivalent") + xlab("% Equity in House") + 
#   theme(legend.title = element_blank()) + 
#   scale_fill_manual("factor", values = c("Better eligibility" = "gray61", "Subsidy-only" = "gray19"))+
#   theme(panel.border = element_rect(colour = "black", fill=NA), legend.key = element_blank()) + 
#   theme(panel.background = element_rect(fill='white', colour='black'),
#         panel.grid.major = element_line(colour = "lightgray"),
#         plot.title = element_text(hjust = 0.5)) + annotate("text", x = "0.4-0.6", y = 0.65, label = "Average: 0.21%", size=4.5, colour="red")
# dev.off()


mydata = data.frame(months=c("-40 to -20", "-20 to 0", "0 to 20", "20 to 40", "40 to 60", "60 to 80", "80 to 100", "Non-owners"), 
                    values=c(2.5, 0.9, 0.05, 0.01, 0.0, -0.11, -0.11, 0.11),
                    values2=c(5.1, 2.4, 0.06, 0.055, 0.05, 0.05, 0.02, -0.025))

mydata$months <-factor(mydata$months, 
                       levels = c("-40 to -20", "-20 to 0", "0 to 20", "20 to 40", "40 to 60", "60 to 80", "80 to 100", "Non-owners"))

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/welfare_gains_dist.pdf",width=6,height=3)
ggplot(mydata, aes(months, values, values2)) +
  geom_bar(stat = "identity") + geom_abline(slope=0, intercept=0.21,  col = "red",lty=2, size=1) + 
  ylab("Consumption Equivalent (%)") + xlab("% Equity in House")+
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom")) + 
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid.major = element_line(colour = "lightgray"),
        plot.title = element_text(hjust = 0.5)) + annotate("text", x = "40 to 60", y = 0.4, label = "Average: 0.21%", size=4.5, colour="red")
dev.off()


# Welfare equivalent gains by wealth

mydata = data.frame(months=c("-0.4 to -0.2", "-0.2 to 0", "0 to 0.2", "0.2 to 0.4", "0.4 to 0.6", "0.6 to 0.8", "0.8 to 1", "Non-owners"), values=c(2.3, 0.9, 0.2, 0.15, 0.15, 0.15, 0.05, 0.05))

mydata$months <-factor(mydata$months, 
                       levels = c("-0.4 to -0.2", "-0.2 to 0", "0 to 0.2", "0.2 to 0.4", "0.4 to 0.6", "0.6 to 0.8", "0.8 to 1", "Non-owners"))

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/welfare_dist.pdf",width=4,height=3)
ggplot(mydata, aes(months, values)) +
  geom_bar(stat = "identity") + geom_abline(slope=0, intercept=0.23,  col = "red",lty=2, size=1) + 
  ylab("% Consumption Equivalent") + xlab("% Equity in House")+
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom")) + 
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid.major = element_line(colour = "lightgray"),
        plot.title = element_text(hjust = 0.5)) + annotate("text", x = "0.4 to 0.6", y = 0.3, label = "Average: 0.23%", size=4.5, colour="red")
dev.off()






# Welfare equivalent gains by age


mydata = data.frame(months=c("20-30", "30-40", "40-50", "50-60", "60-70", "70-80"), 
                    values=c(0.65-0.45, 1.1-0.45, 1.2-0.5, 0.9-0.35, 1-0.2, 0.01))

mydata2 = data.frame(months=c("20-30", "30-40", "40-50", "50-60", "60-70", "70-80"), 
                     values=c(0.45, 0.45, 0.5, 0.35, -0.2, 0.01))

df = data.frame(months=c("20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80"),
                num = c(rep("Better eligibility",6), rep("Subsidy-only",6)), values=c(mydata$values,mydata2$values))

df$months <-factor(mydata$months, 
                   levels = c("20-30", "30-40", "40-50", "50-60", "60-70", "70-80"))


ggplot(df, aes(x=factor(months), y=values, fill=factor(num))) +
  geom_bar(stat="identity") + geom_abline(slope=0, intercept=0.21,  col = "red",lty=2, size=1) + 
  ylab("% Consumption Equivalent") + xlab("% Equity in House") + 
  theme(legend.title = element_blank()) + 
  scale_fill_manual("factor", values = c("Better eligibility" = "gray61", "Subsidy-only" = "gray19"))+
  theme(panel.border = element_rect(colour = "black", fill=NA), legend.key = element_blank()) + 
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid.major = element_line(colour = "lightgray"),
        plot.title = element_text(hjust = 0.5)) + annotate("text", x = "60-70", y = 0.28, label = "Average: ", size=4.5, colour="red") +
  annotate("text", x = "70-80", y = 0.28, label = "0.21%", size=4.5, colour="red")

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/age_gains_dist.pdf",width=5,height=3)
ggplot(mydata, aes(months, values)) +
  geom_bar(stat = "identity") + geom_abline(slope=0, intercept=0.23,  col = "red",lty=2, size=1) + 
  ylab("% Consumption Equivalent") + xlab("Age")+
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom")) + 
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid.major = element_line(colour = "lightgray"),
        plot.title = element_text(hjust = 0.5)) + annotate("text", x = "60-70", y = 0.28, label = "Average: ", size=4.5, colour="red") +
        annotate("text", x = "70-80", y = 0.28, label = "0.21%", size=4.5, colour="red")
dev.off()

