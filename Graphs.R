
# -----------------------------------------------------------#
#                  Prelimiaries                              #
# -----------------------------------------------------------#

setwd("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Model_externalities_19/")

source("Rcpp_Run.R")

# ------------------------------------------------#
#          Graphs for document and beamer         #
# ------------------------------------------------#

# Lifecycles
pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/lifecycle_figs.pdf",width=9.5,height=5.5)
grid.arrange(arrangeGrob(pC + theme(legend.position="none"), pR, pD, ncol = 3, nrow = 1),
             arrangeGrob(pH, pA, pM, ncol = 3, nrow = 1),
             ncol = 1, nrow = 2, heights=c(1,1))
dev.off()


# Lifecycles
pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/lifecycle_figs_all.pdf",width=9.5,height=9.5)
grid.arrange(arrangeGrob(pC + theme(legend.position="none"), pR, pdista, ncol = 3, nrow = 1),
             arrangeGrob(pH, pM, pMort, ncol = 3, nrow = 1),
             arrangeGrob(pA, pD, pdisth, ncol = 3, nrow = 1),
             arrangeGrob(pRR, pequity, pPP, ncol = 3, nrow = 1), 
             ncol = 1, nrow = 4, heights=c(1,1,1,1))
dev.off()


peq_cycle_doc =   ggplot(data = data.frame(dist=eq_cycle, grid=seq(24,80,length.out =11)), aes(x=grid)) + 
  geom_line(aes(y = dist, linetype = "SCF 2007", colour = "SCF 2007"), size=size_line_doc) + ggtitle("% Home Equity by Age") + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot + ylab("%") + xlab("Age") + 
  geom_line(data = data.frame(dist=equity_cycle, grid=seq(24,80,length.out =T)), aes(x=grid, y=dist, linetype="Model", colour="Model"), size=size_line_doc, linetype = "solid") +
  scale_colour_manual("", breaks = c("Model", "SCF 2007"), values = c("Model" = "red", "SCF 2007" = "black")) +
  scale_linetype_manual("", breaks = c("Model", "SCF 2007"), values = c("Model" = "solid", "SCF 2007" = "dashed"))
  

# Lifecycles some
pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/lifecycle_figs_selected.pdf",width=7.5,height=10)
grid.arrange(arrangeGrob(pC + theme(legend.position="none"), pR, ncol = 2, nrow = 1),
             arrangeGrob(peq_cycle, pdebt_cycle, ncol = 2, nrow = 1),
             arrangeGrob(pequity, pPP, ncol = 2, nrow = 1),
             arrangeGrob(mylegend, ncol = 1, nrow = 1),
             ncol = 1, nrow = 4, heights=c(1,1,1,0.2))
dev.off()

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/lifecycle_figs_beamer.pdf",width=8.5,height=5.2)
grid.arrange(arrangeGrob(pR, peq_cycle, ncol = 2, nrow = 1),
             arrangeGrob(pequity, mylegend, ncol = 2, nrow = 1),
             ncol = 1, nrow = 2, heights=c(1,1))
dev.off()


pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/lifecycle_figs_beamer_app.pdf",width=8.5,height=5.2)
grid.arrange(arrangeGrob(pC + theme(legend.position="none"), pD, ncol = 2, nrow = 1),
             arrangeGrob(pRR, pPP, ncol = 2, nrow = 1),
             ncol = 1, nrow = 2, heights=c(1,1))
dev.off()


# pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/moments.pdf",width=8,height=8)
# grid.arrange(tbl2, ncol = 1, nrow = 1)
# dev.off()


# Equity distributions
pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/equity_CDF_tr0.pdf",width=6.5,height=3.5)
ggplot(data = data.frame(dist=dist_equity, grid=seq(0,100,length.out =100)), aes(x=grid)) + 
  geom_line(data = data.frame(dist=home_equity, grid=100*seq(-0.3,1,length.out =130)), aes(x=grid, y=100*dist), size=size_line, linetype = "solid") + ggtitle("Home Equity CDF") + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot + ylab("%") + xlab("Equity (%)")  +
  annotate("text", label = "Pre-2007", x = 0.25, y = 75, size = 4.5, colour = "black")
dev.off()

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/equity_CDF_tr1.pdf",width=6.5,height=3.5)
ggplot(data = data.frame(dist=dist_equity, grid=seq(0,1,length.out =100)), aes(x=grid)) + 
  geom_line(data = data.frame(dist=home_equity, grid=100*seq(-0.3,1,length.out =130)), aes(x=grid, y=100*dist), size=size_line, linetype = "solid") + ggtitle("Home Equity CDF") + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot + ylab("%") + xlab("Equity (%)") + 
  geom_line(data = data.frame(dist=home_equity1, grid=100*seq(-0.3,1,length.out =130)), aes(x=grid, y=100*dist), size=size_line, linetype = "solid", colour="red") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red", size=size_line*0.7)+
  annotate("text", label = "Pre-2007", x = 25, y = 75, size = 4.5, colour = "black")+
  annotate("text", label = "2007-2009", x = 25, y = 65, size = 4.5, colour = "red")+
  annotate("text", label = "Equity < 0:", x = -13, y = 40, size = 4.5, colour = "red")+
  annotate("text", label = "Model: 16.2%", x = -13, y = 30, size = 4.2, colour = "red")+
  annotate("text", label = "Data: 15.0%", x = -13, y = 20, size = 4.2, colour = "red")
dev.off()

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/equity_CDF_tr2.pdf",width=6.5,height=3.5)
ggplot(data = data.frame(dist=dist_equity, grid=seq(0,1,length.out =100)), aes(x=grid)) + 
  geom_line(data = data.frame(dist=home_equity, grid=100*seq(-0.3,1,length.out =130)), aes(x=grid, y=100*dist), size=size_line, linetype = "solid") + ggtitle("Home Equity CDF") + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot + ylab("%") + xlab("Equity") + 
  geom_line(data = data.frame(dist=home_equity1, grid=100*seq(-0.3,1,length.out =130)), aes(x=grid, y=100*dist), size=size_line, linetype = "solid", colour="red") + 
  geom_line(data = data.frame(dist=home_equity2, grid=100*seq(-0.3,1,length.out =130)), aes(x=grid, y=100*dist), size=size_line, linetype = "solid", colour="blue") +
  annotate("text", label = "Pre-2007", x = 25, y = 75, size = 4.5, colour = "black")+
  annotate("text", label = "2007-2009", x = 25, y = 65, size = 4.5, colour = "red")+
  annotate("text", label = "2009-2011", x = 25, y = 55, size = 4.5, colour = "blue")
dev.off()

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/equity_CDF_tr3.pdf",width=6.5,height=3.5)
ggplot(data = data.frame(dist=dist_equity, grid=seq(0,1,length.out =100)), aes(x=grid)) + 
  geom_line(data = data.frame(dist=home_equity, grid=100*seq(-0.3,1,length.out =130)), aes(x=grid, y=100*dist), size=size_line, linetype = "solid") + ggtitle("Home Equity CDF") + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot + ylab("%") + xlab("Equity") + 
  geom_line(data = data.frame(dist=home_equity1, grid=100*seq(-0.3,1,length.out =130)), aes(x=grid, y=100*dist), size=size_line, linetype = "solid", colour="red") + 
  geom_line(data = data.frame(dist=home_equity2, grid=100*seq(-0.3,1,length.out =130)), aes(x=grid, y=100*dist), size=size_line, linetype = "solid", colour="blue") +
  geom_line(data = data.frame(dist=home_equity4, grid=100*seq(-0.3,1,length.out =130)), aes(x=grid, y=100*dist), size=size_line, linetype = "solid", colour="purple") +
  annotate("text", label = "Pre-2007", x = 25, y = 75, size = 4.5, colour = "black")+
  annotate("text", label = "2007-2009", x = 25, y = 65, size = 4.5, colour = "red")+
  annotate("text", label = "2009-2011", x = 25, y = 55, size = 4.5, colour = "blue")+
  annotate("text", label = "2013-2015", x = 25, y = 45, size = 4.5, colour = "purple")
dev.off()

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/equity_CDF_tr4.pdf",width=6.5,height=3.5)
ggplot(data = data.frame(dist=dist_equity, grid=seq(0,1,length.out =100)), aes(x=grid)) + 
  geom_line(data = data.frame(dist=home_equity, grid=100*seq(-0.3,1,length.out =130)), aes(x=grid, y=100*dist), size=size_line, linetype = "solid") + ggtitle("Home Equity CDF") + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot + ylab("%") + xlab("Equity") + 
  geom_line(data = data.frame(dist=home_equity1, grid=100*seq(-0.3,1,length.out =130)), aes(x=grid, y=100*dist), size=size_line, linetype = "solid", colour="red") + 
  geom_line(data = data.frame(dist=home_equity2, grid=100*seq(-0.3,1,length.out =130)), aes(x=grid, y=100*dist), size=size_line, linetype = "solid", colour="blue") +
  geom_line(data = data.frame(dist=home_equity4, grid=100*seq(-0.3,1,length.out =130)), aes(x=grid, y=100*dist), size=size_line, linetype = "solid", colour="purple") +
  geom_line(data = data.frame(dist=home_equity6, grid=100*seq(-0.3,1,length.out =130)), aes(x=grid, y=100*dist), size=size_line, linetype = "solid", colour="orange") +
  annotate("text", label = "Pre-2007", x = 25, y = 75, size = 4.5, colour = "black")+
  annotate("text", label = "2007-2009", x = 25, y = 65, size = 4.5, colour = "red")+
  annotate("text", label = "2009-2011", x = 25, y = 55, size = 4.5, colour = "blue")+
  annotate("text", label = "2013-2015", x = 25, y = 45, size = 4.5, colour = "purple") +
  annotate("text", label = "2017-2019", x = 25, y = 35, size = 4.5, colour = "orange")
dev.off()

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/equity_CDF_doc_trans.pdf",width=7,height=3)
ggplot(data = data.frame(dist=dist_equity, grid=seq(0,1,length.out =100)), aes(x=grid)) + 
  geom_line(data = data.frame(dist=home_equity, grid=100*seq(-0.3,1,length.out =130)), aes(x=grid, y=100*dist, linetype = "Pre-crisis", colour="Pre-crisis"), size=size_line_doc) +
  scale_x_continuous(expand = c(0, 0)) + theme_plot + ylab("%") + xlab("Equity") + 
  geom_line(data = data.frame(dist=home_equity1, grid=100*seq(-0.3,1,length.out =130)), aes(x=grid, y=100*dist, linetype = "2007", colour="2007"), size=size_line_doc) + 
  geom_line(data = data.frame(dist=home_equity4, grid=100*seq(-0.3,1,length.out =130)), aes(x=grid, y=100*dist, linetype = "2013", colour="2013"), size=size_line_doc) +
  geom_line(data = data.frame(dist=home_equity6, grid=100*seq(-0.3,1,length.out =130)), aes(x=grid, y=100*dist, linetype = "2017", colour="2017"), size=size_line_doc) +
  scale_colour_manual("", breaks = c("Pre-crisis", "2007", "2013", "2017"), values = c("Pre-crisis" = "black", "2007" = "red", "2013" = "purple", "2017" = "orange")) +
  scale_linetype_manual("", breaks = c("Pre-crisis", "2007", "2013", "2017"), values = c("Pre-crisis" = "solid", "2007" = "dashed", "2013" = "dotdash", "2017" = "twodash")) + 
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))
dev.off()

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/equity_CDF_tr5.pdf",width=6.5,height=3.5)
ggplot(data = data.frame(dist=dist_equity, grid=seq(0,1,length.out =100)), aes(x=grid)) + 
  geom_line(data = data.frame(dist=home_equity, grid=seq(-0.3,1,length.out =130)), aes(x=grid, y=dist), size=size_line, linetype = "solid") + ggtitle("Home Equity CDF") + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot + ylab("%") + xlab("Equity") + 
  geom_line(data = data.frame(dist=home_equity1, grid=seq(-0.3,1,length.out =130)), aes(x=grid, y=dist), size=size_line, linetype = "solid", colour="red") + 
  geom_line(data = data.frame(dist=home_equity2, grid=seq(-0.3,1,length.out =130)), aes(x=grid, y=dist), size=size_line, linetype = "solid", colour="blue") +
  geom_line(data = data.frame(dist=home_equity4, grid=seq(-0.3,1,length.out =130)), aes(x=grid, y=dist), size=size_line, linetype = "solid", colour="purple") +
  annotate("text", label = "Pre-2007", x = 0.25, y = 0.75, size = 4.5, colour = "black")+
  annotate("text", label = "2007-2009", x = 0.25, y = 0.65, size = 4.5, colour = "red")+
  annotate("text", label = "2009-2011", x = 0.25, y = 0.55, size = 4.5, colour = "blue")+
  annotate("text", label = "2013-2015", x = 0.25, y = 0.45, size = 4.5, colour = "purple")
dev.off()

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/equity_CDF.pdf",width=6.5,height=3.5)
ggplot(data = data.frame(dist=dist_equity, grid=seq(0,1,length.out =100)), aes(x=grid)) + 
  geom_line(aes(y = dist, linetype = "SCF2007", color = "SCF2007"), size=size_line_doc) + 
  geom_line(data = data.frame(dist=home_equity, grid=seq(-0.3,1,length.out =130)), aes(x=grid, y=dist, linetype = "Model", color = "Model"), size=size_line_doc) +
  scale_x_continuous(expand = c(0, 0)) + theme_plot + ylab("%") + xlab("Equity") + 
  scale_linetype_manual("", breaks = c("Model", "SCF2007"), values = c("Model"="solid", "SCF2007"="dashed"))+
  scale_color_manual("", breaks = c("Model", "SCF2007"), values = c("Model"="red", "SCF2007"="black")) +
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))
dev.off()

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/equity_CDF_doc.pdf",width=6,height=3)
ggplot(data = data.frame(dist=dist_equity, grid=seq(0,100,length.out =100)), aes(x=grid)) + 
  geom_line(aes(y = 100*dist, linetype = "SCF2007", color = "SCF2007"), size=size_line_doc) + 
  geom_line(data = data.frame(dist=home_equity, grid=seq(-30,100,length.out =130)), aes(x=grid, y=100*dist, linetype = "Model", color = "Model"), size=size_line_doc) +
  scale_x_continuous(expand = c(0, 0)) + theme_plot + ylab("%") + xlab("Equity (%)") + 
  scale_linetype_manual("", breaks = c("Model", "SCF2007"), values = c("Model"="solid", "SCF2007"="dashed"))+
  scale_color_manual("", breaks = c("Model", "SCF2007"), values = c("Model"="red", "SCF2007"="black")) +
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))
dev.off()

g_eq_cdf = ggplot(data = data.frame(dist=dist_equity, grid=seq(0,100,length.out =100)), aes(x=grid)) + 
  geom_line(aes(y = 100*dist, linetype = "SCF2007", color = "SCF2007"), size=size_line_doc) + 
  geom_line(data = data.frame(dist=home_equity, grid=100*seq(-0.3,1,length.out =130)), aes(x=grid, y=100*dist, linetype = "Model", color = "Model"), size=size_line_doc) +
  scale_x_continuous(expand = c(0, 0)) + theme_plot + ylab("%") + xlab("Equity (%)") + 
  ggtitle("Cummulative Distribution") +
  scale_linetype_manual("", breaks = c("Model", "SCF2007"), values = c("Model"="solid", "SCF2007"="dashed"))+
  scale_color_manual("", breaks = c("Model", "SCF2007"), values = c("Model"="red", "SCF2007"="black")) +
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))

mylegend<-g_legend(g_eq_cdf) 

g_eq_cycle = ggplot(data = data.frame(dist=eq_cycle, grid=seq(24,80,length.out =11)), aes(x=grid)) + 
  geom_line(aes(y = 100*dist, linetype = "Data", color = "Data"), size=size_line_doc) + 
  geom_line(data = data.frame(dist=equity_cycle, grid=seq(24,80,length.out =T)), aes(x=grid, y=100*dist, linetype = "Model", color = "Model"), size=size_line_doc) +
  ggtitle("Life-Cycle") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot + ylab("%") + xlab("Age") + 
  scale_y_continuous(expand = c(0, 0), limit = c(0, 100)) +
  scale_linetype_manual("", breaks = c("Model", "Data"), values = c("Model"="solid", "Data"="dashed"))+
  scale_color_manual("", breaks = c("Model", "Data"), values = c("Model"="red", "Data"="black")) +
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/equity_both.pdf",width=6.5,height=3)
grid.arrange(arrangeGrob(g_eq_cdf + theme(legend.position="none"),
                         g_eq_cycle + theme(legend.position="none"),
                         nrow = 1, ncol = 2, widths=c(1,1)),
             mylegend, nrow = 2, ncol = 1, heights=c(1,0.2))
dev.off()


# Default in steady state
g_def_cycle = ggplot(data = data.frame(D=D, age=age_vector), aes(x=age)) + 
  geom_line(aes(y = 100*D), size=size_line_doc) + 
  ggtitle("Default (%)") + 
  scale_y_continuous(limit = c(0, 30), expand = c(0, 0))  +  xlab("Age")  +  ylab("%") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot

g_def_cycle_tr =  ggplot(data = data.frame(D=D, D1=D1, age=age_vector), aes(x=age)) + 
  geom_line(aes(y = 100*D), size=size_line_doc) + 
  geom_line(aes(y = 100*D1), size=size_line_doc, colour="red") + 
  ggtitle("Default (%) vs. Age") +
  annotate("text", label = "Pre-2007", x = 60, y = 30, size = 4, colour = "black") +
  annotate("text", label = "2007-2009", x = 60, y = 25, size = 4, colour = "red") +
  scale_y_continuous(limit = c(0, 40), expand = c(0, 0))  +  xlab("Age")+  ylab("%") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot


g_eq_tr = ggplot(data = data.frame(Eq=Eq, Eq1=Eq1, age=age_vector), aes(x=age)) + 
  geom_line(aes(y = 100*Eq), size=size_line_doc) + 
  geom_line(aes(y = 100*Eq1), size=size_line_doc, colour="red") + 
  ggtitle("Equity (%) vs. Age") +
  annotate("text", label = "Pre-2007", x = 60, y = 30, size = 4, colour = "black") +
  annotate("text", label = "2007-2009", x = 60, y = 25, size = 4, colour = "red") +
  scale_y_continuous(limit = c(-5, 65), expand = c(0, 0))  +  xlab("Age")+  ylab("%") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot


g_def_sav_tr = ggplot(data = data.frame(Da=Da, Da1=Da1, age=seq(agrid[1],agrid[na],length.out = na)), aes(x=age)) + 
  geom_line(aes(y = Da), size=size_line_doc) + 
  geom_line(aes(y = Da1), size=size_line_doc, colour="red") + 
  ggtitle("Default (%) vs. Assets") + 
  annotate("text", label = "Pre-2007", x = 0.3, y = 0.3, size = 4, colour = "black") +
  annotate("text", label = "2007-2009", x = 0.3, y = 0.25, size = 4, colour = "red") +
  scale_y_continuous(limit = c(0, 0.4), expand = c(0, 0))  +  xlab("Savings (a)")  +  ylab("%") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot

g_pricing_cycle = ggplot(data = data.frame(PP=PP, age=age_vector), aes(x=age)) + 
  geom_line(aes(y = PP), size=size_line_doc) + ggtitle("Mortgage Pricing")   + 
  scale_y_continuous(expand = c(0, 0)) + ylab("$") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/default.pdf",width=6.5,height=2.5)
grid.arrange(arrangeGrob(g_def_cycle + theme(legend.position="none"),
                         g_pricing_cycle + theme(legend.position="none"),
                         nrow = 1, ncol = 2, widths=c(1,1)), nrow = 1, ncol = 1)
dev.off()


pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/default.pdf",width=4,height=2.5)
grid.arrange(arrangeGrob(g_def_cycle + theme(legend.position="none"),
                         nrow = 1, ncol = 1), nrow = 1, ncol = 1)
dev.off()

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/default_tr.pdf",width=7,height=4.5)
grid.arrange(arrangeGrob(g_def_cycle_tr + theme(legend.position="none"),
                         g_def_sav_tr + theme(legend.position="none"),
                         nrow = 1, ncol = 2),
             arrangeGrob(g_eq_tr, rectGrob(gp=gpar(col=NA)),
                         nrow = 1, ncol = 2), 
             nrow = 2, ncol = 1)
dev.off()

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/default_tr_doc.pdf",width=6,height=3)
ggplot(data = data.frame(D=D, D1=D1, age=age_vector), aes(x=age)) + 
  geom_line(aes(y = 100*D, linetype="Steady State", colour="Steady State"), size=size_line_doc) + 
  geom_line(aes(y = 100*D1, linetype="Great Recession", colour="Great Recession"), size=size_line_doc) + 
  scale_y_continuous(limit = c(0, 40), expand = c(0, 0))  +  xlab("Age")+  ylab("%") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot + 
  scale_linetype_manual("", breaks = c("Steady State", "Great Recession"), values = c("Steady State"="solid", "Great Recession"="dashed"))+
  scale_color_manual("", breaks = c("Steady State", "Great Recession"), values = c("Steady State"="black", "Great Recession"="red")) +
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))
dev.off()


pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/default_tr.pdf",width=7,height=2.5)
grid.arrange(arrangeGrob(g_eq_tr + theme(legend.position="none"),
                         g_def_cycle_tr + theme(legend.position="none"),
                         nrow = 1, ncol = 2), nrow = 1, ncol = 1)
dev.off()

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/default_sav_tr.pdf",width=4,height=2.5)
grid.arrange(arrangeGrob(g_def_sav_tr + theme(legend.position="none"),
                         nrow = 1, ncol = 1), nrow = 1, ncol = 1)
dev.off()



g_ref_cycle = ggplot(data = data.frame(R=RR, age=age_vector), aes(x=age)) + 
  geom_line(aes(y = 100*R), size=size_line_doc) + 
  xlab("Age") + 
  scale_y_continuous(limit = c(0, 20), expand = c(0, 0))  +  ylab("%") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/refinancing_doc.pdf",width=3.75,height=2.5)
grid.arrange(arrangeGrob(g_ref_cycle + theme(legend.position="none"),
                         nrow = 1, ncol = 1), nrow = 1, ncol = 1)
dev.off()


g_def_cycle = ggplot(data = data.frame(D=D, age=age_vector), aes(x=age)) + 
  geom_line(aes(y = 100*D), size=size_line_doc) + 
  xlab("Age") + 
  scale_y_continuous(limit = c(0, 15), expand = c(0, 0))  +  ylab("%") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot

pequity_def = ggplot(data = data.frame(dist=dist_equity, grid=seq(0,1,length.out =100)), aes(x=grid)) +
              scale_x_continuous(expand = c(0, 0)) + theme_plot + ylab("%") + xlab("Equity (%)") +
              geom_line(data = data.frame(dist=100*default_dist, grid=100*seq(-0.3,1,length.out =26)), aes(x=grid, y=dist), size=size_line_doc, linetype = "solid")

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/default_doc.pdf",width=7,height=2.5)
grid.arrange(arrangeGrob(g_def_cycle + theme(legend.position="none"),
                         pequity_def + theme(legend.position="none"),
                         nrow = 1, ncol = 2), nrow = 1, ncol = 1)
dev.off()

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/default_equity_doc.pdf",width=3.75,height=2.5)
grid.arrange(arrangeGrob(pequity_def + theme(legend.position="none"),
                         nrow = 1, ncol = 1), nrow = 1, ncol = 1)
dev.off()


pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/equity_cycle.pdf",width=6.5,height=3.5)
ggplot(data = data.frame(dist=eq_cycle, grid=seq(24,80,length.out =11)), aes(x=grid)) + 
  geom_line(aes(y = 100*dist, linetype = "Data", color = "Data"), size=size_line_doc) + 
  geom_line(data = data.frame(dist=equity_cycle, grid=seq(24,80,length.out =T)), aes(x=grid, y=100*dist, linetype = "Model", color = "Model"), size=size_line_doc) +
  scale_x_continuous(expand = c(0, 0)) + theme_plot + ylab("Equity (%)") + xlab("Age") + 
  scale_y_continuous(expand = c(0, 0), limit = c(0, 100)) +
  scale_linetype_manual("", breaks = c("Model", "Data"), values = c("Model"="solid", "Data"="dashed"))+
  scale_color_manual("", breaks = c("Model", "Data"), values = c("Model"="red", "Data"="black")) +
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))
dev.off()



pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/renters_cycle.pdf",width=5.7,height=2.7)
ggplot(data = aggregates, aes(x=age)) + 
  geom_line(aes(y = NR, linetype = "Model", color = "Model"), size=size_line_doc) + 
  geom_line(aes(y = Rp_emp, linetype = "Data", color = "Data"), size=size_line_doc) + 
  ylab("%") + xlab("Age") + 
  scale_y_continuous(limit = c(0, 1), expand = c(0, 0))  + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot +
  scale_linetype_manual("", breaks = c("Model", "Data"), values = c("Model"="solid", "Data"="dashed"))+
  scale_color_manual("", breaks = c("Model", "Data"), values = c("Model"="red", "Data"="black")) +
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))
dev.off()


# Transitional Dynamics
g_price_dyn = ggplot(data = data.frame(Ppathbail=c(1, Ppath_tr), Ppathbase=c(1, Ppath_pol), Ppathdef=c(1, Ppath_defpol), year=seq(2006,2006+yearspp*Ttrans, length.out=Ttrans+1)), aes(x=year)) + 
  geom_line(aes(y = Ppathbase, colour="Baseline Policy", linetype="Baseline Policy"), size=size_line_doc) + ggtitle("House Prices")   + 
  geom_line(aes(y = Ppathdef, colour="Subsidy-only Policy", linetype="Subsidy-only Policy"), size=size_line_doc) + 
  geom_line(aes(y = Ppathbail, colour="Bailout-only Policy", linetype="Bailout-only Policy"), size=size_line_doc) + 
  xlab("Year") + 
  theme(legend.position = c(0.8, 0.3)) +
  scale_colour_manual("", breaks = c("Baseline Policy", "Subsidy-only Policy", "Bailout-only Policy"), values = c("Baseline Policy" = "black", "Subsidy-only Policy" = "blue", "Bailout-only Policy" = "red")) +
  scale_linetype_manual("", breaks = c("Baseline Policy", "Subsidy-only Policy", "Bailout-only Policy"), values = c("Baseline Policy" = "solid", "Subsidy-only Policy" = "solid", "Bailout-only Policy" = "solid")) +
  scale_y_continuous(expand = c(0, 0), limit = c(0.7, 1.05)) + ylab("$") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2006,2006+yearspp*Ttrans, length.out=5)) + theme_plot +
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))


g_price_dyn1 = ggplot(data = data.frame(Ppathbail=c(1, Ppath_tr), Ppathbase=c(1, Ppath_pol), Ppathdef=c(1, Ppath_defpol), year=seq(2006,2006+yearspp*Ttrans, length.out=Ttrans+1)), aes(x=year)) + 
  geom_line(aes(y = Ppathbase, colour="Baseline Policy", linetype="Baseline Policy"), size=size_line_doc) + ggtitle("House Prices")   + 
  xlab("Year") + 
  scale_colour_manual("", breaks = c("Baseline Policy", "Subsidy-only Policy", "Bailout-only Policy"), values = c("Baseline Policy" = "black", "Subsidy-only Policy" = "blue", "Bailout-only Policy" = "red")) +
  scale_linetype_manual("", breaks = c("Baseline Policy", "Subsidy-only Policy", "Bailout-only Policy"), values = c("Baseline Policy" = "solid", "Subsidy-only Policy" = "solid", "Bailout-only Policy" = "solid")) +
  scale_y_continuous(expand = c(0, 0), limit = c(0.7, 1.05)) + ylab("$") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2006,2006+yearspp*Ttrans, length.out=5)) + theme_plot +
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))

g_price_dyn2 = ggplot(data = data.frame(Ppathbail=c(1, Ppath_tr), Ppathbase=c(1, Ppath_pol), Ppathdef=c(1, Ppath_defpol), year=seq(2006,2006+yearspp*Ttrans, length.out=Ttrans+1)), aes(x=year)) + 
  geom_line(aes(y = Ppathbase, colour="Baseline Policy", linetype="Baseline Policy"), size=size_line_doc*0.8) + ggtitle("House Prices")   + 
  geom_line(aes(y = Ppathdef, colour="Subsidy-only Policy", linetype="Subsidy-only Policy"), size=size_line_doc) + 
  xlab("Year") + 
  scale_colour_manual("", breaks = c("Baseline Policy", "Subsidy-only Policy", "Bailout-only Policy"), values = c("Baseline Policy" = "black", "Subsidy-only Policy" = "blue", "Bailout-only Policy" = "red")) +
  scale_linetype_manual("", breaks = c("Baseline Policy", "Subsidy-only Policy", "Bailout-only Policy"), values = c("Baseline Policy" = "solid", "Subsidy-only Policy" = "solid", "Bailout-only Policy" = "solid")) +
  scale_y_continuous(expand = c(0, 0), limit = c(0.7, 1.05)) + ylab("$") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2006,2006+yearspp*Ttrans, length.out=5)) + theme_plot +
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))

g_price_dyn3 = ggplot(data = data.frame(Ppathbail=c(1, Ppath_tr), Ppathbase=c(1, Ppath_pol), Ppathdef=c(1, Ppath_defpol), year=seq(2006,2006+yearspp*Ttrans, length.out=Ttrans+1)), aes(x=year)) + 
  geom_line(aes(y = Ppathbase, colour="Baseline Policy", linetype="Baseline Policy"), size=0.8*size_line_doc) + ggtitle("House Prices")   + 
  geom_line(aes(y = Ppathdef, colour="Subsidy-only Policy", linetype="Subsidy-only Policy"), size=0.8*size_line_doc) + 
  geom_line(aes(y = Ppathbail, colour="Bailout-only Policy", linetype="Bailout-only Policy"), size=size_line_doc) + 
  xlab("Year") + 
  scale_colour_manual("", breaks = c("Baseline Policy", "Subsidy-only Policy", "Bailout-only Policy"), values = c("Baseline Policy" = "black", "Subsidy-only Policy" = "blue", "Bailout-only Policy" = "red")) +
  scale_linetype_manual("", breaks = c("Baseline Policy", "Subsidy-only Policy", "Bailout-only Policy"), values = c("Baseline Policy" = "solid", "Subsidy-only Policy" = "solid", "Bailout-only Policy" = "solid")) +
  scale_y_continuous(expand = c(0, 0), limit = c(0.7, 1.05)) + ylab("$") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2006,2006+yearspp*Ttrans, length.out=5)) + theme_plot +
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))

mylegend_pricedyn1<-g_legend(g_price_dyn1) 
mylegend_pricedyn2<-g_legend(g_price_dyn2) 
mylegend_pricedyn3<-g_legend(g_price_dyn3) 

g_default_dyn = ggplot(data = data.frame(defpathbail=100*c(default_rate, defpath_tr), defpath=100*c(default_rate, defpath_pol), defpathpol=100*c(default_rate, defpath_defpol), year=seq(2006,2006+yearspp*Ttrans, length.out=Ttrans+1)), aes(x=year)) + 
  geom_line(aes(y = defpath, colour="Baseline Policy", linetype="Baseline Policy"), size=size_line_doc) + ggtitle("Foreclosure Rate (%)")   + 
  geom_line(aes(y = defpathpol, colour="Subsidy-only Policy", linetype="Subsidy-only Policy"), size=size_line_doc) + 
  geom_line(aes(y = defpathbail, colour="Bailout-only Policy", linetype="Bailout-only Policy"), size=size_line_doc) + ggtitle("Foreclosure Rate (%)")   + 
  theme(legend.position = c(0.8, 0.8)) +
  scale_colour_manual("", breaks = c("Baseline Policy", "Subsidy-only Policy", "Bailout-only Policy"), values = c("Baseline Policy" = "black", "Subsidy-only Policy" = "blue", "Bailout-only Policy" = "red")) +
  scale_linetype_manual("", breaks = c("Baseline Policy", "Subsidy-only Policy", "Bailout-only Policy"), values = c("Baseline Policy" = "solid", "Subsidy-only Policy" = "solid", "Bailout-only Policy" = "solid")) +
  scale_y_continuous(expand = c(0, 0), limit = c(0.00, 11)) + ylab("%") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2006,2006+yearspp*Ttrans, length.out=5)) + theme_plot + 
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))

g_default_dyn1 = ggplot(data = data.frame(defpathbail=100*c(default_rate, defpath_tr), defpath=100*c(default_rate, defpath_pol), defpathpol=100*c(default_rate, defpath_defpol), year=seq(2006,2006+yearspp*Ttrans, length.out=Ttrans+1)), aes(x=year)) + 
  geom_line(aes(y = defpath, colour="Baseline Policy", linetype="Baseline Policy"), size=size_line_doc) + ggtitle("Foreclosure Rate (%)")   + 
  theme(legend.position = c(0.8, 0.8)) +
  scale_colour_manual("", breaks = c("Baseline Policy", "Subsidy-only Policy", "Bailout-only Policy"), values = c("Baseline Policy" = "black", "Subsidy-only Policy" = "blue", "Bailout-only Policy" = "red")) +
  scale_linetype_manual("", breaks = c("Baseline Policy", "Subsidy-only Policy", "Bailout-only Policy"), values = c("Baseline Policy" = "solid", "Subsidy-only Policy" = "solid", "Bailout-only Policy" = "solid")) +
  scale_y_continuous(expand = c(0, 0), limit = c(0.00, 11)) + ylab("%") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2006,2006+yearspp*Ttrans, length.out=5)) + theme_plot + 
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))

g_default_dyn2 = ggplot(data = data.frame(defpathbail=100*c(default_rate, defpath_tr), defpath=100*c(default_rate, defpath_pol), defpathpol=100*c(default_rate, defpath_defpol), year=seq(2006,2006+yearspp*Ttrans, length.out=Ttrans+1)), aes(x=year)) + 
  geom_line(aes(y = defpath, colour="Baseline Policy", linetype="Baseline Policy"), size=size_line_doc*0.8) + ggtitle("Foreclosure Rate (%)")   + 
  geom_line(aes(y = defpathpol, colour="Subsidy-only Policy", linetype="Subsidy-only Policy"), size=size_line_doc) + 
  theme(legend.position = c(0.8, 0.8)) +
  scale_colour_manual("", breaks = c("Baseline Policy", "Subsidy-only Policy", "Bailout-only Policy"), values = c("Baseline Policy" = "black", "Subsidy-only Policy" = "blue", "Bailout-only Policy" = "red")) +
  scale_linetype_manual("", breaks = c("Baseline Policy", "Subsidy-only Policy", "Bailout-only Policy"), values = c("Baseline Policy" = "solid", "Subsidy-only Policy" = "solid", "Bailout-only Policy" = "solid")) +
  scale_y_continuous(expand = c(0, 0), limit = c(0.00, 11)) + ylab("%") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2006,2006+yearspp*Ttrans, length.out=5)) + theme_plot + 
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))

g_default_dyn3 = ggplot(data = data.frame(defpathbail=100*c(default_rate, defpath_tr), defpath=100*c(default_rate, defpath_pol), defpathpol=100*c(default_rate, defpath_defpol), year=seq(2006,2006+yearspp*Ttrans, length.out=Ttrans+1)), aes(x=year)) + 
  geom_line(aes(y = defpath, colour="Baseline Policy", linetype="Baseline Policy"), size=0.8*size_line_doc) + ggtitle("Foreclosure Rate (%)")   + 
  geom_line(aes(y = defpathpol, colour="Subsidy-only Policy", linetype="Subsidy-only Policy"), size=0.8*size_line_doc) + 
  geom_line(aes(y = defpathbail, colour="Bailout-only Policy", linetype="Bailout-only Policy"), size=size_line_doc) + ggtitle("Foreclosure Rate (%)")   + 
  theme(legend.position = c(0.8, 0.8)) +
  scale_colour_manual("", breaks = c("Baseline Policy", "Subsidy-only Policy", "Bailout-only Policy"), values = c("Baseline Policy" = "black", "Subsidy-only Policy" = "blue", "Bailout-only Policy" = "red")) +
  scale_linetype_manual("", breaks = c("Baseline Policy", "Subsidy-only Policy", "Bailout-only Policy"), values = c("Baseline Policy" = "solid", "Subsidy-only Policy" = "solid", "Bailout-only Policy" = "solid")) +
  scale_y_continuous(expand = c(0, 0), limit = c(0.00, 11)) + ylab("%") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2006,2006+yearspp*Ttrans, length.out=5)) + theme_plot + 
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/transitions.pdf",width=6.5,height=2.5)
grid.arrange(arrangeGrob(g_price_dyn,
                         g_default_dyn,
                         nrow = 1, ncol = 2, widths=c(1,1)), nrow = 1, ncol = 1)
dev.off()

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/transitions0.pdf",width=6.5,height=2.5)
grid.arrange(arrangeGrob(g_price_dyn1 + theme(legend.position="none"),
                         g_default_dyn1 + theme(legend.position="none"),
                         nrow = 1, ncol = 2, widths=c(1,1)),
             nrow = 1, ncol = 1)
dev.off()


pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/transitions1.pdf",width=6.5,height=3.5)
grid.arrange(arrangeGrob(g_price_dyn1 + theme(legend.position="none"),
                         g_default_dyn1 + theme(legend.position="none"),
                         nrow = 1, ncol = 2, widths=c(1,1)), 
             mylegend_pricedyn1,
             nrow = 2, ncol = 1, heights=c(1,0.4))
dev.off()

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/transitions2.pdf",width=6.5,height=3.5)
grid.arrange(arrangeGrob(g_price_dyn2 + theme(legend.position="none"),
                         g_default_dyn2 + theme(legend.position="none"),
                         nrow = 1, ncol = 2, widths=c(1,1)),
             mylegend_pricedyn2,
             nrow = 2, ncol = 1, heights=c(1,0.4))
dev.off()

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/transitions3.pdf",width=6.5,height=3.5)
grid.arrange(arrangeGrob(g_price_dyn3 + theme(legend.position="none"),
                         g_default_dyn3 + theme(legend.position="none"),
                         nrow = 1, ncol = 2, widths=c(1,1)), 
             mylegend_pricedyn3,
             nrow = 2, ncol = 1, heights=c(1,0.4))
dev.off()


g_price_dyn_doc = ggplot(data = data.frame(Ppathbail=c(1, Ppath_tr), Ppathbase=c(1, Ppath_pol), Ppathdef=c(1, Ppath_defpol), year=seq(2006,2006+yearspp*Ttrans, length.out=Ttrans+1)), aes(x=year)) + 
  geom_line(aes(y = Ppathbase, colour="Baseline Policy", linetype="Baseline Policy"), size=0.8*size_line_doc) + 
  geom_line(aes(y = Ppathdef, colour="Subsidy-only Policy", linetype="Subsidy-only Policy"), size=0.8*size_line_doc) + 
  geom_line(aes(y = Ppathbail, colour="Bailout-only Policy", linetype="Bailout-only Policy"), size=size_line_doc) + 
  xlab("Year") + 
  scale_colour_manual("", breaks = c("Baseline Policy", "Subsidy-only Policy", "Bailout-only Policy"), values = c("Baseline Policy" = "black", "Subsidy-only Policy" = "blue", "Bailout-only Policy" = "red")) +
  scale_linetype_manual("", breaks = c("Baseline Policy", "Subsidy-only Policy", "Bailout-only Policy"), values = c("Baseline Policy" = "solid", "Subsidy-only Policy" = "dashed", "Bailout-only Policy" = "twodash")) +
  scale_y_continuous(expand = c(0, 0), limit = c(0.7, 1.05)) + ylab("House Prices ($)") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2006,2006+yearspp*Ttrans, length.out=5)) + theme_plot +
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))


g_default_dyn_doc = ggplot(data = data.frame(defpathbail=100*c(default_rate, defpath_tr), defpath=100*c(default_rate, defpath_pol), defpathpol=100*c(default_rate, defpath_defpol), year=seq(2006,2006+yearspp*Ttrans, length.out=Ttrans+1)), aes(x=year)) + 
  geom_line(aes(y = defpath, colour="Baseline Policy", linetype="Baseline Policy"), size=0.8*size_line_doc) + 
  geom_line(aes(y = defpathpol, colour="Subsidy-only Policy", linetype="Subsidy-only Policy"), size=0.8*size_line_doc) + 
  geom_line(aes(y = defpathbail, colour="Bailout-only Policy", linetype="Bailout-only Policy"), size=size_line_doc) + 
  theme(legend.position = c(0.8, 0.8)) + 
  scale_colour_manual("", breaks = c("Baseline Policy", "Subsidy-only Policy", "Bailout-only Policy"), values = c("Baseline Policy" = "black", "Subsidy-only Policy" = "blue", "Bailout-only Policy" = "red")) +
  scale_linetype_manual("", breaks = c("Baseline Policy", "Subsidy-only Policy", "Bailout-only Policy"), values = c("Baseline Policy" = "solid", "Subsidy-only Policy" = "dashed", "Bailout-only Policy" = "twodash")) +
  scale_y_continuous(expand = c(0, 0), limit = c(0.00, 11)) + ylab("Foreclosure Rate (%)") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2006,2006+yearspp*Ttrans, length.out=5)) + theme_plot + 
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))

mylegend_pricedyn_doc <- g_legend(g_price_dyn_doc) 

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/transitions_doc.pdf",width=7.5,height=3.5)
grid.arrange(arrangeGrob(g_price_dyn_doc + theme(legend.position="none"),
                         g_default_dyn_doc + theme(legend.position="none"),
                         nrow = 1, ncol = 2, widths=c(1,1)), 
             mylegend_pricedyn_doc,
             nrow = 2, ncol = 1, heights=c(1,0.3))
dev.off()

g_default_dyn = ggplot(data = data.frame(defpath=100*c(default_rate, defpath_tr), defpathpol=100*c(default_rate, defpath_pol), year=seq(2006,2006+yearspp*Ttrans, length.out=Ttrans+1)), aes(x=year)) + 
  geom_line(aes(y = defpath), size=size_line_doc) + ggtitle("Foreclosure Rate (%)")   + 
  xlab("Year") + 
  scale_y_continuous(expand = c(0, 0), limit = c(0.00, 10.5)) + ylab("%") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2006,2006+yearspp*Ttrans, length.out=5)) + theme_plot


# Life-cycle patterns
pR = ggplot(data = aggregates, aes(x=age)) + 
  geom_line(aes(y = 100*NR, colour="Model", linetype="Model"), size=size_line_doc) + 
  geom_line(aes(y = 100*Rp_emp, colour="SCF 2007", linetype="SCF 2007"), size=size_line_doc) + 
  ggtitle("Net Renters by Age (%)") + ylab("%") + xlab("Age") + 
  scale_y_continuous(limit = c(0, 110), expand = c(0, 0))  + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot +
  scale_colour_manual("", breaks = c("Model", "SCF 2007"), values = c("Model" = "red", "SCF 2007" = "black")) +
  scale_linetype_manual("", breaks = c("Model", "SCF 2007"), values = c("Model" = "solid", "SCF 2007" = "dashed"))+
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))

pA =  ggplot(data = data.frame(A=A, age=age_vector), aes(x=age)) + 
  geom_line(aes(y = A, colour="Model", linetype="Model"), size=size_line_doc) + 
  ggtitle("Assets")   +  xlab("Age") + ylab("Assets ($)") + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot +
  scale_colour_manual("", breaks = c("Model", "Data"), values = c("Model" = "red", "Data" = "black")) +
  scale_linetype_manual("", breaks = c("Model", "Data"), values = c("Model" = "solid", "Data" = "dashed"))+
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))


pM = ggplot(data = data.frame(M=M, age=age_vector), aes(x=age)) + 
  geom_line(aes(y = M, colour="Model", linetype="Model"), size=size_line_doc) + 
  scale_y_continuous(expand = c(0, 0), limit = c(0, 1.2*max(M))) + 
  scale_x_continuous(expand = c(0, 0)) + 
  ggtitle("Mortgage")  + theme_plot + xlab("Age") + ylab("Mortgages (m)") + 
  scale_colour_manual("", breaks = c("Model", "Data"), values = c("Model" = "red", "Data" = "black")) +
  scale_linetype_manual("", breaks = c("Model", "Data"), values = c("Model" = "solid", "Data" = "dashed"))+
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))


pH = ggplot(data = aggregates, aes(x=age)) + 
  geom_line(aes(y = VU, colour="Model", linetype="Model"), size=size_line_doc) + 
  ggtitle("Owner Occupied Housing") +  xlab("Age") + ylab("Housing (h)") + 
  scale_y_continuous(limit = c(0, 1.9), expand = c(0, 0))  + 
  scale_x_continuous(expand = c(0, 0)) + theme_plot +
  scale_colour_manual("", breaks = c("Model", "Data"), values = c("Model" = "red", "Data" = "black")) +
  scale_linetype_manual("", breaks = c("Model", "Data"), values = c("Model" = "solid", "Data" = "dashed"))+
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))

mylegend_pR <- g_legend(pR) 

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/lifecycle_paths.pdf",width=7,height=8.5)
grid.arrange(arrangeGrob(pH + theme(legend.position="none"), pM + theme(legend.position="none"), ncol = 2, nrow = 1),
             arrangeGrob(peq_cycle_doc + theme(legend.position="none"), pA + theme(legend.position="none"), ncol = 2, nrow = 1),
             arrangeGrob(mylegend_pR, ncol = 1, nrow = 1),
             ncol = 1, nrow = 3, heights=c(1,1,0.2))
dev.off()

# Who are the refinancers?
pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/refinance_paths.pdf",width=6.5,height=3)
ggplot(data = data.frame(RR=RRup, age=age_vector), aes(x=age)) + 
  geom_line(aes(y = RR, colour="Refinance up", linetype = "Refinance up"), size=size_line_doc) + 
  geom_line(data = data.frame(RR=RRdn, age=age_vector), aes(y = RR, colour="Refinance down", linetype = "Refinance down"), size=size_line) + 
  scale_y_continuous(limit = c(0, 0.2), expand = c(0, 0)) + ylab("%") + xlab("Age") +
  scale_colour_manual("", breaks = c("Refinance up", "Refinance down"), values = c("Refinance up" = "blue", "Refinance down" = "red")) +
  scale_linetype_manual("", breaks = c("Refinance up", "Refinance down"), values = c("Refinance up" = "solid", "Refinance down" = "dashed"))+
  scale_x_continuous(expand = c(0, 0)) + theme_plot +
  theme(legend.key.width=unit(3.5,"lines"), panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))
dev.off()


# Who changed default decision after shock?

p_t_change = ggplot(data = data.frame(t_changers=t_changers, age=age_vector), aes(x=age)) +
  geom_line(aes(y = t_changers), size=size_line) +
  scale_y_continuous(limit = c(0, 0.4), expand = c(0, 0))  +  ylab("%") +
  scale_x_continuous(expand = c(0, 0)) + theme_plot




surv = c(0.997698, 0.997331,0.997287,0.997299,0.997232,0.997106,0.996975,0.996818,
         0.996562,0.996189,0.995667,0.994902,0.993824,0.992517,0.990992,0.989251,0.987264,
         0.984991,0.982468,0.980156,0.977663,0.974437,0.970300,0.965084,0.958700,0.950622,
         0.941064,0.930068,0.915904,0.898216,0.104000)

av_rate = array(0,T-2)
pipol = array(0,T-2)

for(it in 1:(T-2)){
  print(it)
  for(iy in 1:ny){
    for(ia in 1:na){
      for(im in 1:nm){
        for(ih in 1:nh){
          for(id in 1:nd){
            
            if(Policym[it, iy, ia, im, ih, id] > 0 & Renew[it, iy, ia, im, ih, id] == 1){
              age = it
              
              iap = Policya[it, iy, ia, im, ih, id]+1
              imp = Policym[it, iy, ia, im, ih, id]+1
              ihp = Policyh[it, iy, ia, im, ih, id]+1
              
              P = Pricing[it, iy, iap, imp, ihp, id]
              
              implied = function(x){
                age = 1
                res = 0
                hasta = min(12.5, 30-age)
                # hasta = 12.5
                for(i in 1:hasta){
                  res = res + surv[age + i - 1]/((1+x)^(i-1))
                }
                
                return(res - P)
              }
              
              x = nleqslv(0.1, implied)
              
              av_rate[it] = av_rate[it] + Pcond[it, iy, ia, im, ih, id]*x$x
              pipol[it] = pipol[it] + Pcond[it, iy, ia, im, ih, id]
              
            }
          }
        }
      }
    }
  }
}

for(it in 1:(T-2)){
  av_rate[it] = av_rate[it]/pipol[it]
}

plot(seq(20,60,length.out=21), av_rate[1:21])