
# -----------------------------------------------------------#
#                  Prelimiaries                              #
# -----------------------------------------------------------#

library("Rcpp")
library("MASS")
library("Rtauchen")
library("ggplot2")
library("foreign")
library("scales")
library("gridExtra")
library("fANCOVA")
library("grid")
library("nleqslv")


# ------------------------#
#     Experiment #1       #
# ------------------------#

# Increase of interest rate of 5% yearly during 2 periods
# Decrease of income of 20% on first period

q = 0.24

price_ext = array(0, 7)
rent_ext  = array(0, 7)
def_ext   = array(0,7)

price_ext[1] = 1.0;
price_ext[2] = 0.755000;
price_ext[3] = 0.955;
price_ext[4] = 0.974;
price_ext[5] = 0.991;
price_ext[6] = 0.998;
price_ext[7] = 1.0;

rent_ext[1] = q;
rent_ext[2] = q-0.044;
rent_ext[3] = q-0.015;
rent_ext[4] = q-0.012;
rent_ext[5] = q-0.003;
rent_ext[6] = q+0.0;
rent_ext[7] = q+0.0;

def_ext[1] = 0.078;
def_ext[2] = 0.311;
def_ext[3] = 0.163;
def_ext[4] = 0.150;
def_ext[5] = 0.106;
def_ext[6] = 0.081;
def_ext[7] = 0.079;

price_noext = array(0, 7)
rent_noext  = array(0, 7)
def_noext   = array(0, 7)

price_noext[1] = 1.0;
price_noext[2] = 0.78;
price_noext[3] = 0.977;
price_noext[4] = 0.987;
price_noext[5] = 0.998;
price_noext[6] = 1.0;
price_noext[7] = 1.0;

rent_noext[1] = q;
rent_noext[2] = q-0.042;
rent_noext[3] = q-0.012;
rent_noext[4] = q-0.01;
rent_noext[5] = q-0.0;
rent_noext[6] = q+0.0;
rent_noext[7] = q+0.0;

def_noext[1] = 0.078;
def_noext[2] = 0.170;
def_noext[3] = 0.105;
def_noext[4] = 0.104;
def_noext[5] = 0.094;
def_noext[6] = 0.072;
def_noext[7] = 0.076;

datos = data.frame(period = seq(0,6,by=1), 
                   price_ext = price_ext, price_noext = price_noext,
                   rent_ext = rent_ext, rent_noext = rent_noext,
                   def_ext = def_ext, def_noext = def_noext)




size_factor = 1
size_line = 1.5
annotate_size = 4
beamer_font = 14
theme_plot = theme(panel.background = element_rect(fill='white', colour='black'),
                   panel.grid.major = element_line(colour = "lightgray"),
                   plot.title = element_text(hjust = 0.5))

pprice = ggplot(data = datos, aes(x=period)) + 
  geom_line(aes(y = price_ext, linetype="Foreclosure Externalities", colour="Foreclosure Externalities"), size = size_line) + 
  geom_line(aes(y = price_noext, linetype="No Externalities", colour="No Externalities"), size = size_line) + 
  ggtitle("House Price") + labs(x = "", y = "") +
  scale_linetype_manual("", breaks = c("Foreclosure Externalities", "No Externalities"), values = c("Foreclosure Externalities"="solid", "No Externalities"="dashed"))+
  scale_colour_manual("", breaks = c("Foreclosure Externalities", "No Externalities"), values = c("Foreclosure Externalities"="black", "No Externalities"="blue"))  +
  scale_x_continuous(breaks = seq(0,6,by=1), expand = c(0, 0)) +
  scale_y_continuous(limit = c(0.72, 1.08), expand = c(0, 0))  + theme_plot +
  geom_hline(yintercept = price_noext[1], linetype = "longdash") +
  annotate("text", label = paste("Price fall:", paste0((price_ext[2]-price_ext[1])*100/price_ext[1]),"%"), x = 4, y = 0.85, size = annotate_size, colour = "black") +
  annotate("text", label = paste("Externality:", paste0(((price_noext[2]-price_noext[1])*100/price_noext[1])-((price_ext[2]-price_ext[1])*100/price_ext[1])),"%"), x = 4, y = 0.82, size = annotate_size, colour = "black") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))


prent = ggplot(data = datos, aes(x=period)) + 
  geom_line(aes(y = rent_ext), colour="black", size = size_line) + 
  geom_line(aes(y = rent_noext), linetype="dashed", colour="blue", size = size_line) + 
  ggtitle("Rental Price") + labs(x = "", y = "") +
  annotate("text", label = paste("Rent fall:", paste0(round((rent_ext[2]-rent_ext[1])*100/rent_ext[1],1)),"%"), x = 4, y = 0.2, size = annotate_size, colour = "black") +
  geom_hline(yintercept = rent_noext[1], linetype = "longdash") +
  scale_x_continuous(breaks = seq(0,6,by=1), expand = c(0, 0)) +
  scale_y_continuous(limit = c(0.17, 0.25), expand = c(0, 0))  + theme_plot

pdefault = ggplot(data = datos, aes(x=period)) + 
  geom_line(aes(y = def_ext), colour="black", size = size_line) + 
  geom_line(aes(y = def_noext), linetype="dashed", colour="blue", size = size_line) + 
  ggtitle("Default Rate") + labs(x = "", y = "") +
  geom_hline(yintercept = def_noext[1], linetype = "longdash") +
  scale_x_continuous(breaks = seq(0,6,by=1), expand = c(0, 0)) +
  scale_y_continuous(limit = c(0.0, 0.35), expand = c(0, 0))  + theme_plot

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(pprice) 

gra = grid.arrange(arrangeGrob(pprice + theme(legend.position="none"), prent + theme(legend.position="none"), ncol = 2, nrow = 1),
                   arrangeGrob(pdefault + theme(legend.position="none"), arrangeGrob(textGrob(""),textGrob("Income shock: -20%"), textGrob("Interest rate shock: +5% yearly"), mylegend, ncol = 1, nrow = 4, heights=c(0.7,0.1,0.1,1)), ncol = 2, nrow = 1),
                   ncol = 1, nrow = 2, heights=c(1,1))
print(gra)

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/Transition_experiment_1.pdf",width=7,height=7)
grid.arrange(arrangeGrob(pprice + theme(legend.position="none"), prent + theme(legend.position="none"), ncol = 2, nrow = 1),
             arrangeGrob(pdefault + theme(legend.position="none"), arrangeGrob(textGrob(""),textGrob("Income shock: -20%"), textGrob("Interest rate shock: +5% yearly"), mylegend, ncol = 1, nrow = 4, heights=c(0.7,0.1,0.1,1)), ncol = 2, nrow = 1),
             ncol = 1, nrow = 2, heights=c(1,1))
print(gra)
dev.off()

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/Transition_experiment_1_beamer.pdf",width=8.5,height=5.2)
grid.arrange(arrangeGrob(pprice + theme(legend.position="none", text = element_text(size=beamer_font)), prent + theme(legend.position="none", text = element_text(size=beamer_font)), ncol = 2, nrow = 1),
             arrangeGrob(pdefault + theme(legend.position="none", text = element_text(size=beamer_font)), arrangeGrob(textGrob(""),textGrob("Income shock: -20%"), textGrob("Interest rate shock: +5% yearly"), mylegend, ncol = 1, nrow = 4, heights=c(0.7,0.1,0.1,1)), ncol = 2, nrow = 1),
             ncol = 1, nrow = 2, heights=c(1,1))
print(gra)
dev.off()


# ------------------------#
#     Experiment #2       #
# ------------------------#

# Taste shock: increase ppsi from 0.65 to 0.73

q = 0.236237

price_ext = array(0, 7)
rent_ext  = array(0, 7)
def_ext   = array(0,7)

price_ext[1] = 1.0;
price_ext[2] = 0.91;
price_ext[3] = 0.987;
price_ext[4] = 1.0;
price_ext[5] = 1.0;
price_ext[6] = 1.0;
price_ext[7] = 1.0;

rent_ext[1] = q;
rent_ext[2] = q-0.066;
rent_ext[3] = q-0.01;
rent_ext[4] = q-0.0;
rent_ext[5] = q-0.0;
rent_ext[6] = q+0.0;
rent_ext[7] = q+0.0;

def_ext[1] = 0.065;
def_ext[2] = 0.268;
def_ext[3] = 0.138;
def_ext[4] = 0.076;
def_ext[5] = 0.057;
def_ext[6] = 0.063;
def_ext[7] = 0.065;

price_noext = array(0, 7)
rent_noext  = array(0, 7)
def_noext   = array(0, 7)

price_noext[1] = 1.0;
price_noext[2] = 0.915;
price_noext[3] = 0.989;
price_noext[4] = 1.0;
price_noext[5] = 1.0;
price_noext[6] = 1.0;
price_noext[7] = 1.0;

rent_noext[1] = q;
rent_noext[2] = q-0.066;
rent_noext[3] = q-0.01;
rent_noext[4] = q-0.0;
rent_noext[5] = q-0.0;
rent_noext[6] = q+0.0;
rent_noext[7] = q+0.0;

def_noext[1] = 0.065;
def_noext[2] = 0.190;
def_noext[3] = 0.118;
def_noext[4] = 0.074;
def_noext[5] = 0.054;
def_noext[6] = 0.062;
def_noext[7] = 0.065;

datos = data.frame(period = seq(0,6,by=1), 
                   price_ext = price_ext, price_noext = price_noext,
                   rent_ext = rent_ext, rent_noext = rent_noext,
                   def_ext = def_ext, def_noext = def_noext)


pprice = ggplot(data = datos, aes(x=period)) + 
  geom_line(aes(y = price_ext, linetype="Foreclosure Externalities", colour="Foreclosure Externalities"), size = size_line) + 
  geom_line(aes(y = price_noext, linetype="No Externalities", colour="No Externalities"), size = size_line) + 
  ggtitle("House Price") + labs(x = "", y = "") +
  scale_linetype_manual("", breaks = c("Foreclosure Externalities", "No Externalities"), values = c("Foreclosure Externalities"="solid", "No Externalities"="dashed"))+
  scale_colour_manual("", breaks = c("Foreclosure Externalities", "No Externalities"), values = c("Foreclosure Externalities"="black", "No Externalities"="blue"))  +
  scale_x_continuous(breaks = seq(0,6,by=1), expand = c(0, 0)) +
  scale_y_continuous(limit = c(0.72, 1.08), expand = c(0, 0))  + theme_plot +
  geom_hline(yintercept = price_noext[1], linetype = "longdash") +
  annotate("text", label = paste("Price fall:", paste0((price_ext[2]-price_ext[1])*100/price_ext[1]),"%"), x = 4, y = 0.85, size = annotate_size, colour = "black") +
  annotate("text", label = paste("Externality:", paste0(((price_noext[2]-price_noext[1])*100/price_noext[1])-((price_ext[2]-price_ext[1])*100/price_ext[1])),"%"), x = 4, y = 0.82, size = annotate_size, colour = "black") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))


prent = ggplot(data = datos, aes(x=period)) + 
  geom_line(aes(y = rent_ext), colour="black", size = size_line) + 
  geom_line(aes(y = rent_noext), linetype="dashed", colour="blue", size = size_line) + 
  ggtitle("Rental Price") + labs(x = "", y = "") +
  annotate("text", label = paste("Rent fall:", paste0(round((rent_ext[2]-rent_ext[1])*100/rent_ext[1],1)),"%"), x = 4, y = 0.2, size = annotate_size, colour = "black") +
  geom_hline(yintercept = rent_noext[1], linetype = "longdash") +
  scale_x_continuous(breaks = seq(0,6,by=1), expand = c(0, 0)) +
  scale_y_continuous(limit = c(0.15, 0.25), expand = c(0, 0))  + theme_plot

pdefault = ggplot(data = datos, aes(x=period)) + 
  geom_line(aes(y = def_ext), colour="black", size = size_line) + 
  geom_line(aes(y = def_noext), linetype="dashed", colour="blue", size = size_line) + 
  ggtitle("Default Rate") + labs(x = "", y = "") +
  geom_hline(yintercept = def_noext[1], linetype = "longdash") +
  scale_x_continuous(breaks = seq(0,6,by=1), expand = c(0, 0)) +
  scale_y_continuous(limit = c(0.0, 0.35), expand = c(0, 0))  + theme_plot

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(pprice) 

gra = grid.arrange(arrangeGrob(pprice + theme(legend.position="none"), prent + theme(legend.position="none"), ncol = 2, nrow = 1),
                   arrangeGrob(pdefault + theme(legend.position="none"), arrangeGrob(textGrob(""),textGrob(""), textGrob("Taste shock"), mylegend, ncol = 1, nrow = 4, heights=c(0.7,0.1,0.1,1)), ncol = 2, nrow = 1),
                   ncol = 1, nrow = 2, heights=c(1,1))
print(gra)

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/Transition_experiment_2.pdf",width=7,height=7)
grid.arrange(arrangeGrob(pprice + theme(legend.position="none"), prent + theme(legend.position="none"), ncol = 2, nrow = 1),
             arrangeGrob(pdefault + theme(legend.position="none"), arrangeGrob(textGrob(""),textGrob(""), textGrob("Taste shock"), mylegend, ncol = 1, nrow = 4, heights=c(0.7,0.1,0.1,1)), ncol = 2, nrow = 1),
             ncol = 1, nrow = 2, heights=c(1,1))
print(gra)
dev.off()

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/Transition_experiment_2.pdf",width=8.5,height=5.2)
grid.arrange(arrangeGrob(pprice + theme(legend.position="none", text = element_text(size=beamer_font)), prent + theme(legend.position="none", text = element_text(size=beamer_font)), ncol = 2, nrow = 1),
             arrangeGrob(pdefault + theme(legend.position="none", text = element_text(size=beamer_font)), arrangeGrob(textGrob(""),textGrob(""), textGrob("Taste shock"), mylegend, ncol = 1, nrow = 4, heights=c(0.7,0.1,0.1,1)), ncol = 2, nrow = 1),
             ncol = 1, nrow = 2, heights=c(1,1))
print(gra)
dev.off()

# ------------------------#
#     Experiment #3       #
# ------------------------#

# Increase of interest rate of 5% yearly during 2 periods
# Decrease of income of 10% on first period

q = 0.24

price_ext = array(0, 7)
rent_ext  = array(0, 7)
def_ext   = array(0,7)

price_ext[1] = 1.0;
price_ext[2] = 0.792;
price_ext[3] = 0.9749;
price_ext[4] = 0.9894;
price_ext[5] = 0.9975;
price_ext[6] = 1.0;
price_ext[7] = 1.0;

rent_ext[1] = q;
rent_ext[2] = q-0.032;
rent_ext[3] = q-0.01;
rent_ext[4] = q-0.006;
rent_ext[5] = q-0.001;
rent_ext[6] = q+0.0;
rent_ext[7] = q+0.0;

def_ext[1] = 0.078;
def_ext[2] = 0.180;
def_ext[3] = 0.110;
def_ext[4] = 0.105;
def_ext[5] = 0.079;
def_ext[6] = 0.078;
def_ext[7] = 0.078;

price_noext = array(0, 7)
rent_noext  = array(0, 7)
def_noext   = array(0, 7)

price_noext[1] = 1.0;
price_noext[2] = 0.815;
price_noext[3] = 0.98;
price_noext[4] = 0.989;
price_noext[5] = 0.997;
price_noext[6] = 1.0;
price_noext[7] = 1.0;

rent_noext[1] = q;
rent_noext[2] = q-0.024;
rent_noext[3] = q-0.01;
rent_noext[4] = q-0.006;
rent_noext[5] = q-0.001;
rent_noext[6] = q+0.0;
rent_noext[7] = q+0.0;

def_noext[1] = 0.078;
def_noext[2] = 0.126;
def_noext[3] = 0.090;
def_noext[4] = 0.104;
def_noext[5] = 0.092;
def_noext[6] = 0.069;
def_noext[7] = 0.075;

datos = data.frame(period = seq(0,6,by=1), 
                   price_ext = price_ext, price_noext = price_noext,
                   rent_ext = rent_ext, rent_noext = rent_noext,
                   def_ext = def_ext, def_noext = def_noext)


pprice = ggplot(data = datos, aes(x=period)) + 
  geom_line(aes(y = price_ext, linetype="Foreclosure Externalities", colour="Foreclosure Externalities"), size = size_line) + 
  geom_line(aes(y = price_noext, linetype="No Externalities", colour="No Externalities"), size = size_line) + 
  ggtitle("House Price") + labs(x = "", y = "") +
  scale_linetype_manual("", breaks = c("Foreclosure Externalities", "No Externalities"), values = c("Foreclosure Externalities"="solid", "No Externalities"="dashed"))+
  scale_colour_manual("", breaks = c("Foreclosure Externalities", "No Externalities"), values = c("Foreclosure Externalities"="black", "No Externalities"="blue"))  +
  scale_x_continuous(breaks = seq(0,6,by=1), expand = c(0, 0)) +
  scale_y_continuous(limit = c(0.72, 1.08), expand = c(0, 0))  + theme_plot +
  geom_hline(yintercept = price_noext[1], linetype = "longdash") +
  annotate("text", label = paste("Price fall:", paste0((price_ext[2]-price_ext[1])*100/price_ext[1]),"%"), x = 4, y = 0.85, size = annotate_size, colour = "black") +
  annotate("text", label = paste("Externality:", paste0(round((price_noext[2]-price_noext[1])*100/price_noext[1])-((price_ext[2]-price_ext[1])*100/price_ext[1]), 1),"%"), x = 4, y = 0.82, size = annotate_size, colour = "black") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.key = element_blank()) + guides(fill = guide_legend(label.position = "bottom"))


prent = ggplot(data = datos, aes(x=period)) + 
  geom_line(aes(y = rent_ext), colour="black", size = size_line) + 
  geom_line(aes(y = rent_noext), linetype="dashed", colour="blue", size = size_line) + 
  ggtitle("Rental Price") + labs(x = "", y = "") +
  annotate("text", label = paste("Rent fall:", paste0(round((rent_ext[2]-rent_ext[1])*100/rent_ext[1],1)),"%"), x = 4, y = 0.2, size = annotate_size, colour = "black") +
  geom_hline(yintercept = rent_noext[1], linetype = "longdash") +
  scale_x_continuous(breaks = seq(0,6,by=1), expand = c(0, 0)) +
  scale_y_continuous(limit = c(0.18, 0.25), expand = c(0, 0))  + theme_plot

pdefault = ggplot(data = datos, aes(x=period)) + 
  geom_line(aes(y = def_ext), colour="black", size = size_line) + 
  geom_line(aes(y = def_noext), linetype="dashed", colour="blue", size = size_line) + 
  ggtitle("Default Rate") + labs(x = "", y = "") +
  geom_hline(yintercept = def_noext[1], linetype = "longdash") +
  scale_x_continuous(breaks = seq(0,6,by=1), expand = c(0, 0)) +
  scale_y_continuous(limit = c(0.0, 0.35), expand = c(0, 0))  + theme_plot

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(pprice) 

gra = grid.arrange(arrangeGrob(pprice + theme(legend.position="none"), prent + theme(legend.position="none"), ncol = 2, nrow = 1),
                   arrangeGrob(pdefault + theme(legend.position="none"), arrangeGrob(textGrob(""),textGrob("Income shock: -10%"), textGrob("Interest rate shock: +5% yearly"), mylegend, ncol = 1, nrow = 4, heights=c(0.7,0.1,0.1,1)), ncol = 2, nrow = 1),
                   ncol = 1, nrow = 2, heights=c(1,1))
print(gra)

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/Transition_experiment_3.pdf",width=7,height=7)
grid.arrange(arrangeGrob(pprice + theme(legend.position="none"), prent + theme(legend.position="none"), ncol = 2, nrow = 1),
             arrangeGrob(pdefault + theme(legend.position="none"), arrangeGrob(textGrob(""),textGrob("Income shock: -10%"), textGrob("Interest rate shock: +5% yearly"), mylegend, ncol = 1, nrow = 4, heights=c(0.7,0.1,0.1,1)), ncol = 2, nrow = 1),
             ncol = 1, nrow = 2, heights=c(1,1))
print(gra)
dev.off()

pdf("/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Document/Figures/Transition_experiment_3.pdf",width=8.5,height=5.2)
grid.arrange(arrangeGrob(pprice + theme(legend.position="none", text = element_text(size=beamer_font)), prent + theme(legend.position="none", text = element_text(size=beamer_font)), ncol = 2, nrow = 1),
             arrangeGrob(pdefault + theme(legend.position="none", text = element_text(size=beamer_font)), arrangeGrob(textGrob(""),textGrob("Income shock: -10%"), textGrob("Interest rate shock: +5% yearly"), mylegend, ncol = 1, nrow = 4, heights=c(0.7,0.1,0.1,1)), ncol = 2, nrow = 1),
             ncol = 1, nrow = 2, heights=c(1,1))
print(gra)
dev.off()
