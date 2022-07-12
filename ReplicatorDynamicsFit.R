library(ggbiplot)
library(ggpubr)
library(tidyr)
library(deSolve)
library(gridExtra)
library(grid)

#data_long is a data frame with patient ID, RTx dosage, Time, Cell Type (MES1 or MES2), and the proportion of this type
data_RTx_4 <- data_long[(data_long$RTxDosage == "0" | data_long$RTxDosage == "4"),]
data_RTx_10 <- data_long[(data_long$RTxDosage == "0" | data_long$RTxDosage == "10"),]

#a,b,c,d fitness matrix parameters
replicator.dynamics.3565.4 <- function(t, x, parms, a=0,b=0.13,c=0.18,d=0) {
  dx1 <- x[1]*(a*x[1]+b*x[2]-(x[1]*x[2]*(c+b)+a*x[1]*x[1]+d*x[2]*x[2]))
  dx2 <- x[2]*(c*x[1]+d*x[2]-(x[1]*x[2]*(c+b)+a*x[1]*x[1]+d*x[2]*x[2]))
  list(c(dx1, dx2))
}

replicator.dynamics.3565.10 <- function(t, x, parms, a=0,b=0.10,c=0.18,d=0) {
  dx1 <- x[1]*(a*x[1]+b*x[2]-(x[1]*x[2]*(c+b)+a*x[1]*x[1]+d*x[2]*x[2]))
  dx2 <- x[2]*(c*x[1]+d*x[2]-(x[1]*x[2]*(c+b)+a*x[1]*x[1]+d*x[2]*x[2]))
  list(c(dx1, dx2))
} 

replicator.dynamics.1919.4 <- function(t, x, parms, a=0,b=0.06,c=1,d=0) {
  dx1 <- x[1]*(a*x[1]+b*x[2]-(x[1]*x[2]*(c+b)+a*x[1]*x[1]+d*x[2]*x[2]))
  dx2 <- x[2]*(c*x[1]+d*x[2]-(x[1]*x[2]*(c+b)+a*x[1]*x[1]+d*x[2]*x[2]))
  list(c(dx1, dx2))
} 

replicator.dynamics.1919.10 <- function(t, x, parms, a=0,b=0.19,c=0.98,d=0) {
  dx1 <- x[1]*(a*x[1]+b*x[2]-(x[1]*x[2]*(c+b)+a*x[1]*x[1]+d*x[2]*x[2]))
  dx2 <- x[2]*(c*x[1]+d*x[2]-(x[1]*x[2]*(c+b)+a*x[1]*x[1]+d*x[2]*x[2]))
  list(c(dx1, dx2))
} 

replicator.dynamics.3128.4 <- function(t, x, parms, a=0,b=0.04,c=0.04,d=0) {
  dx1 <- x[1]*(a*x[1]+b*x[2]-(x[1]*x[2]*(c+b)+a*x[1]*x[1]+d*x[2]*x[2]))
  dx2 <- x[2]*(c*x[1]+d*x[2]-(x[1]*x[2]*(c+b)+a*x[1]*x[1]+d*x[2]*x[2]))
  list(c(dx1, dx2))
} 

replicator.dynamics.3128.10 <- function(t, x, parms, a=0,b=0.74,c=0.92,d=0) {
  dx1 <- x[1]*(a*x[1]+b*x[2]-(x[1]*x[2]*(c+b)+a*x[1]*x[1]+d*x[2]*x[2]))
  dx2 <- x[2]*(c*x[1]+d*x[2]-(x[1]*x[2]*(c+b)+a*x[1]*x[1]+d*x[2]*x[2]))
  list(c(dx1, dx2))
} 


x.ini.3565.4 <- as.numeric(unlist(data[data$TimeAfterRTx==0 &
                                         data$PatientID=="3565" &
                                         data$RTxDosage==4, c("MES1", "MES2")]))
x.ini.3565.10 <- as.numeric(unlist(data[data$TimeAfterRTx==0 &
                                          data$PatientID=="3565" &
                                          data$RTxDosage==10, c("MES1", "MES2")]))
x.ini.1919.4 <- as.numeric(unlist(data[data$TimeAfterRTx==0 &
                                         data$PatientID=="1919" &
                                         data$RTxDosage==4, c("MES1", "MES2")]))
x.ini.1919.10 <- as.numeric(unlist(data[data$TimeAfterRTx==0 &
                                          data$PatientID=="1919" &
                                          data$RTxDosage==10, c("MES1", "MES2")]))
x.ini.3128.4 <- as.numeric(unlist(data[data$TimeAfterRTx==0 &
                                         data$PatientID=="3128" &
                                         data$RTxDosage==4, c("MES1", "MES2")]))
x.ini.3128.10 <- as.numeric(unlist(data[data$TimeAfterRTx==0 &
                                          data$PatientID=="3128" &
                                          data$RTxDosage==10, c("MES1", "MES2")]))

names(x.ini.3565.4) <- c("MESlike1", "MESlike2")
names(x.ini.3565.10) <- c("MESlike1", "MESlike2")
names(x.ini.1919.4) <- c("MESlike1", "MESlike2")
names(x.ini.1919.10) <- c("MESlike1", "MESlike2")
names(x.ini.3128.4) <- c("MESlike1", "MESlike2")
names(x.ini.3128.10) <- c("MESlike1", "MESlike2")


out.3565.10 <- ode (times = seq(from = 0, to = 75, by = 0.01),
                    y = x.ini.1919.10, func = replicator.dynamics.1919.10, parms = NULL, method="ode23")
out.3565.4 <- ode (times = seq(from = 0, to = 75, by = 0.01),
                   y = x.ini.1919.4, func = replicator.dynamics.1919.4, parms = NULL, method="ode23")
out.1919.10 <- ode (times = seq(from = 0, to = 75, by = 0.01),
                    y = x.ini.1919.10, func = replicator.dynamics.1919.10, parms = NULL, method="ode23")
out.1919.4 <- ode (times = seq(from = 0, to = 75, by = 0.01),
                   y = x.ini.1919.4, func = replicator.dynamics.1919.4, parms = NULL, method="ode23")
out.3128.10 <- ode (times = seq(from = 0, to = 75, by = 0.01),
                    y = x.ini.3128.10, func = replicator.dynamics.3128.10, parms = NULL, method="ode23")
out.3128.4 <- ode (times = seq(from = 0, to = 75, by = 0.01),
                   y = x.ini.3128.4, func = replicator.dynamics.3128.4, parms = NULL, method="ode23")


data_3565_RTx_4 <- data_long[(data_long$RTxDosage == "0" | data_long$RTxDosage == "4") & data$PatientID =="3565",]
data_3565_RTx_10 <- data_long[(data_long$RTxDosage == "0" | data_long$RTxDosage == "10") & data$PatientID =="3565",]
data_1919_RTx_4 <- data_long[(data_long$RTxDosage == "0" | data_long$RTxDosage == "4") & data$PatientID =="1919",]
data_1919_RTx_10 <- data_long[(data_long$RTxDosage == "0" | data_long$RTxDosage == "10") & data$PatientID =="1919",]
data_3128_RTx_4 <- data_long[(data_long$RTxDosage == "0" | data_long$RTxDosage == "4") & data$PatientID =="3128",]
data_3128_RTx_10 <- data_long[(data_long$RTxDosage == "0" | data_long$RTxDosage == "10") & data$PatientID =="3128",]


plt1 <- ggplot(data = data_3565_RTx_4) +
  geom_point(data = data_3565_RTx_4,
             aes(y = Proportions, x = as.numeric(TimeAfterRTx), color = CellType),show.legend = FALSE) +
  geom_line(data = as.data.frame(out.3565.4), aes(x=time, y =MESlike2),
            color = "darkgreen", linetype = "dashed",show.legend = FALSE, position = "identity") +
  geom_line(data = as.data.frame(out.3565.4), aes(x=time, y =MESlike1),
            color = "darkred", linetype = "twodash",show.legend = FALSE, position = "identity") +
  ylim(0, 1) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 8),
                     plot.title = element_text(size = 10), 
                     axis.title.x = element_text(size = 9),
                     axis.title.y = element_text(size = 9)) +
  labs(title ="Patient 3 after 4-fold RTx", x = "Time", y = "Proportions", fill = "Cell type")+
  geom_abline(intercept = (13/31), slope = 0, color="red", 
              linetype="solid", size=0.1,alpha=0.2)+
  geom_abline(intercept = (18/31), slope = 0, color="darkgreen", 
              linetype="solid", size=0.1,alpha=0.2)

plt2 <- ggplot(data = data_3565_RTx_10) +
  geom_point(data = data_3565_RTx_10,
             aes(y = Proportions, x = as.numeric(TimeAfterRTx), color = CellType),show.legend = FALSE) +
  geom_line(data = as.data.frame(out.3565.10), aes(x=time, y =MESlike2),
            color = "darkgreen", linetype = "dashed",show.legend = FALSE, position = "identity") +
  geom_line(data = as.data.frame(out.3565.10), aes(x=time, y =MESlike1),
            color = "darkred", linetype = "twodash",show.legend = FALSE, position = "identity") +
  ylim(0, 1) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 8),
                     plot.title = element_text(size = 10), 
                     axis.title.x = element_text(size = 9),
                     axis.title.y = element_text(size = 9)) +
  labs(title ="Patient 3 after 10-fold RTx", x = "Time", y = "Proportions", fill = "Cell type")+
  geom_abline(intercept = (5/14), slope = 0, color="red", 
              linetype="solid", size=0.1,alpha=0.2)+
  geom_abline(intercept = (9/14), slope = 0, color="darkgreen", 
              linetype="solid", size=0.1,alpha=0.2)

plt3 <- ggplot(data = data_1919_RTx_4) +
  geom_point(data = data_1919_RTx_4,
             aes(y = Proportions, x = as.numeric(TimeAfterRTx), color = CellType),show.legend = FALSE) +
  geom_line(data = as.data.frame(out.1919.4), aes(x=time, y =MESlike2),
            color = "darkgreen", linetype = "dashed",show.legend = FALSE, position = "identity") +
  geom_line(data = as.data.frame(out.1919.4), aes(x=time, y =MESlike1),
            color = "darkred", linetype = "twodash",show.legend = FALSE, position = "identity") +
  ylim(0, 1) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 8),
                     plot.title = element_text(size = 10), 
                     axis.title.x = element_text(size = 9),
                     axis.title.y = element_text(size = 9)) +
  labs(title ="Patient 1 after 4-fold RTx", x = "Time", y = "Proportions", fill = "Cell type")+
  geom_abline(intercept = (3/53), slope = 0, color="red", 
              linetype="solid", size=0.1,alpha=0.2)+
  geom_abline(intercept = (50/53), slope = 0, color="darkgreen", 
              linetype="solid", size=0.1,alpha=0.2)

plt4 <- ggplot(data = data_1919_RTx_10) +
  geom_point(data = data_1919_RTx_10,
             aes(y = Proportions, x = as.numeric(TimeAfterRTx), color = CellType),show.legend = FALSE) +
  geom_line(data = as.data.frame(out.1919.10), aes(x=time, y =MESlike2),
            color = "darkgreen", linetype = "dashed",show.legend = FALSE, position = "identity") +
  geom_line(data = as.data.frame(out.1919.10), aes(x=time, y =MESlike1),
            color = "darkred", linetype = "twodash",show.legend = FALSE, position = "identity") +
  ylim(0, 1) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 8),
                     plot.title = element_text(size = 10), 
                     axis.title.x = element_text(size = 9),
                     axis.title.y = element_text(size = 9)) +
  labs(title ="Patient 1 after 10-fold RTx", x = "Time", y = "Proportions", fill = "Cell type")+
  geom_abline(intercept = (19/117), slope = 0, color="red", 
              linetype="solid", size=0.1,alpha=0.2)+
  geom_abline(intercept = (98/117), slope = 0, color="darkgreen", 
              linetype="solid", size=0.1,alpha=0.2)

plt5 <- ggplot(data = data_3128_RTx_4) +
  geom_point(data = data_3128_RTx_4,
             aes(y = Proportions, x = as.numeric(TimeAfterRTx), color = CellType),show.legend = FALSE) +
  geom_line(data = as.data.frame(out.3128.4), aes(x=time, y =MESlike2),
            color = "darkgreen", linetype = "dashed",show.legend = FALSE, position = "identity") +
  geom_line(data = as.data.frame(out.3128.4), aes(x=time, y =MESlike1),
            color = "darkred", linetype = "twodash",show.legend = FALSE, position = "identity") +
  ylim(0, 1) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 8),
                     plot.title = element_text(size = 10), 
                     axis.title.x = element_text(size = 9),
                     axis.title.y = element_text(size = 9)) +  labs(title ="Patient 2 after 4-fold RTx", x = "Time", y = "Proportions", fill = "Cell type")+
  geom_abline(intercept = (1/2), slope = 0, color="red", 
              linetype="solid", size=0.1,alpha=0.2)+
  geom_abline(intercept = (1/2), slope = 0, color="darkgreen", 
              linetype="solid", size=0.1,alpha=0.2)

plt6 <- ggplot(data = data_3128_RTx_10) +
  geom_point(data = data_3128_RTx_10,
             aes(y = Proportions, x = as.numeric(TimeAfterRTx), color = CellType),show.legend = FALSE) +
  geom_line(data = as.data.frame(out.3128.10), aes(x=time, y =MESlike2),
            color = "darkgreen", linetype = "dashed",show.legend = FALSE, position = "identity") +
  geom_line(data = as.data.frame(out.3128.10), aes(x=time, y =MESlike1),
            color = "darkred", linetype = "twodash",show.legend = FALSE, position = "identity") +
  ylim(0, 1) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 8),
                     plot.title = element_text(size = 10), 
                     axis.title.x = element_text(size = 9),
                     axis.title.y = element_text(size = 9)) +
  labs(title ="Patient 2 after 10-fold RTx", x = "Time", y = "Proportions", fill = "Cell type")+
  geom_abline(intercept = (36/83), slope = 0, color="red", 
              linetype="solid", size=0.1,alpha=0.2)+
  geom_abline(intercept = (47/83), slope = 0, color="darkgreen", 
              linetype="solid", size=0.1,alpha=0.2)

get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
  }

plt7 <- ggplot(data = data_3128_RTx_10) +
  geom_point(data = data_3128_RTx_10,
             aes(y = Proportions, x = as.numeric(TimeAfterRTx), color = CellType),show.legend = TRUE)+
  labs(title ="Patient 2 after 10-fold RTx", x = "Time", y = "Proportions", fill = "dupa cycki")


p_legend <- get_legend(plt7)

plot.data<-ggarrange(plt3, plt5, plt1, plt4, plt6, plt2,
                     ncol=3,nrow=2, common.legend = FALSE, legend="bottom")
plot.data

grid.arrange(arrangeGrob(plot.data), 
             p_legend, 
             nrow=2,heights=c(10, 1))
