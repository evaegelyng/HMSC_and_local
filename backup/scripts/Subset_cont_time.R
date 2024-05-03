# load packages
library("ggplot2")
library("reshape2")

#reading files
pc_bs<-read.table("data/merged_metadata_230427.txt", sep="\t", header=T)

#Time
variables_ite<-c("Temperature", "Latitude", "Longitude", "Salinity")

htdj <- melt(pc_bs, id=c("season","habitat","sal_group","cluster","Time"), measure=variables_ite)
htdj$sal_group<-as.character(htdj$sal_group)
htdj$season <- factor(htdj$season, levels = c("spring","autumn"))
htdj$habitat <- factor(htdj$habitat, levels = c("rocks","sand","eelgrass"))
htdj$Time<-as.Date(htdj$Time, format= "%d-%m-%y")

ggplot(htdj, aes(Time, value)) +
geom_point(data=htdj, aes(Time, value), size=0.9, alpha=0.7) +
geom_smooth(data=htdj, aes(group=season), method="lm", size=0.5, se=F, color="red") +
facet_grid(variable~season, scale="free") +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
scale_x_date(date_breaks = "15 days", date_labels = "%d \n %m") +
theme(strip.text = element_text(size=8),
legend.title=element_text(size=4),legend.text=element_text(size=4), axis.title=element_text(size=8), plot.title = element_text(size=8), axis.text.x = element_text(size=8), axis.text.y = element_text(size=8), panel.grid.major = element_line(colour = "gray80", size=0.3, linetype=2))
ggsave("results/Both/Subset_cont_Time_EES.png")