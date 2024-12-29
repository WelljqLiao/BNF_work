windowsFonts(Font = windowsFont("Times New Roman"))
par(family = 'Font')

##Fig.2 global map----
library(terra)
library(tidyterra)
library(ggplot2)
library(rnaturalearth)
coast <- ne_coastline(scale = "medium", returnclass = "sf")#海岸线

SNF_tif <- rast('3_Model result/SNF-0.1.tif')
FNF_tif <- rast('3_Model result/FNF-0.1.tif')

p1 <- ggplot() +
  geom_spatraster(data=SNF_tif)+
  scale_fill_viridis_c(option = 'D',name = 'SNF,kg/(ha*yr)',na.value = 'transparent',begin = 0, end = 1,
                       values = c(0,1),limits = c(0,40)) +  
  geom_sf(data = coast)+
  coord_sf(crs = '+proj=robin')+ 
  # labs(title = expression(a.))+
  theme_minimal()+
  theme(plot.title = element_text(size = 22), 
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "bottom",legend.box = "horizontal",
        legend.margin = margin(-20, 0, 20, 20)
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
p1
# hist(SNF_tif)

p2 <- ggplot() +
  geom_spatraster(data=FNF_tif)+
  scale_fill_viridis_c(option = 'D',name = 'FNF,kg/(ha*yr)',na.value = 'transparent',begin = 0, end = 1,
                       values = c(0,1),limits = c(0,10)) +  
  geom_sf(data = coast)+
  coord_sf(crs = '+proj=robin')+ 
  # labs(title = expression(d.))+
  theme_minimal()+
  theme(plot.title = element_text(size = 22), 
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "bottom",legend.box = "horizontal",
        legend.margin = margin(-20, 0, 20, 20)
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
p2
# hist(FNF_tif)


C_cols <- c('#fce624','#a6d933','#57c463','#26a682','#23868c','#2f678d','#3E4884','#442677','#450056') 
C_cols <- rev(C_cols)

colors <- colorRampPalette(C_cols)(40)
hist(SNF_tif,breaks = 50,freq = F,main = NA,col = colors,xlim = c(0,40),xlab = 'SNF', ylab = NA)

colors <- colorRampPalette(C_cols)(35)
hist(FNF_tif,breaks = 50,freq = F,main = NA,col = colors,xlim = c(0,10),xlab = 'FNF', ylab = NA)

##Fig.S global uncertianty----
SNFun_tif <- rast("3_Model result/SNFcv-0.1.tif")
FNFun_tif <- rast("3_Model result/FNFcv-0.1.tif")

summary(SNFun_tif)
SNFun_tif[SNFun_tif <= 0 ] = 0.14
SNFun_tif[SNFun_tif >= 1 ] = 0.14

summary(FNFun_tif)
FNFun_tif[FNFun_tif <= 0 ] = 0.18
FNFun_tif[FNFun_tif >= 1 ] = 0.18

p5 <- ggplot() +
  geom_spatraster(data=SNFun_tif)+
  scale_fill_viridis_c(option = 'A',name = 'CV',na.value = 'transparent',begin = 0.5, end = 1,
                       values = c(0,1),direction = -1,limits = c(0,0.6)) +  
  geom_sf(data = coast)+
  coord_sf(ylim = c(-56, 90),expand = FALSE)+ 
  labs(title = expression(a.~SNF~uncertainty))+
  theme_minimal()+
  theme(plot.title = element_text(size = 20), 
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15))
p5

p6 <- ggplot() +
  geom_spatraster(data = FNFun_tif)+
  scale_fill_viridis_c(option = 'A', name = 'CV',na.value = 'transparent',begin = 0.5, end = 1,
                       values = c(0,1),direction = -1,limits = c(0,1)) +
  geom_sf(data = coast)+
  coord_sf(ylim = c(-56, 90),expand = FALSE)+
  labs(title = expression(b.~FNF~uncertainty))+
  theme_minimal()+
  theme(plot.title = element_text(size = 20), 
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15))
p6


colors <- colorRampPalette(c("#FDF4B6","#B52F74"))(70)
hist(SNFun_tif,breaks = 150,freq = F,main = NA,col = colors,
     xlim = c(0.05,0.4),xlab = 'SNF uncertainty (CV)', ylab = 'Count',
     cex.axis = 1.5, cex.lab = 1.5)

hist(FNFun_tif,breaks = 100,freq = F,main = NA,col = colors,
     xlim = c(0,1),xlab = 'FNF uncertainty (CV)', ylab = 'Count',
     cex.axis = 1.5, cex.lab = 1.5)


## Fig.2 ecosystem BNF rates----
library(ggplot2) 
library(ggsignif) 
library(gghalves) 
library(dplyr)
C_cols <- c('#fce624','#a6d933','#57c463','#26a682','#23868c','#2f678d','#3E4884','#442677','#450056') 
C_cols <- rev(C_cols)
C_pal <- colorRampPalette(C_cols) #颜色分层函数
df <- readRDS("3_Model result/SNF_df.rds")
df$group <- factor(df$group,levels = c("ENF","EBF","DNF","DBF","MIX","SHB","WAS","SAV","GRS"))
label <- expression("SNF, (kg ha"^{-1} * "yr"^{-1} * ")")


ggplot(df,aes(group,values,fill=group))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.8,color=NA)+
  geom_jitter(data = sample_frac(df,0.1),aes(fill=group),shape=21,size=1,width=0.15,alpha = 0.08)+
  geom_boxplot(width=0.4,size=1.2,outlier.color =NA)+
  geom_hline(yintercept = mean(df$values), linetype = 2, color = "black",linewidth=1)+
  scale_y_continuous(limits = c(0,65),breaks = c(0,20,40,60))+ # sns scale_y
  # scale_y_continuous(limits = c(0,10),breaks = c(0,3,6,9))+  # fnf scale_y
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1),
        axis.text.x = element_text(color = "black", size = 22),
        axis.text.y = element_text(color = "black",size = 22),
        axis.title.x = element_text(size = 24),  # 设置 y 轴标题字体大小
        legend.position = "none",
        axis.ticks = element_line(color="black",linewidth = 1))+
  scale_fill_manual(values = C_pal(9))+
  coord_flip() +
  labs(x = NULL,y=label)


library(agricolae)
variance <- aov(values ~ group, data=df)
MC <- LSD.test(variance,"group", p.adj="none")
MC


##Fige.S BNF和FNFratio----
BNF_tif <- SNF_tif + FNF_tif
FNFratio_tif <- FNF_tif/BNF_tif

s1 <- ggplot() +
  geom_spatraster(data=BNF_tif)+
  scale_fill_viridis_c(option = 'D',name = 'BNF,kg/(ha*yr)',na.value = 'transparent',begin = 0, end = 1,
                       values = c(0,1),limits = c(0,40)) +  
  geom_sf(data = coast)+
  coord_sf(crs = '+proj=robin')+ 
  # labs(title = expression(a.))+
  theme_minimal()+
  theme(plot.title = element_text(size = 22), 
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "bottom",legend.box = "horizontal",
        legend.margin = margin(-20, 0, 20, 20)
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
s1

s2 <- ggplot() +
  geom_spatraster(data=FNFratio_tif)+
  scale_fill_viridis_c(option = 'D',name = 'FNF:BNF',na.value = 'transparent',begin = 0, end = 1,
                       values = c(0,1),limits = c(0,1)) +  
  geom_sf(data = coast)+
  coord_sf(crs = '+proj=robin')+ 
  # labs(title = expression(d.))+
  theme_minimal()+
  theme(plot.title = element_text(size = 22), 
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "bottom",legend.box = "horizontal",
        legend.margin = margin(-20, 0, 20, 20)
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
s2

colors <- colorRampPalette(C_cols)(100)
hist(BNF_tif,breaks = 100,freq = F,main = NA,col = colors,
     xlim = c(0,50),xlab = 'BNF,kg/(ha*yr)', ylab = 'Count',
     cex.axis = 1.5, cex.lab = 1.5)

hist(FNFratio_tif,breaks = 100,freq = F,main = NA,col = colors,
     xlim = c(0,1),xlab = 'FNFratio (%)', ylab = 'Count',
     cex.axis = 1.5, cex.lab = 1.5)

