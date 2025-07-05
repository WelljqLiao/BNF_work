windowsFonts(Font = windowsFont("Times New Roman"))
par(family = 'Font')
library(terra)
library(tidyterra)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(readxl)

##Fig.1a site map----
setwd('../0_Data')
snf_data <- read_excel("BNFdata_use.xlsx", sheet = "Symbiotic")
fnf_data <- read_excel("BNFdata_use.xlsx", sheet = "Free-living")
snf_data$Type <- "SNF"
fnf_data$Type <- "FNF"
all_data <- rbind(snf_data[, c("lat", "lon", "Type")], 
                  fnf_data[, c("lat", "lon", "Type")])

world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_sf(data = world, fill = "grey90", color = "white", size = 0.2) +
  geom_point(data = all_data, 
             aes(x = lon, y = lat, color = Type),
             size = 3, alpha = 0.7) +
  scale_color_manual(values = c("SNF" = "#209D86", "FNF" = "#336990"),
                     labels = c("SNF", "FNF")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12),
    legend.position = "none",
    legend.title = element_blank()
  ) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  labs(x = "", y = "")
# ggsave("study_sites_distribution.png", width = 10, height = 6, dpi = 300)


##Fig.2a and d global map----
setwd('../2_Interim')
coast <- ne_coastline(scale = "medium", returnclass = "sf")
SNF_tif <- rast('SNF_predict_0520.tif')
FNF_tif <- rast('FNF_predict_0520.tif')

# Fig.2A
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
        legend.position = "bottom",legend.box = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
# png(filename = 'Global pattern_SNF.png',width = 2800,height=1400,units='px',bg='white',res=300,family ='Font')   #打开图形窗口
p1
# dev.off()

# Fig.2D
p2 <- ggplot() +
  geom_spatraster(data=FNF_tif)+
  scale_fill_viridis_c(option = 'D',name = 'FNF,kg/(ha*yr)',na.value = 'transparent',begin = 0, end = 1,
                       values = c(0,1),limits = c(0,6)) +  
  geom_sf(data = coast)+
  coord_sf(crs = '+proj=robin')+ 
  # labs(title = expression(d.))+
  theme_minimal()+
  theme(plot.title = element_text(size = 22), 
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "bottom",legend.box = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
# png(filename = 'Global pattern_FNF.png',width = 2800,height=1400,units='px',bg='white',res=300,family ='Font')   #打开图形窗口
p2
# dev.off()


C_cols <- c('#fce624','#a6d933','#57c463','#26a682','#23868c','#2f678d','#3E4884','#442677','#450056') 
C_cols <- rev(C_cols)
colors <- colorRampPalette(C_cols)(75)
# png(filename = 'SNF_hist.png',width = 800,height=800,units='px',bg='NA',res = 300,family ='Font')   #打开图形窗口
hist(SNF_tif,breaks = 100,freq = F,main = NA,col = colors,xlim = c(0,40),xlab = '', ylab = NA)
# dev.off()

colors <- colorRampPalette(C_cols)(65)
# png(filename = 'FNF_hist.png',width = 800,height=800,units='px',bg='NA',res = 300,family ='Font')   #打开图形窗口
hist(FNF_tif,breaks = 100,freq = F,main = NA,col = colors,xlim = c(0,6),xlab = '', ylab = NA)
# dev.off()

##Fig.SX global uncertianty
SNFun_tif <- rast("SNF_uncer_0520.tif")
FNFun_tif <- rast("FNF_uncer_0520.tif")
p5 <- ggplot() +
  geom_spatraster(data=SNFun_tif)+
  scale_fill_viridis_c(option = 'A',name = 'CV',na.value = 'transparent',begin = 0.5, end = 1,
                       values = c(0,1),direction = -1,limits = c(0,0.27)) +  
  geom_sf(data = coast)+
  coord_sf(ylim = c(-56, 90),expand = FALSE)+ 
  labs(title = expression(a.~SNF~uncertainty))+
  theme_minimal()+
  theme(plot.title = element_text(size = 20), 
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15))
p5
# png(filename = 'Global pattern_SNF-CV.png',width = 2800,height=1400,units='px',bg='white',res=300,family ='Font')   #打开图形窗口
# p5
# dev.off()

p6 <- ggplot() +
  geom_spatraster(data = FNFun_tif)+
  scale_fill_viridis_c(option = 'A', name = 'CV',na.value = 'transparent',begin = 0.5, end = 1,
                       values = c(0,1),direction = -1,limits = c(0,0.46)) +
  geom_sf(data = coast)+
  coord_sf(ylim = c(-56, 90),expand = FALSE)+
  labs(title = expression(b.~FNF~uncertainty))+
  theme_minimal()+
  theme(plot.title = element_text(size = 20), 
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15))
p6
# png(filename = 'Global pattern_FNF-CV.png',width = 2800,height=1400,units='px',bg='white',res=300,family ='Font')   #打开图形窗口
# p6
# dev.off()

colors <- colorRampPalette(c("#FDF4B6","#B52F74"))(90)
# png(filename = 'SNF_cv.png',width = 1800,height= 1200,units='px',bg='NA',res = 300,family ='Font')   #打开图形窗口
hist(SNFun_tif,breaks = 150,freq = F,main = NA,col = colors,
     xlim = c(0.0,0.25),xlab = 'SNF uncertainty (CV)', ylab = 'Count',
     cex.axis = 1.5, cex.lab = 1.5)
# dev.off()

colors <- colorRampPalette(c("#FDF4B6","#B52F74"))(50)
# png(filename = 'FNF_cv.png',width = 1800,height= 1200,units='px',bg='NA',res = 300,family ='Font')   #打开图形窗口
hist(FNFun_tif,breaks = 100,freq = F,main = NA,col = colors,
     xlim = c(0,0.46),xlab = 'FNF uncertainty (CV)', ylab = 'Count',
     cex.axis = 1.5, cex.lab = 1.5)
# dev.off()

#### ecosystem summary
SNF <- rast('SNF_predict_0520.tif')
FNF <- rast('FNF_predict_0520.tif')
BNF <- SNF+FNF
FNFratio <- FNF/BNF

SNF_mean = global(SNF, stat = 'mean',na.rm=T)
print(SNF_mean)
FNF_mean = global(FNF, stat = 'mean',na.rm=T)
print(FNF_mean)

Landcover <- rast('Landcover_WGS84.tif')
plot(Landcover)
Landcover[Landcover == 0] = 0 #water
Landcover[Landcover >= 12 & Landcover <= 14] = 0 #Excluding human areas such as farmland
Landcover[Landcover == 15] = 0 #Excluding snow
Landcover[Landcover == 16] = 0 #Excluding bare
Sd_raster <- rast(nrow = 360, ncol = 720,
                    xmin = -180, xmax = 180,
                    ymin = -90, ymax = 90)
Landcover <- crop(Landcover,Sd_raster)  
Land <- resample(Landcover,Sd_raster,method = 'near') 

ENF <- Land == 1
EBF <- Land == 2
DNF <- Land == 3
DBF <- Land == 4
MF <- Land == 5
Shrub <- Land == 7
Sav_w <- Land == 8
Sav <- Land == 9
Grass <- Land == 10
Cropland <- Land >= 12 & Land <=14
process_region <- Land >= 0 & Land <=15

# choose SNF, FNF, or BNF
#SNF
ENF_SNF <- mask(SNF,ENF,maskvalue=0)
EBF_SNF <- mask(SNF,EBF,maskvalue=0)
DNF_SNF <- mask(SNF,DNF,maskvalue=0)
DBF_SNF <- mask(SNF,DBF,maskvalue=0)
MF_SNF <- mask(SNF,MF,maskvalue=0)
Shrub_SNF <- mask(SNF,Shrub,maskvalue=0)
Sav_w_SNF <- mask(SNF,Sav_w,maskvalue=0)
Sav_SNF <- mask(SNF,Sav,maskvalue=0)
Grass_SNF <- mask(SNF,Grass,maskvalue=0)

#FNF
ENF_FNF <- mask(FNF,ENF,maskvalue=0)
EBF_FNF <- mask(FNF,EBF,maskvalue=0)
DNF_FNF <- mask(FNF,DNF,maskvalue=0)
DBF_FNF <- mask(FNF,DBF,maskvalue=0)
MF_FNF <- mask(FNF,MF,maskvalue=0)
Shrub_FNF <- mask(FNF,Shrub,maskvalue=0)
Sav_w_FNF <- mask(FNF,Sav_w,maskvalue=0)
Sav_FNF <- mask(FNF,Sav,maskvalue=0)
Grass_FNF <- mask(FNF,Grass,maskvalue=0)

#BNF
ENF_BNF <- mask(BNF,ENF,maskvalue=0)
EBF_BNF <- mask(BNF,EBF,maskvalue=0)
DNF_BNF <- mask(BNF,DNF,maskvalue=0)
DBF_BNF <- mask(BNF,DBF,maskvalue=0)
MF_BNF <- mask(BNF,MF,maskvalue=0)
Shrub_BNF <- mask(BNF,Shrub,maskvalue=0)
Sav_w_BNF <- mask(BNF,Sav_w,maskvalue=0)
Sav_BNF <- mask(BNF,Sav,maskvalue=0)
Grass_BNF <- mask(BNF,Grass,maskvalue=0)

data_df_A <- as.data.frame(values(ENF_FNF))
data_df_B <- as.data.frame(values(EBF_FNF))
data_df_C <- as.data.frame(values(DNF_FNF))
data_df_D <- as.data.frame(values(DBF_FNF))
data_df_E <- as.data.frame(values(MF_FNF))
data_df_F <- as.data.frame(values(Shrub_FNF))
data_df_G <- as.data.frame(values(Sav_w_FNF))
data_df_H <- as.data.frame(values(Sav_FNF))
data_df_I <- as.data.frame(values(Grass_FNF))

data_df_A <- na.omit(data_df_A)
data_df_B <- na.omit(data_df_B)
data_df_C <- na.omit(data_df_C)
data_df_D <- na.omit(data_df_D)
data_df_E <- na.omit(data_df_E)
data_df_F <- na.omit(data_df_F)
data_df_G <- na.omit(data_df_G)
data_df_H <- na.omit(data_df_H)
data_df_I <- na.omit(data_df_I)

data_df_A$group <- "ENF"
data_df_B$group <- "EBF"
data_df_C$group <- "DNF"
data_df_D$group <- "DBF"
data_df_E$group <- "MIX"
data_df_F$group <- "SHB"
data_df_G$group <- "WAS"
data_df_H$group <- "SAV"
data_df_I$group <- "GRS"

colnames(data_df_A) <- c("values", "group")
colnames(data_df_B) <- c("values", "group")
colnames(data_df_C) <- c("values", "group")
colnames(data_df_D) <- c("values", "group")
colnames(data_df_E) <- c("values", "group")
colnames(data_df_F) <- c("values", "group")
colnames(data_df_G) <- c("values", "group")
colnames(data_df_H) <- c("values", "group")
colnames(data_df_I) <- c("values", "group")

# 合并数据框
data_df <- rbind(data_df_A, data_df_B, data_df_C,
                 data_df_D, data_df_E, data_df_F,
                 data_df_G, data_df_H, data_df_I)

# saveRDS(data_df, "SNF_df.rds")


## Fig.2b and e ecosystem BNF rates----
library(ggplot2) 
library(ggsignif) 
library(gghalves) 
library(dplyr)
C_cols <- c('#fce624','#a6d933','#57c463','#26a682','#23868c','#2f678d','#3E4884','#442677','#450056') 
C_cols <- rev(C_cols)
C_pal <- colorRampPalette(C_cols) 

# read "SNF_df.rds" or "FNF_df.rds"
df <- readRDS("FNF_df.rds")
df$group <- factor(df$group,levels = c("ENF","EBF","DNF","DBF","MIX","SHB","WAS","SAV","GRS"))
label <- expression("FNF, (kg ha"^{-1} * "yr"^{-1} * ")")

# Fig 2b & 2e
# png(filename = 'ecoregionFNF.png',width = 800,height= 1600,units='px',bg= NA,res=200,family ='Font')
ggplot(df,aes(group,values,fill=group))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.8,color=NA)+
  geom_jitter(data = sample_frac(df,0.1),aes(fill=group),shape=21,size=1,width=0.15,alpha = 0.08)+
  geom_boxplot(width=0.4,size=1.2,outlier.color =NA)+
  geom_hline(yintercept = mean(df$values), linetype = 2, color = "black",linewidth=1)+
  # scale_y_continuous(limits = c(0,42),breaks = c(0,20,40))+ # snf scale_y
  scale_y_continuous(limits = c(0,7),breaks = c(0,3,6))+  # fnf scale_y
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1),
        axis.text.x = element_text(color = "black", size = 22),
        axis.text.y = element_text(color = "black",size = 22),
        axis.title.x = element_text(size = 24), 
        legend.position = "none",
        axis.ticks = element_line(color="black",linewidth = 1))+
  scale_fill_manual(values = C_pal(9))+
  coord_flip() +
  labs(x = NULL,y=label)
# dev.off()

# annova
library(agricolae)
variance <- aov(values ~ group, data=df)
MC <- LSD.test(variance,"group", p.adj="none")
MC

##Fige.Support BNF and FNFratio plot----
s1 <- ggplot() +
  geom_spatraster(data=BNF)+
  scale_fill_viridis_c(option = 'D',name = 'BNF,kg/(ha*yr)',na.value = 'transparent',begin = 0, end = 1,
                       values = c(0,1),limits = c(0,40)) +  
  geom_sf(data = coast)+
  coord_sf(crs = '+proj=robin')+ 
  # labs(title = expression(a.))+
  theme_minimal()+
  theme(plot.title = element_text(size = 22), 
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "bottom",legend.box = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
s1
# png(filename = 'Global pattern_BNF.png',width = 2800,height=1400,units='px',bg='white',res=300,family ='Font')   #打开图形窗口
# s1
# dev.off()

s2 <- ggplot() +
  geom_spatraster(data=FNFratio)+
  scale_fill_viridis_c(option = 'D',name = 'FNF:BNF',na.value = 'transparent',begin = 0, end = 1,
                       values = c(0,1),limits = c(0,0.6)) +  
  geom_sf(data = coast)+
  coord_sf(crs = '+proj=robin')+ 
  # labs(title = expression(d.))+
  theme_minimal()+
  theme(plot.title = element_text(size = 22), 
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "bottom",legend.box = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
s2
# png(filename = 'Global pattern_FBratio.png',width = 2800,height=1400,units='px',bg='white',res=300,family ='Font')   #打开图形窗口
# s2
# dev.off()

C_cols <- c('#fce624','#a6d933','#57c463','#26a682','#23868c','#2f678d','#3E4884','#442677','#450056') 
C_cols <- rev(C_cols)
colors <- colorRampPalette(C_cols)(75)
# png(filename = 'BNF_hist.png',width = 1800,height= 1200,units='px',bg='NA',res = 300,family ='Font')   #打开图形窗口
hist(BNF,breaks = 100,freq = F,main = NA,col = colors,
     xlim = c(0,40),xlab = 'BNF,kg/(ha*yr)', ylab = 'Count',
     cex.axis = 1.5, cex.lab = 1.5)
# dev.off()

colors <- colorRampPalette(C_cols)(60)
# png(filename = 'FBratio_hist.png',width = 1800,height= 1200,units='px',bg='NA',res = 300,family ='Font')   #打开图形窗口
hist(FNFratio,breaks = 90,freq = F,main = NA,col = colors,
     xlim = c(0,0.6),xlab = 'FNFratio (%)', ylab = 'Count',
     cex.axis = 1.5, cex.lab = 1.5)
# dev.off()
