library(terra)
library(tidyterra)
library(rnaturalearth)
library(ggplot2)
library(patchwork)
coast <- ne_coastline(scale = "medium", returnclass = "sf")

windowsFonts(Font = windowsFont("Times New Roman"))
par(family = 'Font')
setwd('../2_Interim')

## Fig.C-J
## Future SNF and FNF changes global pattern----
SNFeCO <- rast('SNF_eCO2.tif')
FNFeCO <- rast('FNF_eCO2.tif')
SNFeT <-  rast('SNF_eT.tif')
FNFeT <-  rast('FNF_eT.tif')
SNFPre <-  rast('SNF_Pre.tif')
FNFPre <-  rast('FNF_Pre.tif')
SNFNdep <-  rast('SNF_Ndep.tif')
FNFNdep <-  rast('FNF_Ndep.tif')

Totalchange <- SNFeCO + FNFeCO + SNFeT + FNFeT + SNFPre + FNFPre + SNFNdep + FNFNdep

p1 <- ggplot() +
  geom_spatraster(data=SNFeCO)+
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                       midpoint = 0,  name =expression(kg~ha^{-1}*yr^{-1}), 
                       na.value = 'transparent', limits = c(0,2)) +  
  geom_sf(data = coast)+
  coord_sf(crs = '+proj=robin')+ 
  theme_minimal()+
  theme(axis.text = element_text(size = 15),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 18),
        legend.position = "bottom",legend.box = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.vjust = 0,title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
p1
# hist(SNFeCO)

p2 <- ggplot() +
  geom_spatraster(data = FNFeCO) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                       midpoint = 0, name =expression(kg~ha^{-1}*yr^{-1}), 
                       na.value = 'transparent', limits = c(0, 4)) +  
  geom_sf(data = coast)+
  coord_sf(crs = '+proj=robin')+ 
  theme_minimal()+
  theme(axis.text = element_text(size = 15),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 18),
        legend.position = "bottom",legend.box = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.vjust = 0,title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
# p2
# hist(FNFeCO)

p3 <- ggplot() +
  geom_spatraster(data = SNFeT) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                       midpoint = 0, name =expression(kg~ha^{-1}*yr^{-1}), 
                       na.value = 'transparent', limits = c(0, 25)) +  
  geom_sf(data = coast)+
  coord_sf(crs = '+proj=robin')+ 
  theme_minimal()+
  theme(axis.text = element_text(size = 15),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 18),
        legend.position = "bottom",legend.box = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.vjust = 0,title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
# p3
# hist(SNFeT)

p4 <- ggplot() +
  geom_spatraster(data = FNFeT) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                       midpoint = 0,  name =expression(kg~ha^{-1}*yr^{-1}), 
                       na.value = 'transparent', limits = c(0, 25)) +  
  geom_sf(data = coast)+
  coord_sf(crs = '+proj=robin')+ 
  theme_minimal()+
  theme(axis.text = element_text(size = 15),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 18),
        legend.position = "bottom",legend.box = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.vjust = 0,title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
# p4
# hist(FNFeT)

p5 <- ggplot() +
  geom_spatraster(data = SNFPre) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                       midpoint = 0, name =expression(kg~ha^{-1}*yr^{-1}), 
                       na.value = 'transparent', limits = c(-10, 10)) +  
  geom_sf(data = coast)+
  coord_sf(crs = '+proj=robin')+ 
  theme_minimal()+
  theme(axis.text = element_text(size = 15),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 18),
        legend.position = "bottom",legend.box = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.vjust = 0,title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
# p5
# hist(SNFPre)

p6 <- ggplot() +
  geom_spatraster(data = FNFPre) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                       midpoint = 0, name =expression(kg~ha^{-1}*yr^{-1}), 
                       na.value = 'transparent', limits = c(-2, 2)) +  
  geom_sf(data = coast)+
  coord_sf(crs = '+proj=robin')+ 
  theme_minimal()+
  theme(axis.text = element_text(size = 15),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 18),
        legend.position = "bottom",legend.box = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.vjust = 0,title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
# p6
# hist(FNFPre)

p7 <- ggplot() +
  geom_spatraster(data = SNFNdep) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                       midpoint = 0, name =expression(kg~ha^{-1}*yr^{-1}), 
                       na.value = 'transparent', limits = c(-10, 10)) +  
  geom_sf(data = coast)+
  coord_sf(crs = '+proj=robin')+ 
  theme_minimal()+
  theme(axis.text = element_text(size = 15),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 18),
        legend.position = "bottom",legend.box = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.vjust = 0,title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
# p7
# hist(SNFNdep)

p8 <- ggplot() +
  geom_spatraster(data = FNFNdep) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                       midpoint = 0, name =expression(kg~ha^{-1}*yr^{-1}), 
                       na.value = 'transparent', limits = c(-10, 10)) +  
  geom_sf(data = coast)+
  coord_sf(crs = '+proj=robin')+ 
  theme_minimal()+
  theme(axis.text = element_text(size = 15),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 18),
        legend.position = "bottom",legend.box = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.vjust = 0,title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
# p8
# hist(FNFNdep)

# png(filename = 'Global Change_BNF.png',width = 7200,height=2400,units='px',bg='white',res=300,family ='Font')   #打开图形窗口
(p1|p3|p5|p7)/(p2|p4|p6|p8)
# dev.off()

## Fig.Support
# Histogram of future SNF and FNF frequency changes
# png(filename = 'Global Change_hist.png',width = 2400,height=1200,units='px',bg='white',res=150,family ='Font')   #打开图形窗口
par(mfrow = c(2,4),cex = 1, cex.lab = 1.2, cex.axis = 1.1, cex.main = 1,
    mar = c(4, 4, 1, 1), oma = c(1, 1, 1, 1))
hist(SNFeCO,breaks = 50,freq = F,main = NA,col = "#336990",xlim = c(-0.05,1),xlab = 'ΔSNF-eCO2', ylab = 'Count')
mtext("a", side = 3, line = -1, adj = 0.01, cex = 2)
hist(SNFeT,breaks = 50,freq = F,main = NA,col = "#336990",xlim = c(-3,25),xlab = 'ΔSNF-eT', ylab ='Count')
mtext("b", side = 3, line = -1, adj = 0.01, cex = 2)
hist(SNFPre,breaks = 150,freq = F,main = NA,col ="#336990",xlim = c(-10,10),xlab = 'ΔSNF-Pre', ylab ='Count')
mtext("c", side = 3, line = -1, adj = 0.01, cex = 2)
hist(SNFNdep,breaks = 150,freq = F,main = NA,col = "#336990",xlim = c(-10,10),xlab = 'ΔSNF-Ndep', ylab ='Count')
mtext("d", side = 3, line = -1, adj = 0.01, cex = 2)
hist(FNFeCO,breaks = 50,freq = F,main = NA,col = "#209D86",xlim = c(0,2.5),xlab = 'ΔFNF-eCO2', ylab = 'Count')
mtext("e", side = 3, line = -1, adj = 0.01, cex = 2)
hist(FNFeT,breaks = 50,freq = F,main = NA,col = "#209D86",xlim = c(0,20),xlab = 'ΔFNF-eT', ylab ='Count')
mtext("f", side = 3, line = -1, adj = 0.01, cex = 2)
hist(FNFPre,breaks = 50,freq = F,main = NA,col ="#209D86",xlim = c(-5,10),xlab = 'ΔFNF-Pre', ylab ='Count')
mtext("g", side = 3, line = -1, adj = 0.01, cex = 2)
hist(FNFNdep,breaks = 80,freq = F,main = NA,col = "#209D86",xlim = c(-10,5),xlab = 'ΔFNF-Ndep', ylab ='Count')
mtext("h", side = 3, line = -1, adj = 0.01, cex = 2)
# dev.off()
  
## BNF(sum of all) under all global change 
pall <- ggplot() +
  geom_spatraster(data = Totalchange) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                       midpoint = 0, name =expression(kg~ha^{-1}*yr^{-1}), 
                       na.value = 'transparent', limits = c(-10, 30)) +  
  geom_sf(data = coast)+
  coord_sf(crs = '+proj=robin')+ 
  theme_minimal()+
  theme(axis.text = element_text(size = 15),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 18),
        legend.position = "bottom",legend.box = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.vjust = 0,title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
pall
