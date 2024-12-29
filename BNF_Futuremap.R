
# windowsFonts(Font = windowsFont("Times New Roman"))
# par(family = 'Font')
setwd('4_Future response/')

library(terra)
library(tidyterra)
library(rnaturalearth)
library(ggplot2)
library(patchwork)
coast <- ne_coastline(scale = "medium", returnclass = "sf")

######## Global Change Factor Trend----
change_co2 <- rast('ChangeCO2-0.1.tif')
change_et <- rast('ChangeET-0.1.tif')
change_pre <-  rast('ChangePre-0.1.tif')
change_ndep <-  rast('ChangeNdep-0.1.tif')

g1 <- ggplot() +
  geom_spatraster(data = change_co2) +
  scale_fill_viridis_c(option = 'A',name = 'eCO2 ppm',na.value = 'transparent',begin = 1, end = 0,
                       values = c(0,1),limits = c(110,140)) +  
  geom_sf(data = coast)+
  coord_sf(crs = '+proj=robin')+ 
  theme_minimal()+
  theme(axis.text = element_text(size = 15),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 18),
        legend.position = "bottom",legend.box = "horizontal",
        legend.margin = margin(-20, 0, 20,20)
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.vjust = 0,title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
g1
# hist(change_co2)

g2 <- ggplot() +
  geom_spatraster(data = change_et) +
  scale_fill_viridis_c(option = 'A',name = 'eT ℃',na.value = 'transparent',begin = 1, end = 0,
                       values = c(0,1),limits = c(0,4)) +  
  geom_sf(data = coast)+
  coord_sf(crs = '+proj=robin')+ 
  theme_minimal()+
  theme(axis.text = element_text(size = 15),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 18),
        legend.position = "bottom",legend.box = "horizontal",
        legend.margin = margin(-20, 0, 20,20)
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.vjust = 0,title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
g2
# hist(change_et)

g3 <- ggplot() +
  geom_spatraster(data = change_pre) +
  scale_fill_viridis_c(option = 'H',name = 'ΔPre %',na.value = 'transparent',begin = 1, end = 0,
                       values = c(0,1),limits = c(-100,100)) +  
  geom_sf(data = coast)+
  coord_sf(crs = '+proj=robin')+ 
  theme_minimal()+
  theme(axis.text = element_text(size = 15),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 18),
        legend.position = "bottom",legend.box = "horizontal",
        legend.margin = margin(-20, 0, 20,20)
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.vjust = 0,title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
g3
# hist(change_pre)

g4 <- ggplot() +
  geom_spatraster(data = change_ndep) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                       midpoint = 0, name =expression(ΔNdep~kg~ha^{-1}*yr^{-1}), 
                       na.value = 'transparent', limits = c(-50,100)) +  
  geom_sf(data = coast)+
  coord_sf(crs = '+proj=robin')+ 
  theme_minimal()+
  theme(axis.text = element_text(size = 15),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 18),
        legend.position = "bottom",legend.box = "horizontal",
        legend.margin = margin(-20, 0, 20,20)
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.vjust = 0,title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
g4
# hist(change_ndep)

# png(filename = 'Global Change_4Factors.png',width = 2000,height=1200,units='px',bg='white',res=150,family ='Font')   #打开图形窗口
(g1|g2)/(g3|g4)
# dev.off()

## Future SNF and FNF changes global pattern----
SNFeCO <- rast('SNFeCO2-0.1.tif')
FNFeCO <- rast('FNFeCO2-0.1.tif')
SNFeT <-  rast('SNFeT-0.1.tif')
FNFeT <-  rast('FNFeT-0.1.tif')
SNFPre <-  rast('SNFPre-0.1.tif')
FNFPre <-  rast('FNFPre-0.1.tif')
SNFNdep <-  rast('SNFNdep-0.1.tif')
FNFNdep <-  rast('FNFNdep-0.1.tif')

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
        legend.position = "bottom",legend.box = "horizontal",
        legend.margin = margin(-20, 0, 20,20)
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
        legend.position = "bottom",legend.box = "horizontal",
        legend.margin = margin(-20, 0, 20,20)
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
        legend.position = "bottom",legend.box = "horizontal",
        legend.margin = margin(-20, 0, 20,20)
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
        legend.position = "bottom",legend.box = "horizontal",
        legend.margin = margin(-20, 0, 20,20)
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
        legend.position = "bottom",legend.box = "horizontal",
        legend.margin = margin(-20, 0, 20,20)
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
        legend.position = "bottom",legend.box = "horizontal",
        legend.margin = margin(-20, 0, 20,20)
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
        legend.position = "bottom",legend.box = "horizontal",
        legend.margin = margin(-20, 0, 20,20)
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
        legend.position = "bottom",legend.box = "horizontal",
        legend.margin = margin(-20, 0, 20,20)
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.vjust = 0,title.hjust = 1,
    barwidth = unit(8, "cm"),
    barheight = unit(0.5, "cm")
  ))
# p8
# hist(FNFNdep)

# png(filename = 'Global Change_BNF.png',width = 2000,height=2400,units='px',bg='white',res=150,family ='Font')   #打开图形窗口
(p1|p2)/(p3|p4)/(p5|p6)/(p7|p8)
# dev.off()

# Histogram of future SNF and FNF frequency changes
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

  