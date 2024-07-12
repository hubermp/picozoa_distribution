#title: Global Distribution, Diversity, and Ecological Niche of Picozoa, a widespread and enigmatic marine protist lineage

############################################################################################################## 

##### Script: Picozoa_distribution_and_diversity.R
##### Created: 11/02/2024
##### Author: Daniele De Angelis daniele.deangelis_at_uniroma1.it
##### Objective: to run CANOMI analysis and compute the niche overlap of Picozoas taxa in the sunlit Ocean
##### NOTE: for data processing please refer to the manuscript: Global Distribution, Diversity, and Ecological Niche of Picozoa, a widespread and enigmatic marine protist lineage. This script is designed such that each code chunk can be executed independently. T

############################################################################################################## 


### LOAD LIBRARIES ####

library(dplyr)
library(ggplot2)
library(adehabitatHS)
library(reshape2)
library(raster)
library(adehabitatHR)
library(ENMeval)
library(corrplot)
library(xlsx)



##########################################################################################
####---- LOAD RASTER VARIABLES AND CREATE RASTER STACK SEPARATED BY DEPTH LEVEL ----######

MLD<-raster("cological Niche/env_variables/M02.tif")
names(MLD)<-c("mld")

# Open environmental variables (surface layer)
(names <- list.files(path = pathpath, 
                     pattern = ".tif", full.names = TRUE))
#names <- names[-grep("aux.xml", names)]

index <- grep("epi",ignore.case = T,names) #select only mean variables
stacked_epi <- stack(names[index])
names(stacked_epi) <- gsub(".tif","", gsub(paste0(pathpath,"/"),"", names[index]))
names(stacked_epi)
dim(stacked_epi)

# 

Nstar_epi<-(stacked_epi[["n_epi"]]-16*stacked_epi[["p_epi"]])
names(Nstar_epi) <- "Nstar_epi"
stacked_epi<-stack(stacked_epi, Nstar_epi)


stacked_epi
names(stacked_epi) <- gsub("_epi", "", names(stacked_epi))
stacked_epi <- stacked_epi[[which(!names(stacked_epi) %in% c("si", "O2sat", "Nstar", "chl", "A", "I", "ugo", "vgo"))]]
names(stacked_epi) 
#plot(stacked_epi)



##########################################################################################
####---- LOAD PICOZOA DATA WITH PTU READ COUNTS ----######
dd3 <- read.csv("NEW_PICOZOA_DATA_DDA_AUG23.csv")[,-1]

info <- read.xlsx("pOTU_tax_information_Jan24.xlsx",1)
info <- info[,c("Final_Name","Lat_distribution")]

# pOTU049 – now is considered as a NP
# pOTU054 – now considered as P
# pOTU020 – now considered as P

# OTUs classes
cosmopolitan <- info[info$Lat_distribution == "C","Final_Name"]
polar <- info[info$Lat_distribution == "P","Final_Name"]
nonpolar <- info[info$Lat_distribution == "NP","Final_Name"]

# extract variables
DF_bin <- dd3[!is.na(dd3$Lon),] #0 sampling stations have no coords
DF_bin <- DF_bin[!is.na(DF_bin$depth),]  #129 sampling stations have no indication of the depth
DF_bin <- SpatialPointsDataFrame(data.frame(DF_bin$long,DF_bin$lat), data=DF_bin)
DF_epi <- DF_bin[DF_bin@data$depth_group == "epi",]

ex.epi <- cbind(DF_epi@data, terra::extract(stacked_epi, DF_epi))
pOTU<- read.xlsx("pOTU_tax_information_Jan24.xlsx",1)
sp.c <- ex.epi[,which(colnames(ex.epi) %in% pOTU$Final_Name)] # OTU READS COUNT

LA <- pOTU[,c("Final_Name","Lat_distribution")]
LA <- LA[LA$Lat_distribution == "LA","Final_Name"]

dd4<-cbind(sp.c, ex.epi[,c(184,185,211:216)])
asvs<-dd4[,c(1:137)]
asvs[asvs>=1]<-1
cs1<-colSums(asvs, na.rm=T)

dd4.nona<-na.omit(dd4)
asvs<-dd4.nona[,c(1:137)]
asvs[asvs>=1]<-1
cs2<-colSums(asvs, na.rm=T)

cs1
cs1-cs2
round((cs1-cs2)/cs1, 2)# % lost ASVs because of NA values in raster





##########################################################################################
####---- EXTRACT RASTER VARIABLES ----######
M <- reshape2::melt(dd4.nona, id.vars = names(dd4.nona)[c(138:ncol(dd4.nona))]) # CONVERT TO LONG FORMAT (column "value" indicates READS COUNT)
M <- M[!M$variable %in% LA,]

MM <- M[rep(seq(nrow(M)), M$value),] # each row is the number of reads for each ASV at each sampling station: I repeat each row n times (n = n reads, i.e. column "value")

# EXCLUDE LOW ABUNDANCE OTUs
MM <- MM[!MM$variable %in% LA,]

## --- prepare data for Overlying Mean Index analysis (Doledec et al. 2000) --- #####

# RASTERIZE SAMPLING STATIONS
# create spatial points with sampling stations
sampling_stations <- SpatialPoints(data.frame(dd4.nona$long,dd4.nona$lat))

# rasterize the points and convert to raster
r_pts <- rasterize(sampling_stations, stacked_epi[["t"]])

# mask the original raster with the rasterized points
m_m <- mask(stacked_epi, r_pts)
#reprojected_stack <- projectRaster(m_m, crs = CRS("+init=epsg:3857"))


#### spatial pixel data frame
m = as(m_m, "SpatialPixelsDataFrame")

# image(m)
# points(lcs, col=as.numeric(slot(lcs, "data")[,20]), pch=16)

lcs <- SpatialPointsDataFrame(data.frame(MM$long,MM$lat), 
                              data=data.frame(MM[,"variable"]), 
                              proj4string=CRS("+init=epsg:4326"))

# plot(m_m[[1]])
# plot(lcs, add=T)

#### DUDI PCA #####
X <- slot(m, "data")
dim(X)
dim(na.omit(X))
pc <- dudi.pca(na.omit(X), scannf=FALSE)

#### COUNT NUMBER OF READ COUNTS FOR EACH SPECIS IN EACH PIXEL ####
options("rgdal_show_exportToProj4_warnings"="none") 
cp <- count.points(lcs, m)
U <- slot(cp, "data")
#U <- log(U+1)
# # #remove from cp rows that contain NAs in X
# rowNAS <- which(is.na(X), arr.ind=TRUE)
# U. <- U[-rowNAS[,1],]
# #U. <- U.[,-as.numeric(which(colSums(U.)<5))] #after removing NAs, one OTU has only one loc and we remove it

 
# #### RUN OVERLYING MEAN INDEX ANALYSIS: 
(ni0 <- niche(pc, U, scannf=FALSE)) #average position of each species along the ordination axes
ni0$eig[1]/sum(ni0$eig) # 0.88
ni0$eig[2]/sum(ni0$eig) # 0.06

(ni <- canomi(pc, U, scannf=FALSE)) #average position of each species along the ordination axes
ni$eig[1]/sum(ni$eig) # without currents and masked: 0.5834; without currents: 0.51 #with currents: 0.62
ni$eig[2]/sum(ni$eig) # without currents and masked: 0.3137; without currents: 0.34 #with currents: 0.25

pdf(file="CANOMI_EIGENVAL_Jan2024.pdf", width = 8, height = 6, pointsize = 12)
barplot(ni$eig, col=c("black","black", "darkgrey", "darkgrey", "darkgrey", "darkgrey"),)
dev.off()


## --- compare density of environemntal variables for each ASV with background --- ##
for(i in 1:ncol(U)){
  pdf(file=paste0(colnames(U[i]),"_HistnicheNM.pdf"), width=10)
  histniche(data.frame(m), c(U[,i]))
  #title(paste(colnames(U[i])))
  dev.off()
}

#plot(ni)

#subsss0 <- ni0$ls[sample(1:nrow(ni$ls),10000),] #subsample to get handable numbers
subsss <- ni$ls#[sample(1:nrow(ni$ls),10000),] #subsample to get handable numbers


#pdf(file="PicoZoa_canOMI_Jan24_withoutCurrents.pdf", width=18, height = 12)
par(mfrow=c(2,2))
scatterutil.eigen(ni$eig, wsel = c(1, 2))
s.arrow(ni$c1, boxes=F, sub = "Env Variables")#,xlim=c(-.1,.1))
s.label(subsss, 1, 2, clabel = 0,
        cpoint = 1.5,
        sub = "Sampling stations")
# s.distri(subsss, eval(U., sys.frame(0)), 
#          add.plot = TRUE, cstar=0)
#s.label(ni$li, boxes=T)

s.label(ni$li, boxes=T, sub = "OTUs", clabel=0, xlim=c(-5,5))
points(ni$li, 
       col=ifelse(
         rownames(ni$li) %in% cosmopolitan,"forestgreen", 
         ifelse(rownames(ni$li) %in% polar, "darkblue",
                ifelse(rownames(ni$li) %in% nonpolar, "darkred","lightgrey")))
)
text(data.frame(ni$li[,1]-0.2, ni$li[,2]-0.1), 
  labels = rownames(ni$li),
  col=ifelse(
  rownames(ni$li) %in% cosmopolitan,"forestgreen", 
  ifelse(rownames(ni$li) %in% polar, "darkblue",
         ifelse(rownames(ni$li) %in% nonpolar, "darkred","lightgrey")))
)
# s.arrow(ni$c1*2, boxes=F, sub = "Env Variables and ASVs",  
#         add.plot=T, edge = FALSE, xlim=c(-2,2))
dev.off()


ss <- data.frame(ni$ls) #resource units (sampling stations) coordinates
niche_centers <- data.frame(ni$li) #species niche coordinates (marginality vectors)
vars_arrows <- data.frame(ni$c1) #variables' vertices

library(ggrepel)
pp1 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 1) +
  geom_vline(xintercept = 0, color = "black", linetype = "solid", size = 1) +
  geom_point(data = ss, aes(x = CS1, y = CS2), color = "lightgrey", size = 1) +
  geom_point(data = niche_centers, aes(x = Axis1, y = Axis2),
             color = ifelse(
               rownames(niche_centers) %in% cosmopolitan,"forestgreen",
               ifelse(rownames(niche_centers) %in% polar, "darkblue",
                      ifelse(rownames(niche_centers) %in% nonpolar, "darkred","lightgrey"))),
             size = 2) +
  labs(title = "",#"CANOMI Analysis",
       x = "CANOMI 1 = 58.34%", y = "CANOMI 2 = 31.36%") +
  # geom_label_repel(data = niche_centers, =
  #                  aes(x = Axis1, 
  #                      y = Axis2, label = rownames(niche_centers)),
  #                  box.padding = 0, label.padding = 0,
  #                  label.size = 0, max.overlaps=50,
  #                             color = ifelse(
  #                               rownames(niche_centers) %in% cosmopolitan,"forestgreen",
  #                               ifelse(rownames(niche_centers) %in% polar, "darkblue",
  #                                      ifelse(rownames(niche_centers) %in% nonpolar, "darkred","lightgrey"))),
  #                             size = 3)+
  geom_segment(data = vars_arrows,
               aes(x = 0, y = 0, xend = CS1, yend = CS2),
               arrow = arrow(length = unit(0.1, "inches")),
               color = "black", size = .4) +
  geom_label_repel(data = vars_arrows,
            aes(x = CS1, y = CS2, label = c("Conductivity", "NO3-", "O2", "PO3-4", "Salinity", "Temperature")),
            box.padding = 0) + theme_minimal()
pp1

pdf(file="CANOMI_v1_Jan2024.pdf", width=18, height = 12)
print(pp1)
dev.off()


pp2 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 1) +
  geom_vline(xintercept = 0, color = "black", linetype = "solid", size = 1) +
  geom_point(data = ss, aes(x = CS1, y = CS2), color = "lightgrey", size = 1) +
  # geom_point(data = niche_centers, aes(x = Axis1, y = Axis2),
  #            color = ifelse(
  #              rownames(niche_centers) %in% cosmopolitan,"forestgreen",
  #              ifelse(rownames(niche_centers) %in% polar, "darkblue",
  #                     ifelse(rownames(niche_centers) %in% nonpolar, "darkred","lightgrey"))),
  #            size = 2) +
  labs(title = " ",
       x = "CANOMI 1 = 58.34%", y = "CANOMI 2 = 31.36%") +
  geom_label_repel(data = niche_centers,
                   aes(x = Axis1,
                       y = Axis2, label = rownames(niche_centers)),
                   box.padding = 0, label.padding = 0,
                   label.size = 0, max.overlaps=50,
                              color = ifelse(
                                rownames(niche_centers) %in% cosmopolitan,"forestgreen",
                                ifelse(rownames(niche_centers) %in% polar, "darkblue",
                                       ifelse(rownames(niche_centers) %in% nonpolar, "darkred","lightgrey"))),
                              size = 3)+
  geom_segment(data = vars_arrows,
               aes(x = 0, y = 0, xend = CS1, yend = CS2),
               arrow = arrow(length = unit(0.1, "inches")),
               color = "black", size = .4) +
  geom_label_repel(data = vars_arrows,
                   aes(x = CS1, y = CS2, label = c("Conductivity", "NO3-", "O2", "PO3-4", "Salinity", "Temperature")),
                   box.padding = 0) + theme_minimal()
pp2

pdf(file="CANOMI_v2_Jan2024.pdf", width=18, height = 12)
print(pp2)
dev.off()


pp3 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 1) +
  geom_vline(xintercept = 0, color = "black", linetype = "solid", size = 1) +
  geom_point(data = ss, aes(x = CS1, y = CS2), color = "lightgrey", size = 1) +
  labs(title = " ",
       x = "CANOMI 1 = 58.34%", y = "CANOMI 2 = 31.36%") +
  geom_point(data = niche_centers, aes(x = Axis1, y = Axis2),
             color = ifelse(
               rownames(niche_centers) %in% cosmopolitan,"forestgreen",
               ifelse(rownames(niche_centers) %in% polar, "darkblue",
                      ifelse(rownames(niche_centers) %in% nonpolar, "darkred","lightgrey"))),
             size = .8) +
  geom_label_repel(data = niche_centers,
                   aes(x = Axis1,
                       y = Axis2, label = rownames(niche_centers)),
                   box.padding = 0, label.padding = 0,
                   label.size = 0, max.overlaps=150,
                   color = ifelse(
                     rownames(niche_centers) %in% cosmopolitan,"forestgreen",
                     ifelse(rownames(niche_centers) %in% polar, "darkblue",
                            ifelse(rownames(niche_centers) %in% nonpolar, "darkred","lightgrey"))),
                   size = 3) +theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill="transparent", size=1))
pp3


pdf(file="CANOMI_v5_Jan2024.pdf", width = 8, height = 6, pointsize = 12)
print(pp3)
dev.off()


pp33 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 1) +
  geom_vline(xintercept = 0, color = "black", linetype = "solid", size = 1) +
  geom_point(data = ss, aes(x = CS1, y = CS2), color = "lightgrey", size = 1) +
  labs(title = " ",
       x = "CANOMI 1 = 58.34%", y = "CANOMI 2 = 31.36%") +
  geom_point(data = niche_centers, aes(x = Axis1, y = Axis2),
             color = ifelse(
               rownames(niche_centers) %in% cosmopolitan,"forestgreen",
               ifelse(rownames(niche_centers) %in% polar, "darkblue",
                      ifelse(rownames(niche_centers) %in% nonpolar, "darkred","lightgrey"))),
             size = 2) +
  # geom_label_repel(data = niche_centers,
  #                  aes(x = Axis1,
  #                      y = Axis2, label = rownames(niche_centers)),
  #                  box.padding = 0, label.padding = 0,
  #                  label.size = 0, max.overlaps=150,
  #                  color = ifelse(
  #                    rownames(niche_centers) %in% cosmopolitan,"forestgreen",
  #                    ifelse(rownames(niche_centers) %in% polar, "darkblue",
  #                           ifelse(rownames(niche_centers) %in% nonpolar, "darkred","lightgrey"))),size = 3)+
  # theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill="transparent", size=1))+theme_minimal()
pp33


pdf(file="CANOMI_v6_Jan2024.pdf", width = 8, height = 6, pointsize = 12)
print(pp33)
dev.off()

pp4 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 1) +
  geom_vline(xintercept = 0, color = "black", linetype = "solid", size = 1) +
  labs(title = "",#"CANOMI Analysis",
       x = "CANOMI 1 = 58.34%", y = "CANOMI 2 = 31.36%") +
  geom_segment(data = vars_arrows,
               aes(x = 0, y = 0, xend = CS1, yend = CS2),
               arrow = arrow(length = unit(0.1, "inches")),
               color = "black", size = .4) +
  geom_label_repel(data = vars_arrows,
                   aes(x = CS1, y = CS2, label = c("Conductivity", "NO3-", "O2", "PO3-4", "Salinity", "Temperature")),
                   box.padding = 0) + 
  scale_y_continuous(limits =c(-.5,1)) + 
  scale_x_continuous(limits =c(-.8,.8)) +
                     theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill="transparent", size=1)) 
pp4

pdf(file="CANOMI_variables_Jan2024.pdf", width = 6, height = 4, pointsize = 12)
print(pp4)
dev.off()

# library(gridExtra)
pdf(file="CANOMI_composed_Jan2024.pdf", width = 8, height = 10, pointsize = 12)
print(grid.arrange(pp3,pp4, ncol=1))
dev.off()

pdf(file="EIGEN_PicoZoa_canOMI_Jan24_withoutCurrents.pdf", width=18, height = 12)
scatterutil.eigen(ni$eig, wsel = c(1, 2))
dev.off()

pdf(file="VARS_PicoZoa_canOMI_Jan24_withoutCurrents.pdf", width=18, height = 12)
s.arrow(ni$c1, boxes=F, sub = "Env Variables", clabel=2)#,xlim=c(-.1,.1))
dev.off()



pdf(file="BACKGROUND_PicoZoa_canOMI_Jan24_withoutCurrents.pdf", width=18, height = 12)
s.label(subsss, 1, 2, clabel = 0, 
        cpoint = 1.5, 
        sub = "Available conditions")
dev.off()


pdf(file="OTUs_PicoZoa_canOMI_Jan24_withoutCurrents.pdf", width=10, height = )
s.label(ni$li, boxes=T, sub = "pOTU", 
        clabel=0, xlim=c(-5,5))
points(ni$li,
       cex=.5,
       col=ifelse(
         rownames(ni$li) %in% cosmopolitan,"forestgreen", 
         ifelse(rownames(ni$li) %in% polar, "darkblue",
                ifelse(rownames(ni$li) %in% nonpolar, "darkred","lightgrey")))
)
text(data.frame(ni$li[,1]+0.25, ni$li[,2]), 
     cex=.5,
     labels = rownames(ni$li),
     col=ifelse(
       rownames(ni$li) %in% cosmopolitan,"forestgreen", 
       ifelse(rownames(ni$li) %in% polar, "darkblue",
              ifelse(rownames(ni$li) %in% nonpolar, "darkred","lightgrey")))
)
dev.off()




## ------ Compute kernels ----- ####
##### kude
dls <- cbind(ni$ls, U)

dls_ <- melt(dls, id.vars = names(dls)[1:2])
dls. <- dls_[rep(seq(nrow(dls_)), dls_$value),]
sp.pts <- SpatialPointsDataFrame(data.frame(dls.$CS1,dls.$CS2), data=dls.)
#rsp <- as(r, "SpatialPixels")
#proj4string(sp.pts) <- CRS("+init=epsg:3035")
kud <- kernelUD(sp.pts[,"variable"], 
                h="href", grid=300, same4all = TRUE)
#image(kud)
volud <- getvolumeUD(kud, standardize = TRUE)

#mcp. <- mcp(sp.pts[,"variable"],percent = 50)
#plot(mcp.,add=T)
# image(kud)
# ka <- kernel.area(volud, percent = 95, standardize = T)
# plot(kud)
# plot(ka, add=T)

### RASTERIZE KERNELS TO COMPUTE OVERLAP
raster_list <- list()
for(i in 1:length(names(kud))){
  fud <- kud[[i]]
  hr95 <- as.data.frame(fud)[,1]
  hr95 <- data.frame(hr95)
  coordinates(hr95) <- coordinates(kud[[i]])
  gridded(hr95) <- TRUE
  raster_list[[i]] <- raster(hr95)
  #plot(raster(hr95))
}
ranges<-do.call(stack, raster_list)
names(ranges) <- names(kud)
#plot(ranges)

### KERNEL OVERLAP
homerange <- getverticeshr(volud, percent = 95)

pdf(file="Kernels_canOMI_Jan24_withoutCurrents.pdf", width=10)
s.label(subsss, 1, 2, clabel = 0,
        cpoint = 1.5,
        sub = "ASV kernels", xlim=c(-10,10))
plot(homerange, add=T, border=ifelse(homerange$id %in% polar, "blue",
                                     ifelse(homerange$id %in% nonpolar, "darkred","green"))
     )
dev.off()


par(mfrow=c(1,1))
NO <- ENMeval::calc.niche.overlap(ranges, overlapStat="D")
NO_ <- NO

full_matrix <- t(NO)
full_matrix[lower.tri(full_matrix)] <- t(full_matrix)[lower.tri(full_matrix)]

# Example correlation matrix and group variable
cor_matrix <- full_matrix

mer1 <- merge(data.frame("Final_Name"=rownames(full_matrix)), info, by="Final_Name", all.x=T, sort=F)
names(mer1) <- c("OTU","group")
group_var <- mer1$group

!any(mer1$OTU == rownames(full_matrix))

# Determine the order of the labels based on the group variable
order <- order(match(group_var, c("NP", "P", "C")))

# Reorder the correlation matrix and labels
NO_ordered <- full_matrix[order, order]

# Export Shoener's D matrix
write.csv(NO_ordered, "SchoenersD_Matrix_Ordered_Jan24.csv")

labels_ordered <- rownames(full_matrix)[order]

colors <- ifelse(labels_ordered %in% cosmopolitan,"forestgreen", 
                 ifelse(labels_ordered %in% polar, "darkblue", 
                        ifelse(labels_ordered %in% nonpolar, "red", 
                               "NA")))


pdf(file="SchoenersD_canOMI_Jan24_withoutCurrents.pdf", width=10)
corrplot::corrplot(NO_ordered, method="circle", 
                   is.corr=F, 
                   type='lower', 
                   #order="original",
                   tl.col=colors,
                   tl.cex = .7, cl.cex=.65, na.label="-",
                   col.lim=c(0,1), title = " ")

dev.off()


hr1 <- getverticeshr(kud, percent = 30)
hr2 <- getverticeshr(kud, percent = 50)
hr3 <- getverticeshr(kud, percent = 70)
hr4 <- getverticeshr(kud, percent = 90)

#par(mfrow=c(3,3))
for(i in labels_ordered){
  fud <- kud[[i]]
  volud <- getvolumeUD(fud, standardize = TRUE)
  plot_name <- paste0(i,c(" kernel"))
  
  pdf(file=paste0(plot_name,".pdf"), width=10)
  s.label(subsss*1.8, 1, 2, clabel = 0, 
          cpoint = .01, 
          sub = plot_name)#, xlim=c(-15,15), ylim=c(-15,15))
  # plot(hr1[hr1$id==names(kud)[i],], add=T, border=as.factor(homerange$id))
  # plot(hr2[hr2$id==names(kud)[i],], add=T, border=as.factor(homerange$id))
  # plot(hr3[hr3$id==names(kud)[i],], add=T, border=as.factor(homerange$id))
  # plot(hr4[hr3$id==names(kud)[i],], add=T, border=as.factor(homerange$id))
  
  colorr <- ifelse(i %in% cosmopolitan,"forestgreen", 
                   ifelse(i %in% polar, "darkblue", 
                          ifelse(i %in% nonpolar, "red", 
                                 "NA")))
  graphics::contour(volud, add=TRUE, levels = c(50,95), col=colorr, lwd=1.5, drawlabels=T, labcex=1)
  main=i
  dev.off()
}



############
############
############
############
############
############
############
############
############
############
############
############


## ------ Compute kernels ----- ####
##### kude
dls <- cbind(ni$ls, U)
dls_ <- melt(dls, id.vars = names(dls)[1:2])
dls. <- dls_[rep(seq(nrow(dls_)), dls_$value),] 
sp.pts <- SpatialPointsDataFrame(data.frame(dls.$CS1,dls.$CS2), data=dls.)
#rsp <- as(r, "SpatialPixels")
#proj4string(sp.pts) <- CRS("+init=epsg:3035")
kud <- kernelUD(sp.pts[,"variable"], 
                h="href", grid=300, same4all = TRUE)
#image(kud)
volud <- getvolumeUD(kud, standardize = TRUE)

#mcp. <- mcp(sp.pts[,"variable"],percent = 50)
#plot(mcp.,add=T)
# image(kud)
# ka <- kernel.area(volud, percent = 95, standardize = T)
# plot(kud)
# plot


### RASTERIZE KERNELS TO COMPUTE OVERLAP
raster_list <- list()
for(i in 1:length(names(kud))){
  fud <- kud[[i]]
  hr95 <- as.data.frame(fud)[,1]
  hr95 <- data.frame(hr95)
  coordinates(hr95) <- coordinates(kud[[i]])
  gridded(hr95) <- TRUE
  raster_list[[i]] <- raster(hr95)
  #plot(raster(hr95))
}
ranges<-do.call(stack, raster_list)
names(ranges) <- names(kud)
#plot(ranges)

ranges.sc <- scale(ranges)


# # Plot a kernel for each OTU
# library(rasterVis)
# library(viridis)
# gp1 <- gplot(ranges.sc) + 
#   geom_tile(aes(fill = value)) +
#   facet_wrap(~ variable) +
#   scale_fill_viridis() +
#   coord_equal()
# 
# ranges. <- subset(ranges.sc, rev(order))
# 
# pdf(file="OTU_Kernels_Feb24_vir.pdf", width=15, height = 15)
# ttpp <- terra::plot(ranges., add=TRUE, 
#             axes=FALSE, 
#             plg=list(shrink=.8), 
#             range=c(6.8,8), maxnl=nlayers(ranges.),
#             col=viridis(100),
#             legend=F)
# print(ttpp)
# dev.off()
# 
# library(rasterVis)
# pdf(file="OTU_Kernels_Feb24_v3.pdf", width=15, height = 15)
# gplot(ranges., maxpixels=5000) + 
#   geom_tile(aes(fill = value)) +
#   facet_wrap(~ variable) +
#   scale_fill_viridis_c(direction = -1) +
#   coord_equal()
# dev.off()
# 
# pdf(file="OTU_Kernels_Feb24_v4.pdf", width=15, height = 15)
# gplot(ranges., maxpixels=5000) + 
#   geom_tile(aes(fill = value)) +
#   facet_wrap(~ variable) +
#   scale_fill_gradientn(colours = rev(terrain.colors(225))) +
#   coord_equal() +
#   theme_classic()
# dev.off()
# 
# 
coords <- xyFromCell(ranges., seq_len(ncell(ranges.)))
ndvi <- stack(as.data.frame(getValues(ranges.)))
names(ndvi) <- c('value', 'variable')
ndvi <- cbind(coords, ndvi)
ndvi$colors <- ifelse(ndvi$variable %in% cosmopolitan,"forestgreen", 
                      ifelse(ndvi$variable %in% polar, "darkblue", 
                             ifelse(ndvi$variable %in% nonpolar, "red", 
                                    "NA")))

pdf(file="OTU_Kernels_Feb24_v5.pdf", width=15, height = 15)
ggplot(ndvi) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 1) +
  geom_vline(xintercept = 0, color = "black", linetype = "solid", size = 1) +
  geom_point(data = ss, aes(x = CS1, y = CS2), color = "lightgrey", size = 1) +
  labs(title = " ",
       x = "CANOMI 1", y = "CANOMI 2") +
  geom_contour(aes(x, y, z = value)) +
  facet_wrap(~ variable) +
  #scale_fill_manual(values = c("forestgreen" = "forestgreen", "Sat" = "Red", "Sun" = "Red", "Thur" = "Black"))
  #scale_fill_gradientn(colours = rev(terrain.colors(225))) +
  coord_equal()+
  theme_minimal()
dev.off()

raster_df <- as.data.frame(ranges.[[1]], xy=T)

# Plot contours using ggplot2
ggplot() +
  geom_contour(data = raster_df, aes(x = x, y = y, z = pOTU059), color = "blue") +
  theme_minimal()


### KERNEL OVERLAP
homerange <- getverticeshr(volud, percent = 95)

pdf(file="Kernels_canOMI_Jan24_withoutCurrents.pdf", width=10)
s.label(subsss, 1, 2, clabel = 0,
        cpoint = 1.5,
        sub = "Sampling stations and ASVs", xlim=c(-10,10))
plot(homerange, add=T, border=as.factor(homerange$id))
dev.off()



par(mfrow=c(1,1))
NO <- ENMeval::calc.niche.overlap(ranges, overlapStat="D")
NO_ <- NO

full_matrix <- t(NO)
full_matrix[lower.tri(full_matrix)] <- t(full_matrix)[lower.tri(full_matrix)]

# Example correlation matrix and group variable
cor_matrix <- full_matrix

mer1 <- merge(data.frame("ASV"=rownames(full_matrix)), gr, all.x=T, sort=F)
group_var <- mer1$group

!any(mer1$ASV == rownames(full_matrix))

# Determine the order of the labels based on the group variable
group_var <- factor(group_var, c("NP", "P", "C"))
order <- order(group_var)

# Reorder the correlation matrix and labels
NO_ordered <- full_matrix[order, order]
labels_ordered <- rownames(full_matrix)[order]

#gr <- read.csv("ASV_groups.csv", sep=";")
colors <- ifelse(labels_ordered %in% mer1[mer1$group == "C","OTU"],"forestgreen", 
                 ifelse(labels_ordered %in% mer1[mer1$group == "P","OTU"], "darkblue", 
                        ifelse(labels_ordered%in% mer1[mer1$group == "NP","OTU"], "red", 
                               "NA")))


for(i in 1:length(labels_ordered)){
  ids <- labels_ordered
  pdf(file=paste0("Kernels_Feb_pOTU_", ids[i], ".pdf"), width=10)
  s.label(subsss, 1, 2, clabel = 0,
          cpoint = 1.5,
          sub = "", xlim=c(-10,10))
  plot(homerange[homerange$id==ids[i],], col="darkred", add=T)
  dev.off()
}


pdf(file="SchoenersD_canOMI_Jan24_withoutCurrents.pdf", width=10)
corrplot::corrplot(NO_ordered, method="circle", 
                   is.corr=F, 
                   type='lower', 
                   #order="original",
                   tl.col=colors,
                   tl.cex = .7, cl.cex=.65, na.label="-")
dev.off()

hr1 <- getverticeshr(kud, percent = 30)
hr2 <- getverticeshr(kud, percent = 50)
hr3 <- getverticeshr(kud, percent = 70)
hr4 <- getverticeshr(kud, percent = 90)

dev.off()
#par(mfrow=c(3,3))
for(i in labels_ordered){
  fud <- kud[[i]]
  volud <- getvolumeUD(fud, standardize = TRUE)
  plot_name <- paste0(i,c(" kernel"))
  
  pdf(file=paste0(plot_name,".pdf"), width=10)
  s.label(subsss, 1, 2, clabel = 0, 
          cpoint = .01, 
          sub = plot_name)#, xlim=c(-15,15), ylim=c(-15,15))
  # plot(hr1[hr1$id==names(kud)[i],], add=T, border=as.factor(homerange$id))
  # plot(hr2[hr2$id==names(kud)[i],], add=T, border=as.factor(homerange$id))
  # plot(hr3[hr3$id==names(kud)[i],], add=T, border=as.factor(homerange$id))
  # plot(hr4[hr3$id==names(kud)[i],], add=T, border=as.factor(homerange$id))
  
  colorr <- ifelse(i %in% gr[gr$group == "Ubiquitous","ASV"],"forestgreen", 
                   ifelse(i %in% gr[gr$group == "Polar","ASV"], "darkblue", 
                          ifelse(i%in% gr[gr$group == "Non-polar","ASV"], "red", 
                                 "NA")))
  graphics::contour(volud, add=TRUE, levels = c(50,95), col=colorr, lwd=1.5, drawlabels=T, labcex=1)
  main=i
  dev.off()
}




##### KERNEL PLOT FOR EACH OTU #######
# Determine the order of the labels based on the group variable
group_var2 <- factor(group_var, c("C", "P", "NP"))
order2 <- order(group_var2)

# Reorder the correlation matrix and labels
NO_ordered2 <- full_matrix[order2, order2]
labels_ordered2 <- rownames(full_matrix)[order2]

#gr <- read.csv("ASV_groups.csv", sep=";")
colors2 <- ifelse(labels_ordered2 %in% mer1[mer1$group == "C","OTU"],"forestgreen", 
                 ifelse(labels_ordered2 %in% mer1[mer1$group == "P","OTU"], "darkblue", 
                        ifelse(labels_ordered2%in% mer1[mer1$group == "NP","OTU"], "red", 
                               "NA")))

for(i in 1:length(labels_ordered2)){
  pdf(file=paste0("Kernels_Feb_pOTU_", labels_ordered2[i], ".pdf"), width=10)
  s.label(subsss, 1, 2, clabel = 0,
          cpoint = .8,
          sub = "", xlim=c(-10,10))
  plot(homerange[homerange$id==labels_ordered2[i],], col=colors2[i], add=T, alpha=.5)
  dev.off()
}


