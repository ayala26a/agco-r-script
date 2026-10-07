#AGKISTRODON 3D GEOMETRIC MORPHOMETRICS ANALYSIS

# 19 specimens
# 3000 landmarks
# 3 dimensions (X, Y, Z)
#
# Main file:
# morphologika_unscaled.txt

# 1. CLEAR R ENVIRONMENT

rm(list = ls(all = TRUE))



rm(list=ls(all=T))
library(geomorph)
library(geiger)
library(dplyr)
library(Morpho)
library(nlme)
library(plotrix)
library(ggplot2)
installed.packages("rgl")



library(geomorph)
library(geiger)
library(dplyr)
library(Morpho)
library(nlme)
library(plotrix)
library(ggplot2)
library(rgl)



pal <- colorRampPalette(
  c("cornsilk", "darkorchid")
)



AgkShape <- read.morphologika(
  "~/morphologika_unscaled.txt"
)




dim(AgkShape)

# Expected:
# 3000 3 19




SpecimenNames <- c(
  "Agco004M",
  "Agco005M",
  "Agco010M",
  "Agco011M",
  "Agco012M",
  "Agco013M",
  "Agco015M",
  "Agco016M",
  "Agco018M",
  "Agco019M",
  "Agco020M",
  "Agco021M",
  "Agco024M",
  "Agco026M",
  "Agco027M",
  "Agco029M",
  "Agco030M",
  "Agco031M",
  "Agco22M"
)

SpecimenNames




dimnames(AgkShape)[[3]] <- SpecimenNames





AgkGPA <- gpagen(
  AgkShape,
  ProcD = FALSE
)



AgkGPA

dim(AgkGPA$coords)

# Expected:
# 3000 3 19



AgkCentroidSize <- AgkGPA$Csize

AgkCentroidSize



AgkSizeData <- data.frame(
  Specimen = SpecimenNames,
  CentroidSize = AgkCentroidSize
)

AgkSizeData


write.csv(
  AgkSizeData,
  "Agkistrodon_centroid_sizes.csv",
  row.names = FALSE
)




AgkPCA <- gm.prcomp(
  AgkGPA$coords,
  phy = NULL
)




AgkPCASum <- summary(
  AgkPCA
)

AgkPCASum

head(
  AgkPCA$x
)



AgkPCScores <- data.frame(
  Specimen = SpecimenNames,
  AgkPCA$x
)

write.csv(
  AgkPCScores,
  "Agkistrodon_PCA_scores.csv",
  row.names = FALSE
)


plot(
  AgkPCA,
  pch = 21,
  col = "black",
  bg = "grey",
  cex = 2
)


PC1percent <- round(
  AgkPCASum$PC.summary["Proportion of Variance", 1] * 100,
  2
)

PC2percent <- round(
  AgkPCASum$PC.summary["Proportion of Variance", 2] * 100,
  2
)

plot(
  AgkPCA$x[,1],
  AgkPCA$x[,2],
  
  xlab = paste0(
    "Principal Component 1 (",
    PC1percent,
    "%)"
  ),
  
  ylab = paste0(
    "Principal Component 2 (",
    PC2percent,
    "%)"
  ),
  
  pch = 21,
  col = "black",
  bg = "grey",
  cex = 2,
  
  main = "Agkistrodon Shape PCA"
)




text(
  AgkPCA$x[,1],
  AgkPCA$x[,2],
  
  labels = SpecimenNames,
  
  pos = 3,
  cex = 0.6
)



AgkColors <- pal(100)

SizeRanks <- cut(
  AgkCentroidSize,
  breaks = 100,
  labels = FALSE
)


plot(
  AgkPCA$x[,1],
  AgkPCA$x[,2],
  
  xlab = paste0(
    "Principal Component 1 (",
    PC1percent,
    "%)"
  ),
  
  ylab = paste0(
    "Principal Component 2 (",
    PC2percent,
    "%)"
  ),
  
  pch = 21,
  col = "black",
  bg = AgkColors[SizeRanks],
  cex = 2.5,
  
  main = "Agkistrodon PCA Colored by Centroid Size"
)




text(
  AgkPCA$x[,1],
  AgkPCA$x[,2],
  
  labels = SpecimenNames,
  
  pos = 3,
  cex = 0.55
)




procAgk <- procSym(
  AgkShape
)


summary(
  procAgk
)




plot(
  procAgk$PCscores[,1],
  procAgk$PCscores[,2],
  
  pch = 21,
  bg = "grey",
  col = "black",
  cex = 2,
  
  xlab = "PC1",
  ylab = "PC2",
  
  main = "Agkistrodon Morpho PCA"
)




text(
  procAgk$PCscores[,1],
  procAgk$PCscores[,2],
  
  labels = SpecimenNames,
  
  pos = 3,
  cex = 0.6
)




posPC1Agk <- restoreShapes(
  
  2 * sd(
    procAgk$PCscores[,1]
  ),
  
  procAgk$PCs[,1],
  
  procAgk$mshape
)




negPC1Agk <- restoreShapes(
  
  -2 * sd(
    procAgk$PCscores[,1]
  ),
  
  procAgk$PCs[,1],
  
  procAgk$mshape
)



deformGrid3d(
  posPC1Agk,
  negPC1Agk,
  
  ngrid = 0,
  lines = FALSE,
  pt.size = 0.1
)




posPC2Agk <- restoreShapes(
  
  2 * sd(
    procAgk$PCscores[,2]
  ),
  
  procAgk$PCs[,2],
  
  procAgk$mshape
)




negPC2Agk <- restoreShapes(
  
  -2 * sd(
    procAgk$PCscores[,2]
  ),
  
  procAgk$PCs[,2],
  
  procAgk$mshape
)




deformGrid3d(
  posPC2Agk,
  negPC2Agk,
  
  ngrid = 0,
  lines = FALSE,
  pt.size = 0.1
)




AgkData <- geomorph.data.frame(
  
  shape = AgkGPA$coords,
  
  Csize = AgkGPA$Csize,
  
  specimen = SpecimenNames
)



AgkData






logCsize <- log(
  AgkGPA$Csize
)




plot(
  logCsize,
  AgkPCA$x[,1],
  
  pch = 21,
  bg = "grey",
  col = "black",
  cex = 2,
  
  xlab = "Log Centroid Size",
  
  ylab = "PC1 Score",
  
  main = "PC1 vs Centroid Size"
)




PC1sizeModel <- lm(
  AgkPCA$x[,1] ~ logCsize
)

abline(
  PC1sizeModel,
  lwd = 2
)




summary(
  PC1sizeModel
)




plot(
  logCsize,
  AgkPCA$x[,2],
  
  pch = 21,
  bg = "grey",
  col = "black",
  cex = 2,
  
  xlab = "Log Centroid Size",
  
  ylab = "PC2 Score",
  
  main = "PC2 vs Centroid Size"
)




PC2sizeModel <- lm(
  AgkPCA$x[,2] ~ logCsize
)

abline(
  PC2sizeModel,
  lwd = 2
)



summary(
  PC2sizeModel
)




AgkShapeSizeModel <- procD.lm(
  
  AgkGPA$coords ~
    log(AgkGPA$Csize),
  
  iter = 999
)



summary(
  AgkShapeSizeModel
)




anova(
  AgkShapeSizeModel
)




plot(
  AgkShapeSizeModel,
  
  type = "regression",
  
  predictor = log(
    AgkGPA$Csize
  )
)




AgkShapeDistances <- dist(
  two.d.array(
    AgkGPA$coords
  )
)

AgkShapeDistances



AgkDistanceMatrix <- as.matrix(
  AgkShapeDistances
)

rownames(
  AgkDistanceMatrix
) <- SpecimenNames

colnames(
  AgkDistanceMatrix
) <- SpecimenNames




write.csv(
  AgkDistanceMatrix,
  "Agkistrodon_shape_distance_matrix.csv"
)


PCAvariance <- AgkPCA$sdev^2

PCApercent <- (
  PCAvariance /
    sum(PCAvariance)
) * 100



AgkPCVariance <- data.frame(
  
  PC = paste0(
    "PC",
    seq_along(PCApercent)
  ),
  
  PercentVariance = PCApercent,
  
  CumulativeVariance =
    cumsum(PCApercent)
)




head(
  AgkPCVariance,
  10
)




write.csv(
  AgkPCVariance,
  
  "Agkistrodon_PCA_variance.csv",
  
  row.names = FALSE
)




barplot(
  AgkPCVariance$PercentVariance[1:10],
  names.arg = AgkPCVariance$PC[1:10],
  xlab = "Principal Component",
  ylab = "Percent Variance Explained",
  main = "Agkistrodon PCA Variance",
  las = 2,
  cex.names = 0.8
)



AgkMeanShape <- mshape(
  AgkGPA$coords
)

AgkMeanShape



rgl::plot3d(
  AgkMeanShape,
  type = "s",
  size = 2
)


LowestPC1 <- SpecimenNames[
  which.min(
    AgkPCA$x[,1]
  )
]

LowestPC1




HighestPC1 <- SpecimenNames[
  which.max(
    AgkPCA$x[,1]
  )
]

HighestPC1



LowestPC2 <- SpecimenNames[
  which.min(
    AgkPCA$x[,2]
  )
]

LowestPC2




HighestPC2 <- SpecimenNames[
  which.max(
    AgkPCA$x[,2]
  )
]

HighestPC2




ExtremeSpecimens <- data.frame(
  
  Direction = c(
    "Negative PC1",
    "Positive PC1",
    "Negative PC2",
    "Positive PC2"
  ),
  
  Specimen = c(
    LowestPC1,
    HighestPC1,
    LowestPC2,
    HighestPC2
  )
)

ExtremeSpecimens




write.csv(
  ExtremeSpecimens,
  
  "Agkistrodon_extreme_specimens.csv",
  
  row.names = FALSE
)




AgkResults <- data.frame(
  
  Specimen = SpecimenNames,
  
  CentroidSize =
    AgkGPA$Csize,
  
  PC1 =
    AgkPCA$x[,1],
  
  PC2 =
    AgkPCA$x[,2],
  
  PC3 =
    AgkPCA$x[,3]
)



AgkResults




write.csv(
  
  AgkResults,
  
  "Agkistrodon_main_results.csv",
  
  row.names = FALSE
)



pdf(
  "Agkistrodon_All_Graphs.pdf",
  width = 10,
  height = 8
)
plot(
  AgkPCA$x[,1],
  AgkPCA$x[,2],
  xlab = paste0(
    "Principal Component 1 (",
    PC1percent,
    "%)"
  ),
  ylab = paste0(
    "Principal Component 2 (",
    PC2percent,
    "%)"
  ),
  pch = 21,
  col = "black",
  bg = "grey",
  cex = 2,
  main = "Agkistrodon Shape PCA"
)

text(
  AgkPCA$x[,1],
  AgkPCA$x[,2],
  labels = SpecimenNames,
  pos = 3,
  cex = 0.6
)



pdf(
  "Agkistrodon_All_Graphs.pdf",
  width = 10,
  height = 8
)

############################################################
# GRAPH 1: BASIC PCA
############################################################

plot(
  AgkPCA,
  pch = 21,
  col = "black",
  bg = "grey",
  cex = 2
)


############################################################
# GRAPH 2: PC1 VS PC2 WITH SPECIMEN LABELS
############################################################

plot(
  AgkPCA$x[,1],
  AgkPCA$x[,2],
  
  xlab = paste0(
    "Principal Component 1 (",
    PC1percent,
    "%)"
  ),
  
  ylab = paste0(
    "Principal Component 2 (",
    PC2percent,
    "%)"
  ),
  
  pch = 21,
  col = "black",
  bg = "grey",
  cex = 2,
  
  main = "Agkistrodon Shape PCA"
)

text(
  AgkPCA$x[,1],
  AgkPCA$x[,2],
  
  labels = SpecimenNames,
  pos = 3,
  cex = 0.6
)


############################################################
# GRAPH 3: PCA COLORED BY CENTROID SIZE
############################################################

AgkColors <- pal(100)

SizeRanks <- cut(
  AgkCentroidSize,
  breaks = 100,
  labels = FALSE
)

plot(
  AgkPCA$x[,1],
  AgkPCA$x[,2],
  
  xlab = paste0(
    "Principal Component 1 (",
    PC1percent,
    "%)"
  ),
  
  ylab = paste0(
    "Principal Component 2 (",
    PC2percent,
    "%)"
  ),
  
  pch = 21,
  col = "black",
  bg = AgkColors[SizeRanks],
  cex = 2.5,
  
  main = "Agkistrodon PCA Colored by Centroid Size"
)

text(
  AgkPCA$x[,1],
  AgkPCA$x[,2],
  
  labels = SpecimenNames,
  pos = 3,
  cex = 0.55
)

# Simple legend showing size range
legend(
  "topright",
  legend = c(
    paste("Small:", round(min(AgkCentroidSize), 1)),
    paste("Large:", round(max(AgkCentroidSize), 1))
  ),
  pch = 21,
  pt.bg = c(AgkColors[1], AgkColors[100]),
  pt.cex = 2,
  bty = "n",
  title = "Centroid Size"
)


############################################################
# GRAPH 4: MORPHO PCA
############################################################

plot(
  procAgk$PCscores[,1],
  procAgk$PCscores[,2],
  
  pch = 21,
  bg = "grey",
  col = "black",
  cex = 2,
  
  xlab = "PC1",
  ylab = "PC2",
  
  main = "Agkistrodon Morpho PCA"
)

text(
  procAgk$PCscores[,1],
  procAgk$PCscores[,2],
  
  labels = SpecimenNames,
  pos = 3,
  cex = 0.6
)


############################################################
# GRAPH 5: PC1 VS CENTROID SIZE
############################################################

logCsize <- log(AgkGPA$Csize)

PC1sizeModel <- lm(
  AgkPCA$x[,1] ~ logCsize
)

plot(
  logCsize,
  AgkPCA$x[,1],
  
  pch = 21,
  bg = "grey",
  col = "black",
  cex = 2,
  
  xlab = "Log Centroid Size",
  ylab = "PC1 Score",
  
  main = "PC1 vs Centroid Size"
)

abline(
  PC1sizeModel,
  lwd = 2
)


############################################################
# GRAPH 6: PC2 VS CENTROID SIZE
############################################################

PC2sizeModel <- lm(
  AgkPCA$x[,2] ~ logCsize
)

plot(
  logCsize,
  AgkPCA$x[,2],
  
  pch = 21,
  bg = "grey",
  col = "black",
  cex = 2,
  
  xlab = "Log Centroid Size",
  ylab = "PC2 Score",
  
  main = "PC2 vs Centroid Size"
)

abline(
  PC2sizeModel,
  lwd = 2
)


############################################################
# GRAPH 7: PCA VARIANCE / SCREE PLOT
############################################################

AgkPCVariance <- data.frame(
  PC = paste0(
    "PC",
    seq_len(ncol(AgkPCA$x))
  ),
  
  PercentVariance =
    AgkPCASum$PC.summary[
      "Proportion of Variance",
    ] * 100
)

barplot(
  AgkPCVariance$PercentVariance[1:10],
  
  names.arg =
    AgkPCVariance$PC[1:10],
  
  xlab = "Principal Component",
  ylab = "Percent Variance Explained",
  
  main = "Agkistrodon PCA Variance",
  
  las = 2,
  cex.names = 0.8
)



dev.off()



getwd()
