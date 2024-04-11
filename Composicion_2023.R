# Analsis de Diversidad de especies:
library(vegan)
library(lattice)
library(permute)
library (ggplot2)
###An?lisis de similitud (COMPOSICI?N)
#Archivo: "AnoSim2.txt" ###Deben estar las Spp como columnas###

####DIVIDO LA COMPOSICION POR ESTRATOS ###

#########ESTRA ARBOREO##########

veg<- read.table(file.choose(), header = TRUE)
head(veg) #base de datos veg
str(veg) 
as.matrix(veg)

veg$Nelchi <- as.numeric(veg$Nelchi)
veg$Tepart <- as.numeric(veg$Tepart)
veg$Sempra <- as.numeric(veg$Sempra)
veg$Opuqui <- as.numeric(veg$Opuqui)
veg$Langri <- as.numeric(veg$Langri)
veg$Deiurb <- as.numeric(veg$Deiurb)

fix(veg) #Selgil tiene un ?

veg$Sengil <- as.numeric(veg$Sengil)

#Analysis of Similarities
veg.dist<-vegdist(veg, "bray") # para datos de abundancia se usa bray, la primer celda de la base de datos debe estar VACIA 
summary(veg.dist)
head(veg.dist)

sitios<- read.table(file.choose(), header = TRUE)#cargo la info de los sitios. Archivo"Sitios.txt" o"sitios+alt+DAS"
attach(sitios)

dispersion<-betadisper(veg.dist, group=LU)
permutest(dispersion) 
mod.HSD <- TukeyHSD(dispersion) # no da significativa puedo seguir 
plot(mod.HSD)

#La matriz.de.distancias es la que obtienes con vegdist
#El grupo en tu caso serÃ­a la condic de uso (Lu)
data.ano.v <- anosim(veg.dist, LU)
summary(data.ano.v)

#Dissimilarity ranks between and within classes:
 #0%   25%   50%    75% 100%  N
#Between  2 40.75  65.5  95.50  120 84
#C       36 90.75 100.5 111.75  118  6
#RWC      1  8.50  30.0  62.00   86 15
#WV       4 22.50  34.0  60.50   91 15

plot(data.ano.v)
#R=0.30 P=0.002

#The ANOSIM statistic compares the mean of ranked dissimilarities between groups to 
#the mean of ranked dissimilarities within groups. An R value close to "1.0" suggests
#dissimilarity between groups while an R value close to "0" suggests an even distribution 
#of high and low ranks within and between groups. R values below "0" suggest that dissimilarities are greater within groups
#than between groups. See Clarke and Gorley (2001) for a guide to interpreting ANOSIM R values.

veg.adonis <- adonis2(veg.dist ~ LU,  contr.unordered = "contr.sum", data = veg)
veg.adonis


#Df SumOfSqs      R2      F Pr(>F)  
#LU        2  0.59738 0.21785 1.8104  0.014 * # este
#Residual 13  2.14480 0.78215                
#Total    15  2.74218 1.00000                
---

##############################################################################
#Gr?ficos de ordenaci?n: Non-metric multidimensional scaling NMDS
#Gráficos de ordenación: Non-metric multidimensional scaling NMDS


library(ggplot2)
ord1 <- metaMDS(veg, distance = "bray")
plot(ord1, type = "t")


plot(ord1, disp="sites", type="n")
ordihull(ord1, LU, col=1:2, lwd=2)
ordiellipse(ord1, LU, col=1:2, kind = "ehull", lwd=2)
ordiellipse(ord1, LU, col=1:2, draw="polygon")
ordispider(ord1, LU, col=1:2, label = TRUE)
points(ord1, disp="sites", pch=21, col="red", bg="yellow",col="green", cex=1.3)
text(ord1, display = "spec", cex=0.6, col="blue")


NMDS1 <- ord1$points[,1] ##tiro esto por que no me reconoce los NMDS1 Y 2 en ggplot
NMDS2 <- ord1$points[,2]



## solo los vectores significativos
library(ggplot2)
p1 <-ggplot(veg, aes(NMDS1, NMDS2))+
  geom_point(aes(NMDS1, NMDS2, color=LU),position=position_jitter(.1))+##separates overlapping points
  stat_ellipse(aes(fill=LU), alpha=.2,type='t',size =1, geom="polygon")


p1

# We can see that the 3 land use types shere a great part of their plants species, but the sites with active livestock production systems presente more heterogenity. 
# Maybe the land use management is introducing heterogenity and increasing diversity but disruting the state of plant community. 

#el stress indice permite saber si las diferencias y similitudes de los datos
#ºestan bien expresadas en el graficos. En este caso si. hay que fijrse q tan dispersos 
#estan los puntos de la linea

stressplot(ord2)

#Clarke 1993 suggests the following guidelines for acceptable stress values: 
#<0.05 = excellent, <0.10 = good, <0.20 = usable, >0.20 = not acceptable. The plot shows the border of the 0.20 stress value limit. Solutions with higher 
#stress values should be interpreted with caution and those with stress above 0.30 are highly suspect.

metaMDS(
  veg,
  distance = "bray",
  k = 2,
  trymax = 20,
  autotransform = TRUE) # 0.1889742  el valor de stress es bueno

#bueno


#################################
require(graphics)
install.packages("ggdendro")
library("ggplot2")
library(ggdendro)
head(data)


a <- hclust(veg.dist, "ave") #########mepa q es este 
a <- hclust(dist(veg), "ave")


a <- hclust(veg.dist, "ave")
plot(a)
plot(a, hang = -1) # aca se ve que lo sitios estan re mezclados 
# la condicion no explica ni la variable del suelo ni la composicion de la vegetacion 

ggdendrogram(a)

# Build dendrogram object from hclust results
dend <- as.dendrogram(a)
# Extract the data (for rectangular lines)
# Type can be "rectangle" or "triangle"
dend_data <- dendro_data(dend, type = "rectangle")
# What contains dend_data
names(dend_data)
head(dend_data$segments)
head(dend_data$labels)
