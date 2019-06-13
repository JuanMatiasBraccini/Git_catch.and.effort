################ Estimation des bycatchs de l'OI : 2003-2009
rm(list=ls())
pack<-c("RODBC", "maps", "mapdata", "mapproj", "MASS", "xtable", "fields", "spatstat", "RGeoS", "Hmisc", "mice")
for(i in 1:length(pack)){library(pack[i], character.only = TRUE)}
#  Appel des functions personnelles cr��es pour mes analyses
source("D:/THESE/Algorithmes utiles/Compilation_Fonction_Utiles.txt")
#  Lecture et visualisation de l'ensemble des donn�es par table
source("D:/THESE/Algorithmes utiles/lecture_donnee_dcr_OI.txt")
#  Repertoire de travail
setwd("D:/THESE/Analyse/Estimation_bycatch_OI")
donnee_generale <- subset(donnee_generale, donnee_generale$TYPBANC1 !="IND") 
donnee_generale$TYPBANC1 <- as.factor(as.vector(donnee_generale$TYPBANC1)) ### Pour faire disparaitre MODE DE PECHE == IND


### Nombre de mar�e
length(unique(donnee_generale$NUM_MAREE))
length(donnee_generale$Numlance)
a <- unique(donnee_generale[, c("NUM_MAREE", "D_FIN")])
a$annee <- substr(a$D_FIN, 7,10)
a <- unique(a[, c("NUM_MAREE", "annee")]) ; a
table(a)
colSums(table(a))

with(donnee_generale, tapply(Numlance, list(TYPBANC1), length)) ### Nombre de coups de senne

############################
### NIVEAU DE PRECISION SUR LES ESTIMATIONS DE BYCATCH EN FONCTION DE L'ESPECE ET DU TAUX DE COUVERTURE... Estimation par le ratio .. Stratification FSC/FAD.

donn <- Requete_Calcul_Debarq_Rejet_Totaux_thonides(donnee_generale)
#espece <- c("COH", "ELP", "CFA", "CLO", "WAH", "BCM", "RTY", "LOB", "BLM", "FBA", "BAT", "DVI")
espece <- as.vector(unique(capfauna$C_ESP_F_3L))
# espece <- c("BCM", "ELP", "COH", "CFA", "WAH", "FCR", "CUX", "MZZ", "CLO", "SAP", "COE", "LKE", "RMV", "MCO")
#espece <- c("BCM", "ELP", "COH", "CFA")
donn <- Requete_Calcul_Nombre_Poids_Espece_Associee(donn, espece, "poids")
head(donn, 3); dim(donn)
a <- sort(colSums(donn[, espece])); a
b <- a[a > 0] ; b
espece <-  as.vector(as.factor(names(b))) ; espece
donn$Bycatch <- rowSums(donn[, c(c("tunalike_dbq", "tuna_rej", "Billfish", "Shark", "Rays", "Bonyfish", "Turtles"))])
donn <- subset(donn, is.na(donn$INST)==F)
head(donn)
#write.table(donn, "donn.txt", sep="\t", dec=".", row.names=F)  ## Sauvegarde pour serveur

with(donn, tapply(major_tuna_dbq, list(TYPBANC1), sum)) ### Nombre de coups de senne


# source("script_pour_cluster.r")
getwd()
dir()

#################################################################################
#################  RESUME -- interpr�tation des r�sultats de la simulation --
res <- read.table("res.txt", h=T)
res <- subset(res, is.na(res$rmse)==F & is.na(res$posOccur)==F & res$coverage >=4)
head(res); dim(res)
tab_index <- unique(res[, c("specie", "posOccur", "gini")])

### Modele GLM 
tab <- res
head(tab)
hist(tab$rmse, 1000)
mod <- glm(rmse ~ coverage*posOccur*gini, data = tab, family=Gamma)
summary(mod)
anova(mod, test="Chisq")
round(100*(1 - mod$deviance/mod$null.deviance), 2) ## Part de d�viance expliqu�e par le mod�le

write.table(summary(mod)$coefficients,  file="tableaux/summaryGLM.xls", dec=".", sep="\t")
anov <- anova(mod, test="Chisq"); anov
anov$Explained.Dev <- round(100*anov[, "Deviance"]/anov[, "Resid. Dev"][1], 2)  ; anov
sum(anov$Explained.Dev, na.rm=T)    ## Part de d�viance expliqu�e par le mod�le
write.table(anov,  file="tableaux/anovaGLM.xls", dec=".", sep="\t")


##################### FIGURES ARTICLE  ###################
######## Figure  2: Boxplot of the relative mean square error (RMSE) of all bycatch species estimates conditionally to the observer sampling coverage.

head(res)
a <- unique(res[, c("specie","posOccur", "gini")])
a$coverage <- 100
a$rmse <- 0
res <- rbind(res, a)
#pdf("D:/THESE/Analyse/Estimation_bycatch_OI/figures/figure2.pdf")
#postscript("D:/THESE/Analyse/Estimation_bycatch_OI/figures/figure2.eps")
tiff("D:/THESE/Analyse/Estimation_bycatch_OI/figures/figure2.tif", , bg="transparent", res=1100, width =75, height = 75, units="mm", pointsize = 10)
par(mar=c(3.5, 4, 0.5, 1.05))
pas <- 10
plot(0, xlim=c(0, 100)/pas, ylim=c(0,800), type="n", axes=F, ylab="", xlab="", xaxs="i", yaxs="i")
cover <- c(5, seq(10,100,10))
for(i in 1:length(cover)){
with(subset(res, res$coverage==cover[i]), boxplot(rmse ~ coverage, notch=T, at=cover[i]/pas, lty=1, lwd=0.5, add=T, col="grey", axes=F))}
mtext("RMSE (%)", 2, line=3)
mtext("Coverage rate (%)", 1, line=2)
axis(1, seq(0, 100, 10)/pas, seq(0, 100, 10), pos=0, lwd=0, lwd.ticks=0.5)
axis(2,  seq(0, 900, 100), las=2, pos=0, lwd=0, lwd.ticks=0.5)
box(lwd=0.5)
dev.off()


### Figure 2 avec la somme des co�ts  (PS : AGGRANDIR LA FENETRE AVANT DE FAIRE LE GRAPHE)
couverture <- c(5, seq(10,100,10)); couverture
seuil_rmse <- seq(1, 200, 2)
prop <- matrix(0, nrow=length(seuil_rmse), ncol=length(couverture))
for(i in 1:nrow(prop)){ 
for(j in 1:ncol(prop)) { 
prop[i, j] <- 100*nrow(subset(res, res$coverage==couverture[j] & rmse < seuil_rmse[i]))/nrow(subset(res, res$coverage==couverture[j])) }} 
dimnames(prop) <- list(seuil_rmse, couverture)

pdf("D:/THESE/Analyse/Estimation_bycatch_OI/figures/figure2_cout.pdf")
par(mar=c(7,10,0,1))
image.plot(as.numeric(dimnames(prop)[[1]]), as.numeric(dimnames(prop)[[2]]), prop, ylab="Coverage rate (%)", xlab="RMSE (%)", cex.lab=2, cex.axis=1.5, las=1, horizontal = T, legend.width = 1, legend.args=list(text="Proportion of bycatch species", cex=2, cex.lab=1.5, line=0.5))

spp <- c("CFA", "DVI", "LKE","WHM")
nom <- c("C. falciformis", "D. violacea", "L. kempii", "T. albidus")
couleur <- 1:length(spp)
pchh <- 1:length(spp)
lwdd <- 2
ltyy <- 1:length(spp)
for(i in 1:length(spp)){
aa <-subset(res, res$specie==spp[i]) ; aa <- aa[order(aa$rmse), ]
with(aa, lines(rmse, coverage, col=couleur[i], lwd=lwdd, lty=ltyy[i])); with(aa, points(rmse, coverage, pch=pchh[i], col=couleur[i], lwd=lwdd, lty=ltyy[i])) }
legend("topright", nom, pch=pchh, col=couleur, cex=1.2, lwd=lwdd, lty=ltyy, box.col="white", bg="white")
box()



# Ajout de l'information de co�t
cout_jour <- 242
nb_jr_par_mar <- 40
cout_maree <- cout_jour*nb_jr_par_mar ; cout_maree
nb_tot_maree <- 355 ## Moyenne annuelle sur la p�riode 2003-2009
cout_total <- nb_tot_maree*cout_maree; cout_total
# Cout par unit� de couverture
cout_cov <- cout_total/100; cout_cov
round(seq(10, 90, 10)*cout_cov/1000000, 1)
axis(2, seq(0, 100, 20), round(seq(0, 100, 20)*cout_cov/1000000, 1), cex.axis=1.5, pos=-40, las=1)
mtext("Cost (M�)", 2, line=8.5, cex=2)
lines(c(50, 50), c(0, 100), lty=2, lwd=2)
#mtext(50, 1, at=50, line=0.5, cex=1.5)
dev.off()

## Version NOIR & BLANC
tiff("D:/THESE/Analyse/Estimation_bycatch_OI/figures/figure4.tif", bg="transparent", res=1100, width =75, height = 75, units="mm", pointsize = 5)
par(mar=c(7,10,0,1))
image.plot(as.numeric(dimnames(prop)[[1]]), as.numeric(dimnames(prop)[[2]]), prop, ylab="Coverage rate (%)", xlab="RMSE (%)", cex.lab=2, cex.axis=1.5, las=1, horizontal = T, col=gray(seq(1,0.15,-0.0001)), legend.width = 1, legend.args=list(text="Proportion of bycatch species", cex=1.5, cex.lab=2, line=0.5))

spp <- c("CFA", "DVI", "LKE","WHM")
nom <- c("C. falciformis (shark)", "D. violacea (ray)", "L. kempii (turtle)", "T. albidus (billfish)")
couleur <- c(1,1,"white", "white")
pchh <- 1:length(spp)
lwdd <- 0.5
ltyy <- 1:length(spp)
for(i in 1:length(spp)){
aa <-subset(res, res$specie==spp[i]) ; aa <- aa[order(aa$rmse), ]
with(aa, lines(rmse, coverage, col=couleur[i], lwd=lwdd, lty=ltyy[i])); with(aa, points(rmse, coverage, pch=pchh[i], col=couleur[i], lwd=lwdd, lty=ltyy[i])) }
legend("topright", nom, pch=pchh, col=1, cex=1.2, lwd=lwdd, lty=ltyy, box.col="white", bg="white")
box(lwd=0.5)

# Ajout de l'information de co�t
cout_jour <- 242
nb_jr_par_mar <- 40
cout_maree <- cout_jour*nb_jr_par_mar ; cout_maree
nb_tot_maree <- 355 ## Moyenne annuelle sur la p�riode 2003-2009
cout_total <- nb_tot_maree*cout_maree; cout_total
# Cout par unit� de couverture
cout_cov <- cout_total/100; cout_cov
round(seq(10, 90, 10)*cout_cov/1000000, 1)
axis(2, seq(0, 100, 20), round(seq(0, 100, 20)*cout_cov/1000000, 1), cex.axis=1.5, pos=-40, las=1)
mtext("Cost (M�)", 2, line=8.5, cex=2)
lines(c(50, 50), c(0, 100), lty=2, lwd=1)
#mtext(50, 1, at=50, line=0.5, cex=1.5)
dev.off()




######### Figure 3: Relative mean square error (RMSE) of bycatch estimates as a function of positive occurrence in log-scale (x-axis) and conditional Gini index (y-axis). 

head(res)
res$predd <- predict(mod, type="response", newdata=res )

tab <- subset(res, res$coverage==10 & !(res$gini < 0.05 &  log(res$posOccur, 10) < -3.28))
head(tab); dim(tab)
db <- with(tab, db.create(x1=log(posOccur, 10), x2=gini, z1=predd))
plot(db)    
w <- convexhull.xy(log(res$posOccur,10),  res$gini)
plot(w$bdry[[1]]$x, w$bdry[[1]]$y) ; polygon(w$bdry[[1]]$x, w$bdry[[1]]$y)
x <- c(c(-0.2546528, -0.2714748) + 1, -0.9821456, -1.5075039, -2.4831592, -3.1832698, -3.4844422, -3.0070362, -1.5325054)  ## elargissement de la fen�tre
y <- c(0.6408339-0.03, c(0.7003376, 0.9481568, 0.9546443, 0.7879245, 0.4560299)+ 0.03 , c(0,0), 0.1937891-0.03)               ## elargissement de la fen�tre
lines(x, y, col=2)
Mypolygon <- polygon.create(x , y, polygon=NA)
Mypolygon <- polygon.create(w$bdry[[1]]$x, w$bdry[[1]]$y, polygon=NA)
plot(Mypolygon, col="grey", add=T)
gridd <-  db.create(flag.grid=T, x0=c(-3.5, 0), dx=0.25*c(0.1, 0.05), nx=4*c(35, 20))  
aa <- invdist(db, db.polygon(gridd, Mypolygon), exponent = 3)
aa 
aa <- aa@items    ; head(aa)
aa <- subset(aa, aa$sel==T)
M <- with(aa, tapply(Invdist.z1, INDEX=list(x1, x2), FUN=mean)) 
dev.off()


tiff("D:/THESE/Analyse/Estimation_bycatch_OI/figures/figure3.tif", bg="transparent", res=1100, width =75, height = 75, units="mm", pointsize = 5)
par(mar=c(5,5,0.5, 5))
plot(0, type="n", xlim=c(-3.5, -0.35), ylim=c(0,1), main="", xlab="Percentage of set with positive bycatch occurrence", ylab="Gini index", cex.lab=1.8, axes=F)
i <- seq(-4, 0.5, 0.5) 
axis(1, i , c(round(100*10^seq(-4, 0.5, 0.5), 2)[1:3], round(100*10^seq(-4, 0.5, 0.5), 1)[4], round(100*10^seq(-4, 0.5, 0.5), 0)[5:length(i)]), cex.axis=1.5, lwd=0, lwd.ticks=0.5)
axis(2, seq(0, 1, 0.2), seq(0, 1, 0.2), cex.axis=1.5, las=2, lwd=0, lwd.ticks=0.5)
image.plot(unique(aa$x1), unique(aa$x2), M, col=gray(seq(0, 1, 0.02)), add=T)
#image.plot(unique(aa$x1), unique(aa$x2), M, col=terrain.colors(100), add=T)
#contour(unique(aa$x1), unique(aa$x2), M, add=T)

head(capfauna)
zzz <- unique(capfauna[, c("Group_bis", "C_ESP_F_3L")])
names(zzz)[2] <-"specie"
zzz <- subset(zzz, !(zzz$Group_bis %in% c("WhaleShark", "Cetaceans")))
head(zzz)

group_res <- merge(subset(res, res$coverage==10 & !(res$gini < 0.05 &  log(res$posOccur, 10) < -3.28)), zzz, all=F)
gp <- as.vector(unique(group_res$Group_bis)) ; gp
for(i in 1:5){
zz <- unique(subset(group_res, group_res$coverage > 3 & group_res$Group_bis == gp[i])[, c("specie", "posOccur", "gini")])
w <- convexhull.xy(log(zz$posOccur,10),  zz$gini)
polygon(w$bdry[[1]]$x, w$bdry[[1]]$y, border=i+1, density=NULL, lwd=1)
with(zz, text(log(posOccur, 10), gini, specie, cex=0.8, col=i+1, font=2))
}
gp <- c("Bony fishes","Sharks", "Billfishes", "Turtles", "Rays")  ## Rearrangement des noms des groupes d'esp�ces
legend("topleft", as.vector(gp), lty=1, col=2:7, cex=1, box.col="white")
box(lwd=0.5)
#with(unique(res[, c("specie", "posOccur", "gini")]), points(log(posOccur, 10), gini, cex=0.8, col=i+1, font=2))
dev.off()


######### Figure 4: Proportion of species that gave a relative mean square error (RMSE) lower than a reference of 20, 30, 40 and 50%, given the observer sampling coverage rate.
head(res)

couverture <- c(5, seq(10,90,10)); couverture
seuil_rmse <- seq(20, 50, 10)
prop <- matrix(0, nrow=length(seuil_rmse), ncol=length(couverture))
for(i in 1:nrow(prop)){ 
for(j in 1:ncol(prop)) { 
prop[i, j] <- 100*nrow(subset(res, res$coverage==couverture[j] & rmse < seuil_rmse[i]))/nrow(subset(res, res$coverage==couverture[j])) }} 

par(mar=c(5,5,2,1))
a <- barplot(prop, beside=T, axes=F, ylim=c(0,101), ylab="Proportion of species (%)", xlab="Coverage rate (%)", cex.lab=2.2, col="white")
grid(nx=NA, ny=c(), col="grey")
barplot(prop, beside=T, axes=F, add=T, col=gray(seq(0, 0.85, length=4)))
axis(1, colMeans(a), couverture, cex.axis=1.8)
axis(2, las=2, cex.axis=1.7) 
legend("topleft", paste(paste("RMSE <", seuil_rmse), "%"), col=gray(seq(0, 0.85, length=4)), pch=15, cex=2, box.col="white")
box()



















######### Autre figure : Gain marginal de pr�cision
head(res)
summary(mod)

comb <- unique(res[, c("posOccur", "gini")])

COVERAGE <- seq(1, 90, 1)
pos <- c()
gain_rmse <- c()
#toto <- c()
k <- 0

for(p in 1:nrow(comb)){
posOccur <- comb$posOccur[p]
gini <- comb$gini[p]
		for(co in 1:length(COVERAGE)){
		coverage <- COVERAGE[co]
rmse1 <- predict(mod, type="response", newdata=data.frame(gini, posOccur, coverage=coverage))
rmse2 <- predict(mod, type="response", newdata=data.frame(gini, posOccur, coverage= coverage+1))

k <- k+1
pos[k] <- (2*coverage+1)/2
gain_rmse[k] <- as.numeric(rmse2) - as.numeric(rmse1)
#gain_rmse[k] <- round(100*(as.numeric(rmse2) - as.numeric(rmse1))/as.numeric(rmse1), 2)
#toto[k] <- paste(gini, posOccur, sep=" - ")
}}

#### figure 5_bis
cout_jour <- 242
nb_jr_par_mar <- 40
cout_maree <- cout_jour*nb_jr_par_mar ; cout_maree
nb_tot_maree <- 355 ## Moyenne annuelle sur la p�riode 2003-2009
cout_total <- nb_tot_maree*cout_maree; cout_total
# Cout par unit� de couverture
cout_cov <- cout_total/100; cout_cov


ii <- boxplot(-gain_rmse ~ pos, cex=0.4, pch=16, axes=T, ylab="Value in �",xlab="Coverage rate (%)", cex.lab=2, main="lambda=1000", col="grey")

par(mar=c(5,5,5,5), mfrow=c(2,2))
matplot(t(ii$stats[, 1:30]), type="l", col=1, lty=c(3,2,1,2,3), lwd=c(1,1,2,1,1), ylab="Value in �", xlab="Coverage rate (%)", cex.lab=1.5, cex.axis=1.5, main="lambda=1")
 
matplot(t(ii$stats[, 1:30])*1000, type="l", col=1, lty=c(3,2,1,2,3), lwd=c(1,1,2,1,1), ylab="Value in �", xlab="Coverage rate (%)", cex.lab=1.5, cex.axis=1.5, main="lambda=1000")
lines(seq(0,100,100), cout_cov + 0*seq(0,100,100))

matplot(t(ii$stats[, 1:30])*2000, type="l", col=1, lty=c(3,2,1,2,3), lwd=c(1,1,2,1,1), ylab="Value in �", xlab="Coverage rate (%)", cex.lab=1.5, cex.axis=1.5, main="lambda=2000")
lines(seq(0,100,100), cout_cov + 0*seq(0,100,100))

matplot(t(ii$stats[, 1:30])*5000, type="l", col=1, lty=c(3,2,1,2,3), lwd=c(1,1,2,1,1), ylab="Value in �", xlab="Coverage rate (%)", cex.lab=1.5, cex.axis=1.5, main="lambda=5000")
lines(seq(0,100,100), cout_cov + 0*seq(0,100,100))


