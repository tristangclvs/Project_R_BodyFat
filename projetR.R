#===============================================
# Installation des librairies nécessaires
#===============================================

require("PCAmixdata")
require("FactoMineR")
require("factoextra")
require("Hmisc")
# install.packages("PCAmixdata")
# install.packages("FactoMineR")
# install.packages("factoextra")
# install.packages("Hmisc")
library(corrplot)
library(PCAmixdata)
library(FactoMineR)
library(factoextra)
library(Hmisc)

#===============================================
# Chargement des données 
#===============================================
donneesProjet = load('donneesProjet2A.RData')

#===============================================
# Informations sur les datas
#===============================================

# 250 individus et 14 variables quantitatives
#------------------------------------
# 1. Pct.BF : Pourcentage de graisse corporelle en utilisant l'équation de Siri, 495/Densité – 450,
# 2. Age : Âge (en années),
# 3. Weight : Poids (en livres – 1 livre = 453,592 grammes),
# 4. Height : Taille (en pouces=inches – 1 pouce = 2,54 centimètres),
# 5. Neck : Circonférence du cou (en cm),
# 6. Chest : Circonférence de la poitrine (en cm),
# 7. Abdomen : Tour de l'abdomen (en cm) "à l'ombilic et au niveau de la crête iliaque",
# 8. Hip : Circonférence de la hanche (en cm),
# 9. Thigh : Circonférence de la cuisse (en cm),
# 10. Knee : Circonférence du genou (en cm),
# 11. Ankle : Circonférence de la cheville (en cm),
# 12. Bicep : Circonférence du biceps étendu (en cm),
# 13. Forearm : Circonférence de l'avant-bras (en cm),
# 14. Wrist : Circonférence du poignet (en cm) "distal à l'apophyse styloïde".
#------------------------------------

# Première visualisation du jeu de données
#------------------------------------
plot(donneesProjet) 

# Statistique descriptive
#------------------------------------
head(donneesProjet) 
summary(donneesProjet) #visualisation stat. de base
boxplot(donneesProjet,col = c("yellow"),main = paste("Boxplot"), ylab = "Quantiles")

hist.data.frame(donneesProjet[,2:4])
hist.data.frame(donneesProjet[,5:14])

###############################################################################################
# Début du code recup
###############################################################################################

ACP <- PCA(data.frame(donneesProjet), graph=FALSE)
round(ACP$eig,digit=2)
n = names(donneesProjet)


###########################
# Outliers
###########################

for (i in (1:length(donneesProjet))){
  boxplot(donneesProjet[,i],xlab=n[i])
  print(n[i])
  val = min(max(donneesProjet[,i]),quantile(donneesProjet[,i],0.75)+ 1.5*(quantile(donneesProjet[,i],0.75)-quantile(donneesProjet[,i],0.25)))
  val2 = max(min(donneesProjet[,i]),quantile(donneesProjet[,i],0.25)- 1.5*(quantile(donneesProjet[,i],0.75)-quantile(donneesProjet[,i],0.25)))
  print(which (donneesProjet[,i]> val))
  print(which (donneesProjet[,i]<val2))
  
}
# On retient l'outlier 40 - 43 - 214 - 224
# Car ils reviennent plusieurs fois

plot(donneesProjet$Ankle,donneesProjet$Height)
which(donneesProjet$Ankle>32)

res<-PCAmix(donneesProjet)
val3 = data.frame(res$quanti$contrib.pct)
val4 = data.frame(res$ind$coord)
plot(val4)

# La 5e composante a un problème, on a 2 observations qui posent un problème
# Impactent-elles le modèle ?
# Regardons l'importance de la variable PCT.bf sur la 5e dimension
# cos2 := 6.81 e-04
# C'est très faible, donc peu d'importance de ces outliers, mais ne pas les oublier
# Ce sont les individus 31 et 84

val5 = data.frame(res$quanti$cos2)

round(res$eig,digits=2)
barplot(res$eig[,1],main="Eigenvalues",names.arg=1:nrow(res$eig))
abline(h=1,col=2,lwd=2)

plot(res,axes=c(1,2),choice="ind") # on retrouve ici le graphique des individus (plan 1-2)
plot(res,axes=c(1,2),choice="cor") # on retrouve ici le cercle des corrélations
plot(res,axes=c(1,2),choice="sqload") # on retrouve ici le graphique des "square loadings" (plan 1-2)

res$quanti$cos2

cor(donneesProjet)  # calcul de la matrice des correlations lineaires des donnees
corrplot(cor(donneesProjet),is.corr=TRUE, method="num",type="lower")

# On définie une étude linéaire pour expliquer le modèle de Pct.BF en fonctions des données de donnesProjet
resPct <- lm(Pct.BF~.,data=donneesProjet)
summary(resPct)
# Adjusted R-squared:  0.7368 

step(resPct)

plot(resPct$fitted,resPct$residuals)
plot(resPct$fitted,donneesProjet$Pct.BF)
abline(a=0,b=1)
shapiro.test(resPct$residuals) # p value 0.09 => rejet H0 ou pas ?


# Étude plus fine du modèle
dataTest = donneesProjet[,c(1,2,4,5,7,8,9,13,14)]
resTest = lm(Pct.BF~.,data=dataTest)
resTest2 = lm(Pct.BF~Hip+Forearm+Thigh+Neck+Height+Age+Wrist+Abdomen,data=donneesProjet)
summary(resTest)
summary(resTest2)
# C'est la même chose


step(resTest) # Valide notre choix de variables

plot(resTest$fitted,resTest$residuals)
abline(h=0)
plot(resTest$fitted,dataTest$Pct.BF)
abline(a=0,b=1)
shapiro.test(resPct$residuals) # p value 0.09 => rejet H0 ou pas ?

#
residus.stud<-rstudent(resTest)
plot(residus.stud,ylab="res. student.",ylim=c(-3.5,3.5))
abline(h=c(-2,0,2),lty=c(2,1,2))
which (residus.stud>2)
which(residus.stud<(-2))

# Outliers
# On enlève les observations suivantes: 79 80 133 205 169 202 222 223 229 236
dataRemoved = dataTest[c(-79,-80,-133,-205,-169,-202,-222,-223,-229,-236),]
resRemoved = lm(Pct.BF~.,data=dataRemoved)
summary(resRemoved)
# Adjusted R-squared:  0.7745 --> notre meilleure
residus.stud2<-rstudent(resRemoved)
plot(residus.stud2,ylab="res. student.",ylim=c(-3.5,3.5))
abline(h=c(-2,0,2),lty=c(2,1,2))
which(residus.stud2>2)
shapiro.test(resRemoved$residuals)
# p-value = 0.002063 --> Non distribués normalement
# Revenir au modèle précédent -------------------

###############################################################################################
# fin du code recup
###############################################################################################

# Analyse en Composantes Principales
#------------------------------------
# Eigen values
ACP <- PCA(data.frame(donneesProjet), graph=FALSE)
round(ACP$eig,digit=2)

# Proportion cumulative des Eigen values
res<-PCAmix(donneesProjet)
round(res$eig,digits=2)

# Barplot des 70% d'explication selon les dimensions => retenir dim 1 et 2 (somme = 72.48)
barplot(res$eig[,1],main="Eigenvalues",names.arg=1:nrow(res$eig))
abline(h=1,col=2,lwd=2) # critère de Kaiser  

plot(res,axes=c(1,2),choice="ind")    # on retrouve ici le graphique des individus (plan 1-2)
plot(res,axes=c(1,2),choice="cor")    # on retrouve ici le cercle des corrélations
plot(res,axes=c(1,2),choice="sqload") # on retrouve ici le graphique des "square loadings" (plan 1-2)

# Pct.BF est fortement expliquée par les dimensions 1 et 2 / Age expliqué par les dimensions 2 et 3
# Weight expliqué par la dimension 1 / Height dimension 2 / Neck dimension 1 / Chest dimension 1
# Abdomen dimension 1 / Hip dimension 1 / Thigh dimension 1 / Knee dimension 1 / Ankle faiblement expliqué par dim 1
# Bicep dimension 1 / Forearm dimension 1 / Wrist dimension 1
res$quanti$cos2 # cos2 des dimensions représentés dans la table de chaleur suivante

# Analyse de la corrélation
#------------------------------------
# calcul de la matrice des correlations lineaires des donnees
matCor = cor(donneesProjet)  

# forte corrélation => coeff. > 0.7
matCor2 = matCor
for(i in 1:nrow(matCor)){
  for(j in 1:ncol(matCor)){
    if(matCor2[i,j]<0.7){
      matCor2[i,j]=0
    }
  }
}
corrplot(matCor,is.corr=TRUE, method="shade", type="lower")
corrplot(matCor2,is.corr=TRUE, method="shade", type="lower") # Visualisation des données corrélées exclusivement

#===============================================
# Régression linéaire multiple
#===============================================
# Estimation modèle de régression linéaire multiple
#------------------------------------
resPct <- lm(Pct.BF~.,data=donneesProjet)
summary(resPct)
plot(resPct$fitted,resPct$residuals)
abline(h=0)
plot(resPct$fitted,donneesProjet$Pct.BF)
abline(0,1)
# Test de normalité des résidus
#------------------------------------
shapiro.test(resPct$residuals)
#p-value = 0.094 => non rejet de H0 : population non distribuée normalement

# Etude plus fine du modèle
#------------------------------------

# #Intervalle de confiance
# test = resPct$fitted
# grille.x<-data.frame(x=seq(from=min(test),to=max(test),length.out=5000))
# ICpred<-predict(resPct,new=grille.x,interval="pred",level=0.95)
# ICmoy<-predict(resPct,new=grille.x,interval="conf",level=0.95)
# 
# # Trace de ces intervalles de prevision et de confiance sur le nuage de points
# plot(x,y)
# abline(res)
# matlines(grille.x,cbind(ICpred[,],ICmoy[,-1]),lty=c(1,2,2,3,3),col=c(1,2,2,3,3))
# legend("bottomright",lty=c(2,3),col=c(2,3),c("IC prevision","IC moy"))

# Avec la méthode drop1 ou step descendant
dataTest = donneesProjet[,c(1,2,4,5,7,8,9,13,14)] # viennent du step resPCT
resTest = lm(Pct.BF~Hip+Forearm+Thigh+Neck+Height+Age+Wrist+Abdomen,data=donneesProjet)
summary(resTest)

# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5.04038    8.35881   0.603 0.547074    
# Hip         -0.19488    0.12984  -1.501 0.134689    
# Forearm      0.29550    0.19166   1.542 0.124440    
# Thigh        0.22387    0.12900   1.735 0.083943 .  
# Neck        -0.45133    0.21774  -2.073 0.039252 *  
# Height      -0.26807    0.12612  -2.125 0.034567 *  
# Age          0.07258    0.03030   2.396 0.017361 *  
# Wrist       -1.73072    0.49360  -3.506 0.000542 ***
# Abdomen      0.82271    0.06880  11.958  < 2e-16 ***
# Adjusted R-squared:  0.7385 

# Avec la méthode add1 ou step ascendant
resTest3 = lm(Pct.BF~Abdomen+Weight+Wrist+Bicep+Age+Thigh,data=donneesProjet)
add1(resTest3,~Age+Weight+Height+Neck+Chest+Abdomen+Hip+Thigh+Knee+Ankle+Bicep+Forearm+Wrist)
summary(resTest3)

# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -32.57374    8.36270  -3.895 0.000127 ***
# Abdomen       0.88408    0.06984  12.658  < 2e-16 ***
# Weight       -0.10683    0.03462  -3.086 0.002263 ** 
# Wrist        -1.76369    0.48630  -3.627 0.000350 ***
# Bicep         0.24326    0.15615   1.558 0.120584    
# Age           0.06175    0.03070   2.012 0.045374 *  
# Thigh         0.17838    0.12100   1.474 0.141722    
# Adjusted R-squared:  0.736 



resTest4 = lm(Pct.BF~.,data=dataTest)
summary(resTest4)
plot(resTest4$fitted,resTest4$residuals)
plot(resTest4$fitted,dataTest$Pct.BF)
abline(a=0,b=1)
shapiro.test(resTest4$residuals) # p value 0.038 normalement distribué 


plot(resTest$fitted,resTest$residuals)
plot(resTest$fitted,dataTest$Pct.BF)
abline(a=0,b=1)
shapiro.test(resTest$residuals) # p value 0.13 => non rejet de H0 
