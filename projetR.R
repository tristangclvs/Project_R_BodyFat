#===============================================
# Installation des librairies nécessaires
#===============================================

install.packages("PCAmixdata")
install.packages("FactoMineR")
install.packages("factoextra")
install.packages("Hmisc")
library(corrplot)
library(PCAmixdata)
library(FactoMineR)
library(factoextra)
library(Hmisc)

#===============================================
# Chargement des données 
#===============================================

# Pour charger le dataset donneesProjet:
# double clic sur le fichier doonneesProjet2A.RData dans le dossier sous jacent

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
summary(donneesProjet) # visualisation stat. de base
boxplot(donneesProjet,col = c("yellow"),main = paste("Boxplot"), ylab = "Quantiles")

hist.data.frame(donneesProjet[,2:4])
hist.data.frame(donneesProjet[,6:14])


ACP <- PCA(data.frame(donneesProjet), graph=FALSE)
round(ACP$eig,digit=2)
n = names(donneesProjet)

# Outliers
#------------------------------------

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

# Analyse en Composantes Principales
#------------------------------------

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

# Calcul de la matrice de corrélation
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

# On définie une étude linéaire multiple pour expliquer le modèle de Pct.BF en fonctions des données de donnesProjet
resPct <- lm(Pct.BF~.,data=donneesProjet)
summary(resPct)
# Adjusted R-squared:  0.7368 

step(resPct)

plot(resPct$fitted,resPct$residuals)
plot(resPct$fitted,donneesProjet$Pct.BF)
abline(a=0,b=1)
shapiro.test(resPct$residuals) # p value 0.09 => résidus normalement distribués

# Etude plus fine du modèle
#------------------------------------
# On enlève les variables selon le critère d'AIC - méthode drop1 ou step descendant
dataTest = donneesProjet[,c(1,2,4,5,7,8,9,13,14)]
resTest = lm(Pct.BF~.,data=dataTest)
summary(resTest)
# Adjusted R-squared:  0.7385 

# resTest2 = lm(Pct.BF~Hip+Forearm+Thigh+Neck+Height+Age+Wrist+Abdomen,data=donneesProjet)
# summary(resTest2) 
# Ces 2 lignes servent à montrer que l'on obtient bien le même résultat

# Avec la méthode add1 ou step ascendant
resTest3 = lm(Pct.BF~Abdomen+Weight+Wrist+Bicep+Age+Thigh,data=donneesProjet)
add1(resTest3,~Age+Weight+Height+Neck+Chest+Abdomen+Hip+Thigh+Knee+Ankle+Bicep+Forearm+Wrist)
summary(resTest3)
# Adjusted R-squared:  0.736 

# --> on obtient un meilleur Rsquared pour resTest <--
# On garde resTest pour la suite de l'étude

plot(resTest$fitted,resTest$residuals) # pas de structure dans les résidus
abline(h=0)
plot(resTest$fitted,dataTest$Pct.BF)
abline(a=0,b=1)
shapiro.test(resTest$residuals) 
# p value 0.13 => résidus distribués normalement

# Résidus de student pour trouver de potentiels valeurs aberrantes
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

# Adjusted R-squared:  0.7745
residus.stud2<-rstudent(resRemoved)
plot(residus.stud2,ylab="res. student.",ylim=c(-3.5,3.5))
abline(h=c(-2,0,2),lty=c(2,1,2))
which(residus.stud2>2)
shapiro.test(resRemoved$residuals)
# p-value = 0.002063 --> Non distribués normalement
# Revenir au modèle précédent <-------------------

############################################################################
# Modèle choisi pour prédire le pourcentage de masse graisseuse:
dataTest = donneesProjet[,c(1,2,4,5,7,8,9,13,14)]
resTest = lm(Pct.BF~.,data=dataTest)
# Ainsi, pour expliquer au mieux cette indice on utilise les variables:
# Age, Height, Neck, Abdomen, Hip, Thigh, Forearm, Wrist
# Avec ces dernières on obtient 73,85 % d'explication de la masse graisseuse
############################################################################
