require("PCAmixdata")
require("FactoMineR")
require("factoextra")
# install.packages("PCAmixdata")
# install.packages("FactoMineR")
# install.packages("factoextra")
library(corrplot)
library(PCAmixdata)
library(FactoMineR)
library(factoextra)
donnesProjet = load('donneesProjet2A.RData')

#250 individus et 14 variables quantitatives
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


# Regarder chaque valeur du tableau avant de faire l'ACP
# --> analyse univariée


##### correlation analysis ##### 

# Eigen values
ACP <- PCA(data.frame(donneesProjet), graph=FALSE)
round(ACP$eig,digit=2)

# Proportion cumulative des Eigen values
res<-PCAmix(donneesProjet)
round(res$eig,digits=2)

# Barplot des 70% d'explication selon les dimensions => retenir dim 1 et 2 (somme = 72.48) + critère de Kaiser non ?
barplot(res$eig[,1],main="Eigenvalues",names.arg=1:nrow(res$eig))
abline(h=1,col=2,lwd=2)

plot(res,axes=c(1,2),choice="ind")    # on retrouve ici le graphique des individus (plan 1-2)
plot(res,axes=c(1,2),choice="cor")    # on retrouve ici le cercle des corrélations
plot(res,axes=c(1,2),choice="sqload") # on retrouve ici le graphique des "square loadings" (plan 1-2)

# Pct.BF est fortement expliquée par les dimensions 1 et 2 / Age expliqué par les dimensions 2 et 3
# Weight expliqué par la dimension 1 / Height dimension 2 / Neck dimension 1 / Chest dimension 1
# Abdomen dimension 1 / Hip dimension 1 / Thigh dimension 1 / Knee dimension 1 / Ankle faiblement expliqué par dim 1
# Bicep dimension 1 / Forearm dimension 1 / Wrist dimension 1
res$quanti$cos2 # cos2 des dimensions représentés dans la table de chaleur suivante

# Analyse de la corrélation avec table de chaleur
matCor = cor(donneesProjet)  # calcul de la matrice des correlations lineaires des donnees

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
corrplot(matCor2,is.corr=TRUE, method="shade", type="lower") #Meilleure visualisation des données correlées 
#### end correlation analysis ####

#Analyse des résidus du data frame

resPct <- lm(Pct.BF~.,data=donneesProjet)
summary(resPct)
plot(resPct$fitted,resPct$residuals)
plot(resPct$fitted,donneesProjet$Pct.BF)
abline(a=0,b=1)
shapiro.test(resPct$residuals) # p value 0.09 => non rejet de H0 

##### étude plus fine du modèle

#Intervalle de confiance
test = resPct$fitted
grille.x<-data.frame(x=seq(from=min(test),to=max(test),length.out=5000))

ICpred<-predict(resPct,new=grille.x,interval="pred",level=0.95)
ICmoy<-predict(res,new=grille.x,interval="conf",level=0.95)

# Trace de ces intervalles de prevision et de confiance sur le nuage de points
plot(x,y)
abline(res)
matlines(grille.x,cbind(ICpred[,],ICmoy[,-1]),lty=c(1,2,2,3,3),col=c(1,2,2,3,3))
legend("bottomright",lty=c(2,3),col=c(2,3),c("IC prevision","IC moy"))


# Étude plus fine du modèle
dataTest = donneesProjet[,c(1,2,4,5,7,8,9,13,14)]
resTest = lm(Pct.BF~.,data=dataTest)
resTest2 = lm(Pct.BF~Hip+Forearm+Thigh+Neck+Height+Age+Wrist+Abdomen,data=donneesProjet)
summary(resTest)
summary(resTest2)
# C'est la même chose
step(resTest) # Valide notre choix de variables

plot(resTest$fitted,resTest$residuals)
plot(resTest$fitted,dataTest$Pct.BF)
abline(a=0,b=1)
shapiro.test(resTest$residuals) # p value 0.09 => rejet H0 ou pas ?
