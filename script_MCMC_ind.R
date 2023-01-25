library(MCMCglmm)
library(ggplot2)
library(ape)
library(dplyr)
library(MuMIn)
library(plyr)
library(reshape2)
library(rio)

# Run MCMCGlmm at network level ####



list_posterioriN = list()
list_intervalsN  = list()
list_summariesN  = list()
list_solN        = list()
list_DICN        = list()
list_summaryGeralN = list()
list_vcvN      =list()
list_randomN      =list()

for (i in 1:100){
  
  print(i)
  
  # load the data
  indwebs = read.csv("database_indwebs_v19_01_2023.csv")
   
  # remove species not in the tree
  notintree = setdiff(indwebs$Scientific,tree$tip.label)
  notindata = setdiff(tree$tip.label, indwebs$Scientific)
  toexclude = c(notindata,notintree)
  
  #prune tree
  treeprune2 = drop.tip(tree,toexclude)
  
  #exclude species notintree from the dataset
  indwebs = indwebs[-which(indwebs$Scientific %in% notintree),]
  
  # transform some variables
  summary(indwebs$PC1)
  indwebs$PC1 = indwebs$PC1 + 4.947339 + 1
  
  summary(indwebs$NichePosition)
  indwebs$NichePosition = abs(indwebs$NichePosition)
  
  # scale variables and remove NAs
  toscaleind = c("Mass", "Wing.Length", "Hand.Wing.Index","Beak.Width",
                 "NichePosition", "NicheBreadth","RangeSize_Meters",
                 "ED" , "DR","Nwebs")
  
  dataMI = as.data.frame(scale(indwebs[,toscaleind],center = TRUE, scale = TRUE))
  dataMI = cbind(indwebs[,c("Scientific","PC1","PC2","Family1","Realm","Habitat","Trophic.Level","IUCN")],dataMI)
  
  dataMI$Habitat = as.factor(dataMI$Habitat)
  dataMI$Trophic.Level = as.factor(dataMI$Trophic.Level)
  dataMI$Family1 = as.factor(dataMI$Family1)
  dataMI$Realm = as.factor(dataMI$Realm)
  dataMI$dummy = factor(rep(1, length(dataMI[,1])))
  
  dataMI = na.exclude(dataMI)
  
  
  # Run MCMCglmm models  #
  
  # prepare the covariance matrix
  Ainv<-inverseA(treeprune2)$Ainv
  
  # model wit biogeographic realm
  prior<-list(R=list(V=1, nu=0.002),
              G=list(G1=list(V=1, nu=0.002),
                     G2=list(V=1, nu=0.002)))
  
  MCmodmetaid = MCMCglmm(PC1 ~  RangeSize_Meters + Wing.Length + Mass + NichePosition + 
                           Hand.Wing.Index + Beak.Width + ED + NicheBreadth + DR, 
                         random=~Scientific + Realm,
                         prior = prior,
                         ginverse=list(Scientific=Ainv),
                         nitt=100000, burnin = 25000 ,thin = 50,
                         data=na.omit(dataMI))
  
  
  list_posterioriN[[i]] = posterior.mode(MCmodmetaid$Sol)
  list_intervalsN[[i]] = HPDinterval(MCmodmetaid$Sol)
  list_summariesN[[i]]  = summary(MCmodmetaid$Sol)[[2]]
  list_solN[[i]]        = MCmodmetaid$Sol
  list_DICN[[i]]        = DIC(MCmodmetaid)
  list_summaryGeralN[[i]] = summary(MCmodmetaid)[[5]]
  list_vcvN[[i]]        = MCmodmetaid$VCV
  list_randomN[[i]] =    summary(MCmodmetaid)[[6]]
  
} # fecha looping

export(list_posterioriN,"PosteriorN.xlsx")
export(list_intervalsN,"IntervalsN.xlsx")
export(list_summariesN,"SummariesN.xlsx")
export(list_solN,"SOLN.xlsx")
export(list_DICN,"DICN.xlsx")
export(list_summaryGeralN,"SummaryGeneralN.xlsx")
export(list_randomN, "posterior_randomN.xlsx")



