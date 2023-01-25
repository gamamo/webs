# load packages
library(MCMCglmm)
library(ggplot2)
library(ape)
library(dplyr)
library(MuMIn)
library(plyr)
library(reshape2)
library(rio)


phy=read.tree("AllBirdsHackett1.tre")

##############
list_posterioriH = list()
list_intervalsH  = list()
list_summariesH  = list()
list_solH        = list()
list_DICH        = list()
list_summaryGeralH = list()
list_vcvH          =list()
list_randomH      =list()

for (i in 1:100){
  
  print(i)
  
  # load data metanetwork ####
  data = read.csv(file="birdtraits_spmetrics_phy_v18_01_2023.csv")
  
  ###############################################################################
  # modelling
  
  # transform some variables
  
  data$PC1 = data$PC1 + 3.603434 + 0.001
  data$NichePosition = abs(data$NichePosition)
  
  # scale variables that will be used in the modelling and exclude NAs
  toscale = c("Mass", "Wing.Length" , "Hand.Wing.Index","Beak.Width", "Kipps.Distance",
              "NichePosition", "NicheBreadth","RangeSize_Meters",
              "ED" , "DR","Nwebs")
  
  dataM = as.data.frame(scale(data[,toscale],center = TRUE, scale = TRUE))
  dataM = cbind(data[,c("Scientific","PC1","PC2","Family1","Realm","Habitat","Trophic.Level","PC1_nobet","IUCN")],dataM)
  
  dataM$Habitat = as.factor(dataM$Habitat)
  dataM$Trophic.Level = as.factor(dataM$Trophic.Level)
  dataM$Family1 = as.factor(dataM$Family1)
  dataM$Realm = as.factor(dataM$Realm)
  dataM$dummy = factor(rep(1, length(dataM[,1])))

  
  # Run MCMCglmm models metanetwork ######
  
  #use only one tree
  tree = phy[[i]]
  tree$tip.label = gsub("_"," ",tree$tip.label)
  
  # prepare data #####################################################
  # remove species not in the tree ####
  notintree = setdiff(dataM$Scientific,tree$tip.label)
  notindata = setdiff(tree$tip.label, dataM$Scientific)
  toexclude = c(notindata,notintree)
  
  #prune tree
  treeprune = drop.tip(tree,toexclude)
  
  #exclude species notintree from the dataset
  dataM = dataM[-which(dataM$Scientific %in% notintree),]
  
  # prepare the covariance matrix
  Ainv<-inverseA(treeprune)$Ainv
  
  # model wit biogeographic realm
  prior<-list(R=list(V=1, nu=0.002),
              G=list(G1=list(V=1, nu=0.002),
                     G2=list(V=1, nu=0.001)))
  
  MCmodmetahab = MCMCglmm(PC1 ~  RangeSize_Meters + Wing.Length + Mass + NichePosition + 
                            Hand.Wing.Index + Beak.Width + ED + NicheBreadth + DR, 
                          random=~Scientific + Realm,
                          prior = prior,
                          ginverse=list(Scientific=Ainv),
                          nitt=100000, burnin = 25000 ,thin = 50,
                          data=na.omit(dataM))
  
  list_posterioriH[[i]] = posterior.mode(MCmodmetahab$Sol)
  list_intervalsH[[i]] = HPDinterval(MCmodmetahab$Sol)
  list_summariesH[[i]]  = summary(MCmodmetahab$Sol)[[2]]
  list_solH[[i]]        = MCmodmetahab$Sol
  list_DICH[[i]]        = DIC(MCmodmetahab)
  list_summaryGeralH[[i]] = summary(MCmodmetahab)[[5]]
  list_vcvH[[i]]        = MCmodmetahab$VCV
  list_randomH[[i]] =    summary(MCmodmetahab)[[6]]
  
  
} # fecha looping

#setwd("C:/Users/gabri/Dropbox/postdocINECOL/dadosCompilados/workingfolder_GM/scripts_MCMC_saturn_7abr22")
export(list_posterioriH,"PosteriorH.xlsx")
export(list_intervalsH,"IntervalsH.xlsx")
export(list_summariesH,"SummariesH.xlsx")
export(list_solH,"SOLH.xlsx")
export(list_DICH,"DICH.xlsx")
export(list_summaryGeralH,"SummaryGeneralH.xlsx")
export(list_vcvH, "VCVH.xlsx")
export(list_randomH, "posterior_randomH.xlsx")

#setwd("C:/Users/gabri/Dropbox/postdocINECOL/dadosCompilados/workingfolder_GM/scripts_MCMC_saturn_7abr22")
save.image("workspace_Modglobal.RData")
