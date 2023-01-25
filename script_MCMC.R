library(MCMCglmm)
library(ggplot2)
library(ape)
library(dplyr)
library(MuMIn)
library(plyr)
library(reshape2)
library(ggridges)
library(rio)


phy=read.tree("AllBirdsHackett1.tre")


# make a looping to run MCMC for all 1000 trees

list_posteriori = list()
list_intervals  = list()
list_summaries  = list()
list_sol        = list()
list_DIC        = list()
list_summaryGeral = list()
list_vcv  = list()
list_random = list()

#for (i in 1:length(phy)){
for (i in 1:100){
  
  print(i)
  # load data metanetwork ####
  data = read.csv(file="birdtraits_spmetrics_phy_v18_01_2023.csv")

  ###############################################################################
  # modelling
  
  # transform some variables
  
  data$PC1 = data$PC1 + 3.603434 + 1
  data$NichePosition = abs(data$NichePosition)
  
  # scale variables that will be used in the modelling and exclude NAs
  toscale = c("Mass", "Wing.Length" , "Hand.Wing.Index","Beak.Width", "Kipps.Distance",
              "NichePosition", "NicheBreadth","RangeSize_Meters",
              "ED" , "DR","Nwebs")
  
  dataM = as.data.frame(scale(data[,toscale],center = TRUE, scale = TRUE))
  dataM = cbind(data[,c("Scientific","PC1","PC2","Family1","Realm","Habitat","Trophic.Level","PC1_nobet")],dataM)
  
  dataM$Habitat = as.factor(dataM$Habitat)
  dataM$Trophic.Level = as.factor(dataM$Trophic.Level)
  dataM$Family1 = as.factor(dataM$Family1)
  dataM$Realm = as.factor(dataM$Realm)
  dataM$dummy = factor(rep(1, length(dataM[,1])))
  
  #dataM = na.exclude(dataM)
  
  summary(dataM$PC1)
  
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
  
  # model without biogeographic realm
  prior<-list(R=list(V=1, nu=0.002),
              G=list(G=list(V=1, nu=0.002)))
  
  MCmodmeta = MCMCglmm(PC1 ~  RangeSize_Meters + Wing.Length + Mass + NichePosition + 
                         Hand.Wing.Index + Beak.Width + ED + NicheBreadth, 
                       random=~Scientific,
                       prior = prior,
                       ginverse=list(Scientific=Ainv),
                       nitt=100000, burnin = 25000 ,thin = 50,
                       data=na.omit(dataM))
  
  list_posteriori[[i]] = posterior.mode(MCmodmeta$Sol)
  list_intervals[[i]] = HPDinterval(MCmodmeta$Sol)
  list_summaries[[i]]  = summary(MCmodmeta$Sol)[[2]]
  list_sol[[i]]        = MCmodmeta$Sol
  list_DIC[[i]]        = DIC(MCmodmeta)
  list_summaryGeral[[i]] = summary(MCmodmeta)[[5]]
  list_vcv[[i]]        = MCmodmeta$VCV
  list_random[[i]] = summary(MCmodmeta)[[6]]


  
} # fecha looping

export(list_posteriori,"Posterior.xlsx")
export(list_intervals,"Intervals.xlsx")
export(list_summaries2,"Summaries.xlsx")
export(list_sol,"SOL.xlsx")
export(list_DIC,"DIC.xlsx")
export(list_summaryGeral2,"SummaryGeneral.xlsx")
export(list_vcv, "VCV.xlsx")

###

