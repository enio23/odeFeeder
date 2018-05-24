library(CellNOptR)
library(CNORode2017)
library(readr)
library(igraph)

load(file = "optimalPruned.RData")
load(file = "simData.RData")

temp_pars <- opt_pars

source("computeCorrelation.R")
source("computeMSE.R")

cnolist <- CNOlist(data = "prunedMIDAS.csv")
cnolist <- compatCNOlist(object = cnolist)

model <- readSIF(sifFile = "prunedPKN.sif.txt")

# plotModel(model = model, CNOlist = cnolist)

# indeces <- computeCorrelation(cnolist = cnolist, model = model, ode_pars = opt_pars)
indeces <- computeMSE(cnolist = cnolist, model = model, ode_pars = opt_pars, mseThresh = 0.05, simData = simData)
# indeces <- list()
# indeces[[length(indeces)+1]] <- c(14, 13)

signedInteractions_1 <- read_delim("signedInteractions-1.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
signedInteractions_2 <- read_delim("signedInteractions-6.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

phonemes <- unique(rbind(signedInteractions_1, signedInteractions_2))

prunedPKN_sif <- read_delim("prunedPKN.sif.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
pknAddition <- as.matrix(prunedPKN_sif)

inhibitors <- cnolist$namesInhibitors

inhID <- c("AKT1_HUMAN", "ATM_HUMAN", "ATR_HUMAN", "CDK1_HUMAN", "CDK2_HUMAN", "CDK5_HUMAN", "CHK1_HUMAN", "GSK3B_HUMAN", "MP2K1_HUMAN", "MK01_HUMAN",
           "MK14_HUMAN", "MK03_HUMAN", "MK08_HUMAN", "MTOR_HUMAN", "PDPK1_HUMAN", "PK3CA_HUMAN", "KAPCA_HUMAN", "KPCA_HUMAN", "PRKDC_HUMAN", "KS6B1_HUMAN")

sitesPKN <- cnolist$namesSignals

sites <- c("ATM_HUMAN.S.1981", "CHK1_HUMAN.S.280", "CHK1_HUMAN.S.317", "H2AX_HUMAN.S.140", "H2AX_HUMAN.T.137", "HSPB1_HUMAN.S.82",
           "JUN_HUMAN.S.63", "KAT7_HUMAN.T.85", "MK01_HUMAN.T.185", "MK03_HUMAN.T.202", "NBN_HUMAN.S.343", "NCBP1_HUMAN.S.22", "PARN_HUMAN.S.557",
           "PPM1G_HUMAN.S.183", "PYR1_HUMAN.S.1859", "RAD50_HUMAN.S.635", "RPTOR_HUMAN.S.859", "SMC1A_HUMAN.S.957", "STMN1_HUMAN.S.38")

presentSpecies <- c(sites, inhID)

pknG <- graph_from_data_frame(d = prunedPKN_sif[, c(1, 3)])
adj1 <- get.adjacency(graph = pknG)
phonemesG <- graph_from_data_frame(d = phonemes[, c(1, 3)])
adj2 <- get.adjacency(graph = phonemesG)

phonemesPaths <- list()
cnoPaths <- list()
bindedReactions <- list()
lenList <- c()
sP_all <- list()
for(i in 1:length(indeces)){
  
  idx <- indeces[[i]]
  
  site <- sitesPKN[idx[1]]
  target <- inhibitors[which(cnolist$valueInhibitors[idx[2], ]==1)]
  
  len <- c()
  for(j in 1:length(target)){
    
    sP <- get.all.shortest.paths(graph = pknG, from = which(rownames(adj1)==target[j]), to = which(rownames(adj1)==site))
    len <- c(len, length(sP[[1]]))
    
  }
  
  lenList <- c(lenList, sum(len))
  
  siteID <- sites[which(sitesPKN==site)]
  targetID <- inhID[which(inhibitors%in%target)]
  
  paths <- list()
  counterOfSpecies <- list()
  
  for(j in 1:length(targetID)){
    
    sP <- get.all.shortest.paths(graph = phonemesG, from = which(rownames(adj2)==targetID[j]), to = which(rownames(adj2)==siteID[j]))
    
    if(length(sP[[1]]) > 0){
      
      for(k in 1:length(sP[[1]])){
        
        currPath <- sP[[1]][[k]]
        
        sP_all[[length(sP_all)+1]] <- currPath
        
      }
      
    }
    
  }
  
}

for(i in 1:length(sP_all)){
  
  sP <- sP_all[[i]]
  
  toBind <- matrix(data = , nrow = 1, ncol = 3)
  for(j in 1:(length(sP)-1)){
    
    ss <- rownames(adj2)[sP[j]]
    if(ss%in%inhID){ss <- inhibitors[which(inhID==ss)]}
    if(ss%in%sites){ss <-sitesPKN[which(sites==ss)]}
    tt <- rownames(adj2)[sP[j+1]]
    if(tt%in%inhID){tt <- inhibitors[which(inhID==tt)]}
    if(tt%in%sites){tt <-sitesPKN[which(sites==tt)]}
    
    toBind <- rbind(toBind, t(as.matrix(c(ss, 
                                          as.numeric(as.character(phonemes$Sign[intersect(which(phonemes$Source==rownames(adj2)[sP[j]]), which(phonemes$Target==rownames(adj2)[sP[j+1]]))])),
                                          tt))))
  }
  
  if(nrow(toBind) > 2){
    
    toBind <- toBind[-1, ]
    colnames(toBind) <- colnames(pknAddition)
    
  } else {
    
    toBind <- toBind[-1, ]
    toBind <- t(as.matrix(toBind))
    
  }
  
  toBindSpecies <- unique(c(toBind[, 1], toBind[, 3]))
  toBindRem <- toBindSpecies[grepl(pattern = ".", x = toBindSpecies, fixed = TRUE)]
  if(length(toBindRem) > 0){
    
    connectMissing <- matrix(data = , nrow = 1, ncol = 3)
    for(j in 1:length(toBindRem)){
      
      connectMissing <- rbind(connectMissing, t(as.matrix(c(toBind[which(toBind[, 3]==toBindRem[j]), 1], "1", toBind[which(toBind[, 1]==toBindRem[j]), 3]))))
      
    }
    
    connectMissing <- connectMissing[-1, ]
    toBindRemIdx <- unique(c(which(toBind[, 1]%in%toBindRem), which(toBind[, 3]%in%toBindRem)))
    toBindFinal <- rbind(connectMissing, toBind[-toBindRemIdx, ])
    colnames(toBindFinal) <- colnames(pknAddition)
    pknAddition <- unique(rbind(pknAddition, toBindFinal))
    
    bindedReactions[[length(bindedReactions)+1]] <- toBindFinal
  
  } else {
    
    bindedReactions[[length(bindedReactions)+1]] <- toBind
    
  }
  
}

object <- vector(mode = "list", length = 2)
object[[1]] <- unique(bindedReactions)
object[[2]] <- as.matrix(prunedPKN_sif)
names(object) <- c("SetoOfReactions", "PrunedPKN")

save(object, file = "object.RData")
