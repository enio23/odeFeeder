library(doParallel)
library(CellNOptR)
library(CNORode2017)
library(MEIGOR)
library(readr)

argsJob= commandArgs(trailingOnly=TRUE)
repIndex <- as.numeric(argsJob[1])

cnolist <- CNOlist(data = "prunedMIDAS.csv", verbose = TRUE)
modelSupport <- readSIF(sifFile = "prunedPKN.sif.txt")

load(file = "object.RData")

write.table(x = object$PrunedPKN, file = "pkn.sif.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

sif <- object$PrunedPKN

load(file = "optimalPruned.RData")

optRes <- opt_pars$ssm_results$fbest

library(igraph)

improvements <- c()

####################################
####################################

print(paste0("Step ----- ", repIndex, "/", length(object$SetoOfReactions)))

reactionSet <- object$SetoOfReactions[[repIndex]]

currSIF <- unique(rbind(sif, reactionSet))
rownames(currSIF) <- as.character(1:nrow(currSIF))

write.table(x = currSIF, file = "pkn.sif.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

model <- readSIF(sifFile = "pkn.sif.txt")

ode_parameters_temp=createLBodeContPars(model, LB_n = 1, LB_k = 0,
                                        LB_tau = 0, UB_n = 4, UB_k = 10, UB_tau = 10, default_n = 4,
                                        default_k = 0.5, default_tau = 0.01, opt_n = FALSE, opt_k = TRUE,
                                        opt_tau = TRUE, random = TRUE)

lb <- ode_parameters_temp$LB
ub <- ode_parameters_temp$UB

gg <- graph_from_data_frame(d = currSIF[, c(1, 3)])
adj <- get.adjacency(graph = gg)
# cues <- as.character(reactionSet[1, 1])
cues <- colnames(cnolist@cues)

reacSetSpecies <- unique(c(reactionSet[, 1], reactionSet[, 3]))
measSites <- reacSetSpecies[which(reacSetSpecies%in%colnames(cnolist@signals$`0`))]
allP <- list()
for(ii in 1:length(measSites)){
  
  for(jj in 1:length(cues)){
    
    path <- all_simple_paths(graph = gg, from = which(rownames(adj)==cues[jj]), to = which(rownames(adj)==measSites[ii]))
    
    if(length(path) > 0){
      
      for(kk in 1:length(path)){
        
        allP[[length(allP)+1]] <- path[[kk]]
        
      }
      
    }
    
  }
  
}

edges <- c()
for(ii in 1:length(allP)){
  
  for(jj in 1:(length(allP[[ii]])-1)){
    
    edges <- c(edges, which(grepl(pattern = paste0(rownames(adj)[allP[[ii]][jj]], "=", rownames(adj)[allP[[ii]][jj+1]]), x = model$reacID)))
    
  }
  
}

edges <- unique(edges)
edges <- model$reacID[edges]

setSIF <- matrix(data = , nrow = length(edges), ncol = 3)
setSIF[, 2] <- "1"
for(ii in 1:length(edges)){
  
  setSIF[ii, 1] <- gsub(pattern = "!", replacement = "", x = strsplit(x = edges[ii], split = "=")[[1]][1])
  setSIF[ii, 3] <- gsub(pattern = "!", replacement = "", x = strsplit(x = edges[ii], split = "=")[[1]][2])
  if(grepl(pattern = "!", x = edges[ii], fixed = TRUE)){setSIF[ii, 2] <- "-1"}
  
}

write.table(x = setSIF, file = paste0("setSIF_", repIndex, ".sif.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

setModel <- readSIF(sifFile = paste0("setSIF_", repIndex, ".sif.txt"))

setPars <- createLBodeContPars(setModel, LB_n = 1, LB_k = 0,
                               LB_tau = 0, UB_n = 4, UB_k = 10, UB_tau = 10, default_n = 4,
                               default_k = 0.5, default_tau = 0.01, opt_n = FALSE, opt_k = TRUE,
                               opt_tau = TRUE, random = TRUE)

file.remove(paste0("setSIF_", repIndex, ".sif.txt"))

for(ii in 1:length(ode_parameters_temp$parNames)){
  
  if((ode_parameters_temp$parNames[ii]%in%setPars$parNames)==FALSE){
    
    if(length(which(opt_pars$parNames==ode_parameters_temp$parNames[ii])) > 0){
      ode_parameters_temp$parValues[ii] <- opt_pars$parValues[which(opt_pars$parNames==ode_parameters_temp$parNames[ii])]
      ode_parameters_temp$LB[ii] <- opt_pars$parValues[which(opt_pars$parNames==ode_parameters_temp$parNames[ii])]
      ode_parameters_temp$UB[ii] <- opt_pars$parValues[which(opt_pars$parNames==ode_parameters_temp$parNames[ii])]
    }
    
  }
  
}

if(length(which(opt_pars$parNames==ode_parameters_temp$parNames[ii])) > 0){
  
  ## Parameter Optimization
  # essm
  paramsSSm=defaultParametersSSm()
  paramsSSm$local_solver = "DHC"
  paramsSSm$maxtime = 300*length(edges);
  paramsSSm$maxeval = Inf;
  paramsSSm$atol=1e-6;
  paramsSSm$reltol=1e-6;
  paramsSSm$nan_fac=1000;
  paramsSSm$dim_refset=30;
  paramsSSm$n_diverse=1000;
  paramsSSm$maxStepSize=Inf;
  paramsSSm$maxNumSteps=10000;
  transferFun=4;
  paramsSSm$transfer_function = transferFun;
  
  paramsSSm$lambda_tau=0.1
  paramsSSm$lambda_k=rep(0.01, length(ode_parameters_temp$index_k))
  paramsSSm$lambda_k[which(model$reacID%in%modelSupport$reacID==FALSE)]=0.02
  paramsSSm$bootstrap=F
  paramsSSm$SSpenalty_fac=0.01
  paramsSSm$SScontrolPenalty_fac=0.01
  
  # run the optimisation algorithm
  print(paste0("Fitting ", length(edges)*2, " parameters..."))
  
  opt_pars_temp=parEstimationLBode(cnolist,model, method="essm", ode_parameters=ode_parameters_temp, paramsSSm=paramsSSm)
  
  if(opt_pars_temp$ssm_results$fbest < optRes){
    
    print(paste0("Improvement on step ", repIndex, "/", length(object$SetoOfReactions)))
    improvements <- c(improvements, repIndex)
    write.table(x = improvements, file = paste0("improvements_", repIndex, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    
    save(opt_pars_temp, file = paste0("improvedParameters_", repIndex, ".RData"))
    
  }
  
}