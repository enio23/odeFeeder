computeCorrelation <- function(cnolist = cnolist, model = model, ode_pars = ode_pars){
  
  ##
  cnoSplines <- list()
  
  for(i in 1:ncol(cnolist$valueSignals[[1]])){
    
    temp <- list()
    
    for(j in 1:nrow(cnolist$valueSignals[[1]])){
      
      vals <- c()
      
      for(k in 1:length(cnolist$timeSignals)){
        
        vals <- c(vals, cnolist$valueSignals[[k]][j, i])
        
      }
      
      temp[[length(temp)+1]] <- vals
      
    }
    
    cnoSplines[[length(cnoSplines)+1]] <- temp
    
  }
  
  names(cnoSplines) <- cnolist$namesSignals
  
  ##
  sim_data <- getLBodeDataSim(cnolist = cnolist, model = model, ode_parameters = ode_pars)
  temp=list();
  
  for(i in 1:length(cnolist$timeSignals)){
    
    temp[[i]]=sim_data[[i]];
    
  }
  
  sim_data=temp;
  
  sim1 <- sim_data[c(1, 2)]
  sim2 <- sim_data[c(2, 3)]
  
  ## cno-1
  simSplines <- list()
  for(i in 1:ncol(sim1[[1]])){
    
    temp <- list()
    
    for(j in 1:nrow(sim1[[1]])){
      
      temp[[length(temp)+1]] <- c(sim1[[1]][j, i], sim1[[2]][j, i])
      
    }
    
    simSplines[[length(simSplines)+1]] <- temp
    
  }
  
  names(simSplines) <- cnolist$namesSignals
  
  library(gtools)
  
  rr <- list()
  
  for(i in 1:length(simSplines)){
    
    temp <- list()
    
    for(j in 1:length(simSplines[[i]])){
      
      var1 <- ts(simSplines[[i]][[j]])
      var2 <- ts(cnoSplines[[i]][[j]][c(1, 2)])
      temp[[length(temp)+1]] <- running(X = var1, Y = var2, fun = cor, width = 2)
      
    }
    
    rr[[length(rr)+1]] <- temp
    
  }
  
  for(i in 1:length(rr)){
    
    for(j in 1:length(rr[[1]])){
      
      idx <- which(is.na(rr[[i]][[j]]))
      
      if(length(idx) > 0){
        
        rr[[i]][[j]][idx] <- 0
        
      }
      
    }
    
  }
  
  corr1 <- rr
  
  ## cno-2
  simSplines <- list()
  for(i in 1:ncol(sim2[[1]])){
    
    temp <- list()
    
    for(j in 1:nrow(sim2[[1]])){
      
      temp[[length(temp)+1]] <- c(sim2[[1]][j, i], sim2[[2]][j, i])
      
    }
    
    simSplines[[length(simSplines)+1]] <- temp
    
  }
  
  names(simSplines) <- cnolist$namesSignals
  
  rr <- list()
  
  for(i in 1:length(simSplines)){
    
    temp <- list()
    
    for(j in 1:length(simSplines[[i]])){
      
      var1 <- ts(simSplines[[i]][[j]])
      var2 <- ts(cnoSplines[[i]][[j]][c(2, 3)])
      temp[[length(temp)+1]] <- running(X = var1, Y = var2, fun = cor, width = 2)
      
    }
    
    rr[[length(rr)+1]] <- temp
    
  }
  
  for(i in 1:length(rr)){
    
    for(j in 1:length(rr[[1]])){
      
      idx <- which(is.na(rr[[i]][[j]]))
      
      if(length(idx) > 0){
        
        rr[[i]][[j]][idx] <- 0
        
      }
      
    }
    
  }
  
  corr2 <- rr
  
  ##
  indeces <- list()
  
  for(i in 1:length(corr1)){
    
    for(j in 1:length(corr1[[i]])){
      
      if((as.character(corr1[[i]][[j]])=="-1") || (as.character(corr2[[i]][[j]])=="-1")){
        
        indeces[[length(indeces)+1]] <- c(i, j)
        
      }
      
    }
    
  }
  
  return(indeces)
  
}