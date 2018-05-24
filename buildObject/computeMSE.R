computeMSE <- function(cnolist = cnolist, model = model, ode_pars = ode_pars, mseThresh = mseThresh, simData = simData){
  
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
  # sim_data = plotLBodeFitness(cnolist = cnolist, model = model, ode_parameters = ode_pars)
  sim_data = simData
  
  simSplines <- list()
  
  for(i in 1:ncol(sim_data[[1]])){
    
    temp <- list()
    
    for(j in 1:nrow(sim_data[[1]])){
      
      temp[[length(temp)+1]] <- c(sim_data[[1]][j, i], sim_data[[2]][j, i], sim_data[[3]][j, i])
      
    }
    
    simSplines[[length(simSplines)+1]] <- temp
    
  }
  
  names(simSplines) <- cnolist$namesSignals
  
  ##
  mse <- matrix(data = , nrow = nrow(cnolist$valueSignals$`0`), ncol = ncol(cnolist$valueSignals$`0`))
  for(j in 1:length(simSplines)){
    
    for(i in 1:length(simSplines[[j]])){
      
      ss <- mean(c((cnoSplines[[j]][[i]][1]-simSplines[[j]][[i]][1])^2, 
                   (cnoSplines[[j]][[i]][2]-simSplines[[j]][[i]][2])^2,
                   (cnoSplines[[j]][[i]][3]-simSplines[[j]][[i]][3])^2))
      
      mse[i, j] <- ss
      
    }
    
  }
  
  indeces <- list()
  
  for(i in 1:nrow(mse)){
    
    for(j in 1:ncol(mse)){
      
      if(mse[i, j] > mseThresh){
        
        indeces[[length(indeces)+1]] <- c(j, i)
        
      }
      
    }
    
  }
  
  return(indeces)
  
}