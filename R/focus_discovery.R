ctbab = function(
  y, # response variable
  Cm, # C matrix: WIHZ for the calculation of eigenvalues, then critical value
  Tm, # T vector: all gi used for test statistics calculations
  upnode, # character vecor of(v1, .. ,vm) where t1< ...<tm, 
  #The starting point is the set of all covariates of the upnode (in desecending order), 
  # determines the upper bound of the eigenvalues
  level, # level of interest
  lownode, # character vector of (vm, ..., v1)
  #starting point is the empty set. 
  #determines the lower bound of the eigenvalues
  tmin,  ctrue, lf, ls,
  # tmin, ctrue, lf and ls are used as indicators to avoid their repeated calculations. 
  #For the left branches, "tmin", "ctrue" and "ls" remain the same as its parent node
  #For the right branches, "lf" remains the same as its parent node
  alpha,
  count = 0, ## used for counting the number of iterations of BAB
  maxIt = 0 ){## the maximal number of itereation of each BAB
  set = c(upnode[1:(level-length(lownode))], lownode) ## the worst set at "level", where the full space does not reject its worst case
  if(missing(tmin) || missing(ctrue)){
    tmin = sum(Tm[set])
    ctrue = criticalvalue(round(eigen(tcrossprod(Cm[,set,drop=FALSE]),symmetric = TRUE,only.values = TRUE)$values,8), alpha =alpha)
  }
  
  if(tmin < ctrue){
    rid = "not reject"
  } else{
    if(missing(lf)){
      lf = round(eigen(tcrossprod(Cm[,upnode,drop=FALSE]),symmetric = TRUE,only.values = TRUE)$values,8)
    }
    if(missing(ls)){
      ls = round(eigen(tcrossprod(Cm[,lownode,drop=FALSE]),symmetric = TRUE,only.values = TRUE)$values,8)
    }
    cmax = criticalvalue(getL(lf,ls, level), alpha = alpha)
    
    if(tmin < cmax) rid = "unsure" else rid="reject"
  }
  
  if(count > maxIt) rid="UNSURE" 
  ## CAPITAL UNSURE indicates the chosen maxIt is not enough to get sure outcome.
  
  if(rid == "unsure" && maxIt > 0){
    res1 = ctbab(y=y,Tm=Tm, Cm=Cm, level=level, upnode = upnode[1:(length(upnode)-1)],lownode=lownode,tmin =tmin, ctrue = ctrue, ls=ls,alpha = alpha, count = count+1, maxIt = maxIt)
    rid = res1$res # go with f-, s
    count <- res1$count
    if(rid =="reject"){
      res2 = ctbab(y=y,Tm=Tm, Cm=Cm, level=level, upnode=upnode, lownode = c(lownode, upnode[length(upnode)]), lf = lf, alpha = alpha, count = count+1, maxIt = maxIt) ## go with f, s+
      rid = res2$res 
      count = res2$count
    }
  }
  
  return(list(res = rid, count = count))
  
}




## discov(): shortcut of closed testing to calculate the true discovereis

discov = function(response, alternative, null, data,  maxit = 0, alpha ){
  sqrW = sqrt(mean(response)*(1-mean(response))) 
  designm = as.matrix(data)
  X = as.matrix(designm[,alternative])
  colnames(X) = alternative
  #
  if(missing(null)|| is.null(null)){
    WIHX = sqrW *(sweep(X,2,colMeans(X))) ## W^{1/2}*(I-H)Xs
  } else{
    I = diag(length(response))
    Z = as.matrix(cbind(1,designm[,null]))
    H = Z%*%solve(t(Z)%*%Z)%*%t(Z) ## Z(Z^tZ)^-1Z^t
    WIHX = sqrW *(I-H)%*%X## W^{1/2}*(I-H)Xs
  }
  
  WIHX = scale(WIHX, center = FALSE,scale = TRUE )/sqrt(length(response)-1) ## standardize the WIHX to make sure iw = 1
  IHX = WIHX/sqrW  ##(I-H)X
  
  gi = unlist(sapply(alternative, function(i) sum(colSums(response*IHX[,i,drop=FALSE])^2) ), use.names=TRUE)
  test = sum(gi) 
  lamt = round(eigen(tcrossprod(WIHX),symmetric = TRUE,only.values = TRUE)$values,8)
  crt = criticalvalue(lamt,alpha = alpha)
  
  if(test<crt){
    return(list(discovery = 0,res = "not reject", iterations = 0) ) ## topnode not reject == 0 discovery
  } else{
    layers = names(sort(gi))
    # the alternative model adjusted for null AND sort them in increasing order
    lel = length(layers)
    if(lel==1){
      return(list(discovery = 1,res = "reject", iterations = 0) ) ## topnode not reject == 0 discovery
    } else{
      #starting values
      outcome = list(res = "reject", count = 0) ## at least one discovery because top node is rejected
      iter = 0 ## to store the total number of iterations
      startlevel = lel-1
      while(outcome$res=="reject" && startlevel>0){
        outcome = ctbab(y = response, Cm = WIHX, Tm = gi, upnode = layers, level = startlevel, lownode = NULL,alpha = alpha, maxIt = maxit)
        iter = iter + outcome$count
        startlevel = startlevel-1
      }
      
      if(outcome$res=="reject"){
        return(list(discovery = length(layers)-startlevel,res = outcome$res, iterations = iter) )
      } else{
        return(list(discovery = length(layers)-startlevel-1,res = outcome$res, iterations = iter) )
      }
    }
  }
}



newFocus = function(response, fsets, null, data, maxit = 0, alpha = 0.05 ,adj = 0){
  # results for given focus sets 
  flen = length(fsets)
  alpha = alpha/(flen-adj)
  fmat = matrix(0, flen, 3)
  for(f1 in 1:flen){
    fres = discov(response = response, alternative= fsets[[f1]], null= null, data=data,  maxit = maxit, alpha = alpha )
    fmat[f1,] = c(fres$discovery, fres$res, fres$iterations)
  }
  
  focus_res = as.data.frame(cbind(names(fsets),  sapply(fsets, length), fmat) )
  colnames(focus_res) = c("GO_ID", "size", "td", "rejecttion_ID", "iterations")
  
  
  return(list(focus = focus_res, fsets = fsets))
}


choosepath = function(startingindex = 1, fsets, lowdv){
  if(all(lowdv==0)){
    path = NULL
  } else{
    dvix = order(lowdv, decreasing = TRUE)
    path = dvix[startingindex]
    newdvix = dvix[-startingindex]
    for(k in seq_along(newdvix)){
      cap = sapply(path,function(x) length(intersect(fsets[[newdvix[k]]], fsets[[x]]) )==0 )
      add = all(cap)
      if(add) path = c(path,newdvix[k])
    }
  }
  
  return(path)
}



pick = function(focus_obj, setofinterest){

  dvf = as.numeric(focus_obj$focus$td)
  flsets = focus_obj$fsets
  vk_set = sapply(flsets, function(i) length(setdiff(i,setofinterest)) )
  
  ## calculate d_i -(V_i - setofinterest) for all focus sets but afterwards, we only focus sets with positive discoveries to save computation time
  dr =  dvf- vk_set
  onlypositivedr = which(dr>0)
  dr = dr[onlypositivedr]
  flsets = flsets[onlypositivedr]
  dvf = dvf[onlypositivedr]
  
  ## first count the single discovereis to avoid repetition time during the algorithm
  singlefl = intersect(which(dr==1), which(sapply(flsets,length)==1) ) # subset index for which the individual discovery is claimed
  disc = sum(dr[singlefl]) 
  
  ##new setofinterest after removing the above focus sets with single discoveries
  setofinterest = setdiff(setofinterest,unique(unlist(flsets[singlefl])) )
  
  
  ## algorithm 
 while( length(dr)!=0 ){
   vk_set = sapply(flsets, function(i) length(setdiff(i,setofinterest)) )
   dr =  dvf - vk_set
   pdr_index = which(dr>0)
   dr = dr[pdr_index]
   flsets = flsets[pdr_index]
   dvf = dvf[pdr_index]
   
   dr_choose = choosepath(startingindex = 1,fsets = flsets, lowdv = dr)
   disc = disc + sum(dr[dr_choose])
   
   setofinterest = setdiff(setofinterest,unique(unlist(flsets[dr_choose])) )
  }
  
  return(disc)
}


