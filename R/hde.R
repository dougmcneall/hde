# Package hde source code

emlist = function(X, Y){
  # Create a list of objects to be emulated
  # X             design matrix with (nrow) instances of (ncol) inputs
  # Y             matrix of outputs, with one row
  #               for each row of X
  
  d = ncol(Y)
  em.list = vector(mode='list', length=d)
  
  for(i in 1:d){
    em.obj = NULL
    em.obj$X = X
    em.obj$y = Y[, i]
    em.list[[i]] = em.obj
  }
  em.list
}

# !N.B. need to check that we're passing the formula in correctly.
km.wrap = function(form, em, ...){
  # wrapper function to fit km model to list em.
  # each list element in em contains a design matrix X 
  # and an output vector y
  out = NA
  fit = try(km(form, design=em$X, response=em$y, control=list(trace=FALSE), ...), silent=TRUE)
  #if(class(fit) == "try-error"){
  #  out = NA
  #}
  #else{
  #  out = fit 
  #}
  out = fit
  out
}

# Functions for dimension reduced emulation
pc.project = function(pca,scores.em,Z.em,scale){
  # project principal components
  num.pc <- dim(scores.em)[2]
  if (scale){
    anom = ((pca$rotation[ ,1:num.pc] %*% t(scores.em))*pca$scale)
    anom.sd = ((pca$rotation[ ,1:num.pc] %*% t(Z.em))*pca$scale)          
  }
  else {
    anom = pca$rotation[ ,1:num.pc] %*% t(scores.em)
    anom.sd = pca$rotation[ ,1:num.pc] %*% t(Z.em)   
  }
  tens = t(anom + pca$center)
  return(list(tens = tens, anom.sd = anom.sd))
}

direct.pred = function(form,X,Y,Xnew,...){
  # Directly applies km in parallel to predict each column of an ensemble
  ens.list = emlist(X=X, Y=Y)
  km.list = mclapply(ens.list,FUN=km.wrap, form = form)
  pred.list = mclapply(km.list,
                       FUN = predict,
                       newdata = as.matrix(Xnew, nrow =1),
                       type = 'UK')
  out = sapply(pred.list, function(x) x$mean )
  out
}

finite.cols = function(Y){
  #find where the columns of a matrix have all finite values
  Y.finite = is.finite(Y)
  fincol = apply(Y.finite, 2, all)
  ix = which(fincol)
  return(list(fincol=fincol,ix=ix)) 
}

excise = function(Y){
  # excise non finite columns from a matrix (Y.all) and 
  # index vector of where columns are all finite (all.finite)
  # will this put things back in the right order
  all.finite = finite.cols(Y)
  Y.ex = Y[ , all.finite$ix]
  return(list(Y.ex = Y.ex, all.finite.ix=all.finite$ix))
}

kmpar.pc = function(form = ~., Y, X, newdata, num.pc, scale=FALSE, center=TRUE, type="UK", ...){
  # Base function for emulation of high dimensional data
  # with PCA and Gaussian Process emulator
  
  Ytrunc = excise(Y)
  pca = prcomp(Ytrunc$Y.ex,scale=scale, center=center)
  
  if(is.matrix(newdata)!= TRUE){
    print('matrixifying newdata')
    newdata = matrix(newdata,nrow=1) 
  }
  
  pred = direct.pred(form=form, X=X, Y=pca$x[ ,1:num.pc], Xnew=newdata)
  
  scores.em = matrix(pred$mean, nrow = nrow(newdata))
  Z.em = as.matrix(pred$sd, nrow = nrow(newdata))
  
  # this projection is in the "excised" state
  proj = pc.project(pca, scores.em, Z.em, scale)
  
  # Now insert the projected data back into the original framework
  tens = matrix(NA, nrow=nrow(newdata), ncol=ncol(Y))
  anom.sd = matrix(NA, nrow=nrow(newdata), ncol=ncol(Y))
  
  anom.sd[,Ytrunc$all.finite.ix] = proj$anom.sd
  tens[ , Ytrunc$all.finite.ix] = proj$tens
  
  return(list(tens=tens, scores.em=scores.em, Z.em=Z.em, anom.sd=anom.sd))
}

