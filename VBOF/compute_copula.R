get_copula <- function(Stats, vola_level = NULL, m = 5, type = 'regular', concatenate = FALSE) {
  
  if (!is.null(vola_level) && !concatenate) {
    if (vola_level > 5 || vola_level < 1 || vola_level != as.integer(vola_level)) {
      stop('vola_level must be an integer from 1 to 5') # which level of vola you pick
    }
  } else if (is.null(vola_level) && !concatenate) {
    stop('vola level is not requested. It must be an integer from 1 to 5') # which level of vola you pick
  }
  
  if (m != as.integer(m) || m < 2) {
    stop('m must be an integer larger than 1') #size of the copula
  }
  
  if (type == 'regular') {
    all_stats = Stats$lStats_large
  } else if (tolower(type) == 'fev') {
    all_stats = Stats$lStats_large_fev
  } else if (tolower(type) == 'mom') {
    all_stats = Stats$lStats_large_mom
  } else {
    stop("type must me (i) an empty string, or (ii) 'fev', or (iii) 'mom' ")
  }
  
  if (!concatenate) {
    stats = all_stats[[vola_level]]
  } else {
    stats = all_stats[[1]]
    for (i in 2:5) {
      stats = cbind(stats, all_stats[[i]])
    }
  }
  
  rets = stats[2, ]
  vols = stats[3, ]
  
  cop = compute_copula(rets, vols, m)
  
  return(cop)
}




compute_copula <- function (rets, vols, m){

  N = length(rets)

  I1 = order(rets)
  I2 = order(vols)
  
  Score1 = rep(0, N)
  Score2 = rep(0, N)

  for (j in 1:m){
    range = ((j-1)*floor(N/m)+1):(j*floor(N/m))
    Score1[I1[range]] = j
    Score2[I2[range]] = j
  }

  Copula = matrix(0, m, m)
  for (k in 1:m){
    for (j in 1:m){
      Copula[k, j] = sum(as.integer(Score1==k)*as.integer(Score2==j))
    }
  }
  Copula = Copula/N
  return(Copula)

}
    
  