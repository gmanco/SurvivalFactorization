library(Matrix)
library(ROCR)

optRoc.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = x^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}

logSurvival = function(idxu, idxv, cascade, A, S, K, mode) {
  llk = 0
  
  u = cascade[idxu, 1]
  v = cascade[idxv, 1]
  t_u = cascade[idxu, 2]
  t_v = cascade[idxv, 2]
  
  rate = S[u, K] * A[v, K]
  if (rate == 0)
    rate = 1e-40
  delta = (t_u - t_v)
  
  if (mode == 'exp') {
    llk = -rate * delta
  }
  else if (mode == "ray") {
    llk = -1 / 2 * rate * delta ^ 2
  } #Complete with other functions
  return(llk)
}

hazard = function(idxu, idxv, cascade, A, S, K, mode) {
  h = 0
  u = cascade[idxu, 1]
  v = cascade[idxv, 1]
  t_u = cascade[idxu, 2]
  t_v = cascade[idxv, 2]
  
  rate = S[u, K] * A[v, K]
  if (rate == 0)
    rate = 1e-40
  delta = (t_u - t_v)
  
  
  if (mode == 'exp') {
    h = rate
  }
  else if (mode == "ray") {
    h = rate * delta
  } #Complete with other functions
  return(h)
}


classify_cascade <- function(cascade, A, S, mode = 'exp') {
  K = dim(A)[2]
  
  llk = rep(0, K)
  
  num_active_nodes = nrow(cascade)
  if (num_active_nodes > 1) {
    for (u in 2:num_active_nodes) {
      a_potential = rep(0, K)
      a_bad = rep(0, K)
      for (v in 1:(u - 1)) {
        for (k in 1:K) {
          a_potential[k] = a_potential[k] + logSurvival(u, v, cascade, A, S, k, mode)
          a_bad[k] = a_bad[k] + hazard(u, v, cascade, A, S, k, mode)
        }
      }
      for (k in 1:K) {
        llk[k] = llk[k] + a_potential[k] + log2(a_bad[k])
      }
    }
  }
  return(which.max(llk))
}

get_predictions_for <-
  function(train_cascade, cascade_to_predict, A, S, K,factor_unactive = 2) {
    num_users = dim(A)[1]
    
    users_in_train = train_cascade[, 1]
    all_users = which(cascade_to_predict>0)
    active_users_in_test = setdiff(all_users,users_in_train)
    
    unactive_users = setdiff(1:num_users, all_users)
    num_unactives_to_sample = min(length(active_users_in_test)*factor_unactive,length(unactive_users))
    sampled_unactive_users = sample(unactive_users,num_unactives_to_sample,replace=FALSE)
    users_to_predict = c(sampled_unactive_users,active_users_in_test)
    horizon = max(cascade_to_predict)
    
    predictions = matrix(0, length(users_to_predict), 3)
    
    for (i in 1:length(users_to_predict)) {
      u = users_to_predict[i]
      
      # estimated_times = train_cascade[, 2] + 1 / (S[u, K] * A[users_in_train, K])
      estimated_times = 1 / (S[u, K] * A[users_in_train, K])
      estimated_times_gen = 1 / (A[users_in_train, ] %*% S[u, ])
      predictions[i, ] = c(u, min(estimated_times), min(estimated_times_gen))
    }
    
    predictions = as.data.frame(predictions)
    names(predictions) = c('node', 'predictedts','predictedts_gen')
    return(predictions)
  }



gen_predictions <- function(cascades,test_cascades_idxs, A, S, split_perc,min_train_size, max_cascade_size) {
  num_cascades_to_predict = length(test_cascades_idxs)
  ass_class = rep(0, num_cascades_to_predict)
  
  results <- data.frame(
    cascade = integer(),
    node = integer(),
    topic = integer(),
    predictedts = double(),
    predictedts_gen = double(),
    truets = double(),
    horizon = double()
  )

  for (i in test_cascades_idxs) {
    indices = which(cascades[, i] > 0)
    
    if (length(indices) > max_cascade_size)
      next
    
    #    elems = cbind(indices,cascades_to_predict[indices,i])
    #    elems = elems[order(elems[,2],decreasing=FALSE),]
    elems = data.frame(node = indices, timestamp = cascades[indices, i])
    elems = elems[order(elems[,2],decreasing=FALSE),]
    
    train_size = floor(nrow(elems) * split_perc)
    test_size = nrow(elems) - train_size
    
    if (train_size > min_train_size  & test_size > 0) {
      cascade_to_classify = elems[1:train_size, ]
      #      cascade_to_predict = elems[(train_size+1):(train_size+test_size),]
      
      K = classify_cascade(cascade_to_classify, A, S)
      if (K != 0) {
        predictions = get_predictions_for(cascade_to_classify, cascades[, i], A, S, K)
        offset = nrow(results)
        n_predictions = nrow(predictions)
        current_res = matrix(0, n_predictions, 7)
        
        current_res[, 1] = rep(i, n_predictions)
        current_res[, 2] = predictions[, 1]
        current_res[, 3] = rep(K, n_predictions)
        current_res[, 4] = predictions[, 2]
        current_res[, 5] = predictions[, 3]
        current_res[, 6] = cascades[predictions[, 1], i]
        current_res[, 7] = rep(max(cascades[, i]), n_predictions)
        
        results[(offset + 1):(offset + n_predictions), ] = current_res
      }
    }
  }
  return(results)
}
