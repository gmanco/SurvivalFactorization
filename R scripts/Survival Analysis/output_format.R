require('igraph')
require('Matrix')


generate_relational_format <- function(cascades,filename){
  mat = summary(cascades)
  mat = mat[order(mat[, 1], mat[, 3], decreasing = FALSE), ]
  
  mat = data.frame(User = mat[, 2],
                   Item = mat[, 1],
                   TimeStamp = mat[, 3])
  
  write.table(
    mat,
    file = filename,
    sep = "\t",
    col.names = T,
    row.names = F
  )
}
  

generate_netrate_format <- function(cascades,g,file_cascades,file_network){
  file.create(file_cascades)
  sz = dim(cascades)
  
  n_users = sz[2]

  for (i in 1:n_users){
    x = sprintf('%d,%d',i-1,i);
    write(x,file_cascades,append = TRUE)
  }
  
  cascade_data = summary(cascades)
  n_entries = nrow(cascade_data)
  
  cascade_data = cascade_data[order(cascade_data[, 1], cascade_data[, 3], decreasing = FALSE), ]
  
  n_cascades = 0;
  p_cascade = -1;
  curr_cascade = ""
  write(curr_cascade,file_cascades,append = TRUE)
  
  for (i in 1:n_entries) {
    
    c_cascade = cascade_data[i,1]
    c_user = cascade_data[i, 2]-1
    c_timestamp = cascade_data[i, 3]
    
    
    if (p_cascade != c_cascade) {
      if (p_cascade != -1){
        write(curr_cascade,file_cascades,append=TRUE)
      }
      curr_cascade = sprintf('%d,%f', c_user, c_timestamp)
      p_cascade = c_cascade
      n_cascades = n_cascades + 1
      
    }  else {
      if (i != 1) {
        curr_cascade = paste(curr_cascade,',',sep="")
      }
      curr_cascade = paste(curr_cascade, sprintf('%d,%f', c_user, c_timestamp), sep='')
    }
  }
  write(curr_cascade,file_cascades,append=TRUE)
  
  
  
 # generating network....
  
  file.create(file_network)
  for (i in 1:n_users){
    x = sprintf('%d,%d',i-1,i);
    write(x,file_network,append = TRUE)
  }
  
  write("",file_network,append = TRUE)
  E = get.adjacency(g,sparse = T)
  m = summary(E)
  m[,1] = m[,1] - 1
  m[,2] = m[,2] - 1
  write.table(m, file_network,append = TRUE, col.names = F, row.names = F)
  
}




generate_two_hops <- function(g,file_hops,negative_factor){
  E = get.adjacency(g,sparse = T)
  E_2=E%*%E
  E_2[E_2 > 0,] = 1
  E_2=E_2-E
  E_2[E_2 < 0, ] = 0
  

  pos_indexes=which(E>0,arr.ind = TRUE)

  neg_indexes=which(E_2>0,arr.ind = TRUE);

  n_positives = nrow(pos_indexes)
  
  n_missing_two_hops=nrow(neg_indexes)
  n_negatives=floor(n_missing_two_hops*negative_factor);

  m = n_positives + n_negatives;

  POS = cbind(pos_indexes[],rep(1,n_positives))
  write.table(POS,file_hops,col.names = FALSE,row.names = FALSE,sep = "\t")
  
  
  permutation_missing_two_hops = sample(n_missing_two_hops,n_missing_two_hops,replace=FALSE);
  

  NEG = cbind(neg_indexes[permutation_missing_two_hops[1:n_negatives],],rep(0,n_negatives))
  
  write.table(NEG,file_hops,col.names = FALSE,row.names = FALSE,sep = "\t", append = TRUE)
}




