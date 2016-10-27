match_degree <- function(P, Q, k, h) {
  p_sq <- sum(P[,k]^2)
  q_sq <- sum(Q[,h]^2)
  prod <- sum(P[,k] * Q[,h])^2
  
  if (q_sq == 0)
    return (p_sq)
  else
    return (p_sq - prod / q_sq)
}

chi <- function(P, Q) {
  ncols = ncol(P)
  
  sum <- 0
  j <- 0
  
  for (k in 1:ncols) {
    j <- j + 1
    h_values <- rep(0, times = ncols)
    
    for (h in 1:ncols) {
      h_values[h] <- match_degree(P, Q, k, h)
    }

    print(paste(j,which.min(h_values)))
    sum <- sum + min(h_values)
  }
  
  return(sum)
}



# Main

A_file <- "/Users/ettore/git/SurvivalFactorization/resources/datasets/synth-cwn/s4/s4.communities.nn1000.nc2048.k6.P"
A_model_file <- "/Users/ettore/git/SurvivalFactorization/resources/datasets/synth-cwn/models/full/s4/synth-cwn_s4_6f.model.phi"

A = as.matrix(read.table(A_file, sep = "\t", header = FALSE))
A_tilde = as.matrix(read.table(A_model_file, sep = "\t", header = FALSE))

A = (A - min(A)) / (max(A) -min(A))

sum((A - A_tilde)^2)
chi(A, A_tilde)

A1 = A
A2 = A_tilde


output <- "/Users/ettore/git/SurvivalFactorization/resources/datasets/synth-cwn/s4/s4_P_comparison.pdf"
A1 = cbind(A[,1],A[,2],A[,3]+A[,5],A[,4]+A[,6])
A2 = cbind(A_tilde[,6],A_tilde[,4],A_tilde[,5],A_tilde[,2])

pdf(file=output)
par(mfrow= c(1,2))
cellcol<-color.scale(A1^.4,c(1,0),c(1,0),c(1,0))
color2D.matplot(A1,cellcolors=cellcol,main=expression("True" ~ Phi), border=NA, xlab = NA, ylab = NA)
cellcol<-color.scale(A_star^.4,c(1,0),c(1,0),c(1,0))
color2D.matplot(A2,cellcolors=cellcol,main=expression("Inferred" ~ Phi), border=NA, xlab = NA, ylab = NA)
dev.off()
