root = '/Users/ettore/git/SurvivalFactorization/resources/datasets/twitter2/'
m_data_prefix = ''
m_prefix = 'models/twitter2_'
ma_suffix = 'f.model.A'
ms_suffix = 'f.model.S'
testset = 'test.txt'
output = 'twitter2_result_'

min_train_size = 4
max_cascade_size = 200

topics = c(2, 4, 8, 16, 32, 64)
n = length(topics)

for (p in 5:8) {
  split_perc = p / 10
  i <- 0

  for (k in topics) {
    i <- i + 1
    fileA = paste(root, m_prefix, k, ma_suffix, sep = "")
    fileS = paste(root, m_prefix, k, ms_suffix, sep = "")
    cascades = paste(root, m_data_prefix, testset, sep = "")
  
    A = as.matrix(read.table(fileA, sep = "\t", header = FALSE))
    S = as.matrix(read.table(fileS, sep = "\t", header = FALSE))
    data = read.table(cascades, sep = '\t', header = TRUE)
    
    num_users = max(data[, 1])
    num_cascades = max(data[, 2])
    
    C = sparseMatrix(i=data[,1],j=data[,2],x=data[,3],dims=c(num_users,num_cascades))
    
    set.seed(12345)
    indices = 1:num_cascades
  
    result = gen_predictions(C,indices, A, S, split_perc,min_train_size, max_cascade_size)
    write.table(result, paste(root, output, (p*10), '_', k, 'f.pred', sep = ""), row.names = FALSE, col.names = TRUE)
  }
}
