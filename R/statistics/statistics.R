#####prevalence筛选函数
filter_by_prevalence <- function(matrix, prevalence) {
  # 计算每列的非零比例
  nonzero_proportion <- colSums(matrix != 0) / nrow(matrix)
  
  # 筛选出非零比例大于 prevalence 的列
  filtered_matrix <- matrix[, which(nonzero_proportion > prevalence), drop = FALSE]
  
  # 返回筛选后的矩阵
  return(filtered_matrix)
}