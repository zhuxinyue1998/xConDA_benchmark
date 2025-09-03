
safe_intersect_or_null <- function(target, candidates) {
  if (length(candidates) == 0) return(NULL)
  res <- intersect(target, candidates)
  if (length(res) == 0) return(NULL)
  return(res)
}
remove_low_feature_samples <- function(mat, verbose = TRUE) {
  # 每一行中非零元素的个数
  non_zero_count <- rowSums(mat != 0)
  
  # 保留非零特征数 ≥ 2 的样本
  keep_idx <- non_zero_count >= 2
  
  if (verbose) {
    removed <- sum(!keep_idx)
    message(removed, " samples removed with fewer than 2 non-zero features.")
  }
  
  mat[keep_idx, , drop = FALSE]
}
filter_features_in_multiple_batches <- function(data, min_batch_count = 2, verbose = TRUE) {
  meta <- data$metadata
  abs_mat <- data$simulated_data_abs
  rela_mat <- data$simulated_data_rela
  
  # 确保 SampleID 对齐
  stopifnot(rownames(abs_mat) == meta$SampleID)
  
  # 获取所有唯一的 Batch
  batches <- unique(meta$Batch)
  
  # 构建一个向量，记录每个特征在多少个 batch 中有非零值
  feature_batch_count <- setNames(rep(0, ncol(abs_mat)), colnames(abs_mat))
  
  for (batch in batches) {
    samples <- meta$SampleID[meta$Batch == batch]
    if (length(samples) == 0) next
    
    batch_data <- abs_mat[samples, , drop = FALSE]
    feature_present <- colSums(batch_data != 0) > 0  # 特征在该 batch 中是否出现
    feature_batch_count[feature_present] <- feature_batch_count[feature_present] + 1
  }
  
  # 筛选在至少 min_batch_count 个 batch 中出现的特征
  keep_features <- names(feature_batch_count[feature_batch_count >= min_batch_count])
  
  if (verbose) {
    removed <- length(feature_batch_count) - length(keep_features)
    message(removed, " features removed; ", length(keep_features), " retained.")
  }
  
  # 保留筛选后的特征列
  abs_mat_filtered <- abs_mat[, keep_features, drop = FALSE]
  rela_mat_filtered <- rela_mat[, keep_features, drop = FALSE]
  
  # 返回新的 list
  return(list(
    metadata = meta,
    simulated_data_abs = abs_mat_filtered,
    simulated_data_rela = rela_mat_filtered
  ))
}