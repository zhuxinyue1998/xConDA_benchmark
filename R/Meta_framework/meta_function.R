
####由于多队列的不同设置对部分差异分析方法的函数在这里进行了一些调整
#对协变量部分进行更改，因为目前筛选协变量已经包含在metadata中
fastANCOM_analysis<-function(simulated_data_all,covariates=NULL,categorical_variable_name = NULL,normalization){
  ###covariates以向量形式输入
  library(fastANCOM)
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  # batch <- simulated_data_all$metadata$Batch
  # batch_data <- model.matrix(~ batch -1)
  # simulated_data_all$metadata <- cbind(simulated_data_all$metadata,batch_data)
  # covariates <- c(setdiff(covariates,'Batch'),colnames(batch_data))
  if(setequal(categorical_variable_name,NULL)){}else{
    for(i in 1:length(categorical_variable_name)){
      simulated_data_all$metadata[[categorical_variable_name[i]]] <- as.numeric(as.factor(simulated_data_all$metadata[[categorical_variable_name[i]]]))
      if(length(unique(simulated_data_all$metadata[[categorical_variable_name[i]]])) > 2){
        another_name <- paste0(categorical_variable_name[i], '_temp')
        assign(another_name, simulated_data_all$metadata[[categorical_variable_name[i]]])
        
        # 使用 as.formula 来动态构造公式
        formula <- as.formula(paste("~", another_name, "-1"))
        temp_data <- model.matrix(formula, data = simulated_data_all$metadata)
        
        simulated_data_all$metadata <- cbind(simulated_data_all$metadata, temp_data)
        covariates <- c(setdiff(covariates, categorical_variable_name[i]), colnames(temp_data))
      }
      
    }
  }
  
  
  if(is.null(covariates)){
    fa_out<-fastANCOM(Y=as.matrix(simulated_data),x=simulated_data_all$metadata$Group,zero_cut = 1)  
  }else{
    fa_out<-fastANCOM(Y=as.matrix(simulated_data),x=simulated_data_all$metadata$Group,Z=as.matrix(simulated_data_all$metadata[,covariates]),zero_cut = 1)
  }
  fa_res<-fa_out$results$final
  fa_res$Feature <- rownames(fa_res)
  colnames(fa_res)[3] <- 'pval'
  colnames(fa_res)[4] <- 'p.adj.val'
  return(fa_res)
}
ANCOMBC_analysis<-function(simulated_data_all,covariates=NULL,normalization){
  
  library(ANCOMBC)
  library(tidyverse)
  library(TreeSummarizedExperiment)
  rownames <- rownames(simulated_data_all$simulated_data_abs)
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  rownames(simulated_data) <- rownames
  se <- TreeSummarizedExperiment(
    assays = list(counts = t(simulated_data)), 
    colData = simulated_data_all$metadata
  )
  
  # 捕获错误并处理
  result <- tryCatch({
    if (is.null(covariates)) {
      covariates <- append(covariates, 'Group')
      ancombc_out <- ancombc2(
        data = se, assay_name = "counts",
        fix_formula = covariates,
        p_adj_method = "BH", prv_cut = 0, lib_cut = 100,
        group = "Group", struc_zero = TRUE, alpha = 0.05, global = TRUE, verbose = TRUE
      )
    } else {
      covariates <- append(covariates, 'Group')
      formula_str <- paste(covariates, collapse = ' + ')
      ancombc_out <- ancombc2(
        data = se, assay_name = "counts",
        fix_formula = formula_str,
        p_adj_method = "BH", prv_cut = 0, lib_cut = 100,
        group = "Group", struc_zero = TRUE, alpha = 0.05, global = TRUE, verbose = TRUE
      )
    }
    ancombc_out
  }, error = function(e) {
    # 提取错误信息中的特征名
    error_message <- conditionMessage(e)
    print(error_message)
    features_to_remove <- str_extract_all(error_message, "Feature\\d+")[[1]]
    
    if (length(features_to_remove) > 0) {
      cat("Removing features with zero variance:", paste(features_to_remove, collapse = ", "), "\n")
      # 从原始数据中移除这些特征
      simulated_data_all$simulated_data_abs <- simulated_data_all$simulated_data_abs[, !colnames(simulated_data_all$simulated_data_abs) %in% features_to_remove]
      # 递归调用自身重新运行
      #print(simulated_data_all$simulated_data_abs)
      simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
      se <- TreeSummarizedExperiment(
        assays = list(counts = t(simulated_data)), 
        colData = simulated_data_all$metadata
      )
      if (is.null(covariates)) {
        covariates <- append(covariates, 'Group')
        ancombc_out <- ancombc2(
          data = se, assay_name = "counts",
          fix_formula = covariates,
          p_adj_method = "BH", prv_cut = 0, lib_cut = 100,
          group = "Group", struc_zero = TRUE, alpha = 0.05, global = TRUE, verbose = TRUE
        )
      } else {
        covariates <- append(covariates, 'Group')
        formula_str <- paste(covariates, collapse = ' + ')
        ancombc_out <- ancombc2(
          data = se, assay_name = "counts",
          fix_formula = formula_str,
          p_adj_method = "BH", prv_cut = 0, lib_cut = 100,
          group = "Group", struc_zero = TRUE, alpha = 0.05, global = TRUE, verbose = TRUE
        )
      }
      ancombc_out
    } else {
      stop("Error occurred but no features identified to remove.")
    }
  })
  
  # 处理结果
  ancombc_res <- result$res
  ancombc_res_process <- ancombc_res
  colnames(ancombc_res_process)[1] <- 'Feature'
  #colnames(ancombc_res_process)[((length(covariates) + 1) * 4 + 2 + length(covariates))] <- 'p.adj.val'
  #colnames(ancombc_res_process)[((length(covariates) + 1) * 3 + 2 + length(covariates))] <- 'pval'
  ancombc_res_process <- ancombc_res_process %>%
    rename_with(~ gsub('^q_Group.*', 'p.adj.val', .), everything()) %>%
    rename_with(~ gsub('^p_Group.*', 'pval', .), everything())
  
  rownames(ancombc_res_process) <- ancombc_res_process$Feature
  na_feature <- setdiff(colnames(simulated_data_all$simulated_data_rela), ancombc_res_process$Feature)
  if (length(na_feature) > 0) {
    filter_feature <- as.data.frame(matrix(NA, nrow = length(na_feature), ncol = ncol(ancombc_res_process)))
    rownames(filter_feature) <- na_feature
    colnames(filter_feature) <- colnames(ancombc_res_process)
    ancombc_res_process <- rbind(ancombc_res_process, filter_feature)
  }
  ancombc_res_process$Feature <- rownames(ancombc_res_process)
  return(ancombc_res_process)
}
edgeR_analysis<-function(simulated_data_all,covariates = NULL,normalization){
  library(edgeR)
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  col_names <- covariates
  covariates <- c("Group",covariates)
  if(length(covariates)>1){
    covariates <- paste(covariates,collapse = ' + ')
  }
  formula_str <- paste0('~ ',covariates)
  Group  <- factor(simulated_data_all$metadata$Group)
  y <- DGEList(counts=t(simulated_data),group=Group)
  
  design <- model.matrix(as.formula(formula_str),data = simulated_data_all$metadata)
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef = 2)
  edgeR_res <- topTags(qlf,n = ncol(simulated_data))
  edgeR_res <- edgeR_res$table
  colnames(edgeR_res)[5] <- 'p.adj.val' 
  colnames(edgeR_res)[4] <- 'pval' 
  edgeR_res$Feature <- rownames(edgeR_res)
  return(edgeR_res)
}
limma_analysis<-function(simulated_data_all,covariates=NULL,normalization){
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  col_names <- covariates
  covariates <- c("Group",covariates)
  if(length(covariates)>1){
    covariates <- paste(covariates,collapse = ' + ')
  }
  formula_str <- paste0('~0+',covariates)
  #print(formula_str)
  library(limma)
  library(edgeR)
  simulated_data_all$metadata$Group <- as.factor(simulated_data_all$metadata$Group)
  mm<-model.matrix(as.formula(formula_str),data = simulated_data_all$metadata)
  colnames(mm)[1:2]<- c('group0','group1')
  fit<-lmFit(t(simulated_data),mm)
  contr <- makeContrasts(group1 - group0, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  limma_res <- top.table
  limma_res$Feature <- rownames(limma_res)
  colnames(limma_res)[5] <- 'p.adj.val'
  colnames(limma_res)[4] <- 'pval'
  return(limma_res)
}
LM_random_effect<-function(simulated_data_all,covariates,categorical_variable_name,normalization){
  library(lmerTest)
  library(MuMIn)
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  lmem_res <- data.frame(Feature = rep(NA,ncol(simulated_data_all$simulated_data_rela)),Pvalue=rep(NA,ncol(simulated_data_all$simulated_data_rela)),Effect_size = rep(NA,ncol(simulated_data_all$simulated_data_rela)),stringsAsFactors = F)
  for(k in 1:ncol(simulated_data_all$simulated_data_rela)) {
    # 创建临时数据
    lmem_res$Feature <- as.character(lmem_res$Feature)
    temp_data <- data.frame(cbind(simulated_data[,k], simulated_data_all$metadata))
    colnames(temp_data)[1] <- 'Feature'
    # 使用 paste 构建随机效应的字符串
    if(length(setdiff(covariates,categorical_variable_name))!=0 & !is.null(categorical_variable_name)){
      random_effects_str <- paste0("(1|",categorical_variable_name, ")")
      fixed_effects_str <- paste(setdiff(covariates,categorical_variable_name),collapse = " + ")
      # 构建完整的公式字符串
      formula_str <- paste("Feature ~ Group", paste(random_effects_str,collapse = " + "), fixed_effects_str,sep = " + ")
    }else if(!is.null(categorical_variable_name)){
      random_effects_str <- paste0("(1|",categorical_variable_name, ")")
      # 构建完整的公式字符串
      formula_str <- paste("Feature ~ Group", paste(random_effects_str,collapse = " + "),sep = " + ")
    }else{
      #random_effects_str <- paste0("(1|",categorical_variable_name, ")")
      fixed_effects_str <- paste(setdiff(covariates,categorical_variable_name),collapse = " + ")
      # 构建完整的公式字符串
      formula_str <- paste("Feature ~ Group", fixed_effects_str,sep = " + ")
    }
    #print(formula_str)
    # 使用 tryCatch 捕获错误
    tryCatch({
      # 运行 lmer 模型
      lm_re <- lmer(as.formula(formula_str), data = temp_data)
      #print(r.squaredGLMM(lm_re))
      # 提取 p-value 和估计效应
      # 获取以 "Group" 开头的系数
      group_coeffs <- grep("^Group", rownames(summary(lm_re)$coefficients), value = TRUE)
      
      if (length(group_coeffs) > 0) {
        pval <- summary(lm_re)$coefficients[group_coeffs, 'Pr(>|t|)']
        estimate_effect <- summary(lm_re)$coefficients[group_coeffs, 'Estimate']
        lmem_res[k, ] <- c(colnames(simulated_data)[k], pval, estimate_effect)
      }
      # pval <- summary(lm_re)$coefficients['Group', 'Pr(>|t|)']
      # #print(pval)
      # estimate_effect <- summary(lm_re)$coefficients['Group', 'Estimate']
      # #print(estimate_effect)
      # # 存储结果
      # lmem_res[k,] <- c(colnames(simulated_data)[k], pval, estimate_effect)
    }, error = function(e) {
      # 如果出现错误，则将结果设为 NA
      lmem_res[k,1] <- colnames(simulated_data)[k]
    }, finally = {
      # 确保赋值总是执行
      if (is.na(lmem_res[k,1])) {
        lmem_res[k,1] <- colnames(simulated_data)[k]
      }
    })
  }
  p.adj.val  = p.adjust(lmem_res$Pvalue,method = 'BH')
  lmem_res$p.adj.val = p.adj.val
  colnames(lmem_res)[2] <- 'pval'
  return(lmem_res)
}
LM_fixed_effect<-function(simulated_data_all,covariates,normalization){
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  covariates_str <- paste(covariates, collapse = " + ")
  # 构建完整的公式字符串
  formula_str <- paste("Feature ~ Group +", covariates_str)
  lmf_res <- data.frame(Feature = rep(NA,ncol(simulated_data_all$simulated_data_rela)),Pvalue=rep(NA,ncol(simulated_data_all$simulated_data_rela)),Effect_size = rep(NA,ncol(simulated_data_all$simulated_data_rela)))
  for(k in 1:ncol(simulated_data_all$simulated_data_rela)){
    temp_data <- data.frame(cbind(simulated_data[,k],simulated_data_all$metadata))
    colnames(temp_data)[1]<-'Feature'
    lmf <- lm(as.formula(formula_str),data = temp_data)
    # 获取以 "Group" 开头的系数
    group_coeffs <- grep("^Group", rownames(summary(lmf)$coefficients), value = TRUE)
    
    if (length(group_coeffs) > 0) {
      pval <- summary(lmf)$coefficients[group_coeffs, 'Pr(>|t|)']
      estimate_effect <- summary(lmf)$coefficients[group_coeffs, 'Estimate']
      lmf_res[k, ] <- c(colnames(simulated_data)[k], pval, estimate_effect)
    }
  }
  p.adj.val  = p.adjust(lmf_res$Pvalue,method = 'BH')
  lmf_res$p.adj.val = p.adj.val
  colnames(lmf_res)[2] <- 'pval'
  return(lmf_res)
}
Maaslin2_LM_analysis<-function(simulated_data_all,covariates=NULL,categorical_variable_name = NULL,normalization,transform){
  covariates <- append(covariates,"Group")
  reference = NULL
  if(setequal(categorical_variable_name,NULL)){}else{
    for(i in 1:length(categorical_variable_name)){
      simulated_data_all$metadata[[categorical_variable_name[i]]] <- as.numeric(as.factor(simulated_data_all$metadata[[categorical_variable_name[i]]]))
      if(length(unique(simulated_data_all$metadata[[categorical_variable_name[i]]])) > 2){
        reference <- c(reference,paste(categorical_variable_name[i],unique(simulated_data_all$metadata[[categorical_variable_name[i]]])[1],sep=','))
      }
      
    }
  }
  library(Maaslin2)
  meta <- simulated_data_all$metadata
  data <- simulated_data_all$simulated_data_abs
  filename <- "data_maaslin2_lm"
  rownames(meta)<-meta$Sample
  fit_data = Maaslin2(
    input_data = data, 
    input_metadata = meta, 
    min_prevalence = 0,
    normalization = normalization,analysis_method = "LM",
    output = filename,transform = transform,plot_heatmap = FALSE,plot_scatter = FALSE,
    fixed_effects = covariates,
    random_effects = NULL)
  maa_lm_res<-fit_data$results
  colnames(maa_lm_res)[1] <- 'Feature'
  colnames(maa_lm_res)[8] <- 'p.adj.val'
  return(maa_lm_res)
}
###meta函数
library(metap)
meta_combine_p <- function(result_list, sample_size,method='Stouffer') {
  # 获取所有唯一的 Feature
  all_features <- unique(unlist(lapply(result_list, function(df) df$Feature)))
  
  # 按照 all_features 更新每个 dataframe
  result_list <- lapply(result_list, function(df) {
    # 创建一个包含所有 Feature 的新 dataframe
    updated_df <- data.frame(Feature = all_features)
    # 合并原始数据，缺失部分填充 NA
    updated_df <- merge(updated_df, df, by = "Feature", all.x = TRUE)
    # 返回更新后的 dataframe
    return(updated_df)
  })
  #print(result_list)
  # 计算 cohort 数量
  cohort_num <- length(result_list)
  
  # 初始化合并的结果 dataframe
  combined_p_data <- data.frame(
    Feature = all_features,
    pval = rep(NA, length(all_features)),
    p.adj.val = rep(NA, length(all_features))
  )
  
  # 计算每个样本的权重
  sample_size_weights <- sqrt(sample_size)
  
  # 计算合并的 p 值
  for (i in 1:length(all_features)) {
    # 提取当前 Feature 的所有 p 值
    p <- unlist(lapply(result_list, function(df) {
      df$pval[which(df$Feature == all_features[i])]
    }))
    
    # 删除 NA 值和对应的权重
    valid_indices <- !is.na(p)
    p <- p[valid_indices]
    weights <- sample_size_weights[valid_indices]
    # 检查有效 P 值数量
    valid_p_count <- sum(valid_indices)
    
    # 如果有效 P 值少于 2，跳过计算
    if (valid_p_count < 2) {
      combined_p_data$pval[i] <- NA
    } else {
      # 处理 P 值不全为 NA 的情况
      if(method == 'Fisher'){
        combined_p_val <- sumlog(p = as.numeric(p), log.p = FALSE,log.input = FALSE)$p
        combined_p_data$pval[i] <- combined_p_val
      }else if(method == 'Stouffer'){
        combined_p_val <- sumz(p = p, weights = weights)$p
        combined_p_data$pval[i] <- combined_p_val
      }
    }
  }
  
  # 计算调整后的 p 值
  combined_p_data$p.adj.val <- p.adjust(combined_p_data$pval, method = 'BH')
  
  return(combined_p_data)
}

meta_combine<-function(effect_merge_data,se_merge_data,cohort_num,
                       rma_conv = 1e-10,method_type = method_type,
                       rma_maxit = 1000){
  ind_feature <- !is.na(effect_merge_data) & !is.na(se_merge_data)
  effect_merge_data <- effect_merge_data[rowSums(ind_feature)>1,]
  se_merge_data <- se_merge_data[rowSums(ind_feature)>1,]
  ind_feature<-ind_feature[rowSums(ind_feature)>1,]
  #ate_merge_data<-ate_merge_data[rowSums(ind_feature)==cohort_num,]
  #ate_se_merge_data<-ate_se_merge_data[rowSums(ind_feature)==cohort_num,]
  #ind_feature<-ind_feature[rowSums(ind_feature)==cohort_num,]
  rownames(ind_feature)<-rownames(effect_merge_data[rowSums(ind_feature)>1,])
  result <- data.frame(matrix(NA,
                              nrow = nrow(effect_merge_data),
                              ncol = 9 + cohort_num))
  colnames(result) <- c("coef",
                        "stderr",
                        "pval",
                        "k",
                        "tau2",
                        "stderr.tau2",
                        "pval.tau2",
                        "I2",
                        "H2",
                        paste0("weight_",
                               colnames(effect_merge_data)))
  rownames(result)<- rownames(effect_merge_data)
  for(feature in rownames(effect_merge_data)){
    #print(feature)
    rma_fit <-metafor::rma.uni(yi = unlist(c(effect_merge_data[feature,ind_feature[feature,]])),
                               sei = unlist(c(se_merge_data[feature,ind_feature[feature,]])),
                               method=method_type,verbose = F,
                               control = list(threshold = rma_conv,
                                              maxiter = rma_maxit))
    
    wts <- metafor::weights.rma.uni(rma_fit)
    names(wts) <- colnames(effect_merge_data)[ind_feature[feature,]]
    result[feature, c("coef",
                      "stderr",
                      "pval",
                      "k",
                      "tau2",
                      "stderr.tau2",
                      "pval.tau2",
                      "I2",
                      "H2",
                      paste0("weight_",
                             names(wts))
    )] <- c(unlist(rma_fit[c("beta",
                             "se",
                             "pval",
                             "k",
                             "tau2",
                             "se.tau2",
                             "QEp",
                             "I2",
                             "H2")]),
            wts)
  }
  result$Bonferroni_pval<-p.adjust(result$pval, method = "bonf")
  result$BH_pval<-p.adjust(result$pval, method = "BH")
  result<-cbind(rownames(result),result)
  colnames(result)[1]<-"Feature"
  return(result)
}


run_meta <- function(method_name,normalization,meta_method,data,transform=NULL,data_index){
  #####提取covariates,categorical_variable_name
  covariates <- setdiff(colnames(data$metadata),c('Group','SampleID','Batch'))
  categorical_variable_name <- safe_intersect_or_null(c("gender"), covariates)
  output_file_name <- paste(method_name,normalization,sep ='_')
  method_function_map <- list(
    edgeR = edgeR_analysis,
    lmem = LM_random_effect,
    limma = limma_analysis,
    lfem = LM_fixed_effect,
    ANCOMBC = ANCOMBC_analysis,
    Maaslin2 = Maaslin2_LM_analysis,
    fastANCOM = fastANCOM_analysis
  )
  selected_method_func <- method_function_map[[method_name]]
  data <- filter_features_in_multiple_batches(data, min_batch_count = 2, verbose = TRUE)
  # 提取所有唯一的 Batch 值
  batch_levels <- unique(data$metadata$Batch)
  
  # 创建一个列表保存每个 Batch 对应的样本 ID
  batch_sample_list <- lapply(batch_levels, function(batch_val) {
    data$metadata$SampleID[data$metadata$Batch == batch_val]
  })
  names(batch_sample_list) <- batch_levels
  
  
  batch_data_list <- list()
  
  for (batch_name in names(batch_sample_list)) {
    sample_ids <- batch_sample_list[[batch_name]]
    
    batch_data <- list(
      simulated_data_rela = data$simulated_data_rela[sample_ids, , drop = FALSE],
      simulated_data_abs = data$simulated_data_abs[sample_ids, , drop = FALSE],
      metadata = data$metadata[data$metadata$SampleID %in% sample_ids, , drop = FALSE]
    )
    
    # 保存到列表中（可选命名）
    batch_data_list[[batch_name]] <- batch_data
  }
  data$simulated_data_abs <- remove_low_feature_samples(data$simulated_data_abs)
  data$simulated_data_abs <- filter_by_prevalence(data$simulated_data_abs,0.1)
  data$simulated_data_rela <- remove_low_feature_samples(data$simulated_data_rela)
  data$simulated_data_rela <- filter_by_prevalence(data$simulated_data_rela,0.1)
  data$metadata <- data$metadata[data$metadata$SampleID %in% rownames(data$simulated_data_abs),]
  for (batch_name in names(batch_data_list)) {
    batch_data_list[[batch_name]]$simulated_data_abs <-
      remove_low_feature_samples(batch_data_list[[batch_name]]$simulated_data_abs)
    batch_data_list[[batch_name]]$simulated_data_abs <-
      filter_by_prevalence(batch_data_list[[batch_name]]$simulated_data_abs, 0.1)
    batch_data_list[[batch_name]]$simulated_data_rela <-
      remove_low_feature_samples(batch_data_list[[batch_name]]$simulated_data_rela)
    batch_data_list[[batch_name]]$simulated_data_rela <-
      filter_by_prevalence(batch_data_list[[batch_name]]$simulated_data_rela, 0.1)
    batch_data_list[[batch_name]]$metadata <-
      batch_data_list[[batch_name]]$metadata[batch_data_list[[batch_name]]$metadata$SampleID %in% rownames(batch_data_list[[batch_name]]$simulated_data_abs), ]
  }
  
  if(!method_name %in% c('lmem','Maaslin2','fastANCOM')){
    if(meta_method == 'rma.uni'){
      # 1. 循环调用 selected_method_func，并将结果保存在列表中
      result_list <- lapply(batch_data_list, function(data) {
        selected_method_func(data, covariates, normalization)
      })
      
      # 2. 提取所有 Feature 的并集
      all_feature <- Reduce(union, lapply(result_list, function(res) res$Feature))
      
      # 3. 初始化合并数据框
      effect_merge_data <- data.frame(Feature = all_feature)
      se_merge_data <- data.frame(Feature = all_feature)
      
      # 4. 逐一提取 lfc 和 se，按照 Feature 对齐
      batch_names <- names(batch_data_list)
      for (i in seq_along(result_list)) {
        res <- result_list[[i]]
        batch <- batch_names[i]
        effect_merge_data[[paste0("effect", batch)]] <- res$lfc_Group1[match(all_feature, res$Feature)]
        se_merge_data[[paste0("se", batch)]] <- res$se_Group1[match(all_feature, res$Feature)]
      }
      
      # 5. 设置 rownames，移除 Feature 列
      rownames(effect_merge_data) <- all_feature
      rownames(se_merge_data) <- all_feature
      effect_merge_data <- effect_merge_data[, -1]
      se_merge_data <- se_merge_data[, -1]
      method_res_REML<-meta_combine(effect_merge_data,se_merge_data,length(batch_data_list),rma_conv = 1e-10,method_type = 'REML',rma_maxit = 1000)
      method_res_EB<-meta_combine(effect_merge_data,se_merge_data,length(batch_data_list),rma_conv = 1e-10,method_type = 'EB',rma_maxit = 1000)
      method_res_PM<-meta_combine(effect_merge_data,se_merge_data,length(batch_data_list),rma_conv = 1e-10,method_type = 'PM',rma_maxit = 1000)
      write.table(method_res_REML,file = paste0(output_file_name,'_res_',data_index,'_',meta_method,'_REML','.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
      write.table(method_res_EB,file = paste0(output_file_name,'_res_',data_index,'_',meta_method,'_EB','.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
      write.table(method_res_PM,file = paste0(output_file_name,'_res_',data_index,'_',meta_method,'_PM','.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    }else if(meta_method == 'Stouffer'){
      result_list <- lapply(batch_data_list, function(data) {
        selected_method_func(data, covariates, normalization)
      })
      sample_sizes <- sapply(batch_sample_list, length)
      
      method_res <- meta_combine_p(
        result_list = result_list,
        sample_size = sample_sizes,
        method = "Stouffer"
      )
      stopifnot(identical(names(batch_data_list), names(batch_sample_list)))
      #method_res<-meta_combine_p(result_list = list(data1_res,data2_res,data3_res), sample_size = c(length(data1_sample),length(data2_sample),length(data3_sample)),method = 'Stouffer')
      write.table(method_res,file = paste0(output_file_name,'_res_',data_index,'_',meta_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    }else if(meta_method == 'Fisher'){
      result_list <- lapply(batch_data_list, function(data) {
        selected_method_func(data, covariates, normalization)
      })
      sample_sizes <- sapply(batch_sample_list, length)
      method_res <- meta_combine_p(
        result_list = result_list,
        sample_size = sample_sizes,
        method = "Fisher"
      )
      stopifnot(identical(names(batch_data_list), names(batch_sample_list)))
      #method_res<-meta_combine_p(result_list = list(data1_res,data2_res,data3_res), sample_size = c(length(data1_sample),length(data2_sample),length(data3_sample)),method = 'Fisher')
      write.table(method_res,file = paste0(output_file_name,'_res_',data_index,'_',meta_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    }else if(meta_method == 'not_used_meta'){
      covariates <- append(covariates,'Batch')
      method_res<-selected_method_func(data,covariates,normalization)
      write.table(method_res,file = paste0(output_file_name,'_res_',data_index,'_no_meta.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    }
  }else{
    if(method_name == 'Maaslin2'){
      if(meta_method == 'rma.uni'){
        # 1. 执行方法函数
        result_list <- lapply(batch_data_list, function(data) {
          selected_method_func(data, covariates, categorical_variable_name, normalization, transform)
        })
        
        # 2. 筛选 metadata == "Group"
        result_list <- lapply(result_list, function(res) {
          res[which(res$metadata == "Group"), ]
        })
        
        # 3. 获取所有 feature 并集
        all_feature <- Reduce(union, lapply(result_list, function(res) res$Feature))
        
        # 4. 初始化数据框
        effect_merge_data <- data.frame(Feature = all_feature)
        se_merge_data <- data.frame(Feature = all_feature)
        
        # 5. 添加命名列
        batch_names <- names(batch_data_list)
        
        for (i in seq_along(result_list)) {
          res <- result_list[[i]]
          batch <- batch_names[i]
          effect_merge_data[[paste0("effect_", batch)]] <- res$coef[match(all_feature, res$Feature)]
          se_merge_data[[paste0("se_", batch)]] <- res$stderr[match(all_feature, res$Feature)]
        }
        
        # 6. 设置行名 & 删除 Feature 列
        rownames(effect_merge_data) <- all_feature
        rownames(se_merge_data) <- all_feature
        effect_merge_data <- effect_merge_data[, -1]
        se_merge_data <- se_merge_data[, -1]
        method_res_REML<-meta_combine(effect_merge_data,se_merge_data,length(batch_data_list),rma_conv = 1e-10,method_type = 'REML',rma_maxit = 1000)
        method_res_EB<-meta_combine(effect_merge_data,se_merge_data,length(batch_data_list),rma_conv = 1e-10,method_type = 'EB',rma_maxit = 1000)
        method_res_PM<-meta_combine(effect_merge_data,se_merge_data,length(batch_data_list),rma_conv = 1e-10,method_type = 'PM',rma_maxit = 1000)
        write.table(method_res_REML,file = paste0(output_file_name,'_res_',data_index,'_',meta_method,'_REML','.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
        write.table(method_res_EB,file = paste0(output_file_name,'_res_',data_index,'_',meta_method,'_EB','.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
        write.table(method_res_PM,file = paste0(output_file_name,'_res_',data_index,'_',meta_method,'_PM','.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
      }else if(meta_method == 'Stouffer'){
        # 1. 执行方法函数
        result_list <- lapply(batch_data_list, function(data) {
          selected_method_func(data, covariates, categorical_variable_name, normalization, transform)
        })
        
        # 2. 筛选 metadata == "Group"
        result_list <- lapply(result_list, function(res) {
          res[which(res$metadata == "Group"), ]
        })
        sample_sizes <- sapply(batch_sample_list, length)
        
        method_res <- meta_combine_p(
          result_list = result_list,
          sample_size = sample_sizes,
          method = "Stouffer"
        )
        stopifnot(identical(names(batch_data_list), names(batch_sample_list)))
        
        #method_res<-meta_combine_p(result_list = list(data1_res,data2_res,data3_res), sample_size = c(length(data1_sample),length(data2_sample),length(data3_sample)),method = 'Stouffer')
        write.table(method_res,file = paste0(output_file_name,'_res_',data_index,'_',meta_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
      }else if(meta_method == 'Fisher'){
        # 1. 执行方法函数
        result_list <- lapply(batch_data_list, function(data) {
          selected_method_func(data, covariates, categorical_variable_name, normalization, transform)
        })
        
        # 2. 筛选 metadata == "Group"
        result_list <- lapply(result_list, function(res) {
          res[which(res$metadata == "Group"), ]
        })
        sample_sizes <- sapply(batch_sample_list, length)
        
        method_res <- meta_combine_p(
          result_list = result_list,
          sample_size = sample_sizes,
          method = "Fisher"
        )
        stopifnot(identical(names(batch_data_list), names(batch_sample_list)))
        #method_res<-meta_combine_p(result_list = list(data1_res,data2_res,data3_res), sample_size = c(length(data1_sample),length(data2_sample),length(data3_sample)),method = 'Fisher')
        write.table(method_res,file = paste0(output_file_name,'_res_',data_index,'_',meta_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
      }else if(meta_method == 'not_used_meta'){
        covariates <- append(covariates,'Batch')
        categorical_variable_name <- append(categorical_variable_name,'Batch')
        method_res<-selected_method_func(data,covariates,categorical_variable_name,normalization,transform)
        method_res<-method_res[which(method_res$metadata=='Group'),]
        write.table(method_res,file = paste0(output_file_name,'_res_',data_index,'_no_meta.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
      }
    }else{
      if(meta_method == 'Stouffer'){
        # 1. 执行方法函数
        result_list <- lapply(batch_data_list, function(data) {
          selected_method_func(data, covariates, categorical_variable_name, normalization)
        })
        sample_sizes <- sapply(batch_sample_list, length)
        
        method_res <- meta_combine_p(
          result_list = result_list,
          sample_size = sample_sizes,
          method = "Stouffer"
        )
        stopifnot(identical(names(batch_data_list), names(batch_sample_list)))
        #method_res<-meta_combine_p(result_list = list(data1_res,data2_res,data3_res), sample_size = c(length(data1_sample),length(data2_sample),length(data3_sample)),method = 'Stouffer')
        write.table(method_res,file = paste0(output_file_name,'_res_',data_index,'_',meta_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
      }else if(meta_method == 'Fisher'){
        # 1. 执行方法函数
        result_list <- lapply(batch_data_list, function(data) {
          selected_method_func(data, covariates, categorical_variable_name, normalization)
        })
        sample_sizes <- sapply(batch_sample_list, length)
        
        method_res <- meta_combine_p(
          result_list = result_list,
          sample_size = sample_sizes,
          method = "Fisher"
        )
        stopifnot(identical(names(batch_data_list), names(batch_sample_list)))
        #method_res<-meta_combine_p(result_list = list(data1_res,data2_res,data3_res), sample_size = c(length(data1_sample),length(data2_sample),length(data3_sample)),method = 'Fisher')
        write.table(method_res,file = paste0(output_file_name,'_res_',data_index,'_',meta_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
      }else if(meta_method == 'not_used_meta'){
        covariates <- append(covariates,'Batch')
        categorical_variable_name <- append(categorical_variable_name,'Batch')
        method_res<-selected_method_func(data,covariates,categorical_variable_name,normalization)
        write.table(method_res,file = paste0(output_file_name,'_res_',data_index,'_no_meta.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
      }
    }
  }
  
}