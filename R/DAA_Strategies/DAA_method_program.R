#####所有方法运行的函数
##edgeR
##不能使用CLR，LOG标准化
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
  edgeR_res$Feature <- rownames(edgeR_res)
  return(edgeR_res)
}


##DEseq2
##只能输入count
# DESeq2_analysis<-function(simulated_data_all,covariates=NULL,normalization){
#   simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
#   col_names <- covariates
#   covariates <- c("Group",covariates)
#   if(length(covariates)>1){
#     covariates <- paste(covariates,collapse = ' + ')
#   }
#   formula_str <- paste0('~ ',covariates)
#   library(DESeq2)
#   library(dplyr)
#   simulated_data_all$metadata$Group <- as.factor(simulated_data_all$metadata$Group)
#   simulated_data_all$metadata <- simulated_data_all$metadata %>%
#     mutate(across(contains("Ca"), as.factor))
#   dds <- DESeqDataSetFromMatrix(countData=t(simulated_data) + 1, 
#                                 colData=simulated_data_all$metadata, 
#                                 design=as.formula(formula_str), tidy = FALSE)
#   dds <- DESeq(dds)
#   DESeq2_res <- data.frame(results(dds))
#   DESeq2_res$Feature <- rownames(DESeq2_res)
#   colnames(DESeq2_res)[which(colnames(DESeq2_res)=='padj')] <- 'p.adj.val'
#   return(DESeq2_res)
# }
DESeq2_analysis<-function(simulated_data_all,covariates=NULL,normalization){
  
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  col_names <- covariates
  covariates <- c("Group",covariates)
  if(length(covariates)>1){
    covariates <- paste(covariates,collapse = ' + ')
  }
  formula_str <- paste0('~ ',covariates)
  library(DESeq2)
  library(dplyr)
  simulated_data_all$metadata$Group <- as.factor(simulated_data_all$metadata$Group)
  simulated_data_all$metadata <- simulated_data_all$metadata %>%
    mutate(across(contains("Ca"), as.factor))
  dds <- DESeqDataSetFromMatrix(countData=t(simulated_data) + 1, 
                                colData=simulated_data_all$metadata, 
                                design=as.formula(formula_str), tidy = FALSE)
  dds <- DESeq(dds)
  
  DESeq2_res <- data.frame(results(dds,contrast = c('Group',1,0)))
  DESeq2_res$Feature <- rownames(DESeq2_res)
  colnames(DESeq2_res)[which(colnames(DESeq2_res)=='padj')] <- 'p.adj.val'
  #head(results(dds,contrast = c('Group',1,0)))
  return(DESeq2_res)
}
##Aldex2
##只能用conunt数据
ALDEx2_analysis <- function(simulated_data_all,covariates=NULL,normalization){
  library(ALDEx2)
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  col_names <- covariates
  covariates <- c("Group",covariates)
  if(length(covariates)>1){
    covariates <- paste(covariates,collapse = ' + ')
  }
  formula_str <- paste0('~ ',covariates)
  mm <- model.matrix(as.formula(formula_str), simulated_data_all$metadata)
  x.glm <- aldex.clr(t(simulated_data), mm, denom="all", verbose=F)
  glm.test <- aldex.glm(x.glm, mm, fdr.method='BH')
  aldex2_res <- glm.test
  aldex2_res$Feature  <- rownames(aldex2_res)
  colnames(aldex2_res)[which(colnames(aldex2_res) == 'Group:pval.padj')] <- 'p.adj.val'
  return(aldex2_res)
}


##Zicoseq
##可以输入各种标准化
###行是特征，列是样本

ZicoSeq_analysis <- function(simulated_data_all,covariates=NULL,normalization){
  library(GUniFrac)
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  #count数据输入
  if(normalization == 'count'){
    ZicoSeq.obj <- ZicoSeq(meta.dat = simulated_data_all$metadata, feature.dat = t(simulated_data), 
                           grp.name = 'Group', adj.name = covariates, feature.dat.type = "count",
                           # Winsorization to replace outliers
                           is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                           # Posterior sampling 
                           is.post.sample = TRUE, post.sample.no = 25, 
                           # Use the square-root transformation
                           link.func = list(function (x) x^0.5), stats.combine.func = max,
                           # Permutation-based multiple testing correction
                           perm.no = 999,  strata = NULL, 
                           # Reference-based multiple stage normalization
                           ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                           # Family-wise error rate control
                           is.fwer = FALSE, verbose = FALSE, return.feature.dat = FALSE)
  }
  #相对丰度数据输入
  if(normalization == 'TSS'){
    ZicoSeq.obj <- ZicoSeq(meta.dat = simulated_data_all$metadata, feature.dat = t(simulated_data), 
                           grp.name = 'Group', adj.name = covariates, feature.dat.type = "proportion",
                           # Winsorization to replace outliers
                           is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                           # Posterior sampling 
                           is.post.sample = FALSE, post.sample.no = 25, 
                           # Use the square-root transformation
                           link.func = list(function (x) x^0.5), stats.combine.func = max,
                           # Permutation-based multiple testing correction
                           perm.no = 999,  strata = NULL, 
                           # Reference-based multiple stage normalization
                           ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                           # Family-wise error rate control
                           is.fwer = FALSE, verbose = FALSE, return.feature.dat = FALSE)
  }
  #其他类型数据输入
  if((normalization != 'TSS') && (normalization != 'count')){
    ZicoSeq.obj <- ZicoSeq(meta.dat = simulated_data_all$metadata, feature.dat = t(simulated_data), 
                           grp.name = 'Group', adj.name = covariates, feature.dat.type = "other",
                           # Winsorization to replace outliers
                           is.winsor = FALSE,
                           # Posterior sampling 
                           is.post.sample = FALSE, post.sample.no = 25, 
                           # Use the square-root transformation
                           link.func = list(function (x) x^0.5), stats.combine.func = max,
                           # Permutation-based multiple testing correction
                           perm.no = 999,  strata = NULL, 
                           # Reference-based multiple stage normalization
                           ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                           # Family-wise error rate control
                           is.fwer = FALSE, verbose = FALSE, return.feature.dat = FALSE)
  }
  ZicoSeq_res<-data.frame(Feature = names(ZicoSeq.obj$p.adj.fdr),effect_size = ZicoSeq.obj$R2[,1],p.adj.val = ZicoSeq.obj$p.adj.fdr)
  return(ZicoSeq_res)
}

####要么使用默认设置要么使用推荐设置
####注意不同标准化可能设置的参数不同
##Corncob
##从输入来看应该只能输入count不管其他标准化方式能不能运行
Corncob_analysis<-function(simulated_data_all,covariates=NULL,normalization){
  library(corncob)
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  col_names <- covariates
  covariates <- c("Group",covariates)
  if(length(covariates)>1){
    covariates <- paste(covariates,collapse = ' + ')
  }
  phi.formula_str <- paste0('~ ',covariates)
  formula_str <- paste0('cbind(W, M - W) ~ ',covariates)
  pval_vector <- c()
  Feature <- c()
  effect_size_vector <- c()
  for(p in 1:ncol(simulated_data)){
    one_feature_data <- cbind(simulated_data_all$metadata, 
                              W = unlist(t(simulated_data)[p, ]),
                              M = colSums(t(simulated_data_all$simulated_data_abs)))
    
    # 使用 tryCatch 捕获错误
    result <- tryCatch({
      corncob <- bbdml(formula = as.formula(formula_str),
                       phi.formula = as.formula(phi.formula_str),
                       data = one_feature_data)
      pval <- data.frame(summary(corncob)$coefficients)['mu.Group', 4]
      effect_size <- data.frame(summary(corncob)$coefficients)['mu.Group', 1]
      c(pval,effect_size)  # 返回 pval 作为结果
    }, error = function(e) {
      message(sprintf("Error at p = %d: %s", p, e$message))  # 打印错误信息
      c(NA,NA)  # 如果发生错误，返回 NA
    })
    Feature <- append(Feature,colnames(simulated_data)[p])
    pval_vector <- append(pval_vector,result[1])
    effect_size_vector <- append(effect_size_vector,result[2])
  }
  p.adj.val = p.adjust(pval_vector,method = 'BH')
  Corncob_res <- data.frame(Feature = Feature, p.adj.val = p.adj.val,effect_size = effect_size_vector)
  return(Corncob_res)
}
fastANCOM_analysis<-function(simulated_data_all,covariates=NULL,normalization){
  ###covariates以向量形式输入
  library(fastANCOM)
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  if(is.null(covariates)){
    fa_out<-fastANCOM(Y=as.matrix(simulated_data),x=simulated_data_all$metadata$Group,zero_cut = 1)  
  }else{
    fa_out<-fastANCOM(Y=as.matrix(simulated_data),x=simulated_data_all$metadata$Group,Z=as.matrix(simulated_data_all$metadata[,covariates]),zero_cut = 1)
  }
  fa_res<-fa_out$results$final
  fa_res$Feature <- rownames(fa_res)
  colnames(fa_res)[4] <- 'p.adj.val'
  return(fa_res)
}
ANCOM_analysis<-function(simulated_data_all,covariates=NULL,normalization){
  
  library(nlme)
  library(tidyverse)
  library(compositions)
  source("/home/data/ZXY/meta_causal/benchmark_DA_pipeline_deep_seq/programs/ancom.R")
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  ####step 1 preprocess
  data_abs_abd_ancom_preprocess = feature_table_pre_process(feature_table = t(simulated_data),meta_data = simulated_data_all$metadata,sample_var = "SampleID",
                                                            group_var = "Group",out_cut = 0.05,zero_cut = 0.9,lib_cut = 10,neg_lb = TRUE)
  ####step 2 ANCOM
  tryCatch({
  if(is.null(covariates)){
    ancom_res = ANCOM(feature_table = data_abs_abd_ancom_preprocess$feature_table,meta_data = data_abs_abd_ancom_preprocess$meta_data,struc_zero = data_abs_abd_ancom_preprocess$structure_zeros,
                      main_var = 'Group',p_adj_method = "BH",alpha = 0.05)
  }else{
    covariates <- paste(covariates,collapse = ' + ')
    ancom_res = ANCOM(feature_table = data_abs_abd_ancom_preprocess$feature_table,meta_data = data_abs_abd_ancom_preprocess$meta_data,struc_zero = data_abs_abd_ancom_preprocess$structure_zeros,
                      main_var = 'Group',p_adj_method = "BH",alpha = 0.05,adj_formula = covariates)
  }
  
  ancom_res_process = ancom_res$out
  ###这里显著的标准按照作者github中描述的commonly used detected_0.7
  ancom_res_process = ancom_res_process[ancom_res_process$detected_0.7,]
  ####这里的p.adj.val是为了评估时函数运行方便，实际上并不是p.adj.val值
  ancom_res_process = data.frame(Feature = rownames(t(simulated_data_all$simulated_data_abs)),p.adj.val = ifelse(rownames(t(simulated_data_all$simulated_data_abs)) %in% ancom_res_process$taxa_id,0,1))
  return(ancom_res_process)
  },error = function(e){
    # 捕捉特定错误信息
    if (grepl("grouping factor must have exactly 2 levels", conditionMessage(e))) {
      cat("Error: Grouping factor must have exactly 2 levels. Returning NA DataFrame.\n")
      
      # 构建返回的 DataFrame
      return(data.frame(
        Feature = rownames(t(simulated_data_all$simulated_data_abs)),
        p.adj.val = NA
      ))
    }else{
      # 如果不是特定错误，不输出内容但是仍然返回NA
      return(data.frame(
        Feature = rownames(t(simulated_data_all$simulated_data_abs)),
        p.adj.val = NA
      ))
    }
  })
}
ANCOMBC_analysis<-function(simulated_data_all,covariates=NULL,normalization){
  
  library(ANCOMBC)
  library(tidyverse)
  library(TreeSummarizedExperiment)
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
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
  colnames(ancombc_res_process)[((length(covariates) + 1) * 4 + 2 + length(covariates))] <- 'p.adj.val'
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
Maaslin2_LM_analysis<-function(simulated_data_all,covariates=NULL,normalization,transform){
  covariates <- append(covariates,"Group")
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
Maaslin2_CPLM_analysis<-function(simulated_data_all,covariates=NULL,normalization,transform){
  covariates <- append(covariates,"Group")
  library(Maaslin2)
  meta <- simulated_data_all$metadata
  data <- simulated_data_all$simulated_data_abs
  filename <- "data_maaslin2_cplm"
  rownames(meta)<-meta$Sample
  fit_data = Maaslin2(
    input_data = data, 
    input_metadata = meta, 
    min_prevalence = 0,
    normalization = normalization,analysis_method = "CPLM",
    output = filename,transform = transform,plot_heatmap = FALSE,plot_scatter = FALSE,
    fixed_effects = covariates,
    random_effects = NULL)
  maa_cplm_res<-fit_data$results
  colnames(maa_cplm_res)[1] <- 'Feature'
  colnames(maa_cplm_res)[8] <- 'p.adj.val'
  return(maa_cplm_res)
}
Maaslin2_ZINB_analysis<-function(simulated_data_all,covariates=NULL,normalization,transform){
  covariates <- append(covariates,"Group")
  library(Maaslin2)
  meta <- simulated_data_all$metadata
  data <- simulated_data_all$simulated_data_abs
  filename <- "data_maaslin2_zinb"
  rownames(meta)<-meta$Sample
  fit_data = Maaslin2(
    input_data = data, 
    input_metadata = meta, 
    min_prevalence = 0,
    normalization = normalization,analysis_method = "ZINB",
    output = filename,transform = transform,plot_heatmap = FALSE,plot_scatter = FALSE,
    fixed_effects = covariates,
    random_effects = NULL)
  maa_zinb_res<-fit_data$results
  colnames(maa_zinb_res)[1] <- 'Feature'
  colnames(maa_zinb_res)[8] <- 'p.adj.val'
  return(maa_zinb_res)
}
Maaslin2_NEGBIN_analysis<-function(simulated_data_all,covariates=NULL,normalization,transform){
  covariates <- append(covariates,"Group")
  library(Maaslin2)
  meta <- simulated_data_all$metadata
  data <- simulated_data_all$simulated_data_abs
  filename <- "data_maaslin2_negbin"
  rownames(meta)<-meta$Sample
  fit_data = Maaslin2(
    input_data = data, 
    input_metadata = meta, 
    min_prevalence = 0,
    normalization = normalization,analysis_method = "NEGBIN",
    output = filename,transform = transform,plot_heatmap = FALSE,plot_scatter = FALSE,
    fixed_effects = covariates,
    random_effects = NULL)
  maa_negbin_res<-fit_data$results
  colnames(maa_negbin_res)[1] <- 'Feature'
  colnames(maa_negbin_res)[8] <- 'p.adj.val'
  return(maa_negbin_res)
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
  colnames(mm)<- c(c('group0','group1'),col_names)
  fit<-lmFit(t(simulated_data),mm)
  contr <- makeContrasts(group1 - group0, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  limma_res <- top.table
  limma_res$Feature <- rownames(limma_res)
  colnames(limma_res)[5] <- 'p.adj.val'
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
    random_effects_str <- paste0("(1|",categorical_variable_name, ")")
    fixed_effects_str <- paste(setdiff(covariates,categorical_variable_name),collapse = " + ")
    # 构建完整的公式字符串
    formula_str <- paste("Feature ~ Group", random_effects_str, fixed_effects_str,sep = " + ")
    #print(formula_str)
    # 使用 tryCatch 捕获错误
    tryCatch({
      # 运行 lmer 模型
      lm_re <- lmer(as.formula(formula_str), data = temp_data)
      #print(r.squaredGLMM(lm_re))
      # 提取 p-value 和估计效应
      pval <- summary(lm_re)$coefficients['Group', 'Pr(>|t|)']
      #print(pval)
      estimate_effect <- summary(lm_re)$coefficients['Group', 'Estimate']
      #print(estimate_effect)
      # 存储结果
      lmem_res[k,] <- c(colnames(simulated_data)[k], pval, estimate_effect)
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
  return(lmem_res)
}
LM_interaction_analysis<-function(simulated_data_all,covariates,normalization){
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  covariates_str <- paste(covariates, collapse = " + ")
  interactions_str <- paste("Group:", covariates, collapse = " + ")
  # 构建完整的公式字符串
  formula_str <- paste("Feature ~ Group +", covariates_str, "+", interactions_str)
  lminter_res <- data.frame(Feature = rep(NA,ncol(simulated_data_all$simulated_data_rela)),Pvalue=rep(NA,ncol(simulated_data_all$simulated_data_rela)),Effect_size = rep(NA,ncol(simulated_data_all$simulated_data_rela)))
  for(k in 1:ncol(simulated_data_all$simulated_data_rela)){
    temp_data <- data.frame(cbind(simulated_data[,k],simulated_data_all$metadata))
    colnames(temp_data)[1]<-'Feature'
    lm_inter <- lm(as.formula(formula_str),data = temp_data)
    pval <- summary(lm_inter)$coefficients['Group','Pr(>|t|)']
    estimate_effect <- summary(lm_inter)$coefficients['Group','Estimate']
    lminter_res[k,]<-c(colnames(simulated_data)[k],pval,estimate_effect)
  }
  p.adj.val  = p.adjust(lminter_res$Pvalue,method = 'BH')
  lminter_res$p.adj.val = p.adj.val
  return(lminter_res)
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
    pval <- summary(lmf)$coefficients['Group','Pr(>|t|)']
    estimate_effect <- summary(lmf)$coefficients['Group','Estimate']
    lmf_res[k,]<-c(colnames(simulated_data)[k],pval,estimate_effect)
  }
  p.adj.val  = p.adjust(lmf_res$Pvalue,method = 'BH')
  lmf_res$p.adj.val = p.adj.val
  return(lmf_res)
}
LM_analysis<-function(simulated_data_all,normalization){
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  lm_res <- data.frame(Feature = rep(NA,ncol(simulated_data_all$simulated_data_rela)),Pvalue=rep(NA,ncol(simulated_data_all$simulated_data_rela)),Effect_size = rep(NA,ncol(simulated_data_all$simulated_data_rela)))
  for(k in 1:ncol(simulated_data_all$simulated_data_rela)){
    temp_data <- data.frame(cbind(simulated_data[,k],simulated_data_all$metadata))
    colnames(temp_data)[1]<-'Feature'
    lm_v <- lm(Feature ~ Group,data = temp_data)
    pval <- summary(lm_v)$coefficients['Group','Pr(>|t|)']
    estimate_effect <- summary(lm_v)$coefficients['Group','Estimate']
    lm_res[k,]<-c(colnames(simulated_data)[k],pval,estimate_effect)
  }
  p.adj.val  = p.adjust(lm_res$Pvalue,method = 'BH')
  lm_res$p.adj.val = p.adj.val
  return(lm_res)
}
VTwins_analysis<-function(simulated_data_all,normalization){
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  library(VTwins)
  library(dplyr)
  pheno_data <- data.frame(id = simulated_data_all$metadata$SampleID,grp = simulated_data_all$metadata$Group)
  rownames(pheno_data) <- pheno_data$id
  
  pheno_data[,2]<-gsub(pheno_data[,2],pattern = 1,replacement = "grp2")
  pheno_data[,2]<-gsub(pheno_data[,2],pattern = 0,replacement = "grp1")
  vtwins_res <- pair_find(data=data.frame(simulated_data),
                          phenodata=pheno_data[rownames(simulated_data_all$simulated_data_rela),],
                          k="euclidean",
                          Cut_pair=25, 
                          method_choose="Permutation",
                          SavePath = "./",
                          ShuffleWstat = "ShuffleWstat", 
                          BoundarySample = "BoundarySample",
                          BoundaryPair = "BoundaryPair",
                          ShuffleTime=1000,
                          DownPercent = 0.2,
                          Uppercent=0.8,PvalueCutoff=0.05)
  if(is.null(vtwins_res)){
    vtwins_res <- data.frame(Feature = colnames(simulated_data_all$simulated_data_rela),
                             p.adj.val = rep(NA,ncol(simulated_data_all$simulated_data_rela)))
  }else{
    for(k in 1:nrow(vtwins_res)){
      vtwins_res$p.adj.val[k] <- ifelse(vtwins_res$Decre.aveRank.P.FDR[k] < vtwins_res$Incre.aveRank.P.FDR[k],vtwins_res$Decre.aveRank.P.FDR[k],vtwins_res$Incre.aveRank.P.FDR[k])
    }
    colnames(vtwins_res)[1] <- 'Feature'
  }
  return(vtwins_res)
}
# 函数：根据中位数划分并创建分类变量
median_discretize <- function(data, continuous_vars) {
  for (var in continuous_vars) {
    if (!var %in% names(data)) {
      stop(paste("Variable", var, "not found in the dataset."))
    }
    
    # 计算中位数并划分
    median_value <- median(data[[var]], na.rm = TRUE)
    data <- data %>%
      mutate(
        !!paste0(var, "_cat") := ifelse(
          !!sym(var) <= median_value, "low", "high"
        )
      )
  }
  return(data)
}

# 函数：合并分类变量
merge_categorical_and_continuous <- function(data, continuous_vars) {
  # 对连续变量离散化
  data <- as.data.frame(data)
  if(is.null(continuous_vars)){
    # 找到以 "Ca" 开头的列
    ca_columns <- grep("^Ca", colnames(data))
    # 对这些列逐行拼接
    data$Ca <- apply(data[, ca_columns], 1, paste0, collapse = "")
    return(data)
  }else{
    data <- median_discretize(data, continuous_vars)
    # 动态生成分类变量名
    continuous_cats <- paste0(continuous_vars, "_cat")
    # 找到以 "Ca" 开头的列和分类变量列
    ca_columns <- grep("^Ca", colnames(data))
    cat_columns <- which(colnames(data) %in% continuous_cats)
    selected_columns <- c(ca_columns, cat_columns)
    
    # 对这些列逐行拼接
    data$Ca <- apply(data[, selected_columns], 1, paste0, collapse = "")
    return(data)
  }
  
}

Wilcoxon_analysis <- function(simulated_data_all, covariates, normalization) {
  library(coin)
  library(dplyr)
  
  simulated_data <- data_normalization(simulated_data_all, normalization = normalization)
  wilcox_res <- data.frame(
    Feature = rep(NA, ncol(simulated_data_all$simulated_data_rela)), 
    Pvalue = rep(NA, ncol(simulated_data_all$simulated_data_rela)), 
    Effect_size = rep(NA, ncol(simulated_data_all$simulated_data_rela))
  )
  
  if (is.null(covariates)) {
    for (k in 1:ncol(simulated_data_all$simulated_data_rela)) {
      temp_data <- data.frame(cbind(simulated_data[, k], simulated_data_all$metadata))
      colnames(temp_data)[1] <- 'Feature'
      x <- temp_data[which(temp_data$Group == 1), 'Feature']
      y <- temp_data[which(temp_data$Group == 0), 'Feature']
      lfc <- log2(mean(as.numeric(x)) / mean(as.numeric(y)))
      wilcox <- wilcox.test(x = x, y = y)
      pval <- wilcox$p.value
      wilcox_res[k, ] <- c(colnames(simulated_data_all$simulated_data_rela)[k], pval, lfc)
    }
    p.adj.val <- p.adjust(wilcox_res$Pvalue, method = 'BH')
    wilcox_res$p.adj.val <- p.adj.val
    return(wilcox_res)
  } else {
    continuous_vars <- colnames(simulated_data_all$metadata)[grep('^Con', colnames(simulated_data_all$metadata))]
    simulated_data_all$metadata <- merge_categorical_and_continuous(simulated_data_all$metadata, continuous_vars)
    formula_str <- "Feature ~ Group | Ca"
    
    for (k in 1:ncol(simulated_data_all$simulated_data_rela)) {
      tryCatch({
        temp_data <- data.frame(cbind(simulated_data[, k], simulated_data_all$metadata))
        colnames(temp_data)[1] <- 'Feature'
        temp_data$Group <- as.factor(temp_data$Group)
        temp_data$Ca <- as.factor(temp_data$Ca)
        
        # 过滤掉 `Ca` 变量水平小于 2 的样本
        ca_counts <- table(temp_data$Ca)
        valid_ca_levels <- names(ca_counts[ca_counts >= 2])
        temp_data <- temp_data[temp_data$Ca %in% valid_ca_levels, ]

        # 再次检查 `Ca` 是否还有足够的水平
        if (length(unique(temp_data$Ca)) < 2) {
          stop("Too few Ca levels after filtering")
        }
        
        wilcox <- wilcox_test(as.formula(formula_str), data = temp_data)
        pval <- pvalue(wilcox)
        x <- temp_data[which(temp_data$Group == 1), 'Feature']
        y <- temp_data[which(temp_data$Group == 0), 'Feature']
        lfc <- log2(mean(as.numeric(x)) / mean(as.numeric(y)))
        wilcox_res[k, ] <- c(colnames(simulated_data_all$simulated_data_rela)[k], pval, lfc)
        
      }, error = function(e) {
        cat(sprintf("block var failed: %s\n", e$message))
        
        if (grepl("less than two observations", e$message)) {
          print("Removing samples with rare `Ca` levels")
          
          # 重新计算 `Ca` 变量的分布，删除 `Ca` 过于稀少的样本
          ca_counts <- table(temp_data$Ca)
          valid_ca_levels <- names(ca_counts[ca_counts >= 2])
          temp_data <- temp_data[temp_data$Ca %in% valid_ca_levels, ]
          
          # 如果仍然不满足 Wilcoxon 统计要求，则跳过
          if (length(unique(temp_data$Ca)) < 2) {
            print("Skipping feature due to insufficient Ca levels")
            next
          }
          
          # 重新执行 Wilcoxon 检验
          wilcox <- wilcox_test(as.formula(formula_str), data = temp_data)
          pval <- pvalue(wilcox)
          x <- temp_data[which(temp_data$Group == 1), 'Feature']
          y <- temp_data[which(temp_data$Group == 0), 'Feature']
          lfc <- log2(mean(as.numeric(x)) / mean(as.numeric(y)))
          wilcox_res[k, ] <- c(colnames(simulated_data_all$simulated_data_rela)[k], pval, lfc)
          
        } else {
          stop(e)  # 其他错误直接抛出
        }
      })
    }
    
    p.adj.val <- p.adjust(wilcox_res$Pvalue, method = 'BH')
    wilcox_res$p.adj.val <- p.adj.val
    return(wilcox_res)
  }
}

# 筛选函数：按 prevalence 筛选矩阵
filter_by_prevalence <- function(matrix, prevalence) {
  # 计算每列的非零比例
  nonzero_proportion <- colSums(matrix != 0) / nrow(matrix)
  
  # 筛选出非零比例大于 prevalence 的列
  filtered_matrix <- matrix[, which(nonzero_proportion > prevalence), drop = FALSE]
  
  # 返回筛选后的矩阵
  return(filtered_matrix)
}
merge_all_method_result<-function(simulated_data_all,covariates=NULL,categorical_variable_name=NULL,prevalence = 0.1,data_index){
  simulated_data_all$simulated_data_rela <- filter_by_prevalence(simulated_data_all$simulated_data_rela,prevalence)
  simulated_data_all$simulated_data_abs <- filter_by_prevalence(simulated_data_all$simulated_data_abs,prevalence)
  print(colSums(simulated_data_all$simulated_data_rela != 0))
  #####################################CLR
  print('CLR start!')
  normalization_method = 'CLR'
  #ANCOM_res_name <- paste0("ANCOM_res_",normalization_method)
  #assign(ANCOM_res_name,ANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  #write.table(get(ANCOM_res_name),file = paste0('ANCOM_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  #print('ANCOM FINISHED')
  maaslin2_lm_res_name <- paste0("maaslin2_lm_res_",data_index,"_",normalization_method)
  assign(maaslin2_lm_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_lm_res_name),file = paste0('maaslin2_lm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm FINISHED')
  #maaslin2_lm_LOG_res_name <- paste0("maaslin2_lm_LOG_res_",normalization_method)
  #assign(maaslin2_lm_LOG_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'LOG'))
  #write.table(get(maaslin2_lm_LOG_res_name),file = paste0('maaslin2_lm_LOG_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  #print('maaslin2-lm LOG FINISHED')
  limma_res_name <- paste0("limma_res_",normalization_method)
  assign(limma_res_name,limma_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(limma_res_name),file = paste0('limma_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('limma FINISHED')
  wilcoxon_res_name <- paste0("wilcoxon_res_",normalization_method)
  assign(wilcoxon_res_name,Wilcoxon_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(wilcoxon_res_name),file = paste0('Wilcoxon_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Wilcoxon test FINISHED')
  #ZicoSeq_res_name <- paste0("ZicoSeq_res_",normalization_method)
  #assign(ZicoSeq_res_name,ZicoSeq_analysis(simulated_data_all,covariates,normalization_method))
  #write.table(get(ZicoSeq_res_name),file = paste0('ZicoSeq_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  #print('ZicoSeq FINISHED')
  if(is.null(covariates)){
    lm_res_name <- paste0("lm_res_",normalization_method)
    assign(lm_res_name,LM_analysis(simulated_data_all,normalization_method))
    write.table(get(lm_res_name),file = paste0('LM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM FINISHED')
    all_CLR_result<-list(get(maaslin2_lm_res_name),get(limma_res_name),get(lm_res_name),get(wilcoxon_res_name))
  }else{
    vtwins_res_name <- paste0("vtwins_res_",normalization_method)
    assign(vtwins_res_name,VTwins_analysis(simulated_data_all,normalization_method))
    write.table(get(vtwins_res_name),file = paste0('VTwins_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('VTwins FINISHED')
    lmem_res_name <- paste0("lmem_res_",normalization_method)
    assign(lmem_res_name,LM_random_effect(simulated_data_all,covariates,categorical_variable_name,normalization_method))
    write.table(get(lmem_res_name),file = paste0('lmem_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LMEM FINISHED')
    lm_inter_res_name <- paste0("lm_inter_res_",normalization_method)
    assign(lm_inter_res_name,LM_interaction_analysis(simulated_data_all,covariates,normalization_method))
    write.table(get(lm_inter_res_name),file = paste0('LM-inter_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-interaction FINISHED')
    lfem_res_name <- paste0("lfem_res_",normalization_method)
    assign(lfem_res_name,LM_fixed_effect(simulated_data_all,covariates,normalization_method))
    write.table(get(lfem_res_name),file = paste0('LM-fixed_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-fixed effect FINISHED')
    #all_result<-list(fa_res,ancom_res,ancombc_res,limma_res,lm_re_res,
    #                 lm_inter_res,lm_fe_res,vtwins_res,wilcoxon_res)
    all_CLR_result<-list(get(ANCOM_res_name),get(maaslin2_lm_res_name),get(limma_res_name),get(wilcoxon_res_name),
      get(lmem_res_name),get(lm_inter_res_name),get(lfem_res_name),get(vtwins_res_name))
  }
  ############################################CSS
  print('CSS start!')
  normalization_method = 'CSS'
  fa_res_name <- paste0("fastANCOM_res_",normalization_method)
  assign(fa_res_name,fastANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(fa_res_name),file = paste0('fastANCOM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  ANCOM_res_name <- paste0("ANCOM_res_",normalization_method)
  assign(ANCOM_res_name,ANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(ANCOM_res_name),file = paste0('ANCOM_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ANCOM FINISHED')
  ANCOMBC_res_name <- paste0("ANCOMBC_res_",normalization_method)
  assign(ANCOMBC_res_name,ANCOMBC_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(ANCOMBC_res_name),file = paste0('ANCOMBC_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ANCOMBC FINISHED')
  maaslin2_lm_res_name <- paste0("maaslin2_lm_res_",data_index,"_",normalization_method)
  assign(maaslin2_lm_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_lm_res_name),file = paste0('maaslin2_lm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm FINISHED')
  maaslin2_lm_LOG_res_name <- paste0("maaslin2_lm_LOG_res_",normalization_method)
  assign(maaslin2_lm_LOG_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'LOG'))
  write.table(get(maaslin2_lm_LOG_res_name),file = paste0('maaslin2_lm_LOG_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm LOG FINISHED')
  #maaslin2_cplm_res_name <- paste0("maaslin2_cplm_res_",normalization_method)
  #assign(maaslin2_cplm_res_name,Maaslin2_CPLM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  #write.table(get(maaslin2_cplm_res_name),file = paste0('maaslin2_cplm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  #print('maaslin2-cplm FINISHED')
  maaslin2_cplm_LOG_res_name <- paste0("maaslin2_cplm_LOG_res_",normalization_method)
  assign(maaslin2_cplm_LOG_res_name,Maaslin2_CPLM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'LOG'))
  write.table(get(maaslin2_cplm_LOG_res_name),file = paste0('maaslin2_cplm_LOG_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-cplm LOG FINISHED')
  maaslin2_negbin_res_name <- paste0("maaslin2_negbin_res_",normalization_method)
  assign(maaslin2_negbin_res_name,Maaslin2_NEGBIN_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_negbin_res_name),file = paste0('maaslin2_negbin_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-negbin FINISHED')
  limma_res_name <- paste0("limma_res_",normalization_method)
  assign(limma_res_name,limma_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(limma_res_name),file = paste0('limma_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('limma FINISHED')
  wilcoxon_res_name <- paste0("wilcoxon_res_",normalization_method)
  assign(wilcoxon_res_name,Wilcoxon_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(wilcoxon_res_name),file = paste0('Wilcoxon_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Wilcoxon test FINISHED')
  ZicoSeq_res_name <- paste0("ZicoSeq_res_",normalization_method)
  assign(ZicoSeq_res_name,ZicoSeq_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(ZicoSeq_res_name),file = paste0('ZicoSeq_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ZicoSeq FINISHED')
  edgeR_res_name <- paste0("edgeR_res_",normalization_method)
  assign(edgeR_res_name,edgeR_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(edgeR_res_name),file = paste0('edgeR_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('edgeR FINISHED')
  if(is.null(covariates)){
    lm_res_name <- paste0("lm_res_",normalization_method)
    assign(lm_res_name,LM_analysis(simulated_data_all,normalization_method))
    write.table(get(lm_res_name),file = paste0('LM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM FINISHED')
    all_CSS_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(fa_res_name),get(ANCOM_res_name),get(ANCOMBC_res_name),get(maaslin2_cplm_LOG_res_name),get(maaslin2_negbin_res_name),
                         get(maaslin2_lm_res_name),get(maaslin2_lm_LOG_res_name),get(limma_res_name),get(lm_res_name),get(wilcoxon_res_name))
  }else{
    lmem_res_name <- paste0("lmem_res_",normalization_method)
    assign(lmem_res_name,LM_random_effect(simulated_data_all,covariates,categorical_variable_name,normalization_method))
    write.table(get(lmem_res_name),file = paste0('lmem_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LMEM FINISHED')
    lm_inter_res_name <- paste0("lm_inter_res_",normalization_method)
    assign(lm_inter_res_name,LM_interaction_analysis(simulated_data_all,covariates,normalization_method))
    write.table(get(lm_inter_res_name),file = paste0('LM-inter_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-interaction FINISHED')
    lfem_res_name <- paste0("lfem_res_",normalization_method)
    assign(lfem_res_name,LM_fixed_effect(simulated_data_all,covariates,normalization_method))
    write.table(get(lfem_res_name),file = paste0('LM-fixed_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-fixed effect FINISHED')
    #all_result<-list(fa_res,ancom_res,ancombc_res,limma_res,lm_re_res,
    #                 lm_inter_res,lm_fe_res,vtwins_res,wilcoxon_res)
    all_CSS_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(fa_res_name),get(ANCOM_res_name),get(ANCOMBC_res_name),get(maaslin2_cplm_res_name),get(maaslin2_cplm_LOG_res_name),get(maaslin2_negbin_res_name),get(lmem_res_name),get(lm_inter_res_name),get(lfem_res_name),get(vtwins_res_name),
         get(maaslin2_lm_res_name),get(maaslin2_lm_LOG_res_name),get(limma_res_name),get(wilcoxon_res_name))
  }
  ####################################count
  print('Count start!')
  normalization_method = 'count'
  fa_res_name <- paste0("fastANCOM_res_",normalization_method)
  assign(fa_res_name,fastANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(fa_res_name),file = paste0('fastANCOM_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('fastANCOM FINISHED')
  ANCOM_res_name <- paste0("ANCOM_res_",normalization_method)
  assign(ANCOM_res_name,ANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(ANCOM_res_name),file = paste0('ANCOM_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ANCOM FINISHED')
  ANCOMBC_res_name <- paste0("ANCOMBC_res_",normalization_method)
  assign(ANCOMBC_res_name,ANCOMBC_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(ANCOMBC_res_name),file = paste0('ANCOMBC_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ANCOMBC FINISHED')
  maaslin2_lm_res_name <- paste0("maaslin2_lm_res_",data_index,"_",normalization_method)
  assign(maaslin2_lm_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization='None',transform = 'None'))
  write.table(get(maaslin2_lm_res_name),file = paste0('maaslin2_lm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm FINISHED')
  #maaslin2_cplm_res_name <- paste0("maaslin2_cplm_res_",normalization_method)
  #assign(maaslin2_cplm_res_name,Maaslin2_CPLM_analysis(simulated_data_all,covariates,normalization='None',transform = 'None'))
  #write.table(get(maaslin2_cplm_res_name),file = paste0('maaslin2_cplm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  #print('maaslin2-cplm FINISHED')
  maaslin2_zinb_res_name <- paste0("maaslin2_zinb_res_",normalization_method)
  assign(maaslin2_zinb_res_name,Maaslin2_ZINB_analysis(simulated_data_all,covariates,normalization='None',transform = 'None'))
  write.table(get(maaslin2_zinb_res_name),file = paste0('maaslin2_zinb_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-zinb FINISHED')
  maaslin2_negbin_res_name <- paste0("maaslin2_negbin_res_",normalization_method)
  assign(maaslin2_negbin_res_name,Maaslin2_NEGBIN_analysis(simulated_data_all,covariates,normalization='None',transform = 'None'))
  write.table(get(maaslin2_negbin_res_name),file = paste0('maaslin2_megbin_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-negbin FINISHED')
  limma_res_name <- paste0("limma_res_",normalization_method)
  assign(limma_res_name,limma_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(limma_res_name),file = paste0('limma_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('limma FINISHED')
  wilcoxon_res_name <- paste0("wilcoxon_res_",normalization_method)
  assign(wilcoxon_res_name,Wilcoxon_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(wilcoxon_res_name),file = paste0('Wilcoxon_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Wilcoxon test FINISHED')
  ZicoSeq_res_name <- paste0("ZicoSeq_res_",normalization_method)
  assign(ZicoSeq_res_name,ZicoSeq_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(ZicoSeq_res_name),file = paste0('ZicoSeq_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ZicoSeq FINISHED')
  edgeR_res_name <- paste0("edgeR_res_",normalization_method)
  assign(edgeR_res_name,edgeR_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(edgeR_res_name),file = paste0('edgeR_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('edgeR FINISHED')
  DESeq2_res_name <- paste0("DESeq2_res_",normalization_method)
  assign(DESeq2_res_name,DESeq2_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(DESeq2_res_name),file = paste0('DESeq2_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('DESeq2 FINISHED')
  ALDEx2_res_name <- paste0("ALDEx2_res_",normalization_method)
  assign(ALDEx2_res_name,ALDEx2_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(ALDEx2_res_name),file = paste0('ALDEx2_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ALDEx2 FINISHED')
  Corncob_res_name <- paste0("Corncob_res_",normalization_method)
  assign(Corncob_res_name,Corncob_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(Corncob_res_name),file = paste0('Corncob_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Corncob FINISHED')
  if(is.null(covariates)){
    lm_res_name <- paste0("lm_res_",normalization_method)
    assign(lm_res_name,LM_analysis(simulated_data_all,normalization_method))
    write.table(get(lm_res_name),file = paste0('LM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM FINISHED')
    all_count_result<-list(get(DESeq2_res_name),get(ALDEx2_res_name),get(Corncob_res_name),get(edgeR_res_name),get(ZicoSeq_res_name),get(ANCOM_res_name),get(ANCOMBC_res_name),get(fa_res_name),get(maaslin2_zinb_res_name),get(maaslin2_negbin_res_name),
                         get(maaslin2_lm_res_name),get(limma_res_name),get(lm_res_name),get(wilcoxon_res_name))
  }else{
    vtwins_res_name <- paste0("vtwins_res_",normalization_method)
    assign(vtwins_res_name,VTwins_analysis(simulated_data_all,normalization_method))
    write.table(get(vtwins_res_name),file = paste0('VTwins_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('VTwins FINISHED')
    lmem_res_name <- paste0("lmem_res_",normalization_method)
    assign(lmem_res_name,LM_random_effect(simulated_data_all,covariates,categorical_variable_name,normalization_method))
    write.table(get(lmem_res_name),file = paste0('lmem_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LMEM FINISHED')
    lm_inter_res_name <- paste0("lm_inter_res_",normalization_method)
    assign(lm_inter_res_name,LM_interaction_analysis(simulated_data_all,covariates,normalization_method))
    write.table(get(lm_inter_res_name),file = paste0('LM-inter_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-interaction FINISHED')
    lfem_res_name <- paste0("lfem_res_",normalization_method)
    assign(lfem_res_name,LM_fixed_effect(simulated_data_all,covariates,normalization_method))
    write.table(get(lfem_res_name),file = paste0('LM-fixed_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-fixed effect FINISHED')
    #all_result<-list(fa_res,ancom_res,ancombc_res,limma_res,lm_re_res,
    #                 lm_inter_res,lm_fe_res,vtwins_res,wilcoxon_res)
    all_count_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(lmem_res_name),get(lm_inter_res_name),get(lfem_res_name),get(vtwins_res_name),get(ANCOM_res_name),get(ANCOMBC_res_name),get(fa_res_name),get(maaslin2_cplm_res_name),get(maaslin2_zinb_res_name),get(maaslin2_negbin_res_name),
                         get(maaslin2_lm_res_name),get(limma_res_name),get(wilcoxon_res_name))
  }
  #####################################TMM
  print('TMM start!')
  normalization_method = 'TMM'
  fa_res_name <- paste0("fastANCOM_res_",normalization_method)
  assign(fa_res_name,fastANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(fa_res_name),file = paste0('fastANCOM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  ANCOM_res_name <- paste0("ANCOM_res_",normalization_method)
  assign(ANCOM_res_name,ANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(ANCOM_res_name),file = paste0('ANCOM_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ANCOM FINISHED')
  ANCOMBC_res_name <- paste0("ANCOMBC_res_",normalization_method)
  assign(ANCOMBC_res_name,ANCOMBC_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(ANCOMBC_res_name),file = paste0('ANCOMBC_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ANCOMBC FINISHED')
  maaslin2_lm_res_name <- paste0("maaslin2_lm_res_",data_index,"_",normalization_method)
  assign(maaslin2_lm_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_lm_res_name),file = paste0('maaslin2_lm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm FINISHED')
  maaslin2_lm_LOG_res_name <- paste0("maaslin2_lm_LOG_res_",normalization_method)
  assign(maaslin2_lm_LOG_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'LOG'))
  write.table(get(maaslin2_lm_LOG_res_name),file = paste0('maaslin2_lm_LOG_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm LOG FINISHED')
  #maaslin2_cplm_res_name <- paste0("maaslin2_cplm_res_",normalization_method)
  #assign(maaslin2_cplm_res_name,Maaslin2_CPLM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  #write.table(get(maaslin2_cplm_res_name),file = paste0('maaslin2_cplm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  #print('maaslin2-cplm FINISHED')
  maaslin2_negbin_res_name <- paste0("maaslin2_negbin_res_",normalization_method)
  assign(maaslin2_negbin_res_name,Maaslin2_NEGBIN_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_negbin_res_name),file = paste0('maaslin2_negbin_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-negbin FINISHED')
  limma_res_name <- paste0("limma_res_",normalization_method)
  assign(limma_res_name,limma_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(limma_res_name),file = paste0('limma_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('limma FINISHED')
  wilcoxon_res_name <- paste0("wilcoxon_res_",normalization_method)
  assign(wilcoxon_res_name,Wilcoxon_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(wilcoxon_res_name),file = paste0('Wilcoxon_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Wilcoxon test FINISHED')
  ZicoSeq_res_name <- paste0("ZicoSeq_res_",normalization_method)
  assign(ZicoSeq_res_name,ZicoSeq_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(ZicoSeq_res_name),file = paste0('ZicoSeq_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ZicoSeq FINISHED')
  edgeR_res_name <- paste0("edgeR_res_",normalization_method)
  assign(edgeR_res_name,edgeR_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(edgeR_res_name),file = paste0('edgeR_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('edgeR FINISHED')
  if(is.null(covariates)){
    lm_res_name <- paste0("lm_res_",normalization_method)
    assign(lm_res_name,LM_analysis(simulated_data_all,normalization_method))
    write.table(get(lm_res_name),file = paste0('LM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM FINISHED')
    all_TMM_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(fa_res_name),get(ANCOM_res_name),get(ANCOMBC_res_name),get(maaslin2_negbin_res_name),
                         get(maaslin2_lm_res_name),get(maaslin2_lm_LOG_res_name),get(limma_res_name),get(lm_res_name),get(wilcoxon_res_name))
  }else{
    vtwins_res_name <- paste0("vtwins_res_",normalization_method)
    assign(vtwins_res_name,VTwins_analysis(simulated_data_all,normalization_method))
    write.table(get(vtwins_res_name),file = paste0('VTwins_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('VTwins FINISHED')
    lmem_res_name <- paste0("lmem_res_",normalization_method)
    assign(lmem_res_name,LM_random_effect(simulated_data_all,covariates,categorical_variable_name,normalization_method))
    write.table(get(lmem_res_name),file = paste0('lmem_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LMEM FINISHED')
    lm_inter_res_name <- paste0("lm_inter_res_",normalization_method)
    assign(lm_inter_res_name,LM_interaction_analysis(simulated_data_all,covariates,normalization_method))
    write.table(get(lm_inter_res_name),file = paste0('LM-inter_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-interaction FINISHED')
    lfem_res_name <- paste0("lfem_res_",normalization_method)
    assign(lfem_res_name,LM_fixed_effect(simulated_data_all,covariates,normalization_method))
    write.table(get(lfem_res_name),file = paste0('LM-fixed_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-fixed effect FINISHED')
    #all_result<-list(fa_res,ancom_res,ancombc_res,limma_res,lm_re_res,
    #                 lm_inter_res,lm_fe_res,vtwins_res,wilcoxon_res)
    all_TMM_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(fa_res_name),get(ANCOM_res_name),get(ANCOMBC_res_name),get(maaslin2_cplm_res_name),get(maaslin2_negbin_res_name),get(lmem_res_name),get(lm_inter_res_name),get(lfem_res_name),get(vtwins_res_name),
                         get(maaslin2_lm_res_name),get(maaslin2_lm_LOG_res_name),get(limma_res_name),get(wilcoxon_res_name))
  }
  ################################TSS
  print('TSS start!')
  normalization_method = 'TSS'
  fa_res_name <- paste0("fastANCOM_res_",normalization_method)
  assign(fa_res_name,fastANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(fa_res_name),file = paste0('fastANCOM_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('fastANCOM FINISHED')
  maaslin2_lm_res_name <- paste0("maaslin2_lm_res_",data_index,"_",normalization_method)
  assign(maaslin2_lm_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_lm_res_name),file = paste0('maaslin2_lm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm FINISHED')
  maaslin2_lm_LOG_res_name <- paste0("maaslin2_lm_LOG_res_",normalization_method)
  assign(maaslin2_lm_LOG_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'LOG'))
  write.table(get(maaslin2_lm_LOG_res_name),file = paste0('maaslin2_lm_LOG_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm LOG FINISHED')
  maaslin2_cplm_res_name <- paste0("maaslin2_cplm_res_",normalization_method)
  assign(maaslin2_cplm_res_name,Maaslin2_CPLM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_cplm_res_name),file = paste0('maaslin2_cplm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-cplm FINISHED')
  maaslin2_cplm_LOG_res_name <- paste0("maaslin2_cplm_LOG_res_",normalization_method)
  assign(maaslin2_cplm_LOG_res_name,Maaslin2_CPLM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_cplm_LOG_res_name),file = paste0('maaslin2_cplm_LOG_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-cplm LOG FINISHED')
  limma_res_name <- paste0("limma_res_",normalization_method)
  assign(limma_res_name,limma_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(limma_res_name),file = paste0('limma_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('limma FINISHED')
  wilcoxon_res_name <- paste0("wilcoxon_res_",normalization_method)
  assign(wilcoxon_res_name,Wilcoxon_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(wilcoxon_res_name),file = paste0('Wilcoxon_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Wilcoxon test FINISHED')
  ZicoSeq_res_name <- paste0("ZicoSeq_res_",normalization_method)
  assign(ZicoSeq_res_name,ZicoSeq_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(ZicoSeq_res_name),file = paste0('ZicoSeq_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ZicoSeq FINISHED')
  edgeR_res_name <- paste0("edgeR_res_",normalization_method)
  assign(edgeR_res_name,edgeR_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(edgeR_res_name),file = paste0('edgeR_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('edgeR FINISHED')
  if(is.null(covariates)){
    lm_res_name <- paste0("lm_res_",normalization_method)
    assign(lm_res_name,LM_analysis(simulated_data_all,normalization_method))
    write.table(get(lm_res_name),file = paste0('LM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM FINISHED')
    all_TSS_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(fa_res_name),get(maaslin2_cplm_res_name),get(maaslin2_cplm_LOG_res_name),
                         get(maaslin2_lm_res_name),get(maaslin2_lm_LOG_res_name),get(limma_res_name),get(lm_res_name),get(wilcoxon_res_name))
  }else{
    vtwins_res_name <- paste0("vtwins_res_",normalization_method)
    assign(vtwins_res_name,VTwins_analysis(simulated_data_all,normalization_method))
    write.table(get(vtwins_res_name),file = paste0('VTwins_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('VTwins FINISHED')
    lmem_res_name <- paste0("lmem_res_",normalization_method)
    assign(lmem_res_name,LM_random_effect(simulated_data_all,covariates,categorical_variable_name,normalization_method))
    write.table(get(lmem_res_name),file = paste0('lmem_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LMEM FINISHED')
    lm_inter_res_name <- paste0("lm_inter_res_",normalization_method)
    assign(lm_inter_res_name,LM_interaction_analysis(simulated_data_all,covariates,normalization_method))
    write.table(get(lm_inter_res_name),file = paste0('LM-inter_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-interaction FINISHED')
    lfem_res_name <- paste0("lfem_res_",normalization_method)
    assign(lfem_res_name,LM_fixed_effect(simulated_data_all,covariates,normalization_method))
    write.table(get(lfem_res_name),file = paste0('LM-fixed_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-fixed effect FINISHED')
    #all_result<-list(fa_res,ancom_res,ancombc_res,limma_res,lm_re_res,
    #                 lm_inter_res,lm_fe_res,vtwins_res,wilcoxon_res)
    all_TSS_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(fa_res_name),get(maaslin2_cplm_res_name),get(maaslin2_cplm_LOG_res_name),get(lmem_res_name),get(lm_inter_res_name),get(lfem_res_name),get(vtwins_res_name),
                         get(maaslin2_lm_res_name),get(maaslin2_lm_LOG_res_name),get(limma_res_name),get(wilcoxon_res_name))
  }
  #########################################LOG
  print('LOG start!')
  normalization_method = 'LOG'
  maaslin2_lm_res_name <- paste0("maaslin2_lm_res_",data_index,"_",normalization_method)
  assign(maaslin2_lm_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization='None',transform = normalization_method))
  write.table(get(maaslin2_lm_res_name),file = paste0('maaslin2_lm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm FINISHED')
  
  limma_res_name <- paste0("limma_res_",normalization_method)
  assign(limma_res_name,limma_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(limma_res_name),file = paste0('limma_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('limma FINISHED')
  wilcoxon_res_name <- paste0("wilcoxon_res_",normalization_method)
  assign(wilcoxon_res_name,Wilcoxon_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(wilcoxon_res_name),file = paste0('Wilcoxon_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Wilcoxon test FINISHED')
  ZicoSeq_res_name <- paste0("ZicoSeq_res_",normalization_method)
  assign(ZicoSeq_res_name,ZicoSeq_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(ZicoSeq_res_name),file = paste0('ZicoSeq_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ZicoSeq FINISHED')
  if(is.null(covariates)){
    lm_res_name <- paste0("lm_res_",normalization_method)
    assign(lm_res_name,LM_analysis(simulated_data_all,normalization_method))
    write.table(get(lm_res_name),file = paste0('LM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM FINISHED')
    all_LOG_result<-list(get(ZicoSeq_res_name),get(maaslin2_lm_res_name),get(limma_res_name),get(lm_res_name),get(wilcoxon_res_name))
  }else{
    vtwins_res_name <- paste0("vtwins_res_",normalization_method)
    assign(vtwins_res_name,VTwins_analysis(simulated_data_all,normalization_method))
    write.table(get(vtwins_res_name),file = paste0('VTwins_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('VTwins FINISHED')
    lmem_res_name <- paste0("lmem_res_",normalization_method)
    assign(lmem_res_name,LM_random_effect(simulated_data_all,covariates,categorical_variable_name,normalization_method))
    write.table(get(lmem_res_name),file = paste0('lmem_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LMEM FINISHED')
    lm_inter_res_name <- paste0("lm_inter_res_",normalization_method)
    assign(lm_inter_res_name,LM_interaction_analysis(simulated_data_all,covariates,normalization_method))
    write.table(get(lm_inter_res_name),file = paste0('LM-inter_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-interaction FINISHED')
    lfem_res_name <- paste0("lfem_res_",normalization_method)
    assign(lfem_res_name,LM_fixed_effect(simulated_data_all,covariates,normalization_method))
    write.table(get(lfem_res_name),file = paste0('LM-fixed_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-fixed effect FINISHED')
    #all_result<-list(fa_res,ancom_res,ancombc_res,limma_res,lm_re_res,
    #                 lm_inter_res,lm_fe_res,vtwins_res,wilcoxon_res)
    all_LOG_result<-list(get(ZicoSeq_res_name),get(maaslin2_lm_res_name),get(limma_res_name),get(wilcoxon_res_name),
                         get(lmem_res_name),get(lm_inter_res_name),get(lfem_res_name),get(vtwins_res_name))
  }
  all_result<-list(all_CLR_result=all_CLR_result,all_CSS_result=all_CSS_result,all_TSS_result=all_TSS_result,
  all_TMM_result=all_TMM_result,all_count_result=all_count_result,all_LOG_result=all_LOG_result)
  
  return(all_result)
}
merge_all_fast_method_result<-function(simulated_data_all,covariates=NULL,categorical_variable_name=NULL,prevalence = 0.1,data_index){
  simulated_data_all$simulated_data_rela <- filter_by_prevalence(simulated_data_all$simulated_data_rela,prevalence)
  simulated_data_all$simulated_data_abs <- filter_by_prevalence(simulated_data_all$simulated_data_abs,prevalence)
  print(colSums(simulated_data_all$simulated_data_rela != 0))
  #####################################CLR
  print('CLR start!')
  normalization_method = 'CLR'
  ANCOM_res_name <- paste0("ANCOM_res_",normalization_method)
  assign(ANCOM_res_name,ANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(ANCOM_res_name),file = paste0('ANCOM_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ANCOM FINISHED')
  maaslin2_lm_res_name <- paste0("maaslin2_lm_res_",data_index,"_",normalization_method)
  assign(maaslin2_lm_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_lm_res_name),file = paste0('maaslin2_lm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm FINISHED')
  #maaslin2_lm_LOG_res_name <- paste0("maaslin2_lm_LOG_res_",normalization_method)
  #assign(maaslin2_lm_LOG_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'LOG'))
  #write.table(get(maaslin2_lm_LOG_res_name),file = paste0('maaslin2_lm_LOG_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  #print('maaslin2-lm LOG FINISHED')
  limma_res_name <- paste0("limma_res_",normalization_method)
  assign(limma_res_name,limma_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(limma_res_name),file = paste0('limma_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('limma FINISHED')
  wilcoxon_res_name <- paste0("wilcoxon_res_",normalization_method)
  assign(wilcoxon_res_name,Wilcoxon_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(wilcoxon_res_name),file = paste0('Wilcoxon_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Wilcoxon test FINISHED')
  #ZicoSeq_res_name <- paste0("ZicoSeq_res_",normalization_method)
  #assign(ZicoSeq_res_name,ZicoSeq_analysis(simulated_data_all,covariates,normalization_method))
  #write.table(get(ZicoSeq_res_name),file = paste0('ZicoSeq_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  #print('ZicoSeq FINISHED')
  if(is.null(covariates)){
    lm_res_name <- paste0("lm_res_",normalization_method)
    assign(lm_res_name,LM_analysis(simulated_data_all,normalization_method))
    write.table(get(lm_res_name),file = paste0('LM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM FINISHED')
    all_CLR_result<-list(get(ANCOM_res_name),get(maaslin2_lm_res_name),get(limma_res_name),get(lm_res_name),get(wilcoxon_res_name))
  }else{
    
    lmem_res_name <- paste0("lmem_res_",normalization_method)
    assign(lmem_res_name,LM_random_effect(simulated_data_all,covariates,categorical_variable_name,normalization_method))
    write.table(get(lmem_res_name),file = paste0('lmem_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LMEM FINISHED')
    lm_inter_res_name <- paste0("lm_inter_res_",normalization_method)
    assign(lm_inter_res_name,LM_interaction_analysis(simulated_data_all,covariates,normalization_method))
    write.table(get(lm_inter_res_name),file = paste0('LM-inter_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-interaction FINISHED')
    lfem_res_name <- paste0("lfem_res_",normalization_method)
    assign(lfem_res_name,LM_fixed_effect(simulated_data_all,covariates,normalization_method))
    write.table(get(lfem_res_name),file = paste0('LM-fixed_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-fixed effect FINISHED')
    #all_result<-list(fa_res,ancom_res,ancombc_res,limma_res,lm_re_res,
    #                 lm_inter_res,lm_fe_res,vtwins_res,wilcoxon_res)
    all_CLR_result<-list(get(ANCOM_res_name),get(maaslin2_lm_res_name),get(limma_res_name),get(wilcoxon_res_name),
      get(lmem_res_name),get(lm_inter_res_name),get(lfem_res_name))
  }
  ############################################CSS
  print('CSS start!')
  normalization_method = 'CSS'
  fa_res_name <- paste0("fastANCOM_res_",normalization_method)
  assign(fa_res_name,fastANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(fa_res_name),file = paste0('fastANCOM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  ANCOM_res_name <- paste0("ANCOM_res_",normalization_method)
  assign(ANCOM_res_name,ANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(ANCOM_res_name),file = paste0('ANCOM_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ANCOM FINISHED')
  ANCOMBC_res_name <- paste0("ANCOMBC_res_",normalization_method)
  assign(ANCOMBC_res_name,ANCOMBC_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(ANCOMBC_res_name),file = paste0('ANCOMBC_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ANCOMBC FINISHED')
  maaslin2_lm_res_name <- paste0("maaslin2_lm_res_",data_index,"_",normalization_method)
  assign(maaslin2_lm_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_lm_res_name),file = paste0('maaslin2_lm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm FINISHED')
  maaslin2_lm_LOG_res_name <- paste0("maaslin2_lm_LOG_res_",normalization_method)
  assign(maaslin2_lm_LOG_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'LOG'))
  write.table(get(maaslin2_lm_LOG_res_name),file = paste0('maaslin2_lm_LOG_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm LOG FINISHED')
  maaslin2_cplm_res_name <- paste0("maaslin2_cplm_res_",normalization_method)
  assign(maaslin2_cplm_res_name,Maaslin2_CPLM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_cplm_res_name),file = paste0('maaslin2_cplm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-cplm FINISHED')
  maaslin2_cplm_LOG_res_name <- paste0("maaslin2_cplm_LOG_res_",normalization_method)
  assign(maaslin2_cplm_LOG_res_name,Maaslin2_CPLM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'LOG'))
  write.table(get(maaslin2_cplm_LOG_res_name),file = paste0('maaslin2_cplm_LOG_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-cplm LOG FINISHED')
  maaslin2_negbin_res_name <- paste0("maaslin2_negbin_res_",normalization_method)
  assign(maaslin2_negbin_res_name,Maaslin2_NEGBIN_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_negbin_res_name),file = paste0('maaslin2_negbin_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-negbin FINISHED')
  limma_res_name <- paste0("limma_res_",normalization_method)
  assign(limma_res_name,limma_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(limma_res_name),file = paste0('limma_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('limma FINISHED')
  wilcoxon_res_name <- paste0("wilcoxon_res_",normalization_method)
  assign(wilcoxon_res_name,Wilcoxon_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(wilcoxon_res_name),file = paste0('Wilcoxon_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Wilcoxon test FINISHED')
  ZicoSeq_res_name <- paste0("ZicoSeq_res_",normalization_method)
  assign(ZicoSeq_res_name,ZicoSeq_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(ZicoSeq_res_name),file = paste0('ZicoSeq_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ZicoSeq FINISHED')
  edgeR_res_name <- paste0("edgeR_res_",normalization_method)
  assign(edgeR_res_name,edgeR_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(edgeR_res_name),file = paste0('edgeR_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('edgeR FINISHED')
  if(is.null(covariates)){
    lm_res_name <- paste0("lm_res_",normalization_method)
    assign(lm_res_name,LM_analysis(simulated_data_all,normalization_method))
    write.table(get(lm_res_name),file = paste0('LM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM FINISHED')
    all_CSS_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(fa_res_name),get(ANCOM_res_name),get(ANCOMBC_res_name),get(maaslin2_cplm_res_name),get(maaslin2_cplm_LOG_res_name),get(maaslin2_negbin_res_name),
                         get(maaslin2_lm_res_name),get(maaslin2_lm_LOG_res_name),get(limma_res_name),get(lm_res_name),get(wilcoxon_res_name))
  }else{
    lmem_res_name <- paste0("lmem_res_",normalization_method)
    assign(lmem_res_name,LM_random_effect(simulated_data_all,covariates,categorical_variable_name,normalization_method))
    write.table(get(lmem_res_name),file = paste0('lmem_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LMEM FINISHED')
    lm_inter_res_name <- paste0("lm_inter_res_",normalization_method)
    assign(lm_inter_res_name,LM_interaction_analysis(simulated_data_all,covariates,normalization_method))
    write.table(get(lm_inter_res_name),file = paste0('LM-inter_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-interaction FINISHED')
    lfem_res_name <- paste0("lfem_res_",normalization_method)
    assign(lfem_res_name,LM_fixed_effect(simulated_data_all,covariates,normalization_method))
    write.table(get(lfem_res_name),file = paste0('LM-fixed_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-fixed effect FINISHED')
    #all_result<-list(fa_res,ancom_res,ancombc_res,limma_res,lm_re_res,
    #                 lm_inter_res,lm_fe_res,vtwins_res,wilcoxon_res)
    all_CSS_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(fa_res_name),get(ANCOM_res_name),get(ANCOMBC_res_name),get(maaslin2_cplm_res_name),get(maaslin2_cplm_LOG_res_name),get(maaslin2_negbin_res_name),get(lmem_res_name),get(lm_inter_res_name),get(lfem_res_name),
         get(maaslin2_lm_res_name),get(maaslin2_lm_LOG_res_name),get(limma_res_name),get(wilcoxon_res_name))
  }
  ####################################count
  print('Count start!')
  normalization_method = 'count'
  fa_res_name <- paste0("fastANCOM_res_",normalization_method)
  assign(fa_res_name,fastANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(fa_res_name),file = paste0('fastANCOM_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('fastANCOM FINISHED')
  ANCOM_res_name <- paste0("ANCOM_res_",normalization_method)
  assign(ANCOM_res_name,ANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(ANCOM_res_name),file = paste0('ANCOM_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ANCOM FINISHED')
  ANCOMBC_res_name <- paste0("ANCOMBC_res_",normalization_method)
  assign(ANCOMBC_res_name,ANCOMBC_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(ANCOMBC_res_name),file = paste0('ANCOMBC_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ANCOMBC FINISHED')
  maaslin2_lm_res_name <- paste0("maaslin2_lm_res_",data_index,"_",normalization_method)
  assign(maaslin2_lm_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization='None',transform = 'None'))
  write.table(get(maaslin2_lm_res_name),file = paste0('maaslin2_lm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm FINISHED')
  maaslin2_cplm_res_name <- paste0("maaslin2_cplm_res_",normalization_method)
  assign(maaslin2_cplm_res_name,Maaslin2_CPLM_analysis(simulated_data_all,covariates,normalization='None',transform = 'None'))
  write.table(get(maaslin2_cplm_res_name),file = paste0('maaslin2_cplm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-cplm FINISHED')
  maaslin2_zinb_res_name <- paste0("maaslin2_zinb_res_",normalization_method)
  assign(maaslin2_zinb_res_name,Maaslin2_ZINB_analysis(simulated_data_all,covariates,normalization='None',transform = 'None'))
  write.table(get(maaslin2_zinb_res_name),file = paste0('maaslin2_zinb_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-zinb FINISHED')
  maaslin2_negbin_res_name <- paste0("maaslin2_negbin_res_",normalization_method)
  assign(maaslin2_negbin_res_name,Maaslin2_NEGBIN_analysis(simulated_data_all,covariates,normalization='None',transform = 'None'))
  write.table(get(maaslin2_negbin_res_name),file = paste0('maaslin2_megbin_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-negbin FINISHED')
  limma_res_name <- paste0("limma_res_",normalization_method)
  assign(limma_res_name,limma_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(limma_res_name),file = paste0('limma_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('limma FINISHED')
  wilcoxon_res_name <- paste0("wilcoxon_res_",normalization_method)
  assign(wilcoxon_res_name,Wilcoxon_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(wilcoxon_res_name),file = paste0('Wilcoxon_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Wilcoxon test FINISHED')
  ZicoSeq_res_name <- paste0("ZicoSeq_res_",normalization_method)
  assign(ZicoSeq_res_name,ZicoSeq_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(ZicoSeq_res_name),file = paste0('ZicoSeq_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ZicoSeq FINISHED')
  edgeR_res_name <- paste0("edgeR_res_",normalization_method)
  assign(edgeR_res_name,edgeR_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(edgeR_res_name),file = paste0('edgeR_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('edgeR FINISHED')
  DESeq2_res_name <- paste0("DESeq2_res_",normalization_method)
  assign(DESeq2_res_name,DESeq2_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(DESeq2_res_name),file = paste0('DESeq2_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('DESeq2 FINISHED')
  ALDEx2_res_name <- paste0("ALDEx2_res_",normalization_method)
  assign(ALDEx2_res_name,ALDEx2_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(ALDEx2_res_name),file = paste0('ALDEx2_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ALDEx2 FINISHED')
  Corncob_res_name <- paste0("Corncob_res_",normalization_method)
  assign(Corncob_res_name,Corncob_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(Corncob_res_name),file = paste0('Corncob_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Corncob FINISHED')
  if(is.null(covariates)){
    lm_res_name <- paste0("lm_res_",normalization_method)
    assign(lm_res_name,LM_analysis(simulated_data_all,normalization_method))
    write.table(get(lm_res_name),file = paste0('LM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM FINISHED')
    all_count_result<-list(get(DESeq2_res_name),get(ALDEx2_res_name),get(Corncob_res_name),get(edgeR_res_name),get(ZicoSeq_res_name),get(ANCOM_res_name),get(ANCOMBC_res_name),get(fa_res_name),get(maaslin2_cplm_res_name),get(maaslin2_zinb_res_name),get(maaslin2_negbin_res_name),
                         get(maaslin2_lm_res_name),get(limma_res_name),get(lm_res_name),get(wilcoxon_res_name))
  }else{
    
    lmem_res_name <- paste0("lmem_res_",normalization_method)
    assign(lmem_res_name,LM_random_effect(simulated_data_all,covariates,categorical_variable_name,normalization_method))
    write.table(get(lmem_res_name),file = paste0('lmem_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LMEM FINISHED')
    lm_inter_res_name <- paste0("lm_inter_res_",normalization_method)
    assign(lm_inter_res_name,LM_interaction_analysis(simulated_data_all,covariates,normalization_method))
    write.table(get(lm_inter_res_name),file = paste0('LM-inter_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-interaction FINISHED')
    lfem_res_name <- paste0("lfem_res_",normalization_method)
    assign(lfem_res_name,LM_fixed_effect(simulated_data_all,covariates,normalization_method))
    write.table(get(lfem_res_name),file = paste0('LM-fixed_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-fixed effect FINISHED')
    #all_result<-list(fa_res,ancom_res,ancombc_res,limma_res,lm_re_res,
    #                 lm_inter_res,lm_fe_res,vtwins_res,wilcoxon_res)
    all_count_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(lmem_res_name),get(lm_inter_res_name),get(lfem_res_name),get(ANCOM_res_name),get(ANCOMBC_res_name),get(fa_res_name),get(maaslin2_cplm_res_name),get(maaslin2_zinb_res_name),get(maaslin2_negbin_res_name),
                         get(maaslin2_lm_res_name),get(limma_res_name),get(wilcoxon_res_name))
  }
  #####################################TMM
  print('TMM start!')
  normalization_method = 'TMM'
  fa_res_name <- paste0("fastANCOM_res_",normalization_method)
  assign(fa_res_name,fastANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(fa_res_name),file = paste0('fastANCOM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  ANCOM_res_name <- paste0("ANCOM_res_",normalization_method)
  assign(ANCOM_res_name,ANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(ANCOM_res_name),file = paste0('ANCOM_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ANCOM FINISHED')
  ANCOMBC_res_name <- paste0("ANCOMBC_res_",normalization_method)
  assign(ANCOMBC_res_name,ANCOMBC_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(ANCOMBC_res_name),file = paste0('ANCOMBC_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ANCOMBC FINISHED')
  maaslin2_lm_res_name <- paste0("maaslin2_lm_res_",data_index,"_",normalization_method)
  assign(maaslin2_lm_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_lm_res_name),file = paste0('maaslin2_lm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm FINISHED')
  maaslin2_lm_LOG_res_name <- paste0("maaslin2_lm_LOG_res_",normalization_method)
  assign(maaslin2_lm_LOG_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'LOG'))
  write.table(get(maaslin2_lm_LOG_res_name),file = paste0('maaslin2_lm_LOG_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm LOG FINISHED')
  maaslin2_cplm_res_name <- paste0("maaslin2_cplm_res_",normalization_method)
  assign(maaslin2_cplm_res_name,Maaslin2_CPLM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_cplm_res_name),file = paste0('maaslin2_cplm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-cplm FINISHED')
  maaslin2_negbin_res_name <- paste0("maaslin2_negbin_res_",normalization_method)
  assign(maaslin2_negbin_res_name,Maaslin2_NEGBIN_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_negbin_res_name),file = paste0('maaslin2_negbin_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-negbin FINISHED')
  limma_res_name <- paste0("limma_res_",normalization_method)
  assign(limma_res_name,limma_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(limma_res_name),file = paste0('limma_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('limma FINISHED')
  wilcoxon_res_name <- paste0("wilcoxon_res_",normalization_method)
  assign(wilcoxon_res_name,Wilcoxon_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(wilcoxon_res_name),file = paste0('Wilcoxon_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Wilcoxon test FINISHED')
  ZicoSeq_res_name <- paste0("ZicoSeq_res_",normalization_method)
  assign(ZicoSeq_res_name,ZicoSeq_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(ZicoSeq_res_name),file = paste0('ZicoSeq_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ZicoSeq FINISHED')
  edgeR_res_name <- paste0("edgeR_res_",normalization_method)
  assign(edgeR_res_name,edgeR_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(edgeR_res_name),file = paste0('edgeR_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('edgeR FINISHED')
  if(is.null(covariates)){
    lm_res_name <- paste0("lm_res_",normalization_method)
    assign(lm_res_name,LM_analysis(simulated_data_all,normalization_method))
    write.table(get(lm_res_name),file = paste0('LM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM FINISHED')
    all_TMM_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(fa_res_name),get(ANCOM_res_name),get(ANCOMBC_res_name),get(maaslin2_cplm_res_name),get(maaslin2_negbin_res_name),
                         get(maaslin2_lm_res_name),get(maaslin2_lm_LOG_res_name),get(limma_res_name),get(lm_res_name),get(wilcoxon_res_name))
  }else{
    
    lmem_res_name <- paste0("lmem_res_",normalization_method)
    assign(lmem_res_name,LM_random_effect(simulated_data_all,covariates,categorical_variable_name,normalization_method))
    write.table(get(lmem_res_name),file = paste0('lmem_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LMEM FINISHED')
    lm_inter_res_name <- paste0("lm_inter_res_",normalization_method)
    assign(lm_inter_res_name,LM_interaction_analysis(simulated_data_all,covariates,normalization_method))
    write.table(get(lm_inter_res_name),file = paste0('LM-inter_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-interaction FINISHED')
    lfem_res_name <- paste0("lfem_res_",normalization_method)
    assign(lfem_res_name,LM_fixed_effect(simulated_data_all,covariates,normalization_method))
    write.table(get(lfem_res_name),file = paste0('LM-fixed_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-fixed effect FINISHED')
    #all_result<-list(fa_res,ancom_res,ancombc_res,limma_res,lm_re_res,
    #                 lm_inter_res,lm_fe_res,vtwins_res,wilcoxon_res)
    all_TMM_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(fa_res_name),get(ANCOM_res_name),get(ANCOMBC_res_name),get(maaslin2_cplm_res_name),get(maaslin2_negbin_res_name),get(lmem_res_name),get(lm_inter_res_name),get(lfem_res_name),
                         get(maaslin2_lm_res_name),get(maaslin2_lm_LOG_res_name),get(limma_res_name),get(wilcoxon_res_name))
  }
  ################################TSS
  print('TSS start!')
  normalization_method = 'TSS'
  fa_res_name <- paste0("fastANCOM_res_",normalization_method)
  assign(fa_res_name,fastANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(fa_res_name),file = paste0('fastANCOM_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('fastANCOM FINISHED')
  maaslin2_lm_res_name <- paste0("maaslin2_lm_res_",data_index,"_",normalization_method)
  assign(maaslin2_lm_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_lm_res_name),file = paste0('maaslin2_lm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm FINISHED')
  maaslin2_lm_LOG_res_name <- paste0("maaslin2_lm_LOG_res_",normalization_method)
  assign(maaslin2_lm_LOG_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'LOG'))
  write.table(get(maaslin2_lm_LOG_res_name),file = paste0('maaslin2_lm_LOG_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm LOG FINISHED')
  maaslin2_cplm_res_name <- paste0("maaslin2_cplm_res_",normalization_method)
  assign(maaslin2_cplm_res_name,Maaslin2_CPLM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_cplm_res_name),file = paste0('maaslin2_cplm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-cplm FINISHED')
  maaslin2_cplm_LOG_res_name <- paste0("maaslin2_cplm_LOG_res_",normalization_method)
  assign(maaslin2_cplm_LOG_res_name,Maaslin2_CPLM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_cplm_LOG_res_name),file = paste0('maaslin2_cplm_LOG_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-cplm LOG FINISHED')
  limma_res_name <- paste0("limma_res_",normalization_method)
  assign(limma_res_name,limma_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(limma_res_name),file = paste0('limma_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('limma FINISHED')
  wilcoxon_res_name <- paste0("wilcoxon_res_",normalization_method)
  assign(wilcoxon_res_name,Wilcoxon_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(wilcoxon_res_name),file = paste0('Wilcoxon_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Wilcoxon test FINISHED')
  ZicoSeq_res_name <- paste0("ZicoSeq_res_",normalization_method)
  assign(ZicoSeq_res_name,ZicoSeq_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(ZicoSeq_res_name),file = paste0('ZicoSeq_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ZicoSeq FINISHED')
  edgeR_res_name <- paste0("edgeR_res_",normalization_method)
  assign(edgeR_res_name,edgeR_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(edgeR_res_name),file = paste0('edgeR_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('edgeR FINISHED')
  if(is.null(covariates)){
    lm_res_name <- paste0("lm_res_",normalization_method)
    assign(lm_res_name,LM_analysis(simulated_data_all,normalization_method))
    write.table(get(lm_res_name),file = paste0('LM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM FINISHED')
    all_TSS_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(fa_res_name),get(maaslin2_cplm_res_name),get(maaslin2_cplm_LOG_res_name),
                         get(maaslin2_lm_res_name),get(maaslin2_lm_LOG_res_name),get(limma_res_name),get(lm_res_name),get(wilcoxon_res_name))
  }else{
    
    lmem_res_name <- paste0("lmem_res_",normalization_method)
    assign(lmem_res_name,LM_random_effect(simulated_data_all,covariates,categorical_variable_name,normalization_method))
    write.table(get(lmem_res_name),file = paste0('lmem_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LMEM FINISHED')
    lm_inter_res_name <- paste0("lm_inter_res_",normalization_method)
    assign(lm_inter_res_name,LM_interaction_analysis(simulated_data_all,covariates,normalization_method))
    write.table(get(lm_inter_res_name),file = paste0('LM-inter_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-interaction FINISHED')
    lfem_res_name <- paste0("lfem_res_",normalization_method)
    assign(lfem_res_name,LM_fixed_effect(simulated_data_all,covariates,normalization_method))
    write.table(get(lfem_res_name),file = paste0('LM-fixed_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-fixed effect FINISHED')
    #all_result<-list(fa_res,ancom_res,ancombc_res,limma_res,lm_re_res,
    #                 lm_inter_res,lm_fe_res,vtwins_res,wilcoxon_res)
    all_TSS_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(fa_res_name),get(maaslin2_cplm_res_name),get(maaslin2_cplm_LOG_res_name),get(lmem_res_name),get(lm_inter_res_name),get(lfem_res_name),
                         get(maaslin2_lm_res_name),get(maaslin2_lm_LOG_res_name),get(limma_res_name),get(wilcoxon_res_name))
  }
  #########################################LOG
  print('LOG start!')
  normalization_method = 'LOG'
  maaslin2_lm_res_name <- paste0("maaslin2_lm_res_",data_index,"_",normalization_method)
  assign(maaslin2_lm_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization='None',transform = normalization_method))
  write.table(get(maaslin2_lm_res_name),file = paste0('maaslin2_lm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm FINISHED')
  
  limma_res_name <- paste0("limma_res_",normalization_method)
  assign(limma_res_name,limma_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(limma_res_name),file = paste0('limma_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('limma FINISHED')
  wilcoxon_res_name <- paste0("wilcoxon_res_",normalization_method)
  assign(wilcoxon_res_name,Wilcoxon_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(wilcoxon_res_name),file = paste0('Wilcoxon_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Wilcoxon test FINISHED')
  ZicoSeq_res_name <- paste0("ZicoSeq_res_",normalization_method)
  assign(ZicoSeq_res_name,ZicoSeq_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(ZicoSeq_res_name),file = paste0('ZicoSeq_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ZicoSeq FINISHED')
  if(is.null(covariates)){
    lm_res_name <- paste0("lm_res_",normalization_method)
    assign(lm_res_name,LM_analysis(simulated_data_all,normalization_method))
    write.table(get(lm_res_name),file = paste0('LM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM FINISHED')
    all_LOG_result<-list(get(ZicoSeq_res_name),get(maaslin2_lm_res_name),get(limma_res_name),get(lm_res_name),get(wilcoxon_res_name))
  }else{
    
    lmem_res_name <- paste0("lmem_res_",normalization_method)
    assign(lmem_res_name,LM_random_effect(simulated_data_all,covariates,categorical_variable_name,normalization_method))
    write.table(get(lmem_res_name),file = paste0('lmem_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LMEM FINISHED')
    lm_inter_res_name <- paste0("lm_inter_res_",normalization_method)
    assign(lm_inter_res_name,LM_interaction_analysis(simulated_data_all,covariates,normalization_method))
    write.table(get(lm_inter_res_name),file = paste0('LM-inter_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-interaction FINISHED')
    lfem_res_name <- paste0("lfem_res_",normalization_method)
    assign(lfem_res_name,LM_fixed_effect(simulated_data_all,covariates,normalization_method))
    write.table(get(lfem_res_name),file = paste0('LM-fixed_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-fixed effect FINISHED')
    #all_result<-list(fa_res,ancom_res,ancombc_res,limma_res,lm_re_res,
    #                 lm_inter_res,lm_fe_res,vtwins_res,wilcoxon_res)
    all_LOG_result<-list(get(ZicoSeq_res_name),get(maaslin2_lm_res_name),get(limma_res_name),get(wilcoxon_res_name),
                         get(lmem_res_name),get(lm_inter_res_name),get(lfem_res_name))
  }
  all_result<-list(all_CLR_result=all_CLR_result,all_CSS_result=all_CSS_result,all_TSS_result=all_TSS_result,
  all_TMM_result=all_TMM_result,all_count_result=all_count_result,all_LOG_result=all_LOG_result)
  
  return(all_result)
}
merge_all_fast_method_result2<-function(simulated_data_all,covariates=NULL,categorical_variable_name=NULL,prevalence = 0.1,data_index){
  simulated_data_all$simulated_data_rela <- filter_by_prevalence(simulated_data_all$simulated_data_rela,prevalence)
  simulated_data_all$simulated_data_abs <- filter_by_prevalence(simulated_data_all$simulated_data_abs,prevalence)
  print(colSums(simulated_data_all$simulated_data_rela != 0))
  #####################################CLR
  print('CLR start!')
  normalization_method = 'CLR'
  # ANCOM_res_name <- paste0("ANCOM_res_",normalization_method)
  # assign(ANCOM_res_name,ANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  # write.table(get(ANCOM_res_name),file = paste0('ANCOM_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  # print('ANCOM FINISHED')
  maaslin2_lm_res_name <- paste0("maaslin2_lm_res_",data_index,"_",normalization_method)
  assign(maaslin2_lm_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_lm_res_name),file = paste0('maaslin2_lm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm FINISHED')
  #maaslin2_lm_LOG_res_name <- paste0("maaslin2_lm_LOG_res_",normalization_method)
  #assign(maaslin2_lm_LOG_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'LOG'))
  #write.table(get(maaslin2_lm_LOG_res_name),file = paste0('maaslin2_lm_LOG_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  #print('maaslin2-lm LOG FINISHED')
  limma_res_name <- paste0("limma_res_",normalization_method)
  assign(limma_res_name,limma_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(limma_res_name),file = paste0('limma_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('limma FINISHED')
  wilcoxon_res_name <- paste0("wilcoxon_res_",normalization_method)
  assign(wilcoxon_res_name,Wilcoxon_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(wilcoxon_res_name),file = paste0('Wilcoxon_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Wilcoxon test FINISHED')
  #ZicoSeq_res_name <- paste0("ZicoSeq_res_",normalization_method)
  #assign(ZicoSeq_res_name,ZicoSeq_analysis(simulated_data_all,covariates,normalization_method))
  #write.table(get(ZicoSeq_res_name),file = paste0('ZicoSeq_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  #print('ZicoSeq FINISHED')
  if(is.null(covariates)){
    lm_res_name <- paste0("lm_res_",normalization_method)
    assign(lm_res_name,LM_analysis(simulated_data_all,normalization_method))
    write.table(get(lm_res_name),file = paste0('LM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM FINISHED')
    all_CLR_result<-list(get(maaslin2_lm_res_name),get(limma_res_name),get(lm_res_name),get(wilcoxon_res_name))
  }else{
    
    lmem_res_name <- paste0("lmem_res_",normalization_method)
    assign(lmem_res_name,LM_random_effect(simulated_data_all,covariates,categorical_variable_name,normalization_method))
    write.table(get(lmem_res_name),file = paste0('lmem_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LMEM FINISHED')
    lm_inter_res_name <- paste0("lm_inter_res_",normalization_method)
    assign(lm_inter_res_name,LM_interaction_analysis(simulated_data_all,covariates,normalization_method))
    write.table(get(lm_inter_res_name),file = paste0('LM-inter_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-interaction FINISHED')
    lfem_res_name <- paste0("lfem_res_",normalization_method)
    assign(lfem_res_name,LM_fixed_effect(simulated_data_all,covariates,normalization_method))
    write.table(get(lfem_res_name),file = paste0('LM-fixed_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-fixed effect FINISHED')
    #all_result<-list(fa_res,ancom_res,ancombc_res,limma_res,lm_re_res,
    #                 lm_inter_res,lm_fe_res,vtwins_res,wilcoxon_res)
    all_CLR_result<-list(get(maaslin2_lm_res_name),get(limma_res_name),get(wilcoxon_res_name),
      get(lmem_res_name),get(lm_inter_res_name),get(lfem_res_name))
  }
  ############################################CSS
  print('CSS start!')
  normalization_method = 'CSS'
  fa_res_name <- paste0("fastANCOM_res_",normalization_method)
  assign(fa_res_name,fastANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(fa_res_name),file = paste0('fastANCOM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  # ANCOM_res_name <- paste0("ANCOM_res_",normalization_method)
  # assign(ANCOM_res_name,ANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  # write.table(get(ANCOM_res_name),file = paste0('ANCOM_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  # print('ANCOM FINISHED')
  ANCOMBC_res_name <- paste0("ANCOMBC_res_",normalization_method)
  assign(ANCOMBC_res_name,ANCOMBC_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(ANCOMBC_res_name),file = paste0('ANCOMBC_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ANCOMBC FINISHED')
  maaslin2_lm_res_name <- paste0("maaslin2_lm_res_",data_index,"_",normalization_method)
  assign(maaslin2_lm_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_lm_res_name),file = paste0('maaslin2_lm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm FINISHED')
  maaslin2_lm_LOG_res_name <- paste0("maaslin2_lm_LOG_res_",normalization_method)
  assign(maaslin2_lm_LOG_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'LOG'))
  write.table(get(maaslin2_lm_LOG_res_name),file = paste0('maaslin2_lm_LOG_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm LOG FINISHED')
  maaslin2_cplm_res_name <- paste0("maaslin2_cplm_res_",normalization_method)
  assign(maaslin2_cplm_res_name,Maaslin2_CPLM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_cplm_res_name),file = paste0('maaslin2_cplm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-cplm FINISHED')
  maaslin2_cplm_LOG_res_name <- paste0("maaslin2_cplm_LOG_res_",normalization_method)
  assign(maaslin2_cplm_LOG_res_name,Maaslin2_CPLM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'LOG'))
  write.table(get(maaslin2_cplm_LOG_res_name),file = paste0('maaslin2_cplm_LOG_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-cplm LOG FINISHED')
  maaslin2_negbin_res_name <- paste0("maaslin2_negbin_res_",normalization_method)
  assign(maaslin2_negbin_res_name,Maaslin2_NEGBIN_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_negbin_res_name),file = paste0('maaslin2_negbin_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-negbin FINISHED')
  limma_res_name <- paste0("limma_res_",normalization_method)
  assign(limma_res_name,limma_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(limma_res_name),file = paste0('limma_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('limma FINISHED')
  wilcoxon_res_name <- paste0("wilcoxon_res_",normalization_method)
  assign(wilcoxon_res_name,Wilcoxon_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(wilcoxon_res_name),file = paste0('Wilcoxon_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Wilcoxon test FINISHED')
  ZicoSeq_res_name <- paste0("ZicoSeq_res_",normalization_method)
  assign(ZicoSeq_res_name,ZicoSeq_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(ZicoSeq_res_name),file = paste0('ZicoSeq_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ZicoSeq FINISHED')
  edgeR_res_name <- paste0("edgeR_res_",normalization_method)
  assign(edgeR_res_name,edgeR_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(edgeR_res_name),file = paste0('edgeR_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('edgeR FINISHED')
  if(is.null(covariates)){
    lm_res_name <- paste0("lm_res_",normalization_method)
    assign(lm_res_name,LM_analysis(simulated_data_all,normalization_method))
    write.table(get(lm_res_name),file = paste0('LM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM FINISHED')
    all_CSS_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(fa_res_name),get(ANCOMBC_res_name),get(maaslin2_cplm_res_name),get(maaslin2_cplm_LOG_res_name),get(maaslin2_negbin_res_name),
                         get(maaslin2_lm_res_name),get(maaslin2_lm_LOG_res_name),get(limma_res_name),get(lm_res_name),get(wilcoxon_res_name))
  }else{
    lmem_res_name <- paste0("lmem_res_",normalization_method)
    assign(lmem_res_name,LM_random_effect(simulated_data_all,covariates,categorical_variable_name,normalization_method))
    write.table(get(lmem_res_name),file = paste0('lmem_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LMEM FINISHED')
    lm_inter_res_name <- paste0("lm_inter_res_",normalization_method)
    assign(lm_inter_res_name,LM_interaction_analysis(simulated_data_all,covariates,normalization_method))
    write.table(get(lm_inter_res_name),file = paste0('LM-inter_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-interaction FINISHED')
    lfem_res_name <- paste0("lfem_res_",normalization_method)
    assign(lfem_res_name,LM_fixed_effect(simulated_data_all,covariates,normalization_method))
    write.table(get(lfem_res_name),file = paste0('LM-fixed_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-fixed effect FINISHED')
    #all_result<-list(fa_res,ancom_res,ancombc_res,limma_res,lm_re_res,
    #                 lm_inter_res,lm_fe_res,vtwins_res,wilcoxon_res)
    all_CSS_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(fa_res_name),get(ANCOMBC_res_name),get(maaslin2_cplm_res_name),get(maaslin2_cplm_LOG_res_name),get(maaslin2_negbin_res_name),get(lmem_res_name),get(lm_inter_res_name),get(lfem_res_name),
         get(maaslin2_lm_res_name),get(maaslin2_lm_LOG_res_name),get(limma_res_name),get(wilcoxon_res_name))
  }
  ####################################count
  print('Count start!')
  normalization_method = 'count'
  fa_res_name <- paste0("fastANCOM_res_",normalization_method)
  assign(fa_res_name,fastANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(fa_res_name),file = paste0('fastANCOM_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('fastANCOM FINISHED')
  # ANCOM_res_name <- paste0("ANCOM_res_",normalization_method)
  # assign(ANCOM_res_name,ANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  # write.table(get(ANCOM_res_name),file = paste0('ANCOM_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  # print('ANCOM FINISHED')
  ANCOMBC_res_name <- paste0("ANCOMBC_res_",normalization_method)
  assign(ANCOMBC_res_name,ANCOMBC_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(ANCOMBC_res_name),file = paste0('ANCOMBC_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ANCOMBC FINISHED')
  maaslin2_lm_res_name <- paste0("maaslin2_lm_res_",data_index,"_",normalization_method)
  assign(maaslin2_lm_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization='None',transform = 'None'))
  write.table(get(maaslin2_lm_res_name),file = paste0('maaslin2_lm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm FINISHED')
  maaslin2_cplm_res_name <- paste0("maaslin2_cplm_res_",normalization_method)
  assign(maaslin2_cplm_res_name,Maaslin2_CPLM_analysis(simulated_data_all,covariates,normalization='None',transform = 'None'))
  write.table(get(maaslin2_cplm_res_name),file = paste0('maaslin2_cplm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-cplm FINISHED')
  maaslin2_zinb_res_name <- paste0("maaslin2_zinb_res_",normalization_method)
  assign(maaslin2_zinb_res_name,Maaslin2_ZINB_analysis(simulated_data_all,covariates,normalization='None',transform = 'None'))
  write.table(get(maaslin2_zinb_res_name),file = paste0('maaslin2_zinb_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-zinb FINISHED')
  maaslin2_negbin_res_name <- paste0("maaslin2_negbin_res_",normalization_method)
  assign(maaslin2_negbin_res_name,Maaslin2_NEGBIN_analysis(simulated_data_all,covariates,normalization='None',transform = 'None'))
  write.table(get(maaslin2_negbin_res_name),file = paste0('maaslin2_megbin_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-negbin FINISHED')
  limma_res_name <- paste0("limma_res_",normalization_method)
  assign(limma_res_name,limma_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(limma_res_name),file = paste0('limma_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('limma FINISHED')
  wilcoxon_res_name <- paste0("wilcoxon_res_",normalization_method)
  assign(wilcoxon_res_name,Wilcoxon_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(wilcoxon_res_name),file = paste0('Wilcoxon_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Wilcoxon test FINISHED')
  ZicoSeq_res_name <- paste0("ZicoSeq_res_",normalization_method)
  assign(ZicoSeq_res_name,ZicoSeq_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(ZicoSeq_res_name),file = paste0('ZicoSeq_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ZicoSeq FINISHED')
  edgeR_res_name <- paste0("edgeR_res_",normalization_method)
  assign(edgeR_res_name,edgeR_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(edgeR_res_name),file = paste0('edgeR_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('edgeR FINISHED')
  DESeq2_res_name <- paste0("DESeq2_res_",normalization_method)
  assign(DESeq2_res_name,DESeq2_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(DESeq2_res_name),file = paste0('DESeq2_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('DESeq2 FINISHED')
  ALDEx2_res_name <- paste0("ALDEx2_res_",normalization_method)
  assign(ALDEx2_res_name,ALDEx2_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(ALDEx2_res_name),file = paste0('ALDEx2_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ALDEx2 FINISHED')
  Corncob_res_name <- paste0("Corncob_res_",normalization_method)
  assign(Corncob_res_name,Corncob_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(Corncob_res_name),file = paste0('Corncob_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Corncob FINISHED')
  if(is.null(covariates)){
    lm_res_name <- paste0("lm_res_",normalization_method)
    assign(lm_res_name,LM_analysis(simulated_data_all,normalization_method))
    write.table(get(lm_res_name),file = paste0('LM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM FINISHED')
    all_count_result<-list(get(DESeq2_res_name),get(ALDEx2_res_name),get(Corncob_res_name),get(edgeR_res_name),get(ZicoSeq_res_name),get(ANCOMBC_res_name),get(fa_res_name),get(maaslin2_cplm_res_name),get(maaslin2_zinb_res_name),get(maaslin2_negbin_res_name),
                         get(maaslin2_lm_res_name),get(limma_res_name),get(lm_res_name),get(wilcoxon_res_name))
  }else{
    
    lmem_res_name <- paste0("lmem_res_",normalization_method)
    assign(lmem_res_name,LM_random_effect(simulated_data_all,covariates,categorical_variable_name,normalization_method))
    write.table(get(lmem_res_name),file = paste0('lmem_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LMEM FINISHED')
    lm_inter_res_name <- paste0("lm_inter_res_",normalization_method)
    assign(lm_inter_res_name,LM_interaction_analysis(simulated_data_all,covariates,normalization_method))
    write.table(get(lm_inter_res_name),file = paste0('LM-inter_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-interaction FINISHED')
    lfem_res_name <- paste0("lfem_res_",normalization_method)
    assign(lfem_res_name,LM_fixed_effect(simulated_data_all,covariates,normalization_method))
    write.table(get(lfem_res_name),file = paste0('LM-fixed_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-fixed effect FINISHED')
    #all_result<-list(fa_res,ancom_res,ancombc_res,limma_res,lm_re_res,
    #                 lm_inter_res,lm_fe_res,vtwins_res,wilcoxon_res)
    all_count_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(lmem_res_name),get(lm_inter_res_name),get(lfem_res_name),get(ANCOMBC_res_name),get(fa_res_name),get(maaslin2_cplm_res_name),get(maaslin2_zinb_res_name),get(maaslin2_negbin_res_name),
                         get(maaslin2_lm_res_name),get(limma_res_name),get(wilcoxon_res_name))
  }
  #####################################TMM
  print('TMM start!')
  normalization_method = 'TMM'
  fa_res_name <- paste0("fastANCOM_res_",normalization_method)
  assign(fa_res_name,fastANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(fa_res_name),file = paste0('fastANCOM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  # ANCOM_res_name <- paste0("ANCOM_res_",normalization_method)
  # assign(ANCOM_res_name,ANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  # write.table(get(ANCOM_res_name),file = paste0('ANCOM_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  # print('ANCOM FINISHED')
  ANCOMBC_res_name <- paste0("ANCOMBC_res_",normalization_method)
  assign(ANCOMBC_res_name,ANCOMBC_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(ANCOMBC_res_name),file = paste0('ANCOMBC_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ANCOMBC FINISHED')
  maaslin2_lm_res_name <- paste0("maaslin2_lm_res_",data_index,"_",normalization_method)
  assign(maaslin2_lm_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_lm_res_name),file = paste0('maaslin2_lm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm FINISHED')
  maaslin2_lm_LOG_res_name <- paste0("maaslin2_lm_LOG_res_",normalization_method)
  assign(maaslin2_lm_LOG_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'LOG'))
  write.table(get(maaslin2_lm_LOG_res_name),file = paste0('maaslin2_lm_LOG_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm LOG FINISHED')
  maaslin2_cplm_res_name <- paste0("maaslin2_cplm_res_",normalization_method)
  assign(maaslin2_cplm_res_name,Maaslin2_CPLM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_cplm_res_name),file = paste0('maaslin2_cplm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-cplm FINISHED')
  maaslin2_negbin_res_name <- paste0("maaslin2_negbin_res_",normalization_method)
  assign(maaslin2_negbin_res_name,Maaslin2_NEGBIN_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_negbin_res_name),file = paste0('maaslin2_negbin_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-negbin FINISHED')
  limma_res_name <- paste0("limma_res_",normalization_method)
  assign(limma_res_name,limma_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(limma_res_name),file = paste0('limma_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('limma FINISHED')
  wilcoxon_res_name <- paste0("wilcoxon_res_",normalization_method)
  assign(wilcoxon_res_name,Wilcoxon_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(wilcoxon_res_name),file = paste0('Wilcoxon_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Wilcoxon test FINISHED')
  ZicoSeq_res_name <- paste0("ZicoSeq_res_",normalization_method)
  assign(ZicoSeq_res_name,ZicoSeq_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(ZicoSeq_res_name),file = paste0('ZicoSeq_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ZicoSeq FINISHED')
  edgeR_res_name <- paste0("edgeR_res_",normalization_method)
  assign(edgeR_res_name,edgeR_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(edgeR_res_name),file = paste0('edgeR_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('edgeR FINISHED')
  if(is.null(covariates)){
    lm_res_name <- paste0("lm_res_",normalization_method)
    assign(lm_res_name,LM_analysis(simulated_data_all,normalization_method))
    write.table(get(lm_res_name),file = paste0('LM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM FINISHED')
    all_TMM_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(fa_res_name),get(ANCOMBC_res_name),get(maaslin2_cplm_res_name),get(maaslin2_negbin_res_name),
                         get(maaslin2_lm_res_name),get(maaslin2_lm_LOG_res_name),get(limma_res_name),get(lm_res_name),get(wilcoxon_res_name))
  }else{
    
    lmem_res_name <- paste0("lmem_res_",normalization_method)
    assign(lmem_res_name,LM_random_effect(simulated_data_all,covariates,categorical_variable_name,normalization_method))
    write.table(get(lmem_res_name),file = paste0('lmem_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LMEM FINISHED')
    lm_inter_res_name <- paste0("lm_inter_res_",normalization_method)
    assign(lm_inter_res_name,LM_interaction_analysis(simulated_data_all,covariates,normalization_method))
    write.table(get(lm_inter_res_name),file = paste0('LM-inter_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-interaction FINISHED')
    lfem_res_name <- paste0("lfem_res_",normalization_method)
    assign(lfem_res_name,LM_fixed_effect(simulated_data_all,covariates,normalization_method))
    write.table(get(lfem_res_name),file = paste0('LM-fixed_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-fixed effect FINISHED')
    #all_result<-list(fa_res,ancom_res,ancombc_res,limma_res,lm_re_res,
    #                 lm_inter_res,lm_fe_res,vtwins_res,wilcoxon_res)
    all_TMM_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(fa_res_name),get(ANCOMBC_res_name),get(maaslin2_cplm_res_name),get(maaslin2_negbin_res_name),get(lmem_res_name),get(lm_inter_res_name),get(lfem_res_name),
                         get(maaslin2_lm_res_name),get(maaslin2_lm_LOG_res_name),get(limma_res_name),get(wilcoxon_res_name))
  }
  ################################TSS
  print('TSS start!')
  normalization_method = 'TSS'
  fa_res_name <- paste0("fastANCOM_res_",normalization_method)
  assign(fa_res_name,fastANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(fa_res_name),file = paste0('fastANCOM_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('fastANCOM FINISHED')
  maaslin2_lm_res_name <- paste0("maaslin2_lm_res_",data_index,"_",normalization_method)
  assign(maaslin2_lm_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_lm_res_name),file = paste0('maaslin2_lm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm FINISHED')
  maaslin2_lm_LOG_res_name <- paste0("maaslin2_lm_LOG_res_",normalization_method)
  assign(maaslin2_lm_LOG_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'LOG'))
  write.table(get(maaslin2_lm_LOG_res_name),file = paste0('maaslin2_lm_LOG_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm LOG FINISHED')
  maaslin2_cplm_res_name <- paste0("maaslin2_cplm_res_",normalization_method)
  assign(maaslin2_cplm_res_name,Maaslin2_CPLM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_cplm_res_name),file = paste0('maaslin2_cplm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-cplm FINISHED')
  maaslin2_cplm_LOG_res_name <- paste0("maaslin2_cplm_LOG_res_",normalization_method)
  assign(maaslin2_cplm_LOG_res_name,Maaslin2_CPLM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'None'))
  write.table(get(maaslin2_cplm_LOG_res_name),file = paste0('maaslin2_cplm_LOG_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-cplm LOG FINISHED')
  limma_res_name <- paste0("limma_res_",normalization_method)
  assign(limma_res_name,limma_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(limma_res_name),file = paste0('limma_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('limma FINISHED')
  wilcoxon_res_name <- paste0("wilcoxon_res_",normalization_method)
  assign(wilcoxon_res_name,Wilcoxon_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(wilcoxon_res_name),file = paste0('Wilcoxon_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Wilcoxon test FINISHED')
  ZicoSeq_res_name <- paste0("ZicoSeq_res_",normalization_method)
  assign(ZicoSeq_res_name,ZicoSeq_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(ZicoSeq_res_name),file = paste0('ZicoSeq_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ZicoSeq FINISHED')
  edgeR_res_name <- paste0("edgeR_res_",normalization_method)
  assign(edgeR_res_name,edgeR_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(edgeR_res_name),file = paste0('edgeR_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('edgeR FINISHED')
  if(is.null(covariates)){
    lm_res_name <- paste0("lm_res_",normalization_method)
    assign(lm_res_name,LM_analysis(simulated_data_all,normalization_method))
    write.table(get(lm_res_name),file = paste0('LM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM FINISHED')
    all_TSS_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(fa_res_name),get(maaslin2_cplm_res_name),get(maaslin2_cplm_LOG_res_name),
                         get(maaslin2_lm_res_name),get(maaslin2_lm_LOG_res_name),get(limma_res_name),get(lm_res_name),get(wilcoxon_res_name))
  }else{
    
    lmem_res_name <- paste0("lmem_res_",normalization_method)
    assign(lmem_res_name,LM_random_effect(simulated_data_all,covariates,categorical_variable_name,normalization_method))
    write.table(get(lmem_res_name),file = paste0('lmem_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LMEM FINISHED')
    lm_inter_res_name <- paste0("lm_inter_res_",normalization_method)
    assign(lm_inter_res_name,LM_interaction_analysis(simulated_data_all,covariates,normalization_method))
    write.table(get(lm_inter_res_name),file = paste0('LM-inter_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-interaction FINISHED')
    lfem_res_name <- paste0("lfem_res_",normalization_method)
    assign(lfem_res_name,LM_fixed_effect(simulated_data_all,covariates,normalization_method))
    write.table(get(lfem_res_name),file = paste0('LM-fixed_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-fixed effect FINISHED')
    #all_result<-list(fa_res,ancom_res,ancombc_res,limma_res,lm_re_res,
    #                 lm_inter_res,lm_fe_res,vtwins_res,wilcoxon_res)
    all_TSS_result<-list(get(edgeR_res_name),get(ZicoSeq_res_name),get(fa_res_name),get(maaslin2_cplm_res_name),get(maaslin2_cplm_LOG_res_name),get(lmem_res_name),get(lm_inter_res_name),get(lfem_res_name),
                         get(maaslin2_lm_res_name),get(maaslin2_lm_LOG_res_name),get(limma_res_name),get(wilcoxon_res_name))
  }
  #########################################LOG
  print('LOG start!')
  normalization_method = 'LOG'
  maaslin2_lm_res_name <- paste0("maaslin2_lm_res_",data_index,"_",normalization_method)
  assign(maaslin2_lm_res_name,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization='None',transform = normalization_method))
  write.table(get(maaslin2_lm_res_name),file = paste0('maaslin2_lm_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm FINISHED')
  
  limma_res_name <- paste0("limma_res_",normalization_method)
  assign(limma_res_name,limma_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(limma_res_name),file = paste0('limma_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('limma FINISHED')
  wilcoxon_res_name <- paste0("wilcoxon_res_",normalization_method)
  assign(wilcoxon_res_name,Wilcoxon_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(wilcoxon_res_name),file = paste0('Wilcoxon_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('Wilcoxon test FINISHED')
  ZicoSeq_res_name <- paste0("ZicoSeq_res_",normalization_method)
  assign(ZicoSeq_res_name,ZicoSeq_analysis(simulated_data_all,covariates,normalization_method))
  write.table(get(ZicoSeq_res_name),file = paste0('ZicoSeq_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ZicoSeq FINISHED')
  if(is.null(covariates)){
    lm_res_name <- paste0("lm_res_",normalization_method)
    assign(lm_res_name,LM_analysis(simulated_data_all,normalization_method))
    write.table(get(lm_res_name),file = paste0('LM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM FINISHED')
    all_LOG_result<-list(get(ZicoSeq_res_name),get(maaslin2_lm_res_name),get(limma_res_name),get(lm_res_name),get(wilcoxon_res_name))
  }else{
    
    lmem_res_name <- paste0("lmem_res_",normalization_method)
    assign(lmem_res_name,LM_random_effect(simulated_data_all,covariates,categorical_variable_name,normalization_method))
    write.table(get(lmem_res_name),file = paste0('lmem_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LMEM FINISHED')
    lm_inter_res_name <- paste0("lm_inter_res_",normalization_method)
    assign(lm_inter_res_name,LM_interaction_analysis(simulated_data_all,covariates,normalization_method))
    write.table(get(lm_inter_res_name),file = paste0('LM-inter_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-interaction FINISHED')
    lfem_res_name <- paste0("lfem_res_",normalization_method)
    assign(lfem_res_name,LM_fixed_effect(simulated_data_all,covariates,normalization_method))
    write.table(get(lfem_res_name),file = paste0('LM-fixed_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
    print('LM-fixed effect FINISHED')
    #all_result<-list(fa_res,ancom_res,ancombc_res,limma_res,lm_re_res,
    #                 lm_inter_res,lm_fe_res,vtwins_res,wilcoxon_res)
    all_LOG_result<-list(get(ZicoSeq_res_name),get(maaslin2_lm_res_name),get(limma_res_name),get(wilcoxon_res_name),
                         get(lmem_res_name),get(lm_inter_res_name),get(lfem_res_name))
  }
  all_result<-list(all_CLR_result=all_CLR_result,all_CSS_result=all_CSS_result,all_TSS_result=all_TSS_result,
  all_TMM_result=all_TMM_result,all_count_result=all_count_result,all_LOG_result=all_LOG_result)
  
  return(all_result)
}
merge_all_maaslin2_result <- function(simulated_data_all,covariates=NULL,categorical_variable_name=NULL,data_index,normalization,transform){
  maaslin2_lm_res<-Maaslin2_LM_analysis(simulated_data_all,covariates,normalization,transform)
  print('MaAslin2-LM FINISHED')
  maaslin2_cplm_res<-Maaslin2_CPLM_analysis(simulated_data_all,covariates,normalization,transform)
  print('MaAslin2-CPLM FINISHED')
  maaslin2_zinb_res<-Maaslin2_ZINB_analysis(simulated_data_all,covariates,normalization,transform)
  print('MaAslin2-ZINB FINISHED')
  maaslin2_negbin_res<-Maaslin2_NEGBIN_analysis(simulated_data_all,covariates,normalization,transform)
  print('MaAslin2-NEGBIN FINISHED')
  all_result <- list(maaslin2_lm_res,
                     maaslin2_cplm_res,maaslin2_zinb_res,
                     maaslin2_negbin_res)
  write.table(maaslin2_lm_res,file = paste0('Maaslin2-LM_res',data_index,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  write.table(maaslin2_cplm_res,file = paste0('Maaslin2-CPLM_res',data_index,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  write.table(maaslin2_zinb_res,file = paste0('Maaslin2-ZINB_res',data_index,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  write.table(maaslin2_negbin_res,file = paste0('Maaslin2-NEGBIN_res',data_index,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  return(all_result)
}
merge_all_ancom_result <- function(simulated_data_all,covariates=NULL,categorical_variable_name=NULL,data_index,normalization){
  fa_res<-fastANCOM_analysis(simulated_data_all,covariates,normalization)
  print('fastANCOM FINISHED')
  ancom_res<-ANCOM_analysis(simulated_data_all,covariates,normalization)
  print('ANCOM FINISHED')
  ancombc_res<-ANCOMBC_analysis(simulated_data_all,covariates,normalization)
  print('ANCOM-BC FINISHED')
  all_result <- list(fa_res,ancom_res,ancombc_res)
  write.table(fa_res,file = paste0('fastANCOM_res_',data_index,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  write.table(ancom_res,file = paste0('ANCOM_res_',data_index,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  write.table(ancombc_res,file = paste0('ANCOM-BC_res',data_index,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  return(all_result)
}

merge_better_method_result<-function(simulated_data_all,covariates=NULL,categorical_variable_name=NULL,prevalence = 0.1,data_index){
  simulated_data_all$simulated_data_rela <- filter_by_prevalence(simulated_data_all$simulated_data_rela,prevalence)
  simulated_data_all$simulated_data_abs <- filter_by_prevalence(simulated_data_all$simulated_data_abs,prevalence)
  print(colSums(simulated_data_all$simulated_data_rela != 0))
  ############################################CSS
  print('CSS start!')
  normalization_method = 'CSS'
  fa_res_name1 <- paste0("fastANCOM_res_",normalization_method)
  assign(fa_res_name1,fastANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(fa_res_name1),file = paste0('fastANCOM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  
  ANCOMBC_res_name1 <- paste0("ANCOMBC_res_",normalization_method)
  assign(ANCOMBC_res_name1,ANCOMBC_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(ANCOMBC_res_name1),file = paste0('ANCOMBC_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ANCOMBC FINISHED')

  edgeR_res_name <- paste0("edgeR_res_",normalization_method)
  assign(edgeR_res_name,edgeR_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(edgeR_res_name),file = paste0('edgeR_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('edgeR FINISHED')

  lmem_res_name <- paste0("lmem_res_",normalization_method)
  assign(lmem_res_name,LM_random_effect(simulated_data_all,covariates,categorical_variable_name,normalization=normalization_method))
  write.table(get(lmem_res_name),file = paste0('lmem_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('LM random effect FINISHED')

  limma_res_name <- paste0("limma_res_",normalization_method)
  assign(limma_res_name,limma_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(limma_res_name),file = paste0('limma_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('limma FINISHED')

  lfem_res_name <- paste0("lfem_res_",normalization_method)
  assign(lfem_res_name,LM_fixed_effect(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(lfem_res_name),file = paste0('LM-fixed_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('LM-fixed FINISHED')

  maaslin2_lm_LOG_res_name1 <- paste0("maaslin2_lm_LOG_res_",normalization_method)
  assign(maaslin2_lm_LOG_res_name1,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'LOG'))
  write.table(get(maaslin2_lm_LOG_res_name1),file = paste0('maaslin2_lm_LOG_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm LOG FINISHED')
  ####################################count
  print('Count start!')
  normalization_method = 'count'
  # fa_res_name2 <- paste0("fastANCOM_res_",normalization_method)
  # assign(fa_res_name2,fastANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  # write.table(get(fa_res_name2),file = paste0('fastANCOM_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  # print('fastANCOM FINISHED')
  # ZicoSeq_res_name1 <- paste0("ZicoSeq_res_",normalization_method)
  # assign(ZicoSeq_res_name1,ZicoSeq_analysis(simulated_data_all,covariates,normalization_method))
  # write.table(get(ZicoSeq_res_name1),file = paste0('ZicoSeq_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  # print('ZicoSeq FINISHED')
  # ALDEx2_res_name <- paste0("ALDEx2_res_",normalization_method)
  # assign(ALDEx2_res_name,ALDEx2_analysis(simulated_data_all,covariates,normalization_method))
  # write.table(get(ALDEx2_res_name),file = paste0('ALDEx2_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  # print('ALDEx2 FINISHED')
  ANCOMBC_res_name2 <- paste0("ANCOMBC_res_",normalization_method)
  assign(ANCOMBC_res_name2,ANCOMBC_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(ANCOMBC_res_name2),file = paste0('ANCOMBC_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ANCOMBC FINISHED')
  #####################################TMM
  print('TMM start!')
  normalization_method = 'TMM'
  fa_res_name2 <- paste0("fastANCOM_res_",normalization_method)
  assign(fa_res_name2,fastANCOM_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(fa_res_name2),file = paste0('fastANCOM_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('fastANCOM FINISHED')
  ANCOMBC_res_name3 <- paste0("ANCOMBC_res_",normalization_method)
  assign(ANCOMBC_res_name3,ANCOMBC_analysis(simulated_data_all,covariates,normalization=normalization_method))
  write.table(get(ANCOMBC_res_name3),file = paste0('ANCOMBC_res_',data_index,"_",normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('ANCOMBC FINISHED')
  maaslin2_lm_LOG_res_name2 <- paste0("maaslin2_lm_LOG_res_",normalization_method)
  assign(maaslin2_lm_LOG_res_name2,Maaslin2_LM_analysis(simulated_data_all,covariates,normalization=normalization_method,transform = 'LOG'))
  write.table(get(maaslin2_lm_LOG_res_name2),file = paste0('maaslin2_lm_LOG_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  print('maaslin2-lm LOG FINISHED')
  ################################TSS
  # print('TSS start!')
  # normalization_method = 'TSS'
  # ZicoSeq_res_name2 <- paste0("ZicoSeq_res_",normalization_method)
  # assign(ZicoSeq_res_name2,ZicoSeq_analysis(simulated_data_all,covariates,normalization_method))
  # write.table(get(ZicoSeq_res_name2),file = paste0('ZicoSeq_res_',data_index,'_',normalization_method,'.tsv'),sep = '\t',row.names = F,col.names = T,quote = F)
  # print('ZicoSeq FINISHED')
  
  
  all_result<-list(get(fa_res_name1),get(fa_res_name2),get(ANCOMBC_res_name1),get(ANCOMBC_res_name2),get(ANCOMBC_res_name3),get(lfem_res_name),
                  get(maaslin2_lm_LOG_res_name1),get(maaslin2_lm_LOG_res_name2),get(edgeR_res_name),get(limma_res_name),get(lmem_res_name))
  
  return(all_result)
}

