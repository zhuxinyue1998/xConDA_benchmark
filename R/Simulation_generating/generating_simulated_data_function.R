#####导入训练好的相应疾病模拟数据生成
fitted_ASD <- readRDS("/home/data/ZXY/meta_causal/benchmark_DA_pipeline_deep_ASD/fitted_ASD_WGS.rds")
##绝对丰度转化相对丰度函数
abd_to_rel<-function(feat_abd){
  mat <- apply(feat_abd,2, as.numeric)
  rel_abd<-t(apply(mat,1, function(x){
    if(sum(x) == 0){
      return(x)
    }else{
      return(x/sum(x))
    }
  }))
  rownames(rel_abd)<-rownames(feat_abd)
  return(rel_abd)
}
####提取ground_truth函数
extract_ground_truth<-function(simulation,simulated_data){
  spike_metadata<-simulation$spike_metadata$feature_metadata_spike_df
  spike_metadata <- spike_metadata[spike_metadata$feature_spiked %in% rownames(simulated_data),]
  ground_truth <- spike_metadata
  return(ground_truth)
}
####随机抽样函数
sample_selection<-function(metadata,case_sample_num,control_sample_num){
  ###output selected samples
  case <- rownames(metadata)[metadata[,1]==1]
  control <- rownames(metadata)[metadata[,1]==0]
  case_selected <- sample(case,size = case_sample_num,replace = F)
  control_selected <- sample(control,size = control_sample_num,replace = F)
  selected_sample <- append(case_selected,control_selected)
  
  return(selected_sample)
}
####输出模拟数据函数
output_simulated_data<-function(metadata,simulated_data,ground_truth,data_index,case_sample_num,control_sample_num){

  selected_sample<-sample_selection(metadata,case_sample_num,control_sample_num)
  metadata <- metadata[selected_sample,]
  simulated_data_selected_abs <- simulated_data[,selected_sample]
  simulated_data_selected <- abd_to_rel(t(simulated_data_selected_abs))
  metadata <- as.data.frame(metadata)
  metadata$SampleID <- rownames(metadata)
  write.table(ground_truth,file = paste0('ground_truth',data_index,'.tsv'),row.names = F,col.names = T,quote = F,sep = '\t')
  write.table(simulated_data_selected,file = paste0('simulated_data_',data_index,'.tsv'),row.names = T,col.names = T,sep = "\t",quote = F)
  write.table(t(simulated_data_selected_abs),file = paste0('simulated_data_abs_',data_index,'.tsv'),row.names = T,col.names = T,sep = "\t",quote = F)
  write.table(metadata,file = paste0('simulated_metadata',data_index,'.tsv'),row.names = F,col.names = T,quote = F,sep = '\t')
  result = list(ground_truth = ground_truth,
                simulated_data_rela = simulated_data_selected,
                simulated_data_abs = t(simulated_data_selected_abs),
                metadata = metadata)
  return(result)
}
################生成模拟数据函数包装
##协变量生成的参数
params_df <- data.frame(
  Variable = c("Con1", "Con2", "Con3", "Con4", "Ca1", "Ca2", "Ca3", "Ca4"),
  Group = c(rep("Con", 4), rep("Ca", 4)),
  Mean = c(2, 1.5, 3, 2.5, NA, NA, NA, NA),
  SD = c(0.5, 0.3, 0.8, 0.7, NA, NA, NA, NA),
  Prob0_0 = c(NA, NA, NA, NA, 0.4, 0.3, 0.55, 0.75),
  Prob0_1 = c(NA, NA, NA, NA, 0.6, 0.7, 0.45, 0.25),
  Prob1_0 = c(NA, NA, NA, NA, 0.7, 0.8, 0.35, 0.15),
  Prob1_1 = c(NA, NA, NA, NA, 0.3, 0.2, 0.65, 0.85)
)
####协变量生成函数
generate_covariates_custom <- function(sample_num = 10000, 
                                       group_ratio = 0.5, 
                                       pattern = "1Con+2Ca", 
                                       con_params = list(), 
                                       ca_params = list()) {
  # 检查参数合法性
  if (group_ratio <= 0 || group_ratio >= 1) {
    stop("group_ratio 必须在 (0, 1) 之间")
  }
  
  # 解析 pattern，提取 Con 和 Ca 的数量
  matches <- regmatches(pattern, regexec("(\\d+)Con\\+(\\d+)Ca", pattern))[[1]]
  if (length(matches) < 3) {
    stop("模式输入不正确，正确格式应为如 '1Con+2Ca'")
  }
  num_con <- as.numeric(matches[2])
  num_ca <- as.numeric(matches[3])
  
  # 生成组别信息
  num_group0 <- round(sample_num * group_ratio)
  num_group1 <- sample_num - num_group0
  Group <- c(rep(0, num_group0), rep(1, num_group1))
  
  # 生成数值型协变量（Con）
  Con_list <- list()
  for (i in seq_len(num_con)) {
    param <- con_params[[i]]
    if (is.null(param)) {
      stop(paste("缺少第", i, "个 Con 的参数，请在 con_params 中提供"))
    }
    #set.seed(param$seed)
    Con <- ifelse(Group == 0, 
                  rnorm(num_group0, mean = param$mean0, sd = param$sd0), 
                  rnorm(num_group1, mean = param$mean1, sd = param$sd1))
    Con_list[[paste0("Con", i)]] <- Con
  }
  
  # 生成分类型协变量（Ca）
  Ca_list <- list()
  for (i in seq_len(num_ca)) {
    param <- ca_params[[i]]
    if (is.null(param)) {
      stop(paste("缺少第", i, "个 Ca 的参数，请在 ca_params 中提供"))
    }
    #set.seed(param$seed)
    Ca <- sapply(Group, function(x) {
      if (x == 0) {
        sample(c(0, 1), 1, prob = c(param$prob0[1], param$prob0[2]))
      } else {
        sample(c(0, 1), 1, prob = c(param$prob1[1], param$prob1[2]))
      }
    })
    Ca_list[[paste0("Ca", i)]] <- Ca
  }
  
  # 组合结果为数据框
  metadata <- data.frame(Group = Group, do.call(cbind, Con_list), do.call(cbind, Ca_list))
  return(metadata)
}

# 解析模式并从数据框中选择对应的变量
select_params <- function(params_df, pattern) {
  # 解析模式
  matches <- regmatches(pattern, regexec("(\\d+)Con\\+(\\d+)Ca", pattern))[[1]]
  if (length(matches) < 3) {
    stop("模式输入不正确，正确格式应为如 '1Con+2Ca'")
  }
  num_con <- as.numeric(matches[2])
  num_ca <- as.numeric(matches[3])
  
  # 提取连续变量参数
  con_params <- params_df[params_df$Group == "Con", ][1:num_con, ]
  con_list <- lapply(1:nrow(con_params), function(i) {
    list(
      mean0 = con_params$Mean[i], 
      sd0 = con_params$SD[i], 
      mean1 = con_params$Mean[i] - 1, # 假设 Group 1 的均值为 Group 0 的均值 -1
      sd1 = con_params$SD[i] + 0.2  # 假设 Group 1 的标准差为 Group 0 的标准差 +0.1
    )
  })
  
  # 提取分类变量参数
  ca_params <- params_df[params_df$Group == "Ca", ][1:num_ca, ]
  ca_list <- lapply(1:nrow(ca_params), function(i) {
    list(
      prob0 = c(ca_params$Prob0_0[i], ca_params$Prob0_1[i]),
      prob1 = c(ca_params$Prob1_0[i], ca_params$Prob1_1[i])
    )
  })
  
  list(con_params = con_list, ca_params = ca_list)
}
Simulation<-function(primary_effect,confounding_effect,confounder_characteristics,data_index){
    library(SparseDOSSA2)
    #不指定ground-truth,不对ground-truth进行生成
    sample_num = 10000
    #注意：组别和协变量需要存在关联才是confounder
  
    # 设置全局种子
    set.seed(177)  # 设置全局种子，确保数据生成可重复
    ####大多数种子设置为1015，少数种子设置为1999（比如5rep2,19之后）
    # 使用函数提取参数
    params <- select_params(params_df, confounder_characteristics)
    con_params <- params$con_params
    ca_params <- params$ca_params
    # 使用参数生成协变量
    metadata_matrix <- as.matrix(generate_covariates_custom(
        sample_num = sample_num, 
        group_ratio = 0.5, 
        pattern = confounder_characteristics, 
        con_params = con_params, 
        ca_params = ca_params
    ))
    ASD_simulation <- SparseDOSSA2(template = fitted_ASD, 
                                 n_sample = sample_num, 
                                 new_features = TRUE,n_feature = 1000,median_read_depth = 3000000,metadata_effect_size = primary_effect,perc_feature_spiked_metadata = 0.2,
                                 spike_metadata = "abundance",metadata_matrix = metadata_matrix,
                                 verbose = TRUE)
    spike_metadata <- ASD_simulation$spike_metadata$feature_metadata_spike_df
    ###设置协变量的效应值为2
    spike_metadata$effect_size[spike_metadata$metadata_datum != 1] <- ifelse(runif(sum(spike_metadata$metadata_datum != 1)) > 0.65, confounding_effect, -confounding_effect)
    ASD_simulation <- SparseDOSSA2(template = fitted_ASD, 
                                 n_sample = sample_num, 
                                 new_features = TRUE,n_feature = 1000,median_read_depth = 3000000,spike_metadata = spike_metadata,
                                 metadata_matrix = metadata_matrix,
                                 verbose = TRUE)
    simulated_data <- ASD_simulation$simulated_data
    metadata <- ASD_simulation$spike_metadata$metadata_matrix
    rownames(metadata) <- colnames(simulated_data)
    simulated_data <- simulated_data[which((rowSums(simulated_data!=0)/sample_num)>0.1),]
    simulated_data <- simulated_data[,which(colSums(simulated_data!=0)>10)]
    metadata <- metadata[colnames(simulated_data),]
    
    ground_truth <- extract_ground_truth(ASD_simulation,simulated_data)
    simulation_list <- list(metadata = metadata,
                          simulated_data = simulated_data,
                          ground_truth = ground_truth,
                          data_index = data_index)
    return(simulation_list)
}