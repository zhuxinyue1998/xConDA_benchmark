#CLR标准化函数
library(compositions)
AST_transform <- function(x) {
  y <- sign(x) * asin(sqrt(abs(x)))
  if(any(is.na(y))) {
    logging::logerror(
      paste0("AST transform is only valid for values between -1 and 1. ",
             "Please select an appropriate normalization option or ",
             "normalize your data prior to running."))
    stop()
  }
  return(y)
}
clr_transform<-function(data){
  #要求输入行是样本，列是特征,count
  data_clr <- data.frame(t(compositions::clr(t(data))))
  return(data_clr)
}
#CSS标准化函数
library(metagenomeSeq)
css_transform<-function(data,log_transform=TRUE){
  #要求输入行是样本，列是特征,count
  metaSeqObject = newMRexperiment(t(data))
  metaSeqObject_CSS <- cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) )
  data_CSS = data.frame(t(data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=log_transform))))
  #确定是否返回log转换的数据
  return(data_CSS)
}
#TMM标准化函数
TMM_transform <- function(features) {
  #输入列是特征，行是样本，count
  # Convert to Matrix from Data Frame
  features_norm = as.matrix(features)
  dd <- colnames(features_norm)
  
  # TMM Normalizing the Data
  X <- t(features_norm)
  
  libSize = edgeR::calcNormFactors(X, method = "TMM")
  eff.lib.size = colSums(X) * libSize
  
  ref.lib.size = mean(eff.lib.size)
  #Use the mean of the effective library sizes as a reference library size
  X.output = sweep(X, MARGIN = 2, eff.lib.size, "/") * ref.lib.size
  #Normalized read counts
  
  # Convert back to data frame
  features_TMM <- as.data.frame(t(X.output))
  
  # Rename the True Positive Features - Same Format as Before
  colnames(features_TMM) <- dd
  
  
  # Return as list
  return(features_TMM)
}
##log转换函数
LOG_transform <- function(x) {
  y <- replace(x, x == 0, min(x[x>0]) / 2)
  return(data.frame(log2(y)))
}
##输出指定的标准化数据函数
data_normalization <- function(simulated_data_all,normalization){
  if(normalization == 'CLR'){
  normalized_data<-clr_transform(simulated_data_all$simulated_data_abs)
  }
  if(normalization == 'TMM'){
  normalized_data<-TMM_transform(simulated_data_all$simulated_data_abs)
  }
  if(normalization == 'CSS'){
  normalized_data<-css_transform(simulated_data_all$simulated_data_abs)
  }
  if(normalization == 'LOG'){
  normalized_data<-LOG_transform(simulated_data_all$simulated_data_abs)
  }
  if(normalization == 'count'){
  normalized_data<-simulated_data_all$simulated_data_abs
  }
  if(normalization == 'TSS'){
  normalized_data<-simulated_data_all$simulated_data_rela
  }
  if(normalization == 'TSS+AST'){
    normalized_data<-AST_transform(simulated_data_all$simulated_data_rela)
  }
  return(normalized_data)
}
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
