#####文件原名method_filter.py
import os
import statistics
import pandas as pd
import numpy as np

def weighted_score_with_stability(metrics_evaluation_result_mean, metrics_evaluation_result_sd, gamma=1.0):
    """
    计算方法的加权评分，同时考虑稳定性（标准差的影响）。
    
    - 采用 Min-Max 归一化，使各指标范围归一化到 [0,1]。
    - 对 FPR 和 FDR 进行 (1 - 值) 变换，以统一方向（较好的值趋向 1）。
    - 计算得分时，根据标准差 (SD) 调整权重，标准差越大，该指标的影响越小。
    
    :param metrics_evaluation_result_mean: DataFrame, 各方法的指标均值
    :param metrics_evaluation_result_sd: DataFrame, 各方法的指标标准差
    :param gamma: float, 控制标准差对权重的折损影响, 默认 1.0
    :return: DataFrame, 包含 Method 和 最终 Score
    """

    # 复制数据，避免修改原数据
    df_mean = metrics_evaluation_result_mean.copy()
    df_sd = metrics_evaluation_result_sd.copy()

    # 反转 FPR 和 FDR（1 - 值），使更好情况接近 1
    df_mean['FPR'] = 1 - df_mean['FPR']
    df_mean['FDR'] = 1 - df_mean['FDR']

    # 需要归一化的指标
    metrics = ['FPR', 'FDR', 'Accuracy', 'AUC', 'AUPR', 'Sensitivity']
    
    # 对 mean 进行 min-max 归一化（不处理 sd）
    for metric in metrics:
        min_val = df_mean[metric].min()
        max_val = df_mean[metric].max()
        if max_val > min_val:
            df_mean[metric] = (df_mean[metric] - min_val) / (max_val - min_val)
        else:
            df_mean[metric] = 0.0  # 避免除零错误
    
    # 指标权重
    weights = {
        'AUC': 2,
        'AUPR': 4,
        'Sensitivity': 4,
        'FPR': 1,
        'FDR': 3,
        'Accuracy': 0
    }

    # 计算每个方法的最终得分（折损权重计算）
    final_scores = []
    for i in range(len(df_mean)):
        row_score = 0
        for metric, weight in weights.items():
            weight_adjusted = weight * (1 - gamma * df_sd.loc[i, metric])  # 根据标准差调整权重
            weight_adjusted = max(0, weight_adjusted)  # 确保权重不会变成负数
            row_score += df_mean.loc[i, metric] * weight_adjusted
        final_scores.append(row_score)

    # 构建结果 DataFrame
    result_df = pd.DataFrame({
        'Method': df_mean['Method'],
        'Score': final_scores
    })

    return result_df


def metrics_filter_threshold(metrics_evaluation_result,metrics_evaluation_result_mean,metrics_evaluation_result_sd):
    ###提取数据并计算阈值
    FPR = metrics_evaluation_result['FPR'].tolist()
    FDR = metrics_evaluation_result['FDR'].tolist()
    Accuracy = metrics_evaluation_result['Accuracy'].tolist()
    AUROC = metrics_evaluation_result['AUC'].tolist()
    AUPR = metrics_evaluation_result['AUPR'].tolist()
    Sensitivity = metrics_evaluation_result['Sensitivity'].tolist()
    #F1 = metrics_evaluation_result['F1'].tolist()
    FPR_sd = metrics_evaluation_result_sd['FPR'].tolist()
    FDR_sd = metrics_evaluation_result_sd['FDR'].tolist()
    Accuracy_sd = metrics_evaluation_result_sd['Accuracy'].tolist()
    AUROC_sd = metrics_evaluation_result_sd['AUC'].tolist()
    AUPR_sd = metrics_evaluation_result_sd['AUPR'].tolist()
    Sensitivity_sd = metrics_evaluation_result_sd['Sensitivity'].tolist()
    #F1_sd = metrics_evaluation_result_sd['F1'].tolist()
    #median_FPR_value = statistics.median(FPR)
    #median_FDR_value = statistics.median(FDR)
    #median_Accuracy_value = statistics.median(Accuracy)
    #median_AUROC_value = statistics.median(AUROC)
    #median_Sensitivity_value = statistics.median(Sensitivity)
    q25_FPR_value = statistics.quantiles(FPR, n=4)[0]
    q25_FDR_value = statistics.quantiles(FDR, n=4)[0]
    q75_FPR_value = statistics.quantiles(FPR, n=4)[2]
    q75_FDR_value = statistics.quantiles(FDR, n=4)[2]
    q75_Accuracy_value = statistics.quantiles(Accuracy, n=4)[2]
    q75_AUROC_value = statistics.quantiles(AUROC, n=4)[2]
    q75_AUPR_value = statistics.quantiles(AUPR, n=4)[2]
    q75_Sensitivity_value = statistics.quantiles(Sensitivity, n=4)[2]
    q50_Sensitivity_value = statistics.quantiles(Sensitivity, n=4)[1]
    #q75_F1_value = statistics.quantiles(F1,n=4)[2]

    Method = metrics_evaluation_result_mean['Method'].tolist()
    FPR = metrics_evaluation_result_mean['FPR'].tolist()
    FDR = metrics_evaluation_result_mean['FDR'].tolist()
    Accuracy = metrics_evaluation_result_mean['Accuracy'].tolist()
    AUROC = metrics_evaluation_result_mean['AUC'].tolist()
    Sensitivity = metrics_evaluation_result_mean['Sensitivity'].tolist()
    #F1 = metrics_evaluation_result_mean['F1'].tolist()
    AUPR = metrics_evaluation_result_mean['AUPR'].tolist()

    elimination_FPR_method = [Method[i] for i in range(len(Method)) if FPR[i] > q25_FPR_value]
    elimination_FDR_method = [Method[i] for i in range(len(Method)) if FDR[i] > q25_FDR_value]
    elimination_Accuracy_method = [Method[i] for i in range(len(Method)) if Accuracy[i] < q75_Accuracy_value]
    elimination_AUROC_method = [Method[i] for i in range(len(Method)) if AUROC[i] < q75_AUROC_value]
    elimination_AUPR_method = [Method[i] for i in range(len(Method)) if AUPR[i] < q75_AUPR_value]
    elimination_Sensitivity_method = [Method[i] for i in range(len(Method)) if Sensitivity[i] < q75_Sensitivity_value]
    #elimination_F1_method = [Method[i] for i in range(len(Method)) if F1[i] < q75_F1_value]
    elimination_AUROC_p50_method = [Method[i] for i in range(len(Method)) if AUROC[i] <= 0.5]
    elimination_AUPR_p30_method = [Method[i] for i in range(len(Method)) if AUPR[i] <= 0.3]
    if q50_Sensitivity_value <= 0.2:
        elimination_low_Sensitivity_method = [Method[i] for i in range(len(Method)) if Sensitivity[i] < 0.2]
    else:
        elimination_low_Sensitivity_method = [Method[i] for i in range(len(Method)) if Sensitivity[i] < q50_Sensitivity_value]
    elimination_high_FDR_method = [Method[i] for i in range(len(Method)) if FDR[i] > q75_FDR_value]
    elimination_high_FPR_method = [Method[i] for i in range(len(Method)) if FPR[i] > q75_FPR_value]
    #combined_list = elimination_FPR_method * 2 + elimination_FDR_method * 2 + elimination_Accuracy_method + elimination_AUROC_method + elimination_Sensitivity_method * 2 + elimination_AUPR_method * 2 + elimination_AUROC_p50_method * 10 + elimination_AUPR_p30_method * 10 + elimination_low_Sensitivity_method * 8
    # 计算每个指标的计数
    counts_FPR = [elimination_FPR_method.count(method) for method in Method]
    counts_FDR = [elimination_FDR_method.count(method) for method in Method]
    counts_Accuracy = [elimination_Accuracy_method.count(method) for method in Method]
    counts_AUROC = [elimination_AUROC_method.count(method) for method in Method]
    counts_AUPR = [elimination_AUPR_method.count(method) for method in Method]
    counts_Sensitivity = [elimination_Sensitivity_method.count(method) for method in Method]
    counts_p50_AUROC = [elimination_AUROC_p50_method.count(method) for method in Method]
    counts_p30_AUPR = [elimination_AUPR_p30_method.count(method) for method in Method]
    counts_low_Sensitivity = [elimination_low_Sensitivity_method.count(method) for method in Method]
    counts_high_FDR = [elimination_high_FDR_method.count(method) for method in Method]
    counts_high_FPR = [elimination_high_FPR_method.count(method) for method in Method]
    # 按权重计算每个方法的最终分数
    scores = [
        counts_FPR[i] * 1 + (FPR_sd[i] * counts_FPR[i]) * 1 +
        counts_FDR[i] * 1 + (FDR_sd[i] * counts_FDR[i]) * 1 +
        counts_high_FPR[i] * 4 + (FPR_sd[i] * counts_high_FPR[i]) * 4 +
        counts_high_FDR[i] * 4 + (FDR_sd[i] * counts_high_FDR[i]) * 4 +
        counts_Accuracy[i] * 1 + (Accuracy_sd[i] * counts_Accuracy[i]) * 1 +
        counts_AUROC[i] * 3 + (AUROC_sd[i] * counts_AUROC[i]) * 3 +
        counts_AUPR[i] * 5 + (AUPR_sd[i] * counts_AUPR[i]) * 5 +
        counts_Sensitivity[i] * 2 + (Sensitivity_sd[i] * counts_Sensitivity[i]) * 2 + 
        counts_p50_AUROC[i] * 7 + counts_p30_AUPR[i] * 5 + 
        counts_low_Sensitivity[i] * 8 + (Sensitivity_sd[i] * counts_low_Sensitivity[i]) * 8
        for i in range(len(Method))
    ]
    ###创建dataframe存储
    df = pd.DataFrame(
        {
            'Method':Method,
            'Scores':scores
        }
    )
    return df


def merge_count_result(base_directory):
    metrics_file_path = os.path.join(base_directory,'metrics_evaluation_result_0.05.tsv')
    metrics_mean_file_path = os.path.join(base_directory,'metrics_evaluation_result_0.05_mean.tsv')
    metrics_sd_file_path = os.path.join(base_directory,'metrics_evaluation_result_0.05_sd.tsv')
    metrics_file = pd.read_csv(metrics_file_path,sep='\t')
    metrics_mean_file = pd.read_csv(metrics_mean_file_path,sep='\t')
    metrics_sd_file = pd.read_csv(metrics_sd_file_path,sep='\t')
    count_df = weighted_score_with_stability(metrics_mean_file,metrics_sd_file)
        # if i == 0:
        #     count_all = count_df
        #     count_all.index = count_df['Method']
        #     #count_all = count_all.drop(columns='Method')
        # else:
        #     count_df.index = count_df['Method']
        #     count_df = count_df.drop(columns='Method')
        #     count_all = pd.concat([count_all, count_df], axis=1)
    return count_df

