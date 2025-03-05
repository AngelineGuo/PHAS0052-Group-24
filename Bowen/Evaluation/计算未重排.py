import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def calculate_metrics(TP, TN, FP, FN):
    """ Sensitivity (Sen), Specificity (Spe), Accuracy (Acc), Precision (Pre), MCC """
    Sen = TP / (TP + FN) if (TP + FN) > 0 else 0
    Spe = TN / (TN + FP) if (TN + FP) > 0 else 0
    Acc = (TP + TN) / (TP + FN + TN + FP) if (TP + FN + TN + FP) > 0 else 0
    Pre = TP / (TP + FP) if (TP + FP) > 0 else 0
    MCC_numerator = (TP * TN) - (FP * FN)
    MCC_denominator = ((TP + FN) * (TP + FP) * (TN + FN) * (TN + FP)) ** 0.5
    MCC = MCC_numerator / MCC_denominator if MCC_denominator > 0 else 0

    return [Sen, Spe, Acc, Pre, MCC]


def calculate_brier_score_and_metrics(actual_file, predicted_file, prob_col):
    df_actual = pd.read_csv(actual_file)
    df_predicted = pd.read_csv(predicted_file)
    df_merged = pd.concat([df_actual.reset_index(drop=True), df_predicted.reset_index(drop=True)], axis=1)

    df_merged["ATP_binding_actual"] = df_merged["ATP Binding Site"].map({"B": 1, "N": 0})
    df_merged["pred_prob"] = df_merged[prob_col].astype(float)
    df_merged["pred_label"] = (df_merged["pred_prob"] >= 0.5).astype(int)

    TP = ((df_merged["ATP_binding_actual"] == 1) & (df_merged["pred_label"] == 1)).sum()
    TN = ((df_merged["ATP_binding_actual"] == 0) & (df_merged["pred_label"] == 0)).sum()
    FP = ((df_merged["ATP_binding_actual"] == 0) & (df_merged["pred_label"] == 1)).sum()
    FN = ((df_merged["ATP_binding_actual"] == 1) & (df_merged["pred_label"] == 0)).sum()

    brier_score = ((df_merged["pred_prob"] - df_merged["ATP_binding_actual"]) ** 2).mean()
    metrics = calculate_metrics(TP, TN, FP, FN)

    return [brier_score] + metrics + [TP, FP, TN, FN]


def process_folders(actual_folder, predicted_folder, output_csv):
    results = []
    actual_files = {f[:4]: f for f in os.listdir(actual_folder)}
    predicted_files = {f[:4]: f for f in os.listdir(predicted_folder)}

    for key in actual_files.keys() & predicted_files.keys():
        actual_path = os.path.join(actual_folder, actual_files[key])
        predicted_path = os.path.join(predicted_folder, predicted_files[key])

        if os.path.exists(actual_path) and os.path.exists(predicted_path):
            brier_bind, *metrics_bind, TP_bind, FP_bind, TN_bind, FN_bind = calculate_brier_score_and_metrics(
                actual_path, predicted_path, "ATPbind Probability")
            brier_seq, *metrics_seq, TP_seq, FP_seq, TN_seq, FN_seq = calculate_brier_score_and_metrics(
                actual_path, predicted_path, "ATPseq Probability")

            results.append([actual_files[key]] + [brier_bind] + metrics_bind + [TP_bind, FP_bind, TN_bind, FN_bind] +
                           [brier_seq] + metrics_seq + [TP_seq, FP_seq, TN_seq, FN_seq])

    columns = ["Filename", "Brier Score (ATPbind)", "Sen (ATPbind)", "Spe (ATPbind)", "Acc (ATPbind)", "Pre (ATPbind)",
               "MCC (ATPbind)",
               "TP (ATPbind)", "FP (ATPbind)", "TN (ATPbind)", "FN (ATPbind)", "Brier Score (ATPseq)", "Sen (ATPseq)",
               "Spe (ATPseq)",
               "Acc (ATPseq)", "Pre (ATPseq)", "MCC (ATPseq)", "TP (ATPseq)", "FP (ATPseq)", "TN (ATPseq)",
               "FN (ATPseq)"]
    df_results = pd.DataFrame(results, columns=columns)
    df_results.to_csv(output_csv, index=False)
    return df_results


def plot_metrics_separately(df_results, save_folder="plots"):
    metrics = ["Brier Score", "Sen", "Spe", "Acc", "Pre", "MCC"]

    # 创建保存目录
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)

    for metric in metrics:
        plt.figure(figsize=(7, 4))
        df_melted = df_results.melt(id_vars=["Filename"],
                                    value_vars=[f"{metric} (ATPbind)", f"{metric} (ATPseq)"],
                                    var_name="Model", value_name=metric)
        sns.stripplot(data=df_melted, x="Model", y=metric, hue="Model", jitter=True, alpha=0.7)

        mean_val = df_melted[metric].mean()
        std_val = df_melted[metric].std()
        plt.axhline(mean_val, color='r', linestyle='--', label=f'Mean: {mean_val:.2f}')
        plt.axhline(mean_val + std_val, color='g', linestyle=':', label=f'+1 Std: {mean_val + std_val:.2f}')
        plt.axhline(mean_val - std_val, color='g', linestyle=':', label=f'-1 Std: {mean_val - std_val:.2f}')

        plt.title(f"{metric} Distribution for ATPseq and ATPbind")
        plt.xlabel("Model")
        plt.ylabel(metric)
        plt.legend()

        # 保存图片
        save_path = os.path.join(save_folder, f"{metric.replace(' ', '_')}.png")
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()  # 关闭图形，释放内存
        print(f"Saved plot: {save_path}")



if __name__ == "__main__":
    actual_folder = r"C:\\Users\\ybw14\\Desktop\\第二部分\\28 ATP 结合点位\\结合点位 4A pdb"
    predicted_folder = r"C:\\Users\\ybw14\\Desktop\\大三项目\\predicted_result"
    output_csv = "evaluation 未重排.csv"

    df_results = process_folders(actual_folder, predicted_folder, output_csv)
    print(f"Brier Scores and Metrics saved to {output_csv}")
    plot_metrics_separately(df_results)
