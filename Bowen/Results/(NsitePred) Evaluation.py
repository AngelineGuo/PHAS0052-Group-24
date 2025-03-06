import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns


def calculate_metrics(TP, TN, FP, FN):
    """ Sensitivity (Sen), Specificity (Spe), Accuracy (Acc), Precision (Pre), MCC """
    Sen = TP / (TP + FN) if (TP + FN) > 0 else 0
    Spe = TN / (TN + FP) if (TN + FP) > 0 else 0
    Acc = (TP + TN) / (TP + FN + TN + FP) if (TP + FN + TN + FP) > 0 else 0
    Pre = TP / (TP + FP) if (TP + FP) > 0 else 0
    MCC_numerator = (TP * TN) - (FP * FN)
    MCC_denominator = ((TP + FN) * (TP + FP) * (TN + FN) * (TN + FP)) ** 0.5
    MCC = MCC_numerator / MCC_denominator if MCC_denominator > 0 else 0

    return Sen, Spe, Acc, Pre, MCC


def calculate_brier_score_and_metrics(actual_file, predicted_file):
    df_actual = pd.read_csv(actual_file)
    df_predicted = pd.read_csv(predicted_file)

    # Ensure data is aligned by row index
    df_predicted = df_predicted.iloc[:len(df_actual)].copy()
    df_predicted["NO"] = df_actual["NO"].values  # Assign NO from actual to predicted

    # Merge actual and predicted data based on NO (which is now just row alignment)
    df_merged = pd.merge(df_actual, df_predicted, on="NO", suffixes=("_actual", "_pred"))

    # Standardize column names for ATP probability
    possible_prob_cols = [col for col in df_merged.columns if "ATP prob" in col]
    if possible_prob_cols:
        prob_col = possible_prob_cols[0]  # Select the first matching column
    else:
        raise KeyError("No ATP probability column found in predicted data")

    # Convert ATP binding site labels: "Y" and "B" both represent binding sites
    df_merged["ATP_binding_actual"] = df_merged["ATP Binding Site"].map({"Y": 1, "B": 1, "N": 0})
    df_merged["pred_prob"] = df_merged[prob_col].astype(float)
    df_merged["pred_label"] = (df_merged["pred_prob"] >= 0.5).astype(int)

    # Compute confusion matrix components
    TP = ((df_merged["ATP_binding_actual"] == 1) & (df_merged["pred_label"] == 1)).sum()
    TN = ((df_merged["ATP_binding_actual"] == 0) & (df_merged["pred_label"] == 0)).sum()
    FP = ((df_merged["ATP_binding_actual"] == 0) & (df_merged["pred_label"] == 1)).sum()
    FN = ((df_merged["ATP_binding_actual"] == 1) & (df_merged["pred_label"] == 0)).sum()

    # Compute Brier Score
    brier_score = ((df_merged["pred_prob"] - df_merged["ATP_binding_actual"]) ** 2).mean()

    # Compute other metrics
    Sen, Spe, Acc, Pre, MCC = calculate_metrics(TP, TN, FP, FN)

    return brier_score, Sen, Spe, Acc, Pre, MCC, TP, FP, TN, FN


def process_folders(actual_folder, predicted_folder, output_csv):
    results = []

    actual_files = {f[:4]: f for f in os.listdir(actual_folder)}
    predicted_files = {f[:4]: f for f in os.listdir(predicted_folder)}

    for key in actual_files.keys() & predicted_files.keys():
        actual_path = os.path.join(actual_folder, actual_files[key])
        predicted_path = os.path.join(predicted_folder, predicted_files[key])

        if os.path.exists(actual_path) and os.path.exists(predicted_path):
            brier_score, Sen, Spe, Acc, Pre, MCC, TP, FP, TN, FN = calculate_brier_score_and_metrics(
                actual_path, predicted_path)
            results.append([actual_files[key], brier_score, Sen, Spe, Acc, Pre, MCC, TP, FP, TN, FN])

    # Save results as CSV
    df_results = pd.DataFrame(results,
                              columns=["Filename", "Brier Score", "Sen", "Spe", "Acc", "Pre", "MCC", "TP", "FP", "TN",
                                       "FN"])
    df_results.to_csv(output_csv, index=False)
    return df_results




def plot_metrics(df_results, save_folder="plots"):
    metrics = ["Brier Score", "Sen", "Spe", "Acc", "Pre", "MCC"]  # 需要绘图的所有指标

    if not os.path.exists(save_folder):
        os.makedirs(save_folder)

    for metric in metrics:
        plt.figure(figsize=(7, 5))


        mean_value = df_results[metric].mean()
        std_value = df_results[metric].std()


        sns.stripplot(data=df_results, x=["NsitePred"] * len(df_results), y=metric, jitter=True, alpha=0.7)


        plt.axhline(mean_value, color='red', linestyle='dashed', label=f"Mean: {mean_value:.2f}")
        plt.axhline(mean_value + std_value, color='green', linestyle='dotted', label=f"+1 Std: {mean_value + std_value:.2f}")
        plt.axhline(mean_value - std_value, color='green', linestyle='dotted', label=f"-1 Std: {mean_value - std_value:.2f}")


        plt.xticks(ticks=[0], labels=["NsitePred"])

        plt.title(f"{metric} Distribution")
        plt.xlabel("Model")
        plt.ylabel(metric)
        plt.legend()

        save_path = os.path.join(save_folder, f"{metric.replace(' ', '_')}.png")
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"Saved plot: {save_path}")



if __name__ == "__main__":
    actual_folder = r"C:\Users\ybw14\Desktop\第二部分\28 ATP 结合点位\结合点位 4A cif"  # Path to experimental data
    predicted_folder = r"C:\Users\ybw14\Desktop\第二部分\28 ATP 结合点位\NsitePred result"  # Path to NsitePred results
    output_csv = "NsitePred_evaluation.csv"

    df_results = process_folders(actual_folder, predicted_folder, output_csv)
    print(f"Metrics saved to {output_csv}")

    # Plot
    plot_metrics(df_results)



