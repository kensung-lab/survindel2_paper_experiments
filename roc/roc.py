import argparse
import pysam
import matplotlib.pyplot as plt

def parse_equivalents(file_path):
    equivalents = {}
    with open(file_path, 'r') as file:
        for line in file:
            benchmark_id, caller_id = line.strip().split()
            equivalents[caller_id] = benchmark_id
    return equivalents

def compute_metrics(caller_vcfs, benchmark_vcf, equivalents_list):
    benchmark_sv = {rec.id: rec for rec in pysam.VariantFile(benchmark_vcf)}
    metrics_list = []

    for caller_vcf, equivalents in zip(caller_vcfs, equivalents_list):
        caller_sv = {rec.id: rec for rec in pysam.VariantFile(caller_vcf)}
        inverted_equivalents = {v: k for k, v in equivalents.items()}
        
        thresholds = sorted(set(rec.qual for rec in caller_sv.values() if rec.qual is not None), reverse=True)
        if not thresholds:
            thresholds = [0]
        metrics = []
        
        for threshold in thresholds:
            filtered_caller_sv = {sv_id: rec for sv_id, rec in caller_sv.items() if threshold == 0 or rec.qual >= threshold}
            tp = sum(1 for sv_id in benchmark_sv if sv_id in inverted_equivalents and inverted_equivalents[sv_id] in filtered_caller_sv)
            fp = sum(1 for sv_id in filtered_caller_sv if sv_id not in equivalents)
            fn = sum(1 for sv_id in benchmark_sv if sv_id not in inverted_equivalents or (sv_id in inverted_equivalents and inverted_equivalents[sv_id] not in filtered_caller_sv))
            
            precision = tp / (tp + fp) if tp + fp > 0 else 0
            recall = tp / (tp + fn) if tp + fn > 0 else 0
            
            metrics.append((recall, precision))
        
        metrics_list.append(metrics)
    
    return metrics_list

def plot_curves(metrics_list, labels, output_path):
    plt.figure(figsize=(10, 6))
    for metrics, label in zip(metrics_list, labels):
        recalls, precisions = zip(*metrics)
        plt.plot(recalls, precisions, marker='o', label=label)
    
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve')
    plt.legend()
    plt.grid(True)
    plt.savefig(output_path)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Plot Precision-Recall Curve based on QUAL threshold for multiple variant callers.")
    parser.add_argument("--caller_vcfs", nargs='+', required=True, help="Paths to the VCF files from the callers.")
    parser.add_argument("--benchmark_vcf", required=True, help="Path to the benchmark VCF file.")
    parser.add_argument("--equivalents", nargs='+', required=True, help="Paths to the text files with equivalent SV IDs for each caller.")
    parser.add_argument("--labels", nargs='+', required=True, help="Labels for each caller's curve on the plot.")
    parser.add_argument("--output", required=True, help="Path to save the plot.")
    args = parser.parse_args()
    
    if len(args.caller_vcfs) != len(args.equivalents) or len(args.caller_vcfs) != len(args.labels):
        print("Error: The number of caller VCF files, equivalents files, and labels must match.")
        return
    
    equivalents_list = [parse_equivalents(file) for file in args.equivalents]
    metrics_list = compute_metrics(args.caller_vcfs, args.benchmark_vcf, equivalents_list)
    plot_curves(metrics_list, args.labels, args.output)

if __name__ == "__main__":
    main()
