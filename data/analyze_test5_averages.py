#!/usr/bin/env python3
"""
Visualization script for test5a and test5b results.
Reads all test5a and test5b CSV files and calculates averages.

Test5a: Single fidelity measurement
Test5b: Two fidelity measurements
  - fidelity_full: Full system fidelity
  - fidelity_partial_(qubits0-3): Partial trace fidelity for qubits 0-3
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from pathlib import Path

def read_test_files(pattern):
    """Read all CSV files matching the pattern and return a list of DataFrames."""
    files = glob.glob(pattern)
    dataframes = []
    filenames = []
    
    for file in sorted(files):
        df = pd.read_csv(file)
        dataframes.append(df)
        filenames.append(os.path.basename(file))
    
    return dataframes, filenames

def calculate_averages_test5a(dataframes, filenames):
    """Calculate averages for test5a files (single fidelity column)."""
    results = []
    
    for df, filename in zip(dataframes, filenames):
        avg_fidelity = df['fidelity'].mean()
        std_fidelity = df['fidelity'].std()
        results.append({
            'file': filename,
            'avg_fidelity': avg_fidelity,
            'std_fidelity': std_fidelity,
            'count': len(df)
        })
    
    return pd.DataFrame(results)

def calculate_averages_test5b(dataframes, filenames):
    """
    Calculate averages for test5b files.
    - fidelity_full: Full system fidelity
    - fidelity_partial_(qubits0-3): Partial trace fidelity for qubits 0-3
    """
    results = []
    
    for df, filename in zip(dataframes, filenames):
        avg_full = df['fidelity_full'].mean()
        avg_partial = df['fidelity_partial_(qubits0-3)'].mean()
        std_full = df['fidelity_full'].std()
        std_partial = df['fidelity_partial_(qubits0-3)'].std()
        results.append({
            'file': filename,
            'avg_fidelity_full': avg_full,
            'avg_fidelity_partial_q0_3': avg_partial,
            'std_fidelity_full': std_full,
            'std_fidelity_partial_q0_3': std_partial,
            'count': len(df)
        })
    
    return pd.DataFrame(results)

def extract_parameter(filename):
    """Extract the parameter value from filename (e.g., '0.1' from 'out_test5a_0.1_timestamp.csv')."""
    parts = filename.split('_')
    for i, part in enumerate(parts):
        if part.replace('.', '').replace('e-', '').replace('-', '').isdigit():
            return float(part)
    return None

def plot_test5a_results(results_df):
    """Create visualization for test5a results."""
    # Extract parameters and sort
    results_df['parameter'] = results_df['file'].apply(extract_parameter)
    results_df = results_df.sort_values('parameter')
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.bar(range(len(results_df)), results_df['avg_fidelity'], 
           color='steelblue', alpha=0.7, edgecolor='black',
           yerr=results_df['std_fidelity'], capsize=5, ecolor='black')
    ax.set_xlabel('Test Parameter', fontsize=12, fontweight='bold')
    ax.set_ylabel('Average Fidelity', fontsize=12, fontweight='bold')
    ax.set_title('Test5a: Average Fidelity by Parameter', fontsize=14, fontweight='bold')
    ax.set_xticks(range(len(results_df)))
    ax.set_xticklabels([f"{p:.0e}" if p < 0.01 else f"{p:.2f}" 
                        for p in results_df['parameter']], rotation=45)
    ax.grid(axis='y', alpha=0.3)
    ax.set_ylim(0, 1.05)
    
    # Add value labels on bars
    for i, (idx, row) in enumerate(results_df.iterrows()):
        ax.text(i, row['avg_fidelity'] + row['std_fidelity'] + 0.02, 
                f"{row['avg_fidelity']:.4f}",
                ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    plt.savefig('data/test5a_averages.png', dpi=300, bbox_inches='tight')
    print("Saved: test5a_averages.png")
    
    return fig

def plot_test5b_results(results_df):
    """Create visualization for test5b results."""
    # Extract parameters and sort
    results_df['parameter'] = results_df['file'].apply(extract_parameter)
    results_df = results_df.sort_values('parameter')
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    x = np.arange(len(results_df))
    width = 0.35
    
    bars1 = ax.bar(x - width/2, results_df['avg_fidelity_full'], width,
                   label='Full System', color='coral', alpha=0.7, edgecolor='black',
                   yerr=results_df['std_fidelity_full'], capsize=5, ecolor='black')
    bars2 = ax.bar(x + width/2, results_df['avg_fidelity_partial_q0_3'], width,
                   label='Partial Trace (Qubits 0-3)', color='lightseagreen', 
                   alpha=0.7, edgecolor='black',
                   yerr=results_df['std_fidelity_partial_q0_3'], capsize=5, ecolor='black')
    
    ax.set_xlabel('Test Parameter', fontsize=12, fontweight='bold')
    ax.set_ylabel('Average Fidelity', fontsize=12, fontweight='bold')
    ax.set_title('Test5b: Full System vs Partial Trace (Qubits 0-3) Fidelity', 
                 fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([f"{p:.0e}" if p < 0.01 else f"{p:.2f}" 
                        for p in results_df['parameter']], rotation=45)
    ax.legend(fontsize=10, loc='best')
    ax.grid(axis='y', alpha=0.3)
    ax.set_ylim(0, 1.05)
    
    # Add value labels on bars
    for i, (idx, row) in enumerate(results_df.iterrows()):
        ax.text(i - width/2, row['avg_fidelity_full'] + row['std_fidelity_full'] + 0.02, 
                f"{row['avg_fidelity_full']:.3f}",
                ha='center', va='bottom', fontsize=8)
        ax.text(i + width/2, row['avg_fidelity_partial_q0_3'] + row['std_fidelity_partial_q0_3'] + 0.02, 
                f"{row['avg_fidelity_partial_q0_3']:.3f}",
                ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    plt.savefig('data/test5b_averages.png', dpi=300, bbox_inches='tight')
    print("Saved: test5b_averages.png")
    
    return fig

def plot_comparison(results_5a, results_5b):
    """Create a comparison plot between test5a and test5b."""
    # Sort both by parameter
    results_5a['parameter'] = results_5a['file'].apply(extract_parameter)
    results_5b['parameter'] = results_5b['file'].apply(extract_parameter)
    results_5a = results_5a.sort_values('parameter')
    results_5b = results_5b.sort_values('parameter')
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Plot test5a with blue line
    ax.errorbar(results_5a['parameter'], results_5a['avg_fidelity'], 
                 yerr=results_5a['std_fidelity'],
                 marker='o', markersize=8, linewidth=2.5, label='Test A', 
                 color='blue', capsize=5, capthick=2)
    
    # Plot test5b (partial trace qubits 0-3) with orange line
    ax.errorbar(results_5b['parameter'], results_5b['avg_fidelity_partial_q0_3'], 
                 yerr=results_5b['std_fidelity_partial_q0_3'],
                 marker='s', markersize=8, linewidth=2.5, label='Test B (Qubits 0-3)', 
                 color='orange', capsize=5, capthick=2)
    
    ax.set_xlabel('Probability', fontsize=14, fontweight='bold')
    ax.set_ylabel('Fitness', fontsize=14, fontweight='bold')
    ax.set_title('Fitness vs Probability: Test A and Test B (Qubits 0-3)', 
                  fontsize=15, fontweight='bold')
    ax.legend(fontsize=12, loc='best')
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    ax.set_ylim(0, 1.05)
    
    plt.tight_layout()
    plt.savefig('data/test5_comparison.png', dpi=300, bbox_inches='tight')
    print("Saved: test5_comparison.png")
    
    return fig

def print_summary_statistics(results_5a, results_5b):
    """Print detailed summary statistics."""
    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)
    
    if not results_5a.empty:
        print("\nTest5a Overall:")
        print(f"  Mean Fidelity (across all parameters): {results_5a['avg_fidelity'].mean():.4f}")
        print(f"  Std Dev: {results_5a['avg_fidelity'].std():.4f}")
        print(f"  Min: {results_5a['avg_fidelity'].min():.4f}")
        print(f"  Max: {results_5a['avg_fidelity'].max():.4f}")
    
    if not results_5b.empty:
        print("\nTest5b Overall:")
        print("  Full System:")
        print(f"    Mean Fidelity: {results_5b['avg_fidelity_full'].mean():.4f}")
        print(f"    Std Dev: {results_5b['avg_fidelity_full'].std():.4f}")
        print(f"    Min: {results_5b['avg_fidelity_full'].min():.4f}")
        print(f"    Max: {results_5b['avg_fidelity_full'].max():.4f}")
        print("\n  Partial Trace (Qubits 0-3):")
        print(f"    Mean Fidelity: {results_5b['avg_fidelity_partial_q0_3'].mean():.4f}")
        print(f"    Std Dev: {results_5b['avg_fidelity_partial_q0_3'].std():.4f}")
        print(f"    Min: {results_5b['avg_fidelity_partial_q0_3'].min():.4f}")
        print(f"    Max: {results_5b['avg_fidelity_partial_q0_3'].max():.4f}")
        print("\n  Average Improvement (Partial over Full):")
        improvement = results_5b['avg_fidelity_partial_q0_3'].mean() - results_5b['avg_fidelity_full'].mean()
        print(f"    Δ Fidelity: {improvement:+.4f}")

def main():
    """Main execution function."""
    print("=" * 70)
    print("Quantum Probabilistic Combinator - Test5 Results Analysis")
    print("=" * 70)
    print()
    
    # Get current directory
    script_dir = Path(__file__).parent
    
    # Process test5a files
    print("Processing Test5a files...")
    test5a_pattern = str(script_dir / "out_test5a_*.csv")
    dfs_5a, files_5a = read_test_files(test5a_pattern)
    
    if not dfs_5a:
        print("  ⚠ No test5a files found!")
        results_5a = pd.DataFrame()
    else:
        print(f"  Found {len(dfs_5a)} test5a files:")
        for f in files_5a:
            print(f"    - {f}")
        
        results_5a = calculate_averages_test5a(dfs_5a, files_5a)
        results_5a['parameter'] = results_5a['file'].apply(extract_parameter)
        print("\n  Test5a Averages:")
        display_5a = results_5a[['parameter', 'avg_fidelity']].rename(columns={'parameter': 'prob'})
        display_5a = display_5a.sort_values('prob')
        print(display_5a.to_string(index=False))
        print()
    
    # Process test5b files
    print("Processing Test5b files...")
    test5b_pattern = str(script_dir / "out_test5b_*.csv")
    dfs_5b, files_5b = read_test_files(test5b_pattern)
    
    if not dfs_5b:
        print("  ⚠ No test5b files found!")
        results_5b = pd.DataFrame()
    else:
        print(f"  Found {len(dfs_5b)} test5b files:")
        for f in files_5b:
            print(f"    - {f}")
        
        results_5b = calculate_averages_test5b(dfs_5b, files_5b)
        results_5b['parameter'] = results_5b['file'].apply(extract_parameter)
        print("\n  Test5b Averages:")
        display_5b = results_5b[['parameter', 'avg_fidelity_partial_q0_3']].rename(
            columns={'parameter': 'prob', 'avg_fidelity_partial_q0_3': 'avg_fidelity_partial'}
        )
        display_5b = display_5b.sort_values('prob')
        print(display_5b.to_string(index=False))
        print()
    
    # Create visualizations
    print("Creating visualizations...")
    if not results_5a.empty:
        plot_test5a_results(results_5a)
    if not results_5b.empty:
        plot_test5b_results(results_5b)
    if not results_5a.empty and not results_5b.empty:
        plot_comparison(results_5a, results_5b)
    
    # Print summary statistics
    if not results_5a.empty or not results_5b.empty:
        print_summary_statistics(results_5a, results_5b)
    
    # Save summary statistics to CSV
    if not results_5a.empty:
        results_5a.to_csv('data/test5a_summary.csv', index=False)
        print("\nSaved: test5a_summary.csv")
    if not results_5b.empty:
        results_5b.to_csv('data/test5b_summary.csv', index=False)
        print("Saved: test5b_summary.csv")
    
    print()
    print("=" * 70)
    print("Analysis complete!")
    print("=" * 70)
    
    plt.show()

if __name__ == "__main__":
    main()
