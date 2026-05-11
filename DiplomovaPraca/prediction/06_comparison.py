"""
=============================================================================
Model Comparison and Final Visualization
=============================================================================
Compare ARIMA, SVM, and CNN models across all horizons
=============================================================================
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import mean_squared_error, mean_absolute_error

# Set style
plt.style.use('default')
sns.set_palette("husl")

# =============================================================================
# CONFIGURATION
# =============================================================================
RESULTS_BASE = Path(__file__).parent / "results"
ARIMA_DIR = RESULTS_BASE / "arima"
SVM_DIR = RESULTS_BASE / "svm"
CNN_DIR = RESULTS_BASE / "cnn"
LIGHTGBM_DIR = RESULTS_BASE / "lightgbm"
LSTM_DIR = RESULTS_BASE / "lstm_tcn"
COMPARISON_DIR = RESULTS_BASE / "comparison"
COMPARISON_DIR.mkdir(parents=True, exist_ok=True)

HORIZONS = ['day', 'week', 'month']
MODELS = ['ARIMA', 'SVM', 'CNN', 'LightGBM', 'LSTM', 'TCN']

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
def load_model_results(model_name, horizon):
    """Load results for a specific model and horizon"""
    if model_name == 'ARIMA':
        metrics_file = ARIMA_DIR / 'arima_metrics.csv'
        pred_file = ARIMA_DIR / f'arima_forecast_{horizon}.csv'
    elif model_name == 'SVM':
        metrics_file = SVM_DIR / 'svm_metrics.csv'
        pred_file = SVM_DIR / f'svm_predictions_{horizon}.csv'
    elif model_name == 'CNN':
        metrics_file = CNN_DIR / 'cnn_metrics.csv'
        pred_file = CNN_DIR / f'cnn_predictions_{horizon}.csv'
    elif model_name == 'LightGBM':
        metrics_file = LIGHTGBM_DIR / 'lightgbm_metrics.csv'
        pred_file = LIGHTGBM_DIR / f'lightgbm_predictions_{horizon}.csv'
    elif model_name == 'LSTM':
        metrics_file = LSTM_DIR / 'lstm_metrics.csv'
        pred_file = LSTM_DIR / f'lstm_predictions_{horizon}.csv'
    elif model_name == 'TCN':
        metrics_file = LSTM_DIR / 'tcn_metrics.csv'
        pred_file = LSTM_DIR / f'tcn_predictions_{horizon}.csv'

    try:
        # Load metrics
        if model_name == 'LSTM':
            metrics_df = pd.read_csv(metrics_file, index_col=0)
            if horizon not in metrics_df.index:
                return None
            metrics = metrics_df.loc[horizon].to_dict()
        elif model_name == 'TCN':
            metrics_df = pd.read_csv(metrics_file, index_col=0)
            if horizon not in metrics_df.index:
                return None
            metrics = metrics_df.loc[horizon].to_dict()
        else:
            metrics_df = pd.read_csv(metrics_file, index_col=0)
            if horizon not in metrics_df.index:
                return None
            metrics = metrics_df.loc[horizon].to_dict()
        
        # Load predictions
        pred_df = pd.read_csv(pred_file, index_col=0)
        
        return {
            'metrics': metrics,
            'predictions': pred_df
        }
    except Exception as e:
        print(f"Error loading {model_name} {horizon}: {e}")
        return None

def create_comparison_table():
    """Create comprehensive comparison table"""
    comparison_data = []

    for horizon in HORIZONS:
        horizon_data = {'Horizon': horizon}

        for model in MODELS:
            result = load_model_results(model, horizon)
            if result:
                for metric, value in result['metrics'].items():
                    horizon_data[f'{model}_{metric}'] = value
            else:
                # Use NaN for missing metrics so numeric operations work reliably
                for metric in ['RMSE', 'MAE', 'MAPE']:
                    horizon_data[f'{model}_{metric}'] = np.nan

        comparison_data.append(horizon_data)

    df = pd.DataFrame(comparison_data).set_index('Horizon')

    # Coerce all metric columns to numeric (non-numeric become NaN)
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    return df

def plot_metrics_comparison(comparison_df):
    """Plot metrics comparison across models and horizons"""
    metrics = ['RMSE', 'MAE', 'MAPE']
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    for i, metric in enumerate(metrics):
        ax = axes[i]

        # Prepare data for plotting
        plot_data = []
        for horizon in HORIZONS:
            for model in MODELS:
                col_name = f'{model}_{metric}'
                if col_name in comparison_df.columns:
                    value = comparison_df.loc[horizon, col_name]
                    # Ensure value is numeric and not NaN
                    try:
                        val = float(value)
                        if np.isfinite(val):
                            plot_data.append({
                                'Horizon': horizon.upper(),
                                'Model': model,
                                'Metric': metric,
                                'Value': val
                            })
                    except Exception:
                        # Skip non-numeric or missing values
                        continue

        plot_df = pd.DataFrame(plot_data)
        
        if not plot_df.empty:
            sns.barplot(data=plot_df, x='Horizon', y='Value', hue='Model', ax=ax)
            ax.set_title(f'{metric} Comparison')
            ax.set_ylabel(metric)
            ax.legend(title='Model')
            ax.grid(True, alpha=0.3)
            
            # Rotate x labels if needed
            ax.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(COMPARISON_DIR / 'metrics_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: metrics_comparison.png")

def plot_predictions_comparison():
    """Plot actual vs predicted for all models on each horizon"""
    
    for horizon in HORIZONS:
        fig, axes = plt.subplots(len(MODELS), 1, figsize=(14, 4*len(MODELS)))
        if len(MODELS) == 1:
            axes = [axes]
        
        for i, model in enumerate(MODELS):
            ax = axes[i]
            
            result = load_model_results(model, horizon)
            if result:
                pred_df = result['predictions']
                
                # Plot subset for visibility (first 200 points)
                max_points = min(200, len(pred_df))
                plot_df = pred_df.head(max_points)
                
                ax.plot(plot_df.index, plot_df['actual'], label='Actual', 
                       alpha=0.8, linewidth=1.5, color='black')
                ax.plot(plot_df.index, plot_df['forecast'] if 'forecast' in plot_df.columns 
                       else plot_df['prediction'], 
                       label=f'{model} Prediction', alpha=0.8, linewidth=1)
                
                ax.set_title(f'{model} - {horizon.upper()} Horizon')
                ax.set_xlabel('Time Steps')
                ax.set_ylabel('Demand (MW)')
                ax.legend()
                ax.grid(True, alpha=0.3)
            else:
                ax.text(0.5, 0.5, f'{model} results not available', 
                       transform=ax.transAxes, ha='center', va='center')
                ax.set_title(f'{model} - {horizon.upper()} Horizon (No Data)')
        
        plt.tight_layout()
        plt.savefig(COMPARISON_DIR / f'predictions_{horizon}.png', dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved: predictions_{horizon}.png")

def create_best_model_table(comparison_df):
    """Identify best performing model for each horizon and metric"""
    best_models = {}
    
    for horizon in HORIZONS:
        horizon_best = {}
        
        for metric in ['RMSE', 'MAE', 'MAPE']:
            # For RMSE and MAE, lower is better
            if metric in ['RMSE', 'MAE']:
                best_value = float('inf')
                best_model = None
                for model in MODELS:
                    col_name = f'{model}_{metric}'
                    if col_name in comparison_df.columns:
                        value = comparison_df.loc[horizon, col_name]
                        if not np.isnan(value) and value < best_value:
                            best_value = value
                            best_model = model
            # For MAPE, lower is better
            else:
                best_value = float('inf')
                best_model = None
                for model in MODELS:
                    col_name = f'{model}_{metric}'
                    if col_name in comparison_df.columns:
                        value = comparison_df.loc[horizon, col_name]
                        if not np.isnan(value) and value < best_value:
                            best_value = value
                            best_model = model
            
            horizon_best[metric] = {'model': best_model, 'value': best_value}
        
        best_models[horizon] = horizon_best
    
    return best_models

def plot_model_ranking(best_models):
    """Plot model ranking across horizons"""
    ranking_data = []
    
    for horizon, metrics in best_models.items():
        model_counts = {}
        for metric, info in metrics.items():
            model = info['model']
            if model:
                model_counts[model] = model_counts.get(model, 0) + 1
        
        for model, count in model_counts.items():
            ranking_data.append({
                'Horizon': horizon.upper(),
                'Model': model,
                'Wins': count
            })
    
    if ranking_data:
        ranking_df = pd.DataFrame(ranking_data)
        
        plt.figure(figsize=(10, 6))
        sns.barplot(data=ranking_df, x='Horizon', y='Wins', hue='Model')
        plt.title('Model Performance Ranking (Wins per Horizon)')
        plt.ylabel('Number of Best Performances')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(COMPARISON_DIR / 'model_ranking.png', dpi=150, bbox_inches='tight')
        plt.close()
        print("Saved: model_ranking.png")

# =============================================================================
# MAIN EXECUTION
# =============================================================================
print("=" * 70)
print("MODEL COMPARISON AND FINAL VISUALIZATION")
print("=" * 70)

# 1. Create comparison table
print("\n1. Creating comparison table...")
comparison_df = create_comparison_table()
print("Comparison table created")

# 2. Plot metrics comparison
print("\n2. Plotting metrics comparison...")
plot_metrics_comparison(comparison_df)

# 3. Plot predictions comparison
print("\n3. Plotting predictions comparison...")
plot_predictions_comparison()

# 4. Create best model analysis
print("\n4. Analyzing best performing models...")
best_models = create_best_model_table(comparison_df)
plot_model_ranking(best_models)

# 5. Save results
print("\n5. Saving results...")

# Save comparison table
comparison_df.to_csv(COMPARISON_DIR / 'model_comparison.csv')
print("Saved: model_comparison.csv")

# Save best models summary
with open(COMPARISON_DIR / 'best_models_summary.txt', 'w') as f:
    f.write("Best Performing Models Summary\n")
    f.write("=" * 40 + "\n\n")
    
    for horizon, metrics in best_models.items():
        f.write(f"{horizon.upper()} Horizon:\n")
        f.write("-" * 20 + "\n")
        for metric, info in metrics.items():
            if info['model']:
                f.write(f"  {metric}: {info['model']} ({info['value']:.4f})\n")
            else:
                f.write(f"  {metric}: No data\n")
        f.write("\n")

print("Saved: best_models_summary.txt")

# =============================================================================
# FINAL SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY - ALL MODELS")
print("=" * 70)

print("\nModel Comparison Table:")
print(comparison_df.to_string())

print("\nBest Models by Horizon and Metric:")
for horizon, metrics in best_models.items():
    print(f"\n{horizon.upper()}:")
    for metric, info in metrics.items():
        if info['model']:
            print(f"  {metric}: {info['model']} ({info['value']:.4f})")
        else:
            print(f"  {metric}: No data")

print(f"\nAll results saved in: {COMPARISON_DIR}")
print("\nFiles created:")
print("- model_comparison.csv (detailed metrics)")
print("- metrics_comparison.png (bar charts)")
print("- predictions_{horizon}.png (time series plots)")
print("- model_ranking.png (performance ranking)")
print("- best_models_summary.txt (text summary)")