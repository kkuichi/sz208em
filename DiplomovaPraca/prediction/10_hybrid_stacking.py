"""
Hybrid and simple ensembling of existing model predictions.
Creates average, median, and linear stacking ensembles and evaluates them for each horizon.
Includes visualizations and comprehensive error handling.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import warnings
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error, mean_absolute_error
import logging

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

RESULTS_DIR = Path(__file__).parent / 'results'
COMPARISON_DIR = RESULTS_DIR / 'comparison'
COMPARISON_DIR.mkdir(parents=True, exist_ok=True)

MODEL_DIRS = {
    'ARIMA': RESULTS_DIR / 'arima',
    'SVM': RESULTS_DIR / 'svm',
    'CNN': RESULTS_DIR / 'cnn',
    'LightGBM': RESULTS_DIR / 'lightgbm',
    'LSTM': RESULTS_DIR / 'lstm_tcn',
    'Advanced': RESULTS_DIR / 'advanced'  # For N-BEATS and TFT if available
}

HORIZONS = ['day', 'week', 'month']

def calc_metrics(y_true, y_pred):
    """Calculate RMSE, MAE, and MAPE metrics."""
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    mae = mean_absolute_error(y_true, y_pred)
    mask = y_true != 0
    if np.sum(mask) == 0:
        mape = np.nan
    else:
        mape = np.mean(np.abs((y_true[mask] - y_pred[mask]) / y_true[mask])) * 100
    return {'RMSE': rmse, 'MAE': mae, 'MAPE': mape}

def load_predictions(model_dir, model_name, horizon):
    """Load predictions for a specific model and horizon."""
    possible_files = [
        model_dir / f'{model_name.lower()}_predictions_{horizon}.csv',
        model_dir / f'{model_name.lower()}_prediction_{horizon}.csv',
        model_dir / f'{model_name.lower()}_forecast_{horizon}.csv',
        model_dir / f'{model_name.lower()}_{horizon}_predictions.csv',
        model_dir / f'{model_name.lower()}_{horizon}_forecast.csv'
    ]
    
    for file_path in possible_files:
        if file_path.exists():
            try:
                df = pd.read_csv(file_path, index_col=0)
                # Standardize column names
                if 'forecast' in df.columns:
                    pred_col = 'forecast'
                elif 'prediction' in df.columns:
                    pred_col = 'prediction'
                elif 'predictions' in df.columns:
                    pred_col = 'predictions'
                else:
                    pred_col = df.columns[-1]  # Assume last column is predictions
                
                predictions = df[pred_col].values
                actual = df['actual'].values if 'actual' in df.columns else None
                return predictions, actual
            except Exception as e:
                logger.warning(f"Error loading {file_path}: {e}")
                continue
    
    return None, None

def create_visualizations(actual, predictions_dict, horizon, save_dir):
    """Create and save visualization plots for ensemble predictions."""
    plt.figure(figsize=(15, 10))
    
    # Plot all predictions
    plt.subplot(2, 1, 1)
    plt.plot(actual, label='Actual', color='black', linewidth=2)
    colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown']
    for i, (name, pred) in enumerate(predictions_dict.items()):
        color = colors[i % len(colors)]
        plt.plot(pred, label=name, color=color, alpha=0.7)
    plt.title(f'Ensemble Predictions vs Actual - {horizon.capitalize()} Horizon')
    plt.xlabel('Time Steps')
    plt.ylabel('Demand')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot residuals
    plt.subplot(2, 1, 2)
    for i, (name, pred) in enumerate(predictions_dict.items()):
        color = colors[i % len(colors)]
        residuals = actual - pred
        plt.plot(residuals, label=f'{name} Residuals', color=color, alpha=0.7)
    plt.title(f'Prediction Residuals - {horizon.capitalize()} Horizon')
    plt.xlabel('Time Steps')
    plt.ylabel('Residuals')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(save_dir / f'ensemble_visualization_{horizon}.png', dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved visualization for {horizon} horizon")

ensemble_results = {}

for h in HORIZONS:
    logger.info(f"Processing horizon: {h}")
    preds = {}
    actual = None
    
    for model_name, model_dir in MODEL_DIRS.items():
        predictions, act = load_predictions(model_dir, model_name, h)
        if predictions is not None:
            preds[model_name] = predictions
            if actual is None and act is not None:
                actual = act
            elif actual is not None and act is not None:
                # Verify actual values are consistent (aligning from the end)
                min_a_len = min(len(actual), len(act))
                if not np.allclose(actual[-min_a_len:], act[-min_a_len:], rtol=1e-10):
                    logger.warning(f"Actual values differ for {model_name} in {h} horizon (checked last {min_a_len} steps)")
        else:
            logger.info(f"No predictions found for {model_name} in {h} horizon")

    if actual is None or len(preds) == 0:
        logger.warning(f'No valid predictions found for horizon {h}. Skipping ensemble.')
        continue

    # Align lengths (truncate to shortest, matching from the end since sequence models drop initial steps)
    min_len = min(len(preds[model]) for model in preds)
    actual = actual[-min_len:]
    for model in list(preds.keys()):
        preds[model] = preds[model][-min_len:]

    logger.info(f"Using {len(preds)} models for {h} horizon: {list(preds.keys())}")

    # Create ensembles
    stacked = np.vstack([preds[model] for model in preds]).T
    
    # Simple average ensemble
    avg_pred = stacked.mean(axis=1)
    avg_metrics = calc_metrics(actual, avg_pred)
    
    # Median ensemble
    median_pred = np.median(stacked, axis=1)
    median_metrics = calc_metrics(actual, median_pred)
    
    # Linear stacking (Ridge) - Note: Using same data for simplicity, but ideally should use cross-validation
    ridge_metrics = None
    ridge_pred = None
    try:
        # Split data for meta-model training (use first 70% for training meta-model)
        split_idx = int(0.7 * len(stacked))
        X_train = stacked[:split_idx]
        y_train = actual[:split_idx]
        X_test = stacked[split_idx:]
        y_test = actual[split_idx:]
        
        ridge = Ridge(alpha=1.0)
        ridge.fit(X_train, y_train)
        ridge_pred = ridge.predict(X_test)
        ridge_metrics = calc_metrics(y_test, ridge_pred)
        
        # For full prediction, retrain on all data
        ridge.fit(stacked, actual)
        ridge_pred_full = ridge.predict(stacked)
    except Exception as e:
        logger.error(f"Error in Ridge stacking for {h}: {e}")
        ridge_pred_full = None

    ensemble_results[h] = {
        'models_used': list(preds.keys()),
        'avg_metrics': avg_metrics,
        'median_metrics': median_metrics,
        'ridge_metrics': ridge_metrics
    }

    # Prepare predictions for visualization and saving
    predictions_dict = {
        'Average Ensemble': avg_pred,
        'Median Ensemble': median_pred
    }
    if ridge_pred_full is not None:
        predictions_dict['Ridge Stacking'] = ridge_pred_full

    # Create visualizations
    create_visualizations(actual, predictions_dict, h, COMPARISON_DIR)

    # Save detailed predictions
    out_df = pd.DataFrame({'actual': actual})
    for name, pred in predictions_dict.items():
        out_df[name.replace(' ', '_').lower()] = pred
    out_df.to_csv(COMPARISON_DIR / f'ensemble_predictions_{h}.csv')
    logger.info(f"Saved predictions for {h} horizon")

# Save comprehensive summary
summary_data = []
for h, r in ensemble_results.items():
    base_row = {
        'horizon': h,
        'models_used': ','.join(r['models_used']),
        'num_models': len(r['models_used'])
    }
    
    # Average ensemble metrics
    base_row.update({f'avg_{k}': v for k, v in r['avg_metrics'].items()})
    
    # Median ensemble metrics
    base_row.update({f'median_{k}': v for k, v in r['median_metrics'].items()})
    
    # Ridge metrics (if available)
    if r['ridge_metrics']:
        base_row.update({f'ridge_{k}': v for k, v in r['ridge_metrics'].items()})
    else:
        base_row.update({f'ridge_{k}': np.nan for k in ['RMSE', 'MAE', 'MAPE']})
    
    summary_data.append(base_row)

summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(COMPARISON_DIR / 'ensemble_summary.csv', index=False)

# Print summary
print("\n" + "="*80)
print("ENSEMBLE RESULTS SUMMARY")
print("="*80)
for _, row in summary_df.iterrows():
    print(f"\nHorizon: {row['horizon'].upper()}")
    print(f"Models used: {row['models_used']}")
    print(f"Average Ensemble - RMSE: {row['avg_RMSE']:.4f}, MAE: {row['avg_MAE']:.4f}, MAPE: {row['avg_MAPE']:.4f}")
    print(f"Median Ensemble  - RMSE: {row['median_RMSE']:.4f}, MAE: {row['median_MAE']:.4f}, MAPE: {row['median_MAPE']:.4f}")
    if not np.isnan(row['ridge_RMSE']):
        print(f"Ridge Stacking    - RMSE: {row['ridge_RMSE']:.4f}, MAE: {row['ridge_MAE']:.4f}, MAPE: {row['ridge_MAPE']:.4f}")
    else:
        print("Ridge Stacking    - Not available")

print(f"\nDetailed results saved in: {COMPARISON_DIR}")
print("Ensembling completed successfully!")
logger.info("Ensembling process completed")
