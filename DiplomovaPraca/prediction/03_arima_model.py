"""
=============================================================================
ARIMA Model for Electricity Demand Forecasting
=============================================================================
Three horizons: Day (288), Week (2016), Month (8640)
=============================================================================
"""

import pandas as pd
import numpy as np
from pathlib import Path
import pickle
import warnings
warnings.filterwarnings('ignore')

from statsmodels.tsa.arima.model import ARIMA
from statsmodels.tsa.stattools import adfuller, acf, pacf
from sklearn.metrics import mean_squared_error, mean_absolute_error
import matplotlib.pyplot as plt

# =============================================================================
# CONFIGURATION
# =============================================================================
SPLITS_DIR = Path(__file__).parent / "data/splits"
RESULTS_DIR = Path(__file__).parent / "results/arima"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

TARGET_COL = 'demand_MW'

# ARIMA parameters (p, d, q)
# Will try auto-selection or use predefined
ARIMA_PARAMS = {
    'day': (2, 1, 2),
    'week': (2, 1, 2),
    'month': (2, 1, 2)
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
def calculate_metrics(y_true, y_pred):
    """Calculate RMSE, MAE, MAPE"""
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    mae = mean_absolute_error(y_true, y_pred)
    # MAPE (avoid division by zero)
    mask = y_true != 0
    mape = np.mean(np.abs((y_true[mask] - y_pred[mask]) / y_true[mask])) * 100
    return {'RMSE': rmse, 'MAE': mae, 'MAPE': mape}

def check_stationarity(series, name="Series"):
    """Augmented Dickey-Fuller test"""
    result = adfuller(series.dropna(), autolag='AIC')
    print(f"\n{name} - ADF Test:")
    print(f"  ADF Statistic: {result[0]:.4f}")
    print(f"  p-value: {result[1]:.4f}")
    print(f"  Stationary: {'Yes' if result[1] < 0.05 else 'No (need differencing)'}")
    return result[1] < 0.05

# =============================================================================
# 1. LOAD DATA
# =============================================================================
print("=" * 70)
print("1. LOAD DATA")
print("=" * 70)

with open(SPLITS_DIR / "train_test_splits.pkl", 'rb') as f:
    splits = pickle.load(f)

print(f"Loaded splits for horizons: {list(splits.keys())}")

# =============================================================================
# 2. STATIONARITY CHECK
# =============================================================================
print("\n" + "=" * 70)
print("2. STATIONARITY CHECK")
print("=" * 70)

# Check on day horizon train data
train_series = splits['day']['y_train']
is_stationary = check_stationarity(train_series, "Original Series")

if not is_stationary:
    diff_series = train_series.diff().dropna()
    check_stationarity(diff_series, "First Difference")

# =============================================================================
# 3. TRAIN ARIMA MODELS
# =============================================================================
print("\n" + "=" * 70)
print("3. TRAIN ARIMA MODELS")
print("=" * 70)

results = {}

for horizon_name, data in splits.items():
    print(f"\n{'='*50}")
    print(f"HORIZON: {horizon_name.upper()}")
    print(f"{'='*50}")
    
    y_train = data['y_train']
    y_test = data['y_test']
    
    print(f"Train size: {len(y_train):,}")
    print(f"Test size: {len(y_test):,}")
    
    # Get ARIMA parameters
    p, d, q = ARIMA_PARAMS[horizon_name]
    print(f"ARIMA parameters: p={p}, d={d}, q={q}")
    
    # For large datasets, use a subset for training (last N points)
    # ARIMA can be slow on very large datasets
    max_train_size = 10000  # Use last 10000 points for training
    if len(y_train) > max_train_size:
        y_train_subset = y_train.iloc[-max_train_size:]
        print(f"Using last {max_train_size} points for training")
    else:
        y_train_subset = y_train
    
    # Set frequency to avoid warnings
    y_train_subset.index.freq = '5min'
    
    try:
        print("Fitting ARIMA model...")
        model = ARIMA(y_train_subset, order=(p, d, q))
        model_fit = model.fit()
        
        print(f"AIC: {model_fit.aic:.2f}")
        print(f"BIC: {model_fit.bic:.2f}")
        
        # Forecast
        print("Generating forecast...")
        forecast = model_fit.forecast(steps=len(y_test))
        
        # Calculate metrics
        metrics = calculate_metrics(y_test.values, forecast.values)
        
        print(f"\nMetrics:")
        print(f"  RMSE: {metrics['RMSE']:.4f}")
        print(f"  MAE:  {metrics['MAE']:.4f}")
        print(f"  MAPE: {metrics['MAPE']:.2f}%")
        
        results[horizon_name] = {
            'model': model_fit,
            'forecast': forecast,
            'y_test': y_test,
            'y_train': y_train_subset,
            'metrics': metrics,
            'params': (p, d, q)
        }
        
    except Exception as e:
        print(f"Error fitting ARIMA for {horizon_name}: {e}")
        results[horizon_name] = None

# =============================================================================
# 4. VISUALIZATION
# =============================================================================
print("\n" + "=" * 70)
print("4. VISUALIZATION")
print("=" * 70)

for horizon_name, result in results.items():
    if result is None:
        continue
    
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    
    y_test = result['y_test']
    forecast = result['forecast']
    metrics = result['metrics']
    
    # Plot 1: Actual vs Predicted
    ax1 = axes[0]
    ax1.plot(y_test.index, y_test.values, label='Actual', alpha=0.8, linewidth=0.8)
    ax1.plot(y_test.index, forecast.values, label='ARIMA Forecast', alpha=0.8, linewidth=0.8)
    ax1.set_title(f'ARIMA Forecast - {horizon_name.upper()} Horizon\n'
                  f'RMSE: {metrics["RMSE"]:.2f}, MAE: {metrics["MAE"]:.2f}, MAPE: {metrics["MAPE"]:.2f}%')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Demand (MW)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Residuals
    ax2 = axes[1]
    residuals = y_test.values - forecast.values
    ax2.plot(y_test.index, residuals, alpha=0.7, linewidth=0.5)
    ax2.axhline(y=0, color='r', linestyle='--', alpha=0.5)
    ax2.set_title('Forecast Residuals')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Residual (MW)')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / f'arima_{horizon_name}.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: arima_{horizon_name}.png")

# =============================================================================
# 5. SAVE RESULTS
# =============================================================================
print("\n" + "=" * 70)
print("5. SAVE RESULTS")
print("=" * 70)

# Save metrics summary
metrics_df = pd.DataFrame({
    horizon: result['metrics'] for horizon, result in results.items() if result is not None
}).T
metrics_df.index.name = 'Horizon'
metrics_df.to_csv(RESULTS_DIR / 'arima_metrics.csv')
print(f"Metrics saved: arima_metrics.csv")

# Save forecasts
for horizon_name, result in results.items():
    if result is None:
        continue
    
    forecast_df = pd.DataFrame({
        'actual': result['y_test'],
        'forecast': result['forecast']
    })
    forecast_df.to_csv(RESULTS_DIR / f'arima_forecast_{horizon_name}.csv')
    print(f"Saved: arima_forecast_{horizon_name}.csv")

# Save model summary
with open(RESULTS_DIR / 'arima_summary.txt', 'w') as f:
    f.write("ARIMA Model Results\n")
    f.write("=" * 50 + "\n\n")
    for horizon_name, result in results.items():
        if result is None:
            continue
        f.write(f"\n{horizon_name.upper()} Horizon:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Parameters: ARIMA{result['params']}\n")
        f.write(f"RMSE: {result['metrics']['RMSE']:.4f}\n")
        f.write(f"MAE: {result['metrics']['MAE']:.4f}\n")
        f.write(f"MAPE: {result['metrics']['MAPE']:.2f}%\n")
print(f"Summary saved: arima_summary.txt")

# =============================================================================
# FINAL SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY - ARIMA")
print("=" * 70)

print("\n" + metrics_df.to_string())
