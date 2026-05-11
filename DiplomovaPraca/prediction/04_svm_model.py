"""
=============================================================================
SVM Model for Electricity Demand Forecasting
=============================================================================
Support Vector Regression with RBF kernel
Three horizons: Day (288), Week (2016), Month (8640)
=============================================================================
"""

import pandas as pd
import numpy as np
from pathlib import Path
import pickle
import warnings
warnings.filterwarnings('ignore')

from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt
import joblib

# =============================================================================
# CONFIGURATION
# =============================================================================
SPLITS_DIR = Path(__file__).parent / "data/splits"
RESULTS_DIR = Path(__file__).parent / "results/svm"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

TARGET_COL = 'demand_MW'

# SVM parameters
SVM_PARAMS = {
    'C': [0.1, 1, 10, 100],
    'gamma': ['scale', 'auto', 0.001, 0.01, 0.1],
    'epsilon': [0.01, 0.1, 0.2]
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

def train_svm_with_grid_search(X_train, y_train, cv_folds=3):
    """Train SVM with grid search for hyperparameter tuning"""
    print("Performing grid search for SVM parameters...")
    
    # Use subset for grid search to speed up
    max_samples = min(5000, len(X_train))
    X_subset = X_train.iloc[-max_samples:] if len(X_train) > max_samples else X_train
    y_subset = y_train.iloc[-max_samples:] if len(y_train) > max_samples else y_train
    
    print(f"Using {len(X_subset)} samples for grid search")
    
    svm = SVR(kernel='rbf')
    grid_search = GridSearchCV(
        svm, 
        SVM_PARAMS, 
        cv=cv_folds, 
        scoring='neg_mean_squared_error',
        n_jobs=-1,
        verbose=1
    )
    
    grid_search.fit(X_subset, y_subset)
    
    print(f"Best parameters: {grid_search.best_params_}")
    print(f"Best CV score (neg MSE): {grid_search.best_score_:.4f}")
    
    return grid_search.best_estimator_

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
# 2. TRAIN SVM MODELS
# =============================================================================
print("\n" + "=" * 70)
print("2. TRAIN SVM MODELS")
print("=" * 70)

results = {}

for horizon_name, data in splits.items():
    print(f"\n{'='*50}")
    print(f"HORIZON: {horizon_name.upper()}")
    print(f"{'='*50}")
    
    X_train = data['X_train']
    y_train = data['y_train']
    X_test = data['X_test']
    y_test = data['y_test']
    
    print(f"Train size: {len(X_train):,}")
    print(f"Test size: {len(X_test):,}")
    print(f"Features: {len(X_train.columns)}")
    
    # Scale features
    print("Scaling features...")
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # For large datasets, use subset for training
    max_train_size = 10000
    if len(X_train) > max_train_size:
        indices = np.random.choice(len(X_train), max_train_size, replace=False)
        X_train_scaled = X_train_scaled[indices]
        y_train = y_train.iloc[indices]
        print(f"Using random subset of {max_train_size} samples for training")
    
    try:
        # Train SVM with grid search
        svm_model = train_svm_with_grid_search(
            pd.DataFrame(X_train_scaled, columns=X_train.columns), 
            y_train
        )
        
        # Fit on full training data
        print("Fitting SVM on full training data...")
        svm_model.fit(X_train_scaled, y_train)
        
        # Predict
        print("Generating predictions...")
        y_pred = svm_model.predict(X_test_scaled)
        
        # Calculate metrics
        metrics = calculate_metrics(y_test.values, y_pred)
        
        print(f"\nMetrics:")
        print(f"  RMSE: {metrics['RMSE']:.4f}")
        print(f"  MAE:  {metrics['MAE']:.4f}")
        print(f"  MAPE: {metrics['MAPE']:.2f}%")
        
        results[horizon_name] = {
            'model': svm_model,
            'scaler': scaler,
            'predictions': y_pred,
            'y_test': y_test,
            'y_train': y_train,
            'X_train': X_train,
            'X_test': X_test,
            'metrics': metrics,
            'best_params': svm_model.get_params()
        }
        
    except Exception as e:
        print(f"Error training SVM for {horizon_name}: {e}")
        results[horizon_name] = None

# =============================================================================
# 3. VISUALIZATION
# =============================================================================
print("\n" + "=" * 70)
print("3. VISUALIZATION")
print("=" * 70)

for horizon_name, result in results.items():
    if result is None:
        continue
    
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    
    y_test = result['y_test']
    y_pred = result['predictions']
    metrics = result['metrics']
    
    # Plot 1: Actual vs Predicted
    ax1 = axes[0]
    ax1.plot(y_test.index, y_test.values, label='Actual', alpha=0.8, linewidth=0.8)
    ax1.plot(y_test.index, y_pred, label='SVM Prediction', alpha=0.8, linewidth=0.8)
    ax1.set_title(f'SVM Forecast - {horizon_name.upper()} Horizon\n'
                  f'RMSE: {metrics["RMSE"]:.2f}, MAE: {metrics["MAE"]:.2f}, MAPE: {metrics["MAPE"]:.2f}%')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Demand (MW)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Residuals
    ax2 = axes[1]
    residuals = y_test.values - y_pred
    ax2.plot(y_test.index, residuals, alpha=0.7, linewidth=0.5)
    ax2.axhline(y=0, color='r', linestyle='--', alpha=0.5)
    ax2.set_title('Prediction Residuals')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Residual (MW)')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / f'svm_{horizon_name}.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: svm_{horizon_name}.png")

# =============================================================================
# 4. SAVE RESULTS
# =============================================================================
print("\n" + "=" * 70)
print("4. SAVE RESULTS")
print("=" * 70)

# Save metrics summary
metrics_df = pd.DataFrame({
    horizon: result['metrics'] for horizon, result in results.items() if result is not None
}).T
metrics_df.index.name = 'Horizon'
metrics_df.to_csv(RESULTS_DIR / 'svm_metrics.csv')
print(f"Metrics saved: svm_metrics.csv")

# Save predictions
for horizon_name, result in results.items():
    if result is None:
        continue
    
    pred_df = pd.DataFrame({
        'actual': result['y_test'],
        'prediction': result['predictions']
    })
    pred_df.to_csv(RESULTS_DIR / f'svm_predictions_{horizon_name}.csv')
    print(f"Saved: svm_predictions_{horizon_name}.csv")

# Save models
for horizon_name, result in results.items():
    if result is None:
        continue
    
    model_path = RESULTS_DIR / f'svm_model_{horizon_name}.pkl'
    joblib.dump(result['model'], model_path)
    
    scaler_path = RESULTS_DIR / f'svm_scaler_{horizon_name}.pkl'
    joblib.dump(result['scaler'], scaler_path)
    
    print(f"Saved: svm_model_{horizon_name}.pkl, svm_scaler_{horizon_name}.pkl")

# Save model summary
with open(RESULTS_DIR / 'svm_summary.txt', 'w') as f:
    f.write("SVM Model Results\n")
    f.write("=" * 50 + "\n\n")
    for horizon_name, result in results.items():
        if result is None:
            continue
        f.write(f"\n{horizon_name.upper()} Horizon:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Best parameters: {result['best_params']}\n")
        f.write(f"RMSE: {result['metrics']['RMSE']:.4f}\n")
        f.write(f"MAE: {result['metrics']['MAE']:.4f}\n")
        f.write(f"MAPE: {result['metrics']['MAPE']:.2f}%\n")
print(f"Summary saved: svm_summary.txt")

# =============================================================================
# FINAL SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY - SVM")
print("=" * 70)

print("\n" + metrics_df.to_string())