"""
=============================================================================
LightGBM Baseline for Electricity Demand Forecasting
=============================================================================
Gradient boosting baseline with feature engineering and hyperparameter tuning
=============================================================================
"""

import pandas as pd
import numpy as np
from pathlib import Path
import pickle
import warnings
warnings.filterwarnings('ignore')

try:
    import lightgbm as lgb
    from lightgbm import early_stopping, log_evaluation
except ImportError:
    print("LightGBM not installed. Install with: pip install lightgbm")
    raise

from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import joblib

# =============================================================================
# CONFIGURATION
# =============================================================================
SPLITS_DIR = Path(__file__).parent / "data/splits"
RESULTS_DIR = Path(__file__).parent / "results/lightgbm"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

TARGET_COL = 'demand_MW'

# LightGBM parameters
LGBM_PARAMS = {
    'objective': 'regression',
    'metric': 'rmse',
    'boosting_type': 'gbdt',
    'learning_rate': 0.05,
    'num_leaves': 31,
    'max_depth': -1,
    'min_data_in_leaf': 20,
    'feature_fraction': 0.8,
    'bagging_fraction': 0.8,
    'bagging_freq': 5,
    'lambda_l1': 0.1,
    'lambda_l2': 0.1,
    'verbosity': -1,
    'seed': 42
}

# Training parameters
NUM_BOOST_ROUND = 1000
EARLY_STOPPING_ROUNDS = 50
VALIDATION_SIZE = 0.1

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

def plot_feature_importance(model, cleaned_names, original_names, horizon, save_path):
    """Plot and save feature importance"""
    try:
        importance = model.feature_importance(importance_type='gain')
        indices = np.argsort(importance)[::-1][:20]  # Top 20 features
        
        plt.figure(figsize=(12, 8))
        plt.barh(range(len(indices)), importance[indices])
        plt.yticks(range(len(indices)), [original_names[i] for i in indices])
        plt.xlabel('Feature Importance (Gain)')
        plt.title(f'LightGBM Feature Importance - {horizon.upper()} Horizon')
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Feature importance saved: {save_path}")
    except Exception as e:
        print(f"Could not plot feature importance: {e}")

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
# 2. TRAIN LIGHTGBM MODELS
# =============================================================================
print("\n" + "=" * 70)
print("2. TRAIN LIGHTGBM MODELS")
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
    
    # Clean feature names for LightGBM (remove special JSON characters)
    def clean_feature_names(columns):
        """Clean column names to be compatible with LightGBM"""
        import re
        cleaned = []
        for col in columns:
            # Replace special characters with underscores
            cleaned_col = re.sub(r'[^\w]', '_', str(col))
            # Ensure no leading/trailing underscores and no double underscores
            cleaned_col = re.sub(r'_+', '_', cleaned_col).strip('_')
            cleaned.append(cleaned_col)
        return cleaned

    original_feature_names = list(X_train.columns)
    cleaned_feature_names = clean_feature_names(original_feature_names)

    # Rename columns on training set first, then create validation split (preserve temporal order)
    X_train_renamed = X_train.rename(columns=dict(zip(original_feature_names, cleaned_feature_names)))

    # Create validation split from training data (temporal order preserved)
    X_tr, X_val, y_tr, y_val = train_test_split(
        X_train_renamed, y_train, 
        test_size=VALIDATION_SIZE, 
        shuffle=False,  # Preserve temporal order
        random_state=42
    )

    # Rename test and validation columns (they share the same mapping)
    X_test_renamed = X_test.rename(columns=dict(zip(original_feature_names, cleaned_feature_names)))
    X_tr_renamed = X_tr  # already renamed via X_train_renamed
    X_val_renamed = X_val

    print(f"Original features: {original_feature_names[:5]}...")
    print(f"Cleaned features: {cleaned_feature_names[:5]}...")
    
    print(f"Train subset: {len(X_tr):,}")
    print(f"Validation subset: {len(X_val):,}")
    
    # Create LightGBM datasets
    lgb_train = lgb.Dataset(X_tr, label=y_tr, feature_name=cleaned_feature_names)
    lgb_val = lgb.Dataset(X_val, label=y_val, reference=lgb_train)
    
    print("Training LightGBM model...")
    try:
        # Train with early stopping
        model = lgb.train(
            LGBM_PARAMS,
            lgb_train,
            num_boost_round=NUM_BOOST_ROUND,
            valid_sets=[lgb_train, lgb_val],
            callbacks=[
                early_stopping(EARLY_STOPPING_ROUNDS),
                log_evaluation(0)  # Silent
            ]
        )
        
        print(f"Best iteration: {model.best_iteration}")
        
        # Make predictions
        print("Generating predictions...")
        y_pred = model.predict(X_test_renamed, num_iteration=model.best_iteration)
        
        # Calculate metrics
        metrics = calculate_metrics(y_test.values, y_pred)
        
        print(f"\nMetrics:")
        print(f"  RMSE: {metrics['RMSE']:.4f}")
        print(f"  MAE:  {metrics['MAE']:.4f}")
        print(f"  MAPE: {metrics['MAPE']:.2f}%")
        
        # Plot feature importance
        importance_path = RESULTS_DIR / f'lightgbm_importance_{horizon_name}.png'
        plot_feature_importance(model, cleaned_feature_names, original_feature_names, horizon_name, importance_path)
        
        # Save model and results
        model_path = RESULTS_DIR / f'lightgbm_{horizon_name}.txt'
        model.save_model(str(model_path))
        
        pred_df = pd.DataFrame({
            'actual': y_test,
            'prediction': y_pred
        })
        pred_df.to_csv(RESULTS_DIR / f'lightgbm_predictions_{horizon_name}.csv')
        
        print(f"Saved: lightgbm_{horizon_name}.txt, predictions_{horizon_name}.csv")
        
        results[horizon_name] = {
            'model': model,
            'metrics': metrics,
            'best_iteration': model.best_iteration,
            'feature_names': original_feature_names
        }
        
    except Exception as e:
        print(f"Error training LightGBM for {horizon_name}: {e}")
        results[horizon_name] = None

# =============================================================================
# 3. SAVE SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("3. SAVE SUMMARY")
print("=" * 70)

# Save metrics summary
metrics_df = pd.DataFrame({
    horizon: result['metrics'] for horizon, result in results.items() if result is not None
}).T
metrics_df.index.name = 'Horizon'
metrics_df.to_csv(RESULTS_DIR / 'lightgbm_metrics.csv')
print(f"Metrics saved: lightgbm_metrics.csv")

# Save detailed summary
with open(RESULTS_DIR / 'lightgbm_summary.txt', 'w') as f:
    f.write("LightGBM Model Results\n")
    f.write("=" * 50 + "\n\n")
    for horizon_name, result in results.items():
        if result is None:
            continue
        f.write(f"\n{horizon_name.upper()} Horizon:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Best iteration: {result['best_iteration']}\n")
        f.write(f"RMSE: {result['metrics']['RMSE']:.4f}\n")
        f.write(f"MAE: {result['metrics']['MAE']:.4f}\n")
        f.write(f"MAPE: {result['metrics']['MAPE']:.2f}%\n")
        f.write(f"Features: {len(result['feature_names'])}\n")
print(f"Summary saved: lightgbm_summary.txt")

# =============================================================================
# FINAL SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY - LIGHTGBM")
print("=" * 70)

print("\n" + metrics_df.to_string())
