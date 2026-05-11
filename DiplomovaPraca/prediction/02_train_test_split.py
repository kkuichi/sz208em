"""
=============================================================================
Train/Test Split for Three Prediction Horizons
=============================================================================
Horizons:
- Day:   288 points (24h × 12 intervals per hour at 5-min)
- Week:  2016 points (7 × 288)
- Month: 8640 points (30 × 288)
=============================================================================
"""

import pandas as pd
import numpy as np
from pathlib import Path
import pickle

# =============================================================================
# CONFIGURATION
# =============================================================================
DATA_DIR = Path(__file__).parent / "data/prepared"
SPLITS_DIR = Path(__file__).parent / "data/splits"
SPLITS_DIR.mkdir(parents=True, exist_ok=True)
PREPARED_DATA = DATA_DIR / "prepared_data.csv"


# Prediction horizons (number of 5-minute intervals)
HORIZONS = {
    'day': 288,      # 24 hours
    'week': 2016,    # 7 days
    'month': 8640    # 30 days
}

TARGET_COL = 'demand_MW'

# =============================================================================
# 1. LOAD DATA
# =============================================================================
print("=" * 70)
print("1. LOAD DATA")
print("=" * 70)

df = pd.read_csv(PREPARED_DATA, index_col=0, parse_dates=True)
print(f"Dataset size: {df.shape[0]:,} rows × {df.shape[1]} columns")
print(f"Period: {df.index.min()} — {df.index.max()}")
print(f"Target column: {TARGET_COL}")

# Verify target column exists
if TARGET_COL not in df.columns:
    # Try to find it
    for col in df.columns:
        if 'demand' in col.lower():
            TARGET_COL = col
            print(f"Found target column: {TARGET_COL}")
            break

# =============================================================================
# 2. CREATE TRAIN/TEST SPLITS
# =============================================================================
print("\n" + "=" * 70)
print("2. CREATE TRAIN/TEST SPLITS")
print("=" * 70)

splits = {}

for horizon_name, horizon_size in HORIZONS.items():
    print(f"\n--- {horizon_name.upper()} ({horizon_size} points) ---")
    
    # Split point
    split_idx = len(df) - horizon_size
    
    if split_idx < horizon_size:
        print(f"  WARNING: Not enough data for {horizon_name} horizon!")
        continue
    
    train = df.iloc[:split_idx].copy()
    test = df.iloc[split_idx:].copy()
    
    print(f"  Train: {len(train):,} rows ({train.index.min()} — {train.index.max()})")
    print(f"  Test:  {len(test):,} rows ({test.index.min()} — {test.index.max()})")
    
    splits[horizon_name] = {
        'train': train,
        'test': test,
        'train_target': train[TARGET_COL],
        'test_target': test[TARGET_COL],
        'horizon_size': horizon_size
    }

# =============================================================================
# 3. FEATURE ENGINEERING
# =============================================================================
print("\n" + "=" * 70)
print("3. FEATURE ENGINEERING")
print("=" * 70)

def create_features(df, target_col, lags=[1, 2, 3, 6, 12, 24, 48, 288]):
    """Create lag features and rolling statistics"""
    df_feat = df.copy()
    
    # Lag features
    for lag in lags:
        df_feat[f'lag_{lag}'] = df_feat[target_col].shift(lag)
    
    # Rolling statistics
    for window in [12, 24, 48, 288]:  # 1h, 2h, 4h, 1day
        df_feat[f'rolling_mean_{window}'] = df_feat[target_col].shift(1).rolling(window).mean()
        df_feat[f'rolling_std_{window}'] = df_feat[target_col].shift(1).rolling(window).std()
    
    # Time features (if not already present)
    if 'hour' not in df_feat.columns:
        df_feat['hour'] = df_feat.index.hour
    if 'dayofweek' not in df_feat.columns:
        df_feat['dayofweek'] = df_feat.index.dayofweek
    if 'month' not in df_feat.columns:
        df_feat['month'] = df_feat.index.month
    if 'is_weekend' not in df_feat.columns:
        df_feat['is_weekend'] = df_feat['dayofweek'].isin([5, 6]).astype(int)
    
    # Cyclical encoding for hour
    df_feat['hour_sin'] = np.sin(2 * np.pi * df_feat['hour'] / 24)
    df_feat['hour_cos'] = np.cos(2 * np.pi * df_feat['hour'] / 24)
    
    # Cyclical encoding for day of week
    df_feat['dow_sin'] = np.sin(2 * np.pi * df_feat['dayofweek'] / 7)
    df_feat['dow_cos'] = np.cos(2 * np.pi * df_feat['dayofweek'] / 7)
    
    return df_feat

# Apply feature engineering
print("Creating features...")
df_features = create_features(df, TARGET_COL)

# Drop rows with NaN (from lag/rolling features)
initial_len = len(df_features)
df_features = df_features.dropna()
print(f"Dropped {initial_len - len(df_features)} rows with NaN values")
print(f"Final dataset size: {len(df_features):,} rows")

# Feature columns (excluding target and original time features)
exclude_cols = [TARGET_COL, 'hour', 'dayofweek', 'month']
feature_cols = [col for col in df_features.columns if col not in exclude_cols]
print(f"\nFeature columns ({len(feature_cols)}):")
print(feature_cols)

# =============================================================================
# 4. CREATE FINAL SPLITS WITH FEATURES
# =============================================================================
print("\n" + "=" * 70)
print("4. CREATE FINAL SPLITS WITH FEATURES")
print("=" * 70)

final_splits = {}

for horizon_name, horizon_size in HORIZONS.items():
    print(f"\n--- {horizon_name.upper()} ---")
    
    split_idx = len(df_features) - horizon_size
    
    if split_idx < horizon_size:
        print(f"  WARNING: Not enough data!")
        continue
    
    train_df = df_features.iloc[:split_idx]
    test_df = df_features.iloc[split_idx:]
    
    X_train = train_df[feature_cols]
    y_train = train_df[TARGET_COL]
    X_test = test_df[feature_cols]
    y_test = test_df[TARGET_COL]
    
    print(f"  X_train: {X_train.shape}")
    print(f"  y_train: {y_train.shape}")
    print(f"  X_test:  {X_test.shape}")
    print(f"  y_test:  {y_test.shape}")
    
    final_splits[horizon_name] = {
        'X_train': X_train,
        'y_train': y_train,
        'X_test': X_test,
        'y_test': y_test,
        'train_df': train_df,
        'test_df': test_df,
        'feature_cols': feature_cols,
        'horizon_size': horizon_size
    }

# =============================================================================
# 5. SAVE SPLITS
# =============================================================================
print("\n" + "=" * 70)
print("5. SAVE SPLITS")
print("=" * 70)

# Save as pickle for easy loading
output_file = SPLITS_DIR / "train_test_splits.pkl"
with open(output_file, 'wb') as f:
    pickle.dump(final_splits, f)
print(f"Splits saved: {output_file}")

# Also save individual CSVs for reference
for horizon_name, data in final_splits.items():
    data['train_df'].to_csv(SPLITS_DIR / f"train_{horizon_name}.csv")
    data['test_df'].to_csv(SPLITS_DIR / f"test_{horizon_name}.csv")
    print(f"Saved: train_{horizon_name}.csv, test_{horizon_name}.csv")

# Save feature list
with open(SPLITS_DIR / "feature_cols.pkl", 'wb') as f:
    pickle.dump(feature_cols, f)
print(f"Feature columns saved: {SPLITS_DIR / 'feature_cols.pkl'}")

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

for horizon_name, data in final_splits.items():
    print(f"\n{horizon_name.upper()}:")
    print(f"  Train samples: {len(data['y_train']):,}")
    print(f"  Test samples:  {len(data['y_test']):,}")
    print(f"  Features:      {len(data['feature_cols'])}")
