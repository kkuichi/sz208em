"""
=============================================================================
CNN Model for Electricity Demand Forecasting
=============================================================================
Convolutional Neural Network with Conv1D layers
Three horizons: Day (288), Week (2016), Month (8640)
=============================================================================
"""

import pandas as pd
import numpy as np
from pathlib import Path
import pickle
import warnings
warnings.filterwarnings('ignore')

import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Flatten, Dense, Dropout, BatchNormalization
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import EarlyStopping, ReduceLROnPlateau
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, mean_absolute_error
import matplotlib.pyplot as plt
import joblib

# Set random seeds for reproducibility
tf.random.set_seed(42)
np.random.seed(42)

# Suppress TensorFlow warnings
import os
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
tf.get_logger().setLevel('ERROR')
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

# =============================================================================
# CONFIGURATION
# =============================================================================
SPLITS_DIR = Path(__file__).parent / "data/splits"
RESULTS_DIR = Path(__file__).parent / "results/cnn"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

TARGET_COL = 'demand_MW'

# CNN Architecture parameters
CNN_PARAMS = {
    'filters': [32, 64, 128],
    'kernel_size': [3, 5, 7],
    'pool_size': [2, 3],
    'dense_units': [64, 128],
    'dropout': [0.1, 0.2, 0.3]
}

# Training parameters
EPOCHS = 100
BATCH_SIZE = 32
PATIENCE = 10

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

def create_sequences(X, y, seq_length):
    """Create sequences for CNN input"""
    X_seq, y_seq = [], []
    for i in range(len(X) - seq_length + 1):
        X_seq.append(X[i:i+seq_length])
        y_seq.append(y[i+seq_length-1])
    return np.array(X_seq), np.array(y_seq)

def build_cnn_model(input_shape, params):
    """Build CNN model with Conv1D layers"""
    model = Sequential()
    
    # First Conv1D layer
    model.add(Conv1D(
        filters=params['filters'], 
        kernel_size=params['kernel_size'], 
        activation='relu', 
        input_shape=input_shape
    ))
    model.add(BatchNormalization())
    model.add(MaxPooling1D(pool_size=params['pool_size']))
    model.add(Dropout(params['dropout']))
    
    # Second Conv1D layer
    model.add(Conv1D(
        filters=params['filters']*2, 
        kernel_size=params['kernel_size'], 
        activation='relu'
    ))
    model.add(BatchNormalization())
    model.add(MaxPooling1D(pool_size=params['pool_size']))
    model.add(Dropout(params['dropout']))
    
    # Flatten and Dense layers
    model.add(Flatten())
    model.add(Dense(params['dense_units'], activation='relu'))
    model.add(Dropout(params['dropout']))
    model.add(Dense(1))  # Output layer
    
    # Compile
    model.compile(
        optimizer=Adam(learning_rate=0.001),
        loss='mse',
        metrics=['mae']
    )
    
    return model

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
# 2. TRAIN CNN MODELS
# =============================================================================
print("\n" + "=" * 70)
print("2. TRAIN CNN MODELS")
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
    
    # Create sequences for CNN
    seq_length = 24  # Use 24 time steps (2 hours) as input
    print(f"Creating sequences with length {seq_length}...")
    
    X_train_seq, y_train_seq = create_sequences(X_train_scaled, y_train.values, seq_length)
    X_test_seq, y_test_seq = create_sequences(X_test_scaled, y_test.values, seq_length)
    
    print(f"Train sequences: {X_train_seq.shape}")
    print(f"Test sequences: {X_test_seq.shape}")
    
    # For large datasets, use subset for training
    max_train_size = 5000
    if len(X_train_seq) > max_train_size:
        indices = np.random.choice(len(X_train_seq), max_train_size, replace=False)
        X_train_seq = X_train_seq[indices]
        y_train_seq = y_train_seq[indices]
        print(f"Using random subset of {max_train_size} sequences for training")
    
    try:
        # Build model
        input_shape = (seq_length, X_train.shape[1])
        cnn_params = {
            'filters': 64,
            'kernel_size': 3,
            'pool_size': 2,
            'dense_units': 128,
            'dropout': 0.2
        }
        
        print("Building CNN model...")
        model = build_cnn_model(input_shape, cnn_params)
        
        # Callbacks
        early_stopping = EarlyStopping(
            monitor='val_loss',
            patience=PATIENCE,
            restore_best_weights=True,
            verbose=1
        )
        
        reduce_lr = ReduceLROnPlateau(
            monitor='val_loss',
            factor=0.5,
            patience=5,
            min_lr=1e-6,
            verbose=1
        )
        
        # Train model
        print("Training CNN model...")
        history = model.fit(
            X_train_seq, y_train_seq,
            epochs=EPOCHS,
            batch_size=BATCH_SIZE,
            validation_split=0.2,
            callbacks=[early_stopping, reduce_lr],
            verbose=1
        )
        
        # Predict
        print("Generating predictions...")
        y_pred = model.predict(X_test_seq, verbose=0).flatten()
        
        # Calculate metrics
        metrics = calculate_metrics(y_test_seq, y_pred)
        
        print(f"\nMetrics:")
        print(f"  RMSE: {metrics['RMSE']:.4f}")
        print(f"  MAE:  {metrics['MAE']:.4f}")
        print(f"  MAPE: {metrics['MAPE']:.2f}%")
        
        results[horizon_name] = {
            'model': model,
            'scaler': scaler,
            'predictions': y_pred,
            'y_test': y_test_seq,
            'y_train': y_train_seq,
            'X_train_seq': X_train_seq,
            'X_test_seq': X_test_seq,
            'history': history.history,
            'metrics': metrics,
            'seq_length': seq_length,
            'params': cnn_params
        }
        
    except Exception as e:
        print(f"Error training CNN for {horizon_name}: {e}")
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
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    y_test = result['y_test']
    y_pred = result['predictions']
    history = result['history']
    metrics = result['metrics']
    
    # Plot 1: Actual vs Predicted
    ax1 = axes[0, 0]
    ax1.plot(y_test, label='Actual', alpha=0.8, linewidth=0.8)
    ax1.plot(y_pred, label='CNN Prediction', alpha=0.8, linewidth=0.8)
    ax1.set_title(f'CNN Forecast - {horizon_name.upper()} Horizon\n'
                  f'RMSE: {metrics["RMSE"]:.2f}, MAE: {metrics["MAE"]:.2f}, MAPE: {metrics["MAPE"]:.2f}%')
    ax1.set_xlabel('Time Steps')
    ax1.set_ylabel('Demand (MW)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Residuals
    ax2 = axes[0, 1]
    residuals = y_test - y_pred
    ax2.plot(residuals, alpha=0.7, linewidth=0.5)
    ax2.axhline(y=0, color='r', linestyle='--', alpha=0.5)
    ax2.set_title('Prediction Residuals')
    ax2.set_xlabel('Time Steps')
    ax2.set_ylabel('Residual (MW)')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Training History - Loss
    ax3 = axes[1, 0]
    ax3.plot(history['loss'], label='Train Loss')
    ax3.plot(history['val_loss'], label='Val Loss')
    ax3.set_title('Training History - Loss')
    ax3.set_xlabel('Epoch')
    ax3.set_ylabel('MSE Loss')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Training History - MAE
    ax4 = axes[1, 1]
    ax4.plot(history['mae'], label='Train MAE')
    ax4.plot(history['val_mae'], label='Val MAE')
    ax4.set_title('Training History - MAE')
    ax4.set_xlabel('Epoch')
    ax4.set_ylabel('MAE')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / f'cnn_{horizon_name}.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: cnn_{horizon_name}.png")

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
metrics_df.to_csv(RESULTS_DIR / 'cnn_metrics.csv')
print(f"Metrics saved: cnn_metrics.csv")

# Save predictions
for horizon_name, result in results.items():
    if result is None:
        continue
    
    pred_df = pd.DataFrame({
        'actual': result['y_test'],
        'prediction': result['predictions']
    })
    pred_df.to_csv(RESULTS_DIR / f'cnn_predictions_{horizon_name}.csv')
    print(f"Saved: cnn_predictions_{horizon_name}.csv")

# Save models
for horizon_name, result in results.items():
    if result is None:
        continue
    
    model_path = RESULTS_DIR / f'cnn_model_{horizon_name}.keras'
    result['model'].save(model_path)
    
    scaler_path = RESULTS_DIR / f'cnn_scaler_{horizon_name}.pkl'
    joblib.dump(result['scaler'], scaler_path)
    
    print(f"Saved: cnn_model_{horizon_name}.keras, cnn_scaler_{horizon_name}.pkl")

# Save training history
for horizon_name, result in results.items():
    if result is None:
        continue
    
    history_df = pd.DataFrame(result['history'])
    history_df.to_csv(RESULTS_DIR / f'cnn_history_{horizon_name}.csv')
    print(f"Saved: cnn_history_{horizon_name}.csv")

# Save model summary
with open(RESULTS_DIR / 'cnn_summary.txt', 'w') as f:
    f.write("CNN Model Results\n")
    f.write("=" * 50 + "\n\n")
    for horizon_name, result in results.items():
        if result is None:
            continue
        f.write(f"\n{horizon_name.upper()} Horizon:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Architecture: {result['params']}\n")
        f.write(f"Sequence length: {result['seq_length']}\n")
        f.write(f"RMSE: {result['metrics']['RMSE']:.4f}\n")
        f.write(f"MAE: {result['metrics']['MAE']:.4f}\n")
        f.write(f"MAPE: {result['metrics']['MAPE']:.2f}%\n")
print(f"Summary saved: cnn_summary.txt")

# =============================================================================
# FINAL SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY - CNN")
print("=" * 70)

print("\n" + metrics_df.to_string())