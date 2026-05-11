"""
=============================================================================
LSTM and TCN Models for Electricity Demand Forecasting
=============================================================================
LSTM seq2seq (many-to-one) and Temporal Convolutional Network (TCN)
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
from tensorflow.keras.layers import LSTM, Dense, Conv1D, MaxPooling1D, Flatten, Dropout, BatchNormalization
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
RESULTS_DIR = Path(__file__).parent / "results/lstm_tcn"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

TARGET_COL = 'demand_MW'

# Model parameters
SEQ_LEN = 24  # 2 hours of 5-min intervals
EPOCHS = 100
BATCH_SIZE = 64
PATIENCE = 10
VALIDATION_SPLIT = 0.1

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

def create_sequences(X, y, seq_len):
    """Create sequences for time series forecasting"""
    X_seq, y_seq = [], []
    for i in range(len(X) - seq_len + 1):
        X_seq.append(X[i:i+seq_len])
        y_seq.append(y[i+seq_len-1])
    return np.array(X_seq), np.array(y_seq)

def build_lstm_model(input_shape):
    """Build LSTM seq2seq model"""
    model = Sequential([
        LSTM(128, input_shape=input_shape, return_sequences=True),
        BatchNormalization(),
        Dropout(0.2),
        
        LSTM(64, return_sequences=False),
        BatchNormalization(),
        Dropout(0.2),
        
        Dense(64, activation='relu'),
        Dropout(0.1),
        Dense(1)
    ])
    
    model.compile(
        optimizer=Adam(learning_rate=0.001),
        loss='mse',
        metrics=['mae']
    )
    
    return model

def build_tcn_model(input_shape):
    """Build simplified TCN-like model with dilated convolutions"""
    model = Sequential()
    
    # First dilated conv block
    model.add(Conv1D(64, kernel_size=3, dilation_rate=1, padding='causal', 
                     input_shape=input_shape, activation='relu'))
    model.add(BatchNormalization())
    model.add(Dropout(0.2))
    
    # Second dilated conv block
    model.add(Conv1D(64, kernel_size=3, dilation_rate=2, padding='causal', activation='relu'))
    model.add(BatchNormalization())
    model.add(Dropout(0.2))
    
    # Third dilated conv block
    model.add(Conv1D(64, kernel_size=3, dilation_rate=4, padding='causal', activation='relu'))
    model.add(BatchNormalization())
    model.add(Dropout(0.2))
    
    # Global pooling and dense
    model.add(Flatten())
    model.add(Dense(64, activation='relu'))
    model.add(Dropout(0.1))
    model.add(Dense(1))
    
    model.compile(
        optimizer=Adam(learning_rate=0.001),
        loss='mse',
        metrics=['mae']
    )
    
    return model

def plot_training_history(history, model_name, horizon, save_path):
    """Plot training history"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Loss
    axes[0].plot(history.history['loss'], label='Train Loss')
    axes[0].plot(history.history['val_loss'], label='Val Loss')
    axes[0].set_title(f'{model_name} - Loss')
    axes[0].set_xlabel('Epoch')
    axes[0].set_ylabel('MSE Loss')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # MAE
    axes[1].plot(history.history['mae'], label='Train MAE')
    axes[1].plot(history.history['val_mae'], label='Val MAE')
    axes[1].set_title(f'{model_name} - MAE')
    axes[1].set_xlabel('Epoch')
    axes[1].set_ylabel('MAE')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.suptitle(f'{model_name} Training History - {horizon.upper()} Horizon')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Training history saved: {save_path}")

def plot_predictions(y_true, y_pred, model_name, horizon, metrics, save_path):
    """Plot predictions vs actual"""
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    
    # Predictions plot (subset for visibility)
    max_points = min(200, len(y_true))
    indices = np.arange(max_points)
    
    axes[0].plot(indices, y_true[:max_points], label='Actual', alpha=0.8, linewidth=1.5)
    axes[0].plot(indices, y_pred[:max_points], label='Predicted', alpha=0.8, linewidth=1)
    axes[0].set_title(f'{model_name} Predictions - {horizon.upper()} Horizon\n'
                     f'RMSE: {metrics["RMSE"]:.2f}, MAE: {metrics["MAE"]:.2f}, MAPE: {metrics["MAPE"]:.2f}%')
    axes[0].set_xlabel('Time Steps')
    axes[0].set_ylabel('Demand (MW)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Residuals
    residuals = y_true - y_pred
    axes[1].plot(indices, residuals[:max_points], alpha=0.7, linewidth=0.5)
    axes[1].axhline(y=0, color='r', linestyle='--', alpha=0.5)
    axes[1].set_title('Prediction Residuals')
    axes[1].set_xlabel('Time Steps')
    axes[1].set_ylabel('Residual (MW)')
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Predictions plot saved: {save_path}")

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
# 2. TRAIN MODELS
# =============================================================================
print("\n" + "=" * 70)
print("2. TRAIN LSTM AND TCN MODELS")
print("=" * 70)

models_to_train = ['LSTM', 'TCN']
results = {}

for model_type in models_to_train:
    print(f"\n{'='*60}")
    print(f"MODEL TYPE: {model_type}")
    print(f"{'='*60}")
    
    model_results = {}
    
    for horizon_name, data in splits.items():
        print(f"\n--- {horizon_name.upper()} HORIZON ---")
        
        X_train = data['X_train']
        y_train = data['y_train']
        X_test = data['X_test']
        y_test = data['y_test']
        
        print(f"Train size: {len(X_train):,}")
        print(f"Test size: {len(X_test):,}")
        print(f"Features: {len(X_train.columns)}")
        
        # Scale features
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        # Create sequences
        print(f"Creating sequences with length {SEQ_LEN}...")
        X_train_seq, y_train_seq = create_sequences(X_train_scaled, y_train.values, SEQ_LEN)
        X_test_seq, y_test_seq = create_sequences(X_test_scaled, y_test.values, SEQ_LEN)
        
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
            input_shape = (SEQ_LEN, X_train.shape[1])
            
            if model_type == 'LSTM':
                model = build_lstm_model(input_shape)
            elif model_type == 'TCN':
                model = build_tcn_model(input_shape)
            
            print(f"Built {model_type} model")
            
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
            print(f"Training {model_type}...")
            history = model.fit(
                X_train_seq, y_train_seq,
                epochs=EPOCHS,
                batch_size=BATCH_SIZE,
                validation_split=VALIDATION_SPLIT,
                callbacks=[early_stopping, reduce_lr],
                verbose=1
            )
            
            # Predict
            print("Generating predictions...")
            y_pred = model.predict(X_test_seq, verbose=0).flatten()
            
            # Calculate metrics
            metrics = calculate_metrics(y_test_seq, y_pred)
            
            print(f"\n{model_type} Metrics:")
            print(f"  RMSE: {metrics['RMSE']:.4f}")
            print(f"  MAE:  {metrics['MAE']:.4f}")
            print(f"  MAPE: {metrics['MAPE']:.2f}%")
            
            # Plot training history
            history_path = RESULTS_DIR / f'{model_type.lower()}_history_{horizon_name}.png'
            plot_training_history(history, model_type, horizon_name, history_path)
            
            # Plot predictions
            pred_path = RESULTS_DIR / f'{model_type.lower()}_predictions_{horizon_name}.png'
            plot_predictions(y_test_seq, y_pred, model_type, horizon_name, metrics, pred_path)
            
            # Save model and scaler
            model_path = RESULTS_DIR / f'{model_type.lower()}_{horizon_name}.keras'
            model.save(model_path)
            
            scaler_path = RESULTS_DIR / f'{model_type.lower()}_scaler_{horizon_name}.pkl'
            joblib.dump(scaler, scaler_path)
            
            # Save predictions
            pred_df = pd.DataFrame({
                'actual': y_test_seq,
                'prediction': y_pred
            })
            pred_df.to_csv(RESULTS_DIR / f'{model_type.lower()}_predictions_{horizon_name}.csv')
            
            print(f"Saved: {model_type.lower()}_{horizon_name}.keras, scaler, predictions")
            
            model_results[horizon_name] = {
                'model': model,
                'scaler': scaler,
                'predictions': y_pred,
                'y_test': y_test_seq,
                'history': history.history,
                'metrics': metrics,
                'seq_length': SEQ_LEN
            }
            
        except Exception as e:
            print(f"Error training {model_type} for {horizon_name}: {e}")
            model_results[horizon_name] = None
    
    results[model_type] = model_results

# =============================================================================
# 3. SAVE SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("3. SAVE SUMMARY")
print("=" * 70)

# Save metrics summary for each model type
for model_type, model_results in results.items():
    metrics_df = pd.DataFrame({
        horizon: result['metrics'] for horizon, result in model_results.items() if result is not None
    }).T
    metrics_df.index.name = 'Horizon'
    metrics_df.to_csv(RESULTS_DIR / f'{model_type.lower()}_metrics.csv')
    print(f"Metrics saved: {model_type.lower()}_metrics.csv")

# Save combined summary
with open(RESULTS_DIR / 'lstm_tcn_summary.txt', 'w') as f:
    f.write("LSTM and TCN Model Results\n")
    f.write("=" * 50 + "\n\n")
    
    for model_type, model_results in results.items():
        f.write(f"\n{model_type} MODELS:\n")
        f.write("=" * 30 + "\n")
        
        for horizon_name, result in model_results.items():
            if result is None:
                continue
            f.write(f"\n{horizon_name.upper()} Horizon:\n")
            f.write("-" * 20 + "\n")
            f.write(f"Sequence length: {result['seq_length']}\n")
            f.write(f"RMSE: {result['metrics']['RMSE']:.4f}\n")
            f.write(f"MAE: {result['metrics']['MAE']:.4f}\n")
            f.write(f"MAPE: {result['metrics']['MAPE']:.2f}%\n")

print(f"Summary saved: lstm_tcn_summary.txt")

# =============================================================================
# FINAL SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY - LSTM & TCN")
print("=" * 70)

for model_type, model_results in results.items():
    print(f"\n{model_type} Results:")
    metrics_df = pd.DataFrame({
        horizon: result['metrics'] for horizon, result in model_results.items() if result is not None
    }).T
    print(metrics_df.to_string())
