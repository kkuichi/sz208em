"""
=============================================================================
Advanced SOTA Models for Electricity Demand Forecasting
=============================================================================
N-BEATS and Temporal Fusion Transformer (TFT) implementations
Requires: neuralforecast[torch] and/or pytorch-forecasting
=============================================================================
"""

import os
import logging
import pandas as pd
import numpy as np
from pathlib import Path
import pickle
import warnings
from typing import Tuple

# Suppress Python warnings
warnings.filterwarnings('ignore')

# Reduce TensorFlow/oneDNN/logging noise before importing TF or other libs
# 0 = all messages, 1 = INFO, 2 = WARNING, 3 = ERROR only
os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL', '3')
# Disable oneDNN info messages
os.environ.setdefault('TF_ENABLE_ONEDNN_OPTS', '0')

# Reduce logging verbosity
logging.getLogger('tensorflow').setLevel(logging.ERROR)
logging.getLogger('neuralforecast').setLevel(logging.ERROR)
logger = logging.getLogger(__name__)

from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, mean_absolute_error
import matplotlib.pyplot as plt

# =============================================================================
# CONFIGURATION
# =============================================================================
SPLITS_DIR = Path(__file__).parent / "data/splits"
RESULTS_DIR = Path(__file__).parent / "results/advanced"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

TARGET_COL = 'demand_MW'

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

def prepare_data_for_neuralforecast(splits, horizon_name) -> Tuple[pd.DataFrame, int]:
    """Prepare data in NeuralForecast format"""
    data = splits[horizon_name]
    X_train = data['X_train']
    y_train = data['y_train']
    X_test = data['X_test']
    y_test = data['y_test']

    # Create datetime index if needed
    if not isinstance(y_train.index, pd.DatetimeIndex):
        # Assume 5-min frequency
        start_date = pd.Timestamp('2020-01-01')
        y_train.index = pd.date_range(start=start_date, periods=len(y_train), freq='5min')
        y_test.index = pd.date_range(start=y_train.index[-1] + pd.Timedelta('5min'),
                                   periods=len(y_test), freq='5min')

    # Combine train and test for full series
    y_full = pd.concat([y_train, y_test])
    X_full = pd.concat([X_train, X_test])

    # Create NeuralForecast DataFrame
    df = pd.DataFrame({
        'unique_id': 'demand',
        'ds': y_full.index,
        'y': y_full.values
    })

    # Add exogenous variables
    for col in X_full.columns:
        df[f'exog_{col}'] = X_full[col].values

    return df, len(y_test)

# =============================================================================
# N-BEATS MODEL
# =============================================================================
def train_nbeats_model(splits, horizon_name):
    """Train N-BEATS model using NeuralForecast"""
    print(f"Training N-BEATS for {horizon_name}...")

    try:
        from neuralforecast import NeuralForecast
        from neuralforecast.models import NBEATS
        from neuralforecast.losses.pytorch import MAE

        # Prepare data
        df, horizon_size = prepare_data_for_neuralforecast(splits, horizon_name)

        # Model configuration
        # Determine validation size (integer) for early stopping based on training portion
        train_len = max(1, len(df) - horizon_size)
        val_size_int = max(1, int(0.1 * train_len))
        print(f"Using val_size={val_size_int} for N-BEATS (train_len={train_len}, horizon={horizon_size})")

        model = NBEATS(
            h=horizon_size,  # Forecast horizon
            input_size=2*horizon_size,  # Input size (2x horizon for context)
            loss=MAE(),
            max_steps=100,  # Limited for demo
            early_stop_patience_steps=5,
            val_check_steps=10,
            val_size=val_size_int,  # Integer > 0 required when early stopping is enabled
            learning_rate=1e-3,
            random_seed=42
        )

        # Create and fit NeuralForecast
        nf = NeuralForecast(
            models=[model],
            freq='5min'
        )

        # Fit model (wrapped to capture and report errors)
        try:
            nf.fit(df=df, id_col='unique_id', time_col='ds', target_col='y')
        except KeyboardInterrupt:
            print("Training N-BEATS interrupted by user.")
            return None
        except Exception as e:
            print(f"Error fitting N-BEATS model: {e}")
            return None

        # Predict
        try:
            predictions = nf.predict()
            # neuralforecast may name the column after the model; try common keys
            if 'NBEATS' in predictions.columns:
                y_pred = predictions['NBEATS'].values
            else:
                # fallback to last column
                y_pred = predictions.iloc[:, -1].values
        except Exception as e:
            print(f"Error predicting with N-BEATS: {e}")
            return None

        # Get actual test values
        data = splits[horizon_name]
        y_test = data['y_test'].values[-horizon_size:]  # Last horizon_size values

        # Calculate metrics
        metrics = calculate_metrics(y_test, y_pred)

        print(f"N-BEATS {horizon_name} - RMSE: {metrics['RMSE']:.4f}, MAE: {metrics['MAE']:.4f}")

        # Save results
        results_df = pd.DataFrame({
            'actual': y_test,
            'prediction': y_pred
        })
        results_df.to_csv(RESULTS_DIR / f'nbeats_predictions_{horizon_name}.csv')

        return {
            'model': nf,
            'predictions': y_pred,
            'y_test': y_test,
            'metrics': metrics
        }

    except ImportError:
        print("NeuralForecast not available. Install with: pip install neuralforecast[torch]")
        return None
    except Exception as e:
        print(f"Error training N-BEATS: {e}")
        return None

# =============================================================================
# SIMPLIFIED TFT-LIKE MODEL
# =============================================================================
def train_simple_tft_model(splits, horizon_name):
    """Train a simplified TFT-like model using basic neural networks"""
    print(f"Training Simple TFT-like for {horizon_name}...")

    try:
        import tensorflow as tf
        from tensorflow.keras.models import Model
        from tensorflow.keras.layers import Input, Dense, Dropout, BatchNormalization, Concatenate
        from tensorflow.keras.optimizers import Adam
        from tensorflow.keras.callbacks import EarlyStopping

        # Suppress TF warnings
        tf.get_logger().setLevel('ERROR')

        data = splits[horizon_name]
        X_train = data['X_train']
        y_train = data['y_train']
        X_test = data['X_test']
        y_test = data['y_test']

        # Scale features
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)

        # For TFT-like: separate static and dynamic features
        # Assume first few columns are static, rest are dynamic
        n_static = min(5, X_train.shape[1] // 3)  # Assume some static features
        n_dynamic = X_train.shape[1] - n_static

        # Build simplified TFT-like model
        # Static features input
        static_input = Input(shape=(n_static,), name='static_input')
        static_dense = Dense(32, activation='relu')(static_input)
        static_dense = BatchNormalization()(static_dense)
        static_dense = Dropout(0.1)(static_dense)

        # Dynamic features input
        dynamic_input = Input(shape=(n_dynamic,), name='dynamic_input')
        dynamic_dense = Dense(64, activation='relu')(dynamic_input)
        dynamic_dense = BatchNormalization()(dynamic_dense)
        dynamic_dense = Dropout(0.1)(dynamic_dense)

        # Combine
        combined = Concatenate()([static_dense, dynamic_dense])
        combined = Dense(64, activation='relu')(combined)
        combined = BatchNormalization()(combined)
        combined = Dropout(0.1)(combined)
        output = Dense(1)(combined)

        model = Model(inputs=[static_input, dynamic_input], outputs=output)
        model.compile(optimizer=Adam(learning_rate=0.001), loss='mse', metrics=['mae'])

        # Split data
        X_train_static = X_train_scaled[:, :n_static]
        X_train_dynamic = X_train_scaled[:, n_static:]
        X_test_static = X_test_scaled[:, :n_static]
        X_test_dynamic = X_test_scaled[:, n_static:]

        # Train (catch KeyboardInterrupt and other exceptions)
        early_stop = EarlyStopping(monitor='val_loss', patience=10, restore_best_weights=True, verbose=0)
        try:
            history = model.fit(
                [X_train_static, X_train_dynamic], y_train,
                epochs=50,
                batch_size=32,
                validation_split=0.1,
                callbacks=[early_stop],
                verbose=0
            )
        except KeyboardInterrupt:
            print("Training Simple TFT interrupted by user.")
            return None
        except Exception as e:
            print(f"Error during Simple TFT training: {e}")
            return None

        # Predict
        try:
            y_pred = model.predict([X_test_static, X_test_dynamic], verbose=0).flatten()
        except Exception as e:
            print(f"Error during Simple TFT prediction: {e}")
            return None

        # Calculate metrics
        metrics = calculate_metrics(y_test.values, y_pred)

        print(f"Simple TFT {horizon_name} - RMSE: {metrics['RMSE']:.4f}, MAE: {metrics['MAE']:.4f}")

        # Save results
        results_df = pd.DataFrame({
            'actual': y_test,
            'prediction': y_pred
        })
        results_df.to_csv(RESULTS_DIR / f'simple_tft_predictions_{horizon_name}.csv')

        return {
            'model': model,
            'scaler': scaler,
            'predictions': y_pred,
            'y_test': y_test,
            'metrics': metrics,
            'n_static': n_static,
            'n_dynamic': n_dynamic
        }

    except Exception as e:
        print(f"Error training Simple TFT: {e}")
        return None

# =============================================================================
# MAIN EXECUTION
# =============================================================================
print("=" * 70)
print("ADVANCED SOTA MODELS")
print("=" * 70)

# Load data
with open(SPLITS_DIR / "train_test_splits.pkl", 'rb') as f:
    splits = pickle.load(f)

print(f"Loaded splits for horizons: {list(splits.keys())}")

# Train models
results = {}

try:
    for horizon in ['day', 'week', 'month']:
        print(f"\n{'='*50}")
        print(f"HORIZON: {horizon.upper()}")
        print(f"{'='*50}")

        horizon_results = {}

        # Try N-BEATS
        nbeats_result = train_nbeats_model(splits, horizon)
        if nbeats_result:
            horizon_results['N-BEATS'] = nbeats_result

        # Try Simple TFT
        tft_result = train_simple_tft_model(splits, horizon)
        if tft_result:
            horizon_results['Simple_TFT'] = tft_result

        results[horizon] = horizon_results
except KeyboardInterrupt:
    print('\nExecution interrupted by user. Saving any available results...')
except Exception as e:
    print(f"Unexpected error during model training loop: {e}")
    logger.exception(e)

# =============================================================================
# SAVE SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("SAVE SUMMARY")
print("=" * 70)

# Save metrics
for model_name in ['N-BEATS', 'Simple_TFT']:
    model_metrics = {}
    for horizon, horizon_results in results.items():
        if model_name in horizon_results:
            model_metrics[horizon] = horizon_results[model_name]['metrics']

    if model_metrics:
        metrics_df = pd.DataFrame(model_metrics).T
        metrics_df.index.name = 'Horizon'
        metrics_df.to_csv(RESULTS_DIR / f'{model_name.lower().replace("-", "_")}_metrics.csv')
        print(f"Saved: {model_name.lower().replace('-', '_')}_metrics.csv")

# Save summary
with open(RESULTS_DIR / 'advanced_models_summary.txt', 'w') as f:
    f.write("Advanced SOTA Models Results\n")
    f.write("=" * 50 + "\n\n")

    for horizon, horizon_results in results.items():
        f.write(f"{horizon.upper()} Horizon:\n")
        f.write("-" * 30 + "\n")

        for model_name, result in horizon_results.items():
            f.write(f"\n{model_name}:\n")
            f.write(f"  RMSE: {result['metrics']['RMSE']:.4f}\n")
            f.write(f"  MAE: {result['metrics']['MAE']:.4f}\n")
            f.write(f"  MAPE: {result['metrics']['MAPE']:.2f}%\n")

        f.write("\n")

print("Summary saved: advanced_models_summary.txt")

# =============================================================================
# FINAL SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY - ADVANCED MODELS")
print("=" * 70)

for horizon, horizon_results in results.items():
    print(f"\n{horizon.upper()}:")
    for model_name, result in horizon_results.items():
        metrics = result['metrics']
        print(f"  {model_name}: RMSE={metrics['RMSE']:.4f}, MAE={metrics['MAE']:.4f}, MAPE={metrics['MAPE']:.2f}%")

print(f"\nAll results saved in: {RESULTS_DIR}")
print("\nNote: N-BEATS requires 'pip install neuralforecast[torch]'")
print("Simple TFT is a basic implementation without full TFT complexity")
