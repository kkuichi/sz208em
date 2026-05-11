#!/usr/bin/env python3
"""
Master script to run all prediction models sequentially.
Runs all experiments from EDA to ensemble comparison.
Estimated runtime: 3-5 hours
"""

import subprocess
import sys
import time
from pathlib import Path

# List of scripts to run in order
SCRIPTS = [
    ("01_eda_and_preparation.py", "Exploratory Data Analysis"),
    ("02_train_test_split.py", "Data Preparation & Feature Engineering"),
    ("03_arima_model.py", "ARIMA Model Training"),
    ("04_svm_model.py", "SVM Model Training"),
    ("05_cnn_model.py", "CNN Model Training"),
    ("07_lightgbm.py", "LightGBM Model Training"),
    ("08_lstm_tcn.py", "LSTM & TCN Model Training"),
    ("09_advanced_models.py", "Advanced Models (N-BEATS, TFT)"),
    ("06_comparison.py", "Model Comparison & Visualization"),
    ("10_hybrid_stacking.py", "Ensemble Methods"),
]

def run_script(script_name, description):
    """Run a single Python script and return success status."""
    print("\n" + "="*80)
    print(f"RUNNING: {description}")
    print(f"Script: {script_name}")
    print("="*80)
    
    try:
        result = subprocess.run(
            [sys.executable, script_name],
            cwd=Path(__file__).parent / "prediction",
            check=True,
            text=True
        )
        print(f"✓ SUCCESS: {description} completed!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ ERROR in {description}: {e}")
        return False
    except Exception as e:
        print(f"✗ UNEXPECTED ERROR in {description}: {e}")
        return False

def main():
    """Run all scripts sequentially."""
    print("\n" + "="*80)
    print("ELECTRICITY CONSUMPTION FORECASTING - FULL PIPELINE")
    print("="*80)
    print(f"Total scripts to run: {len(SCRIPTS)}")
    print(f"Estimated runtime: 3-5 hours")
    print("="*80 + "\n")
    
    start_time = time.time()
    results = {}
    
    for script_name, description in SCRIPTS:
        success = run_script(script_name, description)
        results[description] = success
        
        if not success:
            print(f"\n⚠️  Warning: {description} failed. Continuing with next script...\n")
            time.sleep(2)
    
    # Print summary
    elapsed_time = time.time() - start_time
    hours = elapsed_time // 3600
    minutes = (elapsed_time % 3600) // 60
    
    print("\n" + "="*80)
    print("PIPELINE EXECUTION SUMMARY")
    print("="*80)
    print(f"Total Runtime: {int(hours)}h {int(minutes)}m")
    print(f"\nResults:")
    for description, success in results.items():
        status = "✓ PASS" if success else "✗ FAIL"
        print(f"  {status}: {description}")
    
    successful = sum(1 for success in results.values() if success)
    print(f"\nCompleted: {successful}/{len(SCRIPTS)} scripts")
    print("="*80)
    
    if successful == len(SCRIPTS):
        print("\n✓ All experiments completed successfully!")
        print("Check prediction/results/ for outputs and visualizations.")
        return 0
    else:
        print(f"\n⚠️  {len(SCRIPTS) - successful} script(s) failed.")
        print("Check error messages above for details.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
