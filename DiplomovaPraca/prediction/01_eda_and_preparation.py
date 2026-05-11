import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Display settings
plt.style.use('seaborn-v0_8-whitegrid')
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

# =============================================================================
# 1. DATA LOADING
# =============================================================================
print("=" * 70)
print("1. DATA LOADING")
print("=" * 70)

DATA_PATH = Path(__file__).parent.parent / "data" / "Renewable Energy and Electricity Demand Time Series Dataset with Exogenous Variables at 5-minute Interval"
CSV_FILE = DATA_PATH / "renewable_energy_and_electricity_demand_time_series_dataset_with_exogenous_variables_at_5_minute_interval.csv"

df = pd.read_csv(CSV_FILE)
print(f"Dataset size: {df.shape[0]:,} rows × {df.shape[1]} columns")
print(f"\nColumns:\n{df.columns.tolist()}")

# =============================================================================
# 2. DATA STRUCTURE OVERVIEW
# =============================================================================
print("\n" + "=" * 70)
print("2. DATA STRUCTURE OVERVIEW")
print("=" * 70)

print("\nFirst 5 rows:")
print(df.head())

print("\nData types:")
print(df.dtypes)

print("\nNumeric columns statistics:")
print(df.describe().round(2))

# =============================================================================
# 3. TIMESTAMP PROCESSING
# =============================================================================
print("\n" + "=" * 70)
print("3. TIMESTAMP PROCESSING")
print("=" * 70)

# Find datetime column
datetime_col = None
for col in df.columns:
    if 'time' in col.lower() or 'date' in col.lower():
        datetime_col = col
        break

if datetime_col:
    print(f"Datetime column: {datetime_col}")
    df[datetime_col] = pd.to_datetime(df[datetime_col])
    df = df.set_index(datetime_col)
    df = df.sort_index()
    
    print(f"Data period: {df.index.min()} — {df.index.max()}")
    print(f"Duration: {(df.index.max() - df.index.min()).days} days")
    
    # Check interval
    intervals = df.index.to_series().diff().dropna()
    print(f"Median interval: {intervals.median()}")
else:
    print("Datetime column not found!")

# =============================================================================
# 4. MISSING VALUES ANALYSIS
# =============================================================================
print("\n" + "=" * 70)
print("4. MISSING VALUES ANALYSIS")
print("=" * 70)

missing = df.isnull().sum()
missing_pct = (missing / len(df) * 100).round(2)
missing_df = pd.DataFrame({
    'Missing': missing,
    'Percent': missing_pct
})
missing_df = missing_df[missing_df['Missing'] > 0]

if len(missing_df) > 0:
    print(missing_df)
else:
    print("No missing values!")

# =============================================================================
# 5. TARGET VARIABLE IDENTIFICATION
# =============================================================================
print("\n" + "=" * 70)
print("5. TARGET VARIABLE IDENTIFICATION")
print("=" * 70)

# Search for consumption column
target_keywords = ['demand', 'consumption', 'load', 'energy']
target_col = None

for col in df.columns:
    col_lower = col.lower()
    if any(kw in col_lower for kw in target_keywords):
        target_col = col
        print(f"Target column found: {col}")
        break

if target_col:
    print(f"\nStatistics for {target_col}:")
    print(df[target_col].describe())
else:
    print("Target column not found! Select manually.")
    print(f"Available columns: {df.columns.tolist()}")

# =============================================================================
# 6. TIME SERIES VISUALIZATION
# =============================================================================
print("\n" + "=" * 70)
print("6. VISUALIZATION")
print("=" * 70)

# Create plots directory
PLOTS_DIR = Path(r"C:\Users\mnis\Desktop\DP\prediction\plots")
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

if target_col:
    fig, axes = plt.subplots(3, 1, figsize=(14, 10))
    
    # Full series
    axes[0].plot(df.index, df[target_col], linewidth=0.5, alpha=0.7)
    axes[0].set_title(f'Full Time Series: {target_col}')
    axes[0].set_xlabel('Date')
    axes[0].set_ylabel('MW')
    
    # Last week
    last_week = df[target_col].last('7D')
    axes[1].plot(last_week.index, last_week.values, linewidth=1)
    axes[1].set_title('Last Week')
    axes[1].set_xlabel('Date')
    axes[1].set_ylabel('MW')
    
    # Last day
    last_day = df[target_col].last('1D')
    axes[2].plot(last_day.index, last_day.values, linewidth=1.5, marker='o', markersize=3)
    axes[2].set_title('Last Day')
    axes[2].set_xlabel('Time')
    axes[2].set_ylabel('MW')
    
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / '01_time_series_overview.png', dpi=150)
    plt.show()
    print(f"Plot saved: {PLOTS_DIR / '01_time_series_overview.png'}")

# =============================================================================
# 7. CORRELATION ANALYSIS
# =============================================================================
print("\n" + "=" * 70)
print("7. CORRELATION ANALYSIS")
print("=" * 70)

# Numeric columns correlation
numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
print(f"Numeric columns: {numeric_cols}")

if len(numeric_cols) > 1:
    corr_matrix = df[numeric_cols].corr()
    
    # Correlation with target variable
    if target_col in corr_matrix.columns:
        target_corr = corr_matrix[target_col].sort_values(ascending=False)
        print(f"\nCorrelation with {target_col}:")
        print(target_corr)
    
    # Heatmap
    fig, ax = plt.subplots(figsize=(12, 10))
    sns.heatmap(corr_matrix, annot=True, cmap='RdBu_r', center=0, 
                fmt='.2f', ax=ax, square=True)
    ax.set_title('Correlation Matrix')
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / '02_correlation_matrix.png', dpi=150)
    plt.show()
    print(f"Plot saved: {PLOTS_DIR / '02_correlation_matrix.png'}")

# =============================================================================
# 8. SEASONALITY ANALYSIS
# =============================================================================
print("\n" + "=" * 70)
print("8. SEASONALITY ANALYSIS")
print("=" * 70)

if target_col and isinstance(df.index, pd.DatetimeIndex):
    # Create time features
    datetime_index = pd.DatetimeIndex(df.index)
    df['hour'] = datetime_index.hour
    df['dayofweek'] = datetime_index.dayofweek
    df['month'] = datetime_index.month
    df['is_weekend'] = df['dayofweek'].isin([5, 6]).astype(int)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # By hour
    hourly = df.groupby('hour')[target_col].mean()
    axes[0, 0].bar(hourly.index, hourly.values, color='steelblue')
    axes[0, 0].set_title('Average Consumption by Hour')
    axes[0, 0].set_xlabel('Hour')
    axes[0, 0].set_ylabel('MW')
    
    # By day of week
    daily = df.groupby('dayofweek')[target_col].mean()
    days = ['Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun']
    axes[0, 1].bar(range(7), daily.values, color='coral')
    axes[0, 1].set_xticks(range(7))
    axes[0, 1].set_xticklabels(days)
    axes[0, 1].set_title('Average Consumption by Day of Week')
    axes[0, 1].set_ylabel('MW')
    
    # By month
    monthly = df.groupby('month')[target_col].mean()
    axes[1, 0].bar(monthly.index, monthly.values, color='forestgreen')
    axes[1, 0].set_title('Average Consumption by Month')
    axes[1, 0].set_xlabel('Month')
    axes[1, 0].set_ylabel('MW')
    
    # Weekday vs Weekend
    weekend_data = df.groupby('is_weekend')[target_col].mean()
    axes[1, 1].bar(['Weekday', 'Weekend'], weekend_data.values, color=['steelblue', 'coral'])
    axes[1, 1].set_title('Weekday vs Weekend')
    axes[1, 1].set_ylabel('MW')
    
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / '03_seasonality_analysis.png', dpi=150)
    plt.show()
    print(f"Plot saved: {PLOTS_DIR / '03_seasonality_analysis.png'}")

# =============================================================================
# 9. OUTLIERS ANALYSIS
# =============================================================================
print("\n" + "=" * 70)
print("9. OUTLIERS ANALYSIS")
print("=" * 70)

if target_col:
    Q1 = df[target_col].quantile(0.25)
    Q3 = df[target_col].quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    
    outliers = df[(df[target_col] < lower_bound) | (df[target_col] > upper_bound)]
    print(f"Bounds: [{lower_bound:.2f}, {upper_bound:.2f}]")
    print(f"Outliers: {len(outliers)} ({len(outliers)/len(df)*100:.2f}%)")
    
    # Boxplot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.boxplot(df[target_col].dropna())
    ax.set_title(f'Boxplot: {target_col}')
    ax.set_ylabel('MW')
    plt.savefig(PLOTS_DIR / '04_boxplot.png', dpi=150)
    plt.show()

# =============================================================================
# 10. SAVE PREPARED DATA
# =============================================================================
print("\n" + "=" * 70)
print("10. SAVE PREPARED DATA")
print("=" * 70)

OUTPUT_DIR = Path(r"C:\Users\mnis\Desktop\DP\prediction\data")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Save dataset with time features
df.to_csv(OUTPUT_DIR / 'prepared_data.csv')
print(f"Data saved: {OUTPUT_DIR / 'prepared_data.csv'}")

# Summary
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"Size: {df.shape[0]:,} rows × {df.shape[1]} columns")
print(f"Period: {df.index.min()} — {df.index.max()}")
print(f"Target variable: {target_col}")
print(f"Features: {numeric_cols}")
if 'hour' in df.columns:
    print(f"Added time features: hour, dayofweek, month, is_weekend")
