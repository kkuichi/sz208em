# Prediktívne metódy analýzy dát spotreby elektrickej energie

## Prehľad

Tento projekt implementuje a porovnáva viaceré prediktívne modely na prognózovanie spotreby elektrickej energie v rôznych časových horizontoch (deň, týždeň, mesiac). Dataset pozostáva z meraní s rozlíšením 5 minút s exogénnymi premennými počasia a údajmi o generácii obnoviteľnej energie.

**Projekt:** Diplomová práca  
**Téma:** Prediktívne metódy analýzy dát spotreby elektrickej energie  
**Autor:** Bc. Serghei Zabirchenko  
**Rok:** 2026

---

## Štruktúra projektu

```
prediction/
├── 01_eda_and_preparation.py       # Prieskumná analýza údajov (EDA)
├── 02_train_test_split.py          # Príprava údajov a inžinierstvo atribútov
├── 03_arima_model.py               # Klasický ARIMA model
├── 04_svm_model.py                 # Support Vector Machine regresia
├── 05_cnn_model.py                 # 1D Konvolučná neurónová sieť
├── 06_comparison.py                # Porovnanie modelov a vizualizácia
├── 07_lightgbm.py                  # LightGBM gradient boosting
├── 08_lstm_tcn.py                  # LSTM a Temporal Convolutional Networks
├── 09_advanced_models.py           # N-BEATS a zjednodušený Temporal Fusion Transformer
├── 10_hybrid_stacking.py           # Ensemble metódy (priemer, medián, Ridge stacking)
├── plots/                          # Vizualizácie EDA
│   ├── 01_time_series_overview.png
│   ├── 02_correlation_matrix.png
│   ├── 03_seasonality_analysis.png
│   └── 04_boxplot.png
├── results/                        # Predikcie modelov a metriky výkonu
│   ├── arima/
│   ├── svm/
│   ├── cnn/
│   ├── lightgbm/
│   ├── lstm_tcn/
│   ├── advanced/
│   └── comparison/
│       ├── ensemble_visualization_*.png
│       ├── predictions_*.png
│       ├── metrics_comparison.png
│       └── model_ranking.png
└── data/
    └── [pripravené datasety atribútov pre každý horizont]
```

---

## Inštalácia

### 1. Požiadavky

- **Python 3.8 alebo vyšší**
- **Anaconda alebo Miniconda** (odporúčané pre správu prostredia)
- Aspoň 8 GB RAM pre trénovanie modelov hlbokého učenia
- Windows / Linux / macOS

### 2. Klonujte repozitár

```bash
git clone https://github.com/kkuichi/sz208em.git
cd DP
```

### 3. Vytvorte a aktivujte prostredie

**Používanie Conda (odporúčané):**
```bash
conda create -n energy-forecast python=3.10 -y
conda activate energy-forecast
```

**Používanie venv:**
```bash
python -m venv venv
source venv/bin/activate  # Na Windows: venv\Scripts\activate
```

### 4. Nainštalujte závislosti

```bash
pip install pandas numpy matplotlib seaborn scikit-learn tensorflow keras statsmodels lightgbm optuna
pip install neuralforecast  # Voliteľne, pre N-BEATS model
```

**Úplný requirements.txt:**
```
pandas==2.0.0
numpy==1.24.0
matplotlib==3.7.0
seaborn==0.12.0
scikit-learn==1.2.0
tensorflow==2.13.0
keras==2.13.0
statsmodels==0.14.0
lightgbm==4.0.0
optuna==3.0.0
```

Uložte do `requirements.txt` a spustite:
```bash
pip install -r requirements.txt
```

---

## Údaje

### O datasete

Tento projekt používa údaje časových radov o **spotrebe elektrickej energie s exogénnymi premennými počasia**. Modely sú trénované na údajoch s **rozlíšením 5 minút**, ktoré je vhodné na krátkodobú prognózu zaťaženia elektrickej siete (STLF).

### Požadované charakteristiky datasetu

Pre úspešné reprodukovanie týchto experimentov potrebujete dataset s:

**Základné požiadavky:**
- **Údaje časových radov:** Merania spotreby/dopytu elektrickej energie v čase
- **Časové rozlíšenie:** 5-minútové, 15-minútové, 30-minútové alebo hodinové intervaly (5-minútové sú preferované)
- **Trvanie:** Aspoň 6-12 mesiacov súvislých údajov (čím dlhšie, tým lepšie)
- **Cieľová premenná:** Elektrické zaťaženie/dopyt v MW, kW alebo podobných jednotkách
- **Exogénne premenné (Počasie):**
  - Teplota (alebo aspoň vonkajšia teplota)
  - Vlhkosť (voliteľné, ale odporúčané)
  - Rýchlosť vetra (voliteľné)
  - Slnečné žiarenie (DHI, DNI, GHI) alebo oblačnosť (voliteľné, ale zlepšuje modely)

**Dodatočné užitočné atribúty:**
- Deň v týždni, Mesiac, Sezóna
- Údaje o produkcii fotovoltaiky (PV)
- Údaje o produkcii vetra
- Príznaky sviatkov/Špeciálnych udalostí

### Kde nájsť podobné datasety

#### 1. **Google Dataset Search** (Bezplatne)
   - **URL:** https://datasetsearch.research.google.com/
   - **Výrazy na vyhľadávanie:** "electricity demand", "electricity load", "power consumption time series"
   - **Výhoda:** Indexuje datasety z viacerých zdrojov
   - **Ako:** Zadajte vyššie uvedené výrazy, filtrujte podľa "Time Series" a "Weather Data"

#### 2. **Kaggle Datasets** (Bezplatne)
   - **URL:** https://www.kaggle.com/datasets/
   - **Relevantné datasety:**
     - "Electricity Load Forecasting"
     - "Smart Meter Data"
     - "Energy Consumption"
     - "Household Power Consumption"
   - **Výhoda:** Overené komunitou, často obsahujú EDA notebooky
   - **Odporúčané:** [Household Electric Power Consumption](https://www.kaggle.com/datasets/uciml/electric-power-consumption-data) alebo [Smart Meter Data](https://www.kaggle.com/datasets/tamaramerriman/smarter-meter-data)

#### 3. **UCI Machine Learning Repository** (Bezplatne)
   - **URL:** https://archive.ics.uci.edu/
   - **Vyhľadávanie:** "electricity" alebo "energy"
   - **Zbierky:**
     - Datasety na predikciu energetickej účinnosti
     - Datasety spotreby energie budov
     - Údaje z inteligentných meračov

#### 4. **Kaggle Competitions** (Bezplatne & Reálne)
   - Historické súťaže o prognozu zaťaženia často poskytujú datasety
   - Príklad: "Load Forecasting" alebo "Electricity Price Forecasting"

#### 5. **Operátori elektrizačných sústav** (Bezplatne)
   - Mnoho krajín zverejňuje historické údaje o spotrebe elektrickej energie
   - **Veľká Británia:** National Grid ESO ([Demand Forecasting](https://data.nationalgrideso.com/))
   - **USA:** EIA ([U.S. Energy Information Administration](https://www.eia.gov/))
   - **Írsko:** SONI ([System Operator for Northern Ireland](https://www.soni.ltd.uk/))
   - **Slovensko:** SEPS ([Prevádzkovateľ elektrizačnej sústavy na Slovensku](https://www.seps.sk/en/))

#### 6. **OpenEI** (Bezplatne)
   - **URL:** https://openei.org/
   - Veľké úložisko energetických datasetov

#### 7. **zenodo / figshare** (Bezplatne - Výskumné datasety)
   - Často obsahujú recenzované energetické datasety
   - Vyhľadávajte "electricity demand" alebo "smart meter"

### Ako prispôsobiť svoj dataset

Ak má váš dataset **iné časové rozlíšenie** (napr. hodinové namiesto 5-minútového):

1. **Resampling v prípade potreby:** Použite pandas `resample()`
   ```python
   df = df.resample('5min').interpolate()  # Ak máte hodinové údaje
   ```

2. **Úprava inžinierstva atribútov:**
   - Namiesto 288 oneskorení (24h × 5min), použite:
     - Hodinové údaje: 24 oneskorení (1-24 hodín)
     - 30-minútové údaje: 48 oneskorení (24 hodín)

3. **Upravte časové horizonty v skriptoch:**
   ```python
   # Ak máte hodinové údaje, zmeňte:
   HORIZONS = {
       'day': 24,      # 24 hodín
       'week': 168,    # 7 dní
       'month': 720    # 30 dní
   }
   ```

### Príprava údajov

Bez ohľadu na zdroj pripravte vaše údaje s:
- **Oneskorené atribúty:** Historická spotreba (t-1, t-2, ..., t-24, t-48, t-288)
- **Valcové štatistiky:** Priemer a štandardná odchýlka za rôzne okna (1h, 2h, 4h, 1-deň)
- **Cyklické kódovanie:** Sin/Cos transformácia pre hodinu dňa a mesiac
- **Časové rozdelenie:** Train/Test delenie zachováva chronologický poriadok (bez úniku údajov)
- **Škálovanie:** StandardScaler alebo MinMaxScaler podľa typu modelu

Skripty `01_eda_and_preparation.py` a `02_train_test_split.py` to zvládnu automaticky, keď sú údaje umiestnené v `prediction/data/`

---

## Ako spustiť experimenty

### Sprievodca krok za krokom

#### **Krok 1: Prieskumná analýza údajov** (10-15 minút)
```bash
python prediction/01_eda_and_preparation.py
```
**Výstup:** Vizualizačné grafy uložené do `prediction/plots/`
- Prehľad časových radov
- Matica korelácie
- Sezónnymi vzory
- Analýza odľahlých hodnôt

#### **Krok 2: Príprava údajov a inžinierstvo atribútov** (5 minút)
```bash
python prediction/02_train_test_split.py
```
**Výstup:** Pripravené datasety pre každý časový horizont uložené do `prediction/data/`
- Atribúty sú štandardizované podľa horizontu
- Train/test delenie zachováva chronológiu (bez úniku údajov)

#### **Krok 3: Trénovanie ARIMA modelu** (30-60 minút)
```bash
python prediction/03_arima_model.py
```
**Výstup:** ARIMA predikcie a metriky uložené do `prediction/results/arima/`
- Klasický štatistický prístup pre časové rady

#### **Krok 4: Trénovanie SVM modelu** (15-30 minút)
```bash
python prediction/04_svm_model.py
```
**Výstup:** SVM predikcie uložené do `prediction/results/svm/`
- Support Vector Regression s RBF jadrom
- Hyperparametre ladené pomocou grid search

#### **Krok 5: Trénovanie CNN modelu** (20-40 minút)
```bash
python prediction/05_cnn_model.py
```
**Výstup:** CNN predikcie a história trénovania uložené do `prediction/results/cnn/`
- 1D Konvolučná neurónová sieť
- Aplikované early stopping a zníženie miery učenia

#### **Krok 6: Trénovanie LightGBM modelu** (10-20 minút)
```bash
python prediction/07_lightgbm.py
```
**Výstup:** LightGBM predikcie a grafy dôležitosti atribútov uložené do `prediction/results/lightgbm/`
- Gradient Boosting rozhodovacích stromov
- Zahrnutá analýza dôležitosti atribútov

#### **Krok 7: Trénovanie LSTM a TCN modelov** (40-120 minút)
```bash
python prediction/08_lstm_tcn.py
```
**Výstup:** LSTM a TCN predikcie s históriou trénovania uložené do `prediction/results/lstm_tcn/`
- LSTM: Long Short-Term Memory siete
- TCN: Temporal Convolutional Networks s dilatovanými konvolúciami

#### **Krok 8: Trénovanie pokročilých modelov** (Voliteľne, 30-60 minút)
```bash
python prediction/09_advanced_models.py
```
**Výstup:** N-BEATS a TFT predikcie uložené do `prediction/results/advanced/`
- N-BEATS: State-of-the-art Neural Basis Expansion Analysis
- Zjednodušený Temporal Fusion Transformer

#### **Krok 9: Porovnanie modelov** (5-10 minút)
```bash
python prediction/06_comparison.py
```
**Výstup:** Komplexné grafy porovnania a metriky uložené do `prediction/results/comparison/`
- Porovnávacia vizualizácia všetkých modelov
- Tabuľka metrík výkonu (RMSE, MAE, MAPE)
- Žebríček najlepších modelov
- Presadenie predpovedí pre každý horizont

#### **Krok 10: Ensemble metódy** (5-10 minút)
```bash
python prediction/10_hybrid_stacking.py
```
**Výstup:** Výsledky ensemble uložené do `prediction/results/comparison/`
- Priemeraný ensemble
- Medián ensemble
- Ridge regression stacking (meta-learner)
- Vizualizácia výkonu ensemble

### Rýchle spustenie všetkých modelov

```bash
cd prediction
python 01_eda_and_preparation.py
python 02_train_test_split.py
python 03_arima_model.py
python 04_svm_model.py
python 05_cnn_model.py
python 07_lightgbm.py
python 08_lstm_tcn.py
python 09_advanced_models.py
python 06_comparison.py
python 10_hybrid_stacking.py
```

**Odhadovaný celkový čas spustenia:** 3-5 hodín (podľa dostupnosti GPU a špecifikácií systému)

---

## Prehľad modelov

| Model | Typ | Framework | Kľúčové vlastnosti |
|-------|-----|-----------|------------------|
| **ARIMA** | Klasická štatistika | Statsmodels | Testovanie stacionarity, ARIMA(2,1,2) |
| **SVM** | Strojové učenie | Scikit-learn | RBF jadro, ladenie hyperparametrov |
| **CNN** | Hlboké učenie | Keras/TensorFlow | 1D konvolúcie, batch normalizácia |
| **LightGBM** | Gradient Boosting | LightGBM | Dôležitosť atribútov, early stopping |
| **LSTM** | RNN | Keras/TensorFlow | Obojsmerné, mechanizmy pozornosti |
| **TCN** | Temporálny Conv | Keras/TensorFlow | Dilatované konvolúcie, kauzálne padding |
| **N-BEATS** | State-of-the-Art | NeuralForecast | Architektúra založená na stackoch |
| **TFT** | Na základe pozornosti | Keras/TensorFlow | Temporal Fusion s multi-head pozornosťou |
| **Ensemble** | Meta-učenie | Scikit-learn | Priemer, Medián, Ridge Stacking |

---

## Metriky výkonu

Tri časové horizonty sú vyhodnocované:
- **Deň:** 288 časových krokov (24 hodín)
- **Týždeň:** 2,016 časových krokov (7 dní)
- **Mesiac:** 8,640 časových krokov (30 dní)

### Metriky hodnotenia

- **RMSE** (Root Mean Square Error) - meria priemernú veľkosť chyby predikcie
- **MAE** (Mean Absolute Error) - priemerná absolútna odchýlka
- **MAPE** (Mean Absolute Percentage Error) - percentuálna chyba vztiahnutá na skutočné hodnoty

$$RMSE = \sqrt{\frac{1}{n}\sum_{i=1}^{n}(y_i - \hat{y}_i)^2}$$

$$MAE = \frac{1}{n}\sum_{i=1}^{n}|y_i - \hat{y}_i|$$

$$MAPE = \frac{100}{n}\sum_{i=1}^{n}\left|\frac{y_i - \hat{y}_i}{y_i}\right|$$

---

## Výsledky a výstupy

Všetky výsledky sú uložené do `prediction/results/` s nasledujúcou štruktúrou:

### 1. Grafy EDA (prediction/plots/)
- Vizualizácia časových radov v rôznych mierach
- Heatmapa korelácie atribútov
- Rozklad sezónnosti (hodinové, denné, týždenné, mesačné vzory)
- Box ploty na porovnanie pracovných dní a víkendov

### 2. Modely-špecifické výsledky (prediction/results/[model_name]/)
- CSV súbory predpovedí pre každý horizont
- Vizualizačné grafy porovnávajúce predikcie so skutočnosťou
- Grafy histórie trénovania (pre neurónové siete)
- Grafy dôležitosti atribútov (pre modely založené na stromoch)

### 3. Výsledky porovnania (prediction/results/comparison/)
- **ensemble_visualization_day.png, week.png, month.png** - Porovnania ensemble modelov
- **predictions_day.png, week.png, month.png** - Všetky modely naložené na skutočné údaje
- **metrics_comparison.png** - Stĺpcový graf RMSE/MAE/MAPE modelov
- **model_ranking.png** - Vizualizácia žebríčka
- **model_comparison.csv** - Tabuľka podrobných metrík
- **ensemble_summary.csv** - Zhrnutie výkonu ensemble metód

---

## Reprodukovanie výsledkov

Pre reprodukovanie presných výsledkov a experimentov:

### 1. Zabezpečte rovnaké prostredie
```bash
conda create -n energy-forecast python=3.10 -y
conda activate energy-forecast
pip install -r requirements.txt
```

### 2. Nastavte náhodné semienka
Všetky skripty obsahujú `np.random.seed(42)` a `tf.random.set_seed(42)` pre reprodukovateľnosť.

### 3. Spustite postupne kroky
Postupujte podľa vyššie uvedeného sprievodcu v poradí (kroky 1-10).

### 4. Overte výsledky
Skontrolujte, že:
- Všetky grafy sú generované v `prediction/plots/` a `prediction/results/`
- CSV súbory obsahujú predikcie pre všetky tri horizonty (deň, týždeň, mesiac)
- Tabuľka metrík zobrazuje hodnoty RMSE/MAE/MAPE

### Poznámky
- **Prvé spustenie:** Bude pomalšie, pretože modely sa trénujú prvýkrát
- **Pamäť:** Modely hlbokého učenia vyžadujú značné množstvo RAM; znížite veľkosť dávky pri chybách OOM
- **GPU:** GPU s podporou CUDA urýchľuje trénovanie výrazne (5-10x rýchlejšie)
- **Údaje:** Dataset musí ostať v zložke `prediction/data/` na správnu funkčnosť skriptov

---

## Riešenie problémov

### Problém: ModuleNotFoundError

**Riešenie:** Znova nainštalujte závislosti
```bash
pip install --upgrade pandas numpy matplotlib seaborn scikit-learn tensorflow lightgbm statsmodels
```

### Problém: CUDA/GPU sa nezistilo

**Riešenie:** TensorFlow sa automaticky vráti na CPU. Pre podporu GPU:
```bash
pip install tensorflow[and-cuda]  # Vyžaduje CUDA 11.8 a cuDNN 8.6
```

### Problém: Nedostatok pamäte (OOM)

**Riešenie:** Znížte veľkosť dávky alebo dĺžku sekvencie pri trénovaní modelu:
```python
# V skripte modelu znížte:
BATCH_SIZE = 16  # z 32
SEQ_LENGTH = 12  # z 24
```

### Problém: Dataset nenájdený

**Riešenie:** Zabezpečte, že existuje zložka `prediction/data/` s riadne pripravenými CSV súbormi. Najskôr spustite `02_train_test_split.py`.

---

## Licencia

Tento projekt je súčasťou diplomovej práce. Použitie na vzdelávacie a výskumné účely je odporúčané.

---