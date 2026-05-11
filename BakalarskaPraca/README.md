
# Systémová príručka

Systémová príručka na popis a spustenie všetkých kódov v tomto úložisku

# Stiahnutie súborov z PubMed Baseline (`download_pubmed.py`)

Tento skript slúži na stiahnutie súborov z databázy PubMed Baseline a ich uloženie do lokálneho priečinka. Používa knižnice `requests` a `BeautifulSoup` na získanie a analýzu odkazov na súbory. Stiahne prvých 5 súborov s príponou `.gz` z danej URL adresy.

### Požiadavky

Na spustenie tohto skriptu je potrebné mať nainštalované nasledujúce knižnice:
- `requests`
- `beautifulsoup4`

Môžete ich nainštalovať pomocou `pip`:

```sh
pip install requests beautifulsoup4
```

### Použitie
1. Spustite skript pomocou príkazu:
```sh
python download_pubmed.py
```
2. Skript stiahne všetky súbory s príponou `.gz` z URL adresy `https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/` a uloží ich do miestneho priečinka.

### Popis kódu
- Načítanie obsahu HTML stránky pomocou `requests` a `BeautifulSoup`.
- Vyhľadanie všetkých odkazov na súbory s príponou `.gz`.
- Vytvorenie priečinka na ukladanie súborov, ak neexistuje.
- Stiahnutie každého súboru a uloženie do miestneho priečinka.

# Extrakcia súborov z .gz archívov (`extract_gz.py`)

Tento skript slúži na extrakciu súborov z .gz archívov a ich uloženie do samostatného priečinka. Používa knižnice `os`, `gzip` a `shutil` na spracovanie a rozbalenie súborov.

### Požiadavky

Na spustenie tohto skriptu nie sú potrebné žiadne dodatočné knižnice, pretože používa štandardné knižnice Pythonu.

### Použitie

1. Uistite sa, že máte priečinok `downloads`, ktorý obsahuje .gz súbory, ktoré chcete extrahovať.
2. Spustite skript príkazom:

```sh
python extract_gz.py
```

### Funkcie skriptu

1. #### Definovanie ciest:
Skript najprv definuje cestu k priečinku s .gz súbormi (`gz_folder_path`) a cestu k priečinku, kam budú uložené extrahované súbory (`xmls_folder_path`).

2. #### Vytvorenie priečinka pre extrahované súbory:
Skript vytvorí priečinok `xmls`, ak ešte neexistuje, pomocou `os.makedirs`.

3. #### Funkcia na rozbalenie jedného súboru .gz:
Funkcia `decompress_gz_file` rozbalí daný .gz súbor a uloží výsledok do určeného priečinka.

4. #### Spracovanie všetkých .gz súborov v priečinku:
Skript prejde všetky súbory v priečinku `downloads` a ak nájde súbor s príponou `.gz`, použije funkciu `decompress_gz_file` na jeho rozbalenie.

5. #### Výpis dokončenia dekompresie:
Po dokončení dekompresie všetkých súborov skript vypíše správu o úspešnom dokončení.

# Extrakcia textu zo súborov XML pomocou XSLT (`text_extract.py` a `extract_text.xslt`)

Tento skript slúži na extrakciu textu zo súborov XML pomocou XSLT transformácie a uloženie výsledného textu do jediného súboru. Používa knižnicu `lxml` na prácu s XML a XSLT.

### Požiadavky

Na spustenie tohto skriptu je potrebné mať nainštalovanú knižnicu `lxml`. Môžete ju nainštalovať pomocou `pip`:

```sh
pip install lxml
```

### Použitie

1. Uistite sa, že máte priečinok `xmls`, ktorý obsahuje XML súbory, ktoré chcete spracovať.
2. Uistite sa, že máte súbor `extract_text.xslt` v rovnakom priečinku ako skript.
3. Spustite skript príkazom:
```sh
python text_extract.py
```

### Funkcie skriptu

1. #### Definovanie ciest:
Skript najprv definuje cesty k priečinku so XML súbormi (`xml_folder_path`), k súboru XSLT (`xslt_path`) a k priečinku, kam bude uložený konečný textový súbor (`text_folder_path`).

2. #### Vytvorenie priečinka pre konečný textový súbor:
Skript vytvorí priečinok `text`, ak ešte neexistuje, pomocou `os.makedirs`.

3. #### Načítanie šablóny XSLT:
Skript načíta a skompiluje XSLT šablónu pomocou `etree.XSLT`.

4. #### Funkcia na spracovanie jedného súboru XML:
Funkcia `process_xml_file` načíta XML súbor, aplikuje naň XSLT transformáciu a vráti výsledný text.

5. #### Spracovanie všetkých XML súborov v priečinku:
Skript prejde všetky súbory v priečinku `xmls` a ak nájde súbor s príponou `.xml`, použije funkciu `process_xml_file` na jeho spracovanie a výsledný text uloží do konečného textového súboru.

6. #### Výpis dokončenia spracovania:
Po dokončení spracovania všetkých súborov skript vypíše správu o úspešnom dokončení.

# Extrakcia údajov z DrugBank XML (`drugbank.py`)

Tento skript slúži na extrakciu relevantných informácií o liekoch z XML súboru databázy DrugBank a ich uloženie do formátu JSON. Skript používa knižnicu `lxml` na analýzu XML súboru a knižnicu `re` na čistenie textu.

### Požiadavky

Na spustenie tohto skriptu je potrebné mať nainštalovanú knižnicu `lxml`. Môžete ju nainštalovať pomocou `pip`:

```sh
pip install lxml
```

### Použitie
1. Uistite sa, že máte XML súbor s názvom full database.xml v rovnakom priečinku ako skript.
2. Spustite skript príkazom:
```sh
python drugbank.py
```

### Funkcie skriptu
1. #### Definovanie ciest:
Skript najprv definuje cesty k vstupnému XML súboru (`xml_file_path`) a k výstupnému JSON súboru (`json_file_path`).

2. #### Čistenie textu:
Funkcia `clean_text` odstráni nežiaduce znaky, odkazy, HTML entity a značky z textu.

3. #### Analýza relevantných prvkov:
Funkcia `parse_drug_element` analyzuje prvok `drug` a extrahuje z neho relevantné informácie na základe zoznamu `relevant_tags`.

4. #### Bezpečná analýza XML:
Funkcia `safe_iterparse` iteratívne analyzuje XML súbor a ignoruje chyby syntaxe.

5. #### Spracovanie XML súboru:
Skript iteratívne analyzuje XML súbor, spracováva každý prvok `drug` a extrahované informácie ukladá do zoznamu `drug_data`.

6. #### Uloženie údajov do JSON:
Skript uloží extrahované údaje do JSON súboru s názvom `drugbank_corpus.json`.

7. #### Výpis výsledkov:
Skript vypíše správu o úspešnom vytvorení JSON súboru, celkový počet spracovaných článkov a čas spracovania.

# Extrakcia údajov z PubMed (`extract_pm.py`)

Tento skript slúži na extrakciu údajov z databázy PubMed a ich uloženie do formátu JSON. Používa knižnicu `Bio.Entrez` na získanie údajov z PubMed a knižnicu `re` na čistenie textu.

### Požiadavky

Na spustenie tohto skriptu je potrebné mať nainštalovanú knižnicu `biopython`. Môžete ju nainštalovať pomocou `pip`:

```sh
pip install biopython
```
### Použitie

1. Uistite sa, že máte nainštalovanú knižnicu `biopython`.
2. Nahraďte emailovú adresu v skripte svojou emailovou adresou.
3. Spustite skript príkazom:
```sh
python extract_pm.py
```
4. Po zadaní požadovaného počtu článkov skript stiahne údaje z PubMed a uloží ich do súboru `articles_pm.json`.

### Funkcie skriptu

1. #### Čistenie textu:
Funkcia `clean_text` odstráni HTML entity, citácie, nadbytočné medzery a iné nežiaduce znaky z textu.

2. #### Získanie údajov z PubMed:
Funkcia `fetch_pubmed_data` vykonáva vyhľadávanie článkov na PubMed podľa zadaného dotazu a získava ich abstrakty a ďalšie relevantné údaje.

3. #### Uloženie údajov do JSON:
Funkcia `save_to_json` uloží získané údaje do súboru vo formáte JSON.

# Extrakcia údajov z PMC (`extract_pmc.py`)

Tento skript slúži na extrakciu údajov z databázy PubMed Central (PMC) a ich uloženie do formátu JSON. Používa knižnice `requests`, `Bio.Entrez`, `xml.etree.ElementTree`, `re`, `html` a `concurrent.futures` na získanie a spracovanie údajov.

### Požiadavky

Na spustenie tohto skriptu je potrebné mať nainštalované knižnice `biopython` a `requests`. Môžete ich nainštalovať pomocou `pip`:

```sh
pip install biopython requests
```

### Použitie
1. Uistite sa, že máte nainštalované potrebné knižnice.
2. Nahraďte emailovú adresu v skripte svojou emailovou adresou.
3. Spustite skript príkazom:
```sh
python extract_pmc.py
```
4. Po zadaní požadovaného počtu článkov skript stiahne údaje z PMC a uloží ich do súboru `articles_pmc.json`.

### Funkcie skriptu

1. #### Čistenie textu:
Funkcia `clean_text` odstráni HTML entity, citácie, nadbytočné medzery a iné nežiaduce znaky z textu.

2. #### Získanie detailov článku:
Funkcia `fetch_article_detail` získava detaily článku pre daný PMCID s opakovanými pokusmi v prípade chyby.

3. #### Súbežné získanie detailov článkov:
Funkcia `fetch_article_details_concurrent` získava detaily článkov súbežne pomocou `ThreadPoolExecutor`.

4. #### Hlavná funkcia:
Funkcia `main` inicializuje vyhľadávanie článkov, získava PMCIDs, volá súbežnú funkciu na získanie detailov a ukladá výsledky do JSON súboru.

# Zlúčenie článkov z JSON súborov (`merger.py`)

Tento skript slúži na zlúčenie článkov z dvoch JSON súborov, odstránenie duplicitných článkov a uloženie výsledného zoznamu do nového JSON súboru.

### Požiadavky

Na spustenie tohto skriptu nie sú potrebné žiadne dodatočné knižnice, pretože používa štandardnú knižnicu `json` Pythonu.

### Použitie

1. Uistite sa, že máte JSON súbory `articles_pmc.json` a `articles_pm.json` v rovnakom priečinku ako skript.
2. Spustite skript príkazom:
```sh
python merger.py
```
3. Skript zlúči články z oboch súborov, odstráni duplicity a uloží výsledok do súboru `merged_articles.json`.

### Funkcie skriptu

1. #### Načítanie JSON súboru:
Funkcia `load_json` načíta a vráti obsah JSON súboru.

2. #### Uloženie do JSON súboru:
Funkcia `save_json` uloží dáta do JSON súboru s pekným formátovaním.

3. #### Odstránenie duplicitných článkov:
Funkcia `remove_duplicates` odstráni články z `articles_pm`, ktoré sa už nachádzajú v `articles_pmc`, a vráti filtrovaný zoznam a počet odstránených článkov.

4. #### Zlúčenie článkov:
Funkcia `merge_articles` zlúči dva zoznamy článkov do jedného.

5. #### Hlavná funkcia:
Funkcia `main` načíta JSON súbory, odstráni duplicity, zlúči články a uloží výsledok do nového JSON súboru. Na konci vytlačí počet odstránených článkov.

# Zlúčenie dát z JSON súborov (`merger_drug.py`)

Tento skript slúži na zlúčenie dát z dvoch JSON súborov a uloženie výsledného zoznamu do nového JSON súboru.

### Požiadavky

Na spustenie tohto skriptu nie sú potrebné žiadne dodatočné knižnice, pretože používa štandardnú knižnicu `json` Pythonu.

### Použitie

1. Uistite sa, že máte JSON súbory `merged_articles.json` a `drugbank_corpus.json` v rovnakom priečinku ako skript.
2. Spustite skript príkazom:
```sh
python merger_drug.py
```
3. Skript zlúči dáta z oboch súborov a uloží výsledok do súboru `final_merged_data.json`.

### Funkcie skriptu

1. #### Načítanie JSON súboru:
Funkcia `load_json` načíta a vráti obsah JSON súboru.

2. #### Uloženie do JSON súboru:
Funkcia `save_json` uloží dáta do JSON súboru s pekným formátovaním.

3. #### Zlúčenie dát:
Funkcia `merge_files` načíta dáta z dvoch JSON súborov, zlúči ich a vráti výsledok.

4. #### Hlavná funkcia:
Funkcia `main` načíta JSON súbory, zlúči dáta a uloží výsledok do nového JSON súboru. Na konci vytlačí správu o úspešnom uloženom zlúčení.