from Bio import Entrez
import json
import time
import html
import re

# Funkcia na vyčistenie textu
def clean_text(text):
    if not text:
        return None
    # Dekódovanie HTML entít
    text = html.unescape(text)
    # Odstránenie citácií ako [1], [2], atď.
    text = re.sub(r'\[\d+\]', '', text)
    # Zabezpečenie, že kľúčové slová ako ABSTRACT a INTRODUCTION sú nasledované medzerou
    keywords = ["ABSTRACT", "INTRODUCTION", "CONCLUSIONS", "METHODS", "RESULTS", "DISCUSSION"]
    for keyword in keywords:
        text = re.sub(rf'({keyword})([^\s])', rf'\1 \2', text)
    # Nahradenie viacerých nových riadkov a medzier jednou medzerou
    text = re.sub(r'\s+', ' ', text)
    # Odstránenie iných nežiaducich znakov
    text = re.sub(r'[^\x00-\x7F]+', '', text)
    # Odstránenie úvodných a záverečných medzier
    text = text.strip()
    return text

def fetch_pubmed_data(target_count, email):
    Entrez.email = email
    
    # Začiatok merania času
    start_time = time.time()
    
    articles = []
    fetched_count = 0
    batch_size = target_count  # Počiatočná veľkosť dávky
    
    while len(articles) < target_count:
        # Vyhľadávanie článkov podľa všeobecného dotazu s nárastom počiatočnej veľkosti dávky
        search_handle = Entrez.esearch(db="pubmed", term="all[sb]", retstart=fetched_count, retmax=batch_size)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        id_list = search_results["IdList"]
        
        # Získanie informácií o článkoch podľa ich ID
        fetch_handle = Entrez.efetch(db="pubmed", id=",".join(id_list), retmode="xml")
        fetch_data = Entrez.read(fetch_handle)
        fetch_handle.close()
        
        # Parsovanie dát
        for article in fetch_data['PubmedArticle']:
            abstract_text = article['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', '')
            
            # Join list of strings into a single string if AbstractText is a list
            if isinstance(abstract_text, list):
                abstract_text = ' '.join(abstract_text)
            
            article_data = {
                "id": article['MedlineCitation']['PMID'],
                "title": article['MedlineCitation']['Article']['ArticleTitle'],
                "abstract": clean_text(abstract_text)
            }
            
            # Kontrola, či existuje abstrakt
            if not article_data['abstract']:
                continue
            
            articles.append(article_data)
            
            # Prerušenie, ak je dosiahnuté požadované množstvo článkov
            if len(articles) >= target_count:
                break
        
        # Aktualizácia počtu získaných článkov
        fetched_count += len(id_list)
        
        # Zvýšenie veľkosti dávky pre ďalší dotaz
        batch_size = target_count - len(articles)
    
    # Koniec merania času
    end_time = time.time()
    elapsed_time = end_time - start_time
    
    print(f"Čas trvania: {elapsed_time:.2f} sekúnd")
    print(f"Počet získaných článkov: {len(articles)}")
    
    return articles

def save_to_json(data, filename):
    with open(filename, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=4)

# Počet článkov pre získanie zadaný užívateľom
target_count = int(input("Zadajte počet článkov pre získanie: "))
email = "serghei.zabirchenko@student.tuke.sk"  # Nahraďte svojím emailom
data = fetch_pubmed_data(target_count, email)
save_to_json(data, "articles_pm.json")

print("Údaje uložené do súboru articles_pm.json")