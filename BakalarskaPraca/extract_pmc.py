import json
import time
import requests
from Bio import Entrez
import xml.etree.ElementTree as ET
import re
import html
from concurrent.futures import ThreadPoolExecutor, as_completed
from time import sleep
from random import uniform

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

# Funkcia na získanie detailov článku pre jeden PMCID s opakovanými pokusmi
def fetch_article_detail(pmcid):
    retries = 5
    delay = 3
    for _ in range(retries):
        try:
            handle = Entrez.efetch(db="pmc", id=pmcid, rettype="full", retmode="xml")
            records = handle.read()
            root = ET.fromstring(records)
            
            pmid = None
            abstract = None
            full_text = None

            # Extrahovanie PMID
            pmid_elements = root.findall(".//article-id[@pub-id-type='pmid']")
            if pmid_elements:
                pmid = pmid_elements[0].text

            # Extrahovanie abstraktu
            abstract_elements = root.findall(".//abstract")
            if abstract_elements:
                abstract = ET.tostring(abstract_elements[0], encoding='unicode', method='text')
                abstract = clean_text(abstract)

            # Extrahovanie plného textu
            body_elements = root.findall(".//body")
            if body_elements:
                full_text = ET.tostring(body_elements[0], encoding='unicode', method='text')
                full_text = clean_text(full_text)

            return {
                "pmid": pmid,
                "pmcid": pmcid,
                "abstract": abstract,
                "full_text": full_text
            }
        except Exception as e:
            print(f"Chyba pri získavaní detailov pre PMCID {pmcid}: {e}. Opakovaný pokus...")
            sleep(delay)
            delay = min(delay * 2, 60)  # Exponenciálne oneskorenie, maximálne oneskorenie 60 sekúnd
    print(f"Nepodarilo sa získať detaily pre PMCID {pmcid} po {retries} pokusoch")
    return None

# Funkcia na súbežné získanie detailov článkov
def fetch_article_details_concurrent(pmcid_list):
    details = []
    with ThreadPoolExecutor(max_workers=10) as executor:  # Zvýšený počet vlákien
        futures = {executor.submit(fetch_article_detail, pmcid): pmcid for pmcid in pmcid_list}
        for future in as_completed(futures):
            result = future.result()
            if result:
                details.append(result)
                # Spánok na dodržanie obmedzenia rýchlosti
                sleep(uniform(0.2, 0.5))  # Znížený interval náhodného spánku
    return details

def main():
    num_articles = int(input("Zadajte počet článkov na stiahnutie: "))

    # Spustenie časovača
    start_time = time.time()

    # Vyhľadávanie článkov a získanie PMCIDs
    Entrez.email = "serghei.zabirchenko@student.tuke.sk"
    search_handle = Entrez.esearch(db="pmc", term="cancer", retmax=num_articles)
    search_results = Entrez.read(search_handle)
    pmcid_list = search_results['IdList']

    # Súbežné získanie detailov článkov
    articles = fetch_article_details_concurrent(pmcid_list)

    # Uloženie do JSON súboru
    with open('articles_pmc.json', 'w') as f:
        json.dump(articles, f, indent=4)

    # Zastavenie časovača
    end_time = time.time()
    print(f"Čas potrebný na stiahnutie {num_articles} článkov: {end_time - start_time} sekúnd")

if __name__ == "__main__":
    main()
