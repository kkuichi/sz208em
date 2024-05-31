from lxml import etree
import json
import time
import re

# Definovať cesty k súborom
xml_file_path = 'full database.xml'
json_file_path = 'drugbank_corpus.json'

# Nastaviť na držanie extrahovaných údajov
drug_data = []

# Funkcia na vyčistenie textu odstránením nežiaducich znakov a odkazov
def clean_text(text):
    # Odstrániť znaky nového riadku
    text = text.replace('\n', ' ')
    # Odstrániť odkazy ako [L41539]
    text = re.sub(r'\[\w+(,\w+)*\]', '', text)
    # Odstrániť viacnásobné medzery
    text = re.sub(r'\s+', ' ', text)
    # Odstrániť HTML entity (napr. &#13;)
    text = re.sub(r'&#\d+;', '', text)
    # Odstrániť HTML značky (napr. <sub>max</sub>)
    text = re.sub(r'<.*?>', '', text)
    # Odstrániť unicode znaky (napr. \u03b1-thrombin)
    text = re.sub(r'\\u[\da-fA-F]{4}', '', text)
    return text.strip()

# Funkcia na analýzu relevantných prvkov z elementu drug
def parse_drug_element(drug_elem):
    drug_info = {}
    for child in drug_elem:
        tag = etree.QName(child.tag).localname
        if tag in relevant_tags:
            drug_info[tag] = clean_text(child.text) if child.text else None
        elif tag == 'general-references':
            for article in child.findall('.//article'):
                title = article.find('title')
                if title is not None:
                    drug_info['title'] = clean_text(title.text)
        elif tag == 'articles':
            for article in child.findall('article'):
                title = article.find('title')
                if title is not None:
                    drug_info['title'] = clean_text(title.text)
    return drug_info

# Relevantné značky, ktoré sa majú zahrnúť do JSON
relevant_tags = {
    'name', 'description', 'pharmacodynamics', 'mechanism-of-action',
    'indication', 'synthesis-reference', 'toxicity', 'metabolism',
    'absorption', 'half-life', 'protein-binding', 'route-of-elimination',
    'volume-of-distribution', 'clearance', 'subclass'
}

# Funkcia na bezpečnú analýzu XML a ignorovanie chýb
def safe_iterparse(file_path):
    try:
        context = etree.iterparse(file_path, events=('start', 'end'), recover=True)
        return context
    except etree.XMLSyntaxError as e:
        print(f"XMLSyntaxError: {e}")
        return []

# Spustiť časovač
start_time = time.time()

# Iteratívne analyzovať XML súbor
context = safe_iterparse(xml_file_path)

# Inicializovať počítadlo
article_count = 0

# Spracovať každý prvok
if context:
    for event, elem in context:
        if event == 'end' and etree.QName(elem.tag).localname == 'drug':
            drug_info = parse_drug_element(elem)
            # Pridať drug_info iba ak obsahuje viac ako len názov
            if drug_info and any(key in drug_info for key in relevant_tags - {'name'}):
                drug_data.append(drug_info)
                article_count += 1
            elem.clear()  # Vyčistiť spracované prvky na uvoľnenie pamäte

# Uložiť extrahované údaje do JSON súboru
with open(json_file_path, 'w') as json_file:
    json.dump(drug_data, json_file, indent=4)

# Ukončiť časovač
end_time = time.time()

# Vytlačiť výsledky
print(f"JSON súbor s údajmi o liekoch bol úspešne vytvorený na {json_file_path}")
print(f"Celkový počet spracovaných článkov: {article_count}")
print(f"Čas spracovania: {end_time - start_time:.2f} sekúnd")