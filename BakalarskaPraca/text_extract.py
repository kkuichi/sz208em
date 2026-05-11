import os
from lxml import etree

# Cesta k priečinku so súbormi XML
xml_folder_path = './xmls/'
# Cesta k súboru XSLT
xslt_path = './extract_text.xslt'
# Cesta k priečinku na uloženie konečného textového súboru
text_folder_path = './text'
# Názov konečného textového súboru
text_file_name = 'extracted_text.txt'

# Vytvorenie priečinka na uloženie textového súboru, ak neexistuje
os.makedirs(text_folder_path, exist_ok=True)

# Definujte cestu na uloženie konečného textového súboru
text_file_path = os.path.join(text_folder_path, text_file_name)

# Načítanie šablóny XSLT
with open(xslt_path, 'rb') as f:
    xslt_root = etree.XML(f.read())
xslt_transform = etree.XSLT(xslt_root)

# Funkcia na spracovanie jedného súboru XML
def process_xml_file(xml_file_path):
    with open(xml_file_path, 'rb') as f:
        xml_root = etree.XML(f.read())
    result_tree = xslt_transform(xml_root)
    return str(result_tree)

# Otvorenie konečného textového súboru na nahrávanie
with open(text_file_path, 'w', encoding='utf-8') as text_file:
    # Spracovať všetky súbory XML v priečinku
    for filename in os.listdir(xml_folder_path):
        if filename.endswith('.xml'):
            xml_file_path = os.path.join(xml_folder_path, filename)
            extracted_text = process_xml_file(xml_file_path)
            # Zapísať extrahovaný text do súboru
            text_file.write(extracted_text)
            text_file.write('\n')  # Pridanie oddeľovača medzi obsah rôznych súborov
            print(f'Extracted text from {filename} added to {text_file_path}')

print('Processing completed.')
