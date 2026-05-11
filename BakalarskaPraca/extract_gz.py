import os
import gzip
import shutil

# Cesta k priečinku archívu .gz
gz_folder_path = './downloads'
# Nová cesta k priečinku pre rozbalené súbory
xmls_folder_path = './xmls'

# Vytvorenie priečinka pre rozbalené súbory, ak neexistuje
os.makedirs(xmls_folder_path, exist_ok=True)

# Funkcia na rozbalenie jedného súboru .gz
def decompress_gz_file(gz_file_path, output_folder):
    # Získanie názvu súboru bez prípony .gz
    file_name = os.path.basename(gz_file_path).replace('.gz', '')
    # Definujte cestu na uloženie rozbaleného súboru
    output_file_path = os.path.join(output_folder, file_name)
    
    # Rozbalenie súboru
    with gzip.open(gz_file_path, 'rb') as f_in:
        with open(output_file_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    print(f'Decompressed {gz_file_path} to {output_file_path}')

# Spracovať všetky súbory .gz v priečinku
for filename in os.listdir(gz_folder_path):
    if filename.endswith('.gz'):
        gz_file_path = os.path.join(gz_folder_path, filename)
        decompress_gz_file(gz_file_path, xmls_folder_path)

print('Decompression completed.')
