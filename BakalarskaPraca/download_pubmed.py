import requests
from bs4 import BeautifulSoup
import os

# Adresa URL webovej stránky
url = 'https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/'

# Získanie obsahu HTML stránky
response = requests.get(url)
soup = BeautifulSoup(response.content, 'html.parser')

# Vyhľadajte všetky odkazy na súbory s príponou .gz
links = soup.find_all('a', href=True)
gz_links = [url + link['href'] for link in links if link['href'].endswith('.gz')]

# Vytvorenie priečinka na ukladanie súborov, ak neexistuje
os.makedirs('downloads', exist_ok=True)

# Stiahneme prvých 5 súborov
for i in range(min(5, len(gz_links))):
    gz_url = gz_links[i]
    file_name = os.path.join('downloads', gz_url.split('/')[-1])
    print(f'Starting download: {gz_url}')
    with requests.get(gz_url, stream=True) as r:
        r.raise_for_status()
        with open(file_name, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    print(f'Finished download: {file_name}')

print('All downloads completed.')