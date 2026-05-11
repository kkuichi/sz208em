import json

def load_json(filename):
    with open(filename, 'r') as f:
        return json.load(f)

def save_json(filename, data):
    with open(filename, 'w') as f:
        json.dump(data, f, indent=4)

def remove_duplicates(articles_pmc, articles_pm):
    # Vytvorenie množiny pmid z articles_pmc
    pmid_set_pmc = {article['pmid'] for article in articles_pmc}
    # Filtrovanie articles_pm, aby sa odstránili články, ktoré majú pmid v pmid_set_pmc
    filtered_pm = [article for article in articles_pm if article['id'] not in pmid_set_pmc]
    # Počet odstránených článkov
    deleted_count = len(articles_pm) - len(filtered_pm)
    return filtered_pm, deleted_count

def merge_articles(articles_pmc, articles_pm):
    # Zlúčenie dvoch zoznamov článkov
    return articles_pmc + articles_pm

def main():
    # Načítanie JSON súborov
    articles_pmc = load_json('articles_pmc.json')
    articles_pm = load_json('articles_pm.json')
    
    # Odstránenie duplicitných článkov a počítanie odstránených článkov
    filtered_pm, deleted_count = remove_duplicates(articles_pmc, articles_pm)
    
    # Zlúčenie článkov
    merged_articles = merge_articles(articles_pmc, filtered_pm)
    
    # Uloženie zlúčených článkov do nového JSON súboru
    save_json('merged_articles.json', merged_articles)
    
    # Tlač počtu odstránených článkov
    print(f"{deleted_count} článkov bolo odstránených z 'articles_pm.json'")
    print("Zlúčené články uložené do 'merged_articles.json'")

if __name__ == "__main__":
    main()