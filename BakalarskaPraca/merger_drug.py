import json

def load_json(filename):
    with open(filename, 'r') as f:
        return json.load(f)

def save_json(filename, data):
    with open(filename, 'w') as f:
        json.dump(data, f, indent=4)

def merge_files(file1, file2):
    # Načítanie dát zo súborov
    data1 = load_json(file1)
    data2 = load_json(file2)
    
    # Zlúčenie dát
    merged_data = data1 + data2
    
    return merged_data

def main():
    # Názvy súborov
    file1 = 'merged_articles.json'
    file2 = 'drugbank_corpus.json'
    
    # Zlúčenie súborov
    merged_data = merge_files(file1, file2)
    
    # Uloženie zlúčených dát do nového súboru
    save_json('final_merged_data.json', merged_data)
    
    print("Zlúčené dáta uložené do 'final_merged_data.json'")

if __name__ == "__main__":
    main()