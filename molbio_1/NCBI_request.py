import requests
from bs4 import BeautifulSoup
import pandas as pd
import re

def fetch_ncbi_codes():
    url = "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi"
    response = requests.get(url)
    response.raise_for_status()
    soup = BeautifulSoup(response.text, 'html.parser')

    genetic_codes = {}
    all_codons = [f"{b1}{b2}{b3}" for b1 in "TCAG" for b2 in "TCAG" for b3 in "TCAG"]
    pre_tags = soup.find_all('pre')
    for pre in pre_tags:
        text = pre.get_text()
        if 'AAs' not in text or 'Base1' not in text:
            continue

        aa_match = re.search(r'AAs\s*=\s*([A-Z\*]+)', text)
        base1_match = re.search(r'Base1\s*=\s*([TCAG]+)', text)
        base2_match = re.search(r'Base2\s*=\s*([TCAG]+)', text)
        base3_match = re.search(r'Base3\s*=\s*([TCAG]+)', text)

        if not (aa_match and base1_match and base2_match and base3_match):
            continue

        aa_string = aa_match.group(1)
        base1 = base1_match.group(1)
        base2 = base2_match.group(1)
        base3 = base3_match.group(1)

        if not (len(aa_string) == 64 and len(base1) == 64 and len(base2) == 64 and len(base3) == 64):
            continue

        code_name = "Unknown"
        sibling = pre.find_previous('h2')
        if sibling:
            code_name = sibling.get_text(strip=True)
        else:
            pass
        codon_to_aa = {}
        for i in range(64):
            codon = base1[i] + base2[i] + base3[i]
            aa = aa_string[i]
            codon_to_aa[codon] = aa

        series = pd.Series(codon_to_aa, name=code_name)
        series = series.reindex(all_codons)
        genetic_codes[code_name] = series

    return genetic_codes

def save_codes_to_csv(genetic_codes, filename='genetic_codes.csv'):
    df = pd.DataFrame(genetic_codes)
    df.index.name = 'Codon'
    df.to_csv(filename)
    print(f"Сохранено {len(genetic_codes)} вариантов кода в {filename}")

if __name__ == "__main__":
    codes = fetch_ncbi_codes()
    if codes:
        save_codes_to_csv(codes)
        print("complete")
        df = pd.read_csv('genetic_codes.csv', index_col=0)
        print(df.head())
    else:
        print("err")