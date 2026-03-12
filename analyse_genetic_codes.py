import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re

df = pd.read_csv('genetic_codes.csv', index_col=0)
print(f"Загружено {df.shape[1]} вариантов кода и {df.shape[0]} кодонов.")

if df.isnull().any().any():
    df = df.fillna('?')

conservative_codons = []
codon_values = {}

for codon in df.index:
    values = df.loc[codon].unique()
    codon_values[codon] = list(values)
    if len(values) == 1:
        conservative_codons.append(codon)

print(f"\n=== Консервативные кодоны ===")
print(f"Количество: {len(conservative_codons)} из {len(df.index)} ({len(conservative_codons)/len(df.index)*100:.1f}%)")
print("Список:", conservative_codons)

variability = {}
for codon in df.index:
    values = df.loc[codon].unique()
    variability[codon] = len(values)

sorted_variability = sorted(variability.items(), key=lambda x: x[1], reverse=True)

print(f"\n=== Топ-10 самых вариабельных кодонов ===")
for codon, count in sorted_variability[:10]:
    print(f"{codon}: {count} различных значений -> {codon_values[codon]}")

standard_code = df.iloc[:, 0]

changes = {}
for codon in df.index:
    std_val = standard_code[codon]
    other_vals = df.loc[codon, df.columns[1:]]
    changes[codon] = sum(other_vals != std_val)

print(f"\n=== Кодоны, наиболее часто меняющиеся относительно стандарта ===")
sorted_changes = sorted(changes.items(), key=lambda x: x[1], reverse=True)
for codon, count in sorted_changes[:10]:
    print(f"{codon}: меняется в {count} вариантах")

def extract_code_number(col_name):
    match = re.match(r'^(\d+)', col_name)
    return int(match.group(1)) if match else None

code_numbers = [extract_code_number(col) for col in df.columns]

aa_one_to_three = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
    'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
    'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
    '*': 'Stop'}
unique_aa = sorted(set(df.values.flatten()))
aa_to_num = {aa: i for i, aa in enumerate(unique_aa)}
df_numeric = df.applymap(lambda x: aa_to_num[x])
colorbar_ticks = list(range(len(unique_aa)))
colorbar_labels = [aa_one_to_three.get(aa, aa) for aa in unique_aa]

plt.figure(figsize=(16, 12))
ax = sns.heatmap(df_numeric, cmap='tab20', cbar=True,
                 yticklabels=df.index,           
                 xticklabels=code_numbers,       
                 cbar_kws={'ticks': colorbar_ticks, 'label': 'Аминокислота'})
colorbar = ax.collections[0].colorbar
colorbar.set_ticklabels(colorbar_labels)

plt.title('Сравнение вариантов генетического кода')
plt.xlabel('Номер варианта кода (transl_table)')
plt.ylabel('Кодон')
plt.tight_layout()
plt.savefig('genetic_codes_heatmap.png', dpi=150)
plt.show()

with open('analysis_results.txt', 'w', encoding='utf-8') as f:
    f.write(f"РЕЗУЛЬТАТЫ АНАЛИЗА ВАРИАНТОВ ГЕНЕТИЧЕСКОГО КОДА\n")
    f.write(f"=============================================\n\n")
    f.write(f"Всего проанализировано вариантов: {df.shape[1]}\n")
    f.write(f"Всего кодонов: {df.shape[0]}\n\n")

    f.write(f"Консервативные кодоны (одинаковы во всех вариантах):\n")
    f.write(f"Количество: {len(conservative_codons)}\n")
    f.write(f"Список: {', '.join(conservative_codons)}\n\n")

    f.write(f"Топ-10 вариабельных кодонов (по числу различных значений):\n")
    for codon, cnt in sorted_variability:
        f.write(f"{codon}: {cnt} значений -> {codon_values[codon]}\n")

    f.write(f"\nКодоны, наиболее часто изменяющиеся относительно стандартного:\n")
    for codon, cnt in sorted_changes:
        f.write(f"{codon}: изменён в {cnt} вариантах\n")

print("\nРезультаты сохранены в analysis_results.txt")