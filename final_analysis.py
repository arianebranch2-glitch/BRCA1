# Goal: looks up BRCA1’s Ensembl gene ID, downloads all overlapping variants (saved to brca1_variants.json), then computes per-variant GC% and genomic length, writes sequence_summary.csv, and optionally saves two histograms (GC% and length).
# Aurthor: A'Riane Branch 

import requests, json

headers = {"Accept": "application/json"}

# Lookup BRCA1 -> ENSG
url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/BRCA1"
response = requests.get(url, headers=headers)
if response.status_code == 200:
    data = response.json()
    ensg = data["id"]
    print("Ensembl ID:", ensg)

    # All variants for that gene (>=1,000 entries)
    url = f"https://rest.ensembl.org/overlap/id/{ensg}"
    params = {"feature": "variation"}
    response = requests.get(url, params=params, headers=headers)
    if response.status_code == 200:
        variants = response.json()
        print("Entries downloaded:", len(variants))

        # Save dataset
        with open("brca1_variants.json", "w") as f:
            json.dump(variants, f, indent=2)
        print("Saved to brca1_variants.json")
    else:
        print("Variants failed with code:", response.status_code)
else:
    print("Lookup failed with code:", response.status_code)
# Calculate GC-content and sequence length for BRCA1 Variants
import json

def calculate_gc_percent_str(sequence):
    if not sequence:
        return ""
    s = sequence.upper().replace("-", "")
    if not s:
        return ""
    gc = s.count("G") + s.count("C")
    return f"{(100.0 * gc / len(s)):.3f}"

def pick_one_allele_list(alleles):
    # first non-deletion from a list (strings or dicts)
    if not alleles:
        return ""
    for item in alleles:
        a = item.get("allele", "") if isinstance(item, dict) else str(item)
        if a and a != "-":
            return a
    return ""

def variant_length_str(record):
    try:
        start = int(record["start"])
        end = int(record["end"])
        return str(end - start + 1)
    except:
        return ""

def main():
    with open("brca1_variants.json") as f:
        variants = json.load(f)

    lines = ["variant_id,length,gc_percent"]
    count_len_100 = 0

    for v in variants:
        allele = pick_one_allele_list(v.get("alleles"))
        gc_str = calculate_gc_percent_str(allele)
        length_str = variant_length_str(v)
        if length_str == "100":
            count_len_100 += 1
        vid = v.get("id", "")
        lines.append(f"{vid},{length_str},{gc_str}")
# Create sequence_summary.csv 
    with open("sequence_summary.csv", "w") as out:
        out.write("\n".join(lines))

    print("n_variants:", len(variants))
    print("count_length_equals_100:", count_len_100)

if __name__ == "__main__":
    main()
# Make plots 
#Input: sequence_summary.csv
#Output: figure1.png (GC% histogram), figure2.png (length histogram)

import matplotlib.pyplot as plt
gc_vals = []
lengths = []

with open("sequence_summary.csv", "r") as f:
    f.readline()
    for line in f:
        parts = line.split(",")
        if len(parts) < 3:
            continue
        length_str = parts[1].strip()
        gc_str = parts[2].strip()
        if length_str.isdigit():
            lengths.append(int(length_str))
        if gc_str:
            head = gc_str.split(".", 1)[0]
            if head.isdigit():
                gc_vals.append(int(head))
plt.figure()
plt.hist(gc_vals, bins=50)
plt.title("Allele GC% (integar) - BRAC1 variants")
plt.xlabel("Length (bp)")
plt.ylabel("count")
plt.savefig("figure1.png")

plt.figure()
plt.hist(lengths, bins=50)
plt.title("Variants Genomic Length - BRAC1 variants")
plt.xlabel("length (bp)")
plt.ylabel("count")
plt.savefig("figure2.png")
print("Saved figure1.png and figure2.png")

