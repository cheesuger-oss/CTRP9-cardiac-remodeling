#!/usr/bin/env python3
"""Map Ensembl gene IDs to symbols for human/rat/pig via mygene.info batch API.
Saves mapping tables to data/processed/.
"""
import os, time, json, gzip
import requests
import pandas as pd

OUT = "data/processed"
os.makedirs(OUT, exist_ok=True)

def map_ids(ids, species, fields="symbol,name,ensembl.gene", chunk=1000, sleep=0.4):
    """Batch map Ensembl gene IDs to symbols via mygene.info."""
    mapping = {}
    for i in range(0, len(ids), chunk):
        chunk_ids = ids[i:i+chunk]
        for attempt in range(4):
            try:
                r = requests.post("https://mygene.info/v3/query",
                                  json={"q": chunk_ids, "species": species,
                                        "fields": fields, "size": len(chunk_ids)},
                                  timeout=180)
                if r.status_code == 200:
                    break
            except Exception as e:
                print("  retry", attempt, e)
                time.sleep(3)
        else:
            print("  FAILED chunk", i)
            continue
        res = r.json()
        for item in res:
            q = item.get("query")
            sym = item.get("symbol")
            if sym:
                mapping[q] = sym
        time.sleep(sleep)
        print(f"  {species}: {min(i+chunk,len(ids))}/{len(ids)} mapped so far {len(mapping)}")
    return mapping

# ---------- Human: GSE36761 ----------
expr = pd.read_csv("data/raw/GSE36761_gene_expression_levels_normalized.txt.gz",
                   sep="\t", compression="gzip")
human_ids = [str(x) for x in expr["Ensembl Gene ID"] if str(x).startswith("ENSG")]
print("human IDs:", len(human_ids))
human_map = map_ids(human_ids, "human")
pd.Series(human_map).to_csv(f"{OUT}/mapping_human_ensg_symbol.tsv", sep="\t", header=["symbol"])
print("human mapped:", len(human_map))

# ---------- Rat: GSE119754 ENSRNOG ----------
cnt = pd.read_csv("data/raw/GSE119754_all_counts.csv.gz", compression="gzip", index_col=0)
rat_ids = [str(x) for x in cnt.index if str(x).startswith("ENSRNOG")]
print("rat ENSRNOG IDs:", len(rat_ids))
rat_map = map_ids(rat_ids, "rat")
pd.Series(rat_map).to_csv(f"{OUT}/mapping_rat_ensrnog_symbol.tsv", sep="\t", header=["symbol"])
print("rat mapped:", len(rat_map))

# ---------- Pig: GSE224428 ENSSSCG ----------
pig = pd.read_csv("data/raw/GSE224428_subread_counts_RV.txt.gz", sep="\t", compression="gzip", skiprows=1)
pig_ids = [str(x) for x in pig["Geneid"] if str(x).startswith("ENSSSCG")]
print("pig ENSSSCG IDs:", len(pig_ids))
pig_map = map_ids(pig_ids, "pig")
pd.Series(pig_map).to_csv(f"{OUT}/mapping_pig_enssscg_symbol.tsv", sep="\t", header=["symbol"])
print("pig mapped:", len(pig_map))

# sanity checks
for g, m in [("C1QTNF9", human_map), ("C1qtnf9", rat_map)]:
    hits = {k: v for k, v in m.items() if v == g}
    print(f"{g}: {hits}")
