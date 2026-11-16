#!/usr/bin/env python3
"""Fetch selected PubMed records and build a reproducible RIS/TSV/HTML package."""
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import time
import csv, html, json, xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis/citations/final_30"
PMIDS = [x.strip() for x in (ROOT / "analysis/citations/selected_pmids.txt").read_text().splitlines() if x.strip()]

def txt(node, path):
    n = node.find(path)
    return "" if n is None else "".join(n.itertext()).strip()

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows=[]
    articles=[]
    for start in range(0, len(PMIDS), 10):
        batch=PMIDS[start:start+10]
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urlencode({
            "db":"pubmed", "id":",".join(batch), "retmode":"xml", "tool":"ROVOT_citation_audit"})
        last_error=None
        for attempt in range(3):
            try:
                with urlopen(url, timeout=60) as response: payload=response.read()
                articles.extend(ET.fromstring(payload).findall(".//PubmedArticle"))
                break
            except Exception as exc:
                last_error=exc; time.sleep(1.5*(attempt+1))
        else: raise RuntimeError(f"PubMed batch failed after retries: {batch}") from last_error
        time.sleep(0.4)
    for article in articles:
        med = article.find("MedlineCitation")
        art = med.find("Article")
        pmid = txt(med, "PMID")
        title = txt(art, "ArticleTitle")
        journal = txt(art, "Journal/Title")
        year = txt(art, "Journal/JournalIssue/PubDate/Year") or txt(art, "Journal/JournalIssue/PubDate/MedlineDate")[:4]
        volume = txt(art, "Journal/JournalIssue/Volume")
        issue = txt(art, "Journal/JournalIssue/Issue")
        pages = txt(art, "Pagination/MedlinePgn") or txt(art, "ELocationID")
        authors=[]
        for a in art.findall("AuthorList/Author"):
            coll = txt(a, "CollectiveName")
            name = coll or ", ".join(x for x in (txt(a,"LastName"), txt(a,"ForeName")) if x)
            if name: authors.append(name)
        doi=""
        for aid in article.findall("PubmedData/ArticleIdList/ArticleId"):
            if aid.attrib.get("IdType")=="doi": doi=(aid.text or "").strip()
        abstract=" ".join("".join(x.itertext()).strip() for x in art.findall("Abstract/AbstractText"))
        rows.append(dict(pmid=pmid,title=title,journal=journal,year=year,volume=volume,issue=issue,pages=pages,
                         authors=authors,doi=doi,abstract=abstract,url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                         metadata_source="NCBI PubMed E-utilities"))
    # Add the authors' prior RV pressure-overload study supplied by the user.
    rows.append(dict(pmid="", title="Correlative study on the expression of adipokine CTRP9 in a rat model of right ventricular outflow tract obstruction",
                     journal="Archives of Medical Science", year="2025", volume="", issue="", pages="",
                     authors=["Shen, Yanming","Chen, Shaoqin","Fang, Guanhua","Huang, Hongmiao","Chen, Liangwan","Liao, Dongshan"],
                     doi="10.5114/aoms/199929", abstract="",
                     url="https://doi.org/10.5114/aoms/199929", metadata_source="Author-supplied article and DOI record"))
    # Add the multiple-testing method reference, which is not indexed in PubMed.
    rows.append(dict(pmid="", title="Controlling the false discovery rate: a practical and powerful approach to multiple testing",
                     journal="Journal of the Royal Statistical Society: Series B", year="1995", volume="57", issue="1", pages="289-300",
                     authors=["Benjamini, Yoav","Hochberg, Yosef"], doi="10.1111/j.2517-6161.1995.tb02031.x", abstract="",
                     url="https://doi.org/10.1111/j.2517-6161.1995.tb02031.x", metadata_source="DOI record"))
    if len(rows)!=32: raise RuntimeError(f"Expected 32 references, got {len(rows)}")
    with (OUT/"references.ris").open("w",encoding="utf-8") as f:
        for r in rows:
            f.write("TY  - JOUR\n")
            f.write(f"TI  - {r['title']}\n")
            for a in r['authors']: f.write(f"AU  - {a}\n")
            f.write(f"T2  - {r['journal']}\nPY  - {r['year']}\n")
            if r['volume']: f.write(f"VL  - {r['volume']}\n")
            if r['issue']: f.write(f"IS  - {r['issue']}\n")
            if r['pages']: f.write(f"SP  - {r['pages']}\n")
            if r['doi']: f.write(f"DO  - {r['doi']}\n")
            f.write(f"UR  - {r['url']}\nER  -\n\n")
    with (OUT/"references.tsv").open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["index","first_author","year","title","journal","doi","pmid","url","metadata_source"],delimiter="\t")
        w.writeheader()
        for i,r in enumerate(rows,1):
            w.writerow(dict(index=i,first_author=r['authors'][0] if r['authors'] else "",year=r['year'],title=r['title'],journal=r['journal'],doi=r['doi'],pmid=r['pmid'],url=r['url'],metadata_source=r['metadata_source']))
    (OUT/"references.json").write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding="utf-8")
    trs=[]
    for i,r in enumerate(rows,1):
        link=r['url']; trs.append(f"<tr><td>{i}</td><td>{html.escape(r['authors'][0] if r['authors'] else '')}</td><td>{html.escape(r['year'])}</td><td><a href='{html.escape(link)}'>{html.escape(r['title'])}</a></td><td>{html.escape(r['journal'])}</td><td>{html.escape(r['doi'])}</td></tr>")
    page="""<!doctype html><meta charset='utf-8'><title>ROVOT citation set</title><style>body{font-family:Arial;max-width:1500px;margin:30px auto}table{border-collapse:collapse;width:100%}th,td{border:1px solid #ccc;padding:7px;vertical-align:top}th{background:#eee;position:sticky;top:0}input{width:100%;padding:8px;margin-bottom:12px}</style><input id='q' placeholder='Filter references'><table id='t'><thead><tr><th>#</th><th>First author</th><th>Year</th><th>Title</th><th>Journal</th><th>DOI</th></tr></thead><tbody>"""+"".join(trs)+"</tbody></table><script>q.oninput=()=>{for(const r of t.tBodies[0].rows)r.style.display=r.innerText.toLowerCase().includes(q.value.toLowerCase())?'':'none'}</script>"
    (OUT/"references.html").write_text(page,encoding="utf-8")
    print(f"Fetched {len(PMIDS)} PubMed records; total references: {len(rows)}")
    print(OUT/"references.html"); print(OUT/"references.ris")
if __name__=="__main__": main()
