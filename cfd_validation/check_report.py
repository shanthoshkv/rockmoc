"""Lint the LaTeX validation report: py -3.11 cfd_validation/check_report.py [path/to/report.tex]"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
tex_path = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parent / "validation" / "latex" / "rockmoc_CFD_validation_report.tex"
tex = tex_path.read_text(encoding="utf-8", errors="ignore")
log_path = tex_path.with_suffix(".log")
log = log_path.read_text(encoding="utf-8", errors="ignore") if log_path.exists() else ""

labels = set(re.findall(r"\\label\{([^}]*)\}", tex))
refs = set()
for m in re.findall(r"\\(?:[cC]ref|ref|autoref|pageref|vref)\{([^}]*)\}", tex):
    refs.update(x.strip() for x in m.split(","))
unref = sorted(labels - refs)
undef_tex = sorted(refs - labels)
undef_log = len(re.findall(r"undefined", log))
overfull = [float(x) for x in re.findall(r"Overfull \\hbox \(([\d.]+)pt too wide", log)]

# prose only: drop listings, tables, figures, bibliography, math, commands
prose = re.sub(r"\\begin\{(lstlisting|longtable|figure|thebibliography|titlepage)\}.*?\\end\{\1\}", " ", tex, flags=re.S)
prose = prose.split(r"\begin{document}")[-1]
prose = re.sub(r"(?<!\\)%.*", " ", prose)
prose = prose.replace(r"\%", "percent")
prose = re.sub(r"\$[^$]*\$", " X ", prose)
prose = re.sub(r"\\(?:texttt|url|cite|cref|Cref|ref|label)\{[^}]*\}", " X ", prose)
prose = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?", " ", prose)
prose = re.sub(r"[{}\\~]", " ", prose)
prose = re.sub(r"\s+", " ", prose)
sents = [s.strip() for s in re.split(r"(?<=[.!?])\s+(?=[A-Z])", prose) if len(s.split()) >= 4]
mean_len = sum(len(s.split()) for s in sents) / max(len(sents), 1)
passive = sum(1 for s in sents if re.search(r"\b(?:was|were|is|are|been)\s+\w+ed\b", s))

FILLER = ["it is worth noting", "worth noting", "it can be seen", "in order to", "crucial", "comprehensive", "robust", "leverage", "delve"]
low = prose.lower()
filler = {w: low.count(w) for w in FILLER if w in low}

print(f"labels: {len(labels)}  refs: {len(refs)}")
print(f"unreferenced labels: {len(unref)} {unref[:10]}")
print(f"undefined refs (tex): {len(undef_tex)} {undef_tex[:10]}  (log 'undefined' lines: {undef_log})")
print(f"overfull hbox >10pt: {sum(1 for x in overfull if x > 10)} of {len(overfull)}")
print(f"filler hits: {sum(filler.values())} {filler}")
print(f"mean sentence length: {mean_len:.1f} words over {len(sents)} sentences")
print(f"rough passive share: {100 * passive / max(len(sents), 1):.1f} %")
m = re.search(r"Output written on .*? \((\d+) pages", log)
print(f"pages: {m.group(1) if m else 'n/a'}")
