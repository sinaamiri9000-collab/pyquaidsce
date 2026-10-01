"""Extract the user's printed Stata benchmark, preserving coefficient labels.

Scientific notation must be parsed as a complete token; Stata displays some
small Gamma estimates as, for example, -6.07e-06.
"""

from pathlib import Path
import json
import re


ROOT = Path(__file__).resolve().parent


def parse_reference(path):
    text = Path(path).read_text()
    body = text.split("Basic translog model", 1)[1].split("Note: alpha estimates", 1)[0]
    coefficients = []
    block = ""
    demo = ""
    for line in body.splitlines():
        key = re.fullmatch(r"\s*(alpha|Gamma|Nu)\s*\|\s*", line)
        if key:
            block = key[1]
        key = re.search(r"Good#c\.(\w+)\s*\|", line)
        if key:
            demo = key[1]
        row = re.match(r"^\s*(\d+(?:#\d+)?)\s*\|\s*(\S+)\s+(\S+)", line)
        if row:
            label, estimate, se = row.groups()
            name = (f"alpha:alpha_{label}" if block == "alpha" else
                    f"Gamma:gamma_{label.replace('#', '_')}" if block == "Gamma" else
                    f"Nu:nu_{demo}_{label}")
            coefficients.append(dict(name=name, estimate=float(estimate), se=float(se)))
    normalized = text.split("Normalized parameter", 1)[1].split("Note: alpha estimates", 1)[0]
    row = re.search(r"^\s*14\s*\|\s*(\S+)\s+(\S+)", normalized, re.M)
    coefficients.insert(13, dict(name="alpha:alpha_14", estimate=float(row[1]), se=float(row[2])))
    assert len(coefficients) == 217
    elasticities = []
    price = None
    for line in text.split("estat elasticities, uncompensated atmeans", 1)[1].splitlines():
        key = re.match(r"^Good (\d+)\s*\|", line)
        if key:
            price = int(key[1])
        row = re.match(r"^\s*(\d+)\s*\|\s*(\S+)\s+(\S+)", line)
        if row and price is not None:
            elasticities.append(dict(good=int(row[1]), price=price,
                                     estimate=float(row[2]), se=float(row[3])))
    assert len(elasticities) == 196
    r2 = [float(v) for v in re.findall(r"Good \d+\s*=\s*([\d.]+)", body)]
    assert len(r2) == 14
    return dict(model="basic translog", censor=False, method="ifgnls",
                demographics_method="translating", nobs=6573,
                elasticity_display_nobs=6848, llf=124232.06,
                coefficients=coefficients, uncompensated=elasticities, r2=r2)


if __name__ == "__main__":
    reference = parse_reference(ROOT / "reference" / "demandsysresults.txt")
    (ROOT / "reference" / "stata_reference.json").write_text(
        json.dumps(reference, indent=2) + "\n")
