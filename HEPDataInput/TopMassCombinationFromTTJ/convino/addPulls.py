"""Append the CMS post-fit nuisance values (HEPData ins2106483 v2) to a Convino input file.

The Convino file's [correlation matrix] rows are in the header order of the
HEPData `fit_np_correlation` table (verified: constraints and correlations
agree to 1e-16), so the name map is positional. Adds

    [nuisance values]   name = fit_obs_central   (one line per NP)
    [systematics]       CMS_norm_tt0jet = free   (rateTT0Jet has no prior)

and leaves every existing number untouched.

    python3 addPulls.py <convino_input.txt> [more files...]
"""
import json
import re
import sys
from pathlib import Path

HEP = Path(__file__).resolve().parent.parent / "CMS_13TeV_hepdata_v2"
FREE = {"rateTT0Jet"}


def load_hepdata():
    corr = json.load(open(HEP / "fit_np_correlation.json"))
    order = [h["name"] for h in corr["headers"][1:]]
    pulls = json.load(open(HEP / "np_impacts_pulls.json"))
    col = [h["name"] for h in pulls["headers"][1:]].index("fit_obs_central")
    central = {r["x"][0]["value"]: float(r["y"][col]["value"]) for r in pulls["values"]}
    return order, central


def add_pulls(path: Path, order, central):
    text = path.read_text()
    if "[nuisance values]" in text:
        print(f"{path}: already has [nuisance values], skipped")
        return
    block = re.search(r"\[correlation matrix\](.*?)\[end correlation matrix\]", text, re.S)
    rows = [l.split()[0] for l in block.group(1).strip().splitlines() if l.strip()]
    assert len(rows) == len(order), f"{path}: {len(rows)} rows vs {len(order)} HEPData names"

    lines, free = [], []
    for conv_name, hep_name in zip(rows, order):
        if hep_name not in central:          # POI rows (rate_ttj*)
            continue
        lines.append(f"    {conv_name} = {central[hep_name]:.6g}")
        if hep_name in FREE:
            free.append(f"    {conv_name} = free")
    assert len(lines) == len(central), (len(lines), len(central))

    text = text.replace("[systematics]\n", "[systematics]\n" + "\n".join(free) + "\n", 1)
    text = text.rstrip("\n") + "\n\n[nuisance values]\n" + "\n".join(lines) + "\n[end nuisance values]\n"
    path.write_text(text)
    print(f"{path}: added {len(lines)} nuisance values, {len(free)} free")


if __name__ == "__main__":
    order, central = load_hepdata()
    for p in sys.argv[1:]:
        add_pulls(Path(p), order, central)
