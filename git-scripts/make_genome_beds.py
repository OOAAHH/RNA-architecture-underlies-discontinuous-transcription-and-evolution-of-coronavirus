# -*- coding: utf-8 -*-
"""
Generate genome annotation BED files used by this repository.

Writes BED6 files to `data/genome/` by default.
"""

from __future__ import annotations

from collections import OrderedDict
import argparse
import importlib.util
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


def _load_attr_from_py(py_path: Path, attr: str):
    spec = importlib.util.spec_from_file_location(py_path.stem, py_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load module from: {py_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, attr):
        raise AttributeError(f"{py_path} does not define `{attr}`")
    return getattr(module, attr)


def _write_bed6(out_path: Path, chrom: str, gene_sites: Dict[str, Tuple[int, int]]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        f.write("chrom\tstart\tend\tname\tscore\tstrand\n")
        for gene, (start, end) in gene_sites.items():
            f.write(f"{chrom}\t{int(start)}\t{int(end)}\t{gene}\t0\t+\n")


def _split_top_level_commas(s: str) -> List[str]:
    parts: List[str] = []
    depth = 0
    buf: List[str] = []
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == "," and depth == 0:
            part = "".join(buf).strip()
            if part:
                parts.append(part)
            buf = []
            continue
        buf.append(ch)
    tail = "".join(buf).strip()
    if tail:
        parts.append(tail)
    return parts


def _parse_location_to_intervals(loc: str) -> List[Tuple[int, int]]:
    loc = loc.strip().replace("<", "").replace(">", "")
    if loc.startswith("complement(") and loc.endswith(")"):
        return _parse_location_to_intervals(loc[len("complement(") : -1])
    for wrapper in ("join(", "order("):
        if loc.startswith(wrapper) and loc.endswith(")"):
            inner = loc[len(wrapper) : -1]
            intervals: List[Tuple[int, int]] = []
            for part in _split_top_level_commas(inner):
                intervals.extend(_parse_location_to_intervals(part))
            return intervals
    if ".." in loc:
        a, b = loc.split("..", 1)
        return [(int(a), int(b))]
    if "^" in loc:
        a, _b = loc.split("^", 1)
        return [(int(a), int(a))]
    return [(int(loc), int(loc))]


def _intervals_span(intervals_1_based_inclusive: List[Tuple[int, int]]) -> Tuple[int, int]:
    starts = [s for s, _e in intervals_1_based_inclusive]
    ends = [e for _s, e in intervals_1_based_inclusive]
    return min(starts), max(ends)


def _to_bed_coords(start_1_based: int, end_1_based_inclusive: int) -> Tuple[int, int]:
    # BED: 0-based start, end-exclusive. For 1-based inclusive end, BED end is identical.
    return start_1_based - 1, end_1_based_inclusive


def _fetch_ncbi_genbank(accession: str, timeout_s: int = 30) -> str:
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "gb",
        "retmode": "text",
    }
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(
        params
    )
    req = urllib.request.Request(url, headers={"User-Agent": "make_genome_beds.py"})
    with urllib.request.urlopen(req, timeout=timeout_s) as resp:
        return resp.read().decode("utf-8", errors="replace")


def _parse_genbank_features(genbank_text: str):
    lines = genbank_text.splitlines()
    in_features = False
    current = None
    current_quals: Dict[str, str] = {}
    saw_qualifiers = False
    capturing_key: Optional[str] = None
    capturing_val: List[str] = []
    in_quote = False

    def flush_capture():
        nonlocal capturing_key, capturing_val, in_quote
        if capturing_key is None:
            return
        value = "".join(capturing_val).strip()
        if value.startswith('"') and value.endswith('"') and len(value) >= 2:
            value = value[1:-1]
        current_quals[capturing_key] = value
        capturing_key = None
        capturing_val = []
        in_quote = False

    def flush_feature():
        nonlocal current, current_quals, saw_qualifiers
        flush_capture()
        if current is None:
            return None
        feat = {**current, "qualifiers": current_quals}
        current = None
        current_quals = {}
        saw_qualifiers = False
        return feat

    features = []
    for line in lines:
        if line.startswith("FEATURES             Location/Qualifiers"):
            in_features = True
            continue
        if not in_features:
            continue
        if line.startswith("ORIGIN"):
            feat = flush_feature()
            if feat is not None:
                features.append(feat)
            break

        if line.startswith("     ") and not line.startswith("                     "):
            feat = flush_feature()
            if feat is not None:
                features.append(feat)
            key = line[5:21].strip()
            location = line[21:].strip()
            current = {"key": key, "location": location}
            saw_qualifiers = False
            continue

        if current is None:
            continue

        if line.startswith("                     "):
            payload = line[21:].rstrip()

            # Qualifier line
            if payload.startswith("/") and not in_quote:
                saw_qualifiers = True
                flush_capture()
                if "=" in payload:
                    qk, qv = payload[1:].split("=", 1)
                    capturing_key = qk.strip()
                    capturing_val = [qv.strip()]
                    in_quote = (qv.count('"') % 2 == 1)
                    if not in_quote:
                        flush_capture()
                else:
                    current_quals[payload[1:].strip()] = ""
                continue

            # Location continuation (rare, but happens for long joins)
            if not saw_qualifiers and capturing_key is None and not payload.startswith("/"):
                current["location"] += payload.strip()
                continue

            # Qualifier continuation (multi-line quotes)
            if capturing_key is not None:
                capturing_val.append(payload.strip())
                if in_quote and ("".join(capturing_val).count('"') % 2 == 0):
                    in_quote = False
                    flush_capture()

    return features


def _gene_sites_from_ncbi(accession: str, leader_len: int) -> "OrderedDict[str, Tuple[int, int]]":
    gb = _fetch_ncbi_genbank(accession)
    feats = _parse_genbank_features(gb)

    genome_len_bed: Optional[int] = None
    utr5_start_bed: Optional[int] = None
    utr5_end_bed: Optional[int] = None
    utr3: Optional[Tuple[int, int]] = None
    genes: List[Tuple[str, int, int]] = []
    cds: List[Tuple[str, int, int, str]] = []  # (name, start, end, raw_product)

    def normalize_name(raw_name: str) -> str:
        n = raw_name.strip()
        low = n.lower()
        if low in {"s", "e", "m", "n"}:
            return n.upper()
        if low == "polyprotein":
            return "ORF1a"
        if "polyprotein 1ab" in low or "replicase polyprotein 1ab" in low:
            return "ORF1ab"
        if "polyprotein 1a" in low or "replicase polyprotein 1a" in low:
            return "ORF1a"
        if "spike" in low:
            return "S"
        if "envelope" in low:
            return "E"
        if "membrane" in low:
            return "M"
        if "nucleocapsid" in low:
            return "N"
        return n

    for feat in feats:
        key = feat["key"]
        loc = feat["location"]
        qualifiers = feat.get("qualifiers", {})

        if key not in {"source", "5'UTR", "3'UTR", "gene", "CDS"}:
            continue

        intervals = _parse_location_to_intervals(loc)
        s1, e1 = _intervals_span(intervals)
        s_bed, e_bed = _to_bed_coords(s1, e1)

        if key == "source":
            genome_len_bed = e_bed
            continue
        if key == "5'UTR":
            utr5_start_bed, utr5_end_bed = s_bed, e_bed
            continue
        if key == "3'UTR":
            utr3 = (s_bed, e_bed)
            continue

        raw_name = qualifiers.get("gene") or qualifiers.get("product")
        if not raw_name:
            continue
        name = normalize_name(raw_name)

        # ORF1ab frameshift is often represented as join(...,...).
        is_frameshift = (
            "orf1ab" in name.lower()
            or "ribosomal_slippage" in qualifiers
            or "ribosomal_slippage" in qualifiers.keys()
            or "ribosomal frameshift" in (qualifiers.get("note", "").lower())
        )
        if is_frameshift:
            join_intervals = _parse_location_to_intervals(loc)
            if len(join_intervals) == 2:
                (a_s, a_e), (b_s, b_e) = join_intervals
                a0, a1 = _to_bed_coords(a_s, a_e)
                b0, b1 = _to_bed_coords(b_s, b_e)
                genes.append(("ORF1a", a0, a1))
                genes.append(("ORF1b", b0, b1))
                continue

        if key == "gene":
            genes.append((name, s_bed, e_bed))
        else:  # CDS
            cds.append((name, s_bed, e_bed, raw_name))

    # Merge in CDS-derived features (some virus records only annotate CDS, not `gene`).
    if cds:
        seen = {name for (name, _s, _e) in genes}
        for name, s, e, _raw in cds:
            if name not in seen:
                genes.append((name, s, e))
                seen.add(name)

    genes.sort(key=lambda x: (x[1], x[2], x[0]))
    out: "OrderedDict[str, Tuple[int, int]]" = OrderedDict()

    # Determine 5'UTR / 3'UTR boundaries if not annotated.
    if utr5_start_bed is None and genes:
        utr5_start_bed = 0
        utr5_end_bed = genes[0][1]
    if utr3 is None and genes and genome_len_bed is not None:
        utr3 = (genes[-1][2], genome_len_bed)

    # leader + UTR_5 from 5'UTR if available
    if utr5_start_bed is not None and utr5_end_bed is not None:
        leader_start = utr5_start_bed
        leader_end = min(utr5_end_bed, utr5_start_bed + leader_len)
        out["leader"] = (leader_start, leader_end)
        if leader_end < utr5_end_bed:
            out["UTR_5"] = (leader_end, utr5_end_bed)

    for name, s, e in genes:
        if name in {"leader", "UTR_5", "UTR_3"}:
            continue
        out.setdefault(name, (s, e))

    if utr3 is not None:
        out["UTR_3"] = utr3

    return out


def main(argv: Iterable[str] | None = None) -> int:
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent

    known = {
        "NC_045512.2": (script_dir / "genome_382.py", "dict_gene"),
        "MK584552.1": (script_dir / "genome_pedv_deletion.py", "dict_gene"),
        "KT336560.1": (script_dir / "genome_pdcov.py", "dict_gene"),
    }

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=repo_root / "data" / "genome",
        help="Output directory for BED files (default: data/genome).",
    )
    parser.add_argument(
        "--accessions",
        nargs="*",
        default=[
            "NC_045512.2",
            "MK584552.1",
            "KY363867.1",
            "KT336560.1",
            "MK878536.1",
            "NC_004718.3",
            "NC_019843.3",
        ],
        help="Which accessions to generate (default: the 7 accessions referenced by read_genome_sites.py).",
    )
    parser.add_argument(
        "--ncbi",
        action="store_true",
        help="For accessions without a local dict, download GenBank from NCBI and build a best-effort BED.",
    )
    parser.add_argument(
        "--leader-len",
        type=int,
        default=70,
        help="Leader length to split from 5'UTR when using --ncbi (default: 70).",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    out_dir: Path = args.out_dir
    accessions = list(args.accessions)

    for acc in accessions:
        if acc in known:
            py_path, dict_name = known[acc]
            gene_dict = _load_attr_from_py(py_path, dict_name)
        else:
            if not args.ncbi:
                raise SystemExit(
                    f"Accession {acc} has no local coordinate dict.\n"
                    "Re-run with `--ncbi` (requires network access), or add a local dict module."
                )
            gene_dict = _gene_sites_from_ncbi(acc, leader_len=args.leader_len)

        out_path = out_dir / f"{acc}.bed"
        _write_bed6(out_path, chrom=acc, gene_sites=gene_dict)
        print(f"Wrote {out_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
