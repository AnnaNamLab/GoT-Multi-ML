#!/home/mip4018/miniforge3/envs/ase/bin/python
"""
Python replacement for Combine_IronThrone_Parallel_Output.R.

CLI args (same order as legacy R script):
    1) working directory (ironthrone_output)
    2) sample prefix
    3) PCR ratio threshold
    4) Levenshtein distance threshold
    5) duplicate cutoff
    6) worker count
"""

from __future__ import annotations

import csv
import math
import os
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


CONCAT_COLUMNS = ["BC", "UMI", "num.WT.in.dups", "num.MUT.in.dups", "num.amb.in.dups", "call.in.dups"]
COLLAPSED_COLUMNS = CONCAT_COLUMNS + ["WT.calls", "MUT.calls", "amb.calls"]


def split_semicolon(value: str) -> List[str]:
    if value is None:
        return []
    return str(value).split(";")


def parse_float(value: str) -> float:
    try:
        return float(value)
    except Exception:
        return 0.0


def format_num(value: float) -> str:
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    return f"{value:.15g}"


def write_tsv(path: Path, rows: Sequence[Dict[str, str]], columns: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            writer.writerow([row.get(col, "") for col in columns])


def write_empty_outputs(out_parent: Path, prefix: str) -> None:
    write_tsv(out_parent / f"{prefix}.summTable.concat.txt", [], CONCAT_COLUMNS)
    write_tsv(out_parent / f"{prefix}.summTable.concat.umi_collapsed.txt", [], COLLAPSED_COLUMNS)


def maybe_delegate_to_legacy_r(
    wd: Path,
    prefix: str,
    pcr_thresh: float,
    ld: float,
    dupcut: float,
    threads: int,
) -> int | None:
    # Default to the pure-Python implementation. The legacy R script can still
    # be used explicitly for debugging/regression checks if requested.
    if os.environ.get("IRONTHRONE_PURE_PYTHON", "").lower() in {"1", "true", "yes"}:
        return None
    if os.environ.get("IRONTHRONE_USE_LEGACY_R", "").lower() not in {"1", "true", "yes"}:
        return None

    rscript = shutil.which("Rscript")
    legacy_script = Path(__file__).with_suffix(".R")
    if rscript is None or not legacy_script.exists():
        return None

    completed = subprocess.run(
        [
            rscript,
            str(legacy_script),
            str(wd),
            prefix,
            format_num(pcr_thresh),
            format_num(ld),
            format_num(dupcut),
            str(threads),
        ],
        check=False,
    )
    return int(completed.returncode)


def discover_split_tables(wd: Path, prefix: str) -> List[Path]:
    out = []
    for entry in sorted(wd.iterdir()):
        if not entry.is_dir():
            continue
        p = entry / f"{prefix}.summTable.txt"
        if p.exists():
            out.append(p)
    return out


def read_split_rows(paths: Sequence[Path]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for p in paths:
        with p.open("r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                rows.append(dict(row))
    return rows


def classify_call(wt_dup: float, mut_dup: float, pcr_thresh: float) -> str:
    total = wt_dup + mut_dup
    if total == 0:
        return "AMB"
    wt_frac = wt_dup / total
    mut_frac = mut_dup / total
    if wt_frac > pcr_thresh:
        return "WT"
    if mut_frac > pcr_thresh:
        return "MUT"
    return "AMB"


def aggregate_bc_umi(split_rows: Sequence[Dict[str, str]]) -> Dict[str, Dict[str, List[float]]]:
    # Preserve first-seen order to match R unique()/aggregation behavior.
    per_bc: Dict[str, Dict[str, List[float]]] = {}

    for row in split_rows:
        umis = split_semicolon(row.get("UMI", ""))
        if len(umis) == 0:
            continue

        bcs = split_semicolon(row.get("BC", ""))
        wt_dups = split_semicolon(row.get("num.WT.in.dups", ""))
        mut_dups = split_semicolon(row.get("num.MUT.in.dups", ""))
        amb_dups = split_semicolon(row.get("num.amb.in.dups", ""))

        for i in range(len(umis)):
            bc = bcs[i] if i < len(bcs) else (bcs[-1] if bcs else "")
            umi = umis[i] if i < len(umis) else (umis[-1] if umis else "")
            wt = parse_float(wt_dups[i] if i < len(wt_dups) else (wt_dups[-1] if wt_dups else "0"))
            mut = parse_float(mut_dups[i] if i < len(mut_dups) else (mut_dups[-1] if mut_dups else "0"))
            amb = parse_float(amb_dups[i] if i < len(amb_dups) else (amb_dups[-1] if amb_dups else "0"))

            if bc not in per_bc:
                per_bc[bc] = {}
            if umi not in per_bc[bc]:
                per_bc[bc][umi] = [0.0, 0.0, 0.0]
            per_bc[bc][umi][0] += wt
            per_bc[bc][umi][1] += mut
            per_bc[bc][umi][2] += amb

    return per_bc


def build_concat_rows(per_bc: Dict[str, Dict[str, List[float]]], pcr_thresh: float) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for bc, umi_map in per_bc.items():
        umi_list: List[str] = []
        wt_list: List[str] = []
        mut_list: List[str] = []
        amb_list: List[str] = []
        call_list: List[str] = []

        # R aggregate(..., by=Group) returns UMI groups in sorted order.
        for umi in sorted(umi_map.keys()):
            wt, mut, amb = umi_map[umi]
            umi_list.append(umi)
            wt_list.append(format_num(wt))
            mut_list.append(format_num(mut))
            amb_list.append(format_num(amb))
            call_list.append(classify_call(wt, mut, pcr_thresh))

        rows.append(
            {
                "BC": bc,
                "UMI": ";".join(umi_list),
                "num.WT.in.dups": ";".join(wt_list),
                "num.MUT.in.dups": ";".join(mut_list),
                "num.amb.in.dups": ";".join(amb_list),
                "call.in.dups": ";".join(call_list),
            }
        )
    return rows


def max_allowed_distance(seq: str, ld: float) -> int:
    if ld < 0:
        return 0
    if ld < 1:
        return int(math.ceil(ld * len(seq)))
    return int(math.floor(ld))


def agrep_like_distance(pattern: str, text: str) -> int:
    # R's agrep() performs approximate substring matching and is directional:
    # agrep(pattern, text) is not generally symmetric. The collapse behavior
    # depends on that asymmetry.
    if pattern == text:
        return 0
    if len(pattern) == 0:
        return 0
    if len(text) == 0:
        return len(pattern)

    prev = [0] * (len(text) + 1)
    for i, cp in enumerate(pattern, start=1):
        curr = [i] + [0] * len(text)
        for j, ct in enumerate(text, start=1):
            cost = 0 if cp == ct else 1
            curr[j] = min(prev[j] + 1, curr[j - 1] + 1, prev[j - 1] + cost)
        prev = curr
    return min(prev)


def compute_match_list(umis: Sequence[str], ld: float) -> List[List[int]]:
    match_list: List[List[int]] = []
    dist_cache: Dict[Tuple[str, str], int] = {}

    for i, x in enumerate(umis):
        allowed = max_allowed_distance(x, ld)
        matches_i: List[int] = []
        for j, y in enumerate(umis):
            key = (x, y)
            dist = dist_cache.get(key)
            if dist is None:
                dist = agrep_like_distance(x, y)
                dist_cache[key] = dist
            if dist <= allowed:
                matches_i.append(j)
        match_list.append(matches_i)
    return match_list


def split_and_fill(values: str, n: int, default: str = "") -> List[str]:
    arr = split_semicolon(values)
    if len(arr) >= n:
        return arr[:n]
    if len(arr) == 0:
        return [default] * n
    return arr + [arr[-1]] * (n - len(arr))


def collapse_single_row(row: Dict[str, str], pcr_thresh: float, ld: float, dupcut: float) -> Dict[str, str]:
    umis = split_semicolon(row.get("UMI", ""))
    n = len(umis)
    wt = [parse_float(x) for x in split_and_fill(row.get("num.WT.in.dups", ""), n, "0")]
    mut = [parse_float(x) for x in split_and_fill(row.get("num.MUT.in.dups", ""), n, "0")]
    amb = [parse_float(x) for x in split_and_fill(row.get("num.amb.in.dups", ""), n, "0")]
    calls = split_and_fill(row.get("call.in.dups", ""), n, "AMB")

    if n == 0:
        return {
            "BC": row.get("BC", ""),
            "UMI": "",
            "num.WT.in.dups": "",
            "num.MUT.in.dups": "",
            "num.amb.in.dups": "",
            "call.in.dups": "",
            "WT.calls": "0",
            "MUT.calls": "0",
            "amb.calls": "0",
        }

    adjacency = compute_match_list(umis, ld)
    active = [True] * n

    while True:
        max_match = 1
        to_collapse = -1
        active_matches: Dict[int, List[int]] = {}

        for i in range(n):
            if not active[i]:
                continue
            matches_i = [j for j in adjacency[i] if active[j]]
            active_matches[i] = matches_i
            if len(matches_i) > max_match:
                max_match = len(matches_i)
                to_collapse = i

        if to_collapse < 0:
            break

        matches_t0: set[int] = set()
        matches_t1: set[int] = set(active_matches[to_collapse])
        while len(matches_t1) > len(matches_t0):
            to_add: set[int] = set()
            for idx in matches_t1:
                to_add.update(j for j in adjacency[idx] if active[j])
            matches_t0 = set(matches_t1)
            matches_t1 = to_add

        group = sorted(matches_t1)
        wt_sum = sum(wt[i] for i in group)
        mut_sum = sum(mut[i] for i in group)
        amb_sum = sum(amb[i] for i in group)
        wt[to_collapse] = wt_sum
        mut[to_collapse] = mut_sum
        amb[to_collapse] = amb_sum
        calls[to_collapse] = classify_call(wt_sum, mut_sum, pcr_thresh)

        for idx in group:
            if idx != to_collapse:
                active[idx] = False

    active_idx = [i for i in range(n) if active[i]]
    sort_idx = sorted(active_idx, key=lambda i: wt[i] + mut[i] + amb[i], reverse=True)
    keep_idx = [i for i in sort_idx if wt[i] + mut[i] >= dupcut]

    umis = [umis[i] for i in keep_idx]
    wt = [wt[i] for i in keep_idx]
    mut = [mut[i] for i in keep_idx]
    amb = [amb[i] for i in keep_idx]
    calls = [calls[i] for i in keep_idx]

    return {
        "BC": row.get("BC", ""),
        "UMI": ";".join(umis),
        "num.WT.in.dups": ";".join(format_num(x) for x in wt),
        "num.MUT.in.dups": ";".join(format_num(x) for x in mut),
        "num.amb.in.dups": ";".join(format_num(x) for x in amb),
        "call.in.dups": ";".join(calls),
        "WT.calls": str(sum(1 for c in calls if c == "WT")),
        "MUT.calls": str(sum(1 for c in calls if c == "MUT")),
        "amb.calls": str(sum(1 for c in calls if c == "AMB")),
    }


def collapse_worker(args: Tuple[int, Dict[str, str], float, float, float]) -> Tuple[int, Dict[str, str]]:
    idx, row, pcr_thresh, ld, dupcut = args
    return idx, collapse_single_row(row, pcr_thresh, ld, dupcut)


def collapse_rows(rows: Sequence[Dict[str, str]], pcr_thresh: float, ld: float, dupcut: float, threads: int) -> List[Dict[str, str]]:
    if len(rows) == 0:
        return []
    if threads <= 1 or len(rows) == 1:
        return [collapse_single_row(row, pcr_thresh, ld, dupcut) for row in rows]

    tasks = [(i, row, pcr_thresh, ld, dupcut) for i, row in enumerate(rows)]
    chunksize = max(1, len(tasks) // max(threads * 8, 1))
    try:
        with ProcessPoolExecutor(max_workers=threads) as pool:
            out = list(pool.map(collapse_worker, tasks, chunksize=chunksize))
        out.sort(key=lambda x: x[0])
        return [row for _, row in out]
    except Exception:
        # Fallback for restricted environments where multiprocessing semaphores are unavailable.
        return [collapse_single_row(row, pcr_thresh, ld, dupcut) for row in rows]


def main(argv: Sequence[str]) -> int:
    if len(argv) != 7:
        print(
            "Usage: Combine_IronThrone_Parallel_Output.py <wd> <prefix> <pcr_thresh> <ld> <dupcut> <threads>",
            file=sys.stderr,
        )
        return 1

    wd = Path(argv[1]).resolve()
    prefix = argv[2]
    pcr_thresh = float(argv[3])
    ld = float(argv[4])
    dupcut = float(argv[5])
    threads = max(1, int(float(argv[6])))

    delegated_rc = maybe_delegate_to_legacy_r(wd, prefix, pcr_thresh, ld, dupcut, threads)
    if delegated_rc is not None:
        return delegated_rc

    out_parent = wd.parent
    split_files = discover_split_tables(wd, prefix)
    print(f"[combine] discovered {len(split_files)} split tables for {prefix}", flush=True)
    if len(split_files) == 0:
        write_empty_outputs(out_parent, prefix)
        return 0

    split_rows = read_split_rows(split_files)
    print(f"[combine] loaded {len(split_rows)} split rows", flush=True)
    if len(split_rows) == 0:
        write_empty_outputs(out_parent, prefix)
        return 0

    per_bc = aggregate_bc_umi(split_rows)
    print(f"[combine] aggregated into {len(per_bc)} barcode rows", flush=True)
    if len(per_bc) == 0:
        write_empty_outputs(out_parent, prefix)
        return 0

    concat_rows = build_concat_rows(per_bc, pcr_thresh)
    write_tsv(out_parent / f"{prefix}.summTable.concat.txt", concat_rows, CONCAT_COLUMNS)
    print(f"[combine] wrote {len(concat_rows)} concat rows", flush=True)

    collapsed_rows = collapse_rows(concat_rows, pcr_thresh, ld, dupcut, threads)
    write_tsv(out_parent / f"{prefix}.summTable.concat.umi_collapsed.txt", collapsed_rows, COLLAPSED_COLUMNS)
    print(f"[combine] wrote {len(collapsed_rows)} collapsed rows", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
