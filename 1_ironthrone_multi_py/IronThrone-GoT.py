#!/home/mip4018/miniforge3/envs/ase/bin/python
"""
Python port of IronThrone-GoT.

Goals of this port:
- Preserve CLI compatibility with the original Perl entrypoint.
- Preserve output schema and downstream file names expected by IronThroneRunner.sh.
- Keep runtime performant for large FRP/GoT runs.

Notes:
- The core decision logic is intentionally kept close to the Perl implementation.
- Some historical quirks are kept for behavioral compatibility with existing pipelines.
"""

from __future__ import annotations

import argparse
from collections import Counter
import gzip
import math
import os
import shutil
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


WHAT_IS_THIS_CODE = """

     ##### created by the Landau lab (@New York Genome Center)
     ##### for processing GoT amplicon data
"""


def _round_perl(value: float) -> int:
    """Mimic Perl's sprintf("%.0f", x): round-half-away-from-zero."""
    if value >= 0:
        return int(math.floor(value + 0.5))
    return int(math.ceil(value - 0.5))


class IronThroneGoT:
    def __init__(self, args: argparse.Namespace) -> None:
        self.run_mode: str = args.run
        self.fastq_r1: str = args.fastqR1
        self.fastq_r2: str = args.fastqR2
        self.conf_file: str = args.config
        self.whitelist_file: str = args.whitelist

        self.mmtch: float = float(args.mmtch)
        self.post_p: float = float(args.postP)
        self.dupcut: int = int(args.dupcut)
        self.len_barcode: int = int(args.bclen)
        self.len_umi: int = int(args.umilen)
        self.thread_max: int = max(1, int(args.thread))
        self.keep_outs: int = int(args.keepouts)
        self.verbose: int = int(args.verbose)

        self.allow_err: int = 0
        self.out_dir: Path = Path(args.outdir)
        self.samplename: str = args.sample
        self.log_file: str = args.log

        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.log_path = self.out_dir / self.log_file
        self.log_path.write_text("")
        self._log_handle = self.log_path.open("a", encoding="utf-8")

        # Config-derived globals (kept close to Perl naming for parity)
        self.primer_seq: Optional[str] = None
        self.primer_alw: Optional[int] = None
        self.general_seq: Optional[str] = None
        self.general_alw: Optional[int] = None
        self.general_start: Optional[int] = None
        self.general_end: Optional[int] = None
        self.wt_seq: Optional[str] = None
        self.wt_alw: Optional[int] = None
        self.pcr2_fw_wt_seq: Optional[str] = None
        self.pcr2_fw_wt_alw: Optional[int] = None
        self.pcr2_fw_mut_seq: Optional[str] = None
        self.pcr2_fw_mut_alw: Optional[int] = None
        self.pcr2_rv_wt_seq: Optional[str] = None
        self.pcr2_rv_wt_alw: Optional[int] = None
        self.pcr2_rv_mut_seq: Optional[str] = None
        self.pcr2_rv_mut_alw: Optional[int] = None

        # Matrix header registry
        self.matrix_headers: Dict[str, List[str]] = {}

        # Whitelist/cache state
        self.whitelist_barcodes: List[str] = []
        self.whitelist_set: set[str] = set()
        self.whitelist_index: Dict[str, int] = {}
        self.whitelist_alphabet: List[str] = []
        self.barcode_hamm_cache: Dict[str, List[str]] = {}

        self.phred_scores = self.load_phred_scores()

    # ------------------------------------------------------------------
    # CLI helpers
    # ------------------------------------------------------------------

    @staticmethod
    def parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            prog="IronThrone-GoT",
            formatter_class=argparse.RawTextHelpFormatter,
            description=WHAT_IS_THIS_CODE,
        )
        parser.add_argument("-r", "--run", default="linear")
        parser.add_argument("-f1", "--fastqR1", required=True)
        parser.add_argument("-f2", "--fastqR2", required=True)
        parser.add_argument("-c", "--config", required=True)
        parser.add_argument("-w", "--whitelist", default=str(Path(__file__).resolve().parent / "barcodes10X" / "737K-august-2016.txt"))
        parser.add_argument("-m", "--mmtch", default=0.2)
        parser.add_argument("-p", "--postP", default=0.99)
        parser.add_argument("-d", "--dupcut", default=1)
        parser.add_argument("-b", "--bclen", default=16)
        parser.add_argument("-u", "--umilen", default=10)
        parser.add_argument("-s", "--sample", default="myGoT")
        parser.add_argument("-o", "--outdir", default="./out")
        parser.add_argument("-l", "--log", default="myGoT.log")
        parser.add_argument("-t", "--thread", default=4)
        parser.add_argument("-k", "--keepouts", default=0)
        # Keep typo-compatible option from Perl and also support expected spelling.
        parser.add_argument("--vervose", dest="verbose", default=0)
        parser.add_argument("-v", "--verbose", dest="verbose", default=0)
        return parser

    def close(self) -> None:
        try:
            self._log_handle.close()
        except Exception:
            pass

    def log(self, msg: str) -> None:
        now = time.ctime()
        self._log_handle.write(f"[{now}] {msg}\n")
        self._log_handle.flush()

    # ------------------------------------------------------------------
    # Core pipeline
    # ------------------------------------------------------------------

    def run(self) -> None:
        self.pre_process()
        self.post_process()

    def pre_process(self) -> None:
        self.fastq_r1, self.fastq_r2 = self.un_gzip(self.fastq_r1, self.fastq_r2)

        barcode = self.out_dir / f"{self.samplename}.BARCODE"
        barcode_q = self.out_dir / f"{self.samplename}.BARCODE_Q"
        umi = self.out_dir / f"{self.samplename}.UMI"
        seq = self.out_dir / f"{self.samplename}.SEQ"
        seq_q = self.out_dir / f"{self.samplename}.SEQ_Q"

        self.sed_repeat(self.fastq_r1, str(barcode), 2, 4, 1, self.len_barcode)
        self.sed_repeat(self.fastq_r1, str(barcode_q), 4, 4, 1, self.len_barcode)
        self.sed_repeat(self.fastq_r1, str(umi), 2, 4, self.len_barcode + 1, self.len_barcode + self.len_umi)
        self.sed_repeat(self.fastq_r2, str(seq), 2, 4)
        self.sed_repeat(self.fastq_r2, str(seq_q), 4, 4)

        self.convert_phred_bq(str(barcode_q), str(self.out_dir / f"{self.samplename}.BARCODE_Qnum"))
        self.convert_numeric_bq(str(seq_q), str(self.out_dir / f"{self.samplename}.SEQ_Qnum"))
        self.sum_of_bq(str(seq_q))

        self.paste_cols(
            str(barcode),
            str(self.out_dir / f"{self.samplename}.BARCODE_Qnum"),
            str(umi),
            str(seq),
            str(self.out_dir / f"{self.samplename}.SEQ_Qnum"),
            str(self.out_dir / f"{self.samplename}.SEQ_Q.P_errors_SUM"),
            str(self.out_dir / f"{self.samplename}.SEQUENCE"),
        )

        if self.run_mode == "linear":
            self.seqs_looking(str(self.out_dir / f"{self.samplename}.SEQUENCE"))
        elif self.run_mode == "circ":
            self.seqs_looking2(str(self.out_dir / f"{self.samplename}.SEQUENCE"))
        else:
            raise ValueError(f"Unsupported run mode: {self.run_mode}")

    def post_process(self) -> None:
        conf_lines = [line.rstrip("\n\r") for line in Path(self.conf_file).read_text(encoding="utf-8").splitlines()]

        if self.run_mode == "linear":
            self.primer_seq = conf_lines[0].split("\t")[0]
            self.primer_alw = _round_perl(len(self.primer_seq) * self.mmtch)

            general_parts = conf_lines[1].split("\t")
            self.general_seq = general_parts[0]
            self.general_alw = _round_perl(len(self.general_seq) * self.mmtch)
            self.general_start = self._to_int_or_none(general_parts[1] if len(general_parts) > 1 else None)
            self.general_end = self._to_int_or_none(general_parts[2] if len(general_parts) > 2 else None)

            self.wt_seq = conf_lines[2].split("\t")[0]
            self.wt_alw = _round_perl(len(self.wt_seq) * self.mmtch)

        elif self.run_mode == "circ":
            self.primer_seq = conf_lines[0].split("\t")[0]
            self.primer_alw = _round_perl(len(self.primer_seq) * self.mmtch)

            general_parts = conf_lines[1].split("\t")
            self.general_seq = general_parts[0]
            self.general_alw = _round_perl(len(self.general_seq) * self.mmtch)
            self.general_start = self._to_int_or_none(general_parts[1] if len(general_parts) > 1 else None)
            self.general_end = self._to_int_or_none(general_parts[2] if len(general_parts) > 2 else None)

            self.wt_seq = conf_lines[2].split("\t")[0]
            self.wt_alw = _round_perl(len(self.wt_seq) * self.mmtch)

            self.pcr2_fw_wt_seq = conf_lines[4].split("\t")[0]
            self.pcr2_fw_wt_alw = _round_perl(len(self.pcr2_fw_wt_seq) * self.allow_err)

            self.pcr2_rv_wt_seq = conf_lines[6].split("\t")[0]
            self.pcr2_rv_wt_alw = _round_perl(len(self.pcr2_rv_wt_seq) * self.allow_err)

        self.whitelist_barcodes = self.file_read(self.whitelist_file)
        self.whitelist_barcodes = [x.rstrip("\n\r") for x in self.whitelist_barcodes if x.strip()]
        self.whitelist_set = set(self.whitelist_barcodes)
        self.whitelist_index = {bc: i for i, bc in enumerate(self.whitelist_barcodes)}
        alphabet = set()
        for bc in self.whitelist_barcodes:
            alphabet.update(bc)
        self.whitelist_alphabet = sorted(alphabet) if alphabet else ["A", "C", "G", "T", "N"]

        looked_path = self.out_dir / f"{self.samplename}.looked"
        file_count_looked = 0
        primed_looked_done: List[str] = []
        primed_looked_nd: List[str] = []
        filter_stats = {
            "total": 0,
            "primer_pass": 0,
            "shared_pass": 0,
            "both_pass": 0,
            "primer_min": float("inf"),
            "primer_max": float("-inf"),
            "shared_min": float("inf"),
            "shared_max": float("-inf"),
        }

        if self.run_mode == "linear":
            if self.verbose:
                self.log(
                    f"seeking for mismatches | Primer<= {self.primer_alw} | Shared<= {self.general_alw} "
                )
            with looked_path.open("r", encoding="utf-8") as looked_handle:
                for raw_line in looked_handle:
                    line = raw_line.rstrip("\n\r")
                    file_count_looked += 1
                    cols = line.split("\t")
                    if len(cols) < 5:
                        continue
                    mm_primer = self._safe_float(cols[3], None)
                    mm_share = self._safe_float(cols[4], None)
                    if mm_primer is None or mm_share is None:
                        continue

                    barcode = cols[0]
                    filter_stats["total"] += 1
                    filter_stats["primer_min"] = min(filter_stats["primer_min"], mm_primer)
                    filter_stats["primer_max"] = max(filter_stats["primer_max"], mm_primer)
                    filter_stats["shared_min"] = min(filter_stats["shared_min"], mm_share)
                    filter_stats["shared_max"] = max(filter_stats["shared_max"], mm_share)

                    primer_ok = mm_primer <= float(self.primer_alw)
                    shared_ok = mm_share <= float(self.general_alw)
                    if primer_ok:
                        filter_stats["primer_pass"] += 1
                    if shared_ok:
                        filter_stats["shared_pass"] += 1

                    if primer_ok and shared_ok:
                        filter_stats["both_pass"] += 1
                        if barcode in self.whitelist_set:
                            primed_looked_done.append(line)
                        else:
                            primed_looked_nd.append(line)
            if self.verbose:
                self.log(f"   -> {filter_stats['both_pass']} found from {file_count_looked} lines")

        else:
            if self.verbose:
                self.log(
                    "seeking for mismatches | "
                    f"Primer<= {self.primer_alw} | Shared<= {self.general_alw} | "
                    f"PCR2FW_WT/PCR2FW_MUT<= {self.pcr2_fw_wt_alw} | PCR2RV_WT/PCR2RV_MUT<= {self.pcr2_rv_wt_alw} "
                )

            with looked_path.open("r", encoding="utf-8") as looked_handle:
                for raw_line in looked_handle:
                    line = raw_line.rstrip("\n\r")
                    file_count_looked += 1
                    cols = line.split("\t")
                    if len(cols) < 11:
                        continue
                    mm_primer = self._safe_float(cols[3], None)
                    mm_share = self._safe_float(cols[4], None)
                    mm_pcr2fw_wt = self._safe_float(cols[7], None)
                    mm_pcr2fw_mut = self._safe_float(cols[8], None)
                    mm_pcr2rv_wt = self._safe_float(cols[9], None)
                    mm_pcr2rv_mut = self._safe_float(cols[10], None)
                    if (
                        mm_primer is None
                        or mm_share is None
                        or mm_pcr2fw_wt is None
                        or mm_pcr2fw_mut is None
                        or mm_pcr2rv_wt is None
                        or mm_pcr2rv_mut is None
                    ):
                        continue

                    barcode = cols[0]
                    filter_stats["total"] += 1
                    filter_stats["primer_min"] = min(filter_stats["primer_min"], mm_primer)
                    filter_stats["primer_max"] = max(filter_stats["primer_max"], mm_primer)
                    filter_stats["shared_min"] = min(filter_stats["shared_min"], mm_share)
                    filter_stats["shared_max"] = max(filter_stats["shared_max"], mm_share)

                    primer_ok = mm_primer <= float(self.primer_alw)
                    shared_ok = mm_share <= float(self.general_alw)
                    if primer_ok:
                        filter_stats["primer_pass"] += 1
                    if shared_ok:
                        filter_stats["shared_pass"] += 1

                    if (
                        primer_ok
                        and shared_ok
                        and (
                            mm_pcr2fw_wt <= float(self.pcr2_fw_wt_alw)
                            or mm_pcr2fw_mut <= float(self.pcr2_fw_wt_alw)
                        )
                        and (
                            mm_pcr2rv_wt <= float(self.pcr2_rv_wt_alw)
                            or mm_pcr2rv_mut <= float(self.pcr2_rv_wt_alw)
                        )
                    ):
                        filter_stats["both_pass"] += 1
                        if barcode in self.whitelist_set:
                            primed_looked_done.append(line)
                        else:
                            primed_looked_nd.append(line)

            if self.verbose:
                self.log(f"   -> {filter_stats['both_pass']} found from {file_count_looked} lines")

        stats_msg = (
            "filter stats: "
            f"total={filter_stats['total']}, "
            f"primer_pass={filter_stats['primer_pass']}, "
            f"shared_pass={filter_stats['shared_pass']}, "
            f"both_pass={filter_stats['both_pass']}, "
            f"primer_min={self._fmt_stat(filter_stats['primer_min'])}, "
            f"primer_max={self._fmt_stat(filter_stats['primer_max'])}, "
            f"shared_min={self._fmt_stat(filter_stats['shared_min'])}, "
            f"shared_max={self._fmt_stat(filter_stats['shared_max'])}"
        )
        self.log(stats_msg)

        if filter_stats["both_pass"] == 0:
            self._handle_no_reads_after_filter(filter_stats)
            return

        mtx_done = self.array_to_matrix(primed_looked_done)
        mtx_nd = self.array_to_matrix(primed_looked_nd)

        self.matrix_header_create("MTX_done", mtx_done[0] if mtx_done else None)
        self.matrix_header_create("MTX_ND", mtx_nd[0] if mtx_nd else None)

        self.matrix_put_values(mtx_done, "MTX_done", "whitelist", "Y")
        self.matrix_put_values(mtx_nd, "MTX_ND", "whitelist", "N")

        mtx_sub = [row[:] for row in mtx_nd]
        self.matrix_header_create("MTX_sub", None, *self.matrix_header_get("MTX_ND"))
        self.matrix_header_rename("MTX_sub", "V1", "original_BC")
        self.matrix_put_values(mtx_sub, "MTX_sub", "replaced_BC", "NA")

        self.stringdist_hamm_pre_calc(mtx_sub)

        done_rows_by_bc = Counter(self.matrix_get_array_by_col(mtx_done, 0))
        tested: List[List[str]] = []

        for r, row in enumerate(mtx_sub):
            ori_bc = row[0]
            bq = self._split_semicolon(row[1])

            hamm_1 = self.barcode_hamm_cache.get(ori_bc, [])

            replaced_bc_idx = self.matrix_header_get_index("MTX_sub", "replaced_BC")

            if len(hamm_1) < 1:
                row[replaced_bc_idx] = "No_Hamm_1"
            else:
                dat_n_header = ["class", "bc", "qv", "At", "readN", "priorP", "p_edit", "likelihood", "postP"]
                self.matrix_header_create("dat_n", None, *dat_n_header)
                self.matrix_header_create("maxdat_n", None, *dat_n_header)

                dat_n: List[List[float | str]] = []

                for i, hamm_bc in enumerate(hamm_1):
                    diff_pos = self.str_diff_pos_1(ori_bc, hamm_bc)
                    nin_mtx_done_i = done_rows_by_bc.get(hamm_bc, 0)

                    if nin_mtx_done_i == 0:
                        row[replaced_bc_idx] = "NO_in_SAVED"
                        dat_n.append([i, hamm_bc, 0, diff_pos, 0])
                    else:
                        errs_n_i = float(bq[diff_pos]) if diff_pos < len(bq) and bq[diff_pos] != "" else 0.0
                        errs_edit_n_i = self.min_value(33, errs_n_i)
                        dat_n.append([i, hamm_bc, errs_edit_n_i, diff_pos, nin_mtx_done_i])

                self.matrix_calc_prior_p(dat_n)
                self.matrix_calc_p_edit(dat_n)
                self.matrix_calc_likelihood(dat_n)
                self.matrix_calc_post_p(dat_n)

                maxdat_n = self.matrix_calc_maxdat_n(dat_n, "dat_n", "postP")

                if self.matrix_sum_cols(dat_n, self.matrix_header_get_index("dat_n", "readN")) == 0 and len(maxdat_n) >= 2:
                    maxdat_n = [maxdat_n[0]]
                    row[replaced_bc_idx] = "NO_in_SAVED"

                if maxdat_n and maxdat_n[0][self.matrix_header_get_index("maxdat_n", "postP")] != "NA":
                    top_post = float(maxdat_n[0][self.matrix_header_get_index("maxdat_n", "postP")])
                    if top_post >= self.post_p:
                        row[replaced_bc_idx] = str(maxdat_n[0][self.matrix_header_get_index("maxdat_n", "bc")])
                    else:
                        row[replaced_bc_idx] = "Not_Sig"

                if len(maxdat_n) == 2:
                    top_post = float(maxdat_n[0][self.matrix_header_get_index("maxdat_n", "postP")])
                    if top_post == 0.5:
                        row[replaced_bc_idx] = "Not_Sig"
                        maxdat_n = [maxdat_n[0]]

            if self.verbose:
                self.log(f"processing row {r}")

            tested.append(row)

        mtx_changed = tested
        self.matrix_header_create("MTX_changed", None, *self.matrix_header_get("MTX_sub"))
        if self.keep_outs:
            self.matrix_save(mtx_changed, str(self.out_dir / f"{self.samplename}.MTX_changed.txt"), "MTX_changed")

        mtx_changed_extracted = self.matrix_calc_mtx_extract(mtx_changed)
        self.matrix_header_create("MTX_changed_extracted", None, *self.matrix_header_get("MTX_sub"))
        self.matrix_calc_mtx_changed_extracted_put_bc(mtx_changed_extracted)

        if self.run_mode == "linear":
            mtx_changed_extracted = self.matrix_remove_col(mtx_changed_extracted, "MTX_changed_extracted", 13)
            self.matrix_header_rename("MTX_changed_extracted", "original_BC", "V1")
            self.matrix_header_rename("MTX_changed_extracted", "V13", "whitelist")

            mtx_re = mtx_done + mtx_changed_extracted
            self.matrix_header_create("MTX_RE", None, *self.matrix_header_get("MTX_changed_extracted"))
            self.matrix_header_rename("MTX_RE", "V6", "mismatch.WT")
            self.matrix_header_rename("MTX_RE", "V7", "mismatch.MUT")
            self.matrix_header_rename("MTX_RE", "V11", "BQ.in.WT")
            self.matrix_header_rename("MTX_RE", "V12", "BQ.in.MUT")
        else:
            mtx_changed_extracted = self.matrix_remove_col(mtx_changed_extracted, "MTX_changed_extracted", 21)
            self.matrix_header_rename("MTX_changed_extracted", "original_BC", "V1")
            self.matrix_header_rename("MTX_changed_extracted", "V21", "whitelist")

            mtx_re = mtx_done + mtx_changed_extracted
            self.matrix_header_create("MTX_RE", None, *self.matrix_header_get("MTX_changed_extracted"))
            self.matrix_header_rename("MTX_RE", "V6", "mismatch.WT")
            self.matrix_header_rename("MTX_RE", "V7", "mismatch.MUT")
            self.matrix_header_rename("MTX_RE", "V15", "BQ.in.WT")
            self.matrix_header_rename("MTX_RE", "V16", "BQ.in.MUT")

        self.matrix_header_insert("MTX_RE", "BARCODE_UMI")
        self.matrix_calc_mtx_re_barcode_umi(mtx_re, "MTX_RE")

        mtx_re_barcode_umi_list = self.matrix_get_array_by_col(mtx_re, self.matrix_header_get_index("MTX_RE", "BARCODE_UMI"))
        mtx_re_dedup: List[List[str]] = []
        barcode_umi_groups = self.group_indices_by_value(mtx_re_barcode_umi_list)

        for d, idx_rows in enumerate(barcode_umi_groups.values()):
            mtx_re_d = [mtx_re[i] for i in idx_rows]
            nrow_mtx_re_d = len(mtx_re_d)

            self.matrix_header_create("MTX_RE_d", None, *self.matrix_header_get("MTX_RE"))
            self.matrix_put_values(mtx_re_d, "MTX_RE_d", "total_dups", nrow_mtx_re_d)
            self.matrix_put_values(mtx_re_d, "MTX_RE_d", "mismatch.WT.call", "")
            self.matrix_put_values(mtx_re_d, "MTX_RE_d", "mismatch.MUT.call", "")
            self.matrix_put_values(mtx_re_d, "MTX_RE_d", "mismatch.call", "amb")

            self.matrix_calc_mtx_re_d_mismatch_wt_assign(mtx_re_d, "MTX_RE_d")
            self.matrix_calc_mtx_re_d_mismatch_mut_assign(mtx_re_d, "MTX_RE_d")
            self.matrix_calc_mtx_re_d_mismatch_wt_call_assign(mtx_re_d, "MTX_RE_d")
            self.matrix_calc_mtx_re_d_mismatch_mut_call_assign(mtx_re_d, "MTX_RE_d")
            self.matrix_calc_mtx_re_d_mismatch_amb_assign(mtx_re_d, "MTX_RE_d")

            self.matrix_put_values(
                mtx_re_d,
                "MTX_RE_d",
                "WT.inDups",
                len(self.find_in_matrix("in", mtx_re_d, "MTX_RE_d", "mismatch.call", "WT")),
            )
            self.matrix_put_values(
                mtx_re_d,
                "MTX_RE_d",
                "MUT.inDups",
                len(self.find_in_matrix("in", mtx_re_d, "MTX_RE_d", "mismatch.call", "MUT")),
            )
            self.matrix_put_values(
                mtx_re_d,
                "MTX_RE_d",
                "amb.inDups",
                len(self.find_in_matrix("in", mtx_re_d, "MTX_RE_d", "mismatch.call", "amb")),
            )

            self.matrix_header_create("MTX_wrap_d", None, *self.matrix_header_get("MTX_RE_d"))
            self.matrix_header_create("MTX_wrap_d_BEST", None, *self.matrix_header_get("MTX_RE_d"))

            wt_in_dups = float(mtx_re_d[0][self.matrix_header_get_index("MTX_RE_d", "WT.inDups")]) if mtx_re_d else 0
            mut_in_dups = float(mtx_re_d[0][self.matrix_header_get_index("MTX_RE_d", "MUT.inDups")]) if mtx_re_d else 0

            mtx_wrap_d_best: List[List[str]] = []

            if wt_in_dups > mut_in_dups:
                mtx_wrap_d = [mtx_re_d[i] for i in self.find_in_matrix("in", mtx_re_d, "MTX_RE_d", "mismatch.call", "WT")]
                min_bq_wt = self.matrix_stat_col("min", mtx_re_d, self.matrix_header_get_index("MTX_RE_d", "BQ.in.WT"))
                mtx_wrap_d_best = [
                    mtx_re_d[i]
                    for i in self.find_in_matrix("in", mtx_re_d, "MTX_RE_d", "BQ.in.WT", min_bq_wt)
                ]
                if len(mtx_wrap_d_best) > 1 and mtx_wrap_d_best:
                    best_idx = self.matrix_calc_best_quality(
                        mtx_wrap_d_best, self.matrix_header_get_index("MTX_wrap_d_BEST", "V2")
                    )
                    mtx_wrap_d_best = [mtx_wrap_d_best[best_idx]]

            if wt_in_dups < mut_in_dups:
                mtx_wrap_d = [mtx_re_d[i] for i in self.find_in_matrix("in", mtx_re_d, "MTX_RE_d", "mismatch.call", "MUT")]
                min_bq_wt = self.matrix_stat_col("min", mtx_wrap_d, self.matrix_header_get_index("MTX_wrap_d", "BQ.in.WT"))
                mtx_wrap_d_best = [
                    mtx_wrap_d[i]
                    for i in self.find_in_matrix("in", mtx_wrap_d, "MTX_wrap_d", "BQ.in.WT", min_bq_wt)
                ]
                if len(mtx_wrap_d_best) > 1 and mtx_wrap_d_best:
                    best_idx = self.matrix_calc_best_quality(
                        mtx_wrap_d_best, self.matrix_header_get_index("MTX_wrap_d_BEST", "V2")
                    )
                    mtx_wrap_d_best = [mtx_wrap_d_best[best_idx]]

            if wt_in_dups == mut_in_dups:
                min_bq_wt = self.matrix_stat_col("min", mtx_re_d, self.matrix_header_get_index("MTX_RE_d", "BQ.in.WT"))
                mtx_wrap_d_best = [
                    mtx_re_d[i]
                    for i in self.find_in_matrix("in", mtx_re_d, "MTX_RE_d", "BQ.in.WT", min_bq_wt)
                ]
                if len(mtx_wrap_d_best) > 1 and mtx_wrap_d_best:
                    best_idx = self.matrix_calc_best_quality(
                        mtx_wrap_d_best, self.matrix_header_get_index("MTX_wrap_d_BEST", "V2")
                    )
                    mtx_wrap_d_best = [mtx_wrap_d_best[best_idx]]

            if self.run_mode == "linear":
                mtx_wrap_d_best = self.matrix_remove_col(mtx_wrap_d_best, "MTX_wrap_d_BEST", 16)
                mtx_wrap_d_best = self.matrix_remove_col(mtx_wrap_d_best, "MTX_wrap_d_BEST", 15)
                mtx_wrap_d_best = self.matrix_remove_col(mtx_wrap_d_best, "MTX_wrap_d_BEST", 14)
            else:
                mtx_wrap_d_best = self.matrix_remove_col(mtx_wrap_d_best, "MTX_wrap_d_BEST", 25)
                mtx_wrap_d_best = self.matrix_remove_col(mtx_wrap_d_best, "MTX_wrap_d_BEST", 24)
                mtx_wrap_d_best = self.matrix_remove_col(mtx_wrap_d_best, "MTX_wrap_d_BEST", 23)

            if mtx_wrap_d_best:
                mtx_re_dedup.extend(mtx_wrap_d_best)

            if self.verbose:
                self.log(f"processing {d + 1} out of {len(barcode_umi_groups)}")

        self.matrix_header_create("MTX_RE_dedup", None, *self.matrix_header_get("MTX_wrap_d_BEST"))

        if self.keep_outs:
            self.matrix_save(mtx_re_dedup, str(self.out_dir / f"{self.samplename}.MTX_RE.dedup.txt"), "MTX_RE_dedup")

        mtx_merge = [row[:] for row in mtx_re_dedup]
        self.matrix_header_create("MTX_MERGE", None, *self.matrix_header_get("MTX_RE_dedup"))

        if self.run_mode == "linear":
            mtx_merge = self.matrix_remove_col(mtx_merge, "MTX_MERGE", 13)
            col_select = [1, 13, 3, 15, 16, 17, 14, 8, 9, 10, 11, 12, 4, 5, 6, 7]
            col_select = [x - 1 for x in col_select]
            mtx_merge = self.matrix_get_col(mtx_merge, "MTX_MERGE", *col_select)
            self.matrix_header_create(
                "MTX_MERGE",
                None,
                "BC",
                "whitelist",
                "UMI",
                "num.WT.in.dups",
                "num.MUT.in.dups",
                "num.amb.in.dups",
                "call.in.dups",
                "avg.base_error.R2",
                "avg.base_error.primer",
                "avg.base_error.shared",
                "avg.base_error.WT",
                "avg.base_error.MUT",
                "mismatch.primer",
                "mismatch.shared",
                "mismatch.WT",
                "mismatch.MUT",
            )
        else:
            mtx_merge = self.matrix_remove_col(mtx_merge, "MTX_MERGE", 21)
            col_select = [1, 21, 3, 23, 24, 25, 22, 12, 13, 14, 15, 16, 4, 5, 6, 7, 8, 9, 10, 11]
            col_select = [x - 1 for x in col_select]
            mtx_merge = self.matrix_get_col(mtx_merge, "MTX_MERGE", *col_select)
            self.matrix_header_create(
                "MTX_MERGE",
                None,
                "BC",
                "whitelist",
                "UMI",
                "num.WT.in.dups",
                "num.MUT.in.dups",
                "num.amb.in.dups",
                "call.in.dups",
                "avg.base_error.R2",
                "avg.base_error.primer",
                "avg.base_error.shared",
                "avg.base_error.WT",
                "avg.base_error.MUT",
                "mismatch.primer",
                "mismatch.shared",
                "mismatch.WT",
                "mismatch.MUT",
                "mismatchPCR2FW_WT",
                "mismatchPCR2FW_MUT",
                "mismatchPCR2RV_WT",
                "mismatchPCR2RV_MUT",
            )

        if self.keep_outs:
            self.matrix_save(mtx_merge, str(self.out_dir / f"{self.samplename}.MTX_MERGE.txt"), "MTX_MERGE")

        summing: List[List[str]] = []
        mtx_merge_bc = self.matrix_get_array_by_col(mtx_merge, self.matrix_header_get_index("MTX_MERGE", "BC"))
        merge_bc_groups = self.group_indices_by_value(mtx_merge_bc)

        for idx_rows in merge_bc_groups.values():
            mtx_merge_u = [mtx_merge[i] for i in idx_rows]
            self.matrix_header_create("MTX_MERGE_u", None, *self.matrix_header_get("MTX_MERGE"))

            self.matrix_put_values(mtx_merge_u, "MTX_MERGE_u", "mismatch.WT.call", "")
            self.matrix_put_values(mtx_merge_u, "MTX_MERGE_u", "mismatch.MUT.call", "")
            self.matrix_put_values(mtx_merge_u, "MTX_MERGE_u", "mismatch.call", "")

            self.matrix_calc_mtx_merge_u_wt_assign(mtx_merge_u, "MTX_MERGE_u")
            self.matrix_calc_mtx_merge_u_mut_assign(mtx_merge_u, "MTX_MERGE_u")
            self.matrix_calc_mtx_merge_u_wt_call_assign(mtx_merge_u, "MTX_MERGE_u")
            self.matrix_calc_mtx_merge_u_mut_call_assign(mtx_merge_u, "MTX_MERGE_u")
            self.matrix_calc_mtx_merge_u_amb_assign(mtx_merge_u, "MTX_MERGE_u")

            mtx_merge_u = self.matrix_remove_col(
                mtx_merge_u, "MTX_MERGE_u", self.matrix_header_get_index("MTX_MERGE_u", "mismatch.MUT.call")
            )
            mtx_merge_u = self.matrix_remove_col(
                mtx_merge_u, "MTX_MERGE_u", self.matrix_header_get_index("MTX_MERGE_u", "mismatch.WT.call")
            )

            if not mtx_merge_u:
                continue

            arr_ref0_cnt = len(mtx_merge_u[0])
            arr_ref0_tmp = ["0"] * (arr_ref0_cnt - 1 + 3)
            summing_u = [arr_ref0_tmp]

            summing_u_hd = self.matrix_header_get("MTX_MERGE_u")
            summing_u_hd = summing_u_hd[:-1]
            summing_u_hd = summing_u_hd + ["WT.calls", "MUT.calls", "amb.calls"]
            self.matrix_header_create("summing_u", None, *summing_u_hd)

            summing_u[0][0] = mtx_merge_u[0][0]

            cnt_wt, cnt_mut, cnt_amb = self.matrix_calc_get_mismatch_call(mtx_merge_u, "MTX_MERGE_u")
            self.matrix_put_values(summing_u, "summing_u", "WT.calls", cnt_wt)
            self.matrix_put_values(summing_u, "summing_u", "MUT.calls", cnt_mut)
            self.matrix_put_values(summing_u, "summing_u", "amb.calls", cnt_amb)

            for c in range(self.matrix_header_get_ncol("MTX_MERGE_u") - 1):
                tmp = ";".join(str(x) for x in self.matrix_get_array_by_col(mtx_merge_u, c))
                self.matrix_put_values(summing_u, "summing_u", self.matrix_header_get_name("summing_u", c), tmp)

            summing.extend(summing_u)

        final_sum: List[List[str]] = []
        for line_col in summing:
            dup_wt = self.sum_value(self._split_semicolon(line_col[3]))
            dup_mut = self.sum_value(self._split_semicolon(line_col[4]))
            dup_amb = self.sum_value(self._split_semicolon(line_col[5]))
            if dup_wt + dup_mut + dup_amb >= self.dupcut:
                final_sum.append(line_col)

        if not self.keep_outs:
            for p in self.out_dir.glob(f"{self.samplename}.*"):
                try:
                    p.unlink()
                except FileNotFoundError:
                    pass

        self.matrix_save(final_sum, str(self.out_dir / f"{self.samplename}.summTable.txt"), "summing_u")
        print("______________________________END____________________________")

    @staticmethod
    def _to_int_or_none(value) -> Optional[int]:
        try:
            if value is None or str(value).strip() == "":
                return None
            return int(value)
        except Exception:
            return None

    @staticmethod
    def _safe_float(value, default: Optional[float] = 0.0) -> Optional[float]:
        try:
            return float(value)
        except Exception:
            return default

    def _row_float(self, row: List, idx: int) -> Optional[float]:
        if idx < 0 or idx >= len(row):
            return None
        return self._safe_float(row[idx], None)

    @staticmethod
    def _fmt_stat(value: float) -> str:
        if value == float("inf") or value == float("-inf"):
            return "NA"
        return str(value)

    def _summtable_headers(self) -> List[str]:
        if self.run_mode == "circ":
            return [
                "BC",
                "whitelist",
                "UMI",
                "num.WT.in.dups",
                "num.MUT.in.dups",
                "num.amb.in.dups",
                "call.in.dups",
                "avg.base_error.R2",
                "avg.base_error.primer",
                "avg.base_error.shared",
                "avg.base_error.WT",
                "avg.base_error.MUT",
                "mismatch.primer",
                "mismatch.shared",
                "mismatch.WT",
                "mismatch.MUT",
                "mismatchPCR2FW_WT",
                "mismatchPCR2FW_MUT",
                "mismatchPCR2RV_WT",
                "mismatchPCR2RV_MUT",
                "WT.calls",
                "MUT.calls",
                "amb.calls",
            ]
        return [
            "BC",
            "whitelist",
            "UMI",
            "num.WT.in.dups",
            "num.MUT.in.dups",
            "num.amb.in.dups",
            "call.in.dups",
            "avg.base_error.R2",
            "avg.base_error.primer",
            "avg.base_error.shared",
            "avg.base_error.WT",
            "avg.base_error.MUT",
            "mismatch.primer",
            "mismatch.shared",
            "mismatch.WT",
            "mismatch.MUT",
            "WT.calls",
            "MUT.calls",
            "amb.calls",
        ]

    def _write_empty_summ_table(self, reason: str) -> None:
        if not self.keep_outs:
            for p in self.out_dir.glob(f"{self.samplename}.*"):
                try:
                    p.unlink()
                except FileNotFoundError:
                    pass
        self.matrix_header_create("summing_u", None, *self._summtable_headers())
        out_file = str(self.out_dir / f"{self.samplename}.summTable.txt")
        self.matrix_save([], out_file, "summing_u")
        self.log(f"Wrote empty summary table: {out_file} | reason: {reason}")
        print("______________________________END____________________________")

    @staticmethod
    def _best_mismatch_position(ref: str, seq: str) -> Tuple[int, int]:
        if not ref or not seq or len(seq) < len(ref):
            return len(ref), 1
        best_mm = len(ref) + 1
        best_pos = 1
        win = len(ref)
        for start in range(0, len(seq) - win + 1):
            mm = 0
            sub = seq[start : start + win]
            for a, b in zip(ref, sub):
                if a != b:
                    mm += 1
                    if mm >= best_mm:
                        break
            if mm < best_mm:
                best_mm = mm
                best_pos = start + 1
                if best_mm == 0:
                    break
        return best_mm, best_pos

    def _diagnose_general_match_positions(self, max_reads: int = 2000) -> Optional[Dict[str, object]]:
        if not self.general_seq:
            return None
        seq_path = self.out_dir / f"{self.samplename}.SEQUENCE"
        if not seq_path.exists():
            return None

        pos_counter: Counter[int] = Counter()
        mm_counter: Counter[int] = Counter()
        sampled = 0

        with seq_path.open("r", encoding="utf-8") as f:
            for line in f:
                cols = line.rstrip("\n\r").split("\t")
                if len(cols) < 4:
                    continue
                read_seq = cols[3]
                if len(read_seq) < len(self.general_seq):
                    continue
                best_mm, best_pos = self._best_mismatch_position(self.general_seq, read_seq)
                pos_counter[best_pos] += 1
                mm_counter[best_mm] += 1
                sampled += 1
                if sampled >= max_reads:
                    break

        if sampled == 0:
            return None
        return {
            "sampled": sampled,
            "top_positions": pos_counter.most_common(3),
            "top_mismatch": mm_counter.most_common(3),
        }

    def _handle_no_reads_after_filter(self, filter_stats: Dict[str, float]) -> None:
        message = (
            "No reads passed filter criteria: "
            f"Primer<={self.primer_alw}, Shared<={self.general_alw}. "
            f"(total={filter_stats['total']}, both_pass={filter_stats['both_pass']})"
        )
        self.log(message)
        print(message)

        diag = self._diagnose_general_match_positions(max_reads=2000)
        if diag is not None:
            pos_txt = ", ".join([f"pos{p}:{n}" for p, n in diag["top_positions"]])  # type: ignore[index]
            mm_txt = ", ".join([f"mm{m}:{n}" for m, n in diag["top_mismatch"]])  # type: ignore[index]
            diag_msg = f"General-sequence best-match sample (n={diag['sampled']}): positions [{pos_txt}] | mismatches [{mm_txt}]"
            self.log(diag_msg)

            if self.general_start is not None and diag["top_positions"]:  # type: ignore[index]
                top_pos = diag["top_positions"][0][0]  # type: ignore[index]
                shift = abs(int(top_pos) - int(self.general_start))
                if shift >= 8:
                    self.log(
                        "Potential config/read coordinate mismatch: "
                        f"configured shared start={self.general_start}, observed best-match position~{top_pos}."
                    )

        self._write_empty_summ_table("No reads passed primer/shared mismatch filters")

    # ------------------------------------------------------------------
    # Read/sequence helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _split_semicolon(value: str) -> List[str]:
        if value is None:
            return []
        parts = value.split(";")
        if parts and parts[-1] == "":
            parts = parts[:-1]
        return parts

    def find_mm_1(self, ref: Sequence[str], s: str, s_exclude: str) -> List[str]:
        out: List[str] = []
        for e in ref:
            if e == s_exclude:
                continue
            if len(e) != len(s):
                continue
            mm_cnt = sum(1 for a, b in zip(s, e) if a != b)
            if mm_cnt == 1:
                out.append(e)
        return out

    def sum_of_bq(self, file_in: str) -> None:
        p = Path(file_in)
        if not p.exists() or p.stat().st_size <= 0:
            if self.verbose:
                self.log(f"{file_in} does not exist or 0 byte!")
            sys.exit(0)

        out_path = Path(f"{file_in}.P_errors_SUM")
        with p.open("r", encoding="utf-8") as fin, out_path.open("w", encoding="utf-8") as fout:
            for line in fin:
                line = line.rstrip("\n\r")
                total = 0.0
                for base in line:
                    total += float(self.phred_scores.get(base, {"p_err": 0.0})["p_err"])
                fout.write(f"{total}\n")

    def convert_phred_bq(self, file_in: str, file_out: str) -> None:
        p = Path(file_in)
        if not p.exists() or p.stat().st_size <= 0:
            self.log(f"{file_in} does not exist or 0 byte!")
            sys.exit(0)

        with p.open("r", encoding="utf-8") as fin, Path(file_out).open("w", encoding="utf-8") as fout:
            for line in fin:
                line = line.rstrip("\n\r")
                vals = [str(int(self.phred_scores.get(base, {"q_score": 0})["q_score"])) for base in line]
                fout.write(";".join(vals) + ";\n")

    def convert_numeric_bq(self, file_in: str, file_out: str) -> None:
        p = Path(file_in)
        if not p.exists() or p.stat().st_size <= 0:
            self.log(f"{file_in} does not exist or 0 byte!")
            sys.exit(0)

        with p.open("r", encoding="utf-8") as fin, Path(file_out).open("w", encoding="utf-8") as fout:
            for line in fin:
                line = line.rstrip("\n\r")
                vals = [str(float(self.phred_scores.get(base, {"p_err": 0.0})["p_err"])) for base in line]
                fout.write(";".join(vals) + ";\n")

        self.log(f"[{file_out}] created")

    def seqs_looking(self, file_in: str) -> None:
        p = Path(file_in)
        if not p.exists() or p.stat().st_size <= 0:
            self.log(f"{file_in} does not exist or 0 byte!")
            sys.exit(0)

        conf_lines = [line.rstrip("\n\r") for line in Path(self.conf_file).read_text(encoding="utf-8").splitlines()]

        primer_seq, primer_start, primer_end = conf_lines[0].split("\t")
        primer_len = len(primer_seq)
        primer_alw = _round_perl(primer_len * self.allow_err)
        primer_mmtch = _round_perl(primer_len * self.mmtch)

        general_seq, general_start, general_end = conf_lines[1].split("\t")
        wt_seq, wt_start, wt_end = conf_lines[2].split("\t")
        mut_seq, mut_start, mut_end = conf_lines[3].split("\t")

        general_start = int(general_start)
        general_end = int(general_end)
        wt_start = int(wt_start)
        wt_end = int(wt_end)
        mut_start = int(mut_start)
        mut_end = int(mut_end)
        primer_start = int(primer_start)

        looked_file = self.out_dir / f"{self.samplename}.looked"

        with p.open("r", encoding="utf-8") as fin, looked_file.open("w", encoding="utf-8") as fout:
            for raw in fin:
                line_raw = raw.rstrip("\n\r")
                if not line_raw.strip():
                    continue

                line_col = line_raw.split("\t")
                line = line_col[3]

                primer_found: Optional[int] = None
                general_startz = general_start
                general_endz = general_end
                wt_startz = wt_start
                wt_endz = wt_end
                mut_startz = mut_start
                mut_endz = mut_end

                for loc in range(0, 4):
                    seq_in_line = line[(primer_start - 1 + loc) : (primer_start - 1 + loc) + len(primer_seq) + loc]
                    primer_found_check = self.mismatch_check(primer_seq, seq_in_line, primer_alw)
                    if primer_found_check <= primer_mmtch:
                        primer_found = primer_found_check
                        general_startz = general_start + loc
                        general_endz = general_end + loc
                        wt_startz = wt_start + loc
                        wt_endz = wt_end + loc
                        mut_startz = mut_start + loc
                        mut_endz = mut_end + loc
                        break

                if primer_found is None:
                    primer_found = primer_len

                general_found = self.mismatch_check(general_seq, line[general_startz - 1 : general_startz - 1 + len(general_seq)], 0)
                wt_found = self.mismatch_check(wt_seq, line[wt_startz - 1 : wt_startz - 1 + len(wt_seq)], 0)
                mut_found = self.mismatch_check(mut_seq, line[mut_startz - 1 : mut_startz - 1 + len(mut_seq)], 0)

                r2_bq_vals = [float(x) for x in self._split_semicolon(line_col[4]) if x != ""]
                if (
                    (not line)
                    or (not line_col[5])
                    or (len(general_seq) == 0)
                    or (sum(r2_bq_vals[general_startz - 1 : general_endz]) == 0)
                ):
                    sys.exit(0)

                avg_r2_bq = float(line_col[5]) / max(1, len(line))

                bq_primer = r2_bq_vals[0:primer_len]
                bq_primer_avg = sum(bq_primer) / max(1, primer_len)

                bq_general = r2_bq_vals[general_startz - 1 : general_endz]
                bq_general_avg = sum(bq_general) / max(1, len(general_seq))

                bq_wt = r2_bq_vals[wt_startz - 1 : wt_endz]
                bq_wt_avg = sum(bq_wt) / max(1, len(wt_seq))

                bq_mut = r2_bq_vals[mut_startz - 1 : mut_endz]
                bq_mut_avg = sum(bq_mut) / max(1, len(mut_seq))

                fout.write(
                    f"{line_col[0]}\t{line_col[1]}\t{line_col[2]}\t{primer_found}\t{general_found}\t"
                    f"{wt_found}\t{mut_found}\t{avg_r2_bq}\t{bq_primer_avg}\t{bq_general_avg}\t{bq_wt_avg}\t{bq_mut_avg}\n"
                )

    def seqs_looking2(self, file_in: str) -> None:
        p = Path(file_in)
        if not p.exists() or p.stat().st_size <= 0:
            self.log(f"{file_in} does not exist or 0 byte!")
            sys.exit(0)

        conf_lines = [line.rstrip("\n\r") for line in Path(self.conf_file).read_text(encoding="utf-8").splitlines()]

        primer_seq, primer_start, primer_end = conf_lines[0].split("\t")
        primer_len = len(primer_seq)
        primer_alw = _round_perl(primer_len * self.allow_err)
        primer_mmtch = _round_perl(primer_len * self.mmtch)

        general_seq, general_start, general_end = conf_lines[1].split("\t")
        wt_seq, wt_start, wt_end = conf_lines[2].split("\t")
        mut_seq, mut_start, mut_end = conf_lines[3].split("\t")
        pcr2_fw_wt_seq, pcr2_fw_wt_start, pcr2_fw_wt_end = conf_lines[4].split("\t")
        pcr2_fw_mut_seq, pcr2_fw_mut_start, pcr2_fw_mut_end = conf_lines[5].split("\t")
        pcr2_rv_wt_seq, pcr2_rv_wt_start, pcr2_rv_wt_end = conf_lines[6].split("\t")
        pcr2_rv_mut_seq, pcr2_rv_mut_start, pcr2_rv_mut_end = conf_lines[7].split("\t")

        primer_start = int(primer_start)
        general_start = int(general_start)
        general_end = int(general_end)
        wt_start = int(wt_start)
        wt_end = int(wt_end)
        mut_start = int(mut_start)
        mut_end = int(mut_end)
        pcr2_fw_wt_start = int(pcr2_fw_wt_start)
        pcr2_fw_wt_end = int(pcr2_fw_wt_end)
        pcr2_fw_mut_start = int(pcr2_fw_mut_start)
        pcr2_fw_mut_end = int(pcr2_fw_mut_end)
        pcr2_rv_wt_start = int(pcr2_rv_wt_start)
        pcr2_rv_wt_end = int(pcr2_rv_wt_end)
        pcr2_rv_mut_start = int(pcr2_rv_mut_start)
        pcr2_rv_mut_end = int(pcr2_rv_mut_end)

        looked_file = self.out_dir / f"{self.samplename}.looked"

        with p.open("r", encoding="utf-8") as fin, looked_file.open("w", encoding="utf-8") as fout:
            for raw in fin:
                line_raw = raw.rstrip("\n\r")
                if not line_raw.strip():
                    continue

                line_col = line_raw.split("\t")
                line = line_col[3]

                primer_found: Optional[int] = None

                general_startz, general_endz = general_start, general_end
                wt_startz, wt_endz = wt_start, wt_end
                mut_startz, mut_endz = mut_start, mut_end
                pcr2_fw_wt_startz, pcr2_fw_wt_endz = pcr2_fw_wt_start, pcr2_fw_wt_end
                pcr2_fw_mut_startz, pcr2_fw_mut_endz = pcr2_fw_mut_start, pcr2_fw_mut_end
                pcr2_rv_wt_startz, pcr2_rv_wt_endz = pcr2_rv_wt_start, pcr2_rv_wt_end
                pcr2_rv_mut_startz, pcr2_rv_mut_endz = pcr2_rv_mut_start, pcr2_rv_mut_end

                for loc in range(0, 4):
                    seq_in_line = line[(primer_start - 1 + loc) : (primer_start - 1 + loc) + len(primer_seq) + loc]
                    primer_found_check = self.mismatch_check(primer_seq, seq_in_line, primer_alw)
                    if primer_found_check <= primer_mmtch:
                        primer_found = primer_found_check
                        general_startz, general_endz = general_start + loc, general_end + loc
                        wt_startz, wt_endz = wt_start + loc, wt_end + loc
                        mut_startz, mut_endz = mut_start + loc, mut_end + loc
                        pcr2_fw_wt_startz, pcr2_fw_wt_endz = pcr2_fw_wt_start + loc, pcr2_fw_wt_end + loc
                        pcr2_fw_mut_startz, pcr2_fw_mut_endz = pcr2_fw_mut_start + loc, pcr2_fw_mut_end + loc
                        pcr2_rv_wt_startz, pcr2_rv_wt_endz = pcr2_rv_wt_start + loc, pcr2_rv_wt_end + loc
                        pcr2_rv_mut_startz, pcr2_rv_mut_endz = pcr2_rv_mut_start + loc, pcr2_rv_mut_end + loc
                        break

                if primer_found is None:
                    primer_found = primer_len

                general_found = self.mismatch_check(general_seq, line[general_startz - 1 : general_startz - 1 + len(general_seq)], 0)
                wt_found = self.mismatch_check(wt_seq, line[wt_startz - 1 : wt_startz - 1 + len(wt_seq)], 0)
                mut_found = self.mismatch_check(mut_seq, line[mut_startz - 1 : mut_startz - 1 + len(mut_seq)], 0)
                pcr2_fw_wt_found = self.mismatch_check(
                    pcr2_fw_wt_seq,
                    line[pcr2_fw_wt_startz - 1 : pcr2_fw_wt_startz - 1 + len(pcr2_fw_wt_seq)],
                    0,
                )
                pcr2_fw_mut_found = self.mismatch_check(
                    pcr2_fw_mut_seq,
                    line[pcr2_fw_mut_startz - 1 : pcr2_fw_mut_startz - 1 + len(pcr2_fw_mut_seq)],
                    0,
                )
                pcr2_rv_wt_found = self.mismatch_check(
                    pcr2_rv_wt_seq,
                    line[pcr2_rv_wt_startz - 1 : pcr2_rv_wt_startz - 1 + len(pcr2_rv_wt_seq)],
                    0,
                )
                pcr2_rv_mut_found = self.mismatch_check(
                    pcr2_rv_mut_seq,
                    line[pcr2_rv_mut_startz - 1 : pcr2_rv_mut_startz - 1 + len(pcr2_rv_mut_seq)],
                    0,
                )

                r2_bq_vals = [float(x) for x in self._split_semicolon(line_col[4]) if x != ""]
                if (
                    (not line)
                    or (not line_col[5])
                    or (len(general_seq) == 0)
                    or (sum(r2_bq_vals[general_startz - 1 : general_endz]) == 0)
                ):
                    sys.exit(0)

                avg_r2_bq = float(line_col[5]) / max(1, len(line))

                bq_primer = r2_bq_vals[0:primer_len]
                bq_primer_avg = sum(bq_primer) / max(1, primer_len)

                bq_general = r2_bq_vals[general_startz - 1 : general_endz]
                bq_general_avg = sum(bq_general) / max(1, len(general_seq))

                bq_wt = r2_bq_vals[wt_startz - 1 : wt_endz]
                bq_wt_avg = sum(bq_wt) / max(1, len(wt_seq))

                bq_mut = r2_bq_vals[mut_startz - 1 : mut_endz]
                bq_mut_avg = sum(bq_mut) / max(1, len(mut_seq))

                bq_pcr2_fw_wt = r2_bq_vals[pcr2_fw_wt_startz - 1 : pcr2_fw_wt_endz]
                bq_pcr2_fw_wt_avg = sum(bq_pcr2_fw_wt) / max(1, len(pcr2_fw_wt_seq))

                bq_pcr2_fw_mut = r2_bq_vals[pcr2_fw_mut_startz - 1 : pcr2_fw_mut_endz]
                bq_pcr2_fw_mut_avg = sum(bq_pcr2_fw_mut) / max(1, len(pcr2_fw_mut_seq))

                bq_pcr2_rv_wt = r2_bq_vals[pcr2_rv_wt_startz - 1 : pcr2_rv_wt_endz]
                bq_pcr2_rv_wt_avg = sum(bq_pcr2_rv_wt) / max(1, len(pcr2_rv_wt_seq))

                bq_pcr2_rv_mut = r2_bq_vals[pcr2_rv_mut_startz - 1 : pcr2_rv_mut_endz]
                bq_pcr2_rv_mut_avg = sum(bq_pcr2_rv_mut) / max(1, len(pcr2_rv_mut_seq))

                # Keep historical behavior: last field writes WT average twice.
                fout.write(
                    f"{line_col[0]}\t{line_col[1]}\t{line_col[2]}\t{primer_found}\t{general_found}\t{wt_found}\t"
                    f"{mut_found}\t{pcr2_fw_wt_found}\t{pcr2_fw_mut_found}\t{pcr2_rv_wt_found}\t{pcr2_rv_mut_found}\t"
                    f"{avg_r2_bq}\t{bq_primer_avg}\t{bq_general_avg}\t{bq_wt_avg}\t{bq_mut_avg}\t"
                    f"{bq_pcr2_fw_wt_avg}\t{bq_pcr2_fw_mut_avg}\t{bq_pcr2_rv_wt_avg}\t{bq_pcr2_rv_wt_avg}\n"
                )

    @staticmethod
    def mismatch_check(str_ref: str, str_test: str, mismatch_allow_count: int) -> int:
        mismatched = 0
        for i in range(len(str_ref)):
            letter_ref = str_ref[i : i + 1]
            letter_test = str_test[i : i + 1]
            if letter_ref != letter_test:
                mismatched += 1
        return mismatched

    def file_read(self, file: str) -> List[str]:
        p = Path(file)
        if p.suffix == ".gz":
            with gzip.open(p, "rt", encoding="utf-8") as f:
                data = f.readlines()
        else:
            with p.open("r", encoding="utf-8") as f:
                data = f.readlines()
        if self.verbose:
            self.log(f"File read: {file} (size: {p.stat().st_size},   lines: {len(data)})")
        return data

    def unique(self, arr: Sequence[str]) -> List[str]:
        seen = set()
        out: List[str] = []
        for e in arr:
            if e not in seen:
                out.append(e)
                seen.add(e)
        if self.verbose:
            self.log(f"Unique done: {len(arr)} -> {len(out)}")
        return out

    def intersect(self, arr1: Sequence[str], arr2: Sequence[str]) -> List[str]:
        seen = {e.rstrip("\n\r") for e in arr2}
        out: List[str] = []
        for e in arr1:
            ev = e.rstrip("\n\r")
            if ev in seen:
                out.append(ev)
        if self.verbose:
            self.log(f"Intersect done: input count({len(arr1)}, {len(arr2)}),  intersect count: {len(out)}")
        return out

    def find_in_array(self, mode: str, arr1: Sequence[str], arr2: Sequence[str]) -> List[int]:
        seen = {e.rstrip("\n\r") for e in arr2}
        out: List[int] = []
        for i, e in enumerate(arr1):
            ev = e.rstrip("\n\r")
            if mode == "in":
                if ev in seen:
                    out.append(i)
            elif mode == "not in":
                if ev not in seen:
                    out.append(i)
        if self.verbose:
            self.log(f"find {mode} done: input count({len(arr1)}, {len(arr2)}),  found count: {len(out)}")
        return out

    @staticmethod
    def group_indices_by_value(values: Sequence[str]) -> Dict[str, List[int]]:
        groups: Dict[str, List[int]] = {}
        for i, value in enumerate(values):
            if value not in groups:
                groups[value] = []
            groups[value].append(i)
        return groups

    @staticmethod
    def matrix_get_array_by_col(mtx: List[List], idx: int) -> List:
        out = []
        for row in mtx:
            if idx < len(row):
                out.append(row[idx])
            else:
                out.append("")
        return out

    @staticmethod
    def matrix_get_matrix_by_row(mtx: List[List], rows: Sequence[int]) -> List[List]:
        rows_set = set(rows)
        out: List[List] = []
        for i, row in enumerate(mtx):
            if i in rows_set:
                out.append(row)
        return out

    def matrix_get_col(self, mat: List[List], mat_name: str, *cols: int) -> List[List]:
        old_hd = self.matrix_header_get(mat_name)
        self.matrix_headers[mat_name] = [old_hd[idx] for idx in cols]
        new_mat: List[List] = []
        for row in mat:
            new_mat.append([row[idx] if idx < len(row) else "" for idx in cols])
        return new_mat

    def find_in_matrix(self, mode: str, mat: List[List], mat_name: str, col_name: str, find_value) -> List[int]:
        idx = self.matrix_header_get_index(mat_name, col_name)
        out: List[int] = []
        for i, row in enumerate(mat):
            value = row[idx] if idx < len(row) else ""
            if mode == "in":
                if str(value) == str(find_value):
                    out.append(i)
            elif mode == "not in":
                if str(value) != str(find_value):
                    out.append(i)
        if self.verbose:
            self.log(f"find {mode} done: input count({len(mat)}),  found count: {len(out)}")
        return out

    def matrix_header_get(self, matrix_name: str) -> List[str]:
        return list(self.matrix_headers.get(matrix_name, []))

    def matrix_header_create(self, matrix_name: str, matrix_first_row: Optional[Sequence], *names: str) -> None:
        if matrix_first_row is not None:
            matrix_first_row_cols = list(matrix_first_row)
        else:
            matrix_first_row_cols = []

        if not names:
            self.matrix_headers[matrix_name] = [f"V{i + 1}" for i in range(len(matrix_first_row_cols))]
        else:
            self.matrix_headers[matrix_name] = list(names)

        if self.verbose:
            self.log(f"Matrix({matrix_name}) header created: total {len(self.matrix_headers[matrix_name])} columns.")

    def matrix_header_insert(self, matrix_name: str, name: str, idx: Optional[int] = None) -> int:
        if matrix_name not in self.matrix_headers:
            self.matrix_headers[matrix_name] = []
        ref = self.matrix_headers[matrix_name]

        if idx is None:
            idx = len(ref)

        if idx >= len(ref):
            ref.extend([""] * (idx - len(ref) + 1))
            ref[idx] = name
        else:
            ref.insert(idx, name)

        if self.verbose:
            self.log(f"Matrix({matrix_name}) index added: index({idx}), name({name})")
        return idx

    def matrix_header_remove(self, matrix_name: str, idx: int) -> None:
        ref = self.matrix_headers.get(matrix_name, [])
        if 0 <= idx < len(ref):
            ref.pop(idx)
        if self.verbose:
            self.log(f"Matrix({matrix_name}) index removed: index({idx})")

    def matrix_header_rename(self, matrix_name: str, name_old: str, name_new: str) -> None:
        ref = self.matrix_headers.get(matrix_name, [])
        for i, name in enumerate(ref):
            if name == name_old:
                ref[i] = name_new
                break
        if self.verbose:
            self.log(f"Matrix({matrix_name}) header renamed: old({name_old}) new({name_new})")

    @staticmethod
    def array_to_matrix(arr: Sequence[str]) -> List[List[str]]:
        matrix: List[List[str]] = []
        for line in arr:
            line = line.rstrip("\n\r")
            if not line.strip():
                continue
            matrix.append(line.split("\t"))
        return matrix

    def matrix_header_get_index(self, mat_name: str, header_name: str) -> int:
        ref = self.matrix_headers.get(mat_name, [])
        for i, name in enumerate(ref):
            if name == header_name:
                return i
        return -1

    def matrix_header_require_index(self, mat_name: str, header_name: str) -> int:
        idx = self.matrix_header_get_index(mat_name, header_name)
        if idx < 0:
            raise ValueError(f"Missing required header '{header_name}' in matrix '{mat_name}'")
        return idx

    def matrix_header_get_name(self, mat_name: str, header_idx: int) -> str:
        ref = self.matrix_headers.get(mat_name, [])
        if 0 <= header_idx < len(ref):
            return ref[header_idx]
        return "-1"

    def matrix_header_get_ncol(self, mat_name: str) -> int:
        ref = self.matrix_headers.get(mat_name, [])
        return len(ref)

    def matrix_put_values(self, matrix: List[List], matrix_name: str, header_name: str, *values) -> None:
        idx_to_put = self.matrix_header_get_index(matrix_name, header_name)
        if idx_to_put == -1:
            idx_to_put = self.matrix_header_insert(matrix_name, header_name)

        m_row = len(matrix)
        if m_row == 0:
            m_row = len(values)
            for _ in range(m_row):
                matrix.append([])

        values_list = list(values)
        if len(values_list) < m_row and values_list:
            last_value = values_list[-1]
            values_list.extend([last_value] * (m_row - len(values_list)))
        elif len(values_list) == 0:
            values_list = [""] * m_row

        for i in range(m_row):
            row = matrix[i]
            if idx_to_put >= len(row):
                row.extend([""] * (idx_to_put - len(row) + 1))
            row[idx_to_put] = values_list[i]

    @staticmethod
    def matrix_sum_cols(mat: List[List], idx_col: int) -> float:
        s = 0.0
        for row in mat:
            try:
                s += float(row[idx_col])
            except Exception:
                s += 0.0
        return s

    @staticmethod
    def matrix_stat_col(mode: str, mat: List[List], idx_col: int):
        vals: List[float] = []
        for row in mat:
            try:
                vals.append(float(row[idx_col]))
            except Exception:
                vals.append(0.0)
        if not vals:
            if mode == "sum":
                return 0.0
            return 0.0
        if mode == "min":
            return min(vals)
        if mode == "max":
            return max(vals)
        if mode == "sum":
            return sum(vals)
        if mode == "mean":
            return sum(vals) / len(vals)
        return 0.0

    def matrix_calc_best_quality(self, mat: List[List], idx_col: int) -> int:
        max_sum: Optional[float] = None
        max_idx = 0
        for i, row in enumerate(mat):
            q_sum = self.sum_value(self._split_semicolon(str(row[idx_col]))) if idx_col < len(row) else 0.0
            if max_sum is None or q_sum > max_sum:
                max_sum = q_sum
                max_idx = i
        return max_idx

    def matrix_calc_prior_p(self, mat: List[List]) -> None:
        idx_prior = self.matrix_header_get_index("dat_n", "priorP")
        idx_read = self.matrix_header_get_index("dat_n", "readN")
        sum_read = self.matrix_sum_cols(mat, idx_read)

        for row in mat:
            while len(row) <= idx_prior:
                row.append("")
            read_val = float(row[idx_read]) if idx_read < len(row) else 0.0
            row[idx_prior] = (read_val / sum_read) if sum_read else 0.0

    def matrix_calc_p_edit(self, mat: List[List]) -> None:
        idx_p_edit = self.matrix_header_get_index("dat_n", "p_edit")
        idx_qv = self.matrix_header_get_index("dat_n", "qv")

        for row in mat:
            while len(row) <= idx_p_edit:
                row.append("")
            qv = float(row[idx_qv]) if idx_qv < len(row) else 0.0
            row[idx_p_edit] = 10 ** (-(qv / 10.0))

    def matrix_calc_likelihood(self, mat: List[List]) -> None:
        idx_like = self.matrix_header_get_index("dat_n", "likelihood")
        idx_prior = self.matrix_header_get_index("dat_n", "priorP")
        idx_p_edit = self.matrix_header_get_index("dat_n", "p_edit")

        for row in mat:
            while len(row) <= idx_like:
                row.append("")
            prior = float(row[idx_prior]) if idx_prior < len(row) else 0.0
            p_edit = float(row[idx_p_edit]) if idx_p_edit < len(row) else 0.0
            row[idx_like] = prior * p_edit

    def matrix_calc_post_p(self, mat: List[List]) -> None:
        idx_post = self.matrix_header_get_index("dat_n", "postP")
        idx_like = self.matrix_header_get_index("dat_n", "likelihood")
        sum_like = self.matrix_sum_cols(mat, idx_like)

        for row in mat:
            while len(row) <= idx_post:
                row.append("")
            like = float(row[idx_like]) if idx_like < len(row) else 0.0
            row[idx_post] = (like / sum_like) if sum_like else 0.0

    def matrix_calc_mtx_extract(self, mat: List[List]) -> List[List]:
        idx_replaced_bc = self.matrix_header_get_index("MTX_changed", "replaced_BC")
        out = []
        for row in mat:
            val = row[idx_replaced_bc] if idx_replaced_bc < len(row) else ""
            if val not in {"NO_in_SAVED", "Not_Sig", "No_Hamm_1"}:
                out.append(row)
        return out

    def matrix_calc_mtx_changed_extracted_put_bc(self, mat: List[List]) -> None:
        idx_original = self.matrix_header_get_index("MTX_changed_extracted", "original_BC")
        idx_replaced = self.matrix_header_get_index("MTX_changed_extracted", "replaced_BC")
        for row in mat:
            while len(row) <= max(idx_original, idx_replaced):
                row.append("")
            row[idx_original] = row[idx_replaced]

    def matrix_save(self, mtx: List[List], file: str, mat_name: Optional[str] = None) -> None:
        p = Path(file)
        with p.open("w", encoding="utf-8") as f:
            if mat_name is not None:
                headers = self.matrix_header_get(mat_name)
                if headers:
                    f.write("\t".join(str(x) for x in headers) + "\n")

            for row in mtx:
                f.write("\t".join(str(x) for x in row) + "\n")

        col_cnt = len(mtx[0]) if mtx else 0
        self.log(f"Matrix file saved: {file} ({len(mtx)} rows x {col_cnt} cols)")

    def stringdist_hamm(self, s: str, mm_max_cnt: int) -> List[str]:
        # Fast Hamming<=1 lookup using exact set membership; preserve whitelist order.
        if mm_max_cnt != 1:
            out = []
            for bc in self.whitelist_barcodes:
                if len(bc) != len(s):
                    continue
                mm = sum(1 for a, b in zip(s, bc) if a != b)
                if mm <= mm_max_cnt:
                    out.append(bc)
            return out

        candidates = set()
        if s in self.whitelist_set:
            candidates.add(s)

        s_chars = list(s)
        for i, orig in enumerate(s_chars):
            for rep in self.whitelist_alphabet:
                if rep == orig:
                    continue
                mutated = s[:i] + rep + s[i + 1 :]
                if mutated in self.whitelist_set:
                    candidates.add(mutated)

        return sorted(candidates, key=lambda bc: self.whitelist_index.get(bc, 10**12))

    @staticmethod
    def str_diff_pos_1(a: str, b: str) -> int:
        for i, (x, y) in enumerate(zip(a, b)):
            if x != y:
                return i
        return 0

    @staticmethod
    def arr_max(values: Iterable[float]) -> float:
        values = list(values)
        if not values:
            return 0.0
        return max(values)

    @staticmethod
    def min_value(*values: float) -> float:
        return min(values)

    @staticmethod
    def sum_value(data: Iterable[str]) -> float:
        s = 0.0
        for e in data:
            if e is None:
                continue
            e2 = str(e).strip()
            if e2 == "":
                continue
            try:
                s += float(e2)
            except Exception:
                continue
        return s

    def matrix_calc_maxdat_n(self, mat: List[List], mat_name: str, col_name: str) -> List[List]:
        idx_post = self.matrix_header_get_index(mat_name, col_name)
        vals = []
        for row in mat:
            try:
                vals.append(float(row[idx_post]))
            except Exception:
                vals.append(float("-inf"))
        max_val = self.arr_max(vals)

        out = []
        for row in mat:
            try:
                if float(row[idx_post]) == max_val:
                    out.append(row)
            except Exception:
                pass
        return out

    def matrix_remove_col(self, mat: List[List], mat_name: str, col_idx: int) -> List[List]:
        self.matrix_header_remove(mat_name, col_idx)
        out: List[List] = []
        for row in mat:
            out.append([v for i, v in enumerate(row) if i != col_idx])
        return out

    def matrix_calc_mtx_re_barcode_umi(self, mat: List[List], mat_name: str) -> None:
        idx_v1 = self.matrix_header_get_index(mat_name, "V1")
        idx_v3 = self.matrix_header_get_index(mat_name, "V3")
        idx_bu = self.matrix_header_get_index(mat_name, "BARCODE_UMI")

        for row in mat:
            while len(row) <= idx_bu:
                row.append("")
            v1 = row[idx_v1] if idx_v1 < len(row) else ""
            v3 = row[idx_v3] if idx_v3 < len(row) else ""
            row[idx_bu] = f"{v1}_{v3}"

    def matrix_calc_mtx_re_d_mismatch_wt_assign(self, mat: List[List], header_name: str) -> None:
        idx_wt = self.matrix_header_require_index(header_name, "mismatch.WT")
        idx_call = self.matrix_header_require_index(header_name, "mismatch.WT.call")

        wt_allow = self._safe_float(self.wt_alw, 0.0)
        out = []
        for row in mat:
            mm_wt = self._row_float(row, idx_wt)
            if mm_wt is not None and mm_wt <= wt_allow:
                out.append(row)
        if out:
            for row in mat:
                mm_wt = self._row_float(row, idx_wt)
                if mm_wt is not None and mm_wt <= wt_allow:
                    row[idx_call] = "WT"

    def matrix_calc_mtx_merge_u_wt_assign(self, mat: List[List], header_name: str) -> None:
        idx_wt = self.matrix_header_require_index(header_name, "mismatch.WT")
        idx_call = self.matrix_header_require_index(header_name, "mismatch.WT.call")

        wt_allow = self._safe_float(self.wt_alw, 0.0)
        out = []
        for row in mat:
            mm_wt = self._row_float(row, idx_wt)
            if mm_wt is not None and mm_wt <= wt_allow:
                out.append(row)
        if out:
            for row in mat:
                mm_wt = self._row_float(row, idx_wt)
                if mm_wt is not None and mm_wt <= wt_allow:
                    row[idx_call] = "WT"

    def matrix_calc_mtx_merge_u_mut_assign(self, mat: List[List], header_name: str) -> None:
        idx_mut = self.matrix_header_require_index(header_name, "mismatch.MUT")
        idx_call = self.matrix_header_require_index(header_name, "mismatch.MUT.call")

        # Keep historical behavior: compares against WT allowance.
        wt_allow = self._safe_float(self.wt_alw, 0.0)
        out = []
        for row in mat:
            mm_mut = self._row_float(row, idx_mut)
            if mm_mut is not None and mm_mut <= wt_allow:
                out.append(row)
        if out:
            for row in mat:
                mm_mut = self._row_float(row, idx_mut)
                if mm_mut is not None and mm_mut <= wt_allow:
                    row[idx_call] = "MUT"

    def matrix_calc_mtx_merge_u_wt_call_assign(self, mat: List[List], header_name: str) -> None:
        idx_1 = self.matrix_header_require_index(header_name, "mismatch.WT.call")
        idx_2 = self.matrix_header_require_index(header_name, "mismatch.call")

        out = [row for row in mat if row[idx_1] == "WT"]
        if out:
            for row in mat:
                if row[idx_1] == "WT":
                    row[idx_2] = "WT"

    def matrix_calc_mtx_merge_u_mut_call_assign(self, mat: List[List], header_name: str) -> None:
        idx_1 = self.matrix_header_require_index(header_name, "mismatch.MUT.call")
        idx_2 = self.matrix_header_require_index(header_name, "mismatch.call")

        out = [row for row in mat if row[idx_1] == "MUT"]
        if out:
            for row in mat:
                if row[idx_1] == "MUT":
                    row[idx_2] = "MUT"

    def matrix_calc_mtx_merge_u_amb_assign(self, mat: List[List], header_name: str) -> None:
        idx_1 = self.matrix_header_require_index(header_name, "mismatch.WT.call")
        idx_2 = self.matrix_header_require_index(header_name, "mismatch.MUT.call")
        idx_3 = self.matrix_header_require_index(header_name, "mismatch.call")

        for row in mat:
            if row[idx_1] == "" and row[idx_2] == "":
                row[idx_3] = "amb"

    def matrix_calc_mtx_re_d_mismatch_mut_assign(self, mat: List[List], header_name: str) -> None:
        idx_mut = self.matrix_header_require_index(header_name, "mismatch.MUT")
        idx_mut_call = self.matrix_header_require_index(header_name, "mismatch.MUT.call")

        wt_allow = self._safe_float(self.wt_alw, 0.0)
        for row in mat:
            mm_mut = self._row_float(row, idx_mut)
            if mm_mut is not None and mm_mut <= wt_allow:
                row[idx_mut_call] = "MUT"

    def matrix_calc_mtx_re_d_mismatch_wt_call_assign(self, mat: List[List], header_name: str) -> None:
        idx_wt_call = self.matrix_header_require_index(header_name, "mismatch.WT.call")
        idx_mm_call = self.matrix_header_require_index(header_name, "mismatch.call")

        for row in mat:
            if row[idx_wt_call] == "WT":
                row[idx_mm_call] = "WT"

    def matrix_calc_mtx_re_d_mismatch_mut_call_assign(self, mat: List[List], header_name: str) -> None:
        idx_mut_call = self.matrix_header_require_index(header_name, "mismatch.MUT.call")
        idx_mm_call = self.matrix_header_require_index(header_name, "mismatch.call")

        for row in mat:
            if row[idx_mut_call] == "MUT":
                row[idx_mm_call] = "MUT"

    def matrix_calc_mtx_re_d_mismatch_amb_assign(self, mat: List[List], header_name: str) -> None:
        idx_wt_call = self.matrix_header_require_index(header_name, "mismatch.WT.call")
        idx_mut_call = self.matrix_header_require_index(header_name, "mismatch.MUT.call")
        idx_mm_call = self.matrix_header_require_index(header_name, "mismatch.call")

        for row in mat:
            if row[idx_mut_call] == "MUT" and row[idx_wt_call] == "WT":
                row[idx_mm_call] = "amb"

    def matrix_calc_get_mismatch_call(self, mat: List[List], header_name: str) -> Tuple[int, int, int]:
        idx = self.matrix_header_get_index(header_name, "mismatch.call")
        cnt_wt = 0
        cnt_mut = 0
        cnt_amb = 0

        for row in mat:
            if row[idx] == "WT":
                cnt_wt += 1
            if row[idx] == "MUT":
                cnt_mut += 1
            if row[idx] == "amb":
                cnt_amb += 1

        return cnt_wt, cnt_mut, cnt_amb

    @staticmethod
    def sed_repeat(
        file_in: str,
        file_out: Optional[str],
        line_no_select: int,
        lines_repeat: int,
        cut_start: Optional[int] = None,
        cut_end: Optional[int] = None,
    ) -> Optional[List[str]]:
        if cut_start is not None and cut_end is not None:
            cut_start0 = cut_start - 1
            cut_end0 = cut_end - 1
        else:
            cut_start0 = None
            cut_end0 = None

        out: List[str] = []
        line_cnt = 1

        with Path(file_in).open("r", encoding="utf-8") as fin:
            for line in fin:
                line = line.rstrip("\n\r")

                if line_cnt == line_no_select:
                    if cut_start0 is not None and cut_end0 is not None:
                        selected = line[cut_start0 : cut_end0 + 1]
                        out.append(selected)
                    else:
                        out.append(line)

                mod = line_cnt % lines_repeat
                if mod == 0:
                    line_cnt = 1
                else:
                    line_cnt += 1

        if file_out is not None:
            with Path(file_out).open("w", encoding="utf-8") as fout:
                fout.write("\n".join(out))
            return None
        return out

    def paste_cols(self, *args: str) -> None:
        files = list(args[:-1])
        out_file = args[-1]

        handles = [Path(fp).open("r", encoding="utf-8") for fp in files]
        try:
            with Path(out_file).open("w", encoding="utf-8") as fout:
                while True:
                    lines = [h.readline() for h in handles]
                    if any(line == "" for line in lines):
                        break
                    vals = [line.rstrip("\n\r") for line in lines]
                    fout.write("\t".join(vals) + "\n")
        finally:
            for h in handles:
                h.close()

    @staticmethod
    def load_phred_scores() -> Dict[str, Dict[str, float | int]]:
        table = """Q-Score\tP_error\tASCII_Code\tSymbol
0\t1.00000\t33\t!
1\t0.79433\t34\t\"
2\t0.63096\t35\t#
3\t0.50119\t36\t$
4\t0.39811\t37\t%
5\t0.31623\t38\t&
6\t0.25119\t39\t'
7\t0.19953\t40\t(
8\t0.15849\t41\t)
9\t0.12589\t42\t*
10\t0.10000\t43\t+
11\t0.07943\t44\t,
12\t0.06310\t45\t-
13\t0.05012\t46\t.
14\t0.03981\t47\t/
15\t0.03162\t48\t0
16\t0.02512\t49\t1
17\t0.01995\t50\t2
18\t0.01585\t51\t3
19\t0.01259\t52\t4
20\t0.01000\t53\t5
21\t0.00794\t54\t6
22\t0.00631\t55\t7
23\t0.00501\t56\t8
24\t0.00398\t57\t9
25\t0.00316\t58\t:
26\t0.00251\t59\t;
27\t0.00200\t60\t<
28\t0.00158\t61\t=
29\t0.00126\t62\t>
30\t0.00100\t63\t?
31\t0.00079\t64\t@
32\t0.00063\t65\tA
33\t0.00050\t66\tB
34\t0.00040\t67\tC
35\t0.00032\t68\tD
36\t0.00025\t69\tE
37\t0.00020\t70\tF
38\t0.00016\t71\tG
39\t0.00013\t72\tH
40\t0.00010\t73\tI
41\t0.00008\t74\tJ"""

        out: Dict[str, Dict[str, float | int]] = {}
        for line in table.strip().splitlines()[1:]:
            q_score, p_err, ascii_code, symbol = line.split("\t")
            out[symbol] = {"q_score": int(q_score), "p_err": float(p_err)}
        return out

    def un_gzip(self, f1: str, f2: str) -> Tuple[str, str]:
        def maybe_gunzip(path: str) -> str:
            if not path.lower().endswith(".gz"):
                return path
            out = path[:-3]
            try:
                with gzip.open(path, "rb") as fin, open(out, "wb") as fout:
                    shutil.copyfileobj(fin, fout)
            except Exception as exc:
                self.log(f"Error decompressing '{path}': {exc}")
                sys.exit(1)
            return out

        return maybe_gunzip(f1), maybe_gunzip(f2)

    def stringdist_hamm_pre_calc(self, mtx_sub: List[List]) -> None:
        self.barcode_hamm_cache = {}
        n = len(mtx_sub)
        if n == 0:
            return

        for r, row in enumerate(mtx_sub):
            ori_bc = row[0]
            if ori_bc not in self.barcode_hamm_cache:
                self.barcode_hamm_cache[ori_bc] = self.stringdist_hamm(ori_bc, 1)
            if self.verbose and (r + 1) % 50000 == 0:
                self.log(f"hamming-precalc progress: {r + 1}/{n}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = IronThroneGoT.parser()
    args = parser.parse_args(argv)

    runner = IronThroneGoT(args)
    try:
        runner.run()
    finally:
        runner.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
