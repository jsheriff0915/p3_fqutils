#!/usr/bin/env python3
"""
Host removal pipeline.

Required arguments:
  -i / --src-dir            source directory containing fastq files
  -s / --species-common-name   host common name for the reference genome
  -a / --aligner            aligner selection (bowtie2 for short reads; minimap2 for long reads)
  -o / --out-dir            output directory for filtered fastq files

This script will:
  1) search GCF accessions in a pre-computed index database, based on the aligner selection
  2) remotely download reference genome and construct a custom index (if needed)
  3) filter host reads from input fastq files and output filtered fastq file(s)
"""

import argparse
import re
import subprocess
import sys
import zipfile
import csv
from pathlib import Path

BOWTIE2_DB = Path("/home/ac.jsheriff/dev_container/projects/host_removal/reference_indices/bowtie2")
MINIMAP2_DB = Path("/home/ac.jsheriff/dev_container/projects/host_removal/reference_indices/minimap2")
DATASETS_BIN = Path("/home/ac.jsheriff/dev_container/projects/red_fox/datasets")
HOST_GCF_DICT = Path("/home/ac.jsheriff/dev_container/projects/host_removal/files/reference_index_dict_31Aug2026.csv")

BOWTIE2_THREADS = 2

# Matches sample_1.fastq, sample_R1.fastq, and .fastq.gz variants of both
FWD_READ_PATTERN = re.compile(r"_R?1\.fastq(\.gz)?$")
FASTQ_PATTERN = re.compile(r"\.fastq(\.gz)?$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the host removal pipeline.")
    parser.add_argument(
        "-i", "--src-dir", required=True, dest="src_dir",
        help="Source directory containing fastq files",
    )
    parser.add_argument(
        "-s", "--species-common-name", required=True, dest="species_common_name",
        help="Host common name for the reference genome",
    )
    parser.add_argument(
        "-a", "--aligner", required=True, choices=["bowtie2", "minimap2"],
        help="Aligner selection (bowtie2 for short reads; minimap2 for long reads)",
    )
    parser.add_argument(
        "-o", "--out-dir", required=True, dest="out_dir",
        help="Output directory for filtered fastq files",
    )
    return parser.parse_args()


def run(cmd: list, **kwargs) -> None:
    """Run a command, echoing it first, and raise on non-zero exit (mirrors `set -e`)."""
    print("+ " + " ".join(str(part) for part in cmd))
    subprocess.run(cmd, check=True, **kwargs)


def find_fastq_files(src_dir: Path):
    """
    Classifies each .fastq file in src_dir as forward, reverse, or single-end,
    and derives the sample base name (suffix _1/_2 stripped).
    Mirrors the original bash for-loop, including the fact that the last
    matching file in each category "wins" if there are multiple samples.
    """
    fwd_reads = rvs_reads = single_reads = base_name = None

    for file in sorted(src_dir.glob("*.fastq")):
        name = file.name
        stem = re.sub(r"_1$", "", file.stem)
        stem = re.sub(r"_2$", "", stem)
        base_name = stem

        if FWD_READ_PATTERN.search(name):
            fwd_reads = file
        else:
            rvs_reads = file

        if FASTQ_PATTERN.search(name):
            single_reads = file

    return fwd_reads, rvs_reads, single_reads, base_name

def find_precomputed_index(db_dir: Path, species_common_name: str, suffix: str = ""):
    """Looks for a subdirectory under db_dir starting with species common name and
    returns its index prefix path (dir/dir_name + optional suffix), or None."""

    host_gcf_dict = {}

    with open(HOST_GCF_DICT, mode='r', encoding='utf-8') as file:
        csv_reader = csv.DictReader(file)
        
        for row in csv_reader:
            # Use the 'Common_Name' column as the main key
            key = row.pop('Common_Name') 
            # The remaining row data becomes the value
            host_gcf_dict[key] = row 

    gcf_accession = None
    for key, row in host_gcf_dict.items():
        if species_common_name.strip().lower() == key.strip().lower():
            gcf_accession = row['GCF_Accession']  # confirm this matches your actual CSV header exactly
            print(f"'{species_common_name}' found within reference species list. Paired GCF Accession is '{gcf_accession}'. Using this reference index to filer host sequences...")
            break

    if gcf_accession is None:
        print(f"'{species_common_name}' not found in reference species list.")
        return None

    for match in sorted(db_dir.glob(f"{gcf_accession}*")):
        if match.is_dir():
            return match / f"{match.name}{suffix}"
    return None



def run_bowtie2(species_common_name: str, out_dir: Path, fwd_reads: Path, rvs_reads: Path, base_name: str) -> None:
    print("bowtie2 was selected as the aligner; searching for the host species within precomputed index database...")

    precomputed_index = find_precomputed_index(BOWTIE2_DB, species_common_name)

    if precomputed_index is not None:

        run([
            "hostile", "clean", "--force",
            "--fastq1", str(fwd_reads),
            "--fastq2", str(rvs_reads),
            "--aligner", "bowtie2",
            "--index", str(precomputed_index),
            "--out-dir", str(out_dir),
        ])
        return


def run_minimap2(species_common_name: str, out_dir: Path, single_reads: Path, base_name: str) -> None:
    print("minimap2 was selected as the aligner; searching for the GCF Accession within precomputed index database...")

    precomputed_index = find_precomputed_index(MINIMAP2_DB, species_common_name, suffix=".mmi")

    if precomputed_index is not None:

        run([
            "hostile", "clean", "--force",
            "--fastq1", str(single_reads),
            "--aligner", "minimap2",
            "--index", str(precomputed_index),
            "--out-dir", str(out_dir),
        ])
        return


def run_hostile(src_dir: str, species_common_name: str, aligner: str, out_dir: str) -> None:
    """
    Runs the full host-removal pipeline:
      1) classify fastq files in src_dir as forward/reverse/single-end
      2) search for a precomputed index matching gcf_accession, based on aligner
      3) if none exists, download the reference genome and build a custom index
      4) run `hostile clean` against the resolved index

    Can be called directly (e.g. imported into another script) or via the
    command-line entry point below.
    """
    src_dir = Path(src_dir)
    out_dir = Path(out_dir)

    fwd_reads, rvs_reads, single_reads, base_name = find_fastq_files(src_dir)

    if aligner == "bowtie2":
        run_bowtie2(species_common_name, out_dir, fwd_reads, rvs_reads, base_name)
    elif aligner == "minimap2":
        run_minimap2(species_common_name, out_dir, single_reads, base_name)
    else:
        print(f"Error: unrecognized aligner '{aligner}'", file=sys.stderr)
        sys.exit(1)

    print("Host removal: complete")


def main() -> None:
    args = parse_args()
    run_hostile(
        src_dir=args.src_dir,
        species_common_name=args.species_common_name,
        aligner=args.aligner,
        out_dir=args.out_dir,
    )


if __name__ == "__main__":
    main()
