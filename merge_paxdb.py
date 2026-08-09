#!/usr/bin/env python
"""Merge all per-tissue PaxDb integrated protein abundance datasets for a
species into a single wide table (one 'id' column + one column per tissue),
mirroring the layout of the GTEx gene_median_tpm file used elsewhere in the
diseasemodulediscovery pipeline (bin/tissue_specific_filtering.py).
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
import requests

PAXDB_VERSION = "6.1"
LISTING_URL = "https://pax-db.org/downloads/{version}/datasets/{organism}/"
DATASET_URL = "https://pax-db.org/downloads/{version}/datasets/{organism}/{filename}"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Merge all PaxDb tissue datasets for one organism into one wide TSV"
    )
    parser.add_argument(
        "--organism",
        default="9606",
        help="NCBI taxonomy id of the organism (default: 9606, human)",
    )
    parser.add_argument(
        "--version", default=PAXDB_VERSION, help=f"PaxDb dataset version (default: {PAXDB_VERSION})"
    )
    parser.add_argument(
        "--output",
        default="paxdb_all_tissues.tsv",
        help="path to write the merged wide TSV to",
    )
    parser.add_argument(
        "--raw-dir",
        default="paxdb_raw",
        help="directory to cache the raw per-tissue files in, so re-runs skip re-downloading",
    )
    return parser.parse_args(argv)


def list_tissue_files(organism, version):
    url = LISTING_URL.format(version=version, organism=organism)
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    filenames = re.findall(r'href="([^"]+-integrated\.txt)"', response.text)
    if not filenames:
        raise RuntimeError(f"No integrated tissue datasets found at {url}")
    return sorted(set(filenames))


def tissue_name_from_filename(organism, filename):
    match = re.match(rf"^{re.escape(organism)}-(.+)-integrated\.txt$", filename)
    if not match:
        raise ValueError(f"Unexpected PaxDb filename: {filename}")
    return match.group(1)


def download_tissue_file(organism, version, filename, raw_dir):
    raw_dir.mkdir(parents=True, exist_ok=True)
    local_path = raw_dir / filename
    if not local_path.exists():
        url = DATASET_URL.format(version=version, organism=organism, filename=filename)
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        local_path.write_text(response.text)
    return local_path


def load_tissue_expression(local_path, tissue):
    df = pd.read_csv(
        local_path,
        sep="\t",
        comment="#",
        header=None,
        names=["gene_name", "string_external_id", "abundance"],
    )
    # a gene symbol can map to several STRING/Ensembl protein ids; collapse
    # them the same way convert_id_space() does for other id spaces
    df = df.groupby("gene_name", as_index=False)["abundance"].sum()
    return df.rename(columns={"gene_name": "id", "abundance": tissue})


def main(argv=None):
    args = parse_args(argv)
    raw_dir = Path(args.raw_dir)

    filenames = list_tissue_files(args.organism, args.version)
    print(f"Found {len(filenames)} tissue datasets for organism {args.organism}", file=sys.stderr)

    merged = None
    for i, filename in enumerate(filenames, start=1):
        tissue = tissue_name_from_filename(args.organism, filename)
        print(f"[{i}/{len(filenames)}] {tissue}", file=sys.stderr)
        local_path = download_tissue_file(args.organism, args.version, filename, raw_dir)
        tissue_df = load_tissue_expression(local_path, tissue)
        merged = (
            tissue_df
            if merged is None
            else merged.merge(tissue_df, on="id", how="outer")
        )

    merged.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote merged table with {merged.shape[0]} genes x {merged.shape[1] - 1} tissues to {args.output}")


if __name__ == "__main__":
    sys.exit(main())
