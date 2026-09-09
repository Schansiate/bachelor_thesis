#!/usr/bin/env python
"""Download the per-cohort TCGA STAR-TPM expression matrices from the UCSC
Xena GDC hub, mirroring the layout of data/paxdb_raw/: one raw file per
context, cached locally so the diseasemodulediscovery pipeline's TCGA
filtering_source can read them from disk instead of a remote URL.
"""

import argparse
import sys
from pathlib import Path

import requests
from tqdm import tqdm

BASE_URL = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"
FILENAME_TEMPLATE = "TCGA-{cohort}.star_tpm.tsv.gz"

# the 33 standard TCGA project abbreviations hosted on the GDC hub
TCGA_COHORTS = [
    "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM",
    "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
    "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
    "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM",
]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Download raw per-cohort TCGA STAR-TPM matrices into a local cache dir"
    )
    parser.add_argument(
        "--cohorts",
        nargs="+",
        default=TCGA_COHORTS,
        help="TCGA project abbreviations to download (default: all 33 standard cohorts)",
    )
    parser.add_argument(
        "--raw-dir",
        default="data/tcga_raw",
        help="directory to cache the raw per-cohort files in (default: data/tcga_raw)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="re-download files even if they already exist locally",
    )
    return parser.parse_args(argv)


def download_cohort_file(cohort, raw_dir, force):
    raw_dir.mkdir(parents=True, exist_ok=True)
    filename = FILENAME_TEMPLATE.format(cohort=cohort)
    local_path = raw_dir / filename
    if local_path.exists() and not force:
        return local_path, False

    url = f"{BASE_URL}/{filename}"
    with requests.get(url, stream=True, timeout=60) as response:
        response.raise_for_status()
        total = int(response.headers.get("content-length", 0))
        tmp_path = local_path.with_suffix(local_path.suffix + ".part")
        with open(tmp_path, "wb") as f, tqdm(
            total=total, unit="B", unit_scale=True, desc=cohort, leave=False
        ) as bar:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                f.write(chunk)
                bar.update(len(chunk))
        tmp_path.rename(local_path)
    return local_path, True


def main(argv=None):
    args = parse_args(argv)
    raw_dir = Path(args.raw_dir)

    for i, cohort in enumerate(args.cohorts, start=1):
        print(f"[{i}/{len(args.cohorts)}] TCGA-{cohort}", file=sys.stderr)
        try:
            local_path, downloaded = download_cohort_file(cohort, raw_dir, args.force)
        except Exception as e:
            print(f"  Skipping TCGA-{cohort}: {e}", file=sys.stderr)
            continue
        status = "downloaded" if downloaded else "already cached"
        print(f"  {status}: {local_path}", file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main())
