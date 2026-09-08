#!/usr/bin/env python
"""Group the GTEx and PaxDB tissues used elsewhere in this project into
coarse organ-level categories, using the UBERON anatomy ontology's is_a /
part_of hierarchy (e.g. all GTEx "Brain - ..." subregions and all PaxDB
brain-related organs collapse into a single "brain" group). Each organ group
is then further collapsed to an organ-system-level group the same way (e.g.
"brain" and "spinal cord" both roll up to "nervous system"), giving two
nested tiers of grouping.

Each tissue is mapped to a UBERON term (GTEx via its own SMUBRID sample
annotation, PaxDB via the PaxDb organs API), then walked up the ontology to
the nearest ancestor belonging to a curated set of organ-level UBERON terms
(the union of UBERON's own `organ_slim` and `major_organ` subsets, plus a
handful of terms neither subset covers, see SUPPLEMENTARY_GROUPS below).
Each resulting organ group is walked up a second time to the nearest
ancestor in SYSTEM_GROUPS (UBERON's direct is_a children of "anatomical
system", i.e. the standard organ systems: nervous, digestive, ...).

Because many of UBERON's is_a/part_of edges for anatomical subregions (e.g.
brain regions) are only inferable via OWL reasoning and are not asserted in
the plain .obo file, ancestor lookups are delegated to EBI's OLS4 API (which
serves the already-reasoned closure) rather than walking the local .obo
graph directly. Results are cached under data/ontology/ so repeat runs don't
re-hit the network.

Usage:
    python group_tissues_uberon.py
    python group_tissues_uberon.py --force-refresh   # ignore all caches
"""

import argparse
import json
import re
import sys
import time
import urllib.parse
from pathlib import Path

import obonet
import pandas as pd
import requests

DATA_DIR = Path("data")
ONTOLOGY_DIR = DATA_DIR / "ontology"

UBERON_BASIC_URL = "http://purl.obolibrary.org/obo/uberon/basic.obo"
GTEX_SAMPLE_ATTRIBUTES_URL = (
    "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/"
    "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
)
PAXDB_ORGANS_URL = "https://api.pax-db.org/v6/metadata/organs"
OLS4_ANCESTORS_URL = (
    "https://www.ebi.ac.uk/ols4/api/ontologies/uberon/terms/{iri}/hierarchicalAncestors"
)

GCT_PATH = DATA_DIR / "expression_by_tissue.gct"
PAXDB_TABLE_PATH = DATA_DIR / "paxdb_all_tissues.tsv"

GTEX_MAPPING_CACHE = DATA_DIR / "gtex_tissue_uberon_mapping.tsv"
PAXDB_MAPPING_CACHE = DATA_DIR / "paxdb_tissue_uberon_mapping.tsv"
UBERON_OBO_CACHE = ONTOLOGY_DIR / "uberon_basic.obo"
OLS_ANCESTOR_CACHE = ONTOLOGY_DIR / "ols_hierarchical_ancestors_cache.json"

OUTPUT_PATH = DATA_DIR / "tissue_uberon_groups.json"
OUTPUT_SYSTEMS_PATH = DATA_DIR / "tissue_uberon_systems.json"
OUTPUT_LONG_PATH = DATA_DIR / "tissue_uberon_groups.tsv"

# UBERON terms that matter for GTEx/PaxDB tissues but fall outside both
# `organ_slim` ("organs, excluding individual muscles and skeletal elements")
# and `major_organ` (a fuzzy 23-term subset). Verified against the terms
# GTEx/PaxDB actually attach to samples (see fetch_gtex_mapping /
# fetch_paxdb_mapping), not guessed from memory.
SUPPLEMENTARY_GROUPS = {
    "UBERON:0000178": "blood",
    "UBERON:0001013": "adipose tissue",
    "UBERON:0001981": "blood vessel",
    "UBERON:0000029": "lymph node",
    "UBERON:0002371": "bone marrow",
    "UBERON:0001021": "nerve",
    "UBERON:0000002": "uterine cervix",
    "UBERON:0000996": "vagina",
    "UBERON:0003889": "fallopian tube",
    "UBERON:0008367": "breast",  # breast epithelium; sole GTEx breast tissue, untagged in organ_slim/major_organ
    "UBERON:0002372": "tonsil",  # real organ, not covered by organ_slim/major_organ
    "UBERON:0001066": "intervertebral disc",  # sits ambiguously under both cartilage element and skeletal joint
}

# UBERON terms that are technically in organ_slim/major_organ but are
# abstract grouping classes rather than real organs (tagged "upper_level":
# "abstract upper-level terms not directly useful for analysis" in the
# ontology's own subsetdefs). Left in, they cause spurious ambiguity
# (e.g. visceral adipose tissue also being an ancestor-match for "viscus").
EXCLUDE_UPPER_LEVEL_TAG = "upper_level"

# A few GTEx/PaxDB entries are annotated with non-UBERON ids (cell lines,
# cultured cell types) that can't be placed in the UBERON graph at all, but
# have an unambiguous organ of origin.
NON_UBERON_OVERRIDES = {
    "EFO:0002009": "skin of body",  # Cells - Cultured fibroblasts
    "EFO:0000572": "blood",         # Cells - EBV-transformed lymphocytes
    "EFO:0002067": "bone marrow",   # Cells - Leukemia cell line (CML)
    "CL:0000182": "liver",          # Hepatocyte
    "CL:0000127": "brain",          # Astrocyte
    "CL:0000312": "skin of body",   # Keratinocyte
    "CL:0002620": "skin of body",   # Skin fibroblast
    "CL:0000233": "blood",          # Platelet
    "CL:0000025": "ovary",          # Egg cell / oocyte
    "CL:4030032": "heart",          # Valve interstitial cell
}

# A handful of tissues are true multi-parent cases in UBERON (e.g. a
# coronary artery is simultaneously part_of "blood vessel" and part_of
# "heart"); the automatic nearest-ancestor resolution can't break the tie
# on ontology grounds alone, so the call is made explicitly here instead of
# silently picking one. Keyed by (source, tissue_name).
MANUAL_DISAMBIGUATION = {
    ("GTEx", "Artery_Coronary"): "blood vessel",
    ("GTEx", "Adipose_Visceral_Omentum"): "adipose tissue",
    ("GTEx", "Esophagus_Gastroesophageal_Junction"): "esophagus",
    ("GTEx", "Esophagus_Muscularis"): "esophagus",
    ("GTEx", "Minor_Salivary_Gland"): "saliva-secreting gland",
    # GTEx's own SMUBRID for plain "Small Intestine - Terminal Ileum" points
    # at UBERON:0001211 (Peyer's patch, a lymphoid structure) rather than at
    # the ileum itself, so it doesn't ontologically resolve to an intestine
    # group; its sibling "...Lymphoid_Aggregate" sample uses a different,
    # correctly-placed SMUBRID and needs no override.
    ("GTEx", "Small_Intestine_Terminal_Ileum"): "small intestine",
    ("GTEx", "Small_Intestine_Terminal_Ileum_Mixed_Cell"): "small intestine",
    # UBERON models "portal tract" purely as a vascular structure (its
    # ancestors include "blood vessel" but not "liver"), even though GTEx's
    # own SMTS places this sample under Liver.
    ("GTEx", "Liver_Portal_Tract"): "liver",
}

# Organ-system-level target set: UBERON's own direct is_a children of
# "anatomical system" (UBERON:0000467), restricted to the systems actually
# relevant to the organ groups above. Kidney/bladder resolve to "renal
# system", not "genitourinary system" (there's no standalone "urinary
# system" term); eye resolves to "sensory system".
SYSTEM_GROUPS = {
    "UBERON:0001016": "nervous system",
    "UBERON:0001007": "digestive system",
    "UBERON:0001004": "respiratory system",
    "UBERON:0004535": "cardiovascular system",
    "UBERON:0002204": "musculoskeletal system",
    "UBERON:0002416": "integumental system",
    "UBERON:0000990": "reproductive system",
    "UBERON:0000949": "endocrine system",
    "UBERON:0002405": "immune system",
    "UBERON:0001008": "renal system",
    "UBERON:0001032": "sensory system",
}

# Same rationale as MANUAL_DISAMBIGUATION, one tier up: organ groups that are
# either true multi-parent cases among SYSTEM_GROUPS (liver, pituitary
# gland, tonsil - genuinely dual/triple function organs) or aren't linked to
# any SYSTEM_GROUPS ancestor by an asserted/reasoned UBERON edge at all
# (muscle organ, pancreas, breast), despite an obvious system of origin.
# Keyed by organ group label. Adipose tissue is deliberately left
# unresolved here: it doesn't cleanly belong to one classical organ system
# (subcutaneous vs. visceral depots differ), so it's reported as ungrouped
# rather than forced into one.
SYSTEM_MANUAL_DISAMBIGUATION = {
    "liver": "digestive system",
    "pancreas": "digestive system",
    "pituitary gland": "endocrine system",
    "tonsil": "immune system",
    "breast": "integumental system",
    "muscle organ": "musculoskeletal system",
    # UBERON classifies blood under "hemolymphoid system", which isn't one
    # of the SYSTEM_GROUPS above; grouped with the system it physically
    # flows through instead.
    "blood": "cardiovascular system",
}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Group GTEx and PaxDB tissues into UBERON organ-level categories"
    )
    parser.add_argument(
        "--gct", default=str(GCT_PATH), help="GTEx median-TPM .gct file (for tissue column names)"
    )
    parser.add_argument(
        "--paxdb-table",
        default=str(PAXDB_TABLE_PATH),
        help="merged PaxDB wide TSV (for tissue column names)",
    )
    parser.add_argument("--output", default=str(OUTPUT_PATH), help="path to write the grouping JSON to")
    parser.add_argument(
        "--force-refresh",
        action="store_true",
        help="ignore all local caches and re-download/re-query everything",
    )
    return parser.parse_args(argv)


# --------------------------------------------------------------------------
# Downloads / caches
# --------------------------------------------------------------------------

def download_text(url, dest: Path, force=False):
    if dest.exists() and not force:
        return dest.read_text()
    dest.parent.mkdir(parents=True, exist_ok=True)
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    dest.write_text(response.text)
    return response.text


def fetch_gtex_mapping(force=False) -> pd.DataFrame:
    """SMTSD (tissue detail, matches .gct columns after normalization) ->
    SMUBRID (Uberon id), sourced from GTEx's own sample attributes file."""
    if GTEX_MAPPING_CACHE.exists() and not force:
        return pd.read_csv(GTEX_MAPPING_CACHE, sep="\t", dtype=str)

    print("Downloading GTEx sample attributes ...", file=sys.stderr)
    response = requests.get(GTEX_SAMPLE_ATTRIBUTES_URL, timeout=120)
    response.raise_for_status()
    df = pd.read_csv(
        pd.io.common.StringIO(response.text), sep="\t", usecols=["SMTS", "SMTSD", "SMUBRID"], dtype=str
    )
    mapping = df.drop_duplicates(subset=["SMTSD"]).sort_values("SMTSD").reset_index(drop=True)
    mapping.to_csv(GTEX_MAPPING_CACHE, sep="\t", index=False)
    return mapping


def fetch_paxdb_mapping(force=False) -> pd.DataFrame:
    """PaxDb organ_name -> organ_id (Uberon/CL/BTO), via the PaxDb API."""
    if PAXDB_MAPPING_CACHE.exists() and not force:
        return pd.read_csv(PAXDB_MAPPING_CACHE, sep="\t", dtype=str)

    print("Fetching PaxDB organ list ...", file=sys.stderr)
    response = requests.get(PAXDB_ORGANS_URL, params={"species": 9606}, timeout=30)
    response.raise_for_status()
    df = pd.DataFrame(response.json()["data"])
    df = df.rename(columns={"organ_id": "id", "organ_name": "name"})
    df["id"] = df["id"].str.replace("_", ":", n=1, regex=False)
    df.to_csv(PAXDB_MAPPING_CACHE, sep="\t", index=False)
    return df


def fetch_uberon_graph(force=False):
    download_text(UBERON_BASIC_URL, UBERON_OBO_CACHE, force=force)
    return obonet.read_obo(UBERON_OBO_CACHE)


# --------------------------------------------------------------------------
# Tissue name loading (only tissues actually used in this project's data)
# --------------------------------------------------------------------------

def load_gtex_tissue_columns(gct_path) -> list:
    with open(gct_path) as f:
        f.readline()  # "#1.2"
        f.readline()  # dims
        header = f.readline()
    return header.strip().split("\t")[2:]


def load_paxdb_tissue_columns(table_path) -> list:
    header = pd.read_csv(table_path, sep="\t", nrows=0).columns.tolist()
    return [c for c in header if c != "id"]


def normalize(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", name.lower())


# --------------------------------------------------------------------------
# UBERON ancestor resolution via OLS4 (handles reasoned/inferred edges that
# the plain .obo file is missing)
# --------------------------------------------------------------------------

class OntologyResolver:
    """Resolves a UBERON term to the nearest ancestor in a given target-group
    set. target_groups is passed per-call (not fixed at construction) so the
    same ancestor cache can be shared across different grouping tiers
    (organ-level, system-level, ...) without duplicate network calls."""

    def __init__(self, cache_path: Path, force_refresh=False):
        self.cache_path = cache_path
        self.cache = {} if force_refresh else self._load_cache()

    def _load_cache(self):
        if self.cache_path.exists():
            return json.loads(self.cache_path.read_text())
        return {}

    def _save_cache(self):
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        self.cache_path.write_text(json.dumps(self.cache, indent=0))

    def hierarchical_ancestors(self, term_id: str) -> set:
        if term_id in self.cache:
            return set(self.cache[term_id])
        iri = f"http://purl.obolibrary.org/obo/{term_id.replace(':', '_')}"
        encoded = urllib.parse.quote(urllib.parse.quote(iri, safe=""), safe="")
        url = OLS4_ANCESTORS_URL.format(iri=encoded)
        ancestors = set()
        page = 0
        while True:
            response = requests.get(url, params={"size": 200, "page": page}, timeout=30)
            if response.status_code == 404:
                break
            response.raise_for_status()
            body = response.json()
            ancestors.update(t["obo_id"] for t in body.get("_embedded", {}).get("terms", []))
            total_pages = body.get("page", {}).get("totalPages", 1)
            page += 1
            if page >= total_pages:
                break
        self.cache[term_id] = sorted(ancestors)
        time.sleep(0.05)
        return ancestors

    def resolve_group(self, term_id: str, target_groups: dict):
        """Return (group_id, group_label) for the nearest target_groups
        ancestor of term_id, or None if no target group is an ancestor."""
        if term_id in target_groups:
            return term_id, target_groups[term_id]

        ancestors = self.hierarchical_ancestors(term_id)
        candidates = ancestors & target_groups.keys()
        if not candidates:
            return None
        if len(candidates) == 1:
            cid = next(iter(candidates))
            return cid, target_groups[cid]

        # multiple target-group ancestors matched: keep only the most
        # specific ones (those that are not themselves an ancestor of
        # another candidate)
        most_specific = set(candidates)
        for cid in candidates:
            other_ancestors = self.hierarchical_ancestors(cid)
            most_specific -= (other_ancestors & candidates) - {cid}
        if len(most_specific) == 1:
            cid = next(iter(most_specific))
            return cid, target_groups[cid]
        # genuinely ambiguous (e.g. shared between two systems): report all
        # of them rather than silently picking one
        labels = tuple(sorted(target_groups[c] for c in most_specific))
        return tuple(sorted(most_specific)), labels

    def save(self):
        self._save_cache()


def build_target_groups(graph) -> dict:
    groups = {}
    for node, data in graph.nodes(data=True):
        subsets = data.get("subset", [])
        if EXCLUDE_UPPER_LEVEL_TAG in subsets:
            continue
        if "organ_slim" in subsets or "major_organ" in subsets:
            groups[node] = data.get("name", node)
    groups.update(SUPPLEMENTARY_GROUPS)
    return groups


# --------------------------------------------------------------------------
# Main grouping logic
# --------------------------------------------------------------------------

def _resolve_one(source, tissue, uberon_id, resolver, target_groups, assignments, ungrouped, ambiguous):
    override = MANUAL_DISAMBIGUATION.get((source, tissue))
    if override:
        assignments[tissue] = override
        return
    if uberon_id is None or (isinstance(uberon_id, float) and pd.isna(uberon_id)):
        ungrouped[tissue] = uberon_id
        return
    if not str(uberon_id).startswith("UBERON"):
        group = NON_UBERON_OVERRIDES.get(uberon_id)
        if group:
            assignments[tissue] = group
        else:
            ungrouped[tissue] = uberon_id
        return
    result = resolver.resolve_group(uberon_id, target_groups)
    if result is None:
        ungrouped[tissue] = uberon_id
    elif isinstance(result[0], tuple):
        ambiguous[tissue] = result[1]
    else:
        assignments[tissue] = result[1]


def group_gtex(tissue_columns, gtex_mapping, resolver, target_groups):
    lookup = {normalize(row.SMTSD): row.SMUBRID for row in gtex_mapping.itertuples()}
    assignments, ungrouped, ambiguous = {}, {}, {}
    for tissue in tissue_columns:
        uberon_id = lookup.get(normalize(tissue.replace("_", " ")))
        _resolve_one("GTEx", tissue, uberon_id, resolver, target_groups, assignments, ungrouped, ambiguous)
    return assignments, ungrouped, ambiguous


def group_paxdb(tissue_columns, paxdb_mapping, resolver, target_groups):
    lookup = dict(zip(paxdb_mapping["name"], paxdb_mapping["id"]))
    assignments, ungrouped, ambiguous = {}, {}, {}
    for tissue in tissue_columns:
        uberon_id = lookup.get(tissue)
        _resolve_one("PaxDB", tissue, uberon_id, resolver, target_groups, assignments, ungrouped, ambiguous)
    return assignments, ungrouped, ambiguous


def group_organs_into_systems(organ_assignments, organ_name_to_id, resolver):
    """One tier up: map each tissue's already-resolved organ group to the
    nearest organ-system ancestor (e.g. "brain" -> "nervous system")."""
    assignments, ungrouped, ambiguous = {}, {}, {}
    for tissue, organ_group in organ_assignments.items():
        override = SYSTEM_MANUAL_DISAMBIGUATION.get(organ_group)
        if override:
            assignments[tissue] = override
            continue
        organ_id = organ_name_to_id[organ_group]
        result = resolver.resolve_group(organ_id, SYSTEM_GROUPS)
        if result is None:
            ungrouped[tissue] = organ_group
        elif isinstance(result[0], tuple):
            ambiguous[tissue] = result[1]
        else:
            assignments[tissue] = result[1]
    return assignments, ungrouped, ambiguous


def invert_to_group_view(gtex_assignments, paxdb_assignments):
    groups = {}
    for tissue, group in gtex_assignments.items():
        groups.setdefault(group, {"GTEx": [], "PaxDB": []})["GTEx"].append(tissue)
    for tissue, group in paxdb_assignments.items():
        groups.setdefault(group, {"GTEx": [], "PaxDB": []})["PaxDB"].append(tissue)
    return dict(sorted(groups.items()))


def main(argv=None):
    args = parse_args(argv)

    graph = fetch_uberon_graph(force=args.force_refresh)
    target_groups = build_target_groups(graph)
    organ_name_to_id = {name: uid for uid, name in target_groups.items()}
    print(f"Target organ groups: {len(target_groups)}", file=sys.stderr)
    print(f"Target system groups: {len(SYSTEM_GROUPS)}", file=sys.stderr)

    gtex_mapping = fetch_gtex_mapping(force=args.force_refresh)
    paxdb_mapping = fetch_paxdb_mapping(force=args.force_refresh)

    gtex_tissues = load_gtex_tissue_columns(args.gct)
    paxdb_tissues = load_paxdb_tissue_columns(args.paxdb_table)

    resolver = OntologyResolver(OLS_ANCESTOR_CACHE, force_refresh=args.force_refresh)

    gtex_assignments, gtex_ungrouped, gtex_ambiguous = group_gtex(
        gtex_tissues, gtex_mapping, resolver, target_groups
    )
    paxdb_assignments, paxdb_ungrouped, paxdb_ambiguous = group_paxdb(
        paxdb_tissues, paxdb_mapping, resolver, target_groups
    )

    gtex_systems, gtex_sys_ungrouped, gtex_sys_ambiguous = group_organs_into_systems(
        gtex_assignments, organ_name_to_id, resolver
    )
    paxdb_systems, paxdb_sys_ungrouped, paxdb_sys_ambiguous = group_organs_into_systems(
        paxdb_assignments, organ_name_to_id, resolver
    )

    resolver.save()

    group_view = invert_to_group_view(gtex_assignments, paxdb_assignments)
    print(json.dumps(group_view, indent=2))

    Path(args.output).write_text(json.dumps(group_view, indent=2))
    print(f"\nWrote {len(group_view)} organ groups to {args.output}", file=sys.stderr)

    system_view = invert_to_group_view(gtex_systems, paxdb_systems)
    Path(OUTPUT_SYSTEMS_PATH).write_text(json.dumps(system_view, indent=2))
    print(f"Wrote {len(system_view)} system groups to {OUTPUT_SYSTEMS_PATH}", file=sys.stderr)

    long_rows = [
        {
            "source": "GTEx",
            "tissue": t,
            "organ_group": gtex_assignments[t],
            "organ_system": gtex_systems.get(t, gtex_assignments[t]),
        }
        for t in gtex_assignments
    ] + [
        {
            "source": "PaxDB",
            "tissue": t,
            "organ_group": paxdb_assignments[t],
            "organ_system": paxdb_systems.get(t, paxdb_assignments[t]),
        }
        for t in paxdb_assignments
    ]
    pd.DataFrame(long_rows).sort_values(["source", "tissue"]).to_csv(
        OUTPUT_LONG_PATH, sep="\t", index=False
    )
    print(f"Wrote {len(long_rows)} tissue->group rows to {OUTPUT_LONG_PATH}", file=sys.stderr)

    for label, ungrouped in (("GTEx", gtex_ungrouped), ("PaxDB", paxdb_ungrouped)):
        if ungrouped:
            print(f"\n{label} tissues with no organ-group match ({len(ungrouped)}):", file=sys.stderr)
            for tissue, uid in ungrouped.items():
                print(f"  {tissue!r} (uberon={uid})", file=sys.stderr)
    for label, ambiguous in (("GTEx", gtex_ambiguous), ("PaxDB", paxdb_ambiguous)):
        if ambiguous:
            print(f"\n{label} tissues matching >1 organ group ({len(ambiguous)}):", file=sys.stderr)
            for tissue, labels in ambiguous.items():
                print(f"  {tissue!r} -> {labels}", file=sys.stderr)

    for label, ungrouped in (("GTEx", gtex_sys_ungrouped), ("PaxDB", paxdb_sys_ungrouped)):
        if ungrouped:
            print(f"\n{label} organ groups with no system-group match ({len(ungrouped)}):", file=sys.stderr)
            for tissue, organ_group in ungrouped.items():
                print(f"  {tissue!r} (organ group={organ_group!r})", file=sys.stderr)
    for label, ambiguous in (("GTEx", gtex_sys_ambiguous), ("PaxDB", paxdb_sys_ambiguous)):
        if ambiguous:
            print(f"\n{label} organ groups matching >1 system group ({len(ambiguous)}):", file=sys.stderr)
            for tissue, labels in ambiguous.items():
                print(f"  {tissue!r} -> {labels}", file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main())
