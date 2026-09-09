import requests
import pandas as pd 
from functools import reduce 
from tqdm import tqdm



BASE = "https://www.proteomicsdb.org/proteomicsdb/logic/api_v2/api.xsodata"

tissue_list = pd.read_csv("data/proteomicsdb_tissue_list.tsv", sep="\t")

bto_ids = tissue_list["TISSUE_ID"].dropna().unique().tolist()

def get_tissue_expression_v2(bto_id):
    url = f"{BASE}/Tissue('{bto_id}')/ProteinExpression"

    r = requests.get(url, params={"$format": "json"})
    r.raise_for_status()

    data = r.json()
    rows = data.get("d", {}).get("results", data)

    df = pd.DataFrame(rows)
    df = df.drop(columns=["__metadata"], errors="ignore")

    df = df[["PROTEIN_ID", "EXPRESSION"]].copy()
    df["EXPRESSION"] = pd.to_numeric(
        df["EXPRESSION"],
        errors="coerce"
    )
    # make sure one protein only has one value per tissue
    df = (
        df.groupby("PROTEIN_ID", as_index=False)["EXPRESSION"]
        .mean()
    )


    df = df.rename(columns={"EXPRESSION": bto_id})

    return df

all_tissue_dfs = []

for bto_id in tqdm(bto_ids, desc="querying for tissue specific expressions"):
    print(f"Querying {bto_id}")

    try:
        df_tissue = get_tissue_expression_v2(bto_id)
        all_tissue_dfs.append(df_tissue)

    except Exception as e:
        print(f"Skipping {bto_id}: {e}")

expression_matrix = reduce(
    lambda left, right: pd.merge(left, right, on="PROTEIN_ID", how="outer"),
    all_tissue_dfs
)
expression_matrix.fillna(value=0,inplace=True)

def get_protein_mapping(protein_id):
    url = f"{BASE}/Protein({int(protein_id)})"

    r = requests.get(url, params={"$format": "json"})
    r.raise_for_status()

    data = r.json().get("d", r.json())

    return {
        "PROTEIN_ID": int(protein_id),

        # Try common possible UniProt fields
        "UNIPROT_ID": (
            data.get("UNIQUE_IDENTIFIER")
            or data.get("UNIPROT_ID")
            or data.get("UNIPROT_ACCESSION")
            or data.get("ACCESSION")
        ),

        "DATABASE": data.get("DATABASE"),
        "ENTRY_NAME": data.get("ENTRY_NAME"),
        "GENE_NAME": data.get("GENE_NAME"),
        "PROTEIN_DESCRIPTION": data.get("PROTEIN_DESCRIPTION"),
    }

protein_ids = (
    expression_matrix["PROTEIN_ID"]
    .dropna()
    .astype(int)
    .unique()
    .tolist()
)

mapping_rows = []

for protein_id in tqdm(protein_ids):
    try:
        mapping_rows.append(get_protein_mapping(protein_id))
    except Exception as e:
        print(f"Could not map PROTEIN_ID {protein_id}: {e}")

protein_mapping = pd.DataFrame(mapping_rows)


expression_matrix_uniprot = protein_mapping.merge(
    expression_matrix,
    on="PROTEIN_ID",
    how="right"
)


# Optional: if multiple ProteomicsDB IDs map to the same UniProt ID,
# average their expression values
metadata_cols = [
    "UNIPROT_ID",
    "DATABASE",
    "ENTRY_NAME",
    "GENE_NAME",
    "PROTEIN_DESCRIPTION",
    "PROTEIN_ID"
]

expression_cols = [
    col for col in expression_matrix_uniprot.columns
    if col not in metadata_cols
]

expression_matrix_uniprot_final = (
    expression_matrix_uniprot
    .groupby("UNIPROT_ID", as_index=False)[expression_cols]
    .mean()
)

expression_matrix_uniprot_final.to_csv(
    "data/proteomicsdb_expression_matrix_uniprot.tsv",
    sep="\t",
    index=False
)