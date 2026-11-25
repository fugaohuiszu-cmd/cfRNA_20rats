# 用途：细胞注释检查。
# 输入：res0.3标记基因表和细胞类型评分。
# 输出：注释检查表。

from pathlib import Path
import os
import pandas as pd

PROJECT = Path(os.environ.get("CFRNA_GMF_ROOT", ".")).resolve()
REF = (
    PROJECT
    / "CellFreeGMF_work"
    / "09_Set7_GSE169396_human_femoral_head_GMF"
    / "GSE169396_human_femoral_head_reference"
)
ASSESS = REF / "resolution_assessment"
OUT = REF / "curated_annotation"
OUT.mkdir(parents=True, exist_ok=True)

TOP50 = ASSESS / "markers" / "GSE169396_res0p3_top50_markers.csv"
CALLS = ASSESS / "GSE169396_resolution_cluster_marker_type_calls.csv"
if not TOP50.exists():
    raise FileNotFoundError(TOP50)
if not CALLS.exists():
    raise FileNotFoundError(CALLS)

markers = pd.read_csv(TOP50)
calls = pd.read_csv(CALLS)
calls = calls[calls["resolution"].round(3).eq(0.3)].copy()
calls["cluster"] = calls["cluster"].astype(str)
markers["cluster"] = markers["cluster"].astype(str)

canonical_sets = {
    "Granulocytes": {
        "positive": ["S100A8", "S100A9", "FCGR3B", "CSF3R", "CXCR2", "CEACAM8", "MPO", "ELANE", "PRTN3", "AZU1", "CTSG", "DEFA4"],
        "negative": ["CD3D", "CD79A", "PECAM1", "COL1A1", "HBA1"],
    },
    "T_cells": {
        "positive": ["CD3D", "CD3E", "CD3G", "TRAC", "CD2", "IL7R", "CD8A", "CD8B", "GNLY", "KLRD1", "GZMH", "NKG7"],
        "negative": ["MS4A1", "CD79A", "CD14", "LYZ", "PECAM1", "COL1A1", "HBA1"],
    },
    "Monocytes": {
        "positive": ["CD14", "LYZ", "LST1", "FCN1", "S100A8", "S100A9", "MS4A7", "CTSS", "C1QA", "C1QB", "C1QC", "LYVE1", "CD5L"],
        "negative": ["CD3D", "CD79A", "PECAM1", "COL1A1", "HBA1"],
    },
    "B_cells": {
        "positive": ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "CD74", "HLA-DRA", "RAG1", "RAG2", "FCER2", "IGLL1", "VPREB1", "DNTT", "JCHAIN", "TNFRSF17", "IGHG1", "IGHG2", "IGKC", "IGLC2"],
        "negative": ["CD3D", "CD14", "LYZ", "PECAM1", "COL1A1", "HBA1"],
    },
    "Red_blood_cells": {
        "positive": ["HBA1", "HBA2", "HBB", "HBD", "GYPA", "ALAS2", "PKLR", "GATA2", "TRIM58"],
        "negative": ["CD3D", "CD79A", "CD14", "PECAM1", "COL1A1"],
    },
    "Osteoblastic_lineage_cells": {
        "positive": ["COL1A1", "COL1A2", "SPARC", "BGLAP", "RUNX2", "ALPL", "IBSP", "SPP1", "PRSS23", "MXRA8", "TWIST1", "LMOD1"],
        "negative": ["PTPRC", "CD3D", "CD79A", "CD14", "PECAM1", "HBA1"],
    },
    "Endothelial_cells": {
        "positive": ["PECAM1", "VWF", "KDR", "CLDN5", "EMCN", "RAMP2", "ESAM", "GPIHBP1", "CCL14", "BCL6B"],
        "negative": ["PTPRC", "CD3D", "CD79A", "CD14", "COL1A1", "HBA1"],
    },
    "Endothelial_progenitor_cells": {
        "positive": ["KIT", "PROM1", "CD34", "FLT1", "KDR", "ENG"],
        "negative": ["CD3D", "CD79A", "CD14", "HBA1"],
    },
    "Plasmacytoid_dendritic_cells": {
        "positive": ["GZMB", "IL3RA", "CLEC4C", "LILRA4", "TCF4", "IRF7"],
        "negative": ["CD3D", "MS4A1", "CD14", "PECAM1", "COL1A1", "HBA1"],
    },
    "Mast_or_basophil_like": {
        "positive": ["TPSAB1", "TPSB2", "MS4A2", "CPA3", "KIT", "IL1RL1", "HDC"],
        "negative": ["CD3D", "CD79A", "CD14", "PECAM1", "COL1A1", "HBA1"],
    },
}

def genes_for_cluster(cluster: str):
    sub = markers[markers["cluster"].eq(cluster)].copy()
    if sub.empty:
        return []
    sub = sub.sort_values(["p_val_adj", "avg_log2FC"], ascending=[True, False])
    return sub["gene"].astype(str).tolist()

def present(genes, gene_set):
    genes_upper = {g.upper() for g in genes}
    return [g for g in gene_set if g.upper() in genes_upper]

def score_types(genes):
    rows = []
    for cell_type, items in canonical_sets.items():
        pos = present(genes[:50], items["positive"])
        top20_pos = present(genes[:20], items["positive"])
        neg = present(genes[:50], items["negative"])
        rows.append(
            {
                "cell_type": cell_type,
                "positive_hits_top50": pos,
                "positive_hits_top20": top20_pos,
                "negative_hits_top50": neg,
                "score": len(top20_pos) * 2 + len(pos) - len(neg),
            }
        )
    rows = sorted(rows, key=lambda x: x["score"], reverse=True)
    return rows

def curate_cluster(cluster: str, genes):
    scored = score_types(genes)
    best = scored[0]
    second = scored[1]
    top = genes[:15]

    gset = {g.upper() for g in genes[:30]}
    if {"JCHAIN", "TNFRSF17"} & gset and ({"IGHG1", "IGHG2", "IGKC", "IGLC2", "IGLV3-1", "IGKV4-1"} & gset):
        best_type = "Plasma_or_plasmablast_like"
        major = "B_cells"
        confidence = "high"
        note = "Immunoglobulin/JCHAIN/TNFRSF17 program; merge to B_cells for original broad label."
    elif {"IGLL1", "VPREB1", "DNTT"} & gset or ({"IGLL1", "VPREB1"} & gset and {"RAG1", "RAG2", "TMSB15A"} & gset):
        best_type = "B_lineage_immature_like"
        major = "B_cells"
        confidence = "high"
        note = "Immature/pre-B-like markers; merge to B_cells for original broad label."
    elif {"TPSAB1", "TPSB2"} & gset and ({"MS4A2", "IL1RL1", "CPA3", "KIT"} & gset):
        best_type = "Mast_or_basophil_like"
        major = "Other_myeloid_like_or_review"
        confidence = "high"
        note = "Mast/basophil-like markers are strong; this is not one of the original 9 broad FHTC labels and needs curator decision."
    elif {"GNLY", "KLRD1", "GZMH"} & gset and {"CD3D", "CD3E", "TRAC", "CD8A", "CD8B"} & gset:
        best_type = "Cytotoxic_T_or_NK_like"
        major = "T_cells"
        confidence = "high"
        note = "Cytotoxic lymphocyte program; merge to T_cells for original broad label unless NK is separated."
    elif {"CD40LG", "IL7R", "CCR6"} & gset and {"CD3D", "CD3E", "TRAC"} & gset:
        best_type = "CD4_T_like"
        major = "T_cells"
        confidence = "high"
        note = "T-cell markers with CD4/IL7R-like program."
    elif {"RAG1", "RAG2"} & gset and {"CD79A", "CD79B", "MS4A1", "BANK1", "CD74"} & gset:
        best_type = "B_lineage_immature_like"
        major = "B_cells"
        confidence = "high"
        note = "B-lineage markers with RAG genes; likely immature/pre-B-like."
    elif {"C1QA", "C1QB", "C1QC"} & gset and {"LYVE1", "CD5L", "IGSF21"} & gset:
        best_type = "Macrophage_C1Q_LYVE1_like"
        major = "Monocytes"
        confidence = "high"
        note = "Macrophage-like C1Q/LYVE1/CD5L program; merge to monocyte/macrophage broad label."
    else:
        best_type = best["cell_type"]
        major = best_type
        confidence = "high" if best["score"] >= second["score"] + 3 else "medium"
        note = "Curated from top50 marker overlap; review canonical dotplot/feature plots before final GMF reference."

    return {
        "curated_cell_type_fine_initial": best_type,
        "curated_cell_type_major_initial": major,
        "confidence_initial": confidence,
        "curation_note": note,
        "best_rule_type": best["cell_type"],
        "best_rule_score": best["score"],
        "second_rule_type": second["cell_type"],
        "second_rule_score": second["score"],
        "canonical_positive_hits": ";".join(best["positive_hits_top50"]),
        "canonical_positive_hits_top20": ";".join(best["positive_hits_top20"]),
        "canonical_negative_hits": ";".join(best["negative_hits_top50"]),
        "top15_markers": ";".join(top),
        "top50_markers": ";".join(genes[:50]),
    }

rows = []
for _, call in calls.sort_values("cluster", key=lambda s: s.astype(int)).iterrows():
    cluster = str(call["cluster"])
    genes = genes_for_cluster(cluster)
    cur = curate_cluster(cluster, genes)
    rows.append(
        {
            "resolution": 0.3,
            "cluster": cluster,
            "n_cells": int(call["n_cells"]),
            "auto_marker_score_type": call["best_marker_type"],
            "auto_marker_score_margin": call["marker_margin"],
            **cur,
            "review_status": "needs_user_review",
            "final_cell_type_fine": "",
            "final_cell_type_major": "",
            "final_include_in_GMF": "",
            "final_comment": "",
        }
    )

review = pd.DataFrame(rows)
review.to_csv(OUT / "GSE169396_res0p3_cluster_annotation_review_table.csv", index=False)

simple_cols = [
    "cluster",
    "n_cells",
    "top15_markers",
    "curated_cell_type_fine_initial",
    "curated_cell_type_major_initial",
    "confidence_initial",
    "curation_note",
]
review[simple_cols].to_csv(
    OUT / "GSE169396_res0p3_cluster_annotation_review_table_simple.csv",
    index=False,
)

print("Wrote:")
print(OUT / "GSE169396_res0p3_cluster_annotation_review_table.csv")
print(OUT / "GSE169396_res0p3_cluster_annotation_review_table_simple.csv")
print(review[simple_cols].to_string(index=False))
