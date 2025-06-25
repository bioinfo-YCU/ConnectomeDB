# summarize counts depending on filter
score_cols = [
    ligand_human_inparalog_score,
    receptor_human_inparalog_score,
    ligand_mouse_inparalog_score,
    receptor_mouse_inparalog_score,
    ligand_human_seed_score,
    receptor_human_seed_score,
    ligand_mouse_seed_score,
    receptor_mouse_seed_score,
    ligand_bit_score,
    receptor_bit_score
]

for col in score_cols:
    if col in final_result.columns:
        final_result[col] = pd.to_numeric(final_result[col], errors='coerce')

def summarize_orthologs(human_col, species_col, label,
                        confidence_orth_col=None, confidence_orth_threshold=None,
                        GOC_col=None, GOC_threshold=None,
                        perc_identity_col= None, perc_identity_thres= None,
                        ligand_human_inparalog_score_col=None, ligand_human_inparalog_score_threshold=None,
                        receptor_human_inparalog_score_col=None, receptor_human_inparalog_score_threshold=None,
                        ligand_mouse_inparalog_score_col=None, ligand_mouse_inparalog_score_threshold=None,
                        receptor_mouse_inparalog_score_col=None, receptor_mouse_inparalog_score_threshold=None,
                        ligand_human_seed_score_col=None, ligand_human_seed_score_threshold=None,
                        receptor_human_seed_score_col=None, receptor_human_seed_score_threshold=None,
                        ligand_mouse_seed_score_col=None, ligand_mouse_seed_score_threshold=None,
                        receptor_mouse_seed_score_col=None, receptor_mouse_seed_score_threshold=None,
                        ligand_bit_score_col=None, ligand_bit_score_threshold=None,
                        receptor_bit_score_col=None, receptor_bit_score_threshold=None):
    
    df = final_result.copy()

    # Apply filters one by one if thresholds are given
    filters = [
        (confidence_orth_col, lambda x: x == confidence_orth_threshold),
        (GOC_col, lambda x: x >= GOC_threshold),
        (perc_identity_col, lambda x: x >= perc_identity_thres),
        (ligand_human_inparalog_score_col, lambda x: x >= ligand_human_inparalog_score_threshold),
        (receptor_human_inparalog_score_col, lambda x: x >= receptor_human_inparalog_score_threshold),
        (ligand_mouse_inparalog_score_col, lambda x: x >= ligand_mouse_inparalog_score_threshold),
        (receptor_mouse_inparalog_score_col, lambda x: x >= receptor_mouse_inparalog_score_threshold),
        (ligand_human_seed_score_col, lambda x: x >= ligand_human_seed_score_threshold),
        (receptor_human_seed_score_col, lambda x: x >= receptor_human_seed_score_threshold),
        (ligand_mouse_seed_score_col, lambda x: x >= ligand_mouse_seed_score_threshold),
        (receptor_mouse_seed_score_col, lambda x: x >= receptor_mouse_seed_score_threshold),
        (ligand_bit_score_col, lambda x: x >= ligand_bit_score_threshold),
        (receptor_bit_score_col, lambda x: x >= receptor_bit_score_threshold),
    ]


    # Apply filters and track how many rows were removed
    original_rows = df.shape[0]
    for col, condition in filters:
        if col and condition is not None:
            before = df.shape[0]
            df = df[df[col].apply(condition)]
            after = df.shape[0]
            print(f"Filtered {col}: {before - after} rows removed (remaining: {after})")

    # Compute unique ortholog pairs
    unique_pairs = df[[human_col, species_col]].drop_duplicates()

    # Count mouse orthologs per human gene
    counts = (
        unique_pairs
        .groupby(human_col)[species_col]
        .count()
        .sort_values(ascending=False)
        .reset_index(name='count')
    )

    # Build tag from filters
    filter_tag = label.lower()
    if confidence_orth_threshold is not None:
        filter_tag += f"_conf{confidence_orth_threshold}"
    if GOC_threshold is not None:
        filter_tag += f"_GOCge{GOC_threshold}"
    if ligand_human_inparalog_score_threshold is not None:
        filter_tag += f"_LHISge{ligand_human_inparalog_score_threshold}"
    if receptor_bit_score_threshold is not None:
        filter_tag += f"_RBSge{receptor_bit_score_threshold}"

    counts.to_csv(f"data/human_mouse_orth_count_{filter_tag}.csv", index=False)

    summary_counts = counts['count'].value_counts().sort_index()
    total_human_genes = counts.shape[0]

    # Collect active filters for reporting
    active_filters = []
    if confidence_orth_col and confidence_orth_threshold is not None:
        if isinstance(confidence_orth_threshold, (list, set, tuple)):
            active_filters.append(f"{confidence_orth_col} in {sorted(confidence_orth_threshold)}")
        else:
            active_filters.append(f"{confidence_orth_col} equals '{confidence_orth_threshold}'")
    if GOC_col and GOC_threshold is not None:
        active_filters.append(f"{GOC_col} ≥ {GOC_threshold}")
    if perc_identity_col and perc_identity_thres is not None:
        active_filters.append(f"{perc_identity_col} ≥ {perc_identity_thres}")
    if ligand_human_inparalog_score_col and ligand_human_inparalog_score_threshold is not None:
        active_filters.append(f"{ligand_human_inparalog_score_col} ≥ {ligand_human_inparalog_score_threshold}")
    if receptor_human_inparalog_score_col and receptor_human_inparalog_score_threshold is not None:
        active_filters.append(f"{receptor_human_inparalog_score_col} ≥ {receptor_human_inparalog_score_threshold}")
    if ligand_mouse_inparalog_score_col and ligand_mouse_inparalog_score_threshold is not None:
        active_filters.append(f"{ligand_mouse_inparalog_score_col} ≥ {ligand_mouse_inparalog_score_threshold}")
    if receptor_mouse_inparalog_score_col and receptor_mouse_inparalog_score_threshold is not None:
        active_filters.append(f"{receptor_mouse_inparalog_score_col} ≥ {receptor_mouse_inparalog_score_threshold}")
    if ligand_human_seed_score_col and ligand_human_seed_score_threshold is not None:
        active_filters.append(f"{ligand_human_seed_score_col} ≥ {ligand_human_seed_score_threshold}")
    if receptor_human_seed_score_col and receptor_human_seed_score_threshold is not None:
        active_filters.append(f"{receptor_human_seed_score_col} ≥ {receptor_human_seed_score_threshold}")
    if ligand_mouse_seed_score_col and ligand_mouse_seed_score_threshold is not None:
        active_filters.append(f"{ligand_mouse_seed_score_col} ≥ {ligand_mouse_seed_score_threshold}")
    if receptor_mouse_seed_score_col and receptor_mouse_seed_score_threshold is not None:
        active_filters.append(f"{receptor_mouse_seed_score_col} ≥ {receptor_mouse_seed_score_threshold}")
    if ligand_bit_score_col and ligand_bit_score_threshold is not None:
        active_filters.append(f"{ligand_bit_score_col} ≥ {ligand_bit_score_threshold}")
    if receptor_bit_score_col and receptor_bit_score_threshold is not None:
        active_filters.append(f"{receptor_bit_score_col} ≥ {receptor_bit_score_threshold}")
    
    filter_text = "; ".join(active_filters) if active_filters else "No filters applied"
    
    summary_lines = [
        f"Out of {total_human_genes} unique human {label.lower()} genes:",
        f" - Filters applied: {filter_text}"
    ]

    for orth_count, gene_count in summary_counts.items():
        summary_lines.append(
            f" - {gene_count} human {label.lower()} genes had {orth_count} mouse ortholog(s)"
        )

    return "\n".join(summary_lines)


# Detect columns
# Ligand
confidence_orth_ligand = [col for col in final_result.columns if "Ligand Orthology Confidence" in col][0]
GOC_col_ligand = [col for col in final_result.columns if "Ligand GOC" in col][0]
percIdent_col_ligand = [col for col in final_result.columns if "Ligand % Identity" in col][0]
human_ligand_col = [col for col in final_result.columns if "Ligand HGNC ID" in col][0]
ligand_col = [col for col in final_result.columns if "Ligand" in col][0]
ligand_col = [col for col in final_result.columns if "Ligand" in col][0]
ligand_human_inparalog_score = [col for col in final_result.columns if "Ligand_human_inparalog_score" in col][0]
ligand_mouse_inparalog_score = [col for col in final_result.columns if "Ligand_mouse_inparalog_score" in col][0]
ligand_human_seed_score = [col for col in final_result.columns if "Ligand_human_seed_score" in col][0]
ligand_mouse_seed_score = [col for col in final_result.columns if "Ligand_mouse_seed_score" in col][0]
ligand_bit_score = [col for col in final_result.columns if "Ligand_bitscore" in col][0]

#Receptor
human_receptor_col = [col for col in final_result.columns if "Receptor HGNC ID" in col][0]
confidence_orth_receptor = [col for col in final_result.columns if "Receptor Orthology Confidence" in col][0]
GOC_col_receptor = [col for col in final_result.columns if "Receptor GOC" in col][0]
percIdent_col_receptor = [col for col in final_result.columns if "Receptor % Identity" in col][0]
receptor_col = [col for col in final_result.columns if "Receptor" in col][0]
receptor_human_inparalog_score = [col for col in final_result.columns if "Receptor_human_inparalog_score" in col][0]
receptor_mouse_inparalog_score = [col for col in final_result.columns if "Receptor_mouse_inparalog_score" in col][0]
receptor_human_seed_score = [col for col in final_result.columns if "Receptor_human_seed_score" in col][0]
receptor_mouse_seed_score = [col for col in final_result.columns if "Receptor_mouse_seed_score" in col][0]
receptor_bit_score = [col for col in final_result.columns if "Receptor_bitscore" in col][0]


ligand_summary = summarize_orthologs(
    human_col=human_ligand_col,
    species_col=ligand_col,
    label="Ligand",
    # confidence_orth_col=confidence_orth_ligand,
    # confidence_orth_threshold=1,
    # GOC_col=GOC_col_ligand,
    # GOC_threshold=60,
    ligand_human_inparalog_score_col=ligand_human_inparalog_score,
    ligand_human_inparalog_score_threshold=1,
    # ligand_human_seed_score_col=ligand_human_seed_score,
    # ligand_human_seed_score_threshold=1,
    # ligand_mouse_seed_score_col=ligand_mouse_seed_score,
    # ligand_mouse_seed_score_threshold=1,
    ligand_bit_score_col=ligand_bit_score,
    ligand_bit_score_threshold=40
)

print(ligand_summary)

receptor_summary = summarize_orthologs(
    human_col=human_receptor_col,
    species_col=receptor_col,
    label="Receptor",
    # confidence_orth_col=confidence_orth_receptor,
    # confidence_orth_threshold=1,
    # GOC_col=GOC_col_receptor,
    # GOC_threshold=60,
    receptor_human_inparalog_score_col=receptor_human_inparalog_score,
    receptor_human_inparalog_score_threshold=1,
    # receptor_human_seed_score_col=receptor_human_seed_score,
    # receptor_human_seed_score_threshold=1,
    # receptor_mouse_seed_score_col=receptor_mouse_seed_score,
    # receptor_mouse_seed_score_threshold=1,
    receptor_bit_score_col=receptor_bit_score,
    receptor_bit_score_threshold=40
)

print(receptor_summary)
