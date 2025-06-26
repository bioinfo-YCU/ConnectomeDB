# summarize counts depending on filter
species = "horse"  # Change to desired species, e.g., "zebrafish", "sheep"
capital_species = species.capitalize()

# Load final result file
final_result = pd.read_csv(f"data/human_{species}_merged_ensemblBiomaRt_inParanoid.csv")

# Define score columns dynamically
score_cols = [
    f"Ligand_human_inparalog_score",
    f"Receptor_human_inparalog_score",
    f"Ligand_{species}_inparalog_score",
    f"Receptor_{species}_inparalog_score",
    f"Ligand_human_seed_score",
    f"Receptor_human_seed_score",
    f"Ligand_{species}_seed_score",
    f"Receptor_{species}_seed_score",
    f"Ligand_bitscore",
    f"Receptor_bitscore"
]

for col in score_cols:
    if col in final_result.columns:
        final_result[col] = pd.to_numeric(final_result[col], errors='coerce')

# Detect columns
confidence_orth_ligand = [col for col in final_result.columns if "Ligand Orthology Confidence" in col][0]
GOC_col_ligand = [col for col in final_result.columns if "Ligand GOC" in col][0]
percIdent_col_ligand = [col for col in final_result.columns if "Ligand % Identity" in col][0]
human_ligand_col = [col for col in final_result.columns if "Ligand HGNC ID" in col][0]
ligand_col = [col for col in final_result.columns if "Ligand" in col][0]
ligand_human_inparalog_score = [col for col in final_result.columns if "Ligand_human_inparalog_score" in col][0]
ligand_species_inparalog_score = [col for col in final_result.columns if f"Ligand_{species}_inparalog_score" in col][0]
ligand_human_seed_score = [col for col in final_result.columns if "Ligand_human_seed_score" in col][0]
ligand_species_seed_score = [col for col in final_result.columns if f"Ligand_{species}_seed_score" in col][0]
ligand_bit_score = [col for col in final_result.columns if "Ligand_bitscore" in col][0]

human_receptor_col = [col for col in final_result.columns if "Receptor HGNC ID" in col][0]
confidence_orth_receptor = [col for col in final_result.columns if "Receptor Orthology Confidence" in col][0]
GOC_col_receptor = [col for col in final_result.columns if "Receptor GOC" in col][0]
percIdent_col_receptor = [col for col in final_result.columns if "Receptor % Identity" in col][0]
receptor_col = [col for col in final_result.columns if "Receptor" in col][0]
receptor_human_inparalog_score = [col for col in final_result.columns if "Receptor_human_inparalog_score" in col][0]
receptor_species_inparalog_score = [col for col in final_result.columns if f"Receptor_{species}_inparalog_score" in col][0]
receptor_human_seed_score = [col for col in final_result.columns if "Receptor_human_seed_score" in col][0]
receptor_species_seed_score = [col for col in final_result.columns if f"Receptor_{species}_seed_score" in col][0]
receptor_bit_score = [col for col in final_result.columns if "Receptor_bitscore" in col][0]

# Define function

def summarize_orthologs(human_col, species_col, label,
                        confidence_orth_col=None, confidence_orth_threshold=None,
                        GOC_col=None, GOC_threshold=None,
                        perc_identity_col=None, perc_identity_thres=None,
                        ligand_human_inparalog_score_col=None, ligand_human_inparalog_score_threshold=None,
                        receptor_human_inparalog_score_col=None, receptor_human_inparalog_score_threshold=None,
                        ligand_species_inparalog_score_col=None, ligand_species_inparalog_score_threshold=None,
                        receptor_species_inparalog_score_col=None, receptor_species_inparalog_score_threshold=None,
                        ligand_human_seed_score_col=None, ligand_human_seed_score_threshold=None,
                        receptor_human_seed_score_col=None, receptor_human_seed_score_threshold=None,
                        ligand_species_seed_score_col=None, ligand_species_seed_score_threshold=None,
                        receptor_species_seed_score_col=None, receptor_species_seed_score_threshold=None,
                        ligand_bit_score_col=None, ligand_bit_score_threshold=None,
                        receptor_bit_score_col=None, receptor_bit_score_threshold=None):

    df = final_result.copy()

    filters = [
        (confidence_orth_col, lambda x: x == confidence_orth_threshold),
        (GOC_col, lambda x: x >= GOC_threshold),
        (perc_identity_col, lambda x: x >= perc_identity_thres),
        (ligand_human_inparalog_score_col, lambda x: x >= ligand_human_inparalog_score_threshold),
        (receptor_human_inparalog_score_col, lambda x: x >= receptor_human_inparalog_score_threshold),
        (ligand_species_inparalog_score_col, lambda x: x >= ligand_species_inparalog_score_threshold),
        (receptor_species_inparalog_score_col, lambda x: x >= receptor_species_inparalog_score_threshold),
        (ligand_human_seed_score_col, lambda x: x >= ligand_human_seed_score_threshold),
        (receptor_human_seed_score_col, lambda x: x >= receptor_human_seed_score_threshold),
        (ligand_species_seed_score_col, lambda x: x >= ligand_species_seed_score_threshold),
        (receptor_species_seed_score_col, lambda x: x >= receptor_species_seed_score_threshold),
        (ligand_bit_score_col, lambda x: x >= ligand_bit_score_threshold),
        (receptor_bit_score_col, lambda x: x >= receptor_bit_score_threshold),
    ]

    original_rows = df.shape[0]
    for col, condition in filters:
        if col and condition is not None:
            before = df.shape[0]
            df = df[df[col].apply(condition)]
            after = df.shape[0]
            print(f"Filtered {col}: {before - after} rows removed (remaining: {after})")

    unique_pairs = df[[human_col, species_col]].drop_duplicates()

    counts = (
        unique_pairs
        .groupby(human_col)[species_col]
        .count()
        .sort_values(ascending=False)
        .reset_index(name='count')
    )

    filter_tag = f"{label.lower()}"
    if confidence_orth_threshold is not None:
        filter_tag += f"_conf{confidence_orth_threshold}"
    if GOC_threshold is not None:
        filter_tag += f"_GOCge{GOC_threshold}"
    if ligand_human_inparalog_score_threshold is not None:
        filter_tag += f"_LHISge{ligand_human_inparalog_score_threshold}"
    if receptor_bit_score_threshold is not None:
        filter_tag += f"_RBSge{receptor_bit_score_threshold}"

    #counts.to_csv(f"data/human_{species}_orth_count_{filter_tag}.csv", index=False)

    summary_counts = counts['count'].value_counts().sort_index()
    total_human_genes = counts.shape[0]

    summary_lines = [
        f"Out of {total_human_genes} unique human {label.lower()} genes:",
        f" - Filters applied: " + "; ".join([
            f"{col} ≥ {threshold}" for col, threshold in [
                (confidence_orth_col, confidence_orth_threshold),
                (GOC_col, GOC_threshold),
                (perc_identity_col, perc_identity_thres),
                (ligand_human_inparalog_score_col, ligand_human_inparalog_score_threshold),
                (receptor_human_inparalog_score_col, receptor_human_inparalog_score_threshold),
                (ligand_species_inparalog_score_col, ligand_species_inparalog_score_threshold),
                (receptor_species_inparalog_score_col, receptor_species_inparalog_score_threshold),
                (ligand_human_seed_score_col, ligand_human_seed_score_threshold),
                (receptor_human_seed_score_col, receptor_human_seed_score_threshold),
                (ligand_species_seed_score_col, ligand_species_seed_score_threshold),
                (receptor_species_seed_score_col, receptor_species_seed_score_threshold),
                (ligand_bit_score_col, ligand_bit_score_threshold),
                (receptor_bit_score_col, receptor_bit_score_threshold),
            ] if threshold is not None and col is not None
        ])
    ]

    for orth_count, gene_count in summary_counts.items():
        summary_lines.append(
            f" - {gene_count} human {label.lower()} genes had {orth_count} {species} ortholog(s)"
        )

    return "\n".join(summary_lines)

# Ligand
ligand_summary = summarize_orthologs(
    human_col=human_ligand_col,
    species_col=ligand_col,
    label="Ligand",
    # GOC_col=GOC_col_ligand,
    # GOC_threshold=25,
    # confidence_orth_col=confidence_orth_ligand,
    # confidence_orth_threshold=1,
    # perc_identity_col = percIdent_col_ligand,
    # perc_identity_thres = 60
    # ligand_human_inparalog_score_col=ligand_human_inparalog_score,
    # ligand_human_inparalog_score_threshold=1,
    # ligand_bit_score_col=ligand_bit_score,
    # ligand_bit_score_threshold=40,
    # ligand_species_inparalog_score_col=ligand_species_inparalog_score,
    # ligand_species_inparalog_score_threshold=1,
)
print(ligand_summary)

# Receptor
receptor_summary = summarize_orthologs(
    human_col=human_receptor_col,
    species_col=receptor_col,
    label="Receptor",
    # GOC_col=GOC_col_receptor,
    # GOC_threshold=25,
    # confidence_orth_col=confidence_orth_receptor,
    # confidence_orth_threshold=1,
    # perc_identity_col = percIdent_col_receptor,
    # perc_identity_thres = 60
    # receptor_human_inparalog_score_col=receptor_human_inparalog_score,
    # receptor_human_inparalog_score_threshold=1,
    # receptor_bit_score_col=receptor_bit_score,
    # receptor_bit_score_threshold=40,
    # receptor_species_inparalog_score_col=receptor_species_inparalog_score,
    # receptor_species_inparalog_score_threshold=1,
)
print(receptor_summary)

