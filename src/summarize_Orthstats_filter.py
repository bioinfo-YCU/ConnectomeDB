# summarize counts depending on filter
import createDataTable_perSpecies
def summarize_orthologs(human_col, species_col, label,
                        confidence_orth_col=None, confidence_orth_threshold=None,
                        GOC_col=None, GOC_threshold=None):
    df = mouse_gene_pair1.copy()

    if confidence_orth_col and confidence_orth_threshold is not None:
        df = df[df[confidence_orth_col] == confidence_orth_threshold]

    if GOC_col and GOC_threshold is not None:
        df = df[df[GOC_col] >= GOC_threshold]  # Use >= instead of ==

    unique_pairs = df[[human_col, species_col]].drop_duplicates()

    counts = (
        unique_pairs
        .groupby(human_col)[species_col]
        .count()
        .sort_values(ascending=False)
        .reset_index(name='count')
    )

    filter_tag = label.lower()
    if confidence_orth_threshold is not None:
        filter_tag += f"_conf{confidence_orth_threshold}"
    if GOC_threshold is not None:
        filter_tag += f"_GOCge{GOC_threshold}"

    counts.to_csv(f"data/human_mouse_orth_count_{filter_tag}.csv", index=False)

    summary_counts = counts['count'].value_counts().sort_index()
    total_human_genes = counts.shape[0]

    summary_lines = [
        f"Out of {total_human_genes} unique human {label.lower()} genes "
        f"(Orthology Confidence = {confidence_orth_threshold}, GOC ≥ {GOC_threshold}):"
    ]
    for orth_count, gene_count in summary_counts.items():
        summary_lines.append(
            f" - {gene_count} human {label.lower()} genes had {orth_count} mouse ortholog(s)"
        )

    return "\n".join(summary_lines)


# Detect columns
confidence_orth_ligand = [col for col in mouse_gene_pair1.columns if "Ligand Orthology Confidence" in col][0]
GOC_col_ligand = [col for col in mouse_gene_pair1.columns if "Ligand GOC" in col][0]

confidence_orth_receptor = [col for col in mouse_gene_pair1.columns if "Receptor Orthology Confidence" in col][0]
GOC_col_receptor = [col for col in mouse_gene_pair1.columns if "Receptor GOC" in col][0]

# Generate summaries
ligand_summary = summarize_orthologs(
    human_ligand_col, ligand_col, "Ligand",
    confidence_orth_col=confidence_orth_ligand, confidence_orth_threshold=None,
    GOC_col=GOC_col_ligand, GOC_threshold=0
)

receptor_summary = summarize_orthologs(
    human_receptor_col, receptor_col, "Receptor",
    confidence_orth_col=confidence_orth_receptor, confidence_orth_threshold=None,
    GOC_col=GOC_col_receptor, GOC_threshold=0
)

# Print
print(ligand_summary)
print()
print(receptor_summary)