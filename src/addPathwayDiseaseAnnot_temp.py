# Add top pathway per pair
# LR_pairs = gene_pair["Human LR Pair"].unique()
# df= pd.read_csv("data/pathway_annotations_per_pair.csv")

#df = df[df["interaction"].isin(LR_pairs)]

# Sort by absolute value of 'weight', descending (larger abs(weight) first)
# df_sorted = df.reindex(df['weight'].abs().sort_values(ascending=False).index)
# Keep only the first occurrence for each unique 'interaction'
#df_unique = df_sorted.drop_duplicates(subset='interaction', keep='first')
#Keep ALL
# df = df_sorted.reset_index(drop=True)
# top_pathway_df = df[["interaction", "source"]]
# top_pathway_df = top_pathway_df.groupby('interaction')['source'].apply(', '.join).reset_index()
# top_pathway_df = top_pathway_df.rename(columns={
#                                       "source": "PROGENy Pathway"
# })
# top_pathway_df["interaction"] = [value.replace("^", " ") for value in top_pathway_df["interaction"]]
# gene_pair = gene_pair.merge(top_pathway_df, how='left', left_on='Human LR Pair', right_on='interaction')
# gene_pair = gene_pair.drop(columns=["interaction"])
# #df = df_unique.reset_index(drop=True)
# top_pathway_df=fetchGSheet.kegg_pathway_info[["LR Pair", "kegg_pathway_id", "kegg_relationship", "kegg_pathway_name"]].copy()
# # add link to kegg_pathway_name
# top_pathway_df["kegg_pathway_name"] = [
#     f'<a href="https://www.kegg.jp/pathway/{kegg_id}" target="_blank">{name}</a>'
#     for kegg_id, name in zip(top_pathway_df["kegg_pathway_id"], top_pathway_df["kegg_pathway_name"])
# ]

# # link to kegg_pathway_id
# top_pathway_df["kegg_pathway_id"] = [
#     f'<a href="https://www.kegg.jp/pathway/{id}" target="_blank">{id}</a>'
#     for id in top_pathway_df["kegg_pathway_id"]]


# top_pathway_df = top_pathway_df.rename(columns={
#                                       "kegg_pathway_name": "KEGG Pathway",
#                                       "kegg_relationship": "KEGG relationship",
#                                       "kegg_pathway_id": "KEGG Pathway ID"
    
# })
# top_pathway_df1 = top_pathway_df[["LR Pair", "KEGG Pathway"]].drop_duplicates()
# top_pathway_df1 = top_pathway_df1.groupby('LR Pair')['KEGG Pathway'].apply(', '.join).reset_index()
# gene_pair = gene_pair.merge(top_pathway_df1, how='left', left_on='Human LR Pair', right_on='LR Pair')
# gene_pair = gene_pair.drop(columns=["LR Pair"])

# As of latest DB, skip for now
# # Add Disease Category per pair
# df= pd.read_csv("data/disease_annotations_per_pair.csv")
# df_cat=pd.read_csv("data/disease_categories.csv")
# mapping = dict(zip(df_cat['Disease Name'], df_cat['Category']))
# # Replace values in the column based on the mapping
# df["Disease Type"] = df['disease'].replace(mapping)
# df = df[["interaction", "Disease Type"]].drop_duplicates()
# df['Disease Type'] = df['Disease Type'].astype(str)
# df = df.sort_values(by='Disease Type', ascending=True)
# # Group by 'col1' and combine 'col2' values with ', '
# df = df.groupby('interaction')['Disease Type'].apply(', '.join).reset_index()

# # Create "Cancer-related" column based on whether "Cancers & Neoplasms" is in col2
# df['Cancer-related'] = df['Disease Type'].apply(lambda x: 'Yes' if 'Cancer' in x else 'No')
# disease_df = df[df["interaction"].isin(LR_pairs)]
# # Function to update the "Cancer-related" column and modify "col2" if needed

# gene_pair = gene_pair.merge(disease_df, how='left', left_on='Human LR Pair', right_on='interaction')