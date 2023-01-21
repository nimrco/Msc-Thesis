import pandas as pd

random_df = pd.read_csv("random_groups_15.csv")
random_df['category'] = 'random'
random_df = random_df[['CV', 'category']]
tree_df = pd.read_csv("tree_window_15.csv")
tree_df['category'] = 'tree'
tree_df = tree_df[['CV', 'category']]
mlst_df = pd.read_csv("mlst_stats.csv")
mlst_df = mlst_df[mlst_df['MLST count'] >= 5]
mlst_df['category'] = 'mlst'
mlst_df = mlst_df[['MLST CV of pseudogenes', 'category']]
mlst_df.rename(columns={'MLST CV of pseudogenes': 'CV'}, inplace=True)
bioproject_df = pd.read_csv("bioproject_stats.csv")
bioproject_df = bioproject_df[bioproject_df['BioProject count'] >= 5]
bioproject_df['category'] = 'bioproject'
bioproject_df = bioproject_df[['Bio cv', 'category']]
bioproject_df.rename(columns={'Bio cv': 'CV'}, inplace=True)
combined_df = pd.concat([random_df, tree_df, mlst_df, bioproject_df], ignore_index=True)
combined_df.to_csv("combined_cv_15.csv", index=False)
print("script done")
