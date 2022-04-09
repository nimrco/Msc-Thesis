import pandas as pd

window_size = 5


tips = pd.read_csv("tips.csv")
tips["gene state"] = 0
tips["pseudo state"] = 0

cluster_df = pd.read_csv("cluster_2634.csv")
for index, row in tips.iterrows():
    if cluster_df.iloc[row["x"]]["pseudo"]:
        tips.iloc[index]["pseudo state"] = 6
        for i in range(1, window_size + 1):
            forward = index + i
            if forward < len(tips):
                strain = tips.iloc[forward]["x"]
                if cluster_df.iloc[strain]["pseudo"]:
                    tips.iloc[index]["pseudo state"] = 4
            reverse = index - i
            if reverse >= 0:
                strain = tips.iloc[reverse]["x"]
                if cluster_df.iloc[strain]["pseudo"]:
                    tips.iloc[index]["pseudo state"] = 4
    if cluster_df.iloc[row["x"]]["gene"]:
        tips.iloc[index]["gene state"] = 3
        for i in range(1, window_size + 1):
            forward = index + i
            if forward < len(tips):
                strain = tips.iloc[forward]["x"]
                if cluster_df.iloc[strain]["gene"]:
                    tips.iloc[index]["gene state"] = 1
                if cluster_df.iloc[strain]["pseudo"]:
                    tips.iloc[index]["gene state"] = 2
            reverse = index - i
            if reverse >= 0:
                strain = tips.iloc[reverse]["x"]
                if cluster_df.iloc[strain]["gene"]:
                    tips.iloc[index]["gene state"] = 1
                if cluster_df.iloc[strain]["pseudo"]:
                    tips.iloc[index]["gene state"] = 2

tips.to_csv("cluster_2634_counts.csv", index=False)
print("script done")



