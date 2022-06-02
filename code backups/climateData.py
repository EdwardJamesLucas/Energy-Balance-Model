import pandas as pd
import pickle


dfCityNames = pd.read_excel(
    "climate_data.xlsx", index_col=0, header=None, usecols=[2], nrows=1
)

# andreas code
N_ROWS_FOR_EACH_DF = 3
cities_df_first_row = pd.read_excel("climate_data.xlsx", nrows=1)
cities_no_num = [
    x for x in cities_df_first_row.columns if not "Unnamed" in str(x) and type(x) is str
]

cities_dfs_dict = {}
for number, city in enumerate(cities_no_num):
    cities_dfs_dict[city] = pd.read_excel(
        "climate_data.xlsx",
        skiprows=1,
        index_col=0,
        usecols=[
            0,
            (number * N_ROWS_FOR_EACH_DF) + 1,
            (number * N_ROWS_FOR_EACH_DF) + 2,
            (number * N_ROWS_FOR_EACH_DF) + 3,
        ],
    )

for city in cities_dfs_dict:
    cols = cities_dfs_dict[city].columns
    if all("." in col for col in cols):
        cols = [col.split(".")[0] for col in cols]
        cities_dfs_dict[city].columns = cols

# save dictionary to pickle file
with open("climate_data.pickle", "wb") as file:
    pickle.dump(cities_dfs_dict, file, protocol=pickle.HIGHEST_PROTOCOL)
