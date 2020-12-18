import os
import numpy as np
import geopandas as gpd
from tqdm import tqdm
import pulp as p
from shapely.geometry import Point
import matplotlib.pyplot as plt

print("Start")

gdf = gpd.read_file(os.path.join(os.getcwd(), "data", "gemeente_2020_v1_mini_subset.shp"))
gdf["geom_size"] = np.sqrt(gdf.geometry.area) / 2

gdf["_neighbors"] = ""
for index, row in tqdm(gdf.iterrows(), desc="Calculating neighbors"):
    neighbors = gdf[
        ~gdf.geometry.disjoint(row["geometry"])
    ].index.tolist()
    gdf.at[index, "_neighbors"] = ",".join(
        [str(n) for n in neighbors if n != index]
    )
    for n in gdf.at[index, "_neighbors"].split(","):
        print(index, gdf.at[index, "GM_NAAM"], n, gdf.at[int(n), "GM_NAAM"])
    gdf.at[index, "_n_neighbors"] = len(neighbors)

gdf["_n_neighbors"].fillna(0, inplace=True)

print(gdf.head())

problem = p.LpProblem("Problem", p.LpMinimize)

xi = [i for i in range(gdf.shape[0])]
xj = [i for i in range(gdf.shape[0])]


import itertools


def flat_list(t):
    return [item for sublist in t for item in sublist]


hv = []
permutations = [tuple(sorted(x)) for x in flat_list([list(zip(x_i, xj)) for x_i in itertools.permutations(xi,len(xj))]) if x[0] != x[1]]
vars_list = []
for per in permutations:
    i, j = per
    vars_list.append(f"x_{i}")
    vars_list.append(f"x_{j}")
    vars_list.append(f"y_{i}")
    vars_list.append(f"y_{j}")
    vars_list.append(f"h_{i}_{j}")
    vars_list.append(f"v_{i}_{j}")
    vars_list.append(f"dx_{i}")
    vars_list.append(f"dy_{j}")

vars_dict = p.LpVariable.dicts("coords", vars_list)

dist_list = []
for vd in vars_dict:
    vds = vd.split("_")
    if vds[0] == "d":
        _,i = vd.split("_")
        d = vars_dict[vd]

        problem += d >= 0
        if vd[1] == "x":
            dist = gdf.at[int(i), "geometry"].centroid.x - d
        else:
            dist = gdf.at[int(i), "geometry"].centroid.y - d
        dist_list.append(dist)

    elif len(vds) > 2:


problem += p.lpSum(dist_list)

problem.solve()

# print(problem)

for v in problem.variables():
    var_name = v.name.replace("coords_", "").split("_")
    if len(var_name) == 2:
        # All variables with two parts in their name are coordinates. The first part says if it's a x or y
        # coordinate, the second part the geometry it corresponds to (row number in Geodataframe).
        c, i = var_name
        gdf.at[int(i), f"_{c}"] = v.varValue

gdf.geometry = gdf.apply(lambda r: Point(r["_x"], r["_y"]).buffer(r["geom_size"]), axis=1)


gdf.plot()
plt.show()
