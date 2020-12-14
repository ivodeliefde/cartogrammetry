import os
import geopandas as gpd
import pulp as p
import matplotlib.pyplot as plt
from cartogrammetry.create import create_block
from shapely.geometry import Point

shp_path = os.path.join(os.getcwd(), "data", "gemeente_2020_v1_mini_subset.shp")
municipalities = gpd.read_file(shp_path)
gdf = municipalities
print(municipalities.head())

x1 = municipalities.loc[0, "geometry"].centroid.x
y1 = municipalities.loc[0, "geometry"].centroid.y
x2 = municipalities.loc[1, "geometry"].centroid.x
y2 = municipalities.loc[1, "geometry"].centroid.y
x3 = municipalities.loc[2, "geometry"].centroid.x
y3 = municipalities.loc[2, "geometry"].centroid.y

x_1 = p.LpVariable("x_1", lowBound=0)
y_1 = p.LpVariable("y_1", lowBound=0)
x_2 = p.LpVariable("x_2", lowBound=0)
y_2 = p.LpVariable("y_2", lowBound=0)
x_3 = p.LpVariable("x_3", lowBound=0)
y_3 = p.LpVariable("y_3", lowBound=0)

x_dist_1 = x_1 - x1
y_dist_1 = y_1 - y1
x_dist_2 = x_2 - x2
y_dist_2 = y_2 - y2
x_dist_3 = x_3 - x3
y_dist_3 = y_3 - y3

dist_x12 = x_1 - x_2
dist_x21 = x_2 - x_1
dist_y12 = y_1 - y_2
dist_y21 = y_2 - y_1
dist_x13 = x_1 - x_3
dist_x31 = x_3 - x_1
dist_y13 = y_1 - y_3
dist_y31 = y_3 - y_1
dist_x23 = x_2 - x_3
dist_x32 = x_3 - x_2
dist_y23 = y_2 - x_3
dist_y32 = y_3 - x_2

gdf["geom_size"] = [1030, 450, 150]
w12 = gdf.iloc[0, :].geom_size + gdf.iloc[1, :].geom_size
w13 = gdf.iloc[0, :].geom_size + gdf.iloc[2, :].geom_size
w23 = gdf.iloc[1, :].geom_size + gdf.iloc[2, :].geom_size

gap12 = 0
gap13 = 0
gap23 = 0

problem = p.LpProblem("Problem", p.LpMinimize)
problem += (
    # x_dist_1 +
    # y_dist_1 +
    # x_dist_2 +
    # y_dist_2 +
    # x_dist_3 +
    # y_dist_3 +
    dist_x12
    + dist_x13
    + dist_x23
    + dist_y12
    + dist_y13
    + dist_y23
)

problem += dist_x12 >= w12 + gap12
problem += dist_y12 >= w12 + gap12
problem += dist_x13 >= w13 + gap13
problem += dist_y13 >= w13 + gap13
problem += dist_x23 >= w23 + gap23
problem += dist_y23 >= w23 + gap23
problem += dist_x12 >= dist_x12 - w12
problem += dist_y12 >= dist_y12 - w12
problem += dist_x13 >= dist_x13 - w13
problem += dist_y13 >= dist_y13 - w13
problem += dist_x23 >= dist_x23 - w23
problem += dist_y23 >= dist_y23 - w23
problem += dist_x12 >= dist_x21 - w12
problem += dist_y12 >= dist_y21 - w12
problem += dist_x13 >= dist_x31 - w13
problem += dist_y13 >= dist_y31 - w13
problem += dist_x23 >= dist_x32 - w23
problem += dist_y23 >= dist_y32 - w23

problem.solve()
x = []
y = []
for v in problem.variables():
    print(v.name, v.varValue)
    if v.name[0] == "x":
        x.append(v.varValue)
    else:
        y.append(v.varValue)

gdf["_x"] = x
gdf["_y"] = y
gdf.geometry = gdf.apply(lambda r: Point(r["_x"], r["_y"]), axis=1)
gdf.geometry = gdf[["geometry", "geom_size"]].apply(create_block, axis=1)
gdf["i"] = [1, 2, 3]
gdf.plot(column="i", legend=True)
plt.show()
