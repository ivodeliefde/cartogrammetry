import sys
import os
import geopandas as gpd
import matplotlib.pyplot as plt
from cartogrammetry.create import Cartogram
from cartogrammetry.solve import Solver


def main():
    """

    :return:
    """

    shp_path = os.path.join(os.getcwd(), "data", "PROV1980.SHP")
    municipalities = gpd.read_file(shp_path)
    municipalities = municipalities
    print(municipalities.head())
    cg = Cartogram(
        map_type="circle",
        gdf=municipalities,
        size_column="AANT_INW",
        lower_bound_mult=1,
        upper_bound_mult=1,
    )

    f, ax = plt.subplots(1, 2, figsize=(50, 50))
    # cg.gdf["_n_neighbors"] = cg._solver.gdf._n_neighbors
    cg._solver.gdf.plot(column="GM_NAAM", ax=ax[0], alpha=0.8)
    cg.gdf.plot(column="GM_NAAM", ax=ax[1], alpha=0.8)
    # cg.gdf.centroid.plot(color="orange", ax=ax)
    # cg._solver.gdf.centroid.plot(color="white", ax=ax)
    # plt.axis("off")
    plt.show()


if __name__ == "__main__":
    sys.exit(main())
