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

    shp_path = os.path.join(os.getcwd(), "data", "NLD_provinces.geojson")
    # shp_path = os.path.join(os.getcwd(), "data", "gemeente_2020_v1_small_subset.shp")
    municipalities = gpd.read_file(shp_path)
    # municipalities = municipalities.loc[municipalities.geometry.centroid.y > municipalities.geometry.centroid.y.mean(), :].loc[municipalities.H2O.str.contains("NEE"), :].reset_index()

    cg = Cartogram(
        map_type="block",
        gdf=municipalities,
        size_column="AANT_INW",
        lower_bound_mult=1,
        upper_bound_mult=1,
    )

    # print(cg._solver.gdf.bounds.describe())
    # print(cg._solver.gdf.geometry.head(12))
    # print(cg.gdf.geometry.head(12))
    # cg._solver.gdf.geometry.to_file("data/provinces_result.shp")

    f, ax = plt.subplots(1, 2, figsize=(15, 15))
    # # cg.gdf["_n_neighbors"] = cg._solver.gdf._n_neighbors
    cg._solver.gdf.plot(
        # column="GM_NAAM",
        column="name",
        ax=ax[0],
        alpha=0.8
    )
    # cg._solver.gdf.centroid.plot(ax=ax[0], alpha=0.8)
    cg.gdf.plot(
        # column="GM_NAAM",
        column="name",
        ax=ax[1],
        alpha=0.8
    )
    # # cg.gdf.centroid.plot(color="orange", ax=ax)
    # # cg._solver.gdf.centroid.plot(color="white", ax=ax)
    ax[0].axis('off')
    ax[1].axis('off')
    plt.show()


if __name__ == "__main__":
    sys.exit(main())
