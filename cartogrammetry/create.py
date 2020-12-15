#!/usr/bin/env python3

"""Module for creating block or circle style cartograms.
"""
import sys
import os
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from tqdm.auto import tqdm

from .solve import Solver


def create_block(row: pd.Series) -> Polygon:
    """

    :return:

    """
    c = row["geometry"]
    ll = (c.x - row["geom_size"], c.y - row["geom_size"])
    ul = (c.x - row["geom_size"], c.y + row["geom_size"])
    ur = (c.x + row["geom_size"], c.y + row["geom_size"])
    lr = (c.x + row["geom_size"], c.y - row["geom_size"])
    return Polygon([ll, ul, ur, lr])


class Cartogram:
    """ Class to create a cartogram from a geopandas geodataframe.

    """

    def __init__(
        self,
        gdf: gpd.GeoDataFrame,
        map_type: str,
        size_column: str = None,
        lower_bound_mult: float = 1.0,
        upper_bound_mult: float = 1.0,
    ) -> None:
        """

        :param map_type:

        """
        self.map_type = map_type
        self.gdf = gdf.to_crs(28992)
        self.size_column = size_column
        self.lower_bound_mult = lower_bound_mult
        self.upper_bound_mult = upper_bound_mult
        self._lower_bound = 0
        self._upper_bound = 1
        self._scaler = None
        self._solver = None

        # Create the cartogram.
        self.create_cartogram()

    def create_cartogram(self) -> None:
        """

        :return:

        """

        self.calculate_size()
        self.calculate_neighbors()

        gdf = self.gdf.copy()
        gdf.geometry = self.gdf.geometry.centroid

        self._solver = Solver(gdf=gdf)

        if self.map_type.lower() == "circle":
            self.create_circle_map()
        elif self.map_type.lower() == "block":
            self.create_block_map()

    def calculate_size(self) -> None:
        """Calculate the size for each geometry based on a numeric column.

        """

        # Take a lower and upper bound based on the areas and apply the multiplier
        self._lower_bound = (
            np.sqrt(self.gdf.geometry.area.quantile(0.01)) * self.lower_bound_mult
        ) / 2
        self._upper_bound = (
            np.sqrt(self.gdf.geometry.area.quantile(0.99)) * self.upper_bound_mult
        ) / 2

        # Create a scaler and calculate a size value based on the size_column
        print(self._lower_bound, self._upper_bound)
        self._scaler = MinMaxScaler(
            feature_range=(self._lower_bound, self._upper_bound)
        )

        if self.size_column in self.gdf.columns:
            self.gdf["geom_size"] = self._scaler.fit_transform(
                self.gdf.loc[:, self.size_column].values.reshape(-1, 1)
            )
        else:
            # TODO! check if geometry type is polygon
            # If the size_column is None, calculate size from area
            self.gdf["geom_size"] = self._scaler.fit_transform(
                self.gdf.geometry.area.values.reshape(-1, 1)
            )

    def calculate_neighbors(self) -> None:
        """Calculate which geometries are adjacent to each other.

        """
        self.gdf["_neighbors"] = ""
        for index, row in tqdm(self.gdf.iterrows(), desc="Calculating neighbors"):
            neighbors = self.gdf[
                ~self.gdf.geometry.disjoint(row["geometry"])
            ].index.tolist()
            self.gdf.at[index, "_neighbors"] = ",".join([str(n) for n in neighbors if n != index])
            # for n in self.gdf.at[index, "_neighbors"].split(","):
            #     print(index, self.gdf.at[index, "name"], n, self.gdf.at[int(n), "name"])
            self.gdf.at[index, "_n_neighbors"] = len(neighbors)

        self.gdf["_n_neighbors"].fillna(0, inplace=True)

    def create_circle_map(self) -> None:
        """Create a cartogram where every geometry is represented as a circle.

        """

        self._solver.gdf.geometry = self._solver.gdf.geometry.buffer(
            self._solver.gdf["geom_size"]
        )

    def create_block_map(self) -> None:
        """Create a cartogram where every geometry is represented as a block.

        """
        self._solver.gdf.geometry = self._solver.gdf[["geometry", "geom_size"]].apply(
            create_block, axis=1
        )


def main():
    """Main function for testing purposes

    """

    shp_path = os.path.join(Path(os.getcwd()).parent, "data", "gemeente_2020_v1.shp")
    municipalities = gpd.read_file(shp_path)
    municipalities = municipalities.loc[municipalities.H2O.str.contains("NEE"), :]
    # municipalities.plot(column="GM_NAAM")
    # plt.axis("off")
    # plt.show()
    print(municipalities.info())
    cg = Cartogram(
        map_type="block",
        gdf=municipalities,
        size_column="AANT_INW",
        lower_bound_mult=1,
        upper_bound_mult=1,
    )
    cg.gdf.plot(column="geom_size")
    plt.axis("off")
    plt.show()


if __name__ == "__main__":
    sys.exit(main())
