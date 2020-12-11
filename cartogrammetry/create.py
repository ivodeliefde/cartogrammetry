#!/usr/bin/env python3

"""Module for creating a blockstyle cartogram.
"""

import os
from pathlib import Path
import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
import sys


def create_block(row) -> Polygon:
    """

    :return:
    """
    c = row["geometry"].centroid
    ll = (c.x - row["size"], c.y - row["size"])
    ul = (c.x - row["size"], c.y + row["size"])
    ur = (c.x + row["size"], c.y + row["size"])
    lr = (c.x + row["size"], c.y - row["size"])

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
        self.gdf = gdf
        self.size_column = size_column
        self.lower_bound_mult = lower_bound_mult
        self.upper_bound_mult = upper_bound_mult
        # Take a lower and upper bound based on the areas and apply the multiplier
        self._lower_bound = (
            np.sqrt(self.gdf.geometry.area.quantile(0.01)) * self.lower_bound_mult
        )
        self._upper_bound = (
            np.sqrt(self.gdf.geometry.area.quantile(0.99)) * self.upper_bound_mult
        )
        # Create a scaler and calculate a size value based on the size_column
        self._scaler = MinMaxScaler(
            feature_range=(self._lower_bound, self._upper_bound)
        )
        if self.size_column in self.gdf.columns:
            self.gdf["size"] = self._scaler.fit_transform(
                self.gdf.loc[:, self.size_column].values.reshape(-1, 1)
            )
        else:
            # TODO! check if geometry type is polygon
            # If the size_column is None, calculate size from area
            self.gdf["size"] = self._scaler.fit_transform(
                self.gdf.geometry.area.values.reshape(-1, 1)
            )

        if self.map_type.lower() == "circle":
            self.create_circle_map()
        elif self.map_type.lower() == "block":
            self.create_block_map()

    def create_circle_map(self) -> None:
        """

        :return:
        """

        self.gdf.geometry = self.gdf.geometry.centroid.buffer(self.gdf["size"])

    def create_block_map(self) -> None:
        """

        :return:
        """

        self.gdf.geometry = self.gdf[["geometry", "size"]].apply(
            create_block, axis=1
        )


def main():
    """

    :return:
    """

    shp_path = os.path.join(Path(os.getcwd()).parent, "data", "gemeente_2020_v1.shp")
    municipalities = gpd.read_file(shp_path)
    municipalities = municipalities.loc[municipalities.H2O.str.contains("NEE"), :]
    print(municipalities.info())
    cg = Cartogram(map_type="circle", gdf=municipalities, size_column="AANT_INW")
    print(cg.gdf.geometry.head())
    cg.gdf.plot(column="size")
    plt.axis("off")
    plt.show()
    # print(cg._lower_bound)
    # print(cg._upper_bound)


if __name__ == "__main__":
    sys.exit(main())
