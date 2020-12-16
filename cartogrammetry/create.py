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
from shapely.affinity import translate
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from tqdm.auto import tqdm

from .solve import Solver


def create_block(row: pd.Series) -> Polygon:
    """Function to create a polygon feature from a centerpoint and width.

    :return: shapely.geometry.Polygon
    """
    c = row["geometry"]
    ll = (c.x - row["geom_size"], c.y - row["geom_size"])
    ul = (c.x - row["geom_size"], c.y + row["geom_size"])
    ur = (c.x + row["geom_size"], c.y + row["geom_size"])
    lr = (c.x + row["geom_size"], c.y - row["geom_size"])
    return Polygon([ll, ul, ur, lr])


class Cartogram:
    """Class to create a cartogram from a geopandas geodataframe.

    :param gdf: GeoDataFrame with polygon features to create a cartogram from.
    :param map_type: Type of cartogram, either 'block' or 'circle'.
    :param size_column: Numeric dataframe column to use for scaling the cartogram.
    :param lower_bound_mult: multiplier to affect the lower bound for scaling the cartogram.
    :param upper_bound_mult: multiplier to affect the upper bound for scaling the cartogram.
    :param normalize: Boolean value, if true the values in the size_column are divided by the area of the input geometry.
    :param time_limit: The maximum time allowed to find a solution. If exceeded the best solution at that time is returned.
    :param mode: Mode 1, 2 or 3. Corresponds to 1: quick method but suboptimal, 2: quick, optimal, but with overlapping geometries, 3: optimal, no overlapping geometries, but slow.
    """

    def __init__(
        self,
        gdf: gpd.GeoDataFrame,
        map_type: str,
        size_column: str = None,
        lower_bound_mult: float = 1.0,
        upper_bound_mult: float = 1.0,
        normalize: bool = False,
        mode: int = 1,
        time_limit: int = 300,
    ) -> None:
        """Init method.

        :return: None
        """
        self.map_type = map_type
        self.gdf = gdf.to_crs(3857)
        self.size_column = size_column
        self.lower_bound_mult = lower_bound_mult
        self.upper_bound_mult = upper_bound_mult
        self.normalize = normalize
        self.mode = mode
        self.time_limit = time_limit

        self._lower_bound = 0
        self._upper_bound = 1
        self._scaler = None
        self._solver = None
        self._x_offset = -1.0 * self.gdf.bounds.minx.min()
        self._y_offset = -1.0 * self.gdf.bounds.miny.min()

        if self.gdf.bounds.minx.min() > 0:
            self._x_offset += 10000
        else:
            self._x_offset -= 10000
        if self.gdf.bounds.miny.min() > 0:
            self._y_offset += 10000
        else:
            self._y_offset -= 10000

        # Apply the offset.
        self.gdf.geometry = self.gdf.geometry.apply(
            lambda x: translate(x, xoff=self._x_offset, yoff=self._y_offset)
        )

        # Create the cartogram.
        self.create_cartogram()

        # Undo the offset
        self.gdf.geometry = self.gdf.geometry.apply(
            lambda x: translate(x, xoff=-1 * self._x_offset, yoff=-1 * self._y_offset)
        )

        self._solver.gdf.geometry = self._solver.gdf.geometry.apply(
            lambda x: translate(x, xoff=-1 * self._x_offset, yoff=-1 * self._y_offset)
        )

    def create_cartogram(self) -> None:
        """Method for running all the steps required to make a cartogram.

        :return: None
        """

        self.calculate_size()
        self.calculate_neighbors()

        gdf = self.gdf.copy()
        gdf.geometry = self.gdf.geometry.centroid

        self._solver = Solver(gdf=gdf, mode=self.mode, time_limit=self.time_limit)
        self._solver.run()

        if self.map_type.lower() == "circle":
            self.create_circle_map()
        elif self.map_type.lower() == "block":
            self.create_block_map()

    def calculate_size(self) -> None:
        """Calculate the size for each geometry based on a numeric column.

        """

        # Take a lower and upper bound based on the areas and apply the multiplier
        self._lower_bound = 0.1
        self._upper_bound = 1

        # Create a scaler and calculate a size value based on the size_column
        self._scaler = MinMaxScaler(
            feature_range=(self._lower_bound, self._upper_bound)
        )

        if self.size_column in self.gdf.columns:
            self.gdf["geom_size"] = (
                self._scaler.fit_transform(
                    self.gdf.loc[:, self.size_column].values.reshape(-1, 1)
                )
                * np.sqrt(self.gdf.geometry.area.max())
                / 2
            )
        else:
            # TODO! check if geometry type is polygon
            # If the size_column is None, calculate size from area
            self.gdf["geom_size"] = np.sqrt(self.gdf.geometry.area) / 2

        if self.normalize:
            self.gdf.geom_size = self.gdf.geom_size * (
                self.gdf.geometry.area / self.gdf.geometry.area.max()
            )

    def calculate_neighbors(self) -> None:
        """Calculate which geometries are adjacent to each other.

        """
        self.gdf["_neighbors"] = ""
        for index, row in tqdm(self.gdf.iterrows(), desc="Calculating neighbors"):
            neighbors = self.gdf[
                ~self.gdf.geometry.disjoint(row["geometry"])
            ].index.tolist()
            self.gdf.at[index, "_neighbors"] = ",".join(
                [str(n) for n in neighbors if n != index]
            )
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
