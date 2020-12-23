#!/usr/bin/env python3

"""Module for creating block or circle style cartograms.
"""
import logging
import numpy as np
import geopandas as gpd
from shapely.geometry import box, Point, shape, mapping
from shapely.affinity import translate
from tqdm.auto import tqdm

from .solve import Solver


class Cartogram:
    """Class to create a cartogram from a geopandas geodataframe.

    :param gdf: GeoDataFrame with polygon features to create a cartogram from.
    :param size_column: Numeric dataframe column to use for scaling the cartogram.
    :param time_limit: The maximum time allowed to find a solution. If exceeded the best solution at that time is returned.
    :param mode: Mode 1, 2 or 3. Corresponds to 1: quick method but suboptimal, 2: quick, optimal, but with overlapping geometries, 3: optimal, no overlapping geometries, but slow.
    """

    def __init__(
        self,
        gdf: gpd.GeoDataFrame,
        size_column: str = None,
        mode: int = 1,
        time_limit: int = 300,
    ) -> None:
        """

        :param gdf:
        :param map_type:
        :param size_column:
        :param mode:
        :param time_limit:
        """
        self.crs = gdf.crs
        self.gdf_original = gdf.copy()
        self.gdf = self.gdf_original.to_crs(3857)
        self.size_column = size_column
        self.mode = mode
        self.time_limit = time_limit
        logging.debug("Initialized Cartogram")

    def offset(self) -> gpd.geoseries:
        """Offset geometries (required for ILP)

        """
        logging.debug("Method offset called")
        _x_offset = -1.0 * self.gdf.bounds.minx.min()
        _y_offset = -1.0 * self.gdf.bounds.miny.min()

        return self.gdf.geometry.apply(
            lambda x: translate(x, xoff=_x_offset, yoff=_y_offset)
        )

    # def undo_offset(self, geom: gpd.GeoSeries) -> gpd.GeoSeries:
    #     """Return geometries to their original location.
    #
    #     """
    #     orig_center = box(*self.gdf.geometry.total_bounds).centroid
    #     new_center = box(*geom.total_bounds).centroid
    #
    #     _x_offset = orig_center.x - new_center.x
    #     _y_offset = orig_center.y - new_center.y
    #
    #     return geom.apply(
    #         lambda x: translate(x, xoff=_x_offset, yoff=_y_offset)
    #     )
    #
    # def undo_scaling(self, geom: gpd.GeoSeries) -> gpd.geoseries:
    #     """Scale geometries to their origin size.
    #
    #     """
    #     pass

    def solve_cartogram(self) -> None:
        """Method for running all the steps required to make a cartogram.

        """
        logging.debug("Method solve_cartogram called")
        self.calculate_size()
        self.calculate_neighbors()

        gdf = self.gdf.copy()
        gdf.geometry = self.offset().centroid

        solver = Solver(gdf=gdf, mode=self.mode, time_limit=self.time_limit)
        solver.run()

        self.gdf.geometry = solver.gdf.apply(lambda r: Point(r["_x"], r["_y"]), axis=1)

        # self.gdf.geometry = self.undo_offset(gdf.geometry)


    def calculate_size(self) -> None:
        """Calculate the size for each geometry based on a numeric column.

        """
        logging.debug("Method calculate_size called")
        # Take a lower and upper bound based on the areas and apply the multiplier
        lower_bound = np.sqrt(self.gdf.geometry.area.median())/2

        if self.size_column in self.gdf.columns:
            self.gdf["geom_size"] = np.sqrt(self.gdf.loc[:, self.size_column]) / 2 / self.gdf.loc[:, self.size_column].median() * lower_bound
        else:
            # TODO! check if geometry type is polygon
            # If the size_column is None, calculate size from area
            self.gdf["geom_size"] = np.sqrt(self.gdf.geometry.area) / 2

    def calculate_neighbors(self) -> None:
        """Calculate which geometries are adjacent to each other.

        """
        logging.debug("Method calculate_neighbors called")
        self.gdf["_neighbors"] = ""
        for index, row in tqdm(self.gdf.iterrows(), desc="Calculating neighbors"):
            # Find all neighbors that are not disjoint from each geometry
            neighbors = self.gdf[
                ~self.gdf.geometry.disjoint(row["geometry"])
            ].index.tolist()

            perimeter = row["geometry"].length
            neighbors_perimeter = []
            for n in neighbors:
                if n == index:
                    continue
                # Calculate for each neighbor the percentage of overlap with respect to the full perimeter.
                length = shape(mapping(row["geometry"].intersection(self.gdf.at[n, "geometry"]))).length
                neighbors_perimeter.append(round(length / perimeter, 2))


            # Write neighbor numbers to gdf
            self.gdf.at[index, "_neighbors"] = ",".join(
                [str(n) for n in neighbors if n != index]
            )
            self.gdf.at[index, "_neighbors_perimeter"] = ",".join(
                [str(p) for p in neighbors_perimeter]
            )
            self.gdf.at[index, "_n_neighbors"] = len(neighbors) - 1 # subtract 1 because a feature is not counted as its own neighbor.
            # logging.debug(f"{self.gdf.at[index, 'NAME']} has {self.gdf.at[index, '_n_neighbors']} neighbors: {','.join([self.gdf.at[n, 'NAME'] for n in neighbors if n != index])} with perimeters: {self.gdf.at[index, '_neighbors_perimeter']}")

        self.gdf["_n_neighbors"].fillna(0, inplace=True)
