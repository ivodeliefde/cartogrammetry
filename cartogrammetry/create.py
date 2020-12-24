#!/usr/bin/env python3

"""Module for creating square or circle style cartograms.
"""
import logging
import pandas as pd
from shapely.geometry import Polygon

from .cartogram import Cartogram


def create_square(row: pd.Series) -> Polygon:
    """Function to create a square polygon feature from a centerpoint and width.

    :return: shapely.geometry.Polygon
    """
    logging.debug("Function create_square calculate")
    c = row["geometry"]
    ll = (c.x - row["geom_size"], c.y - row["geom_size"])
    ul = (c.x - row["geom_size"], c.y + row["geom_size"])
    ur = (c.x + row["geom_size"], c.y + row["geom_size"])
    lr = (c.x + row["geom_size"], c.y - row["geom_size"])
    return Polygon([ll, ul, ur, lr])


class CircleCartogram(Cartogram):
    """Extends the Cartogram class to create a cartogram where every geometry is represented as a circle.

    """

    def calculate(self) -> None:
        """

        """
        logging.debug("Method calculate from CircleCartogram called")
        self.solve_cartogram()
        self.gdf.geometry = self.gdf.geometry.buffer(self.gdf["geom_size"])
        self.gdf.to_crs(self.crs, inplace=True)
        logging.debug("Finished method calculate")


class SquareCartogram(Cartogram):
    """Extends the Cartogram class to create a cartogram where every geometry is represented as a square.

    """

    def calculate(self) -> None:
        """

        """
        logging.debug("Method calculate from SquareCartogram called")
        self.solve_cartogram()
        self.gdf.geometry = self.gdf[["geometry", "geom_size"]].apply(
            create_square, axis=1
        )
        self.gdf.to_crs(self.crs, inplace=True)
        logging.debug("Finished method calculate")
