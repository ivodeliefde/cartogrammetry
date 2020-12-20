#!/usr/bin/env python3

"""Module for creating block or circle style cartograms.
"""

import pandas as pd
from shapely.geometry import Polygon

from .cartogram import Cartogram


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


class CircleCartogram(Cartogram):
    """Extends the Cartogram class to create a cartogram where every geometry is represented as a circle.

    """
    def calculate(self) -> None:
        """

        """
        self.solve_cartogram()
        self._solver.gdf.geometry = self._solver.gdf.geometry.buffer(
            self._solver.gdf["geom_size"]
        )


class SquareCartogram(Cartogram):
    """Extends the Cartogram class to create a cartogram where every geometry is represented as a square.

    """
    def calculate(self) -> None:
        """

        """
        self.solve_cartogram()
        self._solver.gdf.geometry = self._solver.gdf[["geometry", "geom_size"]].apply(
            create_block, axis=1
        )
