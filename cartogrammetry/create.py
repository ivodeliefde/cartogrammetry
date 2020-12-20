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

from ._cartogram import Cartogram


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
    def calculate(self) -> None:
        """Create a cartogram where every geometry is represented as a circle.

        """
        self.solve_cartogram()
        self._solver.gdf.geometry = self._solver.gdf.geometry.buffer(
            self._solver.gdf["geom_size"]
        )

class BlockCartogram(Cartogram):
    def calculate(self) -> None:
        """Create a cartogram where every geometry is represented as a block.

        """
        self.solve_cartogram()
        self._solver.gdf.geometry = self._solver.gdf[["geometry", "geom_size"]].apply(
            create_block, axis=1
        )


