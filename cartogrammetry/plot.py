#!/usr/bin/env python3

"""Module for creating block or circle style cartograms.
"""
import os
import logging
import math
import tempfile
import requests
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.patheffects as pe


class Map:
    """Class to plot a map with predefined styling from a geopandas geodataframe.

    :param gdfs: List containing one or more GeoDataFrames with polygon features to plot.
    :param title: Title to add to the plot.
    :param column: Name of dataframe column to use for determining fill colours per feature.
    :param labels: Name of dataframe column to use for adding annotations per feature.
    """

    def __init__(
        self,
        gdfs: list,
        title: str = "Cartogram",
        column: str = None,
        labels: str = None,
    ) -> None:
        """

        :param gdfs:
        :param title:
        :param column:
        :param labels:
        """

        self.gdfs = gdfs
        self.title = title
        self.column = column
        self.labels = labels

        # Plot configuration
        self.ncols = 2
        self.nrows = math.ceil(len(gdfs) / self.ncols)
        self.figsize = (20, 25)
        self.alpha = 0.8

        # Title configuration
        self.preferred_title_font_family = "slabo27px"
        self.preferred_title_font = "Slabo27px-Regular.ttf"
        self.title_font = self.preferred_title_font
        self.title_font_path = os.path.join(tempfile.gettempdir(), self.preferred_title_font)
        self.title_font_url = f"https://github.com/google/fonts/blob/master/ofl/{self.preferred_title_font_family}/{self.preferred_title_font}?raw=true"
        self.title_font_size = 35
        self.title_font_properties = fm.FontProperties(fname=self.title_font_path)

        # Label configuration
        self.label_colour = "#EBEBEB"
        self.preferred_label_font_family = "quicksand"
        self.preferred_label_font = "Quicksand-Regular.ttf"
        self.label_font = self.preferred_label_font
        self.label_font_path = os.path.join(tempfile.gettempdir(), self.preferred_label_font)
        self.label_font_url = f"https://github.com/google/fonts/blob/master/ofl/{self.preferred_label_font_family}/static/{self.preferred_label_font}?raw=true"
        self.label_font_size = 6.5
        self.label_font_properties = fm.FontProperties(fname=self.label_font_path)

        # Initialize plot
        self.f, self.ax = plt.subplots(self.nrows, self.ncols, figsize=self.figsize)

    def get_font(self):
        """Get a nice title font.
        """
        logging.debug(f"Downloading font: {self.title_font_url}")
        r = requests.get(self.title_font_url, allow_redirects=True)
        font = self.title_font_url.split("/")[-1].split(".ttf")[0] + ".ttf"
        logging.debug(f"Download {font} to temp directory: {tempfile.gettempdir()}")
        self.title_font_path = os.path.join(tempfile.gettempdir(), font)
        with open(self.title_font_path, 'wb') as f:
            f.write(r.content)

        self.title_font = font

    def check_font(self, title=True):
        """Add specific font to the map.
        """
        if title:
            try:
                if not os.path.exists(self.title_font_path):
                    self.get_font()
            except Exception as e:
                logging.error(e)
                # If we can't download our preferred font we take an already available sans-serif font.
                self.title_font_path = fm.findfont(fm.FontProperties(family=['serif']))
                self.title_font = os.path.basename(self.title_font_path)
                self.title_font_properties = fm.FontProperties(fname=self.title_font_path)

            logging.debug(f"Using title font: {self.title_font}")
        else:
            try:
                if not os.path.exists(self.label_font_path):
                    self.get_font()
            except Exception as e:
                logging.error(e)
                # If we can't download our preferred font we take an already available sans-serif font.
                self.label_font_path = fm.findfont(fm.FontProperties(family=['sans-serif']))
                self.label_font = os.path.basename(self.label_font_path)
                self.label_font_properties = fm.FontProperties(fname=self.label_font_path)

            logging.debug(f"Using label font: {self.label_font}")

    def create_plot(self, **kwargs):
        """Fill the plot with data from the GeoDataframe(s).
        """
        self.check_font()
        self.f.suptitle(self.title, fontsize=self.title_font_size, fontproperties=self.title_font_properties)

        c = 0
        for i in range(self.nrows):
            for j in range(self.ncols):
                # Plot Geodataframe with column for fill if present
                if self.column:
                    self.gdfs[c].plot(
                        column=self.column, ax=self.ax[i][j], alpha=self.alpha, **kwargs
                    )
                else:
                    self.gdfs[c].plot(ax=self.ax[i][j], alpha=self.alpha, **kwargs)

                # Remove the axes
                self.ax[i][j].axis("off")

                c += 1

    def add_labels(self):
        """Add labels to the plot.
        """
        self.check_font(title=False)
        c = 0
        for i in range(self.nrows):
            for j in range(self.ncols):
                if not "geom_size" in self.gdfs[c].columns:
                    try:
                        self.gdfs[c]["geom_size"] = self.gdfs[c].geometry.area
                    except:
                        pass
                self.gdfs[c].apply(
                    lambda x: self.ax[i][j].annotate(
                        text=x[self.labels],
                        xy=x.geometry.centroid.coords[0],
                        ha="center",
                        va="center",
                        color=self.label_colour,
                        fontproperties=self.label_font_properties,
                        fontsize=self.label_font_size + 5 * x["geom_size"] / self.gdfs[c]["geom_size"].median(),
                        path_effects=[pe.withStroke(linewidth=.5, foreground="#939393")]
                    ),
                    axis=1,
                )
                c += 1

    def plot(self, **kwargs):
        """Return the figure to display it.
        """
        self.create_plot(**kwargs)
        if self.labels:
            self.add_labels()
        return self.f
