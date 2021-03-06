import sys
import os
import logging
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

for name in logging.root.manager.loggerDict:
    logging.getLogger(name).setLevel(logging.ERROR)

from cartogrammetry.create import CircleCartogram, SquareCartogram
from cartogrammetry.plot import Map


def main():
    """US state population example

    """

    # Log messages to stdout
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] %(message)s",
        stream=sys.stdout,
    )

    # Load the sample dataset: the US states and their corresponding population number.
    # (data from https://www.census.gov/)
    us_states_path = os.path.join(os.getcwd(), "sample_data", "cb_2018_us_state_5m.shp")
    us_pop_path = os.path.join(os.getcwd(), "sample_data", "nst-est2019-01.xlsx")
    us_states = gpd.read_file(us_states_path)
    us_inhab = pd.read_excel(us_pop_path, skiprows=3, engine="openpyxl").add_prefix(
        "pop_"
    )
    # Tidy up rows and column names
    us_inhab.rename(columns={us_inhab.columns[0]: "NAME"}, inplace=True)
    us_inhab.NAME = us_inhab.NAME.str.replace(".", "")
    # Join population numbers and us state geometries.
    us_states = us_states.merge(us_inhab, on="NAME").reset_index()
    # Inspect the data
    print(us_states.info())

    # Initialize a circle style cartogram for inhabitants per state in 2019.
    circle_cg = CircleCartogram(
        gdf=us_states,
        size_column="pop_2019",
        mode=2,
        time_limit=60,  # The total amount of seconds the model is allowed to run. Useful for working with mode 3.
    )
    square_cg = SquareCartogram(
        gdf=us_states,
        size_column="pop_2019",
        mode=1,
        time_limit=60,  # The total amount of seconds the model is allowed to run. Useful for working with mode 3.
    )
    square2_cg = SquareCartogram(
        gdf=us_states,
        size_column="pop_2019",
        mode=4,
        time_limit=60,  # The total amount of seconds the model is allowed to run. Useful for working with mode 3.
    )

    # Calculate the cartogram geometries.
    circle_cg.calculate()
    square_cg.calculate()
    square2_cg.calculate()

    # Plot both the original map and the cartogram side by side.
    gdfs = [us_states, circle_cg.gdf, square_cg.gdf, square2_cg.gdf]
    m = Map(
        gdfs=gdfs,
        title="Population per US State in 2019",
        column="pop_2019",
        labels="STUSPS",
    )
    m.ax[0][0].set_xlim(-150, -60)
    m.plot()
    plt.show()


if __name__ == "__main__":
    sys.exit(main())
