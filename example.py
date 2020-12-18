import sys
import os
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from cartogrammetry.create import Cartogram


def main():
    """US state population example

    """

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
    us_states = us_states.merge(us_inhab, on="NAME")

    print(us_states.info())

    # Create a block style cartogram for inhabitants per state in 2019.
    cg = Cartogram(
        map_type="block",
        gdf=us_states,
        size_column="pop_2019",
        mode=2,
        time_limit=480, # The total amount of seconds the model is allowed to run. Useful for working with mode 3.
        lower_bound_mult=1
    )

    # Plot both the original map and the cartogram side by side.
    f, ax = plt.subplots(1, 2, figsize=(10, 25))
    cg._solver.gdf.plot(
        column="pop_2019",
        ax=ax[0], alpha=0.8)
    cg.gdf.plot(
        column="pop_2019",
        ax=ax[1], alpha=0.8)
    ax[0].axis("off")
    ax[1].axis("off")
    cg._solver.gdf.apply(
        lambda x: ax[0].annotate(
            text=x.STUSPS,
            xy=x.geometry.centroid.coords[0],
            ha="center",
            color="#B2B2B2",
        ),
        axis=1,
    )
    cg.gdf.apply(
        lambda x: ax[1].annotate(
            text=x.STUSPS,
            xy=x.geometry.centroid.coords[0],
            ha="center",
            color="#B2B2B2",
        ),
        axis=1,
    )
    ax[1].set_xlim(-21000000, -5000000)
    plt.show()


if __name__ == "__main__":
    sys.exit(main())
