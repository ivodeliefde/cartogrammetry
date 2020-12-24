import sys
import os
import logging
import geopandas as gpd
import matplotlib.pyplot as plt

from cartogrammetry.create import CircleCartogram, SquareCartogram
from cartogrammetry.plot import Map


def main():
    """NL municipalities example

    """

    # Log messages to stdout
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        stream=sys.stdout,
    )

    # Load the sample dataset: the NL municipalities.
    # (data from https://www.cbs.nl/nl-nl/dossier/nederland-regionaal/geografische-data)
    nl_muni_path = os.path.join(
        os.getcwd(), "sample_data", "gemeente_2020_v1_subset.shp"
    )
    nl_muni = gpd.read_file(nl_muni_path)

    # Inspect the data
    print(nl_muni.info(""))

    # Initialize a circle style cartogram for inhabitants per state in 2019.
    circle_cg = CircleCartogram(
        gdf=nl_muni,
        mode=2,
        time_limit=60,  # The total amount of seconds the model is allowed to run. Useful for working with mode 3.
    )
    square_cg = SquareCartogram(
        gdf=nl_muni,
        mode=1,
        time_limit=60,  # The total amount of seconds the model is allowed to run. Useful for working with mode 3.
    )
    square2_cg = SquareCartogram(
        gdf=nl_muni,
        mode=4,
        time_limit=60,  # The total amount of seconds the model is allowed to run. Useful for working with mode 3.
    )

    # Calculate the cartogram geometries.
    circle_cg.calculate()
    square_cg.calculate()
    square2_cg.calculate()

    # Plot both the original map and the cartogram side by side.
    gdfs = [nl_muni, circle_cg.gdf, square_cg.gdf, square2_cg.gdf]
    m = Map(
        gdfs=gdfs,
        title="Municipalities in the Netherlands",
        column="GM_NAAM",
        labels="GM_NAAM",
    )
    m.plot()
    plt.show()


if __name__ == "__main__":
    sys.exit(main())
