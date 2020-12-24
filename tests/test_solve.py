import unittest
import pandas as pd
import geopandas as gpd
from shapely import wkt

# import matplotlib.pyplot as plt

from cartogrammetry.solve import Solver

# Create some test data.
FEATURES = {
    "geom_size": [1, 5, 2],
    "_neighbors": ["1", "0,2", "1"],
    "geometry": [
        "POLYGON ((10 10, 10 20, 20 10, 10 10))",
        "POLYGON ((10 20, 20 20, 20 10, 10 20))",
        "POLYGON ((20 10, 20 20, 30 10, 20 10))",
    ],
}


class SolverTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        df = pd.DataFrame(FEATURES)
        df["geometry"] = df.geometry.apply(wkt.loads).apply(lambda x: x.centroid)
        cls.gdf = gpd.GeoDataFrame(df, geometry="geometry")

    def test_solver_generates_output(self):
        s = Solver(gdf=self.gdf, gap_size=0.5, time_limit=300, mode=1, solve_msg=False)
        s.run()
        self.assertTrue("_x" in s.gdf.columns)
        self.assertTrue("_y" in s.gdf.columns)

    def test_mode_one_returns_correct_distances(self):
        s = Solver(gdf=self.gdf, gap_size=0.5, time_limit=300, mode=1, solve_msg=False)
        s.run()
        for i in range(s.gdf.shape[0]):
            for n in s.gdf.at[i, "_neighbors"].split(","):
                dx = abs(s.gdf.at[i, "geometry"].x - s.gdf.at[int(n), "geometry"].x)
                dy = abs(s.gdf.at[i, "geometry"].y - s.gdf.at[int(n), "geometry"].y)
                self.assertTrue(
                    s.gdf.loc[[i, int(n)], "geom_size"].sum() == dx
                    or s.gdf.loc[[i, int(n)], "geom_size"].sum() == dy,
                    f"Solver returns incorrect distance between features {i} and {n}.",
                )

    def test_mode_two_returns_correct_distances(self):
        s = Solver(gdf=self.gdf, gap_size=0.5, time_limit=300, mode=2, solve_msg=False)
        s.run()
        for i in range(s.gdf.shape[0]):
            for n in s.gdf.at[i, "_neighbors"].split(","):
                dist = s.gdf.at[i, "geometry"].distance(s.gdf.at[int(n), "geometry"])
                size = s.gdf.loc[[i, int(n)], "geom_size"].sum()
                self.assertGreaterEqual(
                    dist,
                    size,
                    f"Solver returns incorrect distance between features {i} and {n}.",
                )

    def test_mode_three_returns_correct_distances(self):
        s = Solver(gdf=self.gdf, gap_size=0.5, time_limit=300, mode=3, solve_msg=False)
        s.run()
        for i in range(s.gdf.shape[0]):
            for n in s.gdf.at[i, "_neighbors"].split(","):
                dx = abs(s.gdf.at[i, "geometry"].x - s.gdf.at[int(n), "geometry"].x)
                dy = abs(s.gdf.at[i, "geometry"].y - s.gdf.at[int(n), "geometry"].y)
                self.assertTrue(
                    s.gdf.loc[[i, int(n)], "geom_size"].sum() == dx
                    or s.gdf.loc[[i, int(n)], "geom_size"].sum() == dy,
                    f"Solver returns incorrect distance between features {i} and {n}.",
                )


if __name__ == "__main__":
    unittest.main()
