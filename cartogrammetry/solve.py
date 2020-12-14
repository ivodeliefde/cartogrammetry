#!/usr/bin/env python3

"""Module for solving the location of a blocks or circles in cartogram using linear programming.
"""
from shapely.geometry import Point
import geopandas as gpd
import pulp as p
from tqdm.auto import tqdm


class Solver:
    """ Linear Programming Solver for creating cartograms.

    This class solves the location of geometries in a GeoDataFrame, in order to create a cartogram. The methodology is
    derived from the article 'Computing stable Demers cartograms' by Nickel et al. (2019).

    """

    def __init__(
        self,
        gdf: gpd.GeoDataFrame,
        gap_size: float = 1000.0,
        max_move_dist: float = 10000.0,
    ) -> None:
        """

        :param gdf:
        """

        self.gdf = gdf.copy()
        self.problem = p.LpProblem("Problem", p.LpMinimize)
        self.gap_size = gap_size
        self.max_move_dist = max_move_dist
        self.coord_dict = {}

        self.add_coord_variables()
        self.add_adjacent_variables_constraints()
        self.solve()

    def add_coord_variables(self) -> None:
        """ Method to create LP variables for each coordinate.

        For each geometry in the GeoDataFrame we create 2 variables:
            - X coordinate of centerpoint
            - Y coordinate of centerpoint

        :return:
        """
        for i in range(self.gdf.shape[0]):
            self.coord_dict[i] = {}
            self.coord_dict[i]["x"] = p.LpVariable(f"x_{i}", lowBound=0)
            self.coord_dict[i]["y"] = p.LpVariable(f"y_{i}", lowBound=0)

    def add_adjacent_variables_constraints(self) -> None:
        """ Method to create LP variables for each pair of adjacent geometries.

        For each combination of neighboring geometries we create 2 additional variables:
            - h(i,j) = The non-negative horizontal distance between two squares
            - v(i,j) = The non-negative vertical distance between two squares
        The distance between two squares is calculated as the sum of the width of both
        squares divided by 2.

        The objective is to minimize the sum of the variables h(i,j) en v(i,j). These are
        all the horizontal and vertical distances between geometries.

        The following constraints are added:
            1. The horizontal distance is greater than or equal to the smallest possible width between
            the two centerpoints (p_i and p_j) plus a gap if they are not adjacent.
                x_j - x_i >= (w_j + w_i) / 2 + gap
                where gap = 0 if p_i and p_j are adjacent along the x axis
            2. The vertical distance is greater than or equal to the smallest possible height between
            the two centerpoints (p_i and p_j) plus a gap if they are not adjacent.
                y_j - y_i >= (w_j + w_i) / 2 + gap
                where gap = 0 if i and j are adjacent along the y axis
            3. The horizontal distance between p_i and p_j is greater than or equal to the distance between x_i
            and x_j minus the width of both geometries.
                h(i,j) >= max(x_j - x_i - w(i,j), x_i - x_j - w(i,j))
            4. The vertical distance between p_i and p_j is greater than or equal to the distance between y_i
            and y_j minus the width of both geometries.
                v(i,j) >= max(y_j - y_i - w(i,j), y_i - y_j - w(i,j))

        Additional remarks:
            - For any combination of p_i and p_j either constraint 1 or 2 must apply. Which one is used depends on
            which axis has the greatest range. E.g. if x_j - x_i > y_j - y_i, then constraints 1 applies and p_i and
            p_j are considered horizontally adjacent.

        :return:
        """
        combinations = set()

        for i in tqdm(range(self.gdf.shape[0]), desc="Construction LP"):
            for j in range(self.gdf.shape[0]):
                if j == i:
                    # We don't check the adjacency with itself.
                    continue
                if (i, j) in combinations or (j, i) in combinations:
                    # If the combination already exists we skip it.
                    continue

                # Create a set with all possible combinations to optimize.
                combinations.add((i, j))

        for i, j in combinations:
            # Create the variables for the optimal horizontal en vertical length between p_i en p_j
            h = p.LpVariable(f"h_{i}_{j}", lowBound=0)
            v = p.LpVariable(f"v_{i}_{j}", lowBound=0)

            # Define the objective
            self.problem += h + v

            # Calculate the minimum length between i en j based on their width
            w = self.gdf.at[i, "geom_size"] + self.gdf.at[j, "geom_size"]

            # Add a gap if they are non-adjacent
            if str(j) in self.gdf.iloc[i, :]["_neighbors"].split(","):
                gap = 0
            else:
                gap = self.gap_size

            # Check whether the neighbors should be connected horizontally or vertically.
            geom_dist_x = self.gdf.at[j, "geometry"].x - self.gdf.at[i, "geometry"].x
            geom_dist_y = self.gdf.at[j, "geometry"].y - self.gdf.at[i, "geometry"].y

            # Depending on whether x_i - x_j yields a positive or negative result we decide the order
            # in which they are presented in the constraints.
            if geom_dist_x > 0:
                if abs(geom_dist_x) > abs(geom_dist_y):
                    # Add constraint 1
                    self.problem += (
                        self.coord_dict[j]["x"] - self.coord_dict[i]["x"] >= w + gap
                    )
                # Add constraint 3
                self.problem += (
                    h >= self.coord_dict[j]["x"] - self.coord_dict[i]["x"] - w
                )
            else:
                if abs(geom_dist_x) > abs(geom_dist_y):
                    # Add constraint 1
                    self.problem += (
                        self.coord_dict[i]["x"] - self.coord_dict[j]["x"] >= w + gap
                    )
                # Add constraint 3
                self.problem += (
                    h >= self.coord_dict[i]["x"] - self.coord_dict[j]["x"] - w
                )

            # Depending on whether y_i - y_j yields a positive or negative result we decide the order
            # in which they are presented in the constraints.
            if geom_dist_y > 0:
                if not abs(geom_dist_x) > abs(geom_dist_y):
                    # Add constraint 2
                    self.problem += (
                        self.coord_dict[j]["y"] - self.coord_dict[i]["y"] >= w + gap
                    )
                # Add constraint 4
                self.problem += (
                    v >= self.coord_dict[j]["y"] - self.coord_dict[i]["y"] - w
                )
            else:
                if not abs(geom_dist_x) > abs(geom_dist_y):
                    # Add constraint 2
                    self.problem += (
                        self.coord_dict[i]["y"] - self.coord_dict[j]["y"] >= w + gap
                    )
                # Add constraint 4
                self.problem += (
                    v >= self.coord_dict[i]["y"] - self.coord_dict[j]["y"] - w
                )

    def solve(self):
        """

        :return:
        """

        self.problem.solve()

        for v in self.problem.variables():
            var_name = v.name.split("_")
            if len(var_name) == 2:
                # All variables with two parts in their name are coordinates. The first part says if it's a x or y
                # coordinate, the second part the geometry it corresponds to (row number in Geodataframe).
                c, i = v.name.split("_")
                self.gdf.at[int(i), f"_{c}"] = v.varValue

        self.gdf.geometry = self.gdf.apply(lambda r: Point(r["_x"], r["_y"]), axis=1)