#!/usr/bin/env python3

"""Module for solving the location of a blocks or circles in cartogram using linear programming.

"""
from shapely.geometry import Point
import geopandas as gpd
import pulp as p
from tqdm.auto import tqdm


class Solver:
    """Linear Programming Solver for creating cartograms.

    This class solves the location of geometries in a GeoDataFrame, in order to create a cartogram. The methodology is derived from the article 'Computing stable Demers cartograms' by Nickel et al. (2019).

    :param gdf: Geodataframe to optimize.
    :param gap_size: The minimum length in meters between non-adjacent geometries.
    :param time_limit: The maximum time allowed to find a solution. If exceeded the best solution at that time is returned.
    :param mode: Mode 1, 2 or 3. Corresponds to 1: quick method but suboptimal, 2: quick, optimal, but with overlapping geometries, 3: optimal, no overlapping geometries, but slow.
    :param solve_msg: Boolean value for displaying solver logging.
    """

    def __init__(
        self,
        gdf: gpd.GeoDataFrame,
        gap_size: float = 0.5,
        time_limit: int = 300,
        mode: int = 1,
        solve_msg: bool = True,
    ) -> None:
        """

        :param gdf:
        :param gap_size:
        :param time_limit:
        :param mode:
        :param solve_msg:
        """

        self.gdf = gdf.copy()
        self.problem = p.LpProblem("Problem", p.LpMinimize)
        self.gap_size = gap_size
        self.time_limit = time_limit
        self.mode = mode
        self.solve_msg = solve_msg

        self._coord_dict = {}
        self._distances = []

    def run(self) -> None:
        """Methode to run all required steps to find the optimal centerpoint locations for each geometry in a cartogram.

        """
        self.add_coord_variables()
        self.add_adjacent_variables_constraints()
        self.solve()

    def add_coord_variables(self) -> None:
        """Method to create LP variables for each coordinate.

For each geometry in the GeoDataFrame we create 2 variables:

    - X coordinate of centerpoint
    - Y coordinate of centerpoint

        """
        for i in range(self.gdf.shape[0]):
            self._coord_dict[i] = {}
            self._coord_dict[i]["x"] = p.LpVariable(f"x_{i}")
            self._coord_dict[i]["y"] = p.LpVariable(f"y_{i}")

            self.problem += self._coord_dict[i]["x"] >= 0
            self.problem += self._coord_dict[i]["y"] >= 0

        # diff_x = self._coord_dict[i]["x"] - self.gdf.at[i, "geometry"].x
        # diff_y = self._coord_dict[i]["y"] - self.gdf.at[i, "geometry"].y
        # self._distances.append(diff_x)
        # self._distances.append(diff_y)

    def add_adjacent_variables_constraints(self, s_dict=None) -> None:
        """Method to create LP variables for each pair of adjacent geometries.

For each combination of neighboring geometries we create 2 additional variables:

    - h(i,j) = The non-negative horizontal distance between two squares
    - v(i,j) = The non-negative vertical distance between two squares

The distance between two squares is calculated as the sum of the width of both squares divided by 2.

The objective is to minimize the sum of the variables h(i,j) en v(i,j). These are all the horizontal and vertical distances between geometries.

The following constraints are added:

    1. The horizontal distance is greater than or equal to the smallest possible width between the two centerpoints (p_i and p_j) plus a gap if they are not adjacent.

        x_j - x_i >= (w_j + w_i) / 2 + gap

        where gap = 0 if p_i and p_j are adjacent along the x axis
    2. The vertical distance is greater than or equal to the smallest possible height between the two centerpoints (p_i and p_j) plus a gap if they are not adjacent.

        y_j - y_i >= (w_j + w_i) / 2 + gap

        where gap = 0 if i and j are adjacent along the y axis
    3. The horizontal distance between p_i and p_j is greater than or equal to the distance between x_i and x_j minus the width of both geometries.

        h(i,j) >= max(x_j - x_i - w(i,j), x_i - x_j - w(i,j))

    4. The vertical distance between p_i and p_j is greater than or equal to the distance between y_i and y_j minus the width of both geometries.

        v(i,j) >= max(y_j - y_i - w(i,j), y_i - y_j - w(i,j))

Additional remarks:

    - For any combination of p_i and p_j either constraint 1 or 2 must apply. Which one is used depends on which axis has the greatest range. E.g. if x_j - x_i > y_j - y_i, then constraints 1 applies and p_i and p_j are considered horizontally adjacent.

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
            h = p.LpVariable(f"h_{i}_{j}")
            v = p.LpVariable(f"v_{i}_{j}")
            self.problem += h >= 0
            self.problem += v >= 0

            # Calculate the minimum length between i en j based on their width
            w = self.gdf.at[i, "geom_size"] + self.gdf.at[j, "geom_size"]

            # Add a gap if they are non-adjacent
            if str(j) in self.gdf.iloc[i, :]["_neighbors"].split(","):
                gap = 0
                self._distances.append(h)
                self._distances.append(v)
                # print("no gap:", self.gdf.at[i, "NAME"], "-", self.gdf.at[j, "NAME"])

            else:
                gap = self.gap_size
                # print("gap:", self.gdf.at[i, "NAME"], "-", self.gdf.at[j, "NAME"])

            # Check whether the neighbors should be connected horizontally or vertically.
            geom_dist_x = self.gdf.at[j, "geometry"].x - self.gdf.at[i, "geometry"].x
            geom_dist_y = self.gdf.at[j, "geometry"].y - self.gdf.at[i, "geometry"].y

            # Create a variable to determine whether there is a horizontal or vertical adjacency
            if s_dict is not None:
                # print(i, j)
                # print(s_dict)
                s = s_dict[i][j]
                one_minus_s = 1 - s
            if self.mode == 2 and s_dict is None:
                s_ = p.LpVariable(f"s_{i}_{j}")
            elif self.mode == 3:
                s_ = p.LpVariable(f"s_{i}_{j}", cat="Integer")
            if self.mode in (2, 3) and s_dict is None:
                s = p.LpAffineExpression([(s_, (w + gap))])
                one_minus_s = p.LpAffineExpression(
                    [(s_, -1 * (w + gap))], constant=(w + gap)
                )
                self.problem += s + one_minus_s <= w + gap
                self.problem += s + one_minus_s >= w + gap

            # Depending on whether x_i - x_j yields a positive or negative result we decide the order
            # in which they are presented in the constraints.
            if geom_dist_x > 0:
                if abs(geom_dist_x) > abs(geom_dist_y) and self.mode == 1:
                    # Add constraint 1
                    self.problem += (
                        self._coord_dict[j]["x"] - self._coord_dict[i]["x"] >= w + gap
                    )
                elif self.mode != 1:
                    # Add constraint 1
                    self.problem += (
                        self._coord_dict[j]["x"] - self._coord_dict[i]["x"] >= s
                    )
                # Add constraint 3
                self.problem += (
                    h >= self._coord_dict[j]["x"] - self._coord_dict[i]["x"] - w
                )
            else:
                if abs(geom_dist_x) > abs(geom_dist_y) and self.mode == 1:
                    # Add constraint 1
                    self.problem += (
                        self._coord_dict[i]["x"] - self._coord_dict[j]["x"] >= w + gap
                    )
                elif self.mode != 1:
                    # Add constraint 1
                    self.problem += (
                        self._coord_dict[i]["x"] - self._coord_dict[j]["x"] >= s
                    )
                # Add constraint 3
                self.problem += (
                    h >= self._coord_dict[i]["x"] - self._coord_dict[j]["x"] - w
                )

            # Depending on whether y_i - y_j yields a positive or negative result we decide the order
            # in which they are presented in the constraints.
            if geom_dist_y > 0:
                if not abs(geom_dist_x) > abs(geom_dist_y) and self.mode == 1:
                    # Add constraint 2
                    self.problem += (
                        self._coord_dict[j]["y"] - self._coord_dict[i]["y"] >= w + gap
                    )
                elif self.mode != 1:
                    # Add constraint 2
                    self.problem += (
                        self._coord_dict[j]["y"] - self._coord_dict[i]["y"]
                        >= one_minus_s
                    )
                # Add constraint 4
                self.problem += (
                    v >= self._coord_dict[j]["y"] - self._coord_dict[i]["y"] - w
                )
            else:
                if not abs(geom_dist_x) > abs(geom_dist_y) and self.mode == 1:
                    # Add constraint 2
                    self.problem += (
                        self._coord_dict[i]["y"] - self._coord_dict[j]["y"] >= w + gap
                    )
                elif self.mode != 1:
                    # Add constraint 2
                    self.problem += (
                        self._coord_dict[i]["y"] - self._coord_dict[j]["y"]
                        >= one_minus_s
                    )
                # Add constraint 4
                self.problem += (
                    v >= self._coord_dict[i]["y"] - self._coord_dict[j]["y"] - w
                )

        # Define the objective
        self.problem += p.lpSum(self._distances)

    def solve(self) -> None:
        """Method to solve for the optimal location of geometries given an input Geodataframe with a size column.

        """

        self.problem.solve(
            p.GLPK_CMD(msg=self.solve_msg, options=["--tmlim", f"{self.time_limit}"])
        )

        s_dict = {}
        for v in self.problem.variables():
            var_name = v.name.split("_")
            if self.mode == 2 and var_name[0] == "s":
                s, i, j = v.name.split("_")
                # print(s, i, j)
                if int(i) in s_dict:
                    s_dict[int(i)][int(j)] = round(v.varValue)
                else:
                    s_dict[int(i)] = {int(j): round(v.varValue)}
            if len(var_name) == 2:
                # All variables with two parts in their name are coordinates. The first part says if it's a x or y
                # coordinate, the second part the geometry it corresponds to (row number in Geodataframe).
                c, i = v.name.split("_")
                self.gdf.at[int(i), f"_{c}"] = v.varValue

        # print(s_dict)
        # if self.mode == 2:
        #     self.problem = p.LpProblem("Problem2", p.LpMinimize)
        #     self.add_coord_variables()
        #     self.add_adjacent_variables_constraints(s_dict=s_dict)
        #     self.problem.solve(
        #         p.GLPK_CMD(msg=True, options=["--tmlim", f"{self.time_limit}"])
        #     )
        #     for v in self.problem.variables():
        #         var_name = v.name.split("_")
        #         if len(var_name) == 2:
        #             # All variables with two parts in their name are coordinates. The first part says if it's a x or y
        #             # coordinate, the second part the geometry it corresponds to (row number in Geodataframe).
        #             c, i = v.name.split("_")
        #             self.gdf.at[int(i), f"_{c}"] = v.varValue

        self.gdf.geometry = self.gdf.apply(lambda r: Point(r["_x"], r["_y"]), axis=1)
