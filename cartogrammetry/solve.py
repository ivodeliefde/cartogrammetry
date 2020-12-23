#!/usr/bin/env python3

"""Module for solving the location of a blocks or circles in cartogram using linear programming.

"""
import logging
import geopandas as gpd
import pulp as p
from tqdm.auto import tqdm


class Solver:
    """Linear Programming Solver for creating cartograms.

    This class solves the location of geometries in a GeoDataFrame, in order to create a cartogram. The methodology is derived from the article 'Computing stable Demers cartograms' by Nickel et al. (2019).


**Methodology**

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

There are 4 modes implemented:

    - Mode 1: For any combination of p_i and p_j either constraint 1 or 2 must apply. Which one is used depends on which axis has the greatest range. E.g. if x_j - x_i > y_j - y_i, then constraints 1 applies and p_i and p_j are considered horizontally adjacent.
    - Mode 2: A variable S is created with a value between 0 and 1 to determine the most optimal horizontal or vertical position for adjacent geometries. The result is allowed to have overlapping geometries.
    - Mode 3: Similar to mode 2, but user an integer variable (MIP) instead of a float. This reduces overlap but takes more time to solve.
    - Mode 4: Similar to mode 2, but after the initial solved result the variable S is rounded and added to the problem as a constant.


    :param gdf: Geodataframe to optimize.
    :param gap_size: The minimum length in meters between non-adjacent geometries.
    :param time_limit: The maximum time allowed to find a solution. If exceeded the best solution at that time is returned.
    :param mode: Mode 1, 2 or 3. Corresponds to 1: quick method but suboptimal, 2: quick, optimal, but with overlapping geometries, 3: optimal, no overlapping geometries, but slow.
    :param solve_msg: Boolean value for displaying solver logging.
    """

    def __init__(
        self,
        gdf: gpd.GeoDataFrame,
        gap_size: float = 0.,
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
        logging.debug("Initialized Solver")
        logging.debug("Parameters")
        logging.debug(f"gdf: {type(self.gdf)} with shape {self.gdf.shape}")

    def run(self) -> None:
        """Methode to run all required steps to find the optimal centerpoint locations for each geometry in a cartogram.

        """
        logging.debug("method run called")
        self.add_coord_variables()
        if self.mode == 1:
            self.add_adjacent_variables_constraints_m1()
        elif self.mode == 2:
            self.add_adjacent_variables_constraints_m23()
        elif self.mode == 3:
            self.add_adjacent_variables_constraints_m23(mip=True)
        elif self.mode == 4:
            self.add_adjacent_variables_constraints_m23()
            s_dict = self.solve()

            # Second iteration with rounded S values as constants
            self._coord_dict = {}
            self._distances = []
            self.problem = p.LpProblem("Problem", p.LpMinimize)
            self.add_coord_variables()
            self.add_adjacent_variables_constraints_m4(s_dict)

        self.solve()

    def add_coord_variables(self, original_location=True) -> None:
        """Method to create LP variables for each coordinate.

For each geometry in the GeoDataFrame we create 2 variables:

    - X coordinate of centerpoint
    - Y coordinate of centerpoint

        """
        logging.debug("method add_coord_variables called")
        for i in range(self.gdf.shape[0]):
            self._coord_dict[i] = {}
            self._coord_dict[i]["x"] = p.LpVariable(f"x_{i}")
            self._coord_dict[i]["y"] = p.LpVariable(f"y_{i}")

            self.problem += self._coord_dict[i]["x"] >= 0
            self.problem += self._coord_dict[i]["y"] >= 0

            if original_location:
                # Add to the objective to also minimize the distance to the original location.
                dist_original_x = (self._coord_dict[i]["x"] - self.gdf.at[i, "geometry"].x) / self.gdf.loc[:, "geometry"].x.min()
                dist_original_y = (self._coord_dict[i]["y"] - self.gdf.at[i, "geometry"].y) / self.gdf.loc[:, "geometry"].y.min()

                # self.problem += dist_original_x >= 0
                # self.problem += dist_original_y >= 0
                #
                # self._distances.append(dist_original_x)
                # self._distances.append(dist_original_y)

    def add_adjacent_variables_constraints_m1(self) -> None:
        """Method to create LP variables for each pair of adjacent geometries (using mode 1).

        """
        logging.debug("method add_adjacent_variables_constraints_m1 called")
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
            perimeters = [np for np in self.gdf.iloc[i, :]["_neighbors"].split(",")]
            neighbors = [n for n in self.gdf.iloc[i, :]["_neighbors"].split(",")]
            if str(j) in neighbors or True:
                gap = 0
                # if float(perimeters[neighbors.index(str(j))]) > 0.1 or True:
                if True:
                    self._distances.append(h)
                    self._distances.append(v)
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
                            self._coord_dict[j]["x"] - self._coord_dict[i]["x"] >= w + gap
                    )
                # Add constraint 3
                self.problem += (
                        h >= self._coord_dict[j]["x"] - self._coord_dict[i]["x"] - w
                )
            else:
                if abs(geom_dist_x) > abs(geom_dist_y):
                    # Add constraint 1
                    self.problem += (
                            self._coord_dict[i]["x"] - self._coord_dict[j]["x"] >= w + gap
                    )
                # Add constraint 3
                self.problem += (
                        h >= self._coord_dict[i]["x"] - self._coord_dict[j]["x"] - w
                )

            # Depending on whether y_i - y_j yields a positive or negative result we decide the order
            # in which they are presented in the constraints.
            if geom_dist_y > 0:
                if not abs(geom_dist_x) > abs(geom_dist_y):
                    # Add constraint 2
                    self.problem += (
                            self._coord_dict[j]["y"] - self._coord_dict[i]["y"] >= w + gap
                    )
                # Add constraint 4
                self.problem += (
                        v >= self._coord_dict[j]["y"] - self._coord_dict[i]["y"] - w
                )
            else:
                if not abs(geom_dist_x) > abs(geom_dist_y):
                    # Add constraint 2
                    self.problem += (
                            self._coord_dict[i]["y"] - self._coord_dict[j]["y"] >= w + gap
                    )
                # Add constraint 4
                self.problem += (
                        v >= self._coord_dict[i]["y"] - self._coord_dict[j]["y"] - w
                )

            # Define the objective
        self.problem += p.lpSum(self._distances)

    def add_adjacent_variables_constraints_m23(self, mip: bool = False) -> None:
        """Method to create LP variables for each pair of adjacent geometries (using mode 2 or mode 3).

        """
        logging.debug("method add_adjacent_variables_constraints_m23 called")
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
            else:
                gap = self.gap_size

            # Create the variable s_ to determine horizontal or vertical adjacency in the cartogram.
            if mip:
                s_ = p.LpVariable(f"s_{i}_{j}", cat="Integer")
            else:
                s_ = p.LpVariable(f"s_{i}_{j}")
            s = p.LpAffineExpression(
                [(s_, (w + gap))]
            )
            one_minus_s = p.LpAffineExpression(
                [(s_, -1 * (w + gap))], constant=(w + gap)
            )
            self.problem += s_ >= 0
            self.problem += s_ <= 1
            self.problem += s + one_minus_s <= w + gap
            self.problem += s + one_minus_s >= w + gap

            # Check the distance between j en i to determine in which order to place the coordinates
            geom_dist_x = self.gdf.at[j, "geometry"].x - self.gdf.at[i, "geometry"].x
            geom_dist_y = self.gdf.at[j, "geometry"].y - self.gdf.at[i, "geometry"].y

            # Depending on whether x_i - x_j yields a positive or negative result we decide the order
            # in which they are presented in the constraints.
            if geom_dist_x > 0:
                # Add constraint 1
                self.problem += (
                        self._coord_dict[j]["x"] - self._coord_dict[i]["x"] >= s
                )

                # Add constraint 3
                self.problem += (
                        h >= self._coord_dict[j]["x"] - self._coord_dict[i]["x"] - w
                )
            else:
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
                # Add constraint 2
                self.problem += (
                        self._coord_dict[j]["y"] - self._coord_dict[i]["y"] >= one_minus_s
                )
                # Add constraint 4
                self.problem += (
                        v >= self._coord_dict[j]["y"] - self._coord_dict[i]["y"] - w
                )
            else:
                # Add constraint 2
                self.problem += (
                        self._coord_dict[i]["y"] - self._coord_dict[j]["y"] >= one_minus_s
                )
                # Add constraint 4
                self.problem += (
                        v >= self._coord_dict[i]["y"] - self._coord_dict[j]["y"] - w
                )

            # Define the objective
        self.problem += p.lpSum(self._distances)

    def add_adjacent_variables_constraints_m4(self, s_dict) -> None:
        """Method to create LP variables for each pair of adjacent geometries (using mode 4).

        """
        logging.debug("method add_adjacent_variables_constraints_m4 called")
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

            else:
                gap = self.gap_size

            # Create the S constant based on the results in the first iteration to determine horizontal or vertical adjacency in the cartogram.
            if s_dict[i][j] and str(j) in self.gdf.iloc[i, :]["_neighbors"].split(","):
                s = w
                one_minus_s = 0
            elif not s_dict[i][j] and str(j) in self.gdf.iloc[i, :]["_neighbors"].split(","):
                s = 0
                one_minus_s = w
            else:
                s = 0
                one_minus_s = 0

            # Check the distance between j en i to determine in which order to place the coordinates
            geom_dist_x = self.gdf.at[j, "geometry"].x - self.gdf.at[i, "geometry"].x
            geom_dist_y = self.gdf.at[j, "geometry"].y - self.gdf.at[i, "geometry"].y

            # Depending on whether x_i - x_j yields a positive or negative result we decide the order
            # in which they are presented in the constraints.
            if geom_dist_x > 0:
                if s > 0:
                    # Add constraint 1
                    self.problem += (
                            self._coord_dict[j]["x"] - self._coord_dict[i]["x"] >= s
                    )

                # Add constraint 3
                self.problem += (
                        h >= self._coord_dict[j]["x"] - self._coord_dict[i]["x"] - w
                )
            else:
                if s > 0:
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
                if one_minus_s > 0:
                    # Add constraint 2
                    self.problem += (
                            self._coord_dict[j]["y"] - self._coord_dict[i]["y"] >= one_minus_s
                    )
                # Add constraint 4
                self.problem += (
                        v >= self._coord_dict[j]["y"] - self._coord_dict[i]["y"] - w
                )
            else:
                if one_minus_s > 0:
                    # Add constraint 2
                    self.problem += (
                            self._coord_dict[i]["y"] - self._coord_dict[j]["y"] >= one_minus_s
                    )
                # Add constraint 4
                self.problem += (
                        v >= self._coord_dict[i]["y"] - self._coord_dict[j]["y"] - w
                )

            # Define the objective
        self.problem += p.lpSum(self._distances)

    def solve(self) -> dict:
        """Method to solve for the optimal location of geometries given an input Geodataframe with a size column.

        :return: Dictionary that contains all the values for all S variables.
        """
        logging.debug("method solve called")
        self.problem.solve(
            p.GLPK_CMD(msg=self.solve_msg, options=["--tmlim", f"{self.time_limit}"])
        )

        s_dict = {}
        for v in self.problem.variables():
            var_name = v.name.split("_")
            if var_name[0] == "s":
                s, i, j = v.name.split("_")
                # print(v.varValue)
                if int(i) in s_dict:
                    s_dict[int(i)][int(j)] = v.varValue >= 0.5
                else:
                    s_dict[int(i)] = {int(j): v.varValue >= 0.5}
            if len(var_name) == 2:
                # All variables with two parts in their name are coordinates. The first part says if it's a x or y
                # coordinate, the second part the geometry it corresponds to (row number in Geodataframe).
                c, i = v.name.split("_")
                self.gdf.at[int(i), f"_{c}"] = v.varValue

        return s_dict
