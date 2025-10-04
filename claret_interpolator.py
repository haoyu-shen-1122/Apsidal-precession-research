#!/usr/bin/env python3
"""Define a class that interpolates within a CMD isochrone grid."""

import os.path

from matplotlib import pyplot
import numpy
from astropy import units as u, constants as c

from general_purpose_python_modules import grid_tracks_interpolate


class IsochroneFileIterator:
    """
    Iterate over the sections of the isochrone grid file with line generators.

    Attributes:
        _isochrone:    The opened isochrone file to iterate over.

        _line:    The last line read from the file.

        header:    The collection of comment lines in the beginning of the
            isochrone file.
    """

    def __init__(self, isochrone_fname):
        """Create the iterator."""

        self._line = ""
        self.header = []
        self._isochrone_fname = isochrone_fname
        self._isochrone = None

    def __enter__(self):
        """Just return self."""

        self._isochrone = open(self._isochrone_fname, "r", encoding="utf-8")
        while True:
            self._line = self._isochrone.readline()
            assert self._line[0] == "#" or self._line.startswith("---") or self._line.startswith("yr|") or self._line.startswith('Age')
            if self._line.startswith("---"):
                break
            self.header.append(self._line)

        return self

    def __iter__(self):
        """Just return self."""

        return self

    def __next__(self):
        """Return a generator to the next section in the isochrone grid."""

        def section_line_generator():
            """Generator over lines of a section of the isochrone grid file."""

            assert self._line[0] == "#" or self._line.startswith("Age") or self._line.startswith("---") or self._line.startswith("yr|")

            yield self._line
            self._line = self._isochrone.readline()
            while self._line == "" or self._line[0] != "#":
                if self._line == "" or self._line.startswith("Age") or self._line.startswith("---"):
                    return
                yield self._line
                self._line = self._isochrone.readline()

        while self._line and self._line != "\n":
            return section_line_generator()
        raise StopIteration()

    def __exit__(self, *_):
        """Close the underlying isochrone file."""

        self._isochrone.close()


class CMDInterpolator:
    """
    Use interpolation to estimate quantities at any mass/[Fe/H] from a CMD grid.

    Attributes:
        header([str]):    The comment lines in the beginning of the isochrone
            file.

        isochrone_fname(str):    The name of the file from which interpolation
            data was read.

        grid([str, []):    The names of the indepnedent interpolation variables
            and the values of these variables at the grid points.

        _data([numpy field array]):    A the data contained in the isochrone
            file downloaded from the CMD interface, organized as a list of
            numpy field arrays, one for each section (corresponding to a
            single combination of [Fe/H] and age)
    """

    def _validate_input_data(self):
        """Check if input data satisfies the assumptions for interpolation."""

        valid = True
        message = (
            "Non-monotonic initial mass in isochrone "
        )

        for section_index, section_data in enumerate(self._data):
            invalid = (
                section_data[self.interp_var][1:]
                - section_data[self.interp_var][:-1]
            ) < 0
            if invalid.any():
                for bad_index in numpy.nonzero(invalid):
                    message += (
                        f"Section {section_index} m[{bad_index}] = "
                        f"{section_data[self.interp_var][bad_index]}, "
                        f"m[{bad_index + 1}] = "
                        f"{section_data[self.interp_var][bad_index + 1]}"
                        "\n"
                    )
                valid = False
        if not valid:
            raise ValueError(message)

    def _get_interpolation_grid(self):
        """Sanity check on the input data and prepare the interpolation."""

        for section_index, section_data in enumerate(self._data):
            for quantity in ["MH", "Mass0"]:
                assert numpy.unique(section_data[quantity]).size == 1
            section_label = {
                quantity: float(section_data[quantity][0])
                for quantity in ["MH", "Mass0"]
            }

            if section_index == 0:
                first_label = section_label
            elif section_index == 1:
                # False positive
                # pylint: disable=unsubscriptable-object

                grid = tuple(
                    (
                        quantity,
                        [first_label[quantity]],
                    )
                    for quantity in (
                        ["MH", "Mass0"]
                        if section_label["MH"] == first_label["MH"]
                        else ["Mass0", "MH"]
                    )
                )
                if (first_label["MH"] == section_label["MH"]) == (
                    first_label["Mass0"] == section_label["Mass0"]
                ):
                    raise ValueError(
                        "Age and metallicity of input data must be arranged on "
                        "a grid."
                    )
                # pylint: enable=unsubscriptable-object
            if section_index >= 1:
                # False positive
                # pylint: disable=used-before-assignment
                if section_label[grid[0][0]] == grid[0][1][0]:
                    grid[1][1].append(section_label[grid[1][0]])
                else:
                    if (
                        section_label[grid[1][0]]
                        != grid[1][1][section_index % len(grid[1][1])]
                    ):
                        raise ValueError(
                            "Age and metallicity of input data must be arranged"
                            f" on a grid. Section {section_index} {grid[1][0]} "
                            "should be "
                            f"{grid[1][1][section_index % len(grid[1][1])]} "
                            f"instead of {section_label[grid[1][0]]}."
                        )
                    if section_label[grid[0][0]] != grid[0][1][-1]:
                        grid[0][1].append(section_label[grid[0][0]])
                # pylint: enable=used-before-assignment
        return grid

    def __init__(self, data, interp_var="logAge"):
        """Interpolate within the given isochrone grid."""

        self._data = data
        self.interp_var = interp_var
        self._validate_input_data()
        self._grid = tuple(
            (quantity, numpy.array(values))
            for quantity, values in self._get_interpolation_grid()
            if len(values) > 1
        )

    def __call__(self, quantities, **interpolate_to):
        """
        Estimate the given quantities for a star of given initial mass & [Fe/H].

        Args:
            quantities(iterable of str):    The quantities to estimate. Should
                be a subset of the columns in the CMD isochrone file.

            interpolate_to(dict):    The value of [Z/H], logAge and Mini/Mass of
                the star for which to estimate the given quantities. Values for
                quantities that had only a single value in the input data are
                ignored.

        Return:
            tuple of floats:
                The values of the quantities estimated using linear
                interpolation in all dimensions.
        """

        return grid_tracks_interpolate(
            interpolate_to, quantities, self._grid, self._data
        )

    def get_mini_range(self, meh, mass):
        """Return the available mass range for the given [M/H] and log(age)."""

        return grid_tracks_interpolate(
            {"MH": meh, "Mass0": mass, "logAge": "range"},
            None,
            self._grid,
            self._data,
        )

    def get_log_age_range(self, logage, meh):
        """Return the available log(age) range for the given [M/H] and Mini."""

        if self._grid[0][0] == "MH":
            meh_values = self._grid[0][1]
            age_values = self._grid[1][1]
            meh_step = self._grid[1][1].size
            age_step = 1
        else:
            age_values = self._grid[0][1]
            meh_values = self._grid[1][1]
            meh_step = 1
            age_step = self._grid[1][1].size
        closest_meh = numpy.searchsorted(meh_values, meh)
        if meh_values[closest_meh] == meh:
            closest_meh = [closest_meh]
        elif closest_meh == 0:
            raise ValueError(
                f"Metallicity {meh} is below the range of the grid."
            )
        else:
            closest_meh = [closest_meh - 1, closest_meh]
        min_log_age = -numpy.inf
        max_log_age = numpy.inf
        for meh_index in closest_meh:
            min_log_age = max(
                min_log_age, self._data[meh_index * meh_step]["Mass0"][0]
            )
            age_index = 0
            while (
                self._data[meh_index * meh_step + age_index * age_step]["logAge"][
                    -1
                ]
                > logAge
                and age_index < age_values.size
            ):
                age_index += 1
            if age_index == 0:
                raise ValueError(
                    f"logAge {logage} is above the range of the grid for "
                    f"[M/H] = {meh}."
                )
            max_log_age = min(max_log_age, age_values[age_index - 1])

        return min_log_age, max_log_age

    def get_range(self, quantity):
        """Return the available interpolation range of the given quantity."""

        if quantity == self.interp_var:
            return (
                self._data[0][self.interp_var][0],
                self._data[0][self.interp_var][-1],
            )
        for name, values in self._grid:
            if name == quantity:
                return values[0], values[-1]
        raise ValueError(f"Unknown grid quantity {quantity}!")


def plot_isochrone(cmd_fname):
    """Plot an isochrone read from the CMD interface."""

    interpolator = CMDInterpolator(cmd_fname)
    print(
        "4.7Gyr Sun Teff = "
        + repr(
            10.0
            ** interpolator(
                ("logTe",), Mini=1.0, MH=0.0, logAge=9.0 + numpy.log10(4.7)
            )[0]
        )
    )

    log_age = numpy.linspace(0.0, 1.0, 1000)
    track_values = [
        interpolator(
            ("logTe", "logL", "logg", "Mass"),
            Mini=1.0,
            MH=0.0,
            logAge=9.0 + logt,
        )
        for logt in log_age
    ]
    track = {
        "logTe": numpy.array([v[0] for v in track_values]),
        "logL": numpy.array([v[1] for v in track_values]),
        "logg": numpy.array([v[2] for v in track_values]),
        "Mass": numpy.array([v[3] for v in track_values]),
    }

    pyplot.plot(-track["logTe"], track["logL"], "-k")
    pyplot.xlabel("-Teff")
    pyplot.ylabel("L")
    pyplot.show()

    pyplot.plot(
        log_age,
        numpy.sqrt(
            c.G
            * track["Mass"]
            * c.M_sun
            / (10 ** track["logg"] * u.cm / u.s**2)
        ).to_value(u.R_sun),
        "-r",
    )
    pyplot.plot(
        log_age,
        (
            10.0 ** (track["logL"] / 2.0 - 2.0 * track["logTe"])
            / (2.0 * u.K**2)
            * numpy.sqrt(c.L_sun / (numpy.pi * c.sigma_sb))
        ).to_value(u.R_sun),
        "-g",
    )
    pyplot.show()

    pyplot.plot(
        log_age,
        (
            (
                10.0 ** (track["logL"] / 2.0 - 2.0 * track["logTe"])
                / (2.0 * u.K**2)
                * numpy.sqrt(c.L_sun / (numpy.pi * c.sigma_sb))
            )
            / numpy.sqrt(
                c.G
                * track["Mass"]
                * c.M_sun
                / (10 ** track["logg"] * u.cm / u.s**2)
            )
        ).to_value(""),
        "-k",
    )
    pyplot.show()


def plot_age_range(cmd_fname, meh, mini_range):
    """Plot the available log(age) range for the given [M/H] and Mini."""

    interpolator = CMDInterpolator(cmd_fname)
    plot_x = numpy.linspace(*mini_range, 1000)
    logage_range = numpy.array(
        [list(interpolator.get_log_age_range(mini, meh)) for mini in plot_x]
    ).T
    pyplot.fill_between(
        plot_x,
        logage_range[0],
        logage_range[1],
        color="lightblue",
        label=f"[M/H] = {meh}",
    )
    pyplot.show()

