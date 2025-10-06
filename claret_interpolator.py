from general_purpose_python_modules.cmd_utils.modified_cmd_isochrone_interpolator import IsochroneFileIterator
from general_purpose_python_modules.cmd_utils.modified_cmd_isochrone_interpolator import CMDInterpolator as interpolator
import numpy

from astropy import units as u, constants as c

from general_purpose_python_modules import grid_tracks_interpolate

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
                    "Mass0 and metalicity of input data must be arranged on "
                    "a grid."
                )
        if section_index >= 1:
            if section_label[grid[0][0]] == grid[0][1][0]:
                grid[1][1].append(section_label[grid[1][0]])
            else:
                if (
                    section_label[grid[1][0]]
                    != grid[1][1][section_index % len(grid[1][1])]
                ):
                    raise ValueError(
                        "Metalicity and Mass0 of input data must be arranged"
                        f" on a grid. Section {section_index} {grid[1][0]} "
                        "should be "
                        f"{grid[1][1][section_index % len(grid[1][1])]} "
                        f"instead of {section_label[grid[1][0]]}."
                    )
                if section_label[grid[0][0]] != grid[0][1][-1]:
                    grid[0][1].append(section_label[grid[0][0]])
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


interpolator.__init__ = __init__
interpolator._get_interpolation_grid = _get_interpolation_grid
IsochroneFileIterator.__enter__ = __enter__
claret_interpolator = interpolator