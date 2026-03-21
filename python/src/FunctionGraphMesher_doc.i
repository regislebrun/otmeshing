%feature("docstring") OTMESHING::FunctionGraphMesher
"Function meshing algorithm.

Examples
--------
>>> import openturns as ot
>>> import otmeshing
>>> bbox = ot.Interval([-8] * 2, [8] * 2)
>>> inputDiscretization = [10] * 2
>>> mesher = otmeshing.FunctionGraphMesher(bbox, inputDiscretization)
>>> f = ot.SymbolicFunction(['x', 'y'], ['cos(x) * sin(y)'])
>>> zIndex = 2  # z axis index
>>> zMin, zMax = 0.01, 0.3  # clip f in [0.01, 0.3]
>>> mesh = mesher.build(f, zIndex, zMin, zMax)"

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::FunctionGraphMesher::build
"Generate a mesh from function discretization.

Parameters
----------
function : :py:class:`openturns.Function`
    Surface function, of output dimension 1
outputIndex : int
    Index of the output variable
minOutput : float
    Minimum output value
maxOutput : float
    Maximum output value
outputDiscretization : int, optional
    Number of discretizations along the output axis (default=1)
subGraph : bool, optional
    Whether to mesh below or above the surface (default=True: below)

Returns
-------
mesh : :py:class:`openturns.Mesh`
    The mesh generated."
