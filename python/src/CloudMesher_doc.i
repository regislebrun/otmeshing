%feature("docstring") OTMESHING::CloudMesher
"Mesher from a set of points.

Examples
--------
Triangulate a set of points:

>>> import openturns as ot
>>> import otmeshing
>>> points = ot.JointDistribution([ot.Uniform()] * 2).getSample(15)
>>> mesher = otmeshing.CloudMesher()
>>> triangulation = mesher.build(points)"

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::CloudMesher::build
"Triangulate a set of points.

Parameters
----------
points : :class:`~openturns.Sample`
    A set of points.

Returns
-------
mesh : :class:`~openturns.Mesh`
    The mesh built."
