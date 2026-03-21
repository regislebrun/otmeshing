import openturns as ot
import openturns.testing as ott
import otmeshing as otm

a = [-4.0] * 3
b = [4.0] * 3
f = ot.SymbolicFunction(["x0", "x1"], ["cos(pi_*x0)*sin(pi_*x1)^2"])
f.setName("Paraboloid")
f.setInputDescription([r"$x_0$", r"$x_1$"])
f.setOutputDescription([r"$x_2$"])
inputInterval = ot.Interval([a[0], a[1]], [b[0], b[1]])
inputDiscretization = [100] * 2
mesher = otm.FunctionGraphMesher(inputInterval, inputDiscretization)
print(mesher)
outputIndex = 2
outputDiscretization = 1
subGraph = True
mesh = mesher.build(f, outputIndex, a[2], b[2], outputDiscretization, subGraph)
# mesh.exportToVTKFile("func_graph.vtk")
print(mesh)
assert mesh.getDimension() == 3
assert mesh.getVerticesNumber() == 20402
assert mesh.getSimplicesNumber() == 60000
ott.assert_almost_equal(mesh.getVolume(), 66.2493)
