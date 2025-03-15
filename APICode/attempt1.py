from vapor import session, renderer, dataset

from vapor.animation import Animation

from vapor.dataset import Dataset
ses = session.Session()

print(Dataset.GetDatasetTypes())

data = ses.OpenDataset(dataset.BOV, ["header_file.bov"])

print("Data stuff:")

vars = ["BX"]

for var in data.GetDataVarNames():

    if var in vars:

        print(f" {var}")

        print(f"    Time Varying: False")

        print(f"    Dimensionality:", data.GetVarGeometryDim(var))

        print(f"    Coordinates:", data.GetVarCoordVars(var, True))

        print("     Data Range:", data.GetDataRange(var))


first_2d_var = vars[0]

print(f"Rendering 2D variable {first_2d_var}")



ren = data.NewRenderer(renderer.TwoDDataRenderer)

ren.SetVariableName(first_2d_var)

ren.GetPrimaryTransferFunction().SetMinMapValue(-1)

ren.GetPrimaryTransferFunction().SetMaxMapValue(1)



ses.GetCamera().ViewAll()

output_file = "output.png"

ses.Render(output_file)
