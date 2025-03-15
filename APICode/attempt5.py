from vapor import session, renderer, dataset

from vapor.dataset import Dataset

# Initialize session

ses = session.Session()

# Open dataset with both BOV files

data = ses.OpenDataset(dataset.BOV, ["./binaries/ts_40/BX_040.bov", "./binaries/ts_40/BY_040.bov", "./binaries/ts_40/BZ_040.bov", "./binaries/ts_40/VZ_040.bov", "./binaries/ts_40/CB2_040.bov", "./binaries/ts_40/CT_BT_040.bov"])

print("Data Variables:")
vars = ["BZ", "VZ", "CT_BT", "CB2"]
for var in data.GetDataVarNames():
    if var in vars:
        print(f" {var}")
        print(f"    Time Varying: False")
        print(f"    Dimensionality:", data.GetVarGeometryDim(var))
        print(f"    Coordinates:", data.GetVarCoordVars(var, True))
        print("     Data Range:", data.GetDataRange(var))
        print(f"    Dimension Length (x=y=z):", data.GetDimensionLength("x", 1))

ren3 = data.NewRenderer(renderer.FlowRenderer)

ren3.SetDimensions(3)

ren3.SetFieldVariableNames(["BX", "BY","BZ"])

bounding_box = ren3.GetRakeRegion()
bbox = renderer.BoundingBox(bounding_box)
val = bbox.GetExtents()
print("boundaries:", val)

min_extents = [179, 130.0, 1.0]  # Set your min bounds
max_extents = [181, 230.0, 10.0]  # Set your max bounds

#bbox.SetExtents(min_extents, max_extents)

ren3.SetRenderType(0) #streamlines

ren3.SetSeedGenMode(2) #random with bias

ren3.SetRakeBiasVariable("CT_BT")

ren3.SetRakeBiasStrength(100)

ren3.SetRandomNumOfSeeds(300)

ren3.SetFlowDirection(2) #bidirectional

ren3.SetRenderFadeTailStart(1) #integration steps?
ren3.SetRenderFadeTailStop(100)

a = ren3.GetRenderFadeTailStart()
b = ren3.GetRenderFadeTailStop()
print(f"integration: {a,b} (?)")

c = ren3.GetIsSteady()
print(f"should be true: {c}")
d = ren3.GetRenderFadeTailLength()
print(f"length: {d}")
e = ren3.GetSteadyNumOfSteps()
print(f"num of steps: {e}")

ren3.SetSteadyNumOfSteps(100) #integration steps?

cam = ses.GetCamera()

cam.LookAt((60, -90, 45), ren3.GetTransform().GetOrigin())

cam.Zoom(0.981)

output_file = "output.png"

ses.Render(output_file)

ses.Show()

