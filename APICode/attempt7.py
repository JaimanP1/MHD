from vapor import session, renderer, dataset

from vapor.dataset import Dataset

# Initialize session

ses = session.Session()
ts = "042"
output_file = f"./Test3/Outputs/output_{ts}.png"

# Open dataset with both BOV files
data = ses.OpenDataset(dataset.BOV, [f"./Test1/BrickOfValues/BX_{ts}.bov", f"./Test1/BrickOfValues/BY_{ts}.bov", f"./Test1/BrickOfValues/BZ_{ts}.bov", f"./Test1/BrickOfValues/VZ_{ts}.bov", f"./Test1/BrickOfValues/CB2_{ts}.bov", f"./Test1/BrickOfValues/CT_BT_{ts}.bov"])

print("Data Variables:")
vars = ["BX", "BY", "BZ", "VZ", "CT_BT", "CB2"]
for var in data.GetDataVarNames():
    if var in vars:
        print(f" {var}, {ts}")
        print(f"    Time Varying: False")
        print(f"    Dimensionality:", data.GetVarGeometryDim(var))
        print(f"    Coordinates:", data.GetVarCoordVars(var, True))
        print("     Data Range:", data.GetDataRange(var))
        print(f"    Dimension Length (x=y=z):", data.GetDimensionLength("x", 1))

# Create renderers

ren1 = data.NewRenderer(renderer.SliceRenderer)

ren2 = data.NewRenderer(renderer.SliceRenderer)

ren3 = data.NewRenderer(renderer.FlowRenderer)

# Configure the first renderer 

ren1.SetVariableName("BZ")

ren1.SetZSlicePlaneOrigin(0)

ren1.GetPrimaryTransferFunction().SetMinMapValue(-1)

ren1.GetPrimaryTransferFunction().SetMaxMapValue(1)

# Configure the second renderer

ren2.SetVariableName("CB2")

ren2.SetZSlicePlaneOrigin(0.5) #halfway up the box

ren2.SetYSlicePlaneRotation(-90) #rotated

ren2.GetPrimaryTransferFunction().SetMinMapValue(0) #change these according to data range, CB2 is always >= 0

ren2.GetPrimaryTransferFunction().SetMaxMapValue(1) #check max value, do not need entire range

ren2.GetPrimaryTransferFunction().SetOpacityScale(0.5)

axis = ses.GetAxisAnnotations() #axis labeling

axis.SetAxisAnnotationEnabled(True)

#axis.SetMaxTics((360,360,360))

#axis.SetMinTics((0,0,0))

axis.SetNumTics((9,9,9)) #number of ticks

#axis.SetAxisBackgroundColor((0,0,0)) #black background

#axis.SetAxisFontSize(18)

# Configure the third renderer

ren3.SetDimensions(3)

ren3.SetFieldVariableNames(["","BY","BZ"])  #gives yz cross section, add in BX if want 3d field

bounding_box = ren3.GetRakeRegion()
val = bounding_box.GetExtents()
print("boundaries:", val)

min_extents = [180, 130.0, .0]  
max_extents = [180, 230.0, 360.0]

bounding_box.SetExtents(min_extents, max_extents)

bounding_box = ren3.GetRakeRegion()
val = bounding_box.GetExtents()
print("new boundaries: ", val)

ren3.SetRenderType(0) #streamlines

ren3.SetSeedGenMode(2) #random with bias

ren3.SetRakeBiasVariable("CT_BT")

ren3.GetPrimaryTransferFunction().LoadBuiltinColormap("Sequential/thermal")

ren3.GetPrimaryTransferFunction().SetMinMapValue(-.35)

ren3.GetPrimaryTransferFunction().SetMaxMapValue(.35) #check max value of variable(s). make dynamic, so data.maxvalue(bz)

ren3.SetRakeBiasStrength(1)

ren3.SetRandomNumOfSeeds(100)

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

cam.LookAt(camera_position=(1200, 180, 180), target=(180, 180, 180), up=(0, 0, 1))

cam.Zoom(0.5)

#cam.ViewAll()

ses.Render(output_file)

ses.Show()

