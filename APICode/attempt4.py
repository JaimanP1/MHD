from vapor import session, renderer, dataset

from vapor.dataset import Dataset

# Initialize session

ses = session.Session()

# Open dataset with both BOV files

data = ses.OpenDataset(dataset.BOV, ["./binaries/ts_40/BX_040.bov", "./binaries/ts_40/BY_040.bov", "./binaries/ts_40/BZ_040.bov", "./binaries/ts_40/VZ_040.bov", "./binaries/ts_40/CB2_040.bov", "./binaries/ts_40/CT_BT_040.bov"])

# Create two slice renderers

ren1 = data.NewRenderer(renderer.SliceRenderer)

ren2 = data.NewRenderer(renderer.SliceRenderer)

ren3 = data.NewRenderer(renderer.FlowRenderer)

# Configure the first renderer (BZ)

ren1.SetVariableName("BZ")

ren1.SetZSlicePlaneOrigin(0)

ren1.GetPrimaryTransferFunction().SetMinMapValue(-1)

ren1.GetPrimaryTransferFunction().SetMaxMapValue(1)

# Configure the second renderer (VZ)

ren2.SetVariableName("CB2")

ren2.SetZSlicePlaneOrigin(0.5) #halfway up the box

ren2.SetYSlicePlaneRotation(-90) #rotated

ren2.GetPrimaryTransferFunction().SetMinMapValue(0) #change these according to data range, CB2 is always >= 0

ren2.GetPrimaryTransferFunction().SetMaxMapValue(1) #check max value

# Configure the third renderer (BY)

ren3.SetDimensions(3)

ren3.SetFieldVariableNames(["BX","BY","BZ"])

ren3.SetRenderType(0) #streamlines

ren3.SetSeedGenMode(2) #random with bias

ren3.SetRakeBiasVariable("CT_BT")

ren3.SetRakeBiasStrength(2)

ren3.SetRandomNumOfSeeds(200)

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

#ren3.SetVariableName("BY")

# Set camera

cam = ses.GetCamera()

cam.LookAt((60, -90, 45), ren1.GetTransform().GetOrigin())

cam.Zoom(0.981)

# Add both renderers to the session

#ses.AddRenderer(ren1)

#ses.AddRenderer(ren2)

# Render and display

output_file = "output.png"

ses.Render(output_file)

#ses.Show()

