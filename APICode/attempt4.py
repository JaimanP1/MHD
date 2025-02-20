from vapor import session, renderer, dataset

from vapor.dataset import Dataset

# Initialize session

ses = session.Session()

# Open dataset with both BOV files

data = ses.OpenDataset(dataset.BOV, ["./binaries/BZ_001.bov", "./binaries/VZ_040.bov", "./binaries/CB2_040.bov", "./binaries/CT_BT_040.bov"])

# Create two slice renderers

ren1 = data.NewRenderer(renderer.SliceRenderer)

ren2 = data.NewRenderer(renderer.SliceRenderer)

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

# Set camera

cam = ses.GetCamera()

cam.LookAt((60, -90, 45), ren1.GetTransform().GetOrigin())

cam.Zoom(0.982)

# Add both renderers to the session

#ses.AddRenderer(ren1)

#ses.AddRenderer(ren2)

# Render and display

output_file = "output.png"

ses.Render(output_file)

#ses.Show()

