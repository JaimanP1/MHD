from vapor import session, renderer, dataset
from vapor.animation import Animation
from vapor.dataset import Dataset

ses = session.Session()
print(Dataset.GetDatasetTypes())
data = ses.OpenDataset(dataset.BOV, ["params.bov"])
first_3d_var = "BZ"
print(f"Rendering 3D variable {first_3d_var}")

ren:renderer.SliceRenderer = data.NewRenderer(renderer.SliceRenderer)
#ren = data.NewRenderer(renderer.SliceRenderer)
ren.SetVariableName(first_3d_var)
#ren.SetSliceOffset(0)
ren.SetZSlicePlaneOrigin(0)
ren.GetPrimaryTransferFunction().SetMinMapValue(-1)
ren.GetPrimaryTransferFunction().SetMaxMapValue(1)
cam = ses.GetCamera()
cam.LookAt((32, -100, 100), ren.GetTransform().GetOrigin())
cam.Zoom(0.982)


output_file = "output.png"
ses.Render(output_file)


