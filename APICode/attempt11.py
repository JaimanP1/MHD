import numpy as np
import matplotlib.pyplot as plt
from vapor import session, renderer, dataset
from vapor.dataset import Dataset

# Loop over time steps from 1 to 160 (inclusive)
for ts in range(1, 50): #161  
    ts_str = f"{ts:03d}"  # Ensure consistent zero-padded format if needed (e.g., "001", "002", ..., "160")

    # Initialize a new session for each timestep
    ses = session.Session()

    output_file = f"/project/cstr/Jaiman/sp26/Project1/Test2/API/Outputs/DensityFrontView/output_{ts_str}.png"

    if ts%1 == 0:

        try:
            # Open dataset with corresponding BOV files for the current timestep
            data = ses.OpenDataset(dataset.BOV, [
                f"/project/cstr/Jaiman/sp26/Project1/Test2/API/BOV/RO_{ts_str}.bov"
            ])

            # Create renderers

            ren2 = data.NewRenderer(renderer.SliceRenderer)


            # Configure the second renderer

            data_range = data.GetDataRange("RO")

            max_val = data_range[1]

            ren2.SetVariableName("RO")

            ren2.SetZSlicePlaneOrigin(0.5) #halfway up the box

            ren2.SetYSlicePlaneRotation(-90) #rotated

            ren2.GetPrimaryTransferFunction().SetOpacityScale(0.5)

            #ren2.GetPrimaryTransferFunction().SetMaxMapValue(np.log10(max_val))

            cam = ses.GetCamera()

            cam.LookAt(camera_position=(1400, 180, 180), target=(180, 180, 180), up=(0, 0, 1))

            cam.Zoom(0.5)

            ses.Render(output_file)

        except Exception as e:
            print(f"Error processing timestep {ts_str}: {e}")

        # Clear session to prevent memory buildup
        ses = None  
        
        print("timestep", ts, "done")

    else:
        print("skipping", ts)

print("Finished")
