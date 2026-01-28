
from vapor import session, renderer, dataset
from vapor.dataset import Dataset

# Loop over time steps from 1 to 160 (inclusive)
for ts in range(10, 161): #161  
    ts_str = f"{ts:03d}"  # Ensure consistent zero-padded format if needed (e.g., "001", "002", ..., "160")

    # Initialize a new session for each timestep
    ses = session.Session()

    output_file = f"/project/cstr/Jaiman/fa25/Project1/Test2/Outputs/DensityFrontView/output_{ts_str}.png"

    if ts%1 == 0:

        try:
            # Open dataset with corresponding BOV files for the current timestep
            data = ses.OpenDataset(dataset.BOV, [
                f"/project/cstr/Jaiman/fa25/Project1/Test2/BOV/BX_{ts_str}.bov",
                f"/project/cstr/Jaiman/fa25/Project1/Test2/BOV/BY_{ts_str}.bov",
                f"/project/cstr/Jaiman/fa25/Project1/Test2/BOV/BZ_{ts_str}.bov",
                f"/project/cstr/Jaiman/fa25/Project1/Test2/BOV/VZ_{ts_str}.bov",
                f"/project/cstr/Jaiman/fa25/Project1/Test2/BOV/CB2_{ts_str}.bov",
                f"/project/cstr/Jaiman/fa25/Project1/Test2/BOV/CT_BT_{ts_str}.bov",
                f"/project/cstr/Jaiman/fa25/Project1/Test2/BOV/RO_{ts_str}.bov"
               # f"/project/cstr/Jaiman/sp25/VAPORdata/Test3/Images/BX_{ts_str}.bov",
               # f"/project/cstr/Jaiman/sp25/VAPORdata/Test3/Images/BY_{ts_str}.bov",
               # f"/project/cstr/Jaiman/sp25/VAPORdata/Test3/Images/BZ_{ts_str}.bov",
               # f"/project/cstr/Jaiman/sp25/VAPORdata/Test3/Images/VZ_{ts_str}.bov",
               # f"/project/cstr/Jaiman/sp25/VAPORdata/Test3/Images/CB2_{ts_str}.bov",
               # f"/project/cstr/Jaiman/sp25/VAPORdata/Test3/Images/CT_BT_{ts_str}.bov"
            ])

            # Create renderers
            ren1 = data.NewRenderer(renderer.SliceRenderer)

            ren2 = data.NewRenderer(renderer.SliceRenderer)

            #ren3 = data.NewRenderer(renderer.FlowRenderer)

            ren4 = data.NewRenderer(renderer.FlowRenderer)

            # Configure the first renderer 

            ren1.SetVariableName("BZ")

            ren1.SetZSlicePlaneOrigin(0)

            ren1.GetPrimaryTransferFunction().SetMinMapValue(-1)

            ren1.GetPrimaryTransferFunction().SetMaxMapValue(1)

            # Configure the second renderer

            ren2.SetVariableName("RO")

            ren2.SetZSlicePlaneOrigin(0.5) #halfway up the box

            ren2.SetYSlicePlaneRotation(-90) #rotated

            ren2.GetPrimaryTransferFunction().SetOpacityScale(0.5)

    #ren3 used to be here

            ren4.SetDimensions(3)

            ren4.SetFieldVariableNames(["BX","BY","BZ"])  #gives yz cross section, add in BX if want 3d field

            bounding_box2 = ren4.GetRakeRegion()
            val2 = bounding_box2.GetExtents()

            min_extents2 = [175, 175, 0]  
            max_extents2 = [185, 185, 360]

            bounding_box2.SetExtents(min_extents2, max_extents2)

            bounding_box2 = ren4.GetRakeRegion()
            val2 = bounding_box2.GetExtents()
            
            ren4.SetRenderType(0) #streamlines

            ren4.SetSeedGenMode(0) #gridded

            ren4.SetGridNumOfSeeds([5,5,10])

            ren4.GetPrimaryTransferFunction().LoadBuiltinColormap("Sequential/thermal")

            ren4.SetFlowDirection(2) #bidirectional

            #ren4.SetRenderFadeTailStart(1) #integration steps?
            #ren4.SetRenderFadeTailStop(10000)
            
            print(ren4.GetSteadyNumOfSteps())
            ren4.SetSteadyNumOfSteps(10000)

            ren4.SetRenderRadiusScalar(1) #size of flow line, default is around 1

            ren4.SetDimensions(3)

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
