
from vapor import session, renderer, dataset
from vapor.dataset import Dataset

# Loop over time steps from 1 to 160 (inclusive)
for ts in range(1, 161): #161  
    ts_str = f"{ts:03d}"  # Ensure consistent zero-padded format if needed (e.g., "001", "002", ..., "160")

    # Initialize a new session for each timestep
    ses = session.Session()

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
            ])

            print(f"density range: ", data.GetDataRange("RO")) 

        except Exception as e:
            print(f"Error processing timestep {ts_str}: {e}")

        # Clear session to prevent memory buildup
        ses = None  
        
        print("timestep", ts, "done")

    else:
        print("skipping", ts)

print("Finished")
