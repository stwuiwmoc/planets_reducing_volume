# %%
import importlib

import planets_optimize_myclass as pom

if __name__ == "__main__":
    importlib.reload(pom)

    CONSTS = pom.Constants(physical_radius=925e-3,
                           ignore_radius=25e-3,
                           pixel_number=1024,
                           zernike_max_degree=10,
                           offset_height_percent=2)

    input_filename0 = "../sag_integration_code/" + "mkfolder/psm_test_kagi_data_integration/0117e_height.csv"
    input_filename3 = "../sag_integration_code/" + "mkfolder/psm_test_kagi_data_integration/0117h_height.csv"

    mes03 = pom.CirclePathMeasurementReading(Constants=CONSTS,
                                             circle_path_radius=870e-3,
                                             original_csv_fpath=input_filename0,
                                             deformed_csv_fpath=input_filename3)
