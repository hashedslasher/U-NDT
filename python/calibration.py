import matplotlib.pyplot as plt
from lib.ndt_acquisition import UltrasonicAcquisition

PATH = "calibration-data/calibration-1.h5"
PIEZOID = "Delay 10MHz"

THICKNESS = 0.010 # [m] thickness of the calibration block, used to estimate where the echo should be and to calibrate the acquisition parameters.

for gain in (0, 10, 20, 30, 40, 50):
    a = UltrasonicAcquisition.from_probe(
        piezo_id= PIEZOID,
        Fech=60e6, gain=gain*10,
        pon=70, poff=70, damp=6000,
        target="Calibration 1018 steel block, ~" + str(THICKNESS * 1000) + "mm",
        h5_path=PATH
    )


a = UltrasonicAcquisition.from_probe(
        Fech=60e6, gain=350,                            # Fech is 60Msps, gain values can be from 0 to 500
        piezo_central_freq=5e6,
        piezo_bandwidth=4e6,
        piezo_id= PIEZOID,
        pon=70, poff=70, damp=6000,                     # Pulse parameters, can be tuned to optimize the signal for a given target.
        target="Calibration 1018 steel block, ~10mm",   # Target description, will be saved in the h5 file and can be used to identify the acquisition later.
        h5_path=PATH,                                   # Where data is saved
        overwrite = True                               # If acquisition with these parameters already exists in the h5 file, 
                                                        # it will be loaded instead of acquired again.
    )
a.plot();


result = a.detect_echoes(
    start_us=35, end_us=55,     # window where the main echo is expected
    target_thickness=THICKNESS,    # [m]
    speed_of_sound=5950,       # [m/s]  longitudinal in steel
    smooth=15,
    smooth_kernel=7,
)
print(result["peak_times_us"])


calib = UltrasonicAcquisition.calibrate(
    PATH,
    Fech=60e6,
    start_us=35, end_us=40,        # echo of interest window
    gain=40,
    piezo_central_freq=5e6,
    piezo_bandwidth=4e6,
    piezo_id=PIEZOID,
    target="Calibration 1018 steel block - calibration",
    overwrite=False,
)
print("optimal pon=poff =", calib["best_pon_poff"])


THICKNESS = 0.010 # [m] thickness of the calibration block, used to estimate where the echo should be and to calibrate the acquisition parameters.
a = UltrasonicAcquisition.from_probe(
        Fech=60e6, gain=350,                            # Fech is 60Msps, gain values can be from 0 to 500
        piezo_central_freq=5e6,
        piezo_bandwidth=2e6,
        piezo_id= PIEZOID,
        pon=calib["best_pon_poff"], 
        poff=calib["best_pon_poff"], 
        damp=6000,                                      # Pulse parameters, can be tuned to optimize the signal for a given target.
        target="Calibration 1018 steel block, ~" + str(THICKNESS * 1000) + "mm",   # Target description, will be saved in the h5 file and can be used to identify the acquisition later.
        h5_path=PATH,                                   # Where data is saved
        overwrite = True                               # If acquisition with these parameters already exists in the h5 file, 
                                                        # it will be loaded instead of acquired again.
    )
a.plot();


result = a.detect_echoes(
    start_us=35, end_us=60,     # window where the main echo is expected
    target_thickness=THICKNESS,    # [m]
    speed_of_sound=5950,       # [m/s]  longitudinal in steel
    smooth=15,
    smooth_kernel=7,
)
print(result["peak_times_us"])



THICKNESS = 0.015 # [m] thickness of the calibration block, used to estimate where the echo should be and to calibrate the acquisition parameters.
a = UltrasonicAcquisition.from_probe(
        Fech=60e6, gain=250,                            # Fech is 60Msps, gain values can be from 0 to 500
        piezo_central_freq=5e6,
        piezo_bandwidth=2e6,
        piezo_id= PIEZOID,
        pon=calib["best_pon_poff"], 
        poff=calib["best_pon_poff"], 
        damp=6000,                                      # Pulse parameters, can be tuned to optimize the signal for a given target.
        target="Calibration 1018 steel block, ~" + str(THICKNESS * 1000) + "mm",   # Target description, will be saved in the h5 file and can be used to identify the acquisition later.
        h5_path=PATH,                                   # Where data is saved
        overwrite = True                               # If acquisition with these parameters already exists in the h5 file, 
                                                        # it will be loaded instead of acquired again.
    )
a.plot();



result = a.detect_echoes(
    start_us=25, end_us=60,     # window where the main echo is expected
    target_thickness=THICKNESS,    # [m]
    speed_of_sound=5950,       # [m/s]  longitudinal in steel
    smooth=15,
    smooth_kernel=7,
)
print(result["peak_times_us"])




THICKNESS = 0.025 # [m] thickness of the calibration block, used to estimate where the echo should be and to calibrate the acquisition parameters.
a = UltrasonicAcquisition.from_probe(
        Fech=60e6, gain=50,                            # Fech is 60Msps, gain values can be from 0 to 500
        piezo_central_freq=5e6,
        piezo_bandwidth=2e6,
        piezo_id= PIEZOID,
        pon=calib["best_pon_poff"], 
        poff=calib["best_pon_poff"], 
        damp=6000,                                      # Pulse parameters, can be tuned to optimize the signal for a given target.
        target="Calibration 1018 steel block, ~" + str(THICKNESS * 1000) + "mm",   # Target description, will be saved in the h5 file and can be used to identify the acquisition later.
        h5_path=PATH,                                   # Where data is saved
        overwrite = True                               # If acquisition with these parameters already exists in the h5 file, 
                                                        # it will be loaded instead of acquired again.
    )
a.plot();



result = a.detect_echoes(
    start_us=15, end_us=60,     # window where the main echo is expected
    target_thickness=THICKNESS,    # [m]
    speed_of_sound=5950,       # [m/s]  longitudinal in steel
    smooth=15,
    smooth_kernel=7,
)
print(result["peak_times_us"])
