import numpy as np
import datetime
import matplotlib.pyplot as plt
from lib import device
import lib.ndt_acquisition as ndt

def main():
    probe_state = device.init_device(verbose=True, logging=False)
    PATH = ""
    PIEZOID = "10MHz delay"

    cent = 10e6
    
    gain = 350
    pon = 80
    poff = 80
    damp = 4000
    
    device.dac(probe_state, gain)
    device.pulse_adc_trigger(probe_state, pon=pon, poff=poff, damp=damp)
    raw_response = device.read_device(probe_state)
    
    raw_strings = [x.replace("b'", "") for x in str(raw_response[2]).split(",") if len(x)]
    signal = np.array([(int(x, 16) - 512) / 512.0 for x in raw_strings[:-1]], dtype=np.float32)

    #calib = ndt.calibrate(
    #    PATH,
    #    Fech=60e6,
    #    start_us=35, end_us=40, # echo of interest window
    #    gain=gain,
    #    piezo_central_freq=cent,
    #    piezo_bandwidth=4e6,
    #    piezo_id=PIEZOID,
    #    target="steel block - calibration",
    #    overwrite=False,
    #)
    #print("optimal pon=poff =", calib["best_pon_poff"])
   
    #a = ndt.from_probe(
    #    Fech=60e6, gain=250,                            # Fech is 60Msps, gain values can be from 0 to 500
    #    piezo_central_freq=5e6,
    #    piezo_bandwidth=2e6,
    #    piezo_id= PIEZOID,
    #    pon=calib["best_pon_poff"], 
    #    poff=calib["best_pon_poff"], 
    #    damp=6000,                                      # Pulse parameters, can be tuned to optimize the signal for a given target.
    #    target="Calibration 1018 steel block, ~" + str(THICKNESS * 1000) + "mm",   # Target description, will be saved in the h5 file and can be used to identify the acquisition later.
    #    h5_path=PATH,                                   # Where data is saved
    #    overwrite = True                               # If acquisition with these parameters already exists in the h5 file, 
    #                                                    # it will be loaded instead of acquired again.
    #)
    #a.plot();
    
    acq_data = {
        "signal": signal,
        "Fech": probe_state["Fech"],
        "pon": pon,
        "poff": poff,
        "damp": damp,
        "gain": gain,
        "target": "calibration_block",
        "timestamp": "",
        "piezo_id": "10MHz Delay",
        "piezo_central_freq": cent,
        "piezo_bandwidth": 3.0e6,
    }
    
    echo_results = ndt.detect_echoes(
        acq_data,
        start_us=28,
        end_us=40,
        target_thickness=0.01,
        speed_of_sound=5400.0, # water
        plot=True
    )
    
    #print(f"Calculated Thickness: {echo_results['calculated_thickness'] * 1000:.2f} mm")
    
    #ndt.save_acquisition(acq_data, path="calibration-data/acquisitions.h5")

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    
    filename = f"plots/plot_{timestamp}.png"
    plt.savefig(filename)

if __name__ == "__main__":
    main()
