import numpy as np
import datetime
import matplotlib.pyplot as plt
from lib import device
import lib.ndt_acquisition as ndt

def main():
    probe_state = device.init_device(verbose=True, logging=False)
    
    gain = 25
    pon = 80
    poff = 80
    damp = 6000
    
    device.dac(probe_state, gain)
    device.pulse_adc_trigger(probe_state, pon=pon, poff=poff, damp=damp)
    raw_response = device.read_device(probe_state)
    
    raw_strings = [x.replace("b'", "") for x in str(raw_response[2]).split(",") if len(x)]
    signal = np.array([(int(x, 16) - 512) / 512.0 for x in raw_strings[:-1]], dtype=np.float32)
    
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
        "piezo_central_freq": 10.0e6,
        "piezo_bandwidth": 3.0e6,
    }
    
    echo_results = ndt.detect_echoes(
        acq_data,
        start_us=1.5,
        end_us=15.0,
        target_thickness=0.025,
        speed_of_sound=3560.0, # copper
        plot=True
    )
    
    #print(f"Calculated Thickness: {echo_results['calculated_thickness'] * 1000:.2f} mm")
    
    #ndt.save_acquisition(acq_data, path="calibration-data/acquisitions.h5")

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    
    filename = f"plots/plot_{timestamp}.png"
    plt.savefig(filename)

if __name__ == "__main__":
    main()
