# U-NDT KSU

## How to Run the Acquisition Script

1. **Clone the repo**
   Clone the repo (or download as zip) and cd into u-ndt.

2. **Python**
   - **Windows:** Download from [python.org](https://www.python.org/downloads/windows/) or run:
     ```bash
     winget install -e --id Python.Python.3.14
     ```
   - **MacOS:** Download from [python.org](https://www.python.org/downloads/macos/) or run:
     ```bash
     brew install python3
     ```
   - Open PowerShell or Terminal and run:
     ```bash
     pip install pyserial numpy matplotlib scippy h5py picodev
     ```

3. **Flash the firmware**
   - Hold `BOOTSEL` and plug the board in. It should appear as `RPI-RP2`.
   - Drag and drop `firmware/pico/rp2040.uf2` (Pico) or `firmware/pico/rp2350.uf2` (Pico 2) into the `RPI-RP2` drive.
   - Or with picotool (probably have to run as admin):
   ```bash
   picotool load firmware/pico/rp2040.uf2
   ```

4. **Run the script**
   In PowerShell or Terminal:
   ```bash
   python3 python/{script}.py
   ```

If you use nix for some reason just:
```bash
nix develop .
python3 python/{script}.py
   ```
