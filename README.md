# RotorFlight Gyro Filter Visualization

This interactive web application demonstrates the digital filters used in RotorFlight firmware for processing gyro signals. It provides educational visualizations of how different filter types work and how they affect both time-domain and frequency-domain responses.

## Features

### 1. **Filter Pipeline Overview**
- Visual representation of the complete gyro signal processing pipeline
- Step-by-step breakdown of each filter stage
- Educational explanation of the signal flow

### 2. **Decimation Filter (Bessel 4th Order)**
- Demonstrates the Bessel 4th order low-pass filter used for downsampling
- Shows frequency response and step response
- Interactive cutoff frequency and sample rate controls
- Mathematical explanation of the transfer function

### 3. **RPM Filter (Dynamic Notch)**
- Shows how motor RPM affects notch filter center frequency
- Interactive motor RPM, blade ratio, and Q controls
- Real-time frequency calculation: `f = RPM × ratio / 60`
- Demonstrates motor vibration filtering

### 4. **Low-Pass Filters (LPF1 & LPF2)**
- Multiple filter types: PT1, PT2, PT3, Butterworth, Bessel, Damped
- Side-by-side comparison of different filter characteristics
- Interactive cutoff frequency controls
- Cascaded filter response visualization

### 5. **Static Notch Filters**
- Two configurable notch filters with independent center frequencies
- Q value calculation based on center and cutoff frequencies
- Cascaded response demonstration
- Time-domain response to sine wave inputs

### 6. **Dynamic Notch Filter (FFT-based)**
- FFT spectrum analysis with peak detection
- Automatic notch filter placement at detected frequencies
- Configurable number of notches and Q values
- Real-time adaptation simulation

### 7. **Complete Pipeline Simulation**
- End-to-end simulation of the entire filter chain
- Multiple input signal types:
  - White noise
  - Step response
  - Sine wave
  - Frequency chirp
  - Realistic gyro data
- Before vs after filtering comparison
- Filter stage contribution analysis

## How to Use

1. **Open the Application**
   - Open `gyro_filter_visualization.html` in a modern web browser
   - The application requires an internet connection for Plotly.js and Math.js libraries

2. **Navigate Sections**
   - Click on any filter section header to expand/collapse it
   - Each section contains interactive controls and visualizations

3. **Adjust Parameters**
   - Use sliders to change filter parameters in real-time
   - Watch the plots update automatically as you adjust settings
   - Compare different filter types and configurations

4. **Explore Different Views**
   - Switch between tabs to see different aspects of each filter
   - View frequency responses, time responses, and comparisons
   - Analyze the complete pipeline simulation

## Technical Details

### Filter Types Implemented

1. **PT1 Filter (First Order)**
   - Transfer function: `H(s) = α / (1 - (1-α)z⁻¹)`
   - Where `α = cutoff / (cutoff + Fs/2π)`

2. **PT2 Filter (Second Order)**
   - Cascaded PT1 filters with adjusted cutoff
   - Cutoff frequency multiplied by 1.553773974

3. **PT3 Filter (Third Order)**
   - Three cascaded PT1 filters
   - Cutoff frequency multiplied by 1.961459177

4. **Butterworth Filter**
   - Maximally flat frequency response
   - Q = 0.707106781 (1/√2)

5. **Bessel Filter**
   - Linear phase response
   - Q = 0.577350269 (1/√3)

6. **Notch Filter**
   - Transfer function: `H(s) = (s² + ω₀²) / (s² + ω₀s/Q + ω₀²)`
   - Q = f₀/(f₂ - f₁)

### Constants from RotorFlight Firmware

```javascript
// Bessel 4th order filter constants
const BESSEL_4A_Q = 0.805538282;
const BESSEL_4B_Q = 0.521934582;
const BESSEL_4A_C = 1.603357516;
const BESSEL_4B_C = 1.430171560;

// Standard filter Q values
const BUTTER_Q = 0.707106781;  // 1/√2
const BESSEL_Q = 0.577350269;  // 1/√3
const DAMPED_Q = 0.5;          // 1/√4
```

## Educational Value

This visualization helps understand:

- **Digital Filter Theory**: How different filter types affect signals
- **Frequency Domain Analysis**: Magnitude and phase responses
- **Time Domain Behavior**: Step responses and transient behavior
- **Filter Cascading**: How multiple filters work together
- **Real-world Applications**: How filters clean gyro signals for flight control

## Browser Compatibility

- Modern browsers with ES6 support
- Requires Plotly.js for interactive plots
- Requires Math.js for mathematical operations
- Tested on Chrome, Firefox, Safari, and Edge

## File Structure

```
├── gyro_filter_visualization.html    # Main HTML file
├── gyro_filter_visualization.js      # JavaScript implementation
└── README_GYRO_FILTERS.md           # This documentation
```

## Contributing

This visualization is based on the actual RotorFlight firmware implementation. The filter algorithms and constants are taken directly from the source code to ensure accuracy.

## License

This educational tool is provided as-is for learning purposes. The filter implementations are based on RotorFlight firmware which is licensed under the GNU General Public License v3.
