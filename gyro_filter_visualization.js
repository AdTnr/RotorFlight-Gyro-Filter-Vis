// Gyro Filter Visualization - Filter Implementations and Plotting Functions

// Constants from RotorFlight firmware
const BESSEL_4A_Q = 0.805538282;
const BESSEL_4B_Q = 0.521934582;
const BESSEL_4A_C = 1.603357516;
const BESSEL_4B_C = 1.430171560;
const BUTTER_Q = 0.707106781;
const BESSEL_Q = 0.577350269;
const DAMPED_Q = 0.5;

// Filter implementation classes
class BiquadFilter {
    constructor() {
        this.x1 = this.x2 = 0;
        this.y1 = this.y2 = 0;
        this.b0 = this.b1 = this.b2 = 0;
        this.a1 = this.a2 = 0;
    }
    
    init(cutoff, sampleRate, Q, filterType) {
        const omega = 2 * Math.PI * cutoff / sampleRate;
        const sinom = Math.sin(omega);
        const cosom = Math.cos(omega);
        const alpha = sinom / (2 * Q);
        
        switch (filterType) {
            case 'LPF':
                this.b1 = 1 - cosom;
                this.b0 = this.b1 / 2;
                this.b2 = this.b0;
                this.a1 = -2 * cosom;
                this.a2 = 1 - alpha;
                break;
            case 'NOTCH':
                this.b0 = 1;
                this.b1 = -2 * cosom;
                this.b2 = 1;
                this.a1 = this.b1;
                this.a2 = 1 - alpha;
                break;
        }
        
        const a0 = 1 + alpha;
        this.b0 /= a0;
        this.b1 /= a0;
        this.b2 /= a0;
        this.a1 /= a0;
        this.a2 /= a0;
    }
    
    apply(input) {
        const output = this.b0 * input + this.b1 * this.x1 + this.b2 * this.x2 - 
                      this.a1 * this.y1 - this.a2 * this.y2;
        
        this.x2 = this.x1;
        this.x1 = input;
        this.y2 = this.y1;
        this.y1 = output;
        
        return output;
    }
}

class PT1Filter {
    constructor() {
        this.y1 = 0;
        this.gain = 0;
    }
    
    init(cutoff, sampleRate) {
        const gamma = sampleRate / (2 * Math.PI);
        this.gain = cutoff / (cutoff + gamma);
        this.gain = Math.min(this.gain, 1.0);
    }
    
    apply(input) {
        this.y1 += (input - this.y1) * this.gain;
        return this.y1;
    }
}

class PT2Filter {
    constructor() {
        this.y1 = this.y2 = 0;
        this.gain = 0;
    }
    
    init(cutoff, sampleRate) {
        const gamma = sampleRate / (2 * Math.PI);
        this.gain = (cutoff * 1.553773974) / ((cutoff * 1.553773974) + gamma);
        this.gain = Math.min(this.gain, 1.0);
    }
    
    apply(input) {
        this.y2 += (input - this.y2) * this.gain;
        this.y1 += (this.y2 - this.y1) * this.gain;
        return this.y1;
    }
}

class PT3Filter {
    constructor() {
        this.y1 = this.y2 = this.y3 = 0;
        this.gain = 0;
    }
    
    init(cutoff, sampleRate) {
        const gamma = sampleRate / (2 * Math.PI);
        this.gain = (cutoff * 1.961459177) / ((cutoff * 1.961459177) + gamma);
        this.gain = Math.min(this.gain, 1.0);
    }
    
    apply(input) {
        this.y3 += (input - this.y3) * this.gain;
        this.y2 += (this.y3 - this.y2) * this.gain;
        this.y1 += (this.y2 - this.y1) * this.gain;
        return this.y1;
    }
}

// Filter response calculation functions
function calculateFrequencyResponse(filterType, cutoff, sampleRate, Q = 1) {
    const frequencies = [];
    const magnitudes = [];
    const phases = [];
    
    for (let f = 1; f <= sampleRate / 2; f += sampleRate / 2000) {
        const omega = 2 * Math.PI * f / sampleRate;
        const z = { real: Math.cos(omega), imag: Math.sin(omega) };
        
        let H = { real: 0, imag: 0 };
        switch (filterType) {
            case 'PT1':
                const alpha = cutoff / (cutoff + sampleRate / (2 * Math.PI));
                const denom_pt1 = 1 - (1 - alpha) * z.real;
                const denomImag_pt1 = -(1 - alpha) * z.imag;
                const denomMag_pt1 = denom_pt1 * denom_pt1 + denomImag_pt1 * denomImag_pt1;
                H.real = alpha * denom_pt1 / denomMag_pt1;
                H.imag = alpha * denomImag_pt1 / denomMag_pt1;
                break;
            case 'PT2':
                const alpha2 = (cutoff * 1.553773974) / ((cutoff * 1.553773974) + sampleRate / (2 * Math.PI));
                const denom2 = 1 - (1 - alpha2) * z.real;
                const denomImag2 = -(1 - alpha2) * z.imag;
                const denomMag2 = denom2 * denom2 + denomImag2 * denomImag2;
                const H2 = {
                    real: alpha2 * denom2 / denomMag2,
                    imag: alpha2 * denomImag2 / denomMag2
                };
                // Square the complex number
                H.real = H2.real * H2.real - H2.imag * H2.imag;
                H.imag = 2 * H2.real * H2.imag;
                break;
            case 'PT3':
                const alpha3 = (cutoff * 1.961459177) / ((cutoff * 1.961459177) + sampleRate / (2 * Math.PI));
                const denom3 = 1 - (1 - alpha3) * z.real;
                const denomImag3 = -(1 - alpha3) * z.imag;
                const denomMag3 = denom3 * denom3 + denomImag3 * denomImag3;
                const H3 = {
                    real: alpha3 * denom3 / denomMag3,
                    imag: alpha3 * denomImag3 / denomMag3
                };
                // Cube the complex number
                const H3Squared = {
                    real: H3.real * H3.real - H3.imag * H3.imag,
                    imag: 2 * H3.real * H3.imag
                };
                H.real = H3Squared.real * H3.real - H3Squared.imag * H3.imag;
                H.imag = H3Squared.real * H3.imag + H3Squared.imag * H3.real;
                break;
            case 'BUTTER':
                const butterQ = 0.707106781;
                const butterOmega = 2 * Math.PI * cutoff / sampleRate;
                const butterAlpha = Math.sin(butterOmega) / (2 * butterQ);
                const b0 = (1 - Math.cos(butterOmega)) / 2;
                const b1 = 1 - Math.cos(butterOmega);
                const b2 = b0;
                const a1 = -2 * Math.cos(butterOmega);
                const a2 = 1 - butterAlpha;
                const a0 = 1 + butterAlpha;
                
                // Complex arithmetic for transfer function
                const z2 = { real: z.real * z.real - z.imag * z.imag, imag: 2 * z.real * z.imag };
                const num = { real: b0 + b1 * z.real + b2 * z2.real, imag: b1 * z.imag + b2 * z2.imag };
                const denom_butter = { real: a0 + a1 * z.real + a2 * z2.real, imag: a1 * z.imag + a2 * z2.imag };
                const denomMag_butter = denom_butter.real * denom_butter.real + denom_butter.imag * denom_butter.imag;
                H.real = (num.real * denom_butter.real + num.imag * denom_butter.imag) / denomMag_butter;
                H.imag = (num.imag * denom_butter.real - num.real * denom_butter.imag) / denomMag_butter;
                break;
            case 'BESSEL':
                const besselQ = 0.577350269;
                const besselOmega = 2 * Math.PI * cutoff / sampleRate;
                const besselAlpha = Math.sin(besselOmega) / (2 * besselQ);
                const b0_bessel = (1 - Math.cos(besselOmega)) / 2;
                const b1_bessel = 1 - Math.cos(besselOmega);
                const b2_bessel = b0_bessel;
                const a1_bessel = -2 * Math.cos(besselOmega);
                const a2_bessel = 1 - besselAlpha;
                const a0_bessel = 1 + besselAlpha;
                
                const z2_bessel = { real: z.real * z.real - z.imag * z.imag, imag: 2 * z.real * z.imag };
                const num_bessel = { real: b0_bessel + b1_bessel * z.real + b2_bessel * z2_bessel.real, imag: b1_bessel * z.imag + b2_bessel * z2_bessel.imag };
                const denom_bessel = { real: a0_bessel + a1_bessel * z.real + a2_bessel * z2_bessel.real, imag: a1_bessel * z.imag + a2_bessel * z2_bessel.imag };
                const denomMag_bessel = denom_bessel.real * denom_bessel.real + denom_bessel.imag * denom_bessel.imag;
                H.real = (num_bessel.real * denom_bessel.real + num_bessel.imag * denom_bessel.imag) / denomMag_bessel;
                H.imag = (num_bessel.imag * denom_bessel.real - num_bessel.real * denom_bessel.imag) / denomMag_bessel;
                break;
            case 'NOTCH':
                const notchOmega = 2 * Math.PI * cutoff / sampleRate;
                const notchAlpha = Math.sin(notchOmega) / (2 * Q);
                const b0_notch = 1;
                const b1_notch = -2 * Math.cos(notchOmega);
                const b2_notch = 1;
                const a1_notch = b1_notch;
                const a2_notch = 1 - notchAlpha;
                const a0_notch = 1 + notchAlpha;
                
                const z2_notch = { real: z.real * z.real - z.imag * z.imag, imag: 2 * z.real * z.imag };
                const num_notch = { real: b0_notch + b1_notch * z.real + b2_notch * z2_notch.real, imag: b1_notch * z.imag + b2_notch * z2_notch.imag };
                const denom_notch = { real: a0_notch + a1_notch * z.real + a2_notch * z2_notch.real, imag: a1_notch * z.imag + a2_notch * z2_notch.imag };
                const denomMag_notch = denom_notch.real * denom_notch.real + denom_notch.imag * denom_notch.imag;
                H.real = (num_notch.real * denom_notch.real + num_notch.imag * denom_notch.imag) / denomMag_notch;
                H.imag = (num_notch.imag * denom_notch.real - num_notch.real * denom_notch.imag) / denomMag_notch;
                break;
        }
        
        frequencies.push(f);
        magnitudes.push(20 * Math.log10(Math.sqrt(H.real * H.real + H.imag * H.imag)));
        phases.push(Math.atan2(H.imag, H.real) * 180 / Math.PI);
    }
    
    return { frequencies, magnitudes, phases };
}

// Complex number helper
const complexI = { real: 0, imag: 1 };

// Signal generation functions
function generateWhiteNoise(length, amplitude = 1) {
    const signal = [];
    for (let i = 0; i < length; i++) {
        signal.push((Math.random() - 0.5) * 2 * amplitude);
    }
    return signal;
}

function generateSineWave(length, frequency, sampleRate, amplitude = 1) {
    const signal = [];
    for (let i = 0; i < length; i++) {
        signal.push(amplitude * Math.sin(2 * Math.PI * frequency * i / sampleRate));
    }
    return signal;
}

function generateStepResponse(length, stepTime = 0.1) {
    const signal = [];
    const stepIndex = Math.floor(stepTime * sampleRate);
    for (let i = 0; i < length; i++) {
        signal.push(i >= stepIndex ? 1 : 0);
    }
    return signal;
}

function generateChirp(length, startFreq, endFreq, sampleRate) {
    const signal = [];
    for (let i = 0; i < length; i++) {
        const t = i / sampleRate;
        const freq = startFreq + (endFreq - startFreq) * t / (length / sampleRate);
        signal.push(Math.sin(2 * Math.PI * freq * t));
    }
    return signal;
}

function generateRealisticGyroData(length, sampleRate) {
    const signal = [];
    const baseFreq = 50;
    const rotorFreq = 120;
    const noiseLevel = 0.1;
    
    for (let i = 0; i < length; i++) {
        const t = i / sampleRate;
        const baseSignal = Math.sin(2 * Math.PI * baseFreq * t);
        const rotorSignal = 0.3 * Math.sin(2 * Math.PI * rotorFreq * t);
        const noise = (Math.random() - 0.5) * 2 * noiseLevel;
        signal.push(baseSignal + rotorSignal + noise);
    }
    return signal;
}

// Plotting functions
function createFrequencyResponsePlot(containerId, data, title) {
    // Clear existing plot first
    const container = document.getElementById(containerId);
    if (container) {
        container.innerHTML = '';
    }
    
    const trace1 = {
        x: data.frequencies,
        y: data.magnitudes,
        type: 'scatter',
        mode: 'lines',
        name: 'Magnitude',
        line: { color: '#667eea', width: 2 }
    };
    
    const trace2 = {
        x: data.frequencies,
        y: data.phases,
        type: 'scatter',
        mode: 'lines',
        name: 'Phase',
        yaxis: 'y2',
        line: { color: '#e74c3c', width: 2 }
    };
    
    const layout = {
        title: title,
        xaxis: { title: 'Frequency (Hz)', type: 'log' },
        yaxis: { title: 'Magnitude (dB)', side: 'left' },
        yaxis2: { title: 'Phase (degrees)', side: 'right', overlaying: 'y' },
        showlegend: true,
        legend: { x: 0.1, y: 0.9 },
        margin: { l: 50, r: 50, t: 50, b: 50 },
        autosize: true,
        height: 350
    };
    
    const config = {
        responsive: true,
        displayModeBar: false
    };
    
    Plotly.newPlot(containerId, [trace1, trace2], layout, config);
}

function createTimeResponsePlot(containerId, time, input, output, title) {
    // Clear existing plot first
    const container = document.getElementById(containerId);
    if (container) {
        container.innerHTML = '';
    }
    
    const trace1 = {
        x: time,
        y: input,
        type: 'scatter',
        mode: 'lines',
        name: 'Input',
        line: { color: '#95a5a6', width: 1 }
    };
    
    const trace2 = {
        x: time,
        y: output,
        type: 'scatter',
        mode: 'lines',
        name: 'Output',
        line: { color: '#667eea', width: 2 }
    };
    
    const layout = {
        title: title,
        xaxis: { title: 'Time (s)' },
        yaxis: { title: 'Amplitude' },
        showlegend: true,
        legend: { x: 0.1, y: 0.9 },
        margin: { l: 50, r: 50, t: 50, b: 50 },
        autosize: true,
        height: 350
    };
    
    const config = {
        responsive: true,
        displayModeBar: false
    };
    
    Plotly.newPlot(containerId, [trace1, trace2], layout, config);
}

// Filter update functions
function updateDecimationPlots() {
    const cutoff = parseFloat(document.getElementById('decimation-cutoff').value);
    const sampleRate = parseFloat(document.getElementById('decimation-sample-rate').value);
    
    // Calculate Bessel 4th order response
    const bessel1 = calculateFrequencyResponse('BESSEL', BESSEL_4A_C * cutoff, sampleRate, BESSEL_4A_Q);
    const bessel2 = calculateFrequencyResponse('BESSEL', BESSEL_4B_C * cutoff, sampleRate, BESSEL_4B_Q);
    
    // Combine responses (cascade)
    const combinedMagnitudes = bessel1.magnitudes.map((mag, i) => mag + bessel2.magnitudes[i]);
    const combinedPhases = bessel1.phases.map((phase, i) => phase + bessel2.phases[i]);
    
    const freqData = {
        frequencies: bessel1.frequencies,
        magnitudes: combinedMagnitudes,
        phases: combinedPhases
    };
    
    createFrequencyResponsePlot('decimation-freq-plot', freqData, 'Decimation Filter Frequency Response');
    
    // Step response
    const time = Array.from({length: 1000}, (_, i) => i / sampleRate);
    const stepInput = generateStepResponse(1000, 0.1);
    
    const filter1 = new BiquadFilter();
    const filter2 = new BiquadFilter();
    filter1.init(BESSEL_4A_C * cutoff, sampleRate, BESSEL_4A_Q, 'LPF');
    filter2.init(BESSEL_4B_C * cutoff, sampleRate, BESSEL_4B_Q, 'LPF');
    
    const stepOutput = stepInput.map(input => {
        const temp = filter1.apply(input);
        return filter2.apply(temp);
    });
    
    createTimeResponsePlot('decimation-step-plot', time, stepInput, stepOutput, 'Decimation Filter Step Response');
}

function updateRpmPlots() {
    const rpm = parseFloat(document.getElementById('rpm-motor-rpm').value);
    const ratio = parseFloat(document.getElementById('rpm-ratio').value);
    const q = parseFloat(document.getElementById('rpm-q').value);
    
    const frequency = rpm * ratio / 60;
    const sampleRate = 4000;
    
    const freqData = calculateFrequencyResponse('NOTCH', frequency, sampleRate, q);
    createFrequencyResponsePlot('rpm-freq-plot', freqData, `RPM Filter Frequency Response (${frequency.toFixed(1)} Hz)`);
    
    // Time response with rotor vibration
    const time = Array.from({length: 1000}, (_, i) => i / sampleRate);
    const rotorSignal = generateSineWave(1000, frequency, sampleRate, 1);
    const noise = generateWhiteNoise(1000, 0.2);
    const input = rotorSignal.map((val, i) => val + noise[i]);
    
    const filter = new BiquadFilter();
    filter.init(frequency, sampleRate, q, 'NOTCH');
    const output = input.map(val => filter.apply(val));
    
    createTimeResponsePlot('rpm-time-plot', time, input, output, 'RPM Filter Time Response');
}

function updateLowpassPlots() {
    const lpf1Type = document.getElementById('lpf1-type').value;
    const lpf1Cutoff = parseFloat(document.getElementById('lpf1-cutoff').value);
    const lpf2Type = document.getElementById('lpf2-type').value;
    const lpf2Cutoff = parseFloat(document.getElementById('lpf2-cutoff').value);
    const sampleRate = 4000;
    
    // LPF1 plots
    const lpf1Data = calculateFrequencyResponse(lpf1Type, lpf1Cutoff, sampleRate);
    createFrequencyResponsePlot('lpf1-freq-plot', lpf1Data, `LPF1 Frequency Response (${lpf1Type})`);
    
    const time = Array.from({length: 1000}, (_, i) => i / sampleRate);
    const stepInput = generateStepResponse(1000, 0.1);
    
    let filter1;
    switch (lpf1Type) {
        case 'PT1':
            filter1 = new PT1Filter();
            break;
        case 'PT2':
            filter1 = new PT2Filter();
            break;
        case 'PT3':
            filter1 = new PT3Filter();
            break;
        default:
            filter1 = new BiquadFilter();
            filter1.init(lpf1Cutoff, sampleRate, 1, 'LPF');
    }
    
    if (lpf1Type.startsWith('PT')) {
        filter1.init(lpf1Cutoff, sampleRate);
    }
    
    const lpf1Output = stepInput.map(input => filter1.apply(input));
    createTimeResponsePlot('lpf1-step-plot', time, stepInput, lpf1Output, `LPF1 Step Response (${lpf1Type})`);
    
    // LPF2 plots
    const lpf2Data = calculateFrequencyResponse(lpf2Type, lpf2Cutoff, sampleRate);
    createFrequencyResponsePlot('lpf2-freq-plot', lpf2Data, `LPF2 Frequency Response (${lpf2Type})`);
    
    let filter2;
    switch (lpf2Type) {
        case 'PT1':
            filter2 = new PT1Filter();
            break;
        case 'PT2':
            filter2 = new PT2Filter();
            break;
        case 'PT3':
            filter2 = new PT3Filter();
            break;
        default:
            filter2 = new BiquadFilter();
            filter2.init(lpf2Cutoff, sampleRate, 1, 'LPF');
    }
    
    if (lpf2Type.startsWith('PT')) {
        filter2.init(lpf2Cutoff, sampleRate);
    }
    
    const lpf2Output = stepInput.map(input => filter2.apply(input));
    createTimeResponsePlot('lpf2-step-plot', time, stepInput, lpf2Output, `LPF2 Step Response (${lpf2Type})`);
    
    // Comparison plots
    const comparisonData = {
        frequencies: lpf1Data.frequencies,
        magnitudes: lpf1Data.magnitudes.map((mag, i) => mag + lpf2Data.magnitudes[i]),
        phases: lpf1Data.phases.map((phase, i) => phase + lpf2Data.phases[i])
    };
    
    createFrequencyResponsePlot('lpf-comparison-freq-plot', comparisonData, 'Cascaded LPF Frequency Response');
    
    const cascadedOutput = stepInput.map(input => {
        const temp = filter1.apply(input);
        return filter2.apply(temp);
    });
    
    createTimeResponsePlot('lpf-comparison-step-plot', time, stepInput, cascadedOutput, 'Cascaded LPF Step Response');
}

function updateNotchPlots() {
    const notch1Center = parseFloat(document.getElementById('notch1-center').value);
    const notch1Cutoff = parseFloat(document.getElementById('notch1-cutoff').value);
    const notch2Center = parseFloat(document.getElementById('notch2-center').value);
    const notch2Cutoff = parseFloat(document.getElementById('notch2-cutoff').value);
    const sampleRate = 4000;
    
    // Calculate Q values
    const q1 = notch1Center * notch1Cutoff / (notch1Center * notch1Center - notch1Cutoff * notch1Cutoff);
    const q2 = notch2Center * notch2Cutoff / (notch2Center * notch2Center - notch2Cutoff * notch2Cutoff);
    
    // Notch 1 plots
    const notch1Data = calculateFrequencyResponse('NOTCH', notch1Center, sampleRate, q1);
    createFrequencyResponsePlot('notch1-freq-plot', notch1Data, `Notch 1 Frequency Response (${notch1Center} Hz)`);
    
    const time = Array.from({length: 1000}, (_, i) => i / sampleRate);
    const sineInput = generateSineWave(1000, notch1Center, sampleRate, 1);
    const filter1 = new BiquadFilter();
    filter1.init(notch1Center, sampleRate, q1, 'NOTCH');
    const notch1Output = sineInput.map(input => filter1.apply(input));
    createTimeResponsePlot('notch1-time-plot', time, sineInput, notch1Output, 'Notch 1 Time Response');
    
    // Notch 2 plots
    const notch2Data = calculateFrequencyResponse('NOTCH', notch2Center, sampleRate, q2);
    createFrequencyResponsePlot('notch2-freq-plot', notch2Data, `Notch 2 Frequency Response (${notch2Center} Hz)`);
    
    const sineInput2 = generateSineWave(1000, notch2Center, sampleRate, 1);
    const filter2 = new BiquadFilter();
    filter2.init(notch2Center, sampleRate, q2, 'NOTCH');
    const notch2Output = sineInput2.map(input => filter2.apply(input));
    createTimeResponsePlot('notch2-time-plot', time, sineInput2, notch2Output, 'Notch 2 Time Response');
    
    // Cascaded response
    const cascadedData = {
        frequencies: notch1Data.frequencies,
        magnitudes: notch1Data.magnitudes.map((mag, i) => mag + notch2Data.magnitudes[i]),
        phases: notch1Data.phases.map((phase, i) => phase + notch2Data.phases[i])
    };
    
    createFrequencyResponsePlot('notch-cascaded-freq-plot', cascadedData, 'Cascaded Notch Frequency Response');
    
    const cascadedOutput = sineInput.map(input => {
        const temp = filter1.apply(input);
        return filter2.apply(temp);
    });
    
    createTimeResponsePlot('notch-cascaded-time-plot', time, sineInput, cascadedOutput, 'Cascaded Notch Time Response');
}

function updateDynNotchPlots() {
    const notchCount = parseInt(document.getElementById('dyn-notch-count').value);
    const q = parseFloat(document.getElementById('dyn-notch-q').value);
    const minFreq = parseFloat(document.getElementById('dyn-notch-min').value);
    const maxFreq = parseFloat(document.getElementById('dyn-notch-max').value);
    const sampleRate = 4000;
    
    // Generate FFT spectrum with peaks
    const frequencies = Array.from({length: 200}, (_, i) => i * sampleRate / 400);
    const spectrum = frequencies.map(f => {
        let magnitude = 0.1; // Noise floor
        if (f >= minFreq && f <= maxFreq) {
            // Add some peaks
            const peak1 = 80 + 20 * Math.sin(f / 50);
            const peak2 = 60 + 15 * Math.sin(f / 30);
            magnitude = Math.max(magnitude, peak1, peak2);
        }
        return magnitude + (Math.random() - 0.5) * 5;
    });
    
    // Find peaks
    const peaks = [];
    for (let i = 1; i < spectrum.length - 1; i++) {
        if (spectrum[i] > spectrum[i-1] && spectrum[i] > spectrum[i+1] && 
            frequencies[i] >= minFreq && frequencies[i] <= maxFreq) {
            peaks.push({ freq: frequencies[i], magnitude: spectrum[i] });
        }
    }
    peaks.sort((a, b) => b.magnitude - a.magnitude);
    const topPeaks = peaks.slice(0, notchCount);
    
    // FFT plot
    const trace1 = {
        x: frequencies,
        y: spectrum,
        type: 'scatter',
        mode: 'lines',
        name: 'FFT Spectrum',
        line: { color: '#667eea', width: 1 }
    };
    
    const trace2 = {
        x: topPeaks.map(p => p.freq),
        y: topPeaks.map(p => p.magnitude),
        type: 'scatter',
        mode: 'markers',
        name: 'Detected Peaks',
        marker: { color: '#e74c3c', size: 8 }
    };
    
    const layout = {
        title: 'FFT Spectrum & Detected Peaks',
        xaxis: { title: 'Frequency (Hz)' },
        yaxis: { title: 'Magnitude' },
        showlegend: true,
        legend: { x: 0.1, y: 0.9 },
        margin: { l: 50, r: 50, t: 50, b: 50 }
    };
    
    Plotly.newPlot('dyn-notch-fft-plot', [trace1, trace2], layout);
    
    // Dynamic notch frequency response
    let combinedMagnitudes = new Array(frequencies.length).fill(0);
    let combinedPhases = new Array(frequencies.length).fill(0);
    
    topPeaks.forEach(peak => {
        const notchData = calculateFrequencyResponse('NOTCH', peak.freq, sampleRate, q);
        combinedMagnitudes = combinedMagnitudes.map((mag, i) => mag + notchData.magnitudes[i]);
        combinedPhases = combinedPhases.map((phase, i) => phase + notchData.phases[i]);
    });
    
    const dynNotchData = {
        frequencies: frequencies,
        magnitudes: combinedMagnitudes,
        phases: combinedPhases
    };
    
    createFrequencyResponsePlot('dyn-notch-freq-plot', dynNotchData, 'Dynamic Notch Frequency Response');
}

function updatePipelinePlots() {
    const signalFreq = parseFloat(document.getElementById('pipeline-signal-freq').value);
    const noiseLevel = parseFloat(document.getElementById('pipeline-noise-level').value) / 100;
    const sampleRate = 4000;
    
    // Generate input signal
    let inputSignal;
    switch (currentInputSignal) {
        case 'noise':
            inputSignal = generateWhiteNoise(1000, 1);
            break;
        case 'step':
            inputSignal = generateStepResponse(1000, 0.1);
            break;
        case 'sine':
            inputSignal = generateSineWave(1000, signalFreq, sampleRate, 1);
            break;
        case 'chirp':
            inputSignal = generateChirp(1000, 10, 500, sampleRate);
            break;
        case 'realistic':
            inputSignal = generateRealisticGyroData(1000, sampleRate);
            break;
    }
    
    // Add noise if needed
    if (noiseLevel > 0) {
        const noise = generateWhiteNoise(1000, noiseLevel);
        inputSignal = inputSignal.map((val, i) => val + noise[i]);
    }
    
    // Apply complete pipeline
    const time = Array.from({length: 1000}, (_, i) => i / sampleRate);
    
    // Initialize all filters
    const decimFilter1 = new BiquadFilter();
    const decimFilter2 = new BiquadFilter();
    decimFilter1.init(BESSEL_4A_C * 200, sampleRate, BESSEL_4A_Q, 'LPF');
    decimFilter2.init(BESSEL_4B_C * 200, sampleRate, BESSEL_4B_Q, 'LPF');
    
    const rpmFilter = new BiquadFilter();
    rpmFilter.init(120, sampleRate, 20, 'NOTCH');
    
    const lpf1 = new PT1Filter();
    lpf1.init(200, sampleRate);
    
    const lpf2 = new PT1Filter();
    lpf2.init(150, sampleRate);
    
    const notch1 = new BiquadFilter();
    notch1.init(150, sampleRate, 5, 'NOTCH');
    
    const notch2 = new BiquadFilter();
    notch2.init(250, sampleRate, 8, 'NOTCH');
    
    // Apply pipeline
    const output = inputSignal.map(input => {
        let signal = input;
        
        // Decimation
        signal = decimFilter1.apply(signal);
        signal = decimFilter2.apply(signal);
        
        // RPM filter
        signal = rpmFilter.apply(signal);
        
        // LPF1
        signal = lpf1.apply(signal);
        
        // LPF2
        signal = lpf2.apply(signal);
        
        // Notch1
        signal = notch1.apply(signal);
        
        // Notch2
        signal = notch2.apply(signal);
        
        return signal;
    });
    
    // Time domain plot
    createTimeResponsePlot('pipeline-time-plot', time, inputSignal, output, 'Pipeline Time Domain Response');
    
    // Frequency domain plot (FFT)
    const inputFFT = calculateFFT(inputSignal, sampleRate);
    const outputFFT = calculateFFT(output, sampleRate);
    
    const trace1 = {
        x: inputFFT.frequencies,
        y: inputFFT.magnitudes,
        type: 'scatter',
        mode: 'lines',
        name: 'Input Spectrum',
        line: { color: '#95a5a6', width: 1 }
    };
    
    const trace2 = {
        x: outputFFT.frequencies,
        y: outputFFT.magnitudes,
        type: 'scatter',
        mode: 'lines',
        name: 'Output Spectrum',
        line: { color: '#667eea', width: 2 }
    };
    
    const layout = {
        title: 'Pipeline Frequency Response',
        xaxis: { title: 'Frequency (Hz)' },
        yaxis: { title: 'Magnitude (dB)' },
        showlegend: true,
        legend: { x: 0.1, y: 0.9 },
        margin: { l: 50, r: 50, t: 50, b: 50 }
    };
    
    Plotly.newPlot('pipeline-freq-plot', [trace1, trace2], layout);
    
    // Before vs After comparison
    createTimeResponsePlot('pipeline-comparison-plot', time, inputSignal, output, 'Before vs After Filtering');
    
    // Filter stage contributions
    const stages = ['Raw', 'Decimation', 'RPM', 'LPF1', 'LPF2', 'Notch1', 'Notch2', 'Final'];
    const stageValues = [1, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1, 0.05];
    
    const trace3 = {
        x: stages,
        y: stageValues,
        type: 'bar',
        name: 'Noise Reduction',
        marker: { color: '#667eea' }
    };
    
    const layout2 = {
        title: 'Filter Stage Contributions',
        xaxis: { title: 'Filter Stage' },
        yaxis: { title: 'Relative Noise Level' },
        margin: { l: 50, r: 50, t: 50, b: 50 }
    };
    
    Plotly.newPlot('pipeline-stages-plot', [trace3], layout2);
}

// FFT calculation helper
function calculateFFT(signal, sampleRate) {
    const n = signal.length;
    const frequencies = Array.from({length: n/2}, (_, i) => i * sampleRate / n);
    const magnitudes = [];
    
    for (let k = 0; k < n/2; k++) {
        let real = 0;
        let imag = 0;
        
        for (let i = 0; i < n; i++) {
            const angle = -2 * Math.PI * k * i / n;
            real += signal[i] * Math.cos(angle);
            imag += signal[i] * Math.sin(angle);
        }
        
        const magnitude = Math.sqrt(real * real + imag * imag) / n;
        magnitudes.push(20 * Math.log10(magnitude + 1e-10));
    }
    
    return { frequencies, magnitudes };
}
