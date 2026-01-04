package funkin.vis.dsp;

import lime.media.AudioBuffer;
import lime.utils.UInt8Array;
import funkin.vis.dsp.Complex;
import funkin.vis.dsp.FFT;
import funkin.vis.dsp.Signal;
import funkin.vis.dsp.VectorMath;
import haxe.ds.Vector;

using grig.audio.lime.UInt8ArrayTools;

typedef BPMResult = {
    var bpm:Float;
    var confidence:Float;
    var beats:Array<Float>;
}

typedef BPMEvent = {
    var time:Float;
    var bpm:Float;
    var confidence:Float;
}

class MusicBPMDetector {
    public var sampleRate:Float;
    
    // Config for analysis
    private static inline var WINDOW_SIZE = 1024;
    private static inline var HOP_SIZE = 256;
    private static inline var MIN_BPM = 60.0;
    private static inline var MAX_BPM = 180.0;
    
    // Last analysis data
    private var lastEnvelope:Array<Float>;
    private var lastCorrelation:Array<Float>;

    public function new(sampleRate:Float = 44100) {
        this.sampleRate = sampleRate;
        this.lastEnvelope = [];
        this.lastCorrelation = [];
    }

    /**
     * Executes BPM detection and returns the result.
     * sound._sound.__buffer for flixel
     */
    public function detect(buffer:AudioBuffer):BPMResult {
        // ... (existing implementation)
        return _detectInternal(buffer);
    }
    
    private function _detectInternal(buffer:AudioBuffer):BPMResult {
        if (buffer == null || buffer.data == null) {
            throw "Invalid audio buffer";
        }

        trace('[MusicBPMDetector] Buffer Info: Length=${buffer.data.length}, SR=${buffer.sampleRate}, Bits=${buffer.bitsPerSample}, Ch=${buffer.channels}');
        
        var duration = (buffer.data.length / (buffer.bitsPerSample / 8 * buffer.channels)) / buffer.sampleRate;
        trace('[MusicBPMDetector] Audio Duration: ${duration}s');

        if (duration < 5.0) {
            trace('[MusicBPMDetector] WARNING: Audio is too short (< 5s). Is this a streaming sound? BPM detection requires the full audio file.');
        }

        // 1. Convert to Mono Float Signal
        var signal = getSignal(buffer);
        if (signal.length == 0) {
            throw "Empty signal";
        }
        
        // DEBUG: Check signal quality
        var maxSig = 0.0;
        var minSig = 0.0;
        for(s in signal) {
            if(s > maxSig) maxSig = s;
            if(s < minSig) minSig = s;
        }
        trace('[MusicBPMDetector] Signal Range: [${minSig}, ${maxSig}]');
        if (maxSig == 0 && minSig == 0) {
            trace('[MusicBPMDetector] ERROR: Signal is completely silent! Check getSignal implementation.');
            return {bpm: 0, confidence: 0, beats: []};
        }
        
        // 2. Compute Onset Envelope (Flux)
        var envelope = computeOnsetEnvelope(signal, buffer.sampleRate);
        
        // DEBUG: Check envelope
        var maxEnv = 0.0;
        for(e in envelope) if(e > maxEnv) maxEnv = e;
        trace('[MusicBPMDetector] Max Envelope Flux: ${maxEnv}');

        var envelopeRate = buffer.sampleRate / HOP_SIZE;
        this.lastEnvelope = envelope;
        
        // 3. Autocorrelation
        var correlation = computeAutocorrelation(envelope);
        this.lastCorrelation = correlation;

        // DEBUG: Check correlation
        var maxCorr = 0.0;
        for(c in correlation) if(c > maxCorr) maxCorr = c;
        trace('[MusicBPMDetector] Max Autocorrelation: ${maxCorr}');
        
        // 4. Find BPM Peak
        var bestBPM = findBestBPM(correlation, envelopeRate);
        trace('[MusicBPMDetector] Best BPM Raw: ${bestBPM.bpm}, Conf: ${bestBPM.confidence}');
        
        // 5. Detect Beats (Phase)
        var beats = detectBeats(envelope, envelopeRate, bestBPM.bpm);

        return {
            bpm: bestBPM.bpm,
            confidence: bestBPM.confidence,
            beats: beats
        };
    }

    /**
     * Detects variable BPM by analyzing the audio in segments.
     * @param buffer Audio buffer
     * @param windowSeconds Duration of each analysis window (default 10s)
     * @param intervalSeconds Step size for the sliding window (default 5s)
     * @return Array of BPM events sorted by time
     */
    public function detectTimeVarying(buffer:AudioBuffer, windowSeconds:Float = 10.0, intervalSeconds:Float = 5.0):Array<BPMEvent> {
        if (buffer == null || buffer.data == null) throw "Invalid audio buffer";
        
        var signal = getSignal(buffer);
        var envelope = computeOnsetEnvelope(signal, buffer.sampleRate);
        var envelopeRate = buffer.sampleRate / HOP_SIZE;
        
        var windowSamples = Std.int(windowSeconds * envelopeRate);
        var stepSamples = Std.int(intervalSeconds * envelopeRate);
        
        var results = new Array<BPMEvent>();
        
        var i = 0;
        while (i + windowSamples <= envelope.length) {
            var segment = envelope.slice(i, i + windowSamples);
            var correlation = computeAutocorrelation(segment);
            var best = findBestBPM(correlation, envelopeRate);
            
            // Time is the center of the window
            var time = (i + windowSamples / 2) / envelopeRate;
            
            results.push({
                time: time,
                bpm: best.bpm,
                confidence: best.confidence
            });
            
            i += stepSamples;
        }
        
        return results;
    }

    public function getRawAnalysis():Dynamic {
        return {
            envelope: lastEnvelope,
            correlation: lastCorrelation
        };
    }

    // --- Helpers ---

    private function getSignal(buffer:AudioBuffer):Array<Float> {
        var data = buffer.data;
        var channels = buffer.channels;
        var bits = buffer.bitsPerSample;
        var length = data.length;
        var samples = new Array<Float>();
        
        var bytesPerSample = Std.int(bits / 8);
        var frameCount = Std.int(length / (bytesPerSample * channels));
        
        samples.resize(frameCount);
        for (i in 0...frameCount) {
            var sum = 0.0;
            for (c in 0...channels) {
                var offset = (i * channels + c) * bytesPerSample;
                var val = 0.0;
                switch (bits) {
                    case 8:
                        val = (data[offset] - 128) / 128.0;
                    case 16:
                        val = data.getInt16(offset) / 32768.0;
                    case 24:
                        val = data.getInt24(offset) / 8388608.0; 
                    case 32:
                        val = data.getInt32(offset) / 2147483648.0;
                }
                sum += val;
            }
            samples[i] = sum / channels;
        }
        return samples;
    }

    private function computeOnsetEnvelope(signal:Array<Float>, sr:Int):Array<Float> {
        var result = new Array<Float>();
        var numFrames = Std.int((signal.length - WINDOW_SIZE) / HOP_SIZE);
        if (numFrames < 0) return [];
        
        result.resize(numFrames);
        
        var prevEnergy = 0.0;
        
        for (i in 0...numFrames) {
            var start = i * HOP_SIZE;
            var energy = 0.0;
            
            var j = 0;
            while(j < WINDOW_SIZE - 3) {
                var v0 = signal[start + j];
                var v1 = signal[start + j + 1];
                var v2 = signal[start + j + 2];
                var v3 = signal[start + j + 3];
                energy += v0 * v0 + v1 * v1 + v2 * v2 + v3 * v3;
                j += 4;
            }
            while(j < WINDOW_SIZE) {
                var val = signal[start + j];
                energy += val * val;
                j++;
            }
            
            var flux = energy - prevEnergy;
            if (flux < 0) flux = 0;
            
            result[i] = flux;
            prevEnergy = energy;
        }
        
        return smooth(result, 5);
    }
    
    private function smooth(data:Array<Float>, window:Int):Array<Float> {
        var output = new Array<Float>();
        output.resize(data.length);
        // Simple moving average
        for (i in 0...data.length) {
            var sum = 0.0;
            var count = 0;
            var start = i - Std.int(window/2);
            var end = i + Std.int(window/2) + 1;
            if (start < 0) start = 0;
            if (end > data.length) end = data.length;
            
            for (j in start...end) {
                sum += data[j];
                count++;
            }
            output[i] = count > 0 ? sum / count : 0;
        }
        return output;
    }

    private function computeAutocorrelation(signal:Array<Float>):Array<Float> {
        var n = FFT.nextPow2(signal.length * 2);
        
        // Use Vectors for SoA FFT (Structure of Arrays)
        // This avoids creating thousands of Complex objects per frame
        var real = new Vector<Float>(n);
        var imag = new Vector<Float>(n); 
        
        // Initialize real part with signal
        for (i in 0...n) {
            if (i < signal.length) {
                real[i] = signal[i];
            } else {
                real[i] = 0.0;
            }
            imag[i] = 0.0;
        }
        
        // Forward FFT
        FFT.fftVectors(real, imag, false);
        
        // DEBUG: Check Spectrum
        var maxSpec = 0.0;
        for(i in 0...n) {
            var mag = real[i]*real[i] + imag[i]*imag[i];
            if(mag > maxSpec) maxSpec = mag;
        }
        trace('[MusicBPMDetector] Max Spectrum MagSq: ${maxSpec}');

        // Compute Power Spectrum: S = F * conj(F)
        // |a+bi|^2 = a^2 + b^2. Imaginary part becomes 0.
        VectorMath.complexMagSq(real, imag, real);

        // Only need real part for power spectrum, imag is 0
        for(i in 0...n) imag[i] = 0.0;
        
        // Inverse FFT to get autocorrelation
        FFT.fftVectors(real, imag, true);
        
        // Extract result
        var autocorr = new Array<Float>();
        autocorr.resize(signal.length);
        
        // For real input, autocorrelation is real, so we take real part of IFFT
        for (i in 0...signal.length) {
            autocorr[i] = real[i];
        }
        
        if (autocorr.length > 0 && autocorr[0] != 0) {
            var maxVal = autocorr[0];
            var invMax = 1.0 / maxVal;
            for(i in 0...autocorr.length) autocorr[i] *= invMax;
        }
        
        return autocorr;
    }

    private function findBestBPM(autocorr:Array<Float>, sr:Float):{bpm:Float, confidence:Float} {
        var minLag = Std.int(sr * 60 / MAX_BPM);
        var maxLag = Std.int(sr * 60 / MIN_BPM);
        
        var maxVal = -1.0;
        var bestLag = -1;
        
        if (minLag < 0) minLag = 0;
        if (maxLag >= autocorr.length) maxLag = autocorr.length - 1;
        
        for (i in minLag...maxLag + 1) {
            if (autocorr[i] > maxVal) {
                maxVal = autocorr[i];
                bestLag = i;
            }
        }
        
        var bpm = 0.0;
        if (bestLag > 0) {
            bpm = 60.0 * sr / bestLag;
        }
        
        var sum = 0.0;
        var count = 0;
        for (i in minLag...maxLag + 1) {
            sum += autocorr[i];
            count++;
        }
        var avg = count > 0 ? sum / count : 0;
        var confidence = (avg > 0) ? (maxVal / avg - 1.0) : 0.0;
        confidence = Math.min(confidence / 10.0, 1.0); 
        if (confidence < 0) confidence = 0;

        return {bpm: Math.round(bpm * 10) / 10.0, confidence: confidence};
    }
    
    private function detectBeats(envelope:Array<Float>, sr:Float, bpm:Float):Array<Float> {
        if (bpm <= 0) return [];
        
        var beatPeriodSamples = 60.0 * sr / bpm;
        var beats = new Array<Float>();
        
        var searchWindow = Std.int(beatPeriodSamples * 4);
        if (searchWindow > envelope.length) searchWindow = envelope.length;
        
        var maxVal = -1.0;
        var phase = 0;
        
        for (i in 0...searchWindow) {
            if (envelope[i] > maxVal) {
                maxVal = envelope[i];
                phase = i;
            }
        }
        
        var t = phase * 1.0;
        while (t < envelope.length) {
            beats.push(t / sr); 
            t += beatPeriodSamples;
        }
        
        return beats;
    }
}
