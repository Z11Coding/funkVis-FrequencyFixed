package funkin.vis.dsp;

import funkin.vis.dsp.Complex;
import haxe.ds.Vector;
import funkin.vis.dsp.VectorMath;

// these are only used for testing, down in FFT.main()
using funkin.vis.dsp.OffsetArray;
using funkin.vis.dsp.Signal;


/**
	Fast/Finite Fourier Transforms.
**/
class FFT {
	/**
		Computes the Discrete Fourier Transform (DFT) of a `Complex` sequence.

		If the input has N data points (N should be a power of 2 or padding will be added)
		from a signal sampled at intervals of 1/Fs, the result will be a sequence of N
		samples from the Discrete-Time Fourier Transform (DTFT) - which is Fs-periodic -
		with a spacing of Fs/N Hz between them and a scaling factor of Fs.
	**/
	public static function fft(input:Array<Complex>) : Array<Complex>
	{
		// Compatibility wrapper
		var n = nextPow2(input.length);
		var real = new Vector<Float>(n);
		var imag = new Vector<Float>(n);
		
		for(i in 0...input.length) {
			real[i] = input[i].real;
			imag[i] = input[i].imag;
		}

		for(i in input.length...n) {
			real[i] = 0.0;
			imag[i] = 0.0;
		}

		fftVectors(real, imag, false);

		var res = new Array<Complex>();
		res.resize(n);
		for(i in 0...n) {
			res[i] = new Complex(real[i], imag[i]);
		}
		return res;
	}

	/**
	 * Optimized FFT taking separate real and imaginary vectors (Structure of Arrays).
	 * Performs calculation IN-PLACE.
	 * @param real Real parts
	 * @param imag Imaginary parts
	 * @param inverse If true, performs Inverse FFT
	 */
	public static function fftVectors(real:Vector<Float>, imag:Vector<Float>, inverse:Bool) : Void {
		var n = real.length;
		if (imag.length != n) throw "Real and Imag vectors must be same length";

		if (inverse && twiddleRealInv?.length != n)
			precomputeTwiddleFactors(n, true);
		else if (!inverse && twiddleReal?.length != n)
			precomputeTwiddleFactors(n, false);

        // Bit-reversal permutation
        var j = 0;
        for (i in 0...n - 1) {
            if (i < j) {
                // Swap real[i] and real[j]
                var tr = real[i]; real[i] = real[j]; real[j] = tr;
                // Swap imag[i] and imag[j]
                var ti = imag[i]; imag[i] = imag[j]; imag[j] = ti;
            }
            var m = n >> 1;
            while (m >= 2 && j >= m) {
                j -= m;
                m >>= 1;
            }
            j += m;
        }

        // Iterative Radix-2 FFT
        var tr = inverse ? twiddleRealInv : twiddleReal;
        var ti = inverse ? twiddleImagInv : twiddleImag;
        
        var step = 1;
        while (step < n) {
            var jump = step << 1;
            // The stride for accessing twiddles.
            // When step=1 (N=2), we need W_N^0. stride = N/2.
            // When step=2 (N=4), we need W_N^0, W_N^(N/4). stride = N/4.
            // Twiddle index increments by stride.
            var stride = Std.int(tr.length / jump);
            
            var group = 0;
            while (group < n) {
                for (pair in 0...step) {
                    var idx = pair * stride;
                    var wr = tr[idx];
                    var wi = ti[idx];
                    
                    var match = group + pair + step;
                    var i = group + pair;
                    
                    // t = w * data[match]
                    var d_mr = real[match];
                    var d_mi = imag[match];
                    var tr_val = wr * d_mr - wi * d_mi;
                    var ti_val = wr * d_mi + wi * d_mr;
                    
                    real[match] = real[i] - tr_val;
                    imag[match] = imag[i] - ti_val;
                    real[i] += tr_val;
                    imag[i] += ti_val;
                }
                group += jump;
            }
            step <<= 1;
        }
		
		if (inverse) {
			VectorMath.mulScalar(real, 1.0/n, real);
			VectorMath.mulScalar(imag, 1.0/n, imag);
		}
	}

	/**
		Like `fft`, but for a real (Float) sequence input.

		Since the input time signal is real, its frequency representation is
		Hermitian-symmetric so we only return the positive frequencies.
	**/
	public static function rfft(input:Array<Float>) : Array<Complex> {
		var n = nextPow2(input.length);
		var real = new Vector<Float>(n);
		var imag = new Vector<Float>(n); // inits to 0
		
		for(i in 0...input.length) real[i] = input[i];
		for(i in input.length...n) real[i] = 0.0;
		
		// imag initialized to 0s by default (or we should force it)
		for(i in 0...n) imag[i] = 0.0; 

		fftVectors(real, imag, false);

		var limit = Std.int(n / 2) + 1;
		var res = new Array<Complex>();
		// res.resize(limit); // resize fills with nulls
		for(i in 0...limit) {
			res.push(new Complex(real[i], imag[i]));
		}
		return res;
	}

	/**
		Computes the Inverse DFT of a periodic input sequence.

		If the input contains N (a power of 2) DTFT samples, each spaced Fs/N Hz
		from each other, the result will consist of N data points as sampled
		from a time signal at intervals of 1/Fs with a scaling factor of 1/Fs.
	**/
	public static function ifft(input:Array<Complex>) : Array<Complex>
		return do_fft(input, true);

	// Handles padding and scaling for forwards and inverse FFTs.
	private static function do_fft(input:Array<Complex>, inverse:Bool) : Array<Complex> {
		// Reusing the wrapper logic implemented in fft()
		var n = nextPow2(input.length);
		var real = new Vector<Float>(n);
		var imag = new Vector<Float>(n);
		
		for(i in 0...input.length) {
			real[i] = input[i].real;
			imag[i] = input[i].imag;
		}
		for(i in input.length...n) {
			real[i] = 0.0;
			imag[i] = 0.0;
		}

		fftVectors(real, imag, inverse);

		var res = new Array<Complex>();
		res.resize(n);
		for(i in 0...n) {
			res[i] = new Complex(real[i], imag[i]);
		}
		return res;
	}


	// Radix-2 Decimation-In-Time variant of Cooley–Tukey's FFT, recursive.
	private static function ditfft2(
		time:Array<Complex>, t:Int,
		freq:Array<Complex>, f:Int,
		n:Int, step:Int, inverse: Bool
	) : Void {
        // Deprecated implementation kept if needed, but ditfft4_soa is preferred
	}
    
    // --- SOA Implementation ---

    private static function ditfft4_soa(
        real:Vector<Float>, rStart:Int, 
        imag:Vector<Float>, iStart:Int, 
        n:Int, step:Int, inverse:Bool
    ):Void {
		if (n == 4) {
            // Base case: Compute the 4-point DFT directly
            // Unrolled for performance
            
            // Indices
            var i0 = rStart;
            var i1 = rStart + step;
            var i2 = rStart + 2 * step;
            var i3 = rStart + 3 * step;
            
            // Inputs
            var r0 = real[i0]; var im0 = imag[i0];
            var r1 = real[i1]; var im1 = imag[i1];
            var r2 = real[i2]; var im2 = imag[i2];
            var r3 = real[i3]; var im3 = imag[i3];

            if (!inverse) {
                // Forward transform hardcoded
                // k=0: sum(x[n])
                var t0_r = r0 + r2; var t0_i = im0 + im2;
                var t1_r = r1 + r3; var t1_i = im1 + im3;
                real[rStart] = t0_r + t1_r;
                imag[iStart] = t0_i + t1_i;

                // k=1: sum(x[n]*W^n) -> W^1 = -i, W^2 = -1, W^3 = i
                // x0 - i*x1 - x2 + i*x3 = (r0-r2) + i(im0-im2) + (-i)(r1-r3) + (im1-im3)
                // = (r0-r2) + (im1-im3) + i(im0-im2 - (r1-r3))
                var t2_r = r0 - r2; var t2_i = im0 - im2;
                var t3_r = r1 - r3; var t3_i = im1 - im3;
                real[rStart + 1] = t2_r + t3_i;
                imag[iStart + 1] = t2_i - t3_r;

                // k=2: sum(x[n]*W^2n) -> W^2 = -1
                // x0 - x1 + x2 - x3 = (r0+r2) - (r1+r3)
                real[rStart + 2] = t0_r - t1_r;
                imag[iStart + 2] = t0_i - t1_i;

                // k=3: sum(x[n]*W^3n) -> W^3 = i, W^6 = -1, W^9 = -i
                // x0 + i*x1 - x2 - i*x3
                real[rStart + 3] = t2_r - t3_i;
                imag[iStart + 3] = t2_i + t3_r;
            } else {
                // Inverse transform (conjugate twiddles, or just swap sign of imaginary parts in twiddle math)
                // Actually easier to just use standard formula with sign flip
                // For N=4, twiddles are 1, i, -1, -i
                
                var t0_r = r0 + r2; var t0_i = im0 + im2;
                var t1_r = r1 + r3; var t1_i = im1 + im3;
                real[rStart] = t0_r + t1_r;
                imag[iStart] = t0_i + t1_i;

                var t2_r = r0 - r2; var t2_i = im0 - im2;
                var t3_r = r1 - r3; var t3_i = im1 - im3;
                
                // k=1: W=i => x1*i, x2*-1, x3*-i
                // x0 + i*x1 - x2 - i*x3
                real[rStart + 1] = t2_r - t3_i;
                imag[iStart + 1] = t2_i + t3_r;

                real[rStart + 2] = t0_r - t1_r;
                imag[iStart + 2] = t0_i - t1_i;

                // k=3: W=-i => x1*-i, x2*-1, x3*i
                // x0 - i*x1 - x2 + i*x3
                real[rStart + 3] = t2_r + t3_i;
                imag[iStart + 3] = t2_i - t3_r;
            }

		} else if (n == 2) {
             var r0 = real[rStart]; var im0 = imag[iStart];
             var r1 = real[rStart+step]; var im1 = imag[iStart+step];
             real[rStart] = r0 + r1; imag[iStart] = im0 + im1;
             real[rStart+1] = r0 - r1; imag[iStart+1] = im0 - im1;
        } else {
			final quarterLen = n >> 2;
			ditfft4_soa(real, rStart, imag, iStart, quarterLen, step * 4, inverse);
			ditfft4_soa(real, rStart + step, imag, iStart + step, quarterLen, step * 4, inverse);
			ditfft4_soa(real, rStart + 2 * step, imag, iStart + 2 * step, quarterLen, step * 4, inverse);
			ditfft4_soa(real, rStart + 3 * step, imag, iStart + 3 * step, quarterLen, step * 4, inverse);
	
            var tr = inverse ? twiddleRealInv : twiddleReal;
            var ti = inverse ? twiddleImagInv : twiddleImag;
            
            // Fix: Calculate stride for twiddle factors
            // The global twiddle array is computed for the top-level N (tr.length).
            // At this recursive level 'n', we need W_n^k.
            // Since W_n^k = W_N^(k * (N/n)), we must stride by N/n.
            var twiddleStride = Std.int(tr.length / n);

			for (k in 0...quarterLen) {
                // Twiddle indices
                // k factors for the 3 branches (branch 0 has twiddle=1)
                // W_N^k, W_N^2k, W_N^3k
                var idx1 = k * twiddleStride;
                var idx2 = idx1 * 2;
                var idx3 = idx1 * 3;

                var w1r = tr[idx1]; var w1i = ti[idx1];
                var w2r = tr[idx2]; var w2i = ti[idx2];
                var w3r = tr[idx3]; var w3i = ti[idx3];
                
                // Indices in arrays
                var p0 = rStart + k; // and iStart + k
                var p1 = rStart + k + quarterLen;
                var p2 = rStart + k + 2 * quarterLen;
                var p3 = rStart + k + 3 * quarterLen;
                
                // Load inputs
                var f0r = real[p0]; var f0i = imag[p0];
                var f1r_raw = real[p1]; var f1i_raw = imag[p1];
                var f2r_raw = real[p2]; var f2i_raw = imag[p2];
                var f3r_raw = real[p3]; var f3i_raw = imag[p3];
                
                // Multiply by twiddles (Complex Mul)
                // f1 = f1_raw * w1
                var f1r = f1r_raw * w1r - f1i_raw * w1i;
                var f1i = f1r_raw * w1i + f1i_raw * w1r;
                
                // f2 = f2_raw * w2
                var f2r = f2r_raw * w2r - f2i_raw * w2i;
                var f2i = f2r_raw * w2i + f2i_raw * w2r;

                // f3 = f3_raw * w3
                var f3r = f3r_raw * w3r - f3i_raw * w3i;
                var f3i = f3r_raw * w3i + f3i_raw * w3r;

                // Butterfly operations
                // val0 = f0 + f1 + f2 + f3
                // val1 = f0 + f1 - f2 - f3 (original code: f0 + f1*i ... wait check logic)
                
                // Original Logic:
                // f0 + f1 + f2 + f3
                // f0 + f1*j - f2 - f3*j -> No, check lines 116-119
                // freq[f+k] = f0 + f1 + f2 + f3
                // freq[f+k+q] = f0 + f1*j - f2 - f3*j  <-- Wait, original code says:
                // f1 is already f1*twiddle1. 
                // Ah, the original code lines 116-119:
                // f0 + f1 + f2 + f3
                // f0 + f1 - f2 - f3  <-- This looks like a Radix-2 step?
                // Let's re-verify the original ditfft4 implementation.
                
                /*
                Original:
                freq[f + k] = f0 + f1 +  f2 + f3;
				freq[f + k + quarterLen] = f0 + f1 - f2 - f3;
				freq[f + k + 2 * quarterLen] = f0 -  f1 - f2 + f3;
				freq[f + k + 3 * quarterLen] = f0 -  f1 + f2 - f3;
                
                This seems to be a specific Radix-4 butterfly optimization.
                */
                
                real[p0] = f0r + f1r + f2r + f3r;
                imag[p0] = f0i + f1i + f2i + f3i;
                
                real[p1] = f0r + f1r - f2r - f3r;
                imag[p1] = f0i + f1i - f2i - f3i;
                
                real[p2] = f0r - f1r - f2r + f3r;
                imag[p2] = f0i - f1i - f2i + f3i;
                
                real[p3] = f0r - f1r + f2r - f3r;
                imag[p3] = f0i - f1i + f2i - f3i;
			}
		}
	}

	private static function ditfft4(time:Array<Complex>, t:Int, freq:Array<Complex>, f:Int, n:Int, step:Int, inverse:Bool):Void {
		// Deprecated
	}

	// Naive O(n^2) DFT, used for testing purposes.
	private static function dft(ts:Array<Complex>, ?inverse:Bool) : Array<Complex> {
        // Unchanged for testing
		if (inverse == null) inverse = false;
		final n = ts.length;
		var fs = new Array<Complex>();
		fs.resize(n);
		for (f in 0...n) {
			var sum = Complex.zero;
			for (t in 0...n) {
				sum += ts[t] * Complex.exp((inverse ? 1 : -1) * 2 * Math.PI * f * t / n);
			}
			fs[f] = inverse ? sum.scale(1 / n) : sum;
		}
		return fs;
	}

	private static var twiddleReal:Vector<Float>;
    private static var twiddleImag:Vector<Float>;
    private static var twiddleRealInv:Vector<Float>;
    private static var twiddleImagInv:Vector<Float>;
    
    // Kept for compatibility but unused in optimized path
	private static var twiddleFactorsInversed:Array<Complex>;
	private static var twiddleFactors:Array<Complex>;

	private static function precomputeTwiddleFactors(maxN:Int, inverse:Bool):Void
	{
		var n:Int = maxN;
        // Radix-4 requires factors up to 3*n/4
        // Original loop: 0...n/2.
        // But for Radix-4 we access k*1, k*2, k*3 where k goes up to n/4.
        // So max index is 3*(n/4) = 0.75n.
        // Let's safe allocate n size.
        
        var tr = new Vector<Float>(n);
        var ti = new Vector<Float>(n);
        
		for (k in 0...n) {
			var constant = -2 * Math.PI / n;
		    var angle = constant * k;
            if (inverse) angle = -angle;
            
            tr[k] = Math.cos(angle);
            ti[k] = Math.sin(angle);
		}

		if (inverse) {
			twiddleRealInv = tr;
            twiddleImagInv = ti;
            // Legacy support
            twiddleFactorsInversed = []; // Placeholder to avoid null checks failing if accessed
        } else {
			twiddleReal = tr;
            twiddleImag = ti;
            twiddleFactors = [];
        }
	}

	private static function computeTwiddle(index, fft_len, inverse:Bool = false)
	{
		var constant = -2 * Math.PI / fft_len;
		var angle = constant * index;

		var result:Complex = new Complex(Math.cos(angle), Math.sin(angle));

		if (inverse)
			return result.conj();
		else
			return result;
	}

	private static function useTwiddleFactor(n:Int, k:Int, inverse:Bool = false):Complex {
		// Compute the index adjustment based on the FFT size n
		// var indexAdjustment:Int = Std.int(twiddleFactors.length / (n / 4));
		var twiddlesToUse = inverse ? twiddleFactorsInversed : twiddleFactors;
		return twiddlesToUse[k];
	}

	/**
		Finds the power of 2 that is equal to or greater than the given natural.
	**/
	public static function nextPow2(x:Int) : Int {
		if (x < 2) return 1;
		else if ((x & (x-1)) == 0) return x;
		var pow = 2;
		x--;
		while ((x >>= 1) != 0) pow <<= 1;
		return pow;
	}

	// testing, but also acts like an example
	static function main() {
		// sampling and buffer parameters
		final Fs = 44100.0;
		final N = 512;
		final halfN = Std.int(N / 2);

		// build a time signal as a sum of sinusoids
		final freqs = [5919.911];
		final ts = [for (n in 0...N) freqs.map(f -> Math.sin(2 * Math.PI * f * n / Fs)).sum()];

		// get positive spectrum and use its symmetry to reconstruct negative domain
		final fs_pos = rfft(ts);
		final fs_fft = new OffsetArray(
			[for (k in -(halfN - 1) ... 0) fs_pos[-k].conj()].concat(fs_pos),
			-(halfN - 1)
		);

		// double-check with naive DFT
		final fs_dft = new OffsetArray(
			dft(ts.map(Complex.fromReal)).circShift(halfN - 1),
			-(halfN - 1)
		);
		final fs_err = [for (k in -(halfN - 1) ... halfN) fs_fft[k] - fs_dft[k]];
		final max_fs_err = fs_err.map(z -> z.magnitude).max();
		if (max_fs_err > 1e-6) haxe.Log.trace('FT Error: ${max_fs_err}', null);

		// find spectral peaks to detect signal frequencies
		final freqis = fs_fft.array.map(z -> z.magnitude)
		                           .findPeaks()
		                           .map(k -> (k - (halfN - 1)) * Fs / N)
		                           .filter(f -> f >= 0);
		if (freqis.length != freqs.length) {
			trace('Found frequencies: ${freqis}');
		} else {
			final freqs_err = [for (i in 0...freqs.length) freqis[i] - freqs[i]];
			final max_freqs_err = freqs_err.map(Math.abs).max();
			if (max_freqs_err > Fs / N) trace('Frequency Errors: ${freqs_err}');
		}

		// recover time signal from the frequency domain
		final ts_ifft = ifft(fs_fft.array.circShift(-(halfN - 1)).map(z -> z.scale(1 / Fs)));
		final ts_err = [for (n in 0...N) ts_ifft[n].scale(Fs).real - ts[n]];
		final max_ts_err = ts_err.map(Math.abs).max();
		if (max_ts_err > 1e-6) haxe.Log.trace('IFT Error: ${max_ts_err}', null);
	}
}
