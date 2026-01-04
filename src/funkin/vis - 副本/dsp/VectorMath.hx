package funkin.vis.dsp;

import haxe.ds.Vector;
#if cpp
import cpp.Pointer;
import cpp.NativeArray;
#end

class VectorMath {
    
    // --- Basic Arithmetic ---

    /**
     * out[i] = a[i] + b[i]
     */
    public static inline function add(a:Vector<Float>, b:Vector<Float>, out:Vector<Float>):Void {
        var len = a.length;
        #if cpp
        // C++ specific optimization hint: Use pointers if possible or rely on simple loops for auto-vectorization
        var ptrA = cpp.Pointer.ofArray(a.toData()).raw;
        var ptrB = cpp.Pointer.ofArray(b.toData()).raw;
        var ptrOut = cpp.Pointer.ofArray(out.toData()).raw;
        untyped __cpp__("
            double* pA = (double*){0};
            double* pB = (double*){1};
            double* pOut = (double*){2};
            for(int i = 0; i < {3}; i++) {
                pOut[i] = pA[i] + pB[i];
            }
        ", ptrA, ptrB, ptrOut, len);
        #else
        for (i in 0...len) out[i] = a[i] + b[i];
        #end
    }

    /**
     * out[i] = a[i] * b[i]
     */
    public static inline function mul(a:Vector<Float>, b:Vector<Float>, out:Vector<Float>):Void {
        var len = a.length;
        #if cpp
        var ptrA = cpp.Pointer.ofArray(a.toData()).raw;
        var ptrB = cpp.Pointer.ofArray(b.toData()).raw;
        var ptrOut = cpp.Pointer.ofArray(out.toData()).raw;
        untyped __cpp__("
            double* pA = (double*){0};
            double* pB = (double*){1};
            double* pOut = (double*){2};
            for(int i = 0; i < {3}; i++) {
                pOut[i] = pA[i] * pB[i];
            }
        ", ptrA, ptrB, ptrOut, len);
        #else
        for (i in 0...len) out[i] = a[i] * b[i];
        #end
    }

    /**
     * out[i] = a[i] * scalar
     */
    public static inline function mulScalar(a:Vector<Float>, scalar:Float, out:Vector<Float>):Void {
        var len = a.length;
        #if cpp
        var ptrA = cpp.Pointer.ofArray(a.toData()).raw;
        var ptrOut = cpp.Pointer.ofArray(out.toData()).raw;
        untyped __cpp__("
            double* pA = (double*){0};
            double* pOut = (double*){1};
            double s = {2};
            for(int i = 0; i < {3}; i++) {
                pOut[i] = pA[i] * s;
            }
        ", ptrA, ptrOut, scalar, len);
        #else
        for (i in 0...len) out[i] = a[i] * scalar;
        #end
    }

    /**
     * out[i] += a[i] * scalar
     * (Fused Multiply-Add approximation)
     */
    public static inline function addMulScalar(out:Vector<Float>, a:Vector<Float>, scalar:Float):Void {
        var len = out.length;
        #if cpp
        var ptrA = cpp.Pointer.ofArray(a.toData()).raw;
        var ptrOut = cpp.Pointer.ofArray(out.toData()).raw;
        untyped __cpp__("
            double* pA = (double*){0};
            double* pOut = (double*){1};
            double s = {2};
            for(int i = 0; i < {3}; i++) {
                pOut[i] += pA[i] * s;
            }
        ", ptrA, ptrOut, scalar, len);
        #else
        for (i in 0...len) out[i] += a[i] * scalar;
        #end
    }

    // --- Complex Number Operations (SoA) ---

    /**
     * Complex Multiplication: (ar + i*ai) * (br + i*bi)
     * Real part: ar*br - ai*bi
     * Imag part: ar*bi + ai*br
     */
    public static inline function complexMul(
        ar:Vector<Float>, ai:Vector<Float>, 
        br:Vector<Float>, bi:Vector<Float>, 
        outR:Vector<Float>, outI:Vector<Float>
    ):Void {
        var len = ar.length;
        #if cpp
        var pAr = cpp.Pointer.ofArray(ar.toData()).raw;
        var pAi = cpp.Pointer.ofArray(ai.toData()).raw;
        var pBr = cpp.Pointer.ofArray(br.toData()).raw;
        var pBi = cpp.Pointer.ofArray(bi.toData()).raw;
        var pOutR = cpp.Pointer.ofArray(outR.toData()).raw;
        var pOutI = cpp.Pointer.ofArray(outI.toData()).raw;
        untyped __cpp__("
            double* p_ar = (double*){0};
            double* p_ai = (double*){1};
            double* p_br = (double*){2};
            double* p_bi = (double*){3};
            double* p_outR = (double*){4};
            double* p_outI = (double*){5};
            for(int i = 0; i < {6}; i++) {
                double a_re = p_ar[i];
                double a_im = p_ai[i];
                double b_re = p_br[i];
                double b_im = p_bi[i];
                p_outR[i] = a_re * b_re - a_im * b_im;
                p_outI[i] = a_re * b_im + a_im * b_re;
            }
        ", pAr, pAi, pBr, pBi, pOutR, pOutI, len);
        #else
        for (i in 0...len) {
            var a_re = ar[i];
            var a_im = ai[i];
            var b_re = br[i];
            var b_im = bi[i];
            outR[i] = a_re * b_re - a_im * b_im;
            outI[i] = a_re * b_im + a_im * b_re;
        }
        #end
    }

    /**
     * Calculates magnitude squared for complex numbers: re*re + im*im
     */
    public static inline function complexMagSq(
        re:Vector<Float>, im:Vector<Float>, out:Vector<Float>
    ):Void {
        var len = re.length;
        #if cpp
        var pRe = cpp.Pointer.ofArray(re.toData()).raw;
        var pIm = cpp.Pointer.ofArray(im.toData()).raw;
        var pOut = cpp.Pointer.ofArray(out.toData()).raw;
        untyped __cpp__("
            double* p_re = (double*){0};
            double* p_im = (double*){1};
            double* p_out = (double*){2};
            for(int i = 0; i < {3}; i++) {
                p_out[i] = p_re[i] * p_re[i] + p_im[i] * p_im[i];
            }
        ", pRe, pIm, pOut, len);
        #else
        for (i in 0...len) out[i] = re[i] * re[i] + im[i] * im[i];
        #end
    }
}
