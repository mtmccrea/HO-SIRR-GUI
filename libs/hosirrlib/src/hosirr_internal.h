/*
 *
 * This library depends on
 * Spatial_Audio_Framework
 * (https://github.com/leomccormack/Spatial_Audio_Framework)
 * and follows the structure of HOSIRR, by the same author.
 *
 * Michael McCrea, Tampere Univeristy, 26.10.23
 *
 */

#ifndef __HOSIRR_INTERNAL_H_INCLUDED__
#define __HOSIRR_INTERNAL_H_INCLUDED__

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "saf.h"
#include "saf_externals.h" // to also include saf dependencies (cblas etc.)

#include <complex.h> // For C99 complex types

#include "hosirrlib.h"

#ifdef __cplusplus
#include <cmath> // For std::exp
#include <complex> // For C++ complex types
extern "C"
{
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */
//#ifndef HOSIRR_MAX_SH_ORDER
//#define HOSIRR_MAX_SH_ORDER ( HOSIRR_MAX_SH_ORDER )
//#define MAX_NUM_SH_SIGNALS ( (HOSIRR_MAX_SH_ORDER+1)*(HOSIRR_MAX_SH_ORDER+1) ) /* Maximum number
//of spherical harmonic components */ #endif
#define MAX_NUM_LOUDSPEAKERS                                                                       \
    (HOSIRR_MAX_NUM_OUTPUTS) /* Maximum permitted channels for the VST standard */
#define MIN_NUM_LOUDSPEAKERS (4) /* To help avoid traingulation errors when using, e.g. AllRAD */
#define MAX_NUM_LOUDSPEAKERS_IN_PRESET (MAX_NUM_LOUDSPEAKERS)
#define DEFAULT_WINDOW_LENGTH (128)
#define MAX_WINDOW_LENGTH (256)
#define MIN_WINDOW_LENGTH (32)
#define MAX_DIFF_FREQ_HZ (3000)
#define ALPHA_DIFF_COEFF (0.7165f)

    /* ========================================================================== */
    /*                                 SAF SUBSTITUTES                            */
    /* Conflicting versions of SAF will have discrepancies in these funcs, so     */
    /* define local versions                                                      */
    /* ========================================================================== */

#ifndef HOSIRR_MIN
#define HOSIRR_MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef HOSIRR_MAX
#define HOSIRR_MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef HOSIRR_CLAMP
#define HOSIRR_CLAMP(a, min, max) (HOSIRR_MAX(min, HOSIRR_MIN(max, a)))
#endif

/** Macro to print a warning message along with the filename and line number */
#define hosirr_print_warning(message)                                                              \
    {                                                                                              \
        fprintf(stdout, "\nRIRLIB WARNING: %s [%s LINE %u] \n", message, __FILE__, __LINE__);      \
    }

/** Macro to print a error message along with the filename and line number */
#define hosirr_print_error(message)                                                                \
    {                                                                                              \
        fprintf(stderr, "\nRIRLIB ERROR: %s [%s LINE %u] \n", message, __FILE__, __LINE__);        \
        exit(EXIT_FAILURE);                                                                        \
    }

    /* ========================================================================== */
    /*                                 Structures                                 */
    /* ========================================================================== */

    /**
     * Main structure for hosirrlib
     */
    typedef struct _hosirrlib
    {
        AMBI_RIR_STATUS ambiRIR_status;
        LS_RIR_STATUS lsRIR_status;

        /* new hodecaylib */

        // depend only on nDir
        float** encBeamCoeffs; // nSH x nDir
        float** decBeamCoeffs; // nDir x nSH
        float** dirGainBufDB; // nBand x nDir
        float* t60Buf_omni; // nBand x 1
        float** t60Buf_dir; // nBand x nDir
        float* rdrBuf; // nBand x 1
        int* directOnsetIdx_bnd; // nBand x 1

        // depend on output design (nDir) AND input RIR (nSamp)
        float** rirBuf_sh; // nSH x nSamp (input RIR)
        float*** rirBuf_bnd_sh; // nBand x nSH x nSamp
        float*** rirBuf_bnd_dir; // nBand x nDir x nSamp
        float** edcBufOmn_bnd; // nBand x nSamp
        float*** edcBuf_bnd_dir; // nBand x nDir x nSamp
        float** fdnBuf_dir; // nDir x nSamp
        float*** fdnBuf_bnd_dir; // nBand nDir x nSamp
        float*** edcBufFDN_bnd_dir; // nBand x nDir x nSamp
        float** fdnBuf_sh; // nSH x nSamp (output RIR)

        float** H_bandFilt; // nBand x bandFiltOrder + 1
        float* bandCenterFreqs; // size should match nBand.
        float* bandXOverFreqs; // 1 x nBand-1
        float* srcDirectivity; // 1 x nBand
        float srcPosition[3]; // [X, Y, Z], used for src-rec distance calc
        float recPosition[3]; // [X, Y, Z]

        int nSH, nSamp, shOrder;
        NORMALIZATION_TYPES inputNorm;
        float fs;
        int nDir, nBand, bandFiltOrder;
        int directOnsetIdx_brdbnd; // direct arrival onset index within the input RIR
        int diffuseOnsetIdx; // diffuse onset sample index within the input RIR
        float diffuseOnsetSec; // diffuse onset from t0, in seconds
        float t0; // time-0: when the sound leaves the source (directOnset - sourceDistance/343)
        int t0Idx; // sample index of t0 (can be negative)
        float duration; // seconds
        float diffuseMin; // minimum diffuseness to detect onset, otherwise it's direct onset +10ms
        float diffuseOnsetFallbackDelay; // delay (sec) after direct onset, to substitute
                                         // calculation of diffuseOnset in case of diffuseness
                                         // threshold isn't reached.
        int srcDirectivityFlag; // boolean: 0: omni, 1: "regular" loudspeaker directivity

        ANALYSIS_STAGE analysisStage;
        BEAM_TYPES beamType;

        /* original hosirrlib */

        float** shir; // input SRIR [nSH x length]
        float* lsir; // output LSIR [nLoudpkrs x length]

        /* Misc. */
        int ambiRIRorder;
        int ambiRIRlength_samples;
        float ambiRIRlength_seconds;
        int ambiRIRsampleRate;
        float progress0_1;
        char* progressText;

        /* user parameters */
        int analysisOrder;
        int nLoudpkrs; // number of loudspeakers/virtual loudspeakers
        int windowLength;
        float wetDryBalance;
        int broadBandFirstPeakFLAG;
        float loudpkrs_dirs_deg[MAX_NUM_LOUDSPEAKERS][2];
        CH_ORDERING chOrdering; // only ACN is supported
        NORMALIZATION_TYPES norm; // N3D or SN3D

    } hosirrlib_data;

    /* ========================================================================== */
    /*                             Internal Functions                             */
    /* ========================================================================== */

    /**
     * Returns the loudspeaker directions for a specified loudspeaker array preset
     * (see LS_ARRAY_PRESETS enum)
     *
     * @param[in]  preset   see LS_ARRAY_PRESETS enum
     * @param[out] dirs_deg loudspeaker directions, [azimuth elevation] convention,
     *                      in DEGREES
     * @param[out] nCH      (&) number of loudspeaker directions in the array
     * @param[out] nDims    (&) number of dimensions (2 or 3)
     */
    // void loadLoudspeakerArrayPreset(LS_ARRAY_PRESETS preset,
    //                                float dirs_deg[MAX_NUM_LOUDSPEAKERS_IN_PRESET][2],
    //                                int* nCH);

    void loadSphDesignPreset(SPHDESIGN_PRESETS preset,
        float dirs_deg[MAX_NUM_LOUDSPEAKERS_IN_PRESET]
                      [2], // TODO: SPHDESIGNs coincide with LS arrays for now
        int* nCH);

    /* ========================================================================== */
    /*                            FIR Filter Functions                            */
    /* The three functions                                                        */
    /* applyWindowingFunction, hosirrlib_FIRCoeffs, hosirrlib_FIRFilterbank       */
    /* mirror those of the same name in SAF, though a                             */
    /* patch has been made correcting band crossovers. The patch is made here     */
    /* to avoid having to bump the version of the current SAF dependency.         */
    /* If the SAF dependency is updated in the future, these can be removed,      */
    /* provided this PR is accepted:                                              */
    /* https://github.com/leomccormack/Spatial_Audio_Framework/pull/61            */
    /* ========================================================================== */

    static void hosirrlib_applyWindowingFunction(
        WINDOWING_FUNCTION_TYPES type, int winlength, float* x)
    {
        int i, N;

        /* if winlength is odd -> symmetric window (mid index has value=1) */
        if (!(winlength % 2 == 0))
            N = winlength - 1;
        /* otherwise, if winlength is even (index: winlength/2+1 = 1.0, but first
         * value != last value) */
        else
            N = winlength;

        switch (type) {
        case WINDOWING_FUNCTION_RECTANGULAR:
            break;

        case WINDOWING_FUNCTION_HAMMING:
            for (i = 0; i < winlength; i++)
                x[i] *= 0.54f
                    - 0.46f
                        * (cosf(2.0f * SAF_PI * (float)i
                            / (float)N)); /* more wide-spread coefficient values */
            /* optimal equiripple coefficient values: */
            /*x[i] *= 0.53836f - 0.46164f * (cosf(2.0f*SAF_PI*(float)i/(float)N));*/
            break;

        case WINDOWING_FUNCTION_HANN:
            for (i = 0; i < winlength; i++)
                x[i] *= 0.5f - 0.5f * (cosf(2.0f * SAF_PI * (float)i / (float)N));
            break;

        case WINDOWING_FUNCTION_BARTLETT:
            for (i = 0; i < winlength; i++)
                x[i] *= 1.0f - 2.0f * fabsf((float)i - ((float)N / 2.0f)) / (float)N;
            break;

        case WINDOWING_FUNCTION_BLACKMAN:
            for (i = 0; i < winlength; i++) {
                x[i] *= 0.42659f - 0.49656f * cosf(2.0f * SAF_PI * (float)i / (float)N)
                    + 0.076849f * cosf(4.0f * SAF_PI * (float)i / (float)N);
            }
            break;

        case WINDOWING_FUNCTION_NUTTALL:
            for (i = 0; i < winlength; i++) {
                x[i] *= 0.355768f - 0.487396f * cosf(2.0f * SAF_PI * (float)i / (float)N)
                    + 0.144232f * cosf(4.0f * SAF_PI * (float)i / (float)N)
                    - 0.012604f * cosf(6.0f * SAF_PI * (float)i / (float)N);
            }
            break;

        case WINDOWING_FUNCTION_BLACKMAN_NUTTALL:
            for (i = 0; i < winlength; i++) {
                x[i] *= 0.3635819f - 0.4891775f * cosf(2.0f * SAF_PI * (float)i / (float)N)
                    + 0.1365995f * cosf(4.0f * SAF_PI * (float)i / (float)N)
                    + 0.0106411f * cosf(4.0f * SAF_PI * (float)i / (float)N);
            }
            break;

        case WINDOWING_FUNCTION_BLACKMAN_HARRIS:
            for (i = 0; i < winlength; i++) {
                x[i] *= 0.35875f - 0.48829f * cosf(2.0f * SAF_PI * (float)i / (float)N)
                    + 0.14128f * cosf(4.0f * SAF_PI * (float)i / (float)N)
                    + 0.01168f * cosf(4.0f * SAF_PI * (float)i / (float)N);
            }
            break;
        }
    }

    static void hosirrlib_FIRCoeffs(FIR_FILTER_TYPES filterType, int order, float fc1,
        float fc2, /* only needed for band-pass/stop */
        float fs, WINDOWING_FUNCTION_TYPES windowType, int scalingFLAG, float* h_filt)
    {
        int i, h_len;
        float ft1, ft2, h_sum, f0;
        float_complex h_z_sum;

        h_len = order + 1;
        ft1 = fc1
            / fs; // corrected, was ft1 = fc1/(fs*2.0f); see §8.1.1, J. G. Proakis and D. G.
                  // Manolakis, Digital Signal Processing: Principles, Algorithms, and Applications.

        /* compute filter weights */
        if (order % 2 == 0) {
            /* if order is multiple of 2 */
            switch (filterType) {
            case FIR_FILTER_LPF:
                for (i = 0; i < h_len; i++)
                    h_filt[i] = i == order / 2 ? 2.0f * ft1
                                               : sinf(2.0f * SAF_PI * ft1 * (float)(i - order / 2))
                            / (SAF_PI * (float)(i - order / 2));
                break;

            case FIR_FILTER_HPF:
                for (i = 0; i < h_len; i++)
                    h_filt[i] = i == order / 2 ? 1.0f - 2.0f * ft1
                                               : -sinf(2.0f * ft1 * SAF_PI * (float)(i - order / 2))
                            / (SAF_PI * (float)(i - order / 2));
                break;

            case FIR_FILTER_BPF:
                ft2 = fc2 / fs; // corrected, was fc2/(fs*2.0f);
                for (i = 0; i < h_len; i++) {
                    h_filt[i] = i == order / 2 ? 2.0f * (ft2 - ft1)
                                               : sinf(2.0f * SAF_PI * ft2 * (float)(i - order / 2))
                                / (SAF_PI * (float)(i - order / 2))
                            - sinf(2.0f * SAF_PI * ft1 * (float)(i - order / 2))
                                / (SAF_PI * (float)(i - order / 2));
                }
                break;

            case FIR_FILTER_BSF:
                ft2 = fc2 / fs; // corrected, was fc2/(fs*2.0f);
                for (i = 0; i < h_len; i++) {
                    h_filt[i] = i == order / 2 ? 1.0f - 2.0f * (ft2 - ft1)
                                               : sinf(2.0f * SAF_PI * ft1 * (float)(i - order / 2))
                                / (SAF_PI * (float)(i - order / 2))
                            - sinf(2.0f * SAF_PI * ft2 * (float)(i - order / 2))
                                / (SAF_PI * (float)(i - order / 2));
                }
                break;
            }
        } else
            hosirr_print_error("Please specify an even value for the filter 'order' argument");

        /* Apply windowing function */
        hosirrlib_applyWindowingFunction(windowType, h_len, h_filt);

        /* Scaling, to ensure pass-band is truely at 1 (0dB).
         * [1] "Programs for Digital Signal Processing", IEEE Press John Wiley &
         *     Sons, 1979, pg. 5.2-1.
         */
        if (scalingFLAG) {
            switch (filterType) {
            case FIR_FILTER_LPF:
            case FIR_FILTER_BSF:
                h_sum = 0.0f;
                for (i = 0; i < h_len; i++)
                    h_sum += h_filt[i];
                for (i = 0; i < h_len; i++)
                    h_filt[i] /= h_sum;
                break;
#ifdef __cplusplus // TODO: not awesome, colliding complex types/functions when compiling in C++
                   // project!
            case FIR_FILTER_HPF:
                f0 = 1.0f;
                h_z_sum = cmplxf(0.0f, 0.0f);
                for (i = 0; i < h_len; i++)
                    h_z_sum = ccaddf(h_z_sum,
                        crmulf(std::exp(cmplxf(0.0f, -2.0f * SAF_PI * (float)i * f0 / 2.0f)),
                            h_filt[i]));
                h_sum = std::abs(h_z_sum);
                for (i = 0; i < h_len; i++)
                    h_filt[i] /= h_sum;
                break;

            case FIR_FILTER_BPF:
                f0 = fc1 / fs + fc2 / fs; // correct, was (fc1/fs+fc2/fs)/2.0f;
                h_z_sum = cmplxf(0.0f, 0.0f);
                for (i = 0; i < h_len; i++)
                    h_z_sum = ccaddf(h_z_sum,
                        crmulf(std::exp(cmplxf(0.0f, -2.0f * SAF_PI * (float)i * f0 / 2.0f)),
                            h_filt[i]));
                h_sum = std::abs(h_z_sum);
                for (i = 0; i < h_len; i++)
                    h_filt[i] /= h_sum;
                break;
#else
        case FIR_FILTER_HPF:
            f0 = 1.0f;
            h_z_sum = cmplxf(0.0f, 0.0f);
            for (i = 0; i < h_len; i++)
                h_z_sum = ccaddf(h_z_sum,
                    crmulf(cexpf(cmplxf(0.0f, -2.0f * SAF_PI * (float)i * f0 / 2.0f)), h_filt[i]));
            h_sum = cabsf(h_z_sum);
            for (i = 0; i < h_len; i++)
                h_filt[i] /= h_sum;
            break;

        case FIR_FILTER_BPF:
            f0 = fc1 / fs + fc2 / fs; // correct, was (fc1/fs+fc2/fs)/2.0f;
            h_z_sum = cmplxf(0.0f, 0.0f);
            for (i = 0; i < h_len; i++)
                h_z_sum = ccaddf(h_z_sum,
                    crmulf(cexpf(cmplxf(0.0f, -2.0f * SAF_PI * (float)i * f0 / 2.0f)), h_filt[i]));
            h_sum = cabsf(h_z_sum);
            for (i = 0; i < h_len; i++)
                h_filt[i] /= h_sum;
            break;
#endif
            }
        }
    }

    static void hosirrlib_FIRFilterbank(
        int order, float* fc, /* cut-off frequencies; nCutoffFreq x 1 */
        int nCutoffFreq, float sampleRate, WINDOWING_FUNCTION_TYPES windowType, int scalingFLAG,
        float* filterbank /* (nCutoffFreq+1) x (order+1) */
    )
    {
        int k, nFilt;

        /* Number of filters returned is always one more than the number of cut-off frequencies */
        nFilt = nCutoffFreq + 1;

        /* first and last bands are low-pass and high pass filters, using the first
         * and last cut-off frequencies in vector 'fc', respectively.  */
        hosirrlib_FIRCoeffs(
            FIR_FILTER_LPF, order, fc[0], 0.0f, sampleRate, windowType, scalingFLAG, filterbank);
        hosirrlib_FIRCoeffs(FIR_FILTER_HPF, order, fc[nCutoffFreq - 1], 0.0f, sampleRate,
            windowType, scalingFLAG, &filterbank[(nFilt - 1) * (order + 1)]);

        /* the inbetween bands are then band-pass filters: */
        if (nCutoffFreq > 1) {
            for (k = 1; k < nFilt - 1; k++)
                hosirrlib_FIRCoeffs(FIR_FILTER_BPF, order, fc[k - 1], fc[k], sampleRate, windowType,
                    scalingFLAG, &filterbank[k * (order + 1)]);
        }
    }

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __HOSIRR_INTERNAL_H_INCLUDED__ */
