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

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "saf.h"
#include "saf_externals.h" /* to also include saf dependencies (cblas etc.) */

#include "hosirrlib.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */
//#ifndef HOSIRR_MAX_SH_ORDER
//#define HOSIRR_MAX_SH_ORDER ( HOSIRR_MAX_SH_ORDER )
//#define MAX_NUM_SH_SIGNALS ( (HOSIRR_MAX_SH_ORDER+1)*(HOSIRR_MAX_SH_ORDER+1) ) /* Maximum number of spherical harmonic components */
//#endif
#define MAX_NUM_LOUDSPEAKERS ( HOSIRR_MAX_NUM_OUTPUTS ) /* Maximum permitted channels for the VST standard */
#define MIN_NUM_LOUDSPEAKERS ( 4 )    /* To help avoid traingulation errors when using, e.g. AllRAD */
#define MAX_NUM_LOUDSPEAKERS_IN_PRESET ( MAX_NUM_LOUDSPEAKERS )
#define DEFAULT_WINDOW_LENGTH ( 128 )
#define MAX_WINDOW_LENGTH ( 256 )
#define MIN_WINDOW_LENGTH ( 32 )
#define MAX_DIFF_FREQ_HZ ( 3000 )
#define ALPHA_DIFF_COEFF ( 0.7165f )

/* ========================================================================== */
/*                                 SAF SUBSTITUTES                            */
/* Conflicting versions of SAF will have discrepancies in these funcs, so     */
/* define local versions                                                      */
/* ========================================================================== */

#define HOSIRR_MAX(a,b)      \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b;       \
})

#define HOSIRR_MIN(a,b)      \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b;       \
})

#define HOSIRR_CLAMP(val, min_val, max_val) ({  \
    __typeof__(val) _val = (val);               \
    __typeof__(min_val) _min_val = (min_val);   \
    __typeof__(max_val) _max_val = (max_val);   \
    _val < _min_val ? _min_val : (_val > _max_val ? _max_val : _val); \
})

/** Macro to print a warning message along with the filename and line number */
# define hosirr_print_warning(message) {fprintf(stdout, \
                                    "SAF WARNING: %s [%s LINE %u] \n", message,\
                                    __FILE__, __LINE__);}

/** Macro to print a error message along with the filename and line number */
# define hosirr_print_error(message) {fprintf(stderr, \
                                  "SAF ERROR: %s [%s LINE %u] \n", message, \
                                  __FILE__, __LINE__); \
                                  exit(EXIT_FAILURE);}

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
    float** encBeamCoeffs;      // nSH x nDir
    float** decBeamCoeffs;      // nDir x nSH
    float** dirGainBufDB;       // nBand x nDir
    float*  t60Buf_omni;        // nBand x 1
    float** t60Buf_dir;         // nBand x nDir
    float*  rdrBuf;             // nBand x 1
    int*    directOnsetIdx_bnd; // nBand x 1
    
    // depend on output design (nDir) AND input RIR (nSamp)
    float**  rirBuf_sh;         // nSH x nSamp (input RIR)
    float*** rirBuf_bnd_sh;     // nBand x nSH x nSamp
    float*** rirBuf_bnd_dir;    // nBand x nDir x nSamp
    float**  edcBufOmn_bnd;     // nBand x nSamp
    float*** edcBuf_bnd_dir;    // nBand x nDir x nSamp
    float**  fdnBuf_dir;        // nDir x nSamp
    float*** fdnBuf_bnd_dir;    // nBand nDir x nSamp
    float*** edcBufFDN_bnd_dir; // nBand x nDir x nSamp
    float**  fdnBuf_sh;         // nSH x nSamp (output RIR)
    
    float**  H_bandFilt;        // nBand x bandFiltOrder + 1
    float*   bandCenterFreqs;   // size should match nBand.
    float*   bandXOverFreqs;    // 1 x nBand-1
    float    srcPosition[3];    // [X, Y, Z], used for src-rec distance calc
    float    recPosition[3];    // [X, Y, Z]
    
    int nSH, nSamp, shOrder;
    float fs;
    int nDir, nBand, bandFiltOrder;
    int directOnsetIdx_brdbnd;  // direct arrival onset index within the input RIR
    int diffuseOnsetIdx;        // diffuse onset sample index within the input RIR
    float diffuseOnsetSec;      // diffuse onset from t0, in seconds
    float t0;                   // time-0: when the sound leaves the source (directOnset - sourceDistance/343)
    int t0Idx;                  // sample index of t0 (can be negative)
    float duration;             // seconds
    float diffuseMin;           // minimum diffuseness to detect onset, otherwise it's direct onset +10ms
    float diffuseOnsetFallbackDelay; // delay (sec) after direct onset, to substitute calculation of diffuseOnset in case of diffuseness threshold isn't reached.
    
    ANALYSIS_STAGE analysisStage;
    BEAM_TYPES beamType;
    
    /* original hosirrlib */
    
    float** shir;               // input SRIR [nSH x length]
    float* lsir;                // output LSIR [nLoudpkrs x length]
        
    /* Misc. */
    int ambiRIRorder;
    int ambiRIRlength_samples;
    float ambiRIRlength_seconds;
    int ambiRIRsampleRate; 
    float progress0_1;
    char* progressText;
    
    /* user parameters */
    int analysisOrder;
    int nLoudpkrs;              // number of loudspeakers/virtual loudspeakers
    int windowLength;
    float wetDryBalance;
    int broadBandFirstPeakFLAG;
    float loudpkrs_dirs_deg[MAX_NUM_LOUDSPEAKERS][2];
    CH_ORDERING chOrdering;     // only ACN is supported
    NORMALIZATION_TYPES norm;   // N3D or SN3D

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
//void loadLoudspeakerArrayPreset(LS_ARRAY_PRESETS preset,
//                                float dirs_deg[MAX_NUM_LOUDSPEAKERS_IN_PRESET][2],
//                                int* nCH);

void loadSphDesignPreset(SPHDESIGN_PRESETS preset,
                         float dirs_deg[MAX_NUM_LOUDSPEAKERS_IN_PRESET][2], // TODO: SPHDESIGNs coincide with LS arrays for now
                         int* nCH);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __HOSIRR_INTERNAL_H_INCLUDED__ */
