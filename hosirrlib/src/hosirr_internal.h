/*
 ==============================================================================
 
 This file is part of HOSIRR
 Copyright (c) 2020 - Leo McCormack
 
 HOSIRR is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 HOSIRR is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with HOSIRR.  If not, see <http://www.gnu.org/licenses/>.
 
 ==============================================================================
 */

/*
 * Filename: hosirr_internal.h
 * ---------------------------
 * hosirrlib is a C-port of the Higher-order Spatial Impulse Response Rendering
 * (HO-SIRR) Matlab toolbox: https://github.com/leomccormack/HO-SIRR
 * HO-SIRR is a rendering method, which can synthesise output loudspeaker array
 * room impulse responses (RIRs) using input spherical harmonic (Ambisonic/
 * B-Format) RIRs of arbitrary order. The method makes assumptions regarding
 * the composition of the sound-field and extracts spatial parameters over time,
 * which allows it to map the input to the output in an adaptive and informed
 * manner.
 * The idea is that you then convolve a monophonic source with this loudspeaker
 * array RIR, and it will be reproduced and exhibit the spatial characteristics
 * of the captured space more faithfully (when compared to linear methods such
 * as Ambisonics).
 *
 * Dependencies:
 *     Spatial_Audio_Framework
 * Author, date created:
 *     Leo McCormack, 04.01.2020
 *
 * [1] McCormack, L., Politis, A., Scheuregger, O., and Pulkki, V. (2019).
 *     "Higher-order processing of spatial impulse responses". In Proceedings of
 *     the 23rd International Congress on Acoustics, 9--13 September 2019 in
 *     Aachen, Germany.
 */

#ifndef __HOSIRR_INTERNAL_H_INCLUDED__
#define __HOSIRR_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "hosirrlib.h"
#include "saf.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#define MAX_SH_ORDER ( HOSIRR_MAX_SH_ORDER )
#define MAX_NUM_SH_SIGNALS ( (MAX_SH_ORDER+1)*(MAX_SH_ORDER+1) ) /* Maximum number of spherical harmonic components */
#define MAX_NUM_LOUDSPEAKERS ( HOSIRR_MAX_NUM_OUTPUTS ) /* Maximum permitted channels for the VST standard */
#define MIN_NUM_LOUDSPEAKERS ( 4 )    /* To help avoid traingulation errors when using, e.g. AllRAD */
#define MAX_NUM_LOUDSPEAKERS_IN_PRESET ( MAX_NUM_LOUDSPEAKERS )
#define DEFAULT_WINDOW_LENGTH ( 128 )
#define MAX_WINDOW_LENGTH ( 256 )
#define MIN_WINDOW_LENGTH ( 32 )
#define MAX_DIFF_FREQ_HZ ( 3000 )
#define ALPHA_DIFF_COEFF ( 0.5f )


/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/*
 * Struct: hosirrlib
 * -----------------
 * Main structure for hosirrlib. Contains variables for ....
 */
typedef struct _hosirrlib
{
    AMBI_RIR_STATUS ambiRIR_status;
    LS_RIR_STATUS lsRIR_status;
    float* shir;
    float* lsir;
    
    /* Misc. */
    int ambiRIRorder;
    int ambiRIRlength_samples;
    float ambiRIRlength_seconds;
    int ambiRIRsampleRate; 
    float progress0_1;
    char* progressText;
    
    /* user parameters */
    int analysisOrder;
    int nLoudpkrs;         /* number of loudspeakers/virtual loudspeakers */
    int windowLength;
    float wetDryBalance;
    int broadBandFirstPeakFLAG;
    float loudpkrs_dirs_deg[MAX_NUM_LOUDSPEAKERS][2];
    CH_ORDER chOrdering;   /* only ACN is supported */
    NORM_TYPES norm;       /* N3D or SN3D */
    
    
} hosirrlib_data;


/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/*
 * Function: loadLoudspeakerArrayPreset
 * ------------------------------------
 * Returns the loudspeaker directions for a specified loudspeaker array preset.
 * The function also returns the number of loudspeakers in the array, and a
 * crude estimate of whether it is 2D or 3D (based on whether the sum of the
 * absolute loudspeaker elevation angle is not zero).
 * Note: 2-D setups which are not on the horizontal plane (i.e.
 * all elevations==0) will be registered as a 3-D setup.
 * Further Note: default uniformly distributed points are used to pad the
 * dirs_deg matrix up to the MAX_NUM_LOUDSPEAKERS_IN_PRESET, if nCH is less than
 * this. This can help avoid scenarios of many sources being panned in the same
 * direction, or triangulations errors.
 *
 * Input Arguments:
 *     preset   - see "LOUDSPEAKER_ARRAY_PRESETS" enum
 * Output Arguments:
 *     dirs_deg - loudspeaker directions, [azimuth elevation] convention, in
 *                DEGREES;
 *     nCH      - & number of loudspeaker directions in the array
 *     nDims    - & number of dimensions (2 or 3)
 */
void loadLoudspeakerArrayPreset(LOUDSPEAKER_ARRAY_PRESETS preset,
                                float dirs_deg[MAX_NUM_LOUDSPEAKERS_IN_PRESET][2],
                                int* nCH);
    
    
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __HOSIRR_INTERNAL_H_INCLUDED__ */
