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

/**
 * @file hosirrlib.c
 * @brief A C-port of the Higher-order Spatial Impulse Response Rendering
 *        (HO-SIRR) Matlab toolbox: https://github.com/leomccormack/HO-SIRR
 *
 * HO-SIRR is a rendering method, which can synthesise output loudspeaker array
 * room impulse responses (RIRs) using input spherical harmonic (Ambisonic/
 * B-Format) RIRs of arbitrary order. The method makes assumptions regarding
 * the composition of the sound-field and extracts spatial parameters over time,
 * which allows it to map the input to the output in an adaptive and informed
 * manner.
 *
 * The idea is that you then convolve a monophonic source with this loudspeaker
 * array RIR, and it will be reproduced and exhibit the spatial characteristics
 * of the captured space more faithfully (when compared to linear methods such
 * as Ambisonics).
 *
 * Dependencies: Spatial_Audio_Framework
 * (https://github.com/leomccormack/Spatial_Audio_Framework)
 *
 * @see [1] McCormack, L., Politis, A., Scheuregger, O., and Pulkki, V. (2019).
 *          "Higher-order processing of spatial impulse responses". In
 *          Proceedings of the 23rd International Congress on Acoustics, 9--13
 *          September 2019 in Aachen, Germany.
 *
 * @author Leo McCormack
 * @date 04.01.2020
 */

#include "hosirr_internal.h"

void hosirrlib_create(
                      void ** const phHS
                      )
{
    hosirrlib_data* pData = (hosirrlib_data*)malloc1d(sizeof(hosirrlib_data));
    *phHS = (void*)pData;

    /* new hodecaylib */
    
    // depend only on nDir
    pData->encBeamCoeffs = NULL;    // nSH x nDir
    pData->decBeamCoeffs = NULL;    // nDir x nSH
    pData->dirGainBuf    = NULL;    // nDir x nBand
    pData->t60Buf_omni   = NULL;    // nDir x 1
    pData->t60Buf_dir    = NULL;    // nDir x nBand
    
    // depend on output design (nDir) AND input RIR (nSamp)
    pData->rirBuf_sh        = NULL;    // nSH x nSamp
    pData->rirBuf_bnd_dir   = NULL;    // nDir x nSamp
    pData->fdnBuf_dir       = NULL;    // nDir x nSamp
    pData->fdnBuf_bnd_dir   = NULL;    // nBand x nDir x nSamp
    pData->edcBufOmn_bnd    = NULL;    // nBand x nSamp
    pData->edcBuf_bnd_dir   = NULL;    // nDir x nBand x nSamp
    pData->edcBufFDN_bnd_dir = NULL;   // nDir x nBand x nSamp
    pData->fdnBuf_sh        = NULL;    // nSH x nSamp
    
    pData->H_bandFilt      = NULL;  // nBand x filtOrder+1
    pData->bandCenterFreqs = NULL;
    pData->bandXOverFreqs  = NULL;
    
    /* Constants */
    
    /* If diffuseness never crosses this threshold, diffuse onset
     * defaults to direct onset +5ms. */
    pData->diffuseMin = 0.3f;
    
    /* Note: Don't initialize filters yet... input RIR is required for getting
     * the fsfilterbank constants set in _initBandFilters(). */
    
    /* zero out the state of buffer resources */
    hosirrlib_setUninitialized(pData);
    
    /* Beam shape for decomposition */
    pData->beamType = HYPERCARDIOID_BEAM;
    /* Initialize spherical design for directional analysis */
//    loadLoudspeakerArrayPreset( // default design currently just set by enum LS_ARRAY_PRESET_15PX
//                               LS_ARRAY_PRESET_15PX,
//                               pData->loudpkrs_dirs_deg,
//                               &(pData->nLoudpkrs));
    
    /* Initialize spherical design for directional analysis */
    // TODO: SPHDESIGNs coincide with LS arrays for now, not a clean separation
    loadSphDesignPreset( // default design currently just set by enum SPHDESIGN_ARRAY_PRESET_15PX
                        SPHDESIGN_PRESET_15PX,
                        pData->loudpkrs_dirs_deg,
                        &(pData->nLoudpkrs));
    pData->nDir = pData->nLoudpkrs;
    
    /* ~~~~~~~~~~~~~~~~~~ */
    /* original hosirrlib */
    
    pData->progress0_1 = 0.0f;
    pData->progressText = malloc1d(HOSIRR_PROGRESSTEXT_CHAR_LENGTH * sizeof(char));
    strcpy(pData->progressText,"HOSIRR");
    
    /* input AmbiRIR */
    pData->shir = NULL;
    pData->ambiRIR_status = AMBI_RIR_STATUS_NOT_LOADED;
    pData->ambiRIRorder = -1;
    pData->ambiRIRlength_seconds = pData->ambiRIRlength_samples = 0.0f;
    pData->ambiRIRsampleRate = 0;
    
    /* output Loudspeaker array RIR */
    pData->lsir = NULL;
    pData->lsRIR_status = LS_RIR_STATUS_NOT_RENDERED;
    
    /* default user parameters */
    pData->analysisOrder = 1;
    //loadLoudspeakerArrayPreset(LS_ARRAY_PRESET_T_DESIGN_24, pData->loudpkrs_dirs_deg, &(pData->nLoudpkrs));
    pData->chOrdering = ACN_ORDER;
    pData->norm = SN3D_NORM;
    pData->broadBandFirstPeakFLAG = 1;
    pData->windowLength = DEFAULT_WINDOW_LENGTH;
    pData->wetDryBalance = 1.0f;
}


void hosirrlib_destroy(
                       void ** const phHS
                       )
{
    hosirrlib_data *pData = (hosirrlib_data*)(*phHS);
    
    if (pData != NULL) {
        
        /* new hodecaylib */
        // depend only on output design (nDir)
        free(pData->encBeamCoeffs);
        free(pData->decBeamCoeffs);
        free(pData->dirGainBuf);
        free(pData->t60Buf_omni);
        free(pData->t60Buf_dir);
        // depend on output design (nDir) AND input RIR (nSamp)
        free(pData->rirBuf_sh);
        free(pData->rirBuf_bnd_sh);
        free(pData->rirBuf_bnd_dir);
        free(pData->edcBufOmn_bnd);
        free(pData->edcBuf_bnd_dir);
        free(pData->fdnBuf_dir);
        free(pData->fdnBuf_bnd_dir);
        free(pData->edcBufFDN_bnd_dir);
        free(pData->fdnBuf_sh);
        // unchanging after init
        free(pData->H_bandFilt);
        free(pData->bandCenterFreqs);
        free(pData->bandXOverFreqs);
        
        /* original hosirrlib */
        free(pData->shir);
        free(pData->lsir);
        free(pData->progressText);
        
        free(pData);
        pData = NULL;
    }
}


int hosirrlib_setRIR(
                     void* const hHS,
                     const float** H,
                     int numChannels,
                     int numSamples,
                     int sampleRate
                     )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("setRIR called\n");
    
    /* Check channel count to see if input is actually in the SHD */
    if (fabsf(sqrtf((float)numChannels) - floorf(sqrtf((float)numChannels))) > 0.0001f) {
        
        /* new vars */
        // TODO: only set state to unitialized if there was no previous RIR loaded
        //setRIRState_uninitialized(pData)
        printf("Input file doesn't appear to be a SHD file\n");
        
        /* old vars */
        pData->ambiRIR_status = AMBI_RIR_STATUS_NOT_LOADED;
        pData->ambiRIRorder = -1;
        pData->ambiRIRlength_seconds = pData->ambiRIRlength_samples = 0.0f;
        pData->ambiRIRsampleRate = 0;
        return (int)(pData->ambiRIR_status);
    }
    
    /* new vars */
    pData->nSH      = numChannels;
    pData->nSamp    = numSamples;
    pData->shOrder  = HOSIRR_MIN(sqrt(numChannels-1), HOSIRR_MAX_SH_ORDER);
    pData->fs       = (float)sampleRate;
    pData->duration = numSamples / (float)sampleRate;
    
    /* (Re)Alloc and copy in input RIR */
    pData->rirBuf_sh = (float**)realloc2d((void**)pData->rirBuf_sh, numChannels, numSamples, sizeof(float));
    for(int i = 0; i < numChannels; i++)
        utility_svvcopy(H[i], numSamples, pData->rirBuf_sh[i]);
    
    /* set FLAGS */
    pData->analysisStage = RIR_LOADED;
    
    /* Initialize filters */
    hosirrlib_initBandFilters(pData);
    
    /* Alloc processing resources */
    hosirrlib_allocProcBufs(pData);
    
    /* old vars */
    pData->ambiRIRorder = HOSIRR_MIN(sqrt(numChannels-1), HOSIRR_MAX_SH_ORDER);
    pData->analysisOrder = pData->ambiRIRorder;
    pData->ambiRIRlength_samples = numSamples;
    pData->ambiRIRsampleRate = sampleRate;
    pData->ambiRIRlength_seconds = (float)pData->ambiRIRlength_samples / (float)pData->ambiRIRsampleRate;
    
    // (re)allocate memory and copy in SH RIR data
    pData->shir = (float**)realloc2d((void**)pData->shir, numChannels, numSamples, sizeof(float));
    for(int i = 0; i < numChannels; i++)
        utility_svvcopy(H[i], numSamples, pData->shir[i]);
        // memcpy(&(pData->shir[i*numSamples]), H[i], numSamples * sizeof(float));

    /* set FLAGS */
    pData->ambiRIR_status = AMBI_RIR_STATUS_LOADED;
    pData->lsRIR_status = LS_RIR_STATUS_NOT_RENDERED;
    
    return (int)(pData->ambiRIR_status); // TODO: consider use of returned value
}


void hosirrlib_setUninitialized(
                                void* const hHS
                                )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("setUninitialized called.\n"); // DEBUG
    // pending initializations
    pData->nSH = -1; // input vars (assume for now in params = out params)
    pData->nSamp = -1;
    pData->shOrder = -1;
    pData->fs = -1.f;
    pData->directOnsetIdx = -1;
    pData->diffuseOnsetIdx = -1;
    pData->diffuseOnsetSec = 0;
    pData->sourceDistance = 3.f; // default
    pData->t0 = -1;
    pData->t0Idx = 0;
    pData->duration = 0.0f;
    pData->analysisStage = RIR_NOT_LOADED;
}


void hosirrlib_initBandFilters(
                               void* const hHS
                               )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("initBandFilters called.\n"); // DEBUG
    
    /* Check previous stages are complete */
    if (pData->analysisStage < FILTERS_INTITIALIZED-1) { // TODO: handle fail case
        printf("allocProcBufs called before previous stages were completed: %d\n", pData->analysisStage);
        return;
    }
    
    /* Initialize octave band filters */
    // bandCenterFreqs are "center freqs" used to determine crossover freqs.
    // The lowest and highest bands are LP and HP filters, "center freq" for
    // those is a bit of a misnomer. Xover freqs are center freqs * sqrt(2)
    // (octave bands), see hosirrlib_initBandFilters().
    // LPF cutoff is 125 * sqrt(2) = 177 Hz.
    // HPF cuton is 8k * sqrt(2) = 11313 Hz.
    // Must be the same size as pData->nBand.
    pData->nBand = 8;
    float bandCenterFreqs[8] = {125.f, 250.f, 500.f, 1000.f, 2000.f, 4000.f, 8000.f, 16000.f};
    pData->bandFiltOrder = 400;
    
    pData->bandCenterFreqs = malloc1d(pData->nBand     * sizeof(float));
    pData->bandXOverFreqs  = malloc1d((pData->nBand-1) * sizeof(float));
    
    // populate the bandCenterFreqs member var for external access
    utility_svvcopy(bandCenterFreqs, pData->nBand, pData->bandCenterFreqs);
    
    // set xover freqs
    for (int ib = 0; ib < pData->nBand-1; ib++) {
        pData->bandXOverFreqs[ib] = bandCenterFreqs[ib] * sqrtf(2.f);
        printf("band %d xover %.4f\n", ib, pData->bandXOverFreqs[ib]);
    }
    
    /* Create the filterbank */
    if (pData->H_bandFilt == NULL) { // if this is the first time this function is called...
        
        // Allocate filter coefficients (nBand x bandFiltOrder + 1)
        pData->H_bandFilt = (float**)realloc2d((void**)pData->H_bandFilt,
                                               pData->nBand,
                                               pData->bandFiltOrder + 1,
                                               sizeof(float));
        // Compute FIR Filterbank coefficients
        FIRFilterbank(pData->bandFiltOrder,
                      pData->bandXOverFreqs,
                      pData->nBand-1,
                      pData->fs,
                      WINDOWING_FUNCTION_HAMMING,
                      1,
                      FLATTEN2D(pData->H_bandFilt)
                      );
        // DEBUG
        printf("FS FOR FILTERBANK\n\t%.1f\n", pData->fs);
        //for (int ib = 0; ib<pData->nBand; ib++) {
        //    printf("band %d coeffs\n", ib);
        //    for (int icoeff = 0; icoeff<pData->bandFiltOrder+1; icoeff++)
        //        printf("\t%.9f\n", pData->H_bandFilt[ib][icoeff]);
        //}
    }
    
//    hosirrlib_inspectFilts(pData);
    pData->analysisStage = FILTERS_INTITIALIZED;
}


// (re)allocate the buffers used for storing intermediate processing data
void hosirrlib_allocProcBufs(
                             void * const hHS
                             )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("allocProcBufs called\n");
    /* Check previous stages are complete */
    if (pData->analysisStage < ANALYSIS_BUFS_LOADED-1) { // TODO: handle fail case
        printf("allocProcBufs called before previous stages were completed: %d\n", pData->analysisStage);
        return;
    }
    
    const int nSH   = pData->nSH;
    const int nDir  = pData->nDir;
    const int nBand = pData->nBand;
    const int nSamp = pData->nSamp;

    // These members depend on output design (nDir) AND input RIR (nSH, nSamp)
    pData->rirBuf_bnd_sh = (float***)realloc3d((void***)pData->rirBuf_bnd_sh,
                                             nBand, nSH, nSamp, sizeof(float));
    pData->rirBuf_bnd_dir = (float***)realloc3d((void***)pData->rirBuf_bnd_dir,
                                             nBand, nDir, nSamp, sizeof(float));
    pData->edcBufOmn_bnd   = (float**)realloc2d((void**)pData->edcBufOmn_bnd,
                                              nBand, nSamp, sizeof(float));
    pData->edcBuf_bnd_dir   = (float***)realloc3d((void***)pData->edcBuf_bnd_dir,
                                              nBand, nDir, nSamp, sizeof(float));
    pData->fdnBuf_dir       = (float**)realloc2d((void**)pData->fdnBuf_dir,
                                             nDir, nSamp, sizeof(float));
    pData->fdnBuf_bnd_dir = (float***)realloc3d((void***)pData->fdnBuf_bnd_dir,
                                             nBand, nDir, nSamp, sizeof(float));
    pData->edcBufFDN_bnd_dir   = (float***)realloc3d((void***)pData->edcBufFDN_bnd_dir,
                                              nBand, nDir, nSamp, sizeof(float));
    pData->fdnBuf_sh   = (float**)realloc2d((void**)pData->fdnBuf_sh,
                                             nSH, nSamp, sizeof(float));
    pData->encBeamCoeffs = (float**)realloc2d((void**)pData->encBeamCoeffs,
                                              nSH, nDir, sizeof(float));
    pData->decBeamCoeffs = (float**)realloc2d((void**)pData->decBeamCoeffs,
                                              nDir, nSH, sizeof(float));
    // (re)allocate buffers that depend only on the number of output channels
    // OPTIM: These could alternatively only be reallocated when the out design
    // changes, but it's not heavy and it's more concise
    pData->dirGainBuf    = (float**)realloc2d((void**)pData->dirGainBuf,
                                              nBand, nDir, sizeof(float));
    // NOTE: omni t60 is still a 2D array nBand x 1 for compatibility with calcT60()
    pData->t60Buf_omni   = (float*)realloc1d((void*)pData->t60Buf_omni,
                                              nBand * sizeof(float));
    pData->t60Buf_dir    = (float**)realloc2d((void**)pData->t60Buf_dir,
                                              nBand, nDir, sizeof(float));
    
    pData->analysisStage = ANALYSIS_BUFS_LOADED;
}

void hosirrlib_renderTMP(
                         void  *  const hHS
                         )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("hosirrlib_renderTMP called.\n"); // DEBUG
    /* Check if processing should actually go-ahead */
    if(pData->ambiRIR_status != AMBI_RIR_STATUS_LOADED ||
       pData->lsRIR_status == LS_RIR_STATUS_RENDERED ||
       pData->lsRIR_status == LS_RIR_STATUS_RENDEREDING_ONGOING)
        return;
    else
        pData->lsRIR_status = LS_RIR_STATUS_RENDEREDING_ONGOING;
    
    strcpy(pData->progressText,"Processing");
    
    hosirrlib_processRIR(pData);
    
    /* indicate that rendering is complete */
    pData->progress0_1 = 1.0f;
    pData->lsRIR_status = LS_RIR_STATUS_RENDERED;
}


void hosirrlib_processRIR(
                          void  *  const hHS
                          )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("processRIR called.\n"); // DEBUG
    
    /* Check previous stages are complete */
    if (pData->analysisStage < ANALYSIS_BUFS_LOADED-1) { // TODO: handle fail case
        printf("processRIR called before previous stages were completed: %d\n", pData->analysisStage);
        return;
    }
    
    // TODO: SET CONSTANTS ELSEWHERE
    // lag after direct onset to start measurement
    const float directLagSec = 0.005f; // sec
    // direct arrival onset threashold (dB below omni peak)
    const float directOnsetThreshDb = -6.f;
    // diffuse onset threashold (normalized scalar below diffuseness peak)
    const float diffuseOnsetThresh = 0.707f;
    // starting level of t60 measurement (<= 0)
    const float t60_start_db = -2.f;
    // db falloff span over which t60 is measured (> 0)
    const float t60_span_db = 15.f;
    
    hosirrlib_setDirectOnsetIndex(pData, directOnsetThreshDb);
    hosirrlib_setDiffuseOnsetIndex(pData, diffuseOnsetThresh);
    hosirrlib_splitBands(pData, pData->rirBuf_sh, pData->rirBuf_bnd_sh, 1, RIR_BANDS_SPLIT);
    hosirrlib_beamformRIR(pData, pData->rirBuf_bnd_sh, pData->rirBuf_bnd_dir, BEAMFORMED);
    
    // omni, bandwise EDCs, T60s
    hosirrlib_calcEDC_omni(pData, pData->rirBuf_bnd_sh, pData->edcBufOmn_bnd,
                           pData->nBand, pData->nSamp,
                           RIR_EDC_OMNI_DONE);
    hosirrlib_calcT60_omni(pData, pData->edcBufOmn_bnd, pData->t60Buf_omni,
                           pData->nBand, pData->nSamp,
                           t60_start_db, t60_span_db,
                           // pData->directOnsetIdx + (int)(pData->fs * directLagSec) // measure after this index
                           pData->diffuseOnsetIdx, // measure from diffuseOnset
                           T60_OMNI_DONE);
    
    // directional, bandwise EDCs, T60s
    hosirrlib_calcEDC_beams(pData, pData->rirBuf_bnd_dir, pData->edcBuf_bnd_dir,
                            pData->nBand, pData->nDir, pData->nSamp,
                            RIR_EDC_DIR_DONE);
    hosirrlib_calcT60_beams(pData, pData->edcBuf_bnd_dir, pData->t60Buf_dir,
                            pData->nBand, pData->nDir, pData->nSamp,
                            t60_start_db, t60_span_db,
                            // pData->directOnsetIdx + (int)(pData->fs * directLagSec) // measure after this index
                            pData->diffuseOnsetIdx, // measure from diffuseOnset
                            T60_DIR_DONE
                            );
}


/*
 thresh_dB: threshold (pressure, dB) below the absolute max value in the
            buffer, above which the onset is considered to have occured.

 NOTE: The naive approach is to consider the index of absolute max value
        in the buffer. However, with synthetic RIRs, like from a shoebox,
        coincident reflections may have greater magnitude than the direct
        arrival. So this returns the index of the first value that crosses
        the (max abs value - threshold).
 TODO: Could return true first _peak_.
        This would require a second iteration, searching for a local max
        within, say, 2 ms after the onset index.
*/
void hosirrlib_setDirectOnsetIndex(
                                   void* const hHS,
                                   const float thresh_dB
                                   )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("setDirectOnsetIndex called.\n"); // DEBUG
    
    /* Check previous stages are complete */
    if (pData->analysisStage < DIRECT_ONSET_FOUND-1) { // TODO: handle fail case
        printf("setDirectOnsetIndes called before previous stages were completed: %d\n", pData->analysisStage);
        return;
    }
    
    const int nSamp = pData->nSamp;
    const float fs = pData->fs;
    
    float* vabs_tmp = malloc1d(nSamp * sizeof(float));
    
    /* absolute values of the omni channel */
    int maxIdx;
    utility_svabs(&pData->rirBuf_sh[0][0], nSamp, vabs_tmp); // abs(omni)
    utility_simaxv(vabs_tmp, nSamp, &maxIdx); // index of max(abs(omni))
    
    /* index of first index above threshold */
    float maxVal = vabs_tmp[maxIdx];
    float onsetThresh = maxVal * powf(10.f, thresh_dB / 20.f);
    const int directOnsetIdx = hosirrlib_firstIndexGreaterThan(vabs_tmp, 0, nSamp-1, onsetThresh);
    const float t0 = ((float)directOnsetIdx / fs) - (pData->sourceDistance / 343.f);
    
    pData->directOnsetIdx = directOnsetIdx;
    pData->t0 = t0;
    pData->t0Idx = (int)(t0 * fs); // negative value means t0 is before RIR start time
    
    // DEBUG
    printf("          t0 index: %d (%.3f sec)\n", pData->t0Idx, t0);
    printf("direct onset index: %d (%.3f sec)\n", pData->directOnsetIdx, (float)pData->directOnsetIdx / pData->fs);
    
    free(vabs_tmp);
    pData->analysisStage = DIRECT_ONSET_FOUND;
}


void hosirrlib_setDiffuseOnsetIndex(
                                    void* const hHS,
                                    const float thresh_fac
                                    )
{ /*
   * thresh_fac: threshold (normalized scalar) below the absolute max value in the
   *            buffer, above which the onset is considered to have occured.
   */
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("setDiffuseOnsetIndex called.\n"); // DEBUG
    
    /* Check previous stages are complete */
    if (pData->analysisStage < DIFFUSENESS_ONSET_FOUND-1) { // TODO: handle fail case
        printf("setDiffuseOnsetIndex called before previous stages were completed: %d\n", pData->analysisStage);
        return;
    }
    
    /* Check if processing should actually go-ahead */
    //if(pData->ambiRIR_status != AMBI_RIR_STATUS_LOADED
    // ||  pData->lsRIR_status == LS_RIR_STATUS_RENDERED ||
    // ||  pData->lsRIR_status == LS_RIR_STATUS_RENDEREDING_ONGOING
    //   ) {
    //    printf("bailing in setDiffuseOnsetIndex bc of ambiRIR/lsRIR_status.\n");
    //    return;
    //} else {
    //    pData->lsRIR_status = LS_RIR_STATUS_RENDEREDING_ONGOING;
    //}
    
    /* take a local copy of current configuration to be thread safe */
    const int   nChan   = pData->nDir;
    const int   nBand   = pData->nBand;
    const int   nPV     = 4; // working on WXYZ only
    const float fs      = pData->fs;
    const int   order   = pData->shOrder;
    const int   lSig    = pData->nSamp;
    const int   winsize = pData->windowLength;
    
    /* Inits */
    const int fftsize   = winsize*2;
    const int hopsize   = winsize/2;     /* half the window size time-resolution */
    const int nBins_anl = winsize/2 + 1; /* nBins used for analysis */
    const int nBins_syn = fftsize/2 + 1; /* nBins used for synthesis */
    const int lSig_pad  = winsize/2 + lSig + winsize*2; // winsize/2 at head, sinsize*2 at tail
    float_complex pvCOV[4][4];
    
    const int nwin_smooth = 3; // number of windows to smooth over (similar to a t60 time) TODO: Constant
    
    // Smoothing constant, a refactoring of
    //      tau = (hopsize * nwin_smooth) / fs;                         <-- smoothing coefficient
    //      smooth_const = expf( -1.f / (tau * fs / (float)hopsize));   <-- decay
    const float smooth_const = expf( -1.f / (float)nwin_smooth);
    
    printf("diffuseness smoothing coeff: %.3f vs const %.3f\n", smooth_const, ALPHA_DIFF_COEFF);
    
    /* Max freq bin to calculate diffuseness */
    float nearestVal = 10e5f;
    int   maxDiffFreq_idx = 0;
    for(int i = 0; i < nBins_anl; i++){
        // bin_idx * binwidth (Hz)
        // float tmp = fabsf((float)i * (fs / (float)winsize) - MAX_DIFF_FREQ_HZ); // TODO: should be fs/fftsize, not fs/winsize?
        float tmp = fabsf((float)i * (fs / (float)fftsize) - MAX_DIFF_FREQ_HZ);
        if(tmp < nearestVal){
            nearestVal = tmp;
            maxDiffFreq_idx = i;
        }
    }
    
    // DEBUG
    printf("num diffuseness bins: %d\n", maxDiffFreq_idx+1);
    printf("num    analysis bins: %d\n", nBins_anl);
    
    float* pvir, * pvir_pad, * win, * insig_win, * diff_win;
    float_complex* inspec_anl, * inspec_syn, * pvspec_win;
    void* hFFT_syn; // for FFT
    
    /* make local copy of current Ambi RIR, in WXYZ ordering */
    
    // NOTE: Assumes ACN-N3D
    pvir = malloc1d(nPV * lSig * sizeof(float)); // OPTIM: pvir to 2D ptr? (float**)malloc2d(nPV, lSig, sizeof(float));
    int xyzOrder[4] = {0, 3, 1, 2}; // w y z x -> w x y z
    for(int i = 0; i < nPV; i++) {
        int inChan = xyzOrder[i];
        memcpy(&pvir[i * lSig], &pData->rirBuf_sh[inChan][0], lSig * sizeof(float));
    };
    
    /* scale XYZ to normalized velocity */
    float velScale = 1.f / sqrtf(3.f);
    utility_svsmul(&pvir[1 * lSig], &velScale, (nPV-1) * lSig, &pvir[1 * lSig]);
    
    /* normalise such that peak of the omni is 1 */
    int peakIdx;
    float peakNorm;
    utility_simaxv(pvir, lSig, &peakIdx); /* index of max(abs(omni)) */
    peakNorm = 1.0f / fabsf(pvir[peakIdx]);
    utility_svsmul(pvir, &peakNorm, nPV * lSig, pvir);
    
    /* zero pad the signal's start and end for STFT */
    // winsize/2 at head, winsize*2 at tail
    pvir_pad = calloc1d(nPV * lSig_pad, sizeof(float));
    for(int i = 0; i < nPV; i++) {
        memcpy(&pvir_pad[i * lSig_pad + (winsize/2)],
               &(pvir[i * lSig]),
               lSig * sizeof(float));
    }
    
    /* transform window (symmetric Hann - 'hanning' in matlab) */
    win = malloc1d(winsize * sizeof(float));
    for(int i = 0; i < winsize; i++)
        win[i] = powf( sinf((float)i * (M_PI / (float)winsize)), 2.0f );
    
    /* mem alloc for diffuseness of each window */
    int nDiffFrames = (int)((lSig + (2 * winsize)) / hopsize + 0.5f);
    diff_win    = calloc1d(nDiffFrames, sizeof(float));
    
    // mem alloc for a single window of processing
    insig_win   = calloc1d(fftsize, sizeof(float));
    inspec_syn  = calloc1d(nPV * nBins_syn, sizeof(float_complex)); // ma->ca
    inspec_anl  = calloc1d(nPV * nBins_anl, sizeof(float_complex)); // ma->ca
    pvspec_win  = malloc1d(nPV * nBins_anl * sizeof(float_complex));
    
    saf_rfft_create(&hFFT_syn, fftsize);
    
    /* Initialize 'prev' intensity and energy vars such that initital
     * diffuseness = 0 (appropriate in context of SRIR analysis) */
    float prev_ixyz_smooth[3] = {1.f, 0.f, 0.f};
    float prev_energy_smooth = 1.f;
    
    
    /* Main processing */
    
    strcpy(pData->progressText,"HOSIRR - Rendering");
    
    /* window-hopping loop */
    int irIdx = 0; // sample frame index into pvir_pad, increments by hopsize
    int hopCount = 0;
    while (irIdx + winsize < lSig + 2 * winsize)
    {
        /* update progress */
        pData->progress0_1 = (float)irIdx / (float)(lSig + 2 * winsize);
 
        /* Transform to frequency domain, one channel at a time */
        for(int ipv = 0; ipv < nPV; ipv++){
            
            // window the input, placing it into a zero-padded buffer insig_win
            for(int j = 0; j < winsize; j++)
                insig_win[j] = win[j] * pvir_pad[(ipv * lSig_pad) + irIdx + j];
            
            // full fft, put in inspec_syn
            saf_rfft_forward(hFFT_syn, insig_win, &inspec_syn[ipv * nBins_syn]);
            
            for(int j = 0, k = 0; j < nBins_anl; j++, k += fftsize/winsize)
                inspec_anl[(ipv * nBins_anl) + j] = inspec_syn[(ipv * nBins_syn) + k];
        }
        
        // copy spectrum of windowed input channels into wxyzspec_win (nPV * nBins_anl)
        for(int i = 0; i < nPV; i++)
            memcpy(&pvspec_win[i * nBins_anl],
                   &inspec_anl[i * nBins_anl],
                   nBins_anl * sizeof(float_complex));
                    
        /* Compute broadband covariance matrix */
        const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
                    4, 4, maxDiffFreq_idx+1, &calpha,
                    pvspec_win, nBins_anl,
                    pvspec_win, nBins_anl, &cbeta,
                    FLATTEN2D(pvCOV), 4);
                
        /* Compute broadband intensity and energy */
        float intensity[3], energy;
        float ixyz_smooth[3], energy_smooth, iaMag_smooth;
        
        for(int i = 0; i < 3; i++)
            intensity[i] = crealf(pvCOV[1+i][0]);
        energy = 0.0f;
        for(int i = 0; i < 4; i++)
            energy += crealf(pvCOV[i][i])*0.5f;
        
        // DEBUG
        //printf("intensity %.3f, %.3f, %.3f\n", intensity[0], intensity[1], intensity[2]);
        
        /* Estimating and time averaging of broadband diffuseness */
        iaMag_smooth = 0.0f;
        for(int i = 0; i < 3; i++) {
            // TODO: Constant ALPHA_DIFF_COEFF is ignored
            // ixyz_smooth[i] = ((1.0f-ALPHA_DIFF_COEFF) * intensity[i]) + (ALPHA_DIFF_COEFF * prev_ixyz_smooth[i]);
            ixyz_smooth[i] = ((1.0f-smooth_const) * intensity[i]) + (smooth_const * prev_ixyz_smooth[i]);
            prev_ixyz_smooth[i] = ixyz_smooth[i];
            iaMag_smooth += powf( ixyz_smooth[i], 2.0f ); // NOTE: ixyz is real so fabsf() was removed
        }
        iaMag_smooth = sqrtf(iaMag_smooth);
        // energy_smooth = ((1.0f-ALPHA_DIFF_COEFF) * energy) + (ALPHA_DIFF_COEFF * prev_energy_smooth);
        energy_smooth = ((1.0f-smooth_const) * energy) + (smooth_const * prev_energy_smooth);
        prev_energy_smooth = energy_smooth;
        
        // store broadband diffuseness value
        diff_win[hopCount] = 1.0f - (iaMag_smooth / (energy_smooth + 2.23e-10f));
        
        // DEBUG
        //printf("normIntensity_smooth %.3f\n", intensityMag_smooth);
        //printf("diff %.3f\n", diff_win[hopCount]);
        
        // advance
        irIdx += hopsize;
        hopCount++;
    }
    
    /* diffuse onset */
    
    /* index of max(diffuseness) */
    // begin search for max diffuseness at the first window containing the
    // direct arrival, where we can expect low diffuseness
    int hopIdx_direct, maxWinIdx;
    int directIdx_tmp = pData->directOnsetIdx - winsize;
    
    if (directIdx_tmp <= 0) {
        hopIdx_direct = 0;
    } else {
        hopIdx_direct = (int)ceilf((float)directIdx_tmp / hopsize);
    };
    utility_simaxv(&diff_win[hopIdx_direct],
                   nDiffFrames-hopIdx_direct, &maxWinIdx);
    maxWinIdx += hopIdx_direct; // add the starting point back
    
    /* index of first index above threshold */
    const float maxVal      = diff_win[maxWinIdx];
    const float onsetThresh = maxVal * thresh_fac;
    const int   onsetWinIdx = hosirrlib_firstIndexGreaterThan(diff_win, 0, nDiffFrames-1, onsetThresh);
    
    // NOTE: Normally we'd assume the _sample_ onset index is in the middle
    // of the window. I.e. onsetWinIdx * hopsize + winsize/2.
    // But because pvsig was zero-padded at the head by winsize/2, we can omit
    // that offset here.
    float diffuseOnsetIdx = onsetWinIdx * hopsize;
    
    if (maxVal < pData->diffuseMin) {
        // TODO: better handling of diffuseness fail
        hosirr_print_warning("Diffuseness didn't exceed the valid threshold.\nFalling back on directOnset +5ms.");
        printf("diffuseness: %.2f < %.2f\n", maxVal, pData->diffuseMin);
        diffuseOnsetIdx = pData->directOnsetIdx + (int)(0.005 * fs);
    }
    
    pData->diffuseOnsetIdx = diffuseOnsetIdx; // diffuse onset in the sample buffer
    pData->diffuseOnsetSec = ((float)diffuseOnsetIdx / fs) - pData->t0; // diffuse onset time in the *room*
    
    // DEBUG
    printf("     begin search at hop idx: %d\n", hopIdx_direct);
    printf("       diffuse onset win idx: %d\n", onsetWinIdx);
    printf("    diffuse onset sample idx: %d (%.3f sec)\n", pData->diffuseOnsetIdx, (float)pData->diffuseOnsetIdx / pData->fs);
    printf(" t0-adjusted diff onset time: %.3f sec\n", pData->diffuseOnsetSec);
    printf("         diffuse max win idx: %d\n", maxWinIdx);
    printf("           diffuse max value: %.2f\n", maxVal);
    printf("        diffuse onset thresh: %.2f\n", onsetThresh);
    
    /* Sanity checks */
    // check that time events are properly increasing
    if ((pData->t0Idx > pData->directOnsetIdx) ||
        (pData->directOnsetIdx > pData->diffuseOnsetIdx))
        hosirr_print_error("Order of analyzed events isn't valid. Should be t0Idx < directOnsetIdx < diffuseOnsetIdx.");
    // check that diffuse onset (from t0) is positive
    if (pData->diffuseOnsetSec <= 0)
        hosirr_print_error("Diffuse onset is negative");
    
    pData->analysisStage = DIFFUSENESS_ONSET_FOUND;

    // Cleanup
    free(pvir);
    free(pvir_pad);
    free(win);
    free(diff_win);
    free(insig_win);
    free(inspec_anl);
    free(inspec_syn);
    free(pvspec_win);
    saf_rfft_destroy(&hFFT_syn);
}

/* Freq-domain band filtering via FD convolution */
void hosirrlib_splitBands(
                          void*    const hHS,
                          float**  const inBuf,
                          float*** const bndBuf,
                          int removeFiltDelayFLAG,
                          ANALYSIS_STAGE thisStage
                          )
/*
 removeFiltDelay :  true will truncate the head and tail of the filtered signal,
                    instead of (only) the tail. So either way, the written
                    output will always be pData->nSamp for each channel.
 */
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("splitBands called.\n"); // DEBUG
    
    /* Check previous stages are complete */
    if (pData->analysisStage < thisStage-1) { // TODO: handle fail case
        printf("splitBands called before previous stages were completed: %d\n", pData->analysisStage);
        return;
    }
    
    const int nInSamp   = pData->nSamp; // TODO: assumes nsamp of outbuf == nsamp of RIR
    const int filtOrder = pData->bandFiltOrder;
    
    int startCopyIdx;
    if (removeFiltDelayFLAG) {
        startCopyIdx = (int)(filtOrder / 2.f); // length of filter delay (floor)
    } else {
        startCopyIdx = 0;
    }
    
    /* Apply filterbank to rir_bands.
        Because of fftconv's indexing, it expects filter coefficients
        for each channel. So we do one band/channel at a time.
        The output length will be nSamp + filtOrder.
        See fftconv(): y_len = x_len + h_len - 1;   */
    float* stage = malloc1d((nInSamp + filtOrder) * sizeof(float));

    for(int ish = 0; ish < pData->nSH; ish++) {
        for(int ib = 0; ib < pData->nBand; ib++) {
            fftconv(&inBuf[ish][0],             // input: 1 channel at a time
                    &pData->H_bandFilt[ib][0],  // band filter coeffs
                    nInSamp,                    // input length
                    filtOrder + 1,              // filter length
                    1,                          // 1 channel at a time
                    stage);                     // write to staging buffer

            /* Copy staging buffer to out */
            memcpy((void*)&bndBuf[ib][ish][0],
                   &stage[startCopyIdx],
                   nInSamp * sizeof(float));
        }
    }
    
    pData->analysisStage = thisStage;
    free(stage);
}

    
void hosirrlib_beamformRIR(
                           void*    const hHS,
                           float*** const inBuf,
                           float*** const beamBuf,
                           ANALYSIS_STAGE thisStage)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("beamformRIR called.\n"); // DEBUG
    
    /* Check previous stages are complete */
    if (pData->analysisStage < thisStage-1) { // TODO: handle fail case
        printf("beamformRIR called before previous stages were completed: %d\n", pData->analysisStage);
        return;
    }
    
    const int nDir    = pData->nDir;
    const int nSamp   = pData->nSamp;
    const int nBand   = pData->nBand;
    const int nSH     = pData->nSH;
    const int shOrder = pData->shOrder;
    
    float * c_l = malloc1d((shOrder + 1) * sizeof(float)); // z beam coeffs, only shOrder + 1 are used
    
    /* Calculate beamforming coeffients */
    // OPTIM: These beamforming coeffs only need to be updated if nSH (input)
    //        or the spherical design (output) changes, so could be refactored
    for (int idir = 0; idir < nDir; idir++) {
        switch(pData->beamType){
            case CARDIOID_BEAM:
                beamWeightsCardioid2Spherical(shOrder, c_l); break;
            case HYPERCARDIOID_BEAM:
                beamWeightsHypercardioid2Spherical(shOrder, c_l); break;
            case MAX_EV_BEAM:
                beamWeightsMaxEV(shOrder, c_l); break;
        }
        rotateAxisCoeffsReal(shOrder,
                             (float*)c_l,
                             (SAF_PI / 2.0f) - (pData->loudpkrs_dirs_deg[idir][1] * SAF_PI / 180.0f),
                             pData->loudpkrs_dirs_deg[idir][0] * SAF_PI / 180.0f,
                             &pData->decBeamCoeffs[idir][0]);
    }
    
    /* Apply beam weights
        (ref. beamformer_process())
        A: float** encBeamCoeffs;  // nDir  x nSH
        B: float*** rirBuf_bnd_sh;  // nBand x nSH  x nSamp
        C: float*** rirBuf_bnd_dir;  // nBand x nDir x nSamp  */
    for (int bd = 0; bd < nBand; bd++) {
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    nDir, nSamp, nSH, 1.0f,
                    (const float*)&pData->decBeamCoeffs[0][0], nSH,         // A, 1st dim of A (row-major)
                    (const float*)&pData->rirBuf_bnd_sh[bd][0][0], nSamp,    // B, 1st dim of B
                    0.0f,                                                   // beta scalar for C
                    (float*)&beamBuf[bd][0][0], nSamp);                     // C, 1st dim of C
    }
    
    pData->analysisStage = thisStage;
    free(c_l);
}


/* Calc bandwise EDCs each directional beam */
void hosirrlib_calcEDC_beams(
                       void*    const hHS,
                       float*** const inBuf,
                       float*** const edcBuf,
                       const int nBand,
                       const int nDir,
                       const int nSamp,
                       ANALYSIS_STAGE thisStage)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("calcEDC_beams called.\n"); // DEBUG
    
    /* Check previous stages are complete */
    if (pData->analysisStage < thisStage-1) { // TODO: handle fail case
        printf("calcEDC_beams called before previous stages were completed: %d\n", pData->analysisStage);
        return;
    }
    
    /* Copy the RIR band beams into EDC buffer */
    utility_svvcopy(FLATTEN3D(inBuf),
                    nBand * nDir * nSamp,
                    FLATTEN3D(edcBuf));

    /* EDC, one band/channel at a time */
    for (int ib = 0; ib < nBand; ib++) {
        for (int ich = 0; ich < nDir; ich++) {
            hosirrlib_calcEDC_1ch(&edcBuf[ib][ich][0], nSamp);
        }
    }
    pData->analysisStage = thisStage;
}


/* Calc bandwise EDCs of the omni channel */
void hosirrlib_calcEDC_omni(
                            void*    const hHS,
                            float*** const shInBuf, // nband x nsh x nsamp
                            float**  const edcBuf_omn,  // nband x nsamp
                            const int nBand,
                            const int nSamp,
                            ANALYSIS_STAGE thisStage)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("calcEDC_omni called.\n"); // DEBUG
    
    /* Check previous stages are complete */
    if (pData->analysisStage < thisStage-1) { // TODO: handle fail case
        printf("calcEDC_omni called before previous stages were completed: %d\n", pData->analysisStage);
        return;
    }
    
    /* Copy the RIR's omni bands into EDC buffer */
    for (int ib = 0; ib < nBand; ib++) {
        utility_svvcopy(&shInBuf[ib][0][0], nSamp, &edcBuf_omn[ib][0]);
    }
    
    /* Calculate EDC in-place, one band at a time */
    for (int ib = 0; ib < nBand; ib++) {
        hosirrlib_calcEDC_1ch(&edcBuf_omn[ib][0], nSamp);
    }
    pData->analysisStage = thisStage;
}


/* Calc EDC of a single channel
 Note: in-place operation, so copy data into edcBuf before calling.
 OPTIM: vectorize
 */
void hosirrlib_calcEDC_1ch(
                           float* const edcBuf,
                           const int nSamp)
{
    double sum = 0.0;
    for (int i = nSamp - 1; i >= 0; i--) {
        sum += edcBuf[i] * edcBuf[i];           // energy
        edcBuf[i] = (float)(10.0 * log10(sum)); // store in dB
    }
}


/* Calc T60s of bandwise beams
 start_db : measure the T60 starting at this level of decay
            (after beginIdx), specify as negative (<= 0)
 span_db  : measure the T60 over this dB decay span (specify as positive)
 beginIdx : start the search for measurement bounds at this index onward
            e.g. a sample index after the first arrival
            0 for the beginning of the EDC channel
 */
void hosirrlib_calcT60_beams(
                       void* const hHS,
                       float*** const edcBuf,
                       float** const t60Buf,
                       const int nBand,
                       const int nChan,
                       const int nSamp,
                       const float start_db,
                       const float span_db,
                       const int beginIdx,
                       ANALYSIS_STAGE thisStage)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("calcT60_beams called.\n"); // DEBUG
    
    /* Check previous stages are complete */
    if (pData->analysisStage < thisStage-1) { // TODO: handle fail case
        printf("calcT60_beams called before previous stages were completed: %d\n", pData->analysisStage);
        return;
    }
    
    const int   nSH = pData->nSH;
    const float fs  = pData->fs;
    
    float* const x_slope1 = malloc1d(nSamp * sizeof(float)); // vector with slope of 1/samp
    float* const y_edc0m  = malloc1d(nSamp * sizeof(float)); // zero-mean EDC
    float* const stage    = malloc1d(nSamp * sizeof(float)); // staging buffer
    
    // start and end samples to measure t60 in each channel
    int*** const st_end_meas = (int***)malloc3d(nBand, nChan, 2, sizeof(int));

    /* find the start and end points of the measurement */
    for (int ibnd = 0; ibnd < nBand; ibnd++) {
        for (int ich = 0; ich < nChan; ich++) {
            hosirrlib_findT60_bounds(&edcBuf[ibnd][ich][0],
                                     beginIdx, nSamp, start_db, span_db,
                                     &st_end_meas[ibnd][ich][0] // 2 x 1
                                     );
        }
    }
    /* Measure the t60 by the line of best fit */
    for (int ibnd = 0; ibnd < nBand; ibnd++) {
        for (int ich = 0; ich < nChan; ich++) {
            hosirrlib_t60_lineFit(&edcBuf[ibnd][ich][0], // omni ch = 0
                                  x_slope1, y_edc0m, stage,
                                  st_end_meas[ibnd][ich][0], st_end_meas[ibnd][ich][1],
                                  pData->fs, &t60Buf[ibnd][ich]);
            // DEBUG
            printf("t60: dir %d band %d  %.2f sec\n", ich, ibnd, t60Buf[ibnd][ich]);
        }
    }
        
    free(x_slope1);
    free(y_edc0m);
    free(stage);
    free(st_end_meas);
    pData->analysisStage = thisStage;
}

/* Calc T60s of bandwise omni channel
 See
 */
void hosirrlib_calcT60_omni(
                       void*   const hHS,
                       float** const edcBuf_omn,
                       float*  const t60Buf,
                       const int nBand,
                       const int nSamp,
                       const float start_db,
                       const float span_db,
                       const int beginIdx,
                       ANALYSIS_STAGE thisStage)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("calcT60_omni called.\n"); // DEBUG
    
    /* Check previous stages are complete */
    if (pData->analysisStage < thisStage-1) { // TODO: handle fail case
        printf("calcT60_omni called before previous stages were completed: %d\n", pData->analysisStage);
        return;
    }
    
    float* const x_slope1 = malloc1d(nSamp * sizeof(float)); // vector with slope of 1/samp
    float* const y_edc0m  = malloc1d(nSamp * sizeof(float)); // zero-mean EDC
    float* const stage    = malloc1d(nSamp * sizeof(float)); // staging buffer
    
    // start and end samples to measure t60 in each channel
    int** const st_end_meas = (int**)malloc2d(nBand, 2, sizeof(int));

    /* find the start and end points of the measurement */
    for (int bd = 0; bd < nBand; bd++) {
        hosirrlib_findT60_bounds(&edcBuf_omn[bd][0], // omni ch = 0
                                 beginIdx, nSamp, start_db, span_db,
                                 &st_end_meas[bd][0]);
    }
    
    for (int ibnd = 0; ibnd < nBand; ibnd++) {
        hosirrlib_t60_lineFit(&edcBuf_omn[ibnd][0], // omni ch = 0
                              x_slope1, y_edc0m, stage,
                              st_end_meas[ibnd][0], st_end_meas[ibnd][1],
                              pData->fs, &t60Buf[ibnd]);
        // DEBUG
        printf("t60: dir %d band %d  %.2f sec\n", 0, ibnd, t60Buf[ibnd]);
    }
    
    free(x_slope1);
    free(y_edc0m);
    free(stage);
    free(st_end_meas);
    pData->analysisStage = thisStage;
}


/* Measure the t60 by the line of best fit
    // x0m : zero-mean vector of a line with a slope of 1/sample
    // y0m : vector of edc values (over measurement span), with mean removed
    x0m      = vec_idc - mean(vec_idc);
    y0m      = edc_span - mean(edc_span);
    dc_db    = sum(x0m .* y0m) / sum(x0m.^2); // decay constant (dB/sample)
    t60_meas = (-60 / dc_db) / fs;
 */
void hosirrlib_t60_lineFit(
                           float* const edcBuf,
                           float* x_slopeBuf,
                           float* y_edc0mBuf,
                           float* stageBuf,
                           const int st_meas,
                           const int end_meas,
                           const float fs,
                           float* const t60Buf_wrPtr)
{
    int nSamp_meas = end_meas - st_meas + 1;
    float y_mean = sumf(&edcBuf[st_meas], nSamp_meas) / nSamp_meas;
    
    /* Construct a vector with a slope of 1:samp with zero mean */
    float firstVal = (nSamp_meas - 1) / -2.f; // first value
    for (int i = 0; i < nSamp_meas; i++) {
        x_slopeBuf[i] = firstVal + i;
    };
    // remove mean from EDC, within the measurement span
    utility_svssub(&edcBuf[st_meas], &y_mean, nSamp_meas, y_edc0mBuf);
    
    // covariances: x*y and x*x
    float c_xy, c_xx, dbPerSamp;
    utility_svvmul(x_slopeBuf, y_edc0mBuf, nSamp_meas, stageBuf);
    c_xy = sumf(stageBuf, nSamp_meas);
    utility_svvmul(x_slopeBuf, x_slopeBuf, nSamp_meas, stageBuf);
    c_xx = sumf(stageBuf, nSamp_meas);
    dbPerSamp = c_xy / c_xx; // slope
    
    *t60Buf_wrPtr = (-60.0f / dbPerSamp) / fs; // write out
    //printf("dbPerSamp: %.7f\n", dbPerSamp); // DEBUG
}

void hosirrlib_findT60_bounds(
                              float* const edcBuf, // nsamp x 1
                              const int beginIdx,
                              const int nSamp,
                              const float start_db,
                              const float span_db,
                              int* const st_end_meas) // 2 x 1
{
    /* Check start and end points of T60 measurement */
    float edcMax = edcBuf[beginIdx];
    
    int start_t60 = hosirrlib_firstIndexLessThan(edcBuf, beginIdx, nSamp - 1,
                                                 edcMax + start_db); // start_db is negative
    if (start_t60 < 0) {
        // TODO: warn or error
        // No value found below the start level, fall back to the head of the EDC
        st_end_meas[0] = 0;
    } else {
        st_end_meas[0] = start_t60;
    }
    
    int end_t60 = hosirrlib_firstIndexLessThan(edcBuf, start_t60, nSamp - 1,
                                               edcMax + start_db - span_db); // start_db is negative, span_db is positive
    if (end_t60 < 0) {
        // TODO: warn or error
        // No value found below start_db - span_db, fall back to near the
        // end of the EDC.. this will likely be quite inaccurate for high frequency bands!
        st_end_meas[1] = (int)(0.7f * nSamp);
    } else {
        st_end_meas[1] = beginIdx + end_t60;
    }
    
    // DEBUG
    //printf("edcMax: %.4f\n", edcMax);
    //printf("start idx: %d\n", st_end_meas[bd][ch][0]);
    //printf("\tend idx: %d\n", st_end_meas[bd][ch][1]);
}

// Returns -1 on fail
int hosirrlib_firstIndexLessThan(
                                 float* vec,
                                 int startIdx,
                                 int endIdx,
                                 float thresh)
{
    for (int i = startIdx; i < endIdx+1; i++) {
        if (vec[i] < thresh)
            return i;
    }
    return -1;
}


// Returns -1 on fail
int hosirrlib_firstIndexGreaterThan(
                                    float* vec,
                                    int startIdx,
                                    int endIdx,
                                    float thresh)
{
    for (int i = startIdx; i < endIdx+1; i++) {
        if (vec[i] > thresh)
            return i;
    }
    return -1;
}


// For the GUI to display the directional EDCs
// might modify to view bands for debugging
void hosirrlib_copyNormalizedEDCs_dir(
                                     void* const hHS,
                                     float** edcOut, // ndir x nsamp
                                     float displayRange)
{
    /*
     Note edcBuf_bnd_dir are foat*** nband x ndir x nsamp, but dirEDC pointer
     in the UI is float** ndir x nsamp,
     so for now, just return the first band of each direction
     */
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    const int nBand = pData->nBand;
    const int nSamp = pData->nSamp;
    const int nDir  = pData->nDir;
    float*** const edcIn = pData->edcBuf_bnd_dir;
    
    if (pData->analysisStage >= RIR_EDC_DIR_DONE) { // ensure EDCs are rendered
        
        /* normalise to range [-1 1] for plotting */
        float maxVal, minVal, range, add, scale, sub;

        maxVal = edcIn[0][0][0];        // intializse to first value of first channel
        
        int maxBndIdx = 0, maxDirIdx = 0; // DEBUG vars
        for (int id = 0; id < nDir; id++) {
            for (int ib = 0; ib < nBand; ib++) {
                float val0 = edcIn[ib][id][0]; // first edc value in each channel
                if (val0 > maxVal) {
                    maxVal = edcIn[ib][id][0];
                    maxBndIdx = ib; maxDirIdx = id; // DEBUG vars
                }
            }
        }
        // DEBUG
        printf("\n>> max val found on dir %d, bnd %d: %.2f\n\n", maxDirIdx, maxBndIdx, maxVal);
        
        minVal = maxVal - displayRange; // just display the uper displayRange in dB
        // check the first and last values of every channel
        for (int ib = 1; ib < nBand; ib++)
            for (int ich = 0; ich < nDir; ich++)
                if (edcIn[ib][ich][0] > maxVal)
                    maxVal = edcIn[ib][ich][0];

        range = maxVal - minVal;
        add   = minVal * -1.f;
        scale = 2.0f/fabsf(range);
        sub   = 1.f;
        
        // DEBUG
        //printf("max %.1f, min %.1f, rng %.1f, add %.1f, scl %.1f, sub %.1f, ",
        //       maxVal, minVal, range, add, scale, sub);
        
        printf("TEMP: viewing band channels of a single direction"); // DEBUG
        for(int i = 0; i < nDir; i++) {
            //int bndIdx = 0; // for just lowest band of all directions
            //int chIdx = i;
            int bndIdx = i % nBand; // cycle through the bands of the chIdx
            int chIdx = 0; // idx 8 for single directional decaying pw test

            utility_svsadd(&(edcIn[bndIdx][chIdx][0]), // [bnd][ch][smp]
                           &add, nSamp,
                           &edcOut[i][0]); // [ch][smp]
            utility_svsmul(&edcOut[i][0],  // [ch][smp]
                           &scale, nSamp,
                           &edcOut[i][0]); // [ch][smp]
            utility_svssub(&edcOut[i][0],  // [ch][smp]
                           &sub, nSamp,
                           &edcOut[i][0]); // [ch][smp]
        }
        
        
        /* TESTS: writing out different buffers for inspection */
        
//        // Write out band-filtered signals
//        for(int i = 0; i < pData->nDir; i++) {
//            memcpy(&edcOut[i][0],                    // copy-to channel
//                   &(pData->rirBuf_bnd_sh[i%nBand][0][0]), // [bnd][ch][smp]
//                   nSamp * sizeof(float));
//        }

//        // write out EDCs
//        for(int i = 0; i < pData->nDir; i++) {
//            //            int bndIdx = 0; // for just lowest band of all directions
//            //            int chIdx = i;
//            printf("TEMP: Writing out EDCs for one direction, by band.\n");
//            int bndIdx = i % nBand; // cycle through the bands of the chIdx
//            int chIdx = 0;
//            memcpy(&edcOut[i][0],             // copy-to channel
//                   &edcIn[bndIdx][chIdx][0],
//                   nSamp * sizeof(float));
//        }

//        // write out filter coeffs
//        for(int id = 0; id < pData->nDir; id++) {
//            memcpy(&edcOut[id][0], // write out
//                   &pData->H_bandFilt[id % pData->nBand][0],    // band filter coeffs
//                   (pData->bandFiltOrder + 1) * sizeof(float)   // filter length
//                   );
//        }
        
//        // copy input to output
//        for(int i = 0; i < pData->nDir; i++)
//            memcpy(&edcOut[i][0],             // copy-to channel
//                   &(pData->rirBuf_sh[i][0]),    // [nsh][smp]
//                   nSamp * sizeof(float));

//        // THIS WORKS: write out DC values
//        float** dcs = (float**)calloc2d(nDir, nSamp, sizeof(float));
//        for(int id = 0; id < nDir; id++) {
//            float addThis = 1.f/8 * (id%8);
//            printf("add this: %.2f\n", addThis);
//            utility_svsadd(&dcs[id][0], &addThis, nSamp, &edcOut[id][0]);
//        }
//        free(dcs);
        
//        // run input channels through band filters and write out
//        for(int id = 0; id < pData->nDir; id++) {
//            fftfilt(&pData->rirBuf_sh[id][0],             // input: 1 channel at a time
//                    &pData->H_bandFilt[id % pData->nBand][0],  // band filter coeffs
//                    pData->nSamp,                    // input length
//                    pData->bandFiltOrder + 1,              // filter length
//                    1,                          // 1 channel at a time
//                    &edcOut[id][0]);       // write out
//        }
        
    } else {
        hosirr_print_error("copyNormalizedEDCBufs_dir: EDC hasn't been rendered so can't copy it out.");
    }
}


// For the GUI to display the omni EDCs by band (repeated up to ndir currently for debugging)
void hosirrlib_copyNormalizedEDCs_omni(
                                     void* const hHS,
                                     float** edcOut, // ndir x nsamp
                                     float displayRange)
{
    /*
     Note edcBuf_bnd_dir are foat*** nband x ndir x nsamp, but dirEDC pointer
     in the UI is float** ndir x nsamp,
     so for now, just return the first band of each direction
     */
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    const int nBand = pData->nBand;
    const int nSamp = pData->nSamp;
    const int nDir  = pData->nDir;
    float** const edcIn = pData->edcBufOmn_bnd;
    
    if (pData->analysisStage >= RIR_EDC_DIR_DONE) { // ensure EDCs are rendered
        
        /* normalise to range [-1 1] for plotting */
        float maxVal, minVal, range, add, scale, sub;

        maxVal = edcIn[0][0];        // intializse to first value of first channel
        
        int maxBndIdx = 0; // DEBUG vars
        for (int ib = 0; ib < nBand; ib++) {
            float val0 = edcIn[ib][0]; // first edc value in each channel
            if (val0 > maxVal) {
                maxVal = edcIn[ib][0];
                maxBndIdx = ib; // DEBUG vars
            }
        }
        // DEBUG
        printf("\n>> max val found on omni bnd %d: %.2f\n\n", maxBndIdx, maxVal);
        
        minVal = maxVal - displayRange; // just display the uper displayRange in dB
        // check the first and last values of every channel
        for (int ib = 1; ib < nBand; ib++)
            if (edcIn[ib][0] > maxVal)
                maxVal = edcIn[ib][0];

        range = maxVal - minVal;
        add   = minVal * -1.f;
        scale = 2.0f/fabsf(range);
        sub   = 1.f;
        
        // DEBUG
        //printf("max %.1f, min %.1f, rng %.1f, add %.1f, scl %.1f, sub %.1f, ",
        //       maxVal, minVal, range, add, scale, sub);
        
        printf("TEMP: viewing band channels of a single direction"); // DEBUG
        for(int i = 0; i < nDir; i++) {
            //int bndIdx = 0; // for just lowest band of all directions
            int bndIdx = i % nBand; // cycle through the bands of the chIdx

            utility_svsadd(&(edcIn[bndIdx][0]), // [bnd][ch][smp]
                           &add, nSamp,
                           &edcOut[i][0]); // [ch][smp]
            utility_svsmul(&edcOut[i][0],  // [ch][smp]
                           &scale, nSamp,
                           &edcOut[i][0]); // [ch][smp]
            utility_svssub(&edcOut[i][0],  // [ch][smp]
                           &sub, nSamp,
                           &edcOut[i][0]); // [ch][smp]
        }
    } else {
        hosirr_print_error("copyNormalizedEDCBufs_omni: EDC hasn't been rendered so can't copy it out.");
    }
}


/*  ~~~~~~~~~~~~~~~~ SCRATCH ~~~~~~~~~~~~~~~~ */

//    // TEST: create coefficients from simple spectrum and pass those
//    float y_len = nInSamp + pData->bandFiltOrder+1 - 1;
//    int fftSize =  (int)((float)nextpow2(y_len)+0.5f);
//    int nBins = fftSize/2+1;
//    float* h0 = calloc1d(fftSize, sizeof(float));
//    float_complex* H = malloc1d(nBins * sizeof(float_complex));
//
//    // spectrum ramp DC->nyq 0->1
//    float magstep = 1.f / (float)nBins;
//    for (int ib = 0; ib < nBins; ib++)
//        H[ib] = cmplxf(magstep * ib, 0.f);

//    for(int ish = 0; ish < pData->nSH; ish++) {
//        for(int ibd = 0; ibd < pData->nBand; ibd++) {
//            fftfilt(&inBuf[ish][0],             // input: 1 channel at a time
//                    &pData->H_bandFilt[ibd][0], // band filter coeffs
//                    nInSamp,                    // input length
//                    filtOrder + 1,              // filter length
//                    1,                          // 1 channel at a time
//                    &bndBuf[ibd][ish][0]);      // write out
//        }
//    }
//    free(h0);
//    free(H);

// DEBUG FUNCTION
//void hosirrlib_inspectFilts(void* const hHS)
//{
//    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
//    printf("inspecting filters.\n");
//
//    int filt_len = 50;
//    int fft_len = 128;
//    int fftSize =  (int)((float)nextpow2(fft_len)+0.5f);
//    int nBins = fftSize/2+1;
//
//    //    float* h = &pData->H_bandFilt[0][0]; // just band 1, low pass, for now
//
//    float* h  = calloc1d(filt_len, sizeof(float)); // filter coeffs - zeros
//    h[filt_len/2] = 1.f; // filter is a centered impulse
//
//    float* h0 = calloc1d(fftSize, sizeof(float));  // time-domain filter (zero-padded)
//    float* y0 = calloc1d(nBins,   sizeof(float));  // output buffer for filter mags
//    float_complex* H = malloc1d(nBins * sizeof(float_complex)); // freq domain filter
//    void* hfft;
//    saf_rfft_create(&hfft, fftSize);
//
//    /* copy the time-domain filter kernel into the head of the
//       zero-padding buffer, prior to fft */
//    memcpy(h0, &h, filt_len * sizeof(float));
//
//    // h0 -> H: time-domain filter to freq-domain
//    saf_rfft_forward(hfft, h0, H);
//
// TEST: directly fill spectrum with a ramp of values
//       with magnitudes 0->1 to just observe the filtering works.
//    float magstep = 1.f / (float)nBins;
//    for (int ib = 0; ib < nBins; ib++)
//        H[ib] = cmplxf(0.707f, 0.707f) * (ib * magstep);
//
//    // write out magnitude of the filters
//    for (int ib = 0; ib < nBins; ib++) {
//        y0[ib] = cabsf(H[ib]);
//        printf("re %.6f, imag %.6f\n", crealf(H[ib]), cimagf(H[ib]));
//        printf("\tmag %d,  %.4f\n", ib, y0[ib]);
//    }
//
//    // cleanup
//    saf_rfft_destroy(&hfft);
//    free(h0);
//    free(y0);
//    free(H);
//}
//    /* TEST: bypass beam weighting, just copy W into all "beams" */
//    for (int id = 0; id < nDir; id++) {
//        for (int bd = 0; bd < nBand; bd++) {
//            utility_svvcopy(&inBuf[bd][0][0], // W
//                            nSamp,
//                            &beamBuf[bd][id][0]);
//        }
//    }
    
/*  ~~~~~~~~~~~~~~~~ END SCRATCH ~~~~~~~~~~~~~~~~ */


/* Render */

void hosirrlib_render
    (void  *  const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    /* The leaning tower of variables */
    int i, j, k, n, ch, maxInd, dirwinsize, BB1stPeak, order, nSH, winsize, nLS;
    int N_gtable, N_tri, aziRes, elevRes, N_azi, aziIndex, elevIndex, idx3d;
    int numSec, order_sec, nSH_sec, delay, frameCount;
    int fftsize, hopsize, nBins_anl, nBins_syn, maxDiffFreq_idx, rir_len;
    int o[HOSIRR_MAX_SH_ORDER+2];
    int idx, jl, lSig, lSig_pad;
    float fs, wetDry, peakNorm, normSec, nearestVal, tmp, a2eNorm;
    float intensity[3], ixyz_smoothed[3], energy, energy_smoothed, normSecIntensity_smoothed;
    //float t60[6] = {0.07f, 0.07f, 0.06f, 0.04f, 0.02f, 0.01f};
    float t60[6] = {0.2f, 0.2f, 0.16f, 0.12f, 0.09f, 0.04};
    float fc[6] = {125.0f, 250.0f, 500.0f, 1000.0f, 2000.0f, 4000.0f};
    float IntensityBB[3] = {0.0f};
    float IntensityBB_XYZ[3];
    float* shir, *shir_tmp, *shir_tmp2, *direct_win, *shir_direct, *shir_pad, *gtable, *sec_dirs_deg, *sectorCoeffs_tmp;
    float* lsir_ndiff, *lsir_diff, *win, *Y_enc_tmp, *D_ls_tmp;
    float* prev_energy_smoothed, *prev_ixyz_smooth, *azim, *elev, *diffs;
    float* insig_win, *ndiffs_sqrt, *lsir_win, *M_ifft, *M_ifft_fl, *rir_filt;
    float_complex* A_xyz, *sectorCoeffs, *sectorCoeffs_syn, *Y_enc, *D_ls;
    float_complex* inspec_syn, *inspec_anl, *s_anl, *WXYZ_sec, *z_diff, *z_00;
    float_complex* ndiffgains_interp, *diffgains_interp, *ndiffgains, *diffs_sqrt;
    float_complex* outspec_ndiff, *outspec_diff, *a_diff;
    float_complex pvCOV[4][4];
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    void* hFFT_syn, *hFFT_anl;
    NORMALIZATION_TYPES norm;
    CH_ORDERING chOrdering;
    
    /* Check if processing should actually go-ahead */
    if(pData->ambiRIR_status != AMBI_RIR_STATUS_LOADED ||
       pData->lsRIR_status == LS_RIR_STATUS_RENDERED ||
       pData->lsRIR_status == LS_RIR_STATUS_RENDEREDING_ONGOING)
        return;
    else
        pData->lsRIR_status = LS_RIR_STATUS_RENDEREDING_ONGOING;
    
    /* take a local copy of current configuration to be thread safe */
    fs = pData->ambiRIRsampleRate;
    order = HOSIRR_MIN(pData->analysisOrder, pData->ambiRIRorder);
    nSH = (order+1)*(order+1);
    winsize = pData->windowLength;
    lSig = pData->ambiRIRlength_samples;
    BB1stPeak = pData->broadBandFirstPeakFLAG;
    wetDry = pData->wetDryBalance;
    nLS = pData->nLoudpkrs;
    norm = pData->norm;
    chOrdering = pData->chOrdering;
    
    /* make local copy of current Ambi RIR */
    shir = malloc1d(nSH * lSig * sizeof(float));
    switch(chOrdering){
        case ACN_ORDER:
//            memcpy(shir, pData->shir, nSH * lSig * sizeof(float));
            memcpy(shir, FLATTEN2D(pData->shir), nSH * lSig * sizeof(float));
            break;
        case FUMA_ORDER:
            /* only for first-order, convert to ACN */
            assert(nSH==4);
            memcpy(&shir[0],      &(pData->shir[0][0]), lSig * sizeof(float));
            memcpy(&shir[1*lSig], &(pData->shir[3][0]), lSig * sizeof(float));
            memcpy(&shir[2*lSig], &(pData->shir[1][0]), lSig * sizeof(float));
            memcpy(&shir[3*lSig], &(pData->shir[2][0]), lSig * sizeof(float));
//            memcpy(&shir[0], &(pData->shir[0]), lSig * sizeof(float));
//            memcpy(&shir[1*lSig], &(pData->shir[3*lSig]), lSig * sizeof(float));
//            memcpy(&shir[2*lSig], &(pData->shir[1*lSig]), lSig * sizeof(float));
//            memcpy(&shir[3*lSig], &(pData->shir[2*lSig]), lSig * sizeof(float));
            break;
    }
    
    /* account for input normalisation scheme */
    for(n=0; n<HOSIRR_MAX_SH_ORDER+2; n++){  o[n] = n*n;  }
    switch(norm){
        case N3D_NORM:  /* already in N3D, do nothing */
            break;
        case FUMA_NORM: /* (same as converting SN3D->N3D for first-order) */
        case SN3D_NORM: /* convert to N3D */
            for (n = 0; n<order+1; n++)
                for (ch = o[n]; ch<o[n+1]; ch++)
                    for(i = 0; i<lSig; i++)
                        shir[ch*lSig+i] *= sqrtf(2.0f*(float)n+1.0f);
            break;
    }
    
    /* normalise such that peak of the omni is 1 */
    utility_simaxv(shir, lSig, &maxInd); /* index of max(abs(omni)) */
    peakNorm = 1.0f/fabsf(shir[maxInd]);
    utility_svsmul(shir, &peakNorm, nSH*lSig, shir);
    
    /* isolate first peak */
    if(BB1stPeak){
        utility_simaxv(shir, lSig, &maxInd); /* index of max(abs(omni)) */
        
        /* calculate window and extract peak */
        dirwinsize = 64;
        if(maxInd <dirwinsize/2)
            BB1stPeak=0;
        else{
            shir_tmp = malloc1d(nSH*lSig*sizeof(float));
            shir_tmp2 = malloc1d(nSH*lSig*sizeof(float));
            memcpy(shir_tmp, shir, nSH*lSig*sizeof(float));
            direct_win = calloc1d(nSH*lSig,sizeof(float));
            for(i=0; i<nSH; i++)
                getWindowingFunction(WINDOWING_FUNCTION_HANN, dirwinsize+1 /* force symmetry */, &direct_win[i*lSig + maxInd - dirwinsize/2]);
            shir_direct = malloc1d(nSH*lSig*sizeof(float));
            utility_svvmul(shir_tmp, direct_win, nSH*lSig, shir_direct);
            
            /* flip window and use it to remove peak from the input */
            for(i=0; i<nSH*lSig; i++)
                direct_win[i] = 1.0f-direct_win[i];
            utility_svvmul(shir_tmp, direct_win, nSH*lSig, shir_tmp2);
            memcpy(shir, shir_tmp2, nSH*lSig*sizeof(float));
            
            free(shir_tmp);
            free(shir_tmp2);
            free(direct_win);
        }
    }
    
    /* zero pad the signal's start and end for STFT */
    lSig_pad = winsize/2 + winsize*2 + lSig; // winsize/2 at head, sinsize*2 at tail
    shir_pad = calloc1d(nSH*lSig_pad, sizeof(float));
    for(i=0; i<nSH; i++)
        memcpy(&shir_pad[i*lSig_pad + winsize/2], &(shir[i*lSig]), lSig*sizeof(float));
 
    /* VBAP gain table */
    strcpy(pData->progressText,"Computing VBAP Gain Table");
    pData->progress0_1 = 0.0f;
    gtable = NULL;
    aziRes = 1;
    elevRes = 1;
    generateVBAPgainTable3D((float*)pData->loudpkrs_dirs_deg, nLS, aziRes, elevRes, 0, 0, 0.0f, &gtable, &N_gtable, &N_tri);
    
    /* Sector design */
    strcpy(pData->progressText,"Computing Sector Coefficients");
    numSec = order == 1 ? 1 : __Tdesign_nPoints_per_degree[2*order-1];
    sec_dirs_deg = (float*)__HANDLES_Tdesign_dirs_deg[2*order-1];
    order_sec = order-1;
    nSH_sec = (order_sec+1)*(order_sec+1);
    A_xyz = malloc1d(nSH*nSH_sec*3*sizeof(float_complex));
    computeVelCoeffsMtx(order_sec, A_xyz);
    sectorCoeffs_tmp = malloc1d((numSec*4)*nSH*sizeof(float));
    sectorCoeffs = malloc1d((numSec*4)*nSH*sizeof(float_complex));
    sectorCoeffs_syn = malloc1d((numSec*4)*nSH*sizeof(float_complex));
    normSec = computeSectorCoeffsEP(order_sec, A_xyz, SECTOR_PATTERN_PWD, sec_dirs_deg, numSec, sectorCoeffs_tmp);
    for(i=0; i<numSec*4*nSH; i++){
        sectorCoeffs[i] = cmplxf(sectorCoeffs_tmp[i], 0.0f); /* real->complex data type */
        sectorCoeffs_syn[i] = cmplxf(sectorCoeffs_tmp[i]/sqrtf(4.0f*M_PI), 0.0f);
    }
    free(sectorCoeffs_tmp);
    free(A_xyz);
    
    /* Inits */
    fftsize = winsize*2; 
    hopsize = winsize/2;       /* half the window size time-resolution */
    nBins_anl = winsize/2 + 1; /* nBins used for analysis */
    nBins_syn = fftsize/2 + 1; /* nBins used for synthesis */
    nearestVal = 10e5f;
    for(i=0; i<nBins_anl; i++){
        tmp = fabsf((float)i*(fs/(float)winsize) - MAX_DIFF_FREQ_HZ);
        if(tmp < nearestVal){
            nearestVal = tmp;
            maxDiffFreq_idx = i;
        }
    }
    
    /* transform window (symmetric Hann - 'hanning' in matlab) */
    win = malloc1d(winsize*sizeof(float));
    for(i=0; i<winsize; i++)
        win[i] = powf(sinf((float)i*(M_PI/(float)winsize)), 2.0f);
    
    /* diffuse stream rendering intialisations */
    a2eNorm = (float)nLS/(sqrtf((float)nLS));
    if(order==1){
        D_ls_tmp = malloc1d(nLS*nSH*sizeof(float));
        getLoudspeakerDecoderMtx((float*)pData->loudpkrs_dirs_deg, nLS, LOUDSPEAKER_DECODER_SAD, order, 0, D_ls_tmp);
        utility_svsmul(D_ls_tmp, &a2eNorm, nLS*nSH, D_ls_tmp);
        D_ls = malloc1d(nLS*nSH*sizeof(float_complex));
        for(i=0; i<nLS*nSH; i++)
            D_ls[i] = cmplxf(D_ls_tmp[i], 0.0f);
    }
    else{
        Y_enc_tmp = malloc1d(nSH_sec*numSec*sizeof(float));
        getRSH(order_sec, sec_dirs_deg, numSec, Y_enc_tmp);
        D_ls_tmp = malloc1d(nLS*nSH_sec*sizeof(float));
        getLoudspeakerDecoderMtx((float*)pData->loudpkrs_dirs_deg, nLS, LOUDSPEAKER_DECODER_SAD, order_sec, 0, D_ls_tmp);
        utility_svsmul(D_ls_tmp, &a2eNorm, nLS*nSH_sec, D_ls_tmp);
        Y_enc = malloc1d(nSH_sec*numSec*sizeof(float_complex));
        D_ls = malloc1d(nLS*nSH_sec*sizeof(float_complex));
        for(i=0; i<nSH_sec*numSec; i++)
            Y_enc[i] = cmplxf(Y_enc_tmp[i], 0.0f);
        for(i=0; i<nLS*nSH_sec; i++)
            D_ls[i] = cmplxf(D_ls_tmp[i], 0.0f);
        free(Y_enc_tmp);
    }
    free(D_ls_tmp);
    
    /* mem alloc */
    saf_rfft_create(&hFFT_syn, fftsize);
    saf_rfft_create(&hFFT_anl, fftsize/2);
    lsir_ndiff = calloc1d(nLS * (lSig + (2*fftsize)), sizeof(float));
    lsir_diff  = calloc1d(nLS * (lSig + (2*fftsize)), sizeof(float));
    lsir_win = malloc1d(fftsize * sizeof(float));
    prev_energy_smoothed = calloc1d(numSec, sizeof(float));
    prev_ixyz_smooth = calloc1d(numSec*3, sizeof(float));
    azim = malloc1d(numSec*nBins_anl*sizeof(float));
    elev = malloc1d(numSec*nBins_anl*sizeof(float));
    diffs = malloc1d(numSec*nBins_anl*sizeof(float));
    ndiffs_sqrt = malloc1d(nBins_anl*sizeof(float));
    diffs_sqrt  = malloc1d(nBins_anl*sizeof(float_complex));
    insig_win = calloc1d(fftsize,sizeof(float));
    inspec_syn = calloc1d(nSH*nBins_syn,sizeof(float_complex));      //////ma->ca
    inspec_anl = calloc1d(nSH*nBins_anl,sizeof(float_complex));//////ma->ca
    s_anl = malloc1d(4*numSec*nBins_anl*sizeof(float_complex));
    WXYZ_sec = malloc1d(4*nBins_anl*sizeof(float_complex));
    z_diff = malloc1d(numSec*nBins_syn*sizeof(float_complex));
    z_00 = malloc1d(nBins_syn*sizeof(float_complex));
    a_diff = malloc1d(nSH*nBins_syn*sizeof(float_complex));
    M_ifft = malloc1d(winsize*sizeof(float));
    M_ifft_fl = calloc1d(fftsize, sizeof(float));
    ndiffgains = malloc1d(nLS*nBins_anl*sizeof(float_complex));
    ndiffgains_interp = malloc1d(nLS*nBins_syn*sizeof(float_complex));
    diffgains_interp = malloc1d(nBins_syn*sizeof(float_complex));
    outspec_ndiff = malloc1d(nLS*nBins_syn*sizeof(float_complex));
    outspec_diff = malloc1d(nLS*nBins_syn*sizeof(float_complex));
    
    /* Main processing loop */
    strcpy(pData->progressText,"HOSIRR - Rendering");
    idx = frameCount = 0;
    while (idx + winsize < lSig + 2*winsize){

        /* update progress */
        pData->progress0_1 = (float)idx/(float)(lSig + 2*winsize);
 
        /* Window input and transform to frequency domain */
        for(i=0; i<nSH; i++){
            for(j=0; j<winsize; j++)
                insig_win[j] = win[j] * shir_pad[i*lSig_pad+idx+j];
            saf_rfft_forward(hFFT_syn, insig_win, &inspec_syn[i*nBins_syn]);
            for(j=0, k=0; j<nBins_anl; j++, k+=fftsize/winsize)
                inspec_anl[i*nBins_anl+j] = inspec_syn[i*nBins_syn+k];
        }
        
        /* Form weighted pressure-velocity signals */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, numSec*4, nBins_anl, nSH, &calpha,
                    sectorCoeffs, nSH,
                    inspec_anl, nBins_anl, &cbeta,
                    s_anl, nBins_anl);
        
        /* SIRR analysis for each sector */
        for(n=0; n<numSec; n++){
            for(i=0; i<4; i++)
                memcpy(&WXYZ_sec[i*nBins_anl], &s_anl[n*4*nBins_anl + i*nBins_anl], nBins_anl*sizeof(float_complex));
            
            /* compute Intensity vector for each frequency bin to estimate DoA */
            for(j=0; j<nBins_anl; j++){
                for(i=0; i<3; i++)
                    intensity[i] = crealf( ccmulf(conjf(WXYZ_sec[j]), WXYZ_sec[(i+1)*nBins_anl+j]));
                printf("intensity %.7f, %.7f, %.7f\n", intensity[0], intensity[1], intensity[2]);
                azim[n*nBins_anl+j] = atan2f(intensity[1], intensity[0])*180.0f/M_PI;
                elev[n*nBins_anl+j] = atan2f(intensity[2], sqrtf(powf(intensity[0], 2.0f) + powf(intensity[1], 2.0f)))*180.0f/M_PI;
            }
                
            /* Compute broad-band active-intensity vector */
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, 4, 4, maxDiffFreq_idx+1, &calpha,
                        WXYZ_sec, nBins_anl,
                        WXYZ_sec, nBins_anl, &cbeta,
                        FLATTEN2D(pvCOV), 4);
            for(i=0; i<3; i++)
                intensity[i] = crealf(pvCOV[1+i][0]);
            energy = 0.0f;
            for(i=0; i<4; i++)
                energy += crealf(pvCOV[i][i])*0.5f;
            
//                    printf("intensity %.5f, %.5f, %.5f\n", intensity[0], intensity[1], intensity[2]);
            
            /* Estimating and time averaging of boadband diffuseness */
            normSecIntensity_smoothed = 0.0f;
            for(i=0; i<3; i++){
                ixyz_smoothed[i] = (1.0f-ALPHA_DIFF_COEFF)*intensity[i] + ALPHA_DIFF_COEFF*prev_ixyz_smooth[n*3+i];
                prev_ixyz_smooth[n*3+i] = ixyz_smoothed[i];
                normSecIntensity_smoothed += powf(fabsf(ixyz_smoothed[i]), 2.0f);
            }
            energy_smoothed = (1.0f-ALPHA_DIFF_COEFF)*energy + ALPHA_DIFF_COEFF*prev_energy_smoothed[n];
            prev_energy_smoothed[n] = energy_smoothed;
            normSecIntensity_smoothed = sqrtf(normSecIntensity_smoothed);
            for(i=0; i<nBins_anl; i++)
                diffs[n*nBins_anl+i] = 1.0f - (normSecIntensity_smoothed / (energy_smoothed+2.23e-10f));
        }
        
        /* SIRR Synthesis for each sector */
        memset(outspec_ndiff, 0, nLS*nBins_syn*sizeof(float_complex));
        memset(outspec_diff,  0, nLS*nBins_syn*sizeof(float_complex));
        for(n=0; n<numSec; n++){
            for(i=0; i<nBins_anl; i++){
                ndiffs_sqrt[i] = sqrtf(1.0f-diffs[n*nBins_anl+i]);
                diffs_sqrt[i]  = cmplxf(sqrtf(diffs[n*nBins_anl+i]), 0.0f);
            }
            
            /* Gain factor computation */
            for(i=0; i<nBins_anl; i++){
                N_azi = (int)(360.0f / (float)aziRes + 0.5f) + 1;
                aziIndex = (int)(matlab_fmodf(azim[n*nBins_anl+i] + 180.0f, 360.0f) / (float)aziRes + 0.5f);
                elevIndex = (int)((elev[n*nBins_anl+i] + 90.0f) / (float)elevRes + 0.5f);
                idx3d = elevIndex * N_azi + aziIndex;
                for(j=0; j<nLS; j++)
                    ndiffgains[j*nBins_anl+i] = cmplxf(gtable[idx3d*nLS+j] * ndiffs_sqrt[i], 0.0f);
            }
            
            /* Interpolate panning filters  */
            for(i=0; i<nLS; i++){
                saf_rfft_backward(hFFT_anl, &ndiffgains[i*nBins_anl], M_ifft);
                /* flip */
                for(j=0; j<winsize/2; j++){
                    M_ifft_fl[j] = M_ifft[winsize/2+j];
                    M_ifft_fl[winsize/2+j] = M_ifft[j];
                }
                saf_rfft_forward(hFFT_syn, M_ifft_fl, &ndiffgains_interp[i*nBins_syn]);
            }
            
            /* Generate non-diffuse stream */
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1, nBins_syn, nSH, &calpha,
                        &sectorCoeffs_syn[n*4*nSH], nSH,
                        inspec_syn, nBins_syn, &cbeta,
                        z_00, nBins_syn);
            for(i=0; i<nLS; i++)
                for(j=0; j<nBins_syn; j++)
                    outspec_ndiff[i*nBins_syn+j] = ccaddf(outspec_ndiff[i*nBins_syn+j],
                                                          ccmulf(ndiffgains_interp[i*nBins_syn+j], crmulf(z_00[j], sqrtf(normSec))) );
            
            /* Interpolate diffs  */
            saf_rfft_backward(hFFT_anl, diffs_sqrt, M_ifft);
            /* flip */
            for(j=0; j<winsize/2; j++){
                M_ifft_fl[j] = M_ifft[winsize/2+j];
                M_ifft_fl[winsize/2+j] = M_ifft[j];
            }
            saf_rfft_forward(hFFT_syn, M_ifft_fl, diffgains_interp);
            
            /* Generate diffuse stream */
            if(order==1){
                for(i=0; i<nSH; i++)
                    for(j=0; j<nBins_syn; j++)
                        a_diff[i*nBins_syn+j] = ccmulf(diffgains_interp[j], crmulf(inspec_syn[i*nBins_syn+j], 1.0f/sqrtf((float)nSH)));
            }
            else{
                for(j=0; j<nBins_syn; j++)
                    z_diff[n*nBins_syn+j] = crmulf(ccmulf(diffgains_interp[j], z_00[j]), 1.0f/sqrtf(numSec));
            }
        }
        
        /* Decode diffuse stream to loudspeakers */
        if(order==1){
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nLS, nBins_syn, nSH, &calpha,
                        D_ls, nSH,
                        a_diff, nBins_syn, &cbeta,
                        outspec_diff, nBins_syn);
        }
        else{
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH_sec, nBins_syn, numSec, &calpha,
                        Y_enc, numSec,
                        z_diff, nBins_syn, &cbeta,
                        a_diff, nBins_syn);
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nLS, nBins_syn, nSH_sec, &calpha,
                        D_ls, nSH_sec,
                        a_diff, nBins_syn, &cbeta,
                        outspec_diff, nBins_syn);
        }

        /* Overlap-add synthesis */
        for(i=0; i<nLS; i++){
            saf_rfft_backward(hFFT_syn, &outspec_ndiff[i*nBins_syn], lsir_win);
            for(jl=idx, k=0; jl<fftsize+idx; jl++, k++)
                lsir_ndiff[i*(lSig + (2*fftsize))+jl] += lsir_win[k];
            saf_rfft_backward(hFFT_syn, &outspec_diff[i*nBins_syn], lsir_win);
            for(jl=idx, k=0; jl<fftsize+idx; jl++, k++)
                lsir_diff[i*(lSig + (2*fftsize))+jl] += lsir_win[k];
        }
        
        idx += hopsize;
        frameCount++;
    }
    
    /* Remove delay caused by processing */
    delay = winsize; ///2;    // from winsize/2
    for(i=0; i<nLS; i++){
        memcpy(&lsir_ndiff[i*lSig], &lsir_ndiff[i*(lSig + (2*fftsize))+delay], lSig*sizeof(float));
        memcpy(&lsir_diff[i*lSig], &lsir_diff[i*(lSig + (2*fftsize))+delay], lSig*sizeof(float));
    }
    lsir_ndiff = realloc1d(lsir_ndiff, nLS*lSig*sizeof(float));
    lsir_diff = realloc1d(lsir_diff, nLS*lSig*sizeof(float));
    
    /* Convolution with exponential decaying noise to decorrelate the diffuse stream */
    strcpy(pData->progressText,"Decorrelating Diffuse Stream");
    rir_filt = NULL;
    synthesiseNoiseReverb(nLS, fs, t60, fc, 6, 1, &rir_filt, &rir_len);
    fftfilt(lsir_diff, rir_filt, lSig, rir_len, nLS, lsir_diff);
    
    /* Re-introduce first peak based on broadband analysis  */
    if(BB1stPeak){
        /* Broad-band intensity */
        for(i=0; i<3; i++)
            for(j=0; j<lSig; j++)
                IntensityBB[i] += shir_direct[0*lSig + j] * shir_direct[(i+1)*lSig + j]/sqrtf(3.0f); /* "sqrtf(3.0f)" for N3D->SN3D */
        IntensityBB_XYZ[0] = IntensityBB[2];
        IntensityBB_XYZ[1] = IntensityBB[0];
        IntensityBB_XYZ[2] = IntensityBB[1];
        
        /* Get gtable index for DoA */
        N_azi = (int)(360.0f / aziRes + 0.5f) + 1;
        aziIndex = (int)(matlab_fmodf( (atan2f(IntensityBB_XYZ[1], IntensityBB_XYZ[0]) *180.0f/M_PI) + 180.0f, 360.0f) / (float)aziRes + 0.5f);
        elevIndex = (int)(((atan2f(IntensityBB_XYZ[2], sqrtf(powf(IntensityBB_XYZ[0], 2.0f) + powf(IntensityBB_XYZ[1], 2.0f)))*180.0f/M_PI) + 90.0f) / (float)elevRes + 0.5f);
        idx3d = elevIndex * N_azi + aziIndex;
        
        /* pan omni and sum to non-diffuse stream */
        for(i=0; i<nLS; i++)
            for(j=0; j<lSig; j++)
                lsir_ndiff[i*lSig+j] += gtable[idx3d*nLS+i] * shir_direct[j];
        
        free(shir_direct);
    }
    
    /* Sum the two streams, store */
    pData->lsir = realloc1d(pData->lsir, nLS*lSig*sizeof(float));
    for(i=0; i<nLS; i++){
        for(j=0; j<lSig; j++){
            if(wetDry<1.0f)
                pData->lsir[i*lSig+j] = wetDry * lsir_ndiff[i*lSig+j] + lsir_diff[i*lSig+j];
            else
                pData->lsir[i*lSig+j] = lsir_ndiff[i*lSig+j] + (2.0f-wetDry) * lsir_diff[i*lSig+j];
        }
    }
    
    /* indicate that rendering is complete */
    pData->progress0_1 = 1.0f;
    pData->lsRIR_status = LS_RIR_STATUS_RENDERED;
 
    /* Clean-up */
    free(shir);
    free(shir_pad);
    free(gtable);
    free(sectorCoeffs);
    free(sectorCoeffs_syn);
    free(win);
    if(order>1)
        free(Y_enc);
    free(D_ls);
    saf_rfft_destroy(&hFFT_syn);
    saf_rfft_destroy(&hFFT_anl);
    free(lsir_ndiff);
    free(lsir_diff);
    free(lsir_win);
    free(prev_energy_smoothed);
    free(prev_ixyz_smooth);
    free(azim);
    free(elev);
    free(diffs);
    free(ndiffs_sqrt);
    free(diffs_sqrt);
    free(insig_win);
    free(inspec_syn);
    free(inspec_anl);
    free(s_anl);
    free(WXYZ_sec);
    free(z_diff);
    free(z_00);
    free(a_diff);
    free(M_ifft);
    free(M_ifft_fl);
    free(ndiffgains);
    free(ndiffgains_interp);
    free(diffgains_interp);
    free(outspec_ndiff);
    free(outspec_diff);
    free(rir_filt);
}


/* Set Functions */

int hosirrlib_setAmbiRIR
(
    void* const hHS,
    const float** H,
    int numChannels,
    int numSamples,
    int sampleRate
)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    /* Check channel count to see if input is actually in the SHD */
    if(fabsf(sqrtf((float)numChannels)-floorf(sqrtf((float)numChannels)))>0.0001f){
        pData->ambiRIR_status = AMBI_RIR_STATUS_NOT_LOADED;
        pData->nSH = -1;
        pData->ambiRIRorder = -1;
        pData->ambiRIRlength_seconds = pData->ambiRIRlength_samples = 0.0f;
        pData->ambiRIRsampleRate = 0;
        return (int)(pData->ambiRIR_status);
    }
    
    /* if it is, store RIR data */
    pData->ambiRIRorder = HOSIRR_MIN(sqrt(numChannels-1), HOSIRR_MAX_SH_ORDER);
    pData->nSH = numChannels;
    pData->analysisOrder = pData->ambiRIRorder;
    pData->ambiRIRlength_samples = numSamples;
    pData->ambiRIRsampleRate = sampleRate;
    pData->ambiRIRlength_seconds = (float)pData->ambiRIRlength_samples/(float)pData->ambiRIRsampleRate;
    
    // (re)allocate memory and copy in SH RIR data
    pData->shir = (float**)realloc2d((void**)pData->shir, numChannels, numSamples, sizeof(float));
    for(int i = 0; i < numChannels; i++)
        utility_svvcopy(H[i], numSamples, pData->shir[i]); // mtm
//        memcpy(&(pData->shir[i*numSamples]), H[i], numSamples * sizeof(float));

    /* set FLAGS */
    pData->ambiRIR_status = AMBI_RIR_STATUS_LOADED;
    pData->lsRIR_status = LS_RIR_STATUS_NOT_RENDERED;
    
    return (int)(pData->ambiRIR_status);
}

void hosirrlib_setBroadBandFirstPeakFLAG(void* const hHS, int newState)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    pData->broadBandFirstPeakFLAG = newState;
    pData->lsRIR_status = LS_RIR_STATUS_NOT_RENDERED;
}

void hosirrlib_setWindowLength(void* const hHS, int newValue)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    /* round to nearest multiple of 2 */
    pData->windowLength = newValue % 2 == 0 ? newValue : newValue + 2 - (newValue % 2);
    /* clamp within bounds */
    pData->windowLength = HOSIRR_CLAMP(pData->windowLength, MIN_WINDOW_LENGTH, MAX_WINDOW_LENGTH);
    pData->lsRIR_status = LS_RIR_STATUS_NOT_RENDERED;
}

void hosirrlib_setWetDryBalance(void* const hHS, float newValue)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    pData->wetDryBalance = newValue;
    pData->lsRIR_status = LS_RIR_STATUS_NOT_RENDERED;
}


void hosirrlib_setAnalysisOrder(void* const hHS, int newValue)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    pData->analysisOrder = HOSIRR_MIN(HOSIRR_MAX(newValue,1), HOSIRR_MAX_SH_ORDER);
    /* FUMA only supports 1st order */
    if(pData->analysisOrder!=ANALYSIS_ORDER_FIRST && pData->chOrdering == FUMA_ORDER)
        pData->chOrdering = ACN_ORDER;
    if(pData->analysisOrder!=ANALYSIS_ORDER_FIRST && pData->norm == FUMA_NORM)
        pData->norm = SN3D_NORM;
    pData->lsRIR_status = LS_RIR_STATUS_NOT_RENDERED;
}

void hosirrlib_setLoudspeakerAzi_deg(void* const hHS, int index, float newAzi_deg)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    if(newAzi_deg>180.0f)
        newAzi_deg = -360.0f + newAzi_deg;
    newAzi_deg = HOSIRR_MAX(newAzi_deg, -180.0f);
    newAzi_deg = HOSIRR_MIN(newAzi_deg, 180.0f);
    pData->loudpkrs_dirs_deg[index][0] = newAzi_deg;
    pData->lsRIR_status = LS_RIR_STATUS_NOT_RENDERED;
}

void hosirrlib_setLoudspeakerElev_deg(void* const hHS, int index, float newElev_deg)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    newElev_deg = HOSIRR_MAX(newElev_deg, -90.0f);
    newElev_deg = HOSIRR_MIN(newElev_deg, 90.0f);
    pData->loudpkrs_dirs_deg[index][1] = newElev_deg;
    pData->lsRIR_status = LS_RIR_STATUS_NOT_RENDERED;
}

void hosirrlib_setNumLoudspeakers(void* const hHS, int new_nLoudspeakers)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    pData->nLoudpkrs = new_nLoudspeakers > MAX_NUM_LOUDSPEAKERS ? MAX_NUM_LOUDSPEAKERS : new_nLoudspeakers;
    pData->nLoudpkrs = HOSIRR_MAX(MIN_NUM_LOUDSPEAKERS, pData->nLoudpkrs);
    pData->lsRIR_status = LS_RIR_STATUS_NOT_RENDERED;
}

void  hosirrlib_setOutputConfigPreset(void* const hHS, int newPresetID)
{
    hosirrlib_data *pData = ( hosirrlib_data*)(hHS);
//    loadLoudspeakerArrayPreset(newPresetID, pData->loudpkrs_dirs_deg, &(pData->nLoudpkrs));
    loadSphDesignPreset(newPresetID, pData->loudpkrs_dirs_deg, &(pData->nLoudpkrs)); // TODO: SPHDESIGNs coincide with LS arrays for now
    pData->lsRIR_status = LS_RIR_STATUS_NOT_RENDERED;
}

void hosirrlib_setChOrder(void* const hHS, int newOrder)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    if((CH_ORDERING)newOrder != FUMA_ORDER || pData->analysisOrder==ANALYSIS_ORDER_FIRST) /* FUMA only supports 1st order */
        pData->chOrdering = (CH_ORDERING)newOrder;
    pData->lsRIR_status = LS_RIR_STATUS_NOT_RENDERED;
}

void hosirrlib_setNormType(void* const hHS, int newType)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    if((NORMALIZATION_TYPES)newType != FUMA_NORM || pData->analysisOrder==ANALYSIS_ORDER_FIRST) /* FUMA only supports 1st order */
        pData->norm = (NORMALIZATION_TYPES)newType;
    pData->lsRIR_status = LS_RIR_STATUS_NOT_RENDERED;
}


/* Get Functions */

float hosirrlib_getProgress0_1(void* const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    return pData->progress0_1;
}

void hosirrlib_getProgressText(void* const hHS, char* text)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    memcpy(text, pData->progressText, HOSIRR_PROGRESSTEXT_CHAR_LENGTH*sizeof(char));
}

int hosirrlib_getAmbiRIRstatus(void* const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    return pData->ambiRIR_status;
}

int hosirrlib_getLsRIRstatus(void* const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    return pData->lsRIR_status;
}

int hosirrlib_getAmbiRIRinputOrder(void* const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    return pData->ambiRIRorder;
}

int hosirrlib_getAmbiRIRlength_samples(void* const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    return pData->ambiRIRlength_samples;
}

float hosirrlib_getAmbiRIRlength_seconds(void* const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    return pData->ambiRIRlength_seconds;
}

int hosirrlib_getAmbiRIRsampleRate(void* const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    return pData->ambiRIRsampleRate;
}

int hosirrlib_getBroadBandFirstPeakFLAG(void* const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    return pData->broadBandFirstPeakFLAG;
}

void hosirrlib_getLsRIR(void* const hHS, float** lsRIR)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    int i;
    if(pData->lsRIR_status == LS_RIR_STATUS_RENDERED)
        for(i=0; i< pData->nLoudpkrs; i++)
            memcpy(lsRIR[i], &(pData->lsir[i*(pData->ambiRIRlength_samples)]), pData->ambiRIRlength_samples*sizeof(float));
}

int hosirrlib_getWindowLength (void* const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    return pData->windowLength;
}

float hosirrlib_getWetDryBalance(void* const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    return pData->wetDryBalance;
}

int hosirrlib_getAnalysisOrder(void* const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    return pData->analysisOrder;
}

float hosirrlib_getLoudspeakerAzi_deg(void* const hHS, int index)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    return pData->loudpkrs_dirs_deg[index][0];
}

float hosirrlib_getLoudspeakerElev_deg(void* const hHS, int index)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    return pData->loudpkrs_dirs_deg[index][1];
}

int hosirrlib_getNumLoudspeakers(void* const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    return pData->nLoudpkrs;
}

int hosirrlib_getMaxNumLoudspeakers()
{
    return MAX_NUM_LOUDSPEAKERS;
}

int hosirrlib_getChOrder(void* const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    return (int)pData->chOrdering;
}

int hosirrlib_getNormType(void* const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    return (int)pData->norm;
}

int hosirrlib_getNumDirections(void* const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    return pData->nDir;
}
