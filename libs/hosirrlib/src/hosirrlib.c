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

void hosirrlib_create
(
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
    pData->t60Buf        = NULL;    // nDir x nBand
    
    // depend on output design (nDir) AND input RIR (nSamp)
    pData->rirBuf        = NULL;    // nSH x nSamp
    pData->rirBuf_beams  = NULL;    // nDir x nSamp
    pData->fdnBuf        = NULL;    // nDir x nSamp
    pData->fdnBuf_bands  = NULL;    // nBand x nDir x nSamp
    pData->edcBuf_rir    = NULL;    // nDir x nBand x nSamp
    pData->edcBuf_fdn    = NULL;    // nDir x nBand x nSamp
    pData->fdnBuf_shd    = NULL;    // nSH x nSamp
    
    pData->H_bandFilt    = NULL;    // nBand x filtOrder+1
    pData->bandXOverFreqs = NULL;
    
    /* zero out the state of buffer resources */
    hosirrlib_setUninitialized(pData);
    
    /* Initialize octave band filters */
    pData->nBand = 8;
    pData->bandFiltOrder = 400;
    hosirrlib_initBandFilters(pData);
    
    /* Initialize spherical desing for directional analysis */
    pData->beamType = STATIC_BEAM_TYPE_HYPERCARDIOID;
    // default design currently just set by enum LOUDSPEAKER_ARRAY_PRESET_15PX
    loadLoudspeakerArrayPreset(LOUDSPEAKER_ARRAY_PRESET_15PX, pData->loudpkrs_dirs_deg, &(pData->nLoudpkrs));
    pData->nDir = pData->nLoudpkrs;
    
    
    /* original hosirrlib */
    
    pData->progress0_1 = 0.0f;
    pData->progressText = malloc1d(HOSIRR_PROGRESSTEXT_CHAR_LENGTH*sizeof(char));
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
    //loadLoudspeakerArrayPreset(LOUDSPEAKER_ARRAY_PRESET_T_DESIGN_24, pData->loudpkrs_dirs_deg, &(pData->nLoudpkrs));
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_SN3D;
    pData->broadBandFirstPeakFLAG = 1;
    pData->windowLength = DEFAULT_WINDOW_LENGTH;
    pData->wetDryBalance = 1.0f;
}

void hosirrlib_destroy(void ** const phHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(*phHS);
    
    if (pData != NULL) {
        
        /* new hodecaylib */
        
        // depend only on output design (nDir)
        free(pData->encBeamCoeffs);
        free(pData->decBeamCoeffs);
        free(pData->dirGainBuf);
        free(pData->t60Buf);
        // depend on output design (nDir) AND input RIR (nSamp)
        free(pData->rirBuf);
        free(pData->rirBuf_bands);
        free(pData->rirBuf_beams);
        free(pData->edcBuf_rir);
        free(pData->fdnBuf);
        free(pData->fdnBuf_bands);
        free(pData->edcBuf_fdn);
        free(pData->fdnBuf_shd);
        // unchanging after init
        free(pData->H_bandFilt);
        free(pData->bandXOverFreqs);
        
        /* original hosirrlib */
        
        free(pData->shir);
        free(pData->lsir);
        free(pData->progressText);
        free(pData);
        pData = NULL;
    }
}

void hosirrlib_initBandFilters(void* const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("initBandFilters called.\n");
    // if this is the first time this function is called...
    if (pData->H_bandFilt == NULL) {
        
        pData->bandXOverFreqs = malloc1d((pData->nBand-1) * sizeof(float));
        
        // 'bandFreqs' are basically center freqs used to determine crossover freqs.
        // The lowest and highest bands are LP and HP filters, so 125 Hz isn't really a "center freq" of the LPF. And the HPF "center freq" would implicitly be 16kHz.
        // Must be size pData->nBand-1
        float bandFreqs[7] = {125.f, 250.f, 500.f, 1000.f, 2000.f, 4000.f, 8000.f};
        
        // The xover freqs are the half-octave step above of the center freqs.
        for (int i=0; i < pData->nBand-1; i++)
            pData->bandXOverFreqs[i] = bandFreqs[i] * sqrtf(2.0f);
        
        // Allocate filter coefficients (nBand x bandFiltOrder + 1)
        pData->H_bandFilt = (float**)realloc2d((void**)pData->H_bandFilt,
                                               pData->nBand,
                                               pData->bandFiltOrder + 1,
                                               sizeof(float));
        // Compute FIR Filterbank coefficients
        FIRFilterbank(pData->bandFiltOrder, pData->bandXOverFreqs, pData->nBand-1,
                      pData->fs, WINDOWING_FUNCTION_HAMMING, 1,
                      FLATTEN2D(pData->H_bandFilt));
    }
}

int hosirrlib_setRIR
(
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
    pData->shOrder  = SAF_MIN(sqrt(numChannels-1), MAX_SH_ORDER);
    pData->fs       = sampleRate;
    pData->duration = numSamples / sampleRate;
    
    /* (Re)Alloc and copy in input RIR */
    pData->rirBuf = (float**)realloc2d((void**)pData->rirBuf, numChannels, numSamples, sizeof(float));
    for(int i = 0; i < numChannels; i++)
        utility_svvcopy(H[i], numSamples, pData->rirBuf[i]);
    
    /* set FLAGS */
    pData->analysisStage = RIR_LOADED;
    
    /* Alloc processing resources */
    hosirrlib_allocProcBufs(pData);
    
    /* old vars */
    pData->ambiRIRorder = SAF_MIN(sqrt(numChannels-1), MAX_SH_ORDER);
    pData->analysisOrder = pData->ambiRIRorder;
    pData->ambiRIRlength_samples = numSamples;
    pData->ambiRIRsampleRate = sampleRate;
    pData->ambiRIRlength_seconds = (float)pData->ambiRIRlength_samples / (float)pData->ambiRIRsampleRate;
    
    // (re)allocate memory and copy in SH RIR data
    pData->shir = (float**)realloc2d((void**)pData->shir, numChannels, numSamples, sizeof(float));
    for(int i = 0; i < numChannels; i++)
        utility_svvcopy(H[i], numSamples, pData->shir[i]);
//        memcpy(&(pData->shir[i*numSamples]), H[i], numSamples * sizeof(float));

    /* set FLAGS */
    pData->ambiRIR_status = AMBI_RIR_STATUS_LOADED;
    pData->lsRIR_status = LS_RIR_STATUS_NOT_RENDERED;
    
    return (int)(pData->ambiRIR_status); // TODO: consider use of returned value
}

void hosirrlib_setUninitialized(void* const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("setUninitialized called.\n");
    // pending initializations
    pData->nSH = -1; // input vars (assume for now in params = out params)
    pData->nSamp = -1;
    pData->shOrder = -1;
    pData->fs = -1;
    pData->directOnsetIdx = -1;
    pData->diffuseOnsetIdx = -1;
    pData->duration = 0.0f;
    pData->analysisStage = RIR_NOT_LOADED;
}

// (re)allocate the buffers used for storing intermediate processing data
void hosirrlib_allocProcBufs(void * const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("allocProcBufs called\n");
    /* Check previous stages are complete */
    if (pData->analysisStage < ANALYSIS_BUFS_LOADED-1) {
        printf("allocProcBufs called before previous stages were completed: %d\n", pData->analysisStage);
        return; // TODO: handle fail case
    }
    
    const int nSH   = pData->nSH;
    const int nDir  = pData->nDir;
    const int nBand = pData->nBand;
    const int nSamp = pData->nSamp;

    // TODO: remember to zero out any bufs that might need it
    // These members depend on output design (nDir) AND input RIR (nSH, nSamp)
    pData->rirBuf_bands = (float***)realloc3d((void***)pData->rirBuf_bands,
                                             nBand, nSH, nSamp, sizeof(float));
    pData->rirBuf_beams = (float***)realloc3d((void***)pData->rirBuf_beams,
                                             nBand, nDir, nSamp, sizeof(float));
    pData->edcBuf_rir   = (float***)realloc3d((void***)pData->edcBuf_rir,
                                              nBand, nDir, nSamp, sizeof(float));
    pData->edcBuf_fdn   = (float***)realloc3d((void***)pData->edcBuf_fdn,
                                              nBand, nDir, nSamp, sizeof(float));
    pData->fdnBuf_bands = (float***)realloc3d((void***)pData->fdnBuf_bands,
                                             nBand, nDir, nSamp, sizeof(float));
    pData->fdnBuf       = (float**)realloc2d((void**)pData->fdnBuf,
                                             nDir, nSamp, sizeof(float));
    pData->fdnBuf_shd   = (float**)realloc2d((void**)pData->fdnBuf_shd,
                                             nSH, nSamp, sizeof(float));
    pData->encBeamCoeffs = (float**)realloc2d((void**)pData->encBeamCoeffs,
                                              nSH, nDir, sizeof(float));
    pData->decBeamCoeffs = (float**)realloc2d((void**)pData->decBeamCoeffs,
                                              nDir, nSH, sizeof(float));
    // (re)allocate buffers that depend only on the number of output channels
    // These could alternatively only be reallocated when the out design
    // changes, but it's not heavy and it's more concise
    pData->dirGainBuf    = (float**)realloc2d((void**)pData->dirGainBuf,
                                              nBand, nDir, sizeof(float));
    pData->t60Buf        = (float**)realloc2d((void**)pData->t60Buf,
                                              nBand, nDir, sizeof(float));
    pData->analysisStage = ANALYSIS_BUFS_LOADED;
}

void hosirrlib_renderTMP
    (void  *  const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("hosirrlib_renderTMP called.\n");
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

void hosirrlib_processRIR
    (void  *  const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("processRIR called.\n");
    
    /* Check previous stages are complete */
    if (pData->analysisStage < ANALYSIS_BUFS_LOADED-1) {
        printf("processRIR called before previous stages were completed: %d\n", pData->analysisStage);
        return; // TODO: handle fail case
    }
    
    // lag after direct onset to start measurement
    const float directLagSec = 0.002f; // 2 ms
    // direct arrival onset threashold (dB below omni peak)
    const float directOnsetThreshDb = -12.f;
    // diffuse onset threashold (factor below diffuseness peak)
    const float diffuseOnsetThreshDb = 0.707f;
    
    hosirrlib_setDirectOnsetIndex(pData, directOnsetThreshDb);
    hosirrlib_setDiffuseOnsetIndex(pData, diffuseOnsetThreshDb);
    hosirrlib_splitBands(pData, pData->rirBuf, pData->rirBuf_bands, 1, RIR_BANDS_SPLIT);
    hosirrlib_beamformRIR(pData);
    hosirrlib_calcEDC(pData, pData->rirBuf_beams, pData->edcBuf_rir, RIR_EDC_DONE);
    hosirrlib_calcT60(pData,
                      -3.f, 18.f, // startDb <=0, spanDb > 0
                      pData->directOnsetIdx + (int)(pData->fs * directLagSec) // measure after this index
                      );

}

void hosirrlib_setDirectOnsetIndex(void* const hHS, const float thresh_dB)
/*
    NOTE:   The naive approach is to consider the index of absolute max value
            in the buffer. However, with synthetic RIRs, like from a shoebox,
            coincident reflections may have greater magnitude than the direct
            arrival. So this returns the index of the first value that crosses
            the (max abs value - threshold).
    TODO:   Could return true first _peak_.
            This would require a second iteration, searching for a local max
            within, say, 2 ms after the onset index.
*/
{ /*
   * thresh_dB: threshold (pressure, dB) below the absolute max value in the
   *            buffer, above which the onset is considered to have occured.
   */
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("setDirectOnsetIndex called.\n");
    
    /* Check previous stages are complete */
    if (pData->analysisStage < DIRECT_ONSET_FOUND-1) {
        printf("setDirectOnsetIndes called before previous stages were completed: %d\n", pData->analysisStage);
        return; // TODO: handle fail case
    }
    
    const int nSamp = pData->nSamp;
    
    float* vabs_tmp = malloc1d(nSamp * sizeof(float));
    
    /* absolute values of the omni channel */
    int maxIdx;
    utility_svabs(&pData->rirBuf[0][0], nSamp, vabs_tmp); // abs(omni)
    utility_simaxv(vabs_tmp, nSamp, &maxIdx); // index of max(abs(omni))
    
    /* index of first index above threshold */
    float maxVal        = vabs_tmp[maxIdx];
    float onsetThresh   = maxVal * powf(10.f, thresh_dB / 20.f);
    pData->directOnsetIdx = hosirrlib_firstIndexGreaterThan(vabs_tmp, 0, nSamp-1,
                                                          onsetThresh);
    
    pData->analysisStage = DIRECT_ONSET_FOUND;
    free(vabs_tmp);
}

void hosirrlib_setDiffuseOnsetIndex(void* const hHS, const float thresh_fac)
{ /*
   * thresh_fac: threshold (energy, dB) below the absolute max value in the
   *            buffer, above which the onset is considered to have occured.
   */
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("setDiffuseOnsetIndex called.\n");
    
    /* Check previous stages are complete */
    if (pData->analysisStage < DIFFUSENESS_ONSET_FOUND-1) {
        printf("setDiffuseOnsetIndex called before previous stages were completed: %d\n", pData->analysisStage);
        return; // TODO: handle fail case
    }
    /* Check if processing should actually go-ahead */
//    if(pData->ambiRIR_status != AMBI_RIR_STATUS_LOADED
//     ||  pData->lsRIR_status == LS_RIR_STATUS_RENDERED ||
//     ||  pData->lsRIR_status == LS_RIR_STATUS_RENDEREDING_ONGOING
//       ) {
//        printf("bailing in setDiffuseOnsetIndex bc of ambiRIR/lsRIR_status.\n");
//        return;
//    } else {
//        pData->lsRIR_status = LS_RIR_STATUS_RENDEREDING_ONGOING;
//    }
    
    /* take a local copy of current configuration to be thread safe */
    const int   nChan   = pData->nDir;
    const int   nBand   = pData->nBand;
    // const int   nSH     = pData->nSH;
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
    const int lSig_pad  = winsize/2 + winsize*2 + lSig; // winsize/2 at head, sinsize*2 at tail
    float_complex pvCOV[4][4];
    
    /* Max freq bin to calculate diffuseness */
    float nearestVal = 10e5f;
    int   maxDiffFreq_idx = 0;
    for(int i = 0; i < nBins_anl; i++){
        // bin_idx * binwidth (Hz)
        float tmp = fabsf((float)i * (fs / (float)winsize) - MAX_DIFF_FREQ_HZ); // TODO: should be fs/fftsize, not fs/winsize?
        if(tmp < nearestVal){
            nearestVal = tmp;
            maxDiffFreq_idx = i;
        }
    }
    
    float* pvir, * pvir_pad, * win, * insig_win, * diff_win;
    float_complex* inspec_anl, * inspec_syn, * wxyzspec_win;
    void* hFFT_syn, *hFFT_anl; // for FFT
    
    /* make local copy of current Ambi RIR, in WXYZ ordering */
    
    // NOTE: Assumes ACN-N3D
    pvir = malloc1d(nPV * lSig * sizeof(float)); // TODO: pvir to 2D ptr? (float**)malloc2d(nPV, lSig, sizeof(float));
    int xyzOrder[4] = {0, 3, 1, 2}; // w y z x -> w x y z
    for(int i = 0; i < nPV; i++) {
        int inChan = xyzOrder[i];
        memcpy(&pvir[i * lSig], &pData->rirBuf[inChan][0], lSig * sizeof(float));
    };
    float velScale = 1.f / sqrtf(3.f);
    /* scale XYZ to normalized velocity */
    utility_svsmul(&pvir[1 * lSig], &velScale, 3 * lSig, &pvir[1 * lSig]);
    
    /* freq domain diffuseness */
    
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
        memcpy(&pvir_pad[i * lSig_pad + winsize/2],
               &(pvir[i * lSig]),
               lSig * sizeof(float));
    }
    
    /* transform window (symmetric Hann - 'hanning' in matlab) */
    win = malloc1d(winsize * sizeof(float));
    for(int i = 0; i < winsize; i++)
        win[i] = powf( sinf((float)i * (M_PI / (float)winsize)), 2.0f );
    
    /* mem alloc for diffuseness of each window */
    int nDiffFrames = (int)((lSig + 2 * winsize) / hopsize + 0.5f);
    diff_win = calloc1d(nDiffFrames, sizeof(float));
    
    /* Main processing */
    
    strcpy(pData->progressText,"HOSIRR - Rendering");
    
    // mem alloc for a single frame of processing
    insig_win    = calloc1d(fftsize, sizeof(float));
    inspec_anl   = calloc1d(nPV * nBins_anl, sizeof(float_complex)); // ma->ca
    inspec_syn   = calloc1d(nPV * nBins_syn, sizeof(float_complex)); // ma->ca
    wxyzspec_win = malloc1d(4 * nBins_anl * sizeof(float_complex));
    saf_rfft_create(&hFFT_syn, fftsize);
    saf_rfft_create(&hFFT_anl, fftsize/2);
    
    /* window-hopping loop */
    int idx = 0;
    int frameCount = 0;
    while (idx + winsize < lSig + 2 * winsize)
    {
        /* update progress */
        pData->progress0_1 = (float)idx / (float)(lSig + 2 * winsize);
 
        /* Window input and transform to frequency domain */
        for(int i = 0; i < nPV; i++){
            for(int j = 0; j < winsize; j++)
                insig_win[j] = win[j] * pvir_pad[i*lSig_pad + idx + j];
            // fft
            saf_rfft_forward(hFFT_syn, insig_win, &inspec_syn[i * nBins_syn]);
            // trim to just analysis bins
            for(int j = 0, k = 0; j < nBins_anl; j++, k += fftsize/winsize)
                inspec_anl[i * nBins_anl + j] = inspec_syn[i * nBins_syn + k];
        }
        
        for(int i = 0; i < nPV; i++)
            memcpy(&wxyzspec_win[i * nBins_anl],
                   &inspec_anl[i * nBins_anl],
                   nBins_anl * sizeof(float_complex));
                    
        /* Compute broad-band active-intensity vector */
        // TODO: This could perhaps be simplified, since we're not really using the full covariance matrix
        const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
        float intensity[3];
        float energy = 0.0f;
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, 4, 4, maxDiffFreq_idx+1, &calpha,
                    wxyzspec_win, nBins_anl,
                    wxyzspec_win, nBins_anl, &cbeta,
                    FLATTEN2D(pvCOV), 4);
        for(int i = 0; i < 3; i++)
            intensity[i] = crealf(pvCOV[1+i][0]);
        for(int i = 0; i < 4; i++)
            energy += crealf(pvCOV[i][i])*0.5f;
        
        /* Estimating and time averaging of boadband diffuseness */
        float ixyz_smooth[3];
        float energy_smooth;
        float prev_ixyz_smooth[3] = {0.f, 0.f, 0.f};
        float prev_energy_smooth = 0.f;
        float normIntensity_smooth = 0.0f;
        
        for(int i = 0; i < 3; i++) {
            ixyz_smooth[i] = ((1.0f-ALPHA_DIFF_COEFF) * intensity[i]) + (ALPHA_DIFF_COEFF * prev_ixyz_smooth[i]);
            prev_ixyz_smooth[i] = ixyz_smooth[i];
            normIntensity_smooth += powf( fabsf(ixyz_smooth[i]), 2.0f ); // TODO: abs unnecessary?
        }
        energy_smooth = ((1.0f-ALPHA_DIFF_COEFF) * energy) + (ALPHA_DIFF_COEFF * prev_energy_smooth);
        prev_energy_smooth = energy_smooth;
        normIntensity_smooth = sqrtf(normIntensity_smooth);
        // store broadband diffuseness value
        diff_win[frameCount] = 1.0f - (normIntensity_smooth / (energy_smooth + 2.23e-10f));
        
        idx += hopsize;
        frameCount++;
    }
    
    /* diffuse onset */
    
    /* index of max(diffuseness) */
    int maxIdx;
    utility_simaxv(diff_win, nDiffFrames, &maxIdx);
    
    /* index of first index above threshold */
    float maxVal = diff_win[maxIdx];
    float onsetThresh = maxVal * thresh_fac;
    int onsetWinIdx = hosirrlib_firstIndexGreaterThan(diff_win, 0, nDiffFrames-1, onsetThresh);
    pData->diffuseOnsetIdx = onsetWinIdx * hopsize + ((int)hopsize/2); // place onset in the middle of the analysis window
    
    pData->analysisStage = DIFFUSENESS_ONSET_FOUND;

    // TODO: free things
    free(pvir);
    free(pvir_pad);
    free(win);
    free(insig_win);
    free(diff_win);
    free(inspec_anl);
    free(inspec_syn);
    free(wxyzspec_win);
    saf_rfft_destroy(&hFFT_syn);
    saf_rfft_destroy(&hFFT_anl);
}


void hosirrlib_splitBands(
                          void* const hHS,
                          float** const inBuf,
                          float*** const bndBuf,
                          int removeFiltDelayFLAG,
                          ANALYSIS_STAGE thisStage)
{
    /* removeFiltDelay = true will truncate the head and tail of the
     * filtered signal, instead of (only) the tail. So either way,
     * the written output will always be pData->nSamp for each channel. */
    
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("splitBands called.\n");
    
    /* Check previous stages are complete */
    if (pData->analysisStage < thisStage-1) {
        printf("splitBands called before previous stages were completed: %d\n", pData->analysisStage);
        return; // TODO: handle fail case
    }
    
    const int nInSamp   = pData->nSamp; // TODO: todo assumes nsamp of outbuf == nsamp of RIR
    const int filtOrder = pData->bandFiltOrder;
    
    int startCopyIdx;
    if (removeFiltDelayFLAG) {
        startCopyIdx = (int)(filtOrder / 2.f); // length of filter delay (floor)
    } else {
        startCopyIdx = 0;
    }
    
    /* Apply filterbank to rir_bands.
     * Because of fftconv's indexing, it expects filter coefficients
     * for each channel. So we do one channel and band at a time.
     * The output size will be nSamp + filtOrder.
     * See fftconv(): y_len = x_len + h_len - 1; */
    float* temp = malloc1d((nInSamp + filtOrder) * sizeof(float));
    
    for(int ish = 0; ish < pData->nSH; ish++) {
        for(int ib = 0; ib < pData->nBand; ib++) {
            fftconv(&inBuf[ish][0],            // input: 1 channel at a time
                    &pData->H_bandFilt[ib][0],  // band filter coeffs
                    nInSamp,                    // input length
                    filtOrder + 1,              // filter length
                                                // add 1 for internal rounding and compensate for out length calc y_len = x_len + h_len - 1; (?)
                    1,                          // 1 channel at a time
                    temp);
            
            /* Copy temp to out */
            memcpy((void*)&bndBuf[ib][ish][0],
                   &temp[startCopyIdx],
                   nInSamp * sizeof(float));
        }
    }
    pData->analysisStage = thisStage;
    free(temp);
}


void hosirrlib_beamformRIR(
                           void* const hHS)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("beamformRIR called.\n");
    
    /* Check previous stages are complete */
    if (pData->analysisStage < BEAMFORMED-1){
        printf("beamformRIR called before previous stages were completed: %d\n", pData->analysisStage);
        return; // TODO: handle fail case
    }
    
    const int nBeam     = pData->nDir;
    const int nSamp     = pData->nSamp;
    const int nBand     = pData->nBand;
    const int nSH       = pData->nSH;
    const int shOrder   = pData->shOrder;
    
    float * c_l = malloc1d((shOrder + 1) * sizeof(float)); // z beam coeffs, only shOrder + 1 are used
    
    /* Apply beamformer */
    // TODO: These beamforming coeffs only need to be updated if nSH (input) or the spherical design (output) changes, so could be refactored
    /* Calculate beamforming coeffients */
    for (int ib = 0; ib < nBeam; ib++) {
        switch(pData->beamType){
            case STATIC_BEAM_TYPE_CARDIOID:
                beamWeightsCardioid2Spherical(shOrder,  c_l); break;
            case STATIC_BEAM_TYPE_HYPERCARDIOID:
                beamWeightsHypercardioid2Spherical(shOrder,  c_l); break;
            case STATIC_BEAM_TYPE_MAX_EV:
                beamWeightsMaxEV(shOrder,  c_l); break;
        }
        rotateAxisCoeffsReal(shOrder,
                             (float*)c_l,
                             (SAF_PI / 2.0f) - (pData->loudpkrs_dirs_deg[ib][1] * SAF_PI / 180.0f),
                             pData->loudpkrs_dirs_deg[ib][0] * SAF_PI / 180.0f,
                             &pData->decBeamCoeffs[ib][0]);
    }
    /* Apply beam weights
    // float beamWeights[MAX_NUM_BEAMS][MAX_NUM_SH_SIGNALS]; A
    // float prev_SHFrameTD[MAX_NUM_SH_SIGNALS][BEAMFORMER_FRAME_SIZE]; B
    // float outputFrameTD[MAX_NUM_BEAMS][BEAMFORMER_FRAME_SIZE]; C
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nBeams, BEAMFORMER_FRAME_SIZE, nSH, 1.0f,
                (const float*)pData->beamWeights, MAX_NUM_SH_SIGNALS,       // A, 1st dim of A (row-major)
                (const float*)pData->prev_SHFrameTD, BEAMFORMER_FRAME_SIZE, // B, 1st dim of B
                0.0f,                                                       // beta scalar for C
                (float*)pData->outputFrameTD, BEAMFORMER_FRAME_SIZE);       // C, 1st dim of C
     */
    /* Apply beam weights
     // (ref. beamformer_process())
     // A: float** encBeamCoeffs;  // nDir  x nSH
     // B: float*** rirBuf_bands;  // nBand x nSH  x nSamp
     // C: float*** rirBuf_beams;  // nBand x nDir x nSamp */
    for (int bd = 0; bd < nBand; bd++) {
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    nBeam, nSamp, nSH, 1.0f,
                    (const float*)&pData->decBeamCoeffs[0][0], nSH, // A, "1st" dim of A (row-major)
                    (const float*)&pData->rirBuf_bands[bd][0][0], nSamp, // B, 1st dim of B
                    0.0f, // beta scalar for C
                    (float*)&pData->rirBuf_beams[bd][0][0], nSamp); // C, 1st dim of C
    }
    pData->analysisStage = BEAMFORMED;
    
    free(c_l);
}


void hosirrlib_calcEDC(
                       void* const hHS,
                       float*** const inBuf,
                       float*** const edcBuf,
                       ANALYSIS_STAGE thisStage)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("calcEDC called.\n");
    
    /* Check previous stages are complete */
    if (pData->analysisStage < thisStage-1){
        printf("calcEDC called before previous stages were completed: %d\n", pData->analysisStage);
        return; // TODO: handle fail case
    }
    
    const int nSamp = pData->nSamp;
    const int nDir  = pData->nDir;
    const int nBand = pData->nBand;
    const int nSH   = pData->nSH;
    
    /* Copy the RIR band beams into EDC buffer */
    utility_svvcopy(FLATTEN3D(inBuf),
                    nBand * nDir * nSamp,
                    FLATTEN3D(edcBuf));
    
    float*** const edc = edcBuf;
    double sum = 0.0; // TODO: double?

    /* EDC: reverse cummulative sum of signal energy, one channel at a time */
    for (int bd = 0; bd < nBand; bd++) {
        for (int ch = 0; ch < nDir; ch++) {
            sum = 0.0;
            // reverse iterate for backwards cumulative sum
            for (int i = nSamp - 1; i > -1; i--) {
                // in-place: replace signal with its EDC
                // TODO: optim by vectorizing
                sum += edc[bd][ch][i] * edc[bd][ch][i]; // energy
                edc[bd][ch][i] = (float)(10.0 * log10(sum)); // store in dB
            }
        }
    }
    pData->analysisStage = thisStage;
}


void hosirrlib_calcT60(void* const hHS, const float startDb, const float spanDb, const int beginIdx)
{ /*
   startDb  : measure the T60 starting at this level of decay
              (after beginIdx), specify as negative (<= 0)
   spanDb   : measure the T60 over this dB decay span (specify as positive)
   beginIdx : start the search for measurement bounds at this point onward
              e.g. a sample index after the first arrival
              0 for the beginning of the EDC buffer
   */
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("calcT60 called.\n");
    
    /* Check previous stages are complete */
    if (pData->analysisStage < T60_DONE-1){
        printf("calcT60 called before previous stages were completed: %d\n", pData->analysisStage);
        return; // TODO: handle fail case
    }
    
    const int nSamp = pData->nSamp;
    const int nChan = pData->nDir;
    const int nBand = pData->nBand;
    const int nSH   = pData->nSH;
    const float fs  = pData->fs;
    
    float ***edc = pData->edcBuf_rir; // TODO: this will become an arg when/if FDN needs to be measured
    
    float* const x_slope1 = malloc1d(nSamp * sizeof(float)); // vector with slope of 1/samp
    float* const y_edc0m  = malloc1d(nSamp * sizeof(float)); // zero-mean EDC
    float* const stage    = malloc1d(nSamp * sizeof(float)); // staging buffer
    int*** const st_end_meas = (int***)malloc3d(nBand, nChan, 2, sizeof(int)); // start and end samples to measure t60 in each channel

    /* find the start and end points of the measurement */
    for (int bd = 0; bd < nBand; bd++) {
        for (int ch = 0; ch < nChan; ch++) {
            float edcMax = edc[bd][ch][beginIdx];
            
            int start_t60 = hosirrlib_firstIndexLessThan(&edc[bd][ch][0],
                                                         beginIdx, nSamp - 1,
                                                         edcMax + startDb); // startDb is negative
            if (start_t60 < 0) {
                // no value found below the start level // TODO: warn or error
                st_end_meas[bd][ch][0] = 0; // fall back to the head of the EDC
            } else {
                st_end_meas[bd][ch][0] = start_t60;
            }
            
            int end_t60 = hosirrlib_firstIndexLessThan(&edc[bd][ch][0],
                                                       start_t60, nSamp - 1,
                                                       edcMax + startDb - spanDb); // startDb is negative, spanDb is positive
            if (end_t60 < 0) {
                // no value found below startDb - spanDb // TODO: warn or error
                // fall back to near the end of the EDC.. this will likely be quite inaccurate for high frequency bands!
                st_end_meas[bd][ch][1] = (int)(0.7f * nSamp);
            } else {
                st_end_meas[bd][ch][1] = beginIdx + end_t60;
            }
        }
    }
    
    /* Calculate the line of best fit
        // x0m : zero-mean vector of a line with a slope of 1/sample
        // y0m : vector of edc values (over measurement span), with mean removed
        x0m      = vec_idc - mean(vec_idc);
        y0m      = edc_span - mean(edc_span);
        dpc_db   = sum(x0m .* y0m) / sum(x0m.^2);   // decay per sample (dB)
        t60_meas = (-60 / dpc_db) / fs;             //  T60
     */
    
    /* Measure the t60 */
    for (int ib = 0; ib < nBand; ib++) {
        for (int ich = 0; ich < nChan; ich++) {
            int st_meas    = st_end_meas[ib][ich][0];
            int end_meas   = st_end_meas[ib][ich][1];
            int nSamp_meas = end_meas - st_meas + 1;
            
            float y_mean = sumf(&edc[ib][ich][st_meas], nSamp_meas) / nSamp_meas;
            
            /* Construct a vector with a slope of 1:samp with zero mean */
            float first_val = (nSamp_meas - 1) / 2.f; // first value
            for (int i = 0; i < nSamp_meas; i++) {
                x_slope1[i] = first_val + i;
            };
            // remove mean from EDC, within the measurement span
            utility_svssub(&edc[ib][ich][st_meas], &y_mean, nSamp_meas, y_edc0m);
            // covariance x * y and x * x
            utility_svvmul(x_slope1, y_edc0m, nSamp_meas, stage);
            float c_xy = sumf(stage, nSamp_meas);
            utility_svvmul(x_slope1, x_slope1, nSamp_meas, stage);
            float c_xx = sumf(stage, nSamp_meas);
            // slope dB/samp
            float dbPerSamp = c_xy / c_xx;
            // write out
            pData->t60Buf[ib][ich] = -60.f / dbPerSamp / fs;
            // debug
            printf("t60: dir %d band %d  %.2f sec\n", ich, ib, pData->t60Buf[ib][ich]);
        }
    }
        
    free(x_slope1);
    free(y_edc0m);
    free(stage);
    free(st_end_meas);
    pData->analysisStage = T60_DONE;
}

int hosirrlib_firstIndexLessThan(float* vec, int startIdx, int endIdx, float thresh)
{
//    printf("irstIndexLessThan called.\n");
    
    for (int i = startIdx; i < endIdx+1; i++) {
        if (vec[i] < thresh)
            return i;
    }
    return -1;
}

int hosirrlib_firstIndexGreaterThan(float* vec, int startIdx, int endIdx, float thresh)
{
//    printf("firstIndexGreaterThan called.\n");
    
    for (int i = startIdx; i < endIdx+1; i++) {
        if (vec[i] > thresh)
            return i;
    }
    return -1;
}

// for the GUI to display the EDCs
void hosirrlib_copyNormalizedEDCBufs(void* const hHS, float** edcCopy, float displayRange)
{
    /*
     Note edcBuf_rir are foat*** nband x ndir x nsamp, but dirEDC pointer
     in the UI is float** ndir x nsamp,
     so for now, just return the first band of each direction
     */
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    const int nBand = pData->nBand;
    const int nSamp = pData->nSamp;
    const int nDir = pData->nDir;
    float*** const edc = pData->edcBuf_rir;
    float* const edc_flat = (float*)FLATTEN3D(pData->edcBuf_rir);
    const int flatLen = nBand * nDir * nSamp;
    
    if(pData->analysisStage >= RIR_EDC_DONE) {    // confirm EDCs are rendered
        
        /* normalise to range [-1 1] for plotting */
        float maxVal, minVal, range, add, scale, sub;
        // intializse to first/last value of first channel
        maxVal = edc[0][0][0];
        minVal = maxVal - displayRange; // just display the top displayRange in dB
        // minVal = edc[0][0][nSamp-1];
        // check the first and last values of every channel
        for (int ib = 1; ib < nBand; ib++)
            for (int ich = 0; ich < nDir; ich++) {
                if (edc[ib][ich][0] > maxVal)
                    maxVal = edc[ib][ich][0];
                // if (edc[ib][ich][nSamp-1] < minVal)
                //     minVal = edc[ib][ich][nSamp-1];
            }
        
        range = maxVal - minVal;
        add = minVal * -1.f;
        scale = 2.0f/fabsf(range);
        sub = 1.f;
        printf("max %.1f, min %.1f, rng %.1f, add %.1f, scl %.1f, sub %.1f, ",
               maxVal, minVal, range, add, scale, sub);
        
        for(int i = 0; i < pData->nDir; i++) {
            utility_svsadd(&(pData->edcBuf_rir[0][i][0]), // [bnd][ch][smp]
                           &add, nSamp,
                           &edcCopy[i][0]); // [ch][smp]
            utility_svsmul(&edcCopy[i][0],  // [ch][smp]
                           &scale, nSamp,
                           &edcCopy[i][0]); // [ch][smp]
            utility_svssub(&edcCopy[i][0],  // [ch][smp]
                           &sub, nSamp,
                           &edcCopy[i][0]); // [ch][smp]
        }
        // simple copy
        //        for(int i = 0; i < pData->nDir; i++)
        //            memcpy(&edcCopy[i][0],                    // copy-to channel
        //                   &(pData->edcBuf_rir[0][i][0]), // [bnd][ch][smp]
        //                   nSamp * sizeof(float));
    }
}

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
    int o[MAX_SH_ORDER+2];
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
    NORM_TYPES norm;
    CH_ORDER chOrdering;
    
    /* Check if processing should actually go-ahead */
    if(pData->ambiRIR_status != AMBI_RIR_STATUS_LOADED ||
       pData->lsRIR_status == LS_RIR_STATUS_RENDERED ||
       pData->lsRIR_status == LS_RIR_STATUS_RENDEREDING_ONGOING)
        return;
    else
        pData->lsRIR_status = LS_RIR_STATUS_RENDEREDING_ONGOING;
    
    /* take a local copy of current configuration to be thread safe */
    fs = pData->ambiRIRsampleRate;
    order = SAF_MIN(pData->analysisOrder, pData->ambiRIRorder);
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
        case CH_ACN:
//            memcpy(shir, pData->shir, nSH * lSig * sizeof(float));
            memcpy(shir, FLATTEN2D(pData->shir), nSH * lSig * sizeof(float));
            break;
        case CH_FUMA:
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
    for(n=0; n<MAX_SH_ORDER+2; n++){  o[n] = n*n;  }
    switch(norm){
        case NORM_N3D:  /* already in N3D, do nothing */
            break;
        case NORM_FUMA: /* (same as converting SN3D->N3D for first-order) */
        case NORM_SN3D: /* convert to N3D */
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
    pData->ambiRIRorder = SAF_MIN(sqrt(numChannels-1), MAX_SH_ORDER);
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
    pData->windowLength = SAF_CLAMP(pData->windowLength, MIN_WINDOW_LENGTH, MAX_WINDOW_LENGTH);
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
    pData->analysisOrder = SAF_MIN(SAF_MAX(newValue,1), MAX_SH_ORDER);
    /* FUMA only supports 1st order */
    if(pData->analysisOrder!=ANALYSIS_ORDER_FIRST && pData->chOrdering == CH_FUMA)
        pData->chOrdering = CH_ACN;
    if(pData->analysisOrder!=ANALYSIS_ORDER_FIRST && pData->norm == NORM_FUMA)
        pData->norm = NORM_SN3D;
    pData->lsRIR_status = LS_RIR_STATUS_NOT_RENDERED;
}

void hosirrlib_setLoudspeakerAzi_deg(void* const hHS, int index, float newAzi_deg)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    if(newAzi_deg>180.0f)
        newAzi_deg = -360.0f + newAzi_deg;
    newAzi_deg = SAF_MAX(newAzi_deg, -180.0f);
    newAzi_deg = SAF_MIN(newAzi_deg, 180.0f);
    pData->loudpkrs_dirs_deg[index][0] = newAzi_deg;
    pData->lsRIR_status = LS_RIR_STATUS_NOT_RENDERED;
}

void hosirrlib_setLoudspeakerElev_deg(void* const hHS, int index, float newElev_deg)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    newElev_deg = SAF_MAX(newElev_deg, -90.0f);
    newElev_deg = SAF_MIN(newElev_deg, 90.0f);
    pData->loudpkrs_dirs_deg[index][1] = newElev_deg;
    pData->lsRIR_status = LS_RIR_STATUS_NOT_RENDERED;
}

void hosirrlib_setNumLoudspeakers(void* const hHS, int new_nLoudspeakers)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    pData->nLoudpkrs = new_nLoudspeakers > MAX_NUM_LOUDSPEAKERS ? MAX_NUM_LOUDSPEAKERS : new_nLoudspeakers;
    pData->nLoudpkrs = SAF_MAX(MIN_NUM_LOUDSPEAKERS, pData->nLoudpkrs);
    pData->lsRIR_status = LS_RIR_STATUS_NOT_RENDERED;
}

void  hosirrlib_setOutputConfigPreset(void* const hHS, int newPresetID)
{
    hosirrlib_data *pData = ( hosirrlib_data*)(hHS);
    loadLoudspeakerArrayPreset(newPresetID, pData->loudpkrs_dirs_deg, &(pData->nLoudpkrs));
    pData->lsRIR_status = LS_RIR_STATUS_NOT_RENDERED;
}

void hosirrlib_setChOrder(void* const hHS, int newOrder)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    if((CH_ORDER)newOrder != CH_FUMA || pData->analysisOrder==ANALYSIS_ORDER_FIRST) /* FUMA only supports 1st order */
        pData->chOrdering = (CH_ORDER)newOrder;
    pData->lsRIR_status = LS_RIR_STATUS_NOT_RENDERED;
}

void hosirrlib_setNormType(void* const hHS, int newType)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    if((NORM_TYPES)newType != NORM_FUMA || pData->analysisOrder==ANALYSIS_ORDER_FIRST) /* FUMA only supports 1st order */
        pData->norm = (NORM_TYPES)newType;
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
