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

#include "hosirr_internal.h"

void hosirrlib_create(
                      void ** const phHS
                      )
{
    hosirrlib_data* pData = (hosirrlib_data*)malloc1d(sizeof(hosirrlib_data));
    *phHS = (void*)pData;

    /* new hodecaylib */
    
    // depend only on nDir
    pData->encBeamCoeffs        = NULL;     // nSH x nDir
    pData->decBeamCoeffs        = NULL;     // nDir x nSH
    pData->dirGainBufDB         = NULL;     // nDir x nBand
    pData->t60Buf_omni          = NULL;     // nBand x 1
    pData->t60Buf_dir           = NULL;     // nDir x nBand
    pData->rdrBuf               = NULL;     // nBand x 1
    pData->directOnsetIdx_bnd   = NULL;     // nBand x 1
    
    // depend on output design (nDir) AND input RIR (nSamp)
    pData->rirBuf_sh            = NULL;     // nSH x nSamp
    pData->rirBuf_bnd_dir       = NULL;     // nDir x nSamp
    pData->fdnBuf_dir           = NULL;     // nDir x nSamp
    pData->fdnBuf_bnd_dir       = NULL;     // nBand x nDir x nSamp
    pData->edcBufOmn_bnd        = NULL;     // nBand x nSamp
    pData->edcBuf_bnd_dir       = NULL;     // nDir x nBand x nSamp
    pData->edcBufFDN_bnd_dir    = NULL;     // nDir x nBand x nSamp
    pData->fdnBuf_sh            = NULL;     // nSH x nSamp
    
    pData->H_bandFilt           = NULL;     // nBand x filtOrder+1
    pData->bandCenterFreqs      = NULL;     // 1 x nBand
    pData->bandXOverFreqs       = NULL;     // 1 x nBand-1
    pData->srcDirectivity       = NULL;     // 1 x nBand
    
    /* Constants */
    
    /* If diffuseness never crosses this threshold, diffuse onset
     * defaults to direct onset +5ms. */
    pData->diffuseMin = 0.3f;
    pData->diffuseOnsetFallbackDelay = 0.005;
    
    /* Defaults */
    for (int i = 0; i < 3; i++) {
        pData->srcPosition[i] = 0.f;
        pData->recPosition[i] = 0.f;
    }
    pData->srcPosition[0] = 1.f; // source defaults to 1 m in front of receiver
    
    pData->srcDirectivityFlag = 1; // default to "regular" loudspeaker directivity
        
    /* Zero out the state of buffer resources */
    hosirrlib_setUninitialized(pData);
    
    /* Beam shape for decomposition */
    pData->beamType = HYPERCARDIOID_BEAM;
    
    /* Initialize spherical design for directional analysis */
    // TODO: SPHDESIGNs coincide with LS arrays for now, not a clean separation
    // default design currently just set by enum SPHDESIGN_ARRAY_PRESET_15PX
    loadSphDesignPreset(
                        SPHDESIGN_PRESET_15PX,
                        pData->loudpkrs_dirs_deg,
                        &(pData->nLoudpkrs));

    pData->nDir = pData->nLoudpkrs; // TODO: update
    
    /* NOTE: Don't initialize filters yet, input RIR is required for getting
     * the fsfilterbank constants (fs) set in _initBandFilters(). */
    
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
        
        // depend only on output design (nDir)
        free(pData->encBeamCoeffs);
        free(pData->decBeamCoeffs);
        free(pData->dirGainBufDB);
        free(pData->t60Buf_omni);
        free(pData->t60Buf_dir);
        free(pData->rdrBuf);
        
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
        free(pData->srcDirectivity);
        
        /* original hosirrlib */
        free(pData->shir);
        free(pData->lsir);
        free(pData->progressText);
        
        free(pData);
        pData = NULL;
    }
}

void hosirrlib_setInputNorm(void* const hHS, int newType)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    if((NORMALIZATION_TYPES)newType != FUMA_NORM)    /* FUMA not currently supported */
    {
        if (newType < 4) {                               /* only 3 input norm types */
            pData->inputNorm = (NORMALIZATION_TYPES)newType;
        } else {
            hosirr_print_error("Input norm specification is outside the valid range (1..3).");
        }
    } else {
        hosirr_print_warning("FUMA is not currently supported as an input format.");
        return;
    }
}

/* Note: pData->inNorm can't be inferred from the should be set prior to calling this. */
int hosirrlib_setRIR(
                     void* const hHS,
                     const float** H,
                     int numChannels,
                     int numSamples,
                     int sampleRate,
                     int inNormInt // N3D, SN3D, FUMA : 1, 2, 3 
                     )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    /* Check channel count to see if input is actually in the SHD */
    if (fabsf(sqrtf((float)numChannels) - floorf(sqrtf((float)numChannels))) > 0.0001f) {
        
        // TODO: only set state to unitialized if there was no previous RIR loaded
        // setRIRState_uninitialized(pData)
        printf("setRIR: Based on numChannels, input file doesn't appear to be a SHD file\n");
        
        /* old vars */
        pData->ambiRIR_status = AMBI_RIR_STATUS_NOT_LOADED;
        pData->ambiRIRorder = -1;
        pData->ambiRIRlength_seconds = pData->ambiRIRlength_samples = 0.0f;
        pData->ambiRIRsampleRate = 0;
        return (int)(pData->ambiRIR_status);
    }
    
    pData->nSH      = numChannels;
    pData->nSamp    = numSamples;
    pData->shOrder  = HOSIRR_MIN(sqrt(numChannels-1), HOSIRR_MAX_SH_ORDER);
    pData->fs       = (float)sampleRate;
    pData->duration = numSamples / (float)sampleRate;
    
    NORMALIZATION_TYPES inNorm = (NORMALIZATION_TYPES)inNormInt;
    hosirrlib_setInputNorm(pData, inNorm);
    
    /* (Re)alloc and copy in input RIR */
    pData->rirBuf_sh = (float**)realloc2d((void**)pData->rirBuf_sh, numChannels, numSamples, sizeof(float));
    
    /* convert to N3D if needed */
    switch (pData->inputNorm) {
        case N3D_NORM:  /* already in N3D, just copy it in */
            printf("\n\tProcessing N3D\n"); // dbg
            for(int i = 0; i < numChannels; i++)
                utility_svvcopy(H[i], numSamples, pData->rirBuf_sh[i]);
            break;
        case SN3D_NORM: /* convert to N3D */
            printf("\n\tProcessing SN3D\n"); // dbg
            for (int n = 0; n < pData->shOrder+1; n++) {
                int numOrderChans = n*2 + 1;
                int orderBaseIdx = n * n;
                for (int i = 0; i < numOrderChans; i++) {
                    float orderScale = sqrtf(2.0f * (float)n + 1.0f);
                    utility_svsmul((float *)H[orderBaseIdx+i],
                                   &orderScale, numSamples,
                                   pData->rirBuf_sh[orderBaseIdx+i]);
                }
            }
            break;
        case FUMA_NORM: /* Not implemented */
            hosirr_print_error("FUMA_NORM isn't currently supported as an input format");
            break;
    }
    
    pData->analysisStage = RIR_LOADED;
    
    /* Initialize band-processing filters */
    hosirrlib_initBandProcessing(pData, FILTERS_INTITIALIZED);
    
    /* Alloc processing resources */
    hosirrlib_allocProcBufs(pData, ANALYSIS_BUFS_LOADED);
    
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


void hosirrlib_setSrcPosition(
                              void* const hHS,
                              const float x,
                              const float y,
                              const float z
                              )
{
    printf("\n\tsetSrcPosition called.\n"); // dbg
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    float * sPos = pData->srcPosition;
    sPos[0] = x;
    sPos[1] = y;
    sPos[2] = z;
}

void hosirrlib_setRecPosition(
                              void* const hHS,
                              const float x,
                              const float y,
                              const float z
                              )
{
    printf("\n\tsetRecPosition called.\n"); // dbg
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    float * rPos = pData->recPosition;
    rPos[0] = x;
    rPos[1] = y;
    rPos[2] = z;
}

float hosirrlib_getSrcRecDistance(
                                  void* const hHS
                                  )
{
    //printf("\n\tgetSrcRecDistance called.\n"); // dbg
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    float * rPos = pData->recPosition;
    float * sPos = pData->srcPosition;
    float sqDiff[3];
    for (int i = 0; i < 3; i++) {
        float diff = sPos[i] - rPos[i];
        sqDiff[i] = diff * diff;
    }
    return sqrtf(sumf(sqDiff, 3));
}

void hosirrlib_setUninitialized(
                                void* const hHS
                                )
{
    printf("\nsetUninitialized called.\n"); // dbg
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    // pending initializations
    pData->nSH = -1; // input vars (assume for now in params = out params)
    pData->nSamp = -1;
    pData->shOrder = -1;
    pData->fs = -1.f;
    pData->directOnsetIdx_brdbnd = -1;
    pData->diffuseOnsetIdx = -1;
    pData->diffuseOnsetSec = 0;
    pData->t0 = -1;
    pData->t0Idx = 0;
    pData->duration = 0.0f;
    pData->analysisStage = RIR_NOT_LOADED;
}


void hosirrlib_initBandProcessing(
                                  void* const hHS,
                                  ANALYSIS_STAGE thisStage
                                  )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    checkProperProcessingOrder(pData, thisStage, __func__);

    /* Initialize octave band filters
     * bandCenterFreqs are "center freqs" used to determine crossover freqs.
     * The lowest and highest bands are LP and HP filters, "center freq" for
     * those is a bit of a misnomer. Xover freqs are center freqs * sqrt(2)
     * (octave bands), see hosirrlib_initBandProcessing().
     *      LPF cutoff is 125 * sqrt(2) = 177 Hz.
     *      HPF cuton is 8k * sqrt(2) = 11313 Hz.
     *      Must be the same size as pData->nBand. */
    pData->nBand = 8;
    pData->bandFiltOrder = 400;
    float bandCenterFreqs[8] = { 125.f, 250.f, 500.f, 1000.f, 2000.f, 4000.f, 8000.f, 16000.f };
    
    /* Bandwise directivity factor of a loudspeaker (Genelec_8351A)
     * Calculated with https://github.com/AppliedAcousticsChalmers/sound-source-directivities */
    float srcLSDirectivity[8] = { 1.1216f, 1.8731f, 3.6854f, 5.4908f, 5.1385f, 6.8701f, 5.1735f, 6.3951f };
    
    pData->bandCenterFreqs = malloc1d(pData->nBand * sizeof(float));
    pData->bandXOverFreqs  = malloc1d((pData->nBand-1) * sizeof(float));
    pData->srcDirectivity  = malloc1d(pData->nBand * sizeof(float));
    
    /* Populate the member vars for external access */
    utility_svvcopy(bandCenterFreqs, pData->nBand, pData->bandCenterFreqs);
    utility_svvcopy(srcLSDirectivity, pData->nBand, pData->srcDirectivity);
    
    /* Set xover freqs */
    for (int ib = 0; ib < pData->nBand-1; ib++) {
        pData->bandXOverFreqs[ib] = bandCenterFreqs[ib] * sqrtf(2.f);
        //printf("\tXOver band %d %.1f\n", ib, pData->bandXOverFreqs[ib]); // dbg
    }
    
    /* Create the filterbank */
    if (pData->H_bandFilt == NULL) { // if this is the first time this function is called ...
        
        /* Allocate filter coefficients (nBand x bandFiltOrder + 1) */
        pData->H_bandFilt = (float**)realloc2d((void**)pData->H_bandFilt,
                                               pData->nBand,
                                               pData->bandFiltOrder + 1,
                                               sizeof(float));
        
        /* Compute FIR Filterbank coefficients */
        // FIRFilterbank(pData->bandFiltOrder,          // SAF version: bug at currently linked SAF version
        hosirrlib_FIRFilterbank(pData->bandFiltOrder,   // locally patched version
                                pData->bandXOverFreqs,
                                pData->nBand-1,
                                pData->fs,
                                WINDOWING_FUNCTION_HAMMING,
                                1,
                                FLATTEN2D(pData->H_bandFilt));
    }
    
    // hosirrlib_inspectFilts(pData); // dbg func
    
    pData->analysisStage = thisStage;
}


// (re)allocate the buffers used for storing intermediate processing data
void hosirrlib_allocProcBufs(
                             void * const hHS,
                             ANALYSIS_STAGE thisStage
                             )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    checkProperProcessingOrder(pData, thisStage, __func__);
    
    const int nSH   = pData->nSH;
    const int nDir  = pData->nDir;
    const int nBand = pData->nBand;
    const int nSamp = pData->nSamp;

    // TODO: these don't all need to be allocated in the case of a mono RIR
    
    // These members depend on output design (nDir) AND input RIR (nSH, nSamp)
    pData->rirBuf_bnd_sh    = (float***)realloc3d((void***)pData->rirBuf_bnd_sh,
                                                  nBand, nSH, nSamp, sizeof(float));
    pData->rirBuf_bnd_dir   = (float***)realloc3d((void***)pData->rirBuf_bnd_dir,
                                                  nBand, nDir, nSamp, sizeof(float));
    pData->edcBufOmn_bnd    = (float**)realloc2d((void**)pData->edcBufOmn_bnd,
                                                 nBand, nSamp, sizeof(float));
    pData->edcBuf_bnd_dir   = (float***)realloc3d((void***)pData->edcBuf_bnd_dir,
                                                  nBand, nDir, nSamp, sizeof(float));
    pData->fdnBuf_dir       = (float**)realloc2d((void**)pData->fdnBuf_dir,
                                                 nDir, nSamp, sizeof(float));
    pData->fdnBuf_bnd_dir   = (float***)realloc3d((void***)pData->fdnBuf_bnd_dir,
                                                  nBand, nDir, nSamp, sizeof(float));
    pData->edcBufFDN_bnd_dir = (float***)realloc3d((void***)pData->edcBufFDN_bnd_dir,// TODO: remove?
                                                   nBand, nDir, nSamp, sizeof(float));
    pData->fdnBuf_sh        = (float**)realloc2d((void**)pData->fdnBuf_sh,
                                                 nSH, nSamp, sizeof(float));
    pData->encBeamCoeffs    = (float**)realloc2d((void**)pData->encBeamCoeffs,
                                                 nSH, nDir, sizeof(float));
    pData->decBeamCoeffs    = (float**)realloc2d((void**)pData->decBeamCoeffs,
                                                 nDir, nSH, sizeof(float));
    
    /* (Re)allocate buffers that depend only on the number of output channels
     * OPTIM: These could alternatively only be reallocated when the out design
     * changes, but it's not heavy and it's more concise */
    pData->dirGainBufDB     = (float**)realloc2d((void**)pData->dirGainBufDB,
                                                 nBand, nDir, sizeof(float));
    pData->t60Buf_omni      = (float*)realloc1d((void*)pData->t60Buf_omni,
                                                nBand * sizeof(float));
    pData->t60Buf_dir       = (float**)realloc2d((void**)pData->t60Buf_dir,
                                                 nBand, nDir, sizeof(float));
    pData->rdrBuf           = (float*)realloc1d((void*)pData->rdrBuf,
                                                nBand * sizeof(float));
    pData->directOnsetIdx_bnd  = (int*)realloc1d((void*)pData->directOnsetIdx_bnd,
                                                 nBand * sizeof(int));

    pData->analysisStage = thisStage;
}


void hosirrlib_renderTMP(
                         void * const hHS
                         )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    printf("\nrenderTMP called.\n"); // dbg
    
    /* Check if processing should actually go-ahead */
    if(pData->ambiRIR_status != AMBI_RIR_STATUS_LOADED ||
       pData->lsRIR_status == LS_RIR_STATUS_RENDERED ||
       pData->lsRIR_status == LS_RIR_STATUS_RENDEREDING_ONGOING)
        return;
    else
        pData->lsRIR_status = LS_RIR_STATUS_RENDEREDING_ONGOING;
    
    strcpy(pData->progressText,"Processing");
    
    hosirrlib_processRIR(pData, DIRGAIN_DONE);
    
    /* indicate that rendering is complete */
    pData->progress0_1 = 1.0f;
    pData->lsRIR_status = LS_RIR_STATUS_RENDERED;
}


void hosirrlib_processRIR(
                          void * const hHS,
                          ANALYSIS_STAGE endStage
                          )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    ANALYSIS_STAGE requiredPreviousStage = ANALYSIS_BUFS_LOADED;
    checkProperProcessingOrder(pData, requiredPreviousStage, __func__);
    
    // TODO: set constants elsewhere
    const float directOnsetThreshDb     = -3.f;     // direct arrival onset threashold (dB below omni peak)
    const float diffuseOnsetThresh_shd  = 0.707f;   // diffuse onset threashold (normalized scalar below diffuseness peak)
    const float diffuseOnsetThresh_mono = 0.707f;   // diffuse onset threashold (normalized scalar below diffuseness peak)
    const float t60_start_db            = -2.f;     // starting level of t60 measurement (<= 0)
    const float t60_span_db             = 15.f;     // db falloff span over which t60 is measured (> 0)
    const int   nWin_smooth             = 3;        // number of analysis frames to smooth diffuseness value (framesize: 128, hop: 64)
    const int   maxDirGainAdjustment    = 12;       // maximum ajustment of directional gain allowed (absolute value)
    
    const int nBand = pData->nBand;
    const int nDir  = pData->nDir;
    const int nSamp = pData->nSamp;
    
    
    /* Omni, bandwise analysis */
    
    if (pData->analysisStage+1 > endStage)
        return;
    
    hosirrlib_splitBands(pData, pData->rirBuf_sh, pData->rirBuf_bnd_sh,
                         1,
                         BANDS_SPLIT);
    if (pData->analysisStage+1 > endStage)
        return;
    
    hosirrlib_setDirectOnsetIndices(pData, &pData->rirBuf_sh[0][0], pData->rirBuf_bnd_sh,
                                    directOnsetThreshDb,
                                    DIRECT_ONSETS_FOUND);
    if (pData->analysisStage+1 > endStage)
        return;
    
    if (pData->shOrder > 0) {
        // SHD SRIR
        hosirrlib_setDiffuseOnsetIndex_shd(pData,
                                           pData->rirBuf_sh,
                                           diffuseOnsetThresh_shd,
                                           nWin_smooth,
                                           DIFFUSENESS_ONSET_FOUND);
    } else {
        // Mono RIR
        hosirrlib_setDiffuseOnsetIndex_mono(pData,
                                            &pData->rirBuf_sh[0][0],        // ptr to head of first channel (0th-order shd (2D) signal)
                                            diffuseOnsetThresh_mono,
                                            1024,                           // window size
                                            64,                             // hop size
                                            pData->directOnsetIdx_brdbnd,   // start index (center of first window)
                                            DIFFUSENESS_ONSET_FOUND);
    }
    if (pData->analysisStage+1 > endStage)
        return;
  
    hosirrlib_calcEDC_omni(pData, pData->rirBuf_bnd_sh, pData->edcBufOmn_bnd,
                           nBand, nSamp,
                           EDC_OMNI_DONE);
    if (pData->analysisStage+1 > endStage)
        return;
    
    hosirrlib_calcT60_omni(pData, pData->edcBufOmn_bnd, pData->t60Buf_omni,
                           nBand, nSamp,
                           t60_start_db, t60_span_db, pData->diffuseOnsetIdx,
                           T60_OMNI_DONE);
    if (pData->analysisStage+1 > endStage)
        return;
    
    // Requires: omni broadband diffuse onset, omni bandwise direct onsets, omni bandwise T60s
    hosirrlib_calcRDR(pData, pData->rirBuf_bnd_sh, pData->rdrBuf,
                      nBand, nSamp,
                      pData->diffuseOnsetIdx,
                      pData->directOnsetIdx_bnd,
                      pData->t60Buf_omni,
                      pData->srcDirectivityFlag,
                      RDR_DONE);
    if (pData->analysisStage+1 > endStage)
        return;
    
    
    /* Directional bandwise analysis */
    
    if (pData->shOrder > 0) { // Only for shd signals, not mono
        
        hosirrlib_beamformRIR(pData, pData->rirBuf_bnd_sh, pData->rirBuf_bnd_dir,
                              BEAMFORMED);
        if (pData->analysisStage+1 > endStage)
            return;
        
        hosirrlib_calcEDC_beams(pData, pData->rirBuf_bnd_dir, pData->edcBuf_bnd_dir,
                                nBand, nDir, nSamp,
                                EDC_DIR_DONE);
        if (pData->analysisStage+1 > endStage)
            return;
        
        hosirrlib_calcT60_beams(pData, pData->edcBuf_bnd_dir, pData->t60Buf_dir,
                                nBand, nDir, nSamp,
                                t60_start_db, t60_span_db, pData->diffuseOnsetIdx,
                                T60_DIR_DONE);
        if (pData->analysisStage+1 > endStage)
            return;
        
        hosirrlib_calcDirectionalGainDB(pData, pData->dirGainBufDB,
                                        t60_start_db, t60_span_db,
                                        pData->diffuseOnsetIdx,
                                        maxDirGainAdjustment,
                                        DIRGAIN_DONE);
        if (pData->analysisStage+1 > endStage)
            return;
    }
}


/* Find the direct onset
 * NOTE: The naive approach is to consider the index of absolute max value
 * in the buffer. However, with synthetic RIRs, like from a shoebox, coincident
 * reflections may have greater magnitude than the direct arrival.
 * So this returns the index of the first value that crosses the
 * (max abs value - threshold).
 *
 * thresh_dB: threshold (pressure, dB) below the absolute max value in the
 *            buffer, above which the onset is considered to have occured. */
void hosirrlib_setDirectOnsetIndices(
                                     void* const hHS,
                                     float* const brdbndBuf,    // 1 x nSamp
                                     float*** const bndBuf,     // nBand x nSH x nSamp
                                     const float thresh_dB,
                                     ANALYSIS_STAGE thisStage
                                     )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    checkProperProcessingOrder(pData, thisStage, __func__);
    
    const int nSamp = pData->nSamp;
    const int nBand = pData->nBand;
    const float  fs = pData->fs;
    
    /* Buffer to store signal absolute values (1ch at a time) */
    float * const vabs_tmp = malloc1d(nSamp * sizeof(float));
    
    const float srcRecDist = hosirrlib_getSrcRecDistance(hHS);
    int directOnsetIdx = getDirectOnset_1ch(brdbndBuf, vabs_tmp, thresh_dB, nSamp);
    
    /* T0: set based on broadband direct onset */
    const float t0 = ((float)directOnsetIdx / fs) - (srcRecDist / 343.f);
    pData->directOnsetIdx_brdbnd = HOSIRR_MAX(directOnsetIdx, 0);
    pData->t0 = t0; // sec
    pData->t0Idx = (int)(t0 * fs); // A negative t0Idx means t0 is before the file begins (RIR can be trimmed ahead of the direct arrival)
    
    // dbg
    printf("\n     src->rec dist: %.3f m", srcRecDist);
    printf("\ndirect onset index: %d\t(%.4f sec)\n\n", pData->directOnsetIdx_brdbnd, (float)pData->directOnsetIdx_brdbnd / pData->fs);
    printf("\n          t0 index: %d\t(%.4f sec)", pData->t0Idx, t0);
    
    /* Bandwise direct onsets */
    printf("Bandwise direct onset indices:\n");
    for (int ib = 0; ib < nBand; ib++) {
        directOnsetIdx = getDirectOnset_1ch(&bndBuf[ib][0][0], // omni band: nBand x nSH x nSamp
                                            vabs_tmp, thresh_dB, nSamp);
        pData->directOnsetIdx_bnd[ib] = HOSIRR_MAX(directOnsetIdx, 0);
        printf("\t(b%d): %d \t(%.4f sec)\n", ib,  pData->directOnsetIdx_bnd[ib], (float)pData->directOnsetIdx_bnd[ib] / pData->fs); // dbg
    }
    
    pData->analysisStage = thisStage;
    free(vabs_tmp);
}

// returns -1 on fail
int getDirectOnset_1ch(
                       const float * const chan,
                       float * tmp, // buffer to hold copied channel data, NULL uses local memory
                       const float thresh_dB,
                       const int nSamp
                       )
{
    int freeLocalBuf = 0;
    if (tmp == NULL) {
        freeLocalBuf = 1;
        tmp = malloc1d(nSamp * sizeof(float));
    }
    /* Absolute values of the omni channel */
    int maxIdx;
    utility_svabs(chan, nSamp, tmp);     // abs
    utility_simaxv(tmp, nSamp, &maxIdx); // index of max(abs(omni))
    
    /* Index of first index above threshold */
    float maxVal = tmp[maxIdx];
    float onsetThresh = maxVal * powf(10.f, thresh_dB / 20.f);
    const int directOnsetIdx = hosirrlib_firstIndexGreaterThan(tmp, 0, nSamp-1, onsetThresh);
    
    if (freeLocalBuf)
        free(tmp);
    
    return directOnsetIdx;
}


/* thresh_fac: threshold (normalized scalar) below the absolute max value in the
 *             buffer, above which the onset is considered to have occured.
 * NOTE: rirBuf_sh is expected to be ACN-N3D
 */
void hosirrlib_setDiffuseOnsetIndex_shd(
                                        void* const hHS,
                                        float ** const rirBuf_sh,
                                        const float thresh_fac,
                                        const int nWin_smooth,
                                        ANALYSIS_STAGE thisStage
                                        )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    checkProperProcessingOrder(pData, thisStage, __func__);
    
    /* Take a local copy of current configuration to be thread safe */
    const int nChan    = pData->nDir;
    const int nBand    = pData->nBand;
    const float fs     = pData->fs;
    const int order    = pData->shOrder;
    const int lenInSig = pData->nSamp;
    const int winSize  = pData->windowLength;
    const int nPV      = 4;                         // P-V channels: WXYZ only
    
    /* Inits */
    const int hopSize    = winSize/2;               // time-domain hop (for input)
    const int fftSize    = winSize*2;               // fftsize = 2*winsize (input is zero-padded)
    const int nBins_anl  = winSize/2 + 1;           // nBins used for analysis
    const int nBins_syn  = fftSize/2 + 1;           // nBins used for synthesis
    const int lenSig_pad = winSize/2 + lenInSig + winSize*2; // winsize/2 at head, sinsize*2 at tail
    float_complex pvCOV[4][4];
    
    /* Smoothing constant, a refactoring of:
     * tau          = (hopsize * nwin_smooth) / fs; <-- smoothing coefficient
     * smooth_const = expf( -1.f / (tau * fs / (float)hopsize)); <-- "decay" */
    const float smooth_const = expf( -1.f / (float)nWin_smooth);
    
    printf("diffuseness smoothing coeff: %.3f vs const %.3f\n", smooth_const, ALPHA_DIFF_COEFF); //dbg
    
    /* Max freq bin to calculate diffuseness */
    float nearestVal = 10e5f;
    int   maxDiffFreq_idx = 0;
    for(int i = 0; i < nBins_anl; i++){
        // bin_idx * binwidth (Hz)
        // float tmp = fabsf((float)i * (fs / (float)winsize) - MAX_DIFF_FREQ_HZ); // TODO: should be fs/fftsize, not fs/winsize?
        float tmp = fabsf( (float)i * ( fs / (float)fftSize ) - MAX_DIFF_FREQ_HZ );
        if(tmp < nearestVal){
            nearestVal = tmp;
            maxDiffFreq_idx = i;
        }
    }
    
    printf("num diffuseness / analysis bins: %d / %d\n", maxDiffFreq_idx+1, nBins_anl); // dbg
    
    float * pvir, * pvir_pad, * win, * insig_win, * diff_frames;
    float_complex * inspec_anl, * inspec_syn, * pvspec_win;
    void * hFFT_syn; // for FFT
    
    /* Make local copy of the Ambi RIR, from ACN to pressure-velocity (WXYZ) ordering */
    pvir = malloc1d(nPV * lenInSig * sizeof(float));
    int acnOrder[4] = {0, 3, 1, 2};
    for(int i = 0; i < nPV; i++) {
        int inIdx = acnOrder[i]; // w y z x -> w x y z
        memcpy(&pvir[i * lenInSig],
               &rirBuf_sh[inIdx][0],
               lenInSig * sizeof(float));
    };
    
    /* Scale XYZ to normalized velocity from N3D */
    float velScale = 1.f / sqrtf(3.f);
    utility_svsmul(&pvir[1 * lenInSig], &velScale,
                   (nPV-1) * lenInSig, &pvir[1 * lenInSig]);
    
    /* Normalise so the peak of the omni is 1 */
    int peakIdx;
    utility_simaxv(pvir, lenInSig, &peakIdx);       // index of max(abs(omni))
    float peakNorm = 1.0f / fabsf(pvir[peakIdx]);
    utility_svsmul(pvir, &peakNorm, nPV * lenInSig, pvir);
    
    /* Zero pad the signal's start and end for STFT
     * winsize/2 at head, winsize*2 at tail */
    pvir_pad = calloc1d(nPV * lenSig_pad, sizeof(float));
    for(int i = 0; i < nPV; i++) {
        memcpy(&pvir_pad[i * lenSig_pad + (winSize/2)],
               &(pvir[i * lenInSig]),
               lenInSig * sizeof(float));
    }
    
    /* Transform window (symmetric Hann - 'hanning' in matlab) */
    win = malloc1d(winSize * sizeof(float));
    for(int i = 0; i < winSize; i++)
        win[i] = powf( sinf((float)i * (M_PI / (float)winSize)), 2.0f );
    
    /* Mem alloc for diffuseness of each window */
    int nDiffFrames = (int)((lenInSig + (2 * winSize)) / hopSize + 0.5f);
    diff_frames = calloc1d(nDiffFrames, sizeof(float));
    
    /* Mem alloc for a single window of processing */
    insig_win  = calloc1d(fftSize, sizeof(float));                  // for single-channel fft
    inspec_syn = calloc1d(nPV * nBins_syn, sizeof(float_complex));  // store 4-ch fft
    inspec_anl = calloc1d(nPV * nBins_anl, sizeof(float_complex));
    pvspec_win = malloc1d(nPV * nBins_anl * sizeof(float_complex));
    
    saf_rfft_create(&hFFT_syn, fftSize);
    
    /* Initialize 'prev' intensity and energy vars such that initital
     * diffuseness = 0 (appropriate in context of SRIR analysis) */
    float prev_ixyz_smooth[3] = { 1.f, 0.f, 0.f };
    float prev_energy_smooth  = 1.f;
    
    /* Window processing */
    
    int irIdx = 0; // sample frame index into pvir_pad, increments by hopsize
    int hopCount = 0;
    while (irIdx + winSize < lenInSig + 2 * winSize)
    {
        /* update progress */
        pData->progress0_1 = (float)irIdx / (float)(lenInSig + 2 * winSize);
 
        /* Transform to frequency domain, one channel at a time */
        for(int ipv = 0; ipv < nPV; ipv++){
            // window the input, placing it into a zero-padded buffer insig_win
            for(int j = 0; j < winSize; j++)
                insig_win[j] = win[j] * pvir_pad[(ipv * lenSig_pad) + irIdx + j];
            
            // full fft, put in inspec_syn
            saf_rfft_forward(hFFT_syn, insig_win, &inspec_syn[ipv * nBins_syn]);
            
            for(int j = 0, k = 0; j < nBins_anl; j++, k += fftSize/winSize)
                inspec_anl[(ipv * nBins_anl) + j] = inspec_syn[(ipv * nBins_syn) + k];
        }
        
        /* Copy spectrum of windowed input channels into pvspec_win (nPV * nBins_anl) */
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
        for(int i = 0; i < 3; i++)
            intensity[i] = crealf(pvCOV[1+i][0]);
        energy = 0.0f;
        for(int i = 0; i < 4; i++)
            energy += crealf(pvCOV[i][i])*0.5f;
        
        /* Estimating and time averaging of broadband diffuseness */
        float ixyz_smooth[3], energy_smooth, iaMag_smooth;
        iaMag_smooth = 0.0f;
        for(int i = 0; i < 3; i++) {
            ixyz_smooth[i] = ( (1.0f-smooth_const) * intensity[i] ) + ( smooth_const * prev_ixyz_smooth[i] );
            prev_ixyz_smooth[i] = ixyz_smooth[i];
            iaMag_smooth += powf( ixyz_smooth[i], 2.0f );
        }
        iaMag_smooth = sqrtf(iaMag_smooth);
        energy_smooth = ( (1.0f-smooth_const) * energy ) + ( smooth_const * prev_energy_smooth );
        prev_energy_smooth = energy_smooth;
        
        /* Store broadband diffuseness value and advance */
        diff_frames[hopCount] = 1.0f - ( iaMag_smooth / (energy_smooth + 2.23e-10f) );
        irIdx += hopSize;
        hopCount++;
    }
    
    /* Diffuse onset */
    
    /* Begin search for max diffuseness at the first window containing the
     * direct arrival, where we can expect low diffuseness */
    int directStartIdx = HOSIRR_MAX(pData->directOnsetIdx_brdbnd - hopSize, 0);
    int directHopIdx = (int)ceilf((float)directStartIdx / hopSize);
    int maxHopIdx;
    utility_simaxv(&diff_frames[directHopIdx],
                   nDiffFrames - directHopIdx,
                   &maxHopIdx);
    maxHopIdx += directHopIdx; // add the starting position back
    
    /* Index of first window above threshold */
    const float diffuseMax  = diff_frames[maxHopIdx];
    const float onsetThresh = diffuseMax * thresh_fac;
    const int   onsetHopIdx = hosirrlib_firstIndexGreaterThan(diff_frames, directHopIdx, nDiffFrames-1, onsetThresh);
    
    /* NOTE: Normally we'd assume the sample onset index is in the middle of the
     * window. I.e. onsetWinIdx * hopsize + winsize/2. But because pvsig was
     * zero-padded at the head by winsize/2, we can omit that offset here. */
    int diffuseOnsetIdx = onsetHopIdx * hopSize;
    
    if (diffuseMax < pData->diffuseMin) {
        hosirr_print_warning("Diffuseness didn't exceed the valid threshold.\nFalling back on directOnset + directOnsetFallbackDelay.");
        diffuseOnsetIdx = pData->directOnsetIdx_brdbnd + (int)(pData->diffuseOnsetFallbackDelay * fs);
    }
    /* Diffuse onset sample index.
     * This is the integer sample position in the the analyzed buffer. */
    pData->diffuseOnsetIdx = diffuseOnsetIdx;
    /* Diffuse onset time.
     * This is the acoustic property of the room. It is relative to t0, the time
     * at which the sound is emitted from the source (which can be before
     * beginning of the provided buffer, e.g. if the silence preceding the
     * direct arrival was trimmed from the recording). */
    pData->diffuseOnsetSec = ((float)diffuseOnsetIdx / fs) - pData->t0;
    
    // dbg
    printf("     begin search at hop idx: %d\n",            directHopIdx);
    printf("       diffuse onset win idx: %d\n",            onsetHopIdx);
    printf("    diffuse onset sample idx: %d (%.3f sec)\n", pData->diffuseOnsetIdx, (float)pData->diffuseOnsetIdx / pData->fs);
    printf(" t0-adjusted diff onset time: %.3f sec\n",      pData->diffuseOnsetSec);
    printf("         diffuse max hop idx: %d\n",            maxHopIdx);
    printf("           diffuse max value: %.3f\n",          diffuseMax);
    printf("        diffuse onset thresh: %.3f\n",          onsetThresh);
    
    /* Sanity checks */
    // Time events should be properly increasing
    if ((pData->t0Idx >= pData->directOnsetIdx_brdbnd) ||
        (pData->directOnsetIdx_brdbnd > pData->diffuseOnsetIdx)) {
        printf("\nt0Idx = %d \ndirectOnsetIdx_brdbnd = %d \ndiffuseOnsetIdx = %d\n",
               pData->t0Idx, pData->directOnsetIdx_brdbnd, pData->diffuseOnsetIdx);
        hosirr_print_error("! Order of analyzed events isn't valid. Should be t0Idx < directOnsetIdx < diffuseOnsetIdx.");
    }
    // Diffuse onset from t0 should be positive
    if (pData->diffuseOnsetSec < 0)
        hosirr_print_error("! Diffuse onset is negative");
    
    pData->analysisStage = thisStage;

    // Cleanup
    free(pvir);
    free(pvir_pad);
    free(win);
    free(diff_frames);
    free(insig_win);
    free(inspec_anl);
    free(inspec_syn);
    free(pvspec_win);
    saf_rfft_destroy(&hFFT_syn);
}


/* Diffuseness based on echo density:
 * Abel & Huang 2006, "A simple, robust measure of reverberation echo
 * density", In: Proc. of the 121st AES Convention, San Francisco */
void hosirrlib_setDiffuseOnsetIndex_mono(
                                         void* const hHS,
                                         float * const monoBuf,
                                         const float thresh_fac,
                                         const int winSize,         // will be made even, rounding down
                                         const int hopSize,
                                         const int startSamp,       // center of first window
                                         ANALYSIS_STAGE thisStage
                                         )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);

    checkProperProcessingOrder(pData, thisStage, __func__);

    const int nBand = pData->nBand;
    const float fs  = pData->fs;
    const int nSamp = pData->nSamp;
    const int wsize = winSize - (winSize % 2);  // make sure it's even
    const int halfWsize = wsize / 2;
    const float erfcNorm = 1.f / 0.317310507862914f; // 1 / erfc(1 / sqrt(2))
    const int nHops = (int)((nSamp - startSamp - wsize - halfWsize) / hopSize);

    if (nHops < 1) {
        hosirr_print_error("[setDiffuseOnsetIndex_mono] buffer is smaller than the window size.");
    }

    float * const window = malloc1d(wsize * sizeof(float));
    getWindowingFunction(WINDOWING_FUNCTION_HANN, wsize, window);
    float winnorm = 1.f / sumf(window, wsize);
    utility_svsmul(window, &winnorm, wsize, NULL);

    int foundOnset = 0;
    int hopIdx = 0;
    int winStart = startSamp - halfWsize;   // center the window at startSamp (OK to start at 0)
    int winOffset;
    double echoDens;
    for (int ih = 0; ih < nHops; ih++) {
        winStart = winStart + hopSize;
        
        int thisWsize, bufOffset;
        if (winStart < 0) {                 // window starts at negative index
            thisWsize = wsize + winStart;
            bufOffset = 0;
            winOffset = abs(winStart);
        } else {
            thisWsize = wsize;
            bufOffset = winStart;
            winOffset = 0;
        }

        // accumulate windowed, squared signal for standard deviation
        double sum_sd = 0.f;
        for (int iw = 0; iw < thisWsize; iw++) {
            const float sig = monoBuf[bufOffset+iw];
//            printf("sig[%d, %d]:  %.3f\n", iw, bufOffset+iw, monoBuf[bufOffset+iw]); // dbg
            sum_sd += sig * sig * window[winOffset+iw];
        }
        double sd = sqrt(sum_sd);

        // sum window scaled by tip count for echo density
        double sum_ed = 0.f;
        for (int iw = 0; iw < thisWsize; iw++) {
            if (fabsf(monoBuf[bufOffset+iw]) > sd)
                sum_ed += window[winOffset+iw];
        }
//        printf("win data\n\twinOffset: %d\n\tbufOffset: %d\n\tsum_sd: %.3f\n", winOffset, bufOffset, sum_sd); // dbg

        // normalized echo densty
        echoDens = sum_ed * erfcNorm;
        printf("[%d] winStart: %d, sampAt: %d, echoDens: %.3f\n", ih, winStart, winStart + halfWsize, echoDens); // dbg
        
        if (echoDens >= thresh_fac) {
            printf("FOUND onset\n\t hopidx: %d\n\t sum_sd: %.3f\n\t sd: %.3f \n\t sum_ed: %.3f\n\t echoDens: %.3f\n\t thresh_fac: %.3f\n\n",
                   hopIdx, sum_sd, sd, sum_ed, echoDens, thresh_fac); // dbg
            foundOnset = 1;
            break;
        }
            
        hopIdx++;
    }
    
    int diffuseOnsetIdx;
    if (foundOnset) {
        // We're calling the onset the center of the window where the measure is
        // most "sensitive". Could alternatively linearly interp between prev
        // window peak and this one.
        diffuseOnsetIdx = winStart + halfWsize;
    } else {
        hosirr_print_warning("Echo density (diffuseness) didn't exceed the threshold.\nFalling back on directOnset + directOnsetFallbackDelay.");
        diffuseOnsetIdx = pData->directOnsetIdx_brdbnd + (int)(pData->diffuseOnsetFallbackDelay * fs);
    }
    // Diffuse onset sample index in the _input buffer_
    pData->diffuseOnsetIdx = diffuseOnsetIdx;
    // Diffuse onset time in the _room_, from t0, which can be before beginning
    // of input the buffer
    pData->diffuseOnsetSec = ((float)diffuseOnsetIdx / fs) - pData->t0;

    // dbg
    printf("       diffuse onset win idx: %d\n",            hopIdx);
    printf("    diffuse onset sample idx: %d (%.3f sec)\n", pData->diffuseOnsetIdx, (float)pData->diffuseOnsetIdx / pData->fs);
    printf(" t0-adjusted diff onset time: %.3f sec\n",      pData->diffuseOnsetSec);
    printf("           diffuse max value: %.3f\n",          echoDens);
    printf("        diffuse onset thresh: %.3f\n",          thresh_fac);

    /* Sanity checks */ // TODO: Move to helper func
    //  Time events should be properly increasing
    if ((pData->t0Idx >= pData->directOnsetIdx_brdbnd) ||
        (pData->directOnsetIdx_brdbnd >= pData->diffuseOnsetIdx)) {
        printf("\nt0Idx = %d \ndirectOnsetIdx_brdbnd = %d \ndiffuseOnsetIdx = %d\n",
               pData->t0Idx, pData->directOnsetIdx_brdbnd, pData->diffuseOnsetIdx);
        hosirr_print_error("Order of analyzed events isn't valid. Should be t0Idx < directOnsetIdx < diffuseOnsetIdx.");
    }
    //  Diffuse onset from t0 should be positive
    if (pData->diffuseOnsetSec < 0)
        hosirr_print_error("Diffuse onset is negative");

    pData->analysisStage = thisStage;

    free(window);
}


/* Freq-domain band filtering via FD convolution
 * removeFiltDelay : true will truncate the head and tail of the filtered
 *                   signal, instead of (only) the tail. So either way, the
 *                   written output will always be pData->nSamp for each channel. */
void hosirrlib_splitBands(
                          void* const hHS,
                          float** const inBuf,
                          float*** const bndBuf,
                          int removeFiltDelay, // flag: 0/1
                          ANALYSIS_STAGE thisStage
                          )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    checkProperProcessingOrder(pData, thisStage, __func__);
    
    const int nInSamp   = pData->nSamp; // TODO: assumes nsamp of outbuf == nsamp of RIR
    const int filtOrder = pData->bandFiltOrder;
    
    int startCopyIdx;
    if (removeFiltDelay) {
        startCopyIdx = (int)(filtOrder / 2.f); // length of filter delay (floor)
    } else {
        startCopyIdx = 0;
    }
    
    /* Apply filterbank to rir_bands
     * Because of fftconv's indexing, it expects filter coefficients
     * for each channel. So we do one band/channel at a time.
     * The output length will be nSamp + filtOrder.
     * See fftconv(): y_len = x_len + h_len - 1; */
    float* stage = malloc1d((nInSamp + filtOrder) * sizeof(float));

    for(int ish = 0; ish < pData->nSH; ish++) {
        for(int ib = 0; ib < pData->nBand; ib++) {
            fftconv(&inBuf[ish][0],                     // input: 1 channel at a time
                    &pData->H_bandFilt[ib][0],          // band filter coeffs
                    nInSamp,                            // input length
                    filtOrder + 1,                      // filter length
                    1,                                  // 1 channel at a time
                    stage);                             // write to staging buffer

            /* Copy staging buffer to out */
            memcpy((void*)&bndBuf[ib][ish][0],
                   &stage[startCopyIdx],
                   nInSamp * sizeof(float));
        }
    }
    
    pData->analysisStage = thisStage;
    free(stage);
}

    
/* Calc bandwise RDR of the omni channel */
void hosirrlib_calcRDR(
                       void* const hHS,
                       float*** const shInBuf,          // nband x nsh x nsamp
                       float* const rdrBuf,             // nband x 1
                       const int nBand,
                       const int nSamp,
                       const int diffuseOnsetIdx,
                       int * const directOnsetIdx_bnd,  // bandwise direct onsets
                       float * const t60Buf_omni,
                       const int srcDirectivityFlag,
                       ANALYSIS_STAGE thisStage)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    checkProperProcessingOrder(pData, thisStage, __func__);
    
    // TODO: make these config constants?
    const int nCyclePerWin = 2;                 // size of the direct onset window as multiple of wavelength of band center freq
    const float minWinSizeSec = 0.0025f;        // limit minimum window size (esp. high frequencies)
    const float winAlignFac = 0.4f;             // where to align the window relative to the direct onset. 0.0 window begins at onset, 0.5 window centered around the onset
    const float srcRecDist = HOSIRR_MAX(        // source distance clipped at 25cm
                                        hosirrlib_getSrcRecDistance(hHS), 0.25);
    
    printf("Source directivity FLAG: %d\n", srcDirectivityFlag); // dbg
    
    for (int ib = 0; ib < nBand; ib++) {
        
        /* Direct energy */

        float winSizeSec = 1.f / pData->bandCenterFreqs[ib];
        winSizeSec = winSizeSec * nCyclePerWin;                     // winsize to n cycles
        winSizeSec = HOSIRR_MAX(winSizeSec, minWinSizeSec);         // clip winsize to at least minWinSizeSec
        int winSize_direct = (int)(winSizeSec * pData->fs + 0.5f);
        
        const int winSize_preOnset = (int)(winSize_direct * winAlignFac);
        
        int winStart_direct = directOnsetIdx_bnd[ib] - winSize_preOnset;
        if (winStart_direct < 0) {
            printf("! Clipping direct window [%d], windowStart: %d \n", ib, winStart_direct); // dbg
            winSize_direct = winSize_direct + winStart_direct;      // window shrinks to maintain alignment with onset
            winStart_direct = 0;
        }
        
        // Window end sample (exclusive)
        int winEnd_direct = winStart_direct + winSize_direct;
        if (winEnd_direct > nSamp)
            winEnd_direct = nSamp;
        
        double directEnergy_sum = 0.f;
        for (int is = winStart_direct; is < winEnd_direct; is++) {
            // Omni pressure from this band, scaled to bring source to 1m
            float directPress = shInBuf[ib][0][is] * srcRecDist;
            directEnergy_sum += directPress * directPress;
        }
        printf("direct onset [%d] \twinsize: %d->%d [%d %d]\n",
               ib, (int)(winSizeSec * pData->fs + 0.5f), winEnd_direct-winStart_direct, winStart_direct, winEnd_direct);
        
        /* Diffuse energy */
        
        // Diffuse energy measured up to T30 time after the diffuse onset to
        // avoid noise floor
        float t30 = t60Buf_omni[ib] * 0.5f;
        int winSize_diffuse = (int)(t30 * pData->fs);
        // Window start (inclusive) and end (exclusive) sample
        int winStart_diffuse = HOSIRR_MAX(diffuseOnsetIdx, winEnd_direct);
        int winEnd_diffuse = HOSIRR_MIN(winStart_diffuse + winSize_diffuse, nSamp);
        double diffuseEnergy_sum = 0.f;
        for (int is = winStart_diffuse; is < winEnd_diffuse; is++) {
            diffuseEnergy_sum += shInBuf[ib][0][is] * shInBuf[ib][0][is]; // omni energy from this band
        }
        
        /* RDR */
        
        // Source directivity not accounted for (assumed omni)
        // For DDR: ddr_bnd_db(ib) = rdr_bnd_db(ib) - 41;
        rdrBuf[ib] = diffuseEnergy_sum / directEnergy_sum;
        
        // printf("FDN measured late energy (TD): %.5f (%.4f dB)\n", diffuseEnergy_sum, 10.f * log10f(diffuseEnergy_sum)); // dbg
        if (srcDirectivityFlag) {
            printf("Source directivity: %.5f (%.4f dB)\n", pData->srcDirectivity[ib], 10.f * log10f(pData->srcDirectivity[ib])); // dbg
            rdrBuf[ib] = rdrBuf[ib] * pData->srcDirectivity[ib];
        }

        // dbg
        printf("   diff / dir sum [%d]: %.5f / %.5f\n", ib, diffuseEnergy_sum, directEnergy_sum);
        printf("              rdr [%d]: %.1f (%.1f dB)\n", ib, rdrBuf[ib], 10 * log10f(rdrBuf[ib]));
//        printf("  direct st/en : size [%d]: %d / %d : %d\n", ib, winStart_direct, winEnd_direct, winEnd_direct-winStart_direct);
//        printf("  direct winsize [%d]: %d (%.1f ms)\n", ib, winEnd_direct - winStart_direct, (winEnd_direct - winStart_direct)/pData->fs*1000);
//        printf("   diffuse start [%d]: %d (%.1f ms)\n", ib, winStart_diffuse, winStart_diffuse/pData->fs*1000);
//        printf(" diffuse winsize [%d]: %d (%.1f ms)\n", ib, winEnd_diffuse - winStart_diffuse, (winEnd_diffuse - winStart_diffuse)/pData->fs*1000);
    }
        
    pData->analysisStage = thisStage;
}


void hosirrlib_beamformRIR(
                           void* const hHS,
                           float*** const inBuf,
                           float*** const beamBuf,
                           ANALYSIS_STAGE thisStage)
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    checkProperProcessingOrder(pData, thisStage, __func__);
    
    const int nDir    = pData->nDir;
    const int nSamp   = pData->nSamp;
    const int nBand   = pData->nBand;
    const int nSH     = pData->nSH;
    const int shOrder = pData->shOrder;
    
    // m = 0 beam coeffs
    float * c_l = malloc1d((shOrder + 1) * sizeof(float));
    
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
     *   (ref. beamformer_process())
     *   A: float** decBeamCoeffs;      // nDir  x nSH
     *   B: float*** rirBuf_bnd_sh;     // nBand x nSH  x nSamp
     *   C: float*** rirBuf_bnd_dir;    // nBand x nDir x nSamp  */
    for (int bd = 0; bd < nBand; bd++) {
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    nDir, nSamp, nSH, 1.0f,
                    (const float*)&pData->decBeamCoeffs[0][0], nSH,         // A, 1st dim of A (row-major)
                    (const float*)&pData->rirBuf_bnd_sh[bd][0][0], nSamp,   // B, 1st dim of B
                    0.0f,                                                   // beta scalar for C
                    (float*)&beamBuf[bd][0][0], nSamp);                     // C, 1st dim of C
    }
    
    pData->analysisStage = thisStage;
    free(c_l);
}


/* Calc bandwise EDCs each directional beam */
void hosirrlib_calcEDC_beams(
                             void* const hHS,
                             float*** const inBuf,
                             float*** const edcBuf,
                             const int nBand,
                             const int nDir,
                             const int nSamp,
                             ANALYSIS_STAGE thisStage
                             )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    checkProperProcessingOrder(pData, thisStage, __func__);
    
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
                            void* const hHS,
                            float*** const shInBuf,     // nband x nsh x nsamp
                            float** const edcBuf_omn,   // nband x nsamp
                            const int nBand,
                            const int nSamp,
                            ANALYSIS_STAGE thisStage
                            )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    checkProperProcessingOrder(pData, thisStage, __func__);
    
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
 * NOTE: in-place operation, so copy data into edcBuf before calling.
 * OPTIM: vectorize */
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

/* Calc the gain offset between a source and target EDC
 * returns pressure gain in dB */
float hosirrlib_gainOffsetDB_1ch(
                                 float* const srcEDC,
                                 float* const targetEDC,
                                 const int startIdx,
                                 const int endIdx
                                 )
{
    /* Because EDCs are in dB, and assuming linear late decays,
     * the error minimization factor (gain offset) is just is just the mean
     * difference between the target and source over the chosen span, i.e.
     *    g_db = mean( edc_target_db - edc_src_db);
     *    g    = 10.^( Asd_db / 10); % to linear gain
     *    g    = sqrt( Asd); // energy -> pressure */
    int nSamp_meas = endIdx - startIdx + 1;
    float src_mn = sumf(&srcEDC[startIdx], nSamp_meas) / nSamp_meas;
    float tar_mn = sumf(&targetEDC[startIdx], nSamp_meas) / nSamp_meas;
    float gain_db = tar_mn - src_mn;
    
    // Scalar _pressure_ gain factor (linear: powf(10.f, gain_db / 20.f))
    return gain_db;
}


/* Directional gain.
 * The synthesis model assumes a fixed decay for each band for all directions,
 * so this decay constant comes form analysis of the omni response.
 * In reality decays will vary by direction, so we approximate this decay rate
 * variation with a gain variation applied to each directional decay channel
 * which have identical decay rates. The gain will be matched in the timespan
 * determine by start_db and span_db.
 */
void hosirrlib_calcDirectionalGainDB(
                                     void* const hHS,
                                     float** const dirGainBuf, // nBand x nChan
                                     const float start_db,
                                     const float span_db,
                                     const int beginIdx,
                                     const int maxGainAdjustment, // absolute value
                                     ANALYSIS_STAGE thisStage
                                     )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    checkProperProcessingOrder(pData, thisStage, __func__);
    
    const int nBand = pData->nBand;
    const int nDir = pData->nDir;
    const int nSamp = pData->nSamp;
    
    /* Here the reference for determining the gain is the omni EDC for each
     * band. So an overall energy matching will be needed.
     * An alternative would be to make the directional EDCs the target,
     * and the EDCs of the rendered FDN outputs the sources.
     */
    float**  const edcOmn_bnd = pData->edcBufOmn_bnd;  // nbnd x nsamp
    float*** const edcDir_bnd = pData->edcBuf_bnd_dir; // nbnd x ndir x nsamp
    
    /* Store start and end indices of range used for the gain calc for each
     * band, as observed in the omni (target) */
    int** const st_end_meas = (int**)malloc2d(nBand, 2, sizeof(int));
    
    for (int ib = 0; ib < nBand; ib++) {
        printf("\tband %d\n", ib); // dbg
        
        hosirrlib_findDecayBounds(&edcOmn_bnd[ib][0],
                                  beginIdx, nSamp, start_db, span_db,
                                  &st_end_meas[ib][0]);
        
        // First pass: find directional gain and accumulate for mean
        float gOffset, gOffset_mean;
        float gOffset_sum = 0.f;
        for (int id = 0; id < nDir; id++) {
            
            /* Note, omni EDC serves as a reference level, so it is passed as
             * the _source_ argument. As such, the returned gain represents the
             * gain to be applied to match the target energy distribution */
            gOffset = hosirrlib_gainOffsetDB_1ch(&edcOmn_bnd[ib][0],
                                                       &edcDir_bnd[ib][id][0],
                                                       st_end_meas[ib][0],
                                                       st_end_meas[ib][1]);
            dirGainBuf[ib][id] = gOffset;
            gOffset_sum += gOffset; // for mean
        }
        // Mean gain offset across directions
        gOffset_mean = gOffset_sum / nDir;
        printf("\tdirGain mean: %.2f dB\n", gOffset_mean); // dbg
        
        // Second pass: remove the mean, clip to +/- maxGainAdjustment, update mean
        gOffset_sum = 0.f;
        
        for (int id = 0; id < nDir; id++) {
            
            // Gain offset made to be zero-mean across directions
            gOffset = dirGainBuf[ib][id];
            printf("\t\tdir %d pre-norm dirGain: (%.2f)\t%.2f dB\n", id, powf(10, gOffset / 20.f), gOffset); // dbg
            gOffset -= gOffset_mean;
            // Clip gain offset to +/- maxGainAdjustment
            gOffset = HOSIRR_MAX(HOSIRR_MIN(gOffset, maxGainAdjustment), -maxGainAdjustment);
//            printf("\t\t       post-norm/clip dirGain: (%.2f)\t%.2f dB\n", powf(10, gOffset / 20.f), gOffset); // dbg
            
            dirGainBuf[ib][id] = gOffset;
            gOffset_sum += gOffset; // for mean
        }
        gOffset_mean = gOffset_sum / nDir;
        printf("\tpost-norm/clip dirGain mean: %.2f dB\n\n", gOffset_mean); // dbg
        
        // Second pass: remove updated mean
        for (int id = 0; id < nDir; id++) {
            printf("\t\tdir %d dirGain: \t%.2f", id, dirGainBuf[ib][id]); // dbg
            dirGainBuf[ib][id] -= gOffset_mean;   // for zero-mean normalization (preserve omni energy)
            printf(" / %.2f dB (post-clip)\n", dirGainBuf[ib][id]); // dbg
        }
        
    }
    
    free(st_end_meas);
    
    pData->analysisStage = thisStage;
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
                             float*** const edcBuf, // nBand x nChan x nSamp
                             float** const t60Buf,  // nBand x nChan
                             const int nBand,
                             const int nChan,
                             const int nSamp,
                             const float start_db,
                             const float span_db,
                             const int beginIdx,
                             ANALYSIS_STAGE thisStage
                             )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    checkProperProcessingOrder(pData, thisStage, __func__);
        
    /* find the start and end indices of the T60 measurement for each channel */
    int*** const st_end_meas = (int***)malloc3d(nBand, nChan, 2, sizeof(int));
    
    for (int ibnd = 0; ibnd < nBand; ibnd++) {
        for (int ich = 0; ich < nChan; ich++) {
            hosirrlib_findDecayBounds(&edcBuf[ibnd][ich][0],
                                      beginIdx, nSamp, start_db, span_db,
                                      &st_end_meas[ibnd][ich][0] // 2 x 1
                                      );
        }
    }
    
    /* Measure the t60 by the line of best fit */
    float* const x_slope1 = malloc1d(nSamp * sizeof(float)); // vector with slope of 1/samp
    float* const y_edc0m  = malloc1d(nSamp * sizeof(float)); // zero-mean EDC
    float* const stage    = malloc1d(nSamp * sizeof(float)); // staging buffer
    
    for (int ibnd = 0; ibnd < nBand; ibnd++) {
        for (int ich = 0; ich < nChan; ich++) {
            t60Buf[ibnd][ich] = hosirrlib_T60_lineFit(&edcBuf[ibnd][ich][0], // omni ch = 0
                                                      x_slope1, y_edc0m, stage,
                                                      st_end_meas[ibnd][ich][0], st_end_meas[ibnd][ich][1],
                                                      pData->fs);
            
            printf("t60: dir %d band %d  %.2f sec\n", ich, ibnd, t60Buf[ibnd][ich]); // dbg
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
                            void* const hHS,
                            float** const edcBuf_omn,
                            float* const t60Buf,
                            const int nBand,
                            const int nSamp,
                            const float start_db,
                            const float span_db,
                            const int beginIdx,
                            ANALYSIS_STAGE thisStage
                            )
{
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    const float fs = pData->fs;
    
    checkProperProcessingOrder(pData, thisStage, __func__);
    
    float* const x_slope1 = malloc1d(nSamp * sizeof(float)); // vector with slope of 1/samp
    float* const y_edc0m  = malloc1d(nSamp * sizeof(float)); // zero-mean EDC
    float* const stage    = malloc1d(nSamp * sizeof(float)); // staging buffer
    
    // start and end samples to measure t60 in each channel
    int** const st_end_meas = (int**)malloc2d(nBand, 2, sizeof(int));

    /* find the start and end points of the measurement */
    for (int bd = 0; bd < nBand; bd++) {
        hosirrlib_findDecayBounds(&edcBuf_omn[bd][0], // omni ch = 0
                                  beginIdx, nSamp, start_db, span_db,
                                  &st_end_meas[bd][0]);
    }
    
    for (int ibnd = 0; ibnd < nBand; ibnd++) {
        t60Buf[ibnd] = hosirrlib_T60_lineFit(&edcBuf_omn[ibnd][0], // omni ch = 0
                                             x_slope1, y_edc0m, stage,
                                             st_end_meas[ibnd][0], st_end_meas[ibnd][1],
                                             pData->fs);
        
        printf("t60 (omni): band %d  %.2f sec\n", ibnd, t60Buf[ibnd]); // dbg
        if (t60Buf[ibnd] * fs > nSamp) {
            hosirr_print_warning("Measured T60 for this band is longer than the duration of the buffer. If this is unexpected, double check that your start_db and span_db are appropriate.")
            printf("\t(Band %d, t60 %.1f sec)\n", ibnd, t60Buf[ibnd]);
        }
    }
    
    free(x_slope1);
    free(y_edc0m);
    free(stage);
    free(st_end_meas);
    pData->analysisStage = thisStage;
}


/* Measure the T60 by the line of best fit
 * x0m: zero-mean vector of a line with a slope of 1/sample
 * y0m: vector of edc values (over measurement span), with mean removed
 *      x0m      = vec_idc - mean(vec_idc);
 *      y0m      = edc_span - mean(edc_span);
 *      dc_db    = sum(x0m .* y0m) / sum(x0m.^2); // decay constant (dB/sample)
 *      t60_meas = (-60 / dc_db) / fs;
 */
float hosirrlib_T60_lineFit(
                            float* const edcBuf,
                            float* x_slopeBuf,
                            float* y_edc0mBuf,
                            float* stageBuf,
                            const int startIdx,
                            const int endIdx,
                            const float fs
                            )
{
    int nSamp_meas = endIdx - startIdx + 1;
    float y_mean = sumf(&edcBuf[startIdx], nSamp_meas) / nSamp_meas;
    
    /* Construct a vector with a slope of 1:samp with zero mean */
    float firstVal = (nSamp_meas - 1) / -2.f;
    for (int i = 0; i < nSamp_meas; i++) {
        x_slopeBuf[i] = firstVal + i;
    };
    // remove mean from EDC, within the measurement span
    utility_svssub(&edcBuf[startIdx], &y_mean, nSamp_meas, y_edc0mBuf);
    
    // covariances: x*y and x*x
    float c_xy, c_xx, dbPerSamp;
    utility_svvmul(x_slopeBuf, y_edc0mBuf, nSamp_meas, stageBuf);
    c_xy = sumf(stageBuf, nSamp_meas);
    utility_svvmul(x_slopeBuf, x_slopeBuf, nSamp_meas, stageBuf);
    c_xx = sumf(stageBuf, nSamp_meas);
    dbPerSamp = c_xy / c_xx; // slope
    
    return (-60.f / dbPerSamp) / fs;
}


void hosirrlib_findDecayBounds(
                               float* const edcBuf,     // 1 x nSamp
                               const int beginIdx,
                               const int bufLength,     // number of samples
                               const float start_db,    // <= 0
                               const float span_db,     // > 0
                               int* const st_end_meas   // output: 1 x 2
                               )
{
    /* Check start and end points of T60 measurement */
    float edcMax = edcBuf[beginIdx];
    
    int start_t60 = hosirrlib_firstIndexLessThan(edcBuf, beginIdx, bufLength - 1,
                                                 edcMax + start_db); // start_db is negative
    if (start_t60 < 0) {
        st_end_meas[0] = beginIdx;
        printf("! [findDecayBounds] No value found below the start level, returning start index as the provided beginIdx.\n");
    } else {
        st_end_meas[0] = start_t60;
    }
    
    int end_t60 = hosirrlib_firstIndexLessThan(edcBuf, start_t60, bufLength - 1,
                                               edcMax + start_db - span_db); // start_db is negative, span_db is positive
    if (end_t60 < 0) {
        // No value found below start_db - span_db, fall back to near the end of
        // the EDC.. this will likely be quite inaccurate for high frequency bands!
        st_end_meas[1] = (int)(0.7f * bufLength);
        printf("! [findDecayBounds] No value found below the decay span, returning end index 0.7 * bufLength.\n");
    } else {
        st_end_meas[1] = beginIdx + end_t60;
    }
    
    // printf("edcMax: %.4f\n", edcMax); // dbg
    // printf("start idx: %d\n", st_end_meas[bd][ch][0]); // dbg
    // printf("\tend idx: %d\n", st_end_meas[bd][ch][1]); // dbg
}

// Returns -1 on fail
int hosirrlib_firstIndexLessThan(
                                 float* vec,
                                 int startIdx,
                                 int endIdx,
                                 float thresh
                                 )
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
                                    float thresh
                                    )
{
    for (int i = startIdx; i < endIdx+1; i++) {
        if (vec[i] > thresh)
            return i;
    }
    return -1;
}

void checkProperProcessingOrder(
                          void* const hHS,
                          ANALYSIS_STAGE currentStage,
                          const char *funcName
                          )
{
    printf("\n~~ %s called. ~~\n", funcName); // dbg
    
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    if (pData->analysisStage < currentStage - 1) {
        fprintf(stderr, "ERROR [srirlib]");
        fprintf(stderr, " %s (%d)", funcName, (int)currentStage);
        fprintf(stderr, " was called before required processing stages were completed (currently at %d).\n",
                (int)pData->analysisStage);
        exit(EXIT_FAILURE);
    }
}

/* For TESTING. For the GUI to display the directional EDCs
 * might modify to view bands for debugging */
void hosirrlib_copyNormalizedEDCs_dir(
                                      void* const hHS,
                                      float** edcOut,       // ndir x nsamp
                                      float displayRange
                                      )
{
    /*
     * Note edcBuf_bnd_dir are foat*** nband x ndir x nsamp, but dirEDC pointer
     * in the UI is float** ndir x nsamp,
     * so for now, just return the first band of each direction */
    hosirrlib_data *pData = (hosirrlib_data*)(hHS);
    
    const int nBand = pData->nBand;
    const int nSamp = pData->nSamp;
    const int nDir  = pData->nDir;
    float*** const edcIn = pData->edcBuf_bnd_dir;
    
    if (pData->analysisStage >= EDC_DIR_DONE) {     // ensure EDCs are rendered
        
        /* normalise to range [-1 1] for plotting */
        float maxVal, minVal, range, add, scale, sub;

        maxVal = edcIn[0][0][0];                        // intializse to first value of first channel
        
        int maxBndIdx = 0, maxDirIdx = 0;               // dbg vars
        for (int id = 0; id < nDir; id++) {
            for (int ib = 0; ib < nBand; ib++) {
                float val0 = edcIn[ib][id][0];          // first edc value in each channel
                if (val0 > maxVal) {
                    maxVal = edcIn[ib][id][0];
                    maxBndIdx = ib; maxDirIdx = id;     // dbg vars
                }
            }
        }
        
        printf("\n>> max val found on dir %d, bnd %d: %.2f\n\n", maxDirIdx, maxBndIdx, maxVal); // dbg
        
        minVal = maxVal - displayRange;                 // just display the uper displayRange in dB
        
        // Check the first and last values of every channel
        for (int ib = 1; ib < nBand; ib++)
            for (int ich = 0; ich < nDir; ich++)
                if (edcIn[ib][ich][0] > maxVal)
                    maxVal = edcIn[ib][ich][0];

        range = maxVal - minVal;
        add   = minVal * -1.f;
        scale = 2.0f / fabsf(range);
        sub   = 1.f;
        
        // printf("max %.1f, min %.1f, rng %.1f, add %.1f, scl %.1f, sub %.1f, ",
        //       maxVal, minVal, range, add, scale, sub); // dbg
        
        printf("TEMP: viewing band channels of a single direction"); // dbg
        for(int i = 0; i < nDir; i++) {
            //int bndIdx = 0;       // for just lowest band of all directions
            //int chIdx = i;
            int bndIdx = i % nBand; // cycle through the bands of the chIdx
            int chIdx = 0;          // idx 8 for single directional decaying pw test

            utility_svsadd(&(edcIn[bndIdx][chIdx][0]),
                           &add, nSamp,
                           &edcOut[i][0]);
            utility_svsmul(&edcOut[i][0],
                           &scale, nSamp,
                           &edcOut[i][0]);
            utility_svssub(&edcOut[i][0],
                           &sub, nSamp,
                           &edcOut[i][0]);
        }
        
        /* TESTS: writing out different buffers for inspection */
        
//        // Write out band-filtered signals
//        for(int i = 0; i < pData->nDir; i++) {
//            memcpy(&edcOut[i][0],                    // copy-to channel
//                   &(pData->rirBuf_bnd_sh[i%nBand][0][0]), // bnd x ch x samp
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
        
        // copy input to output
        for(int i = 0; i < pData->nDir; i++)
            memcpy(&edcOut[i][0],             // copy-to channel
                   &(pData->rirBuf_sh[i][0]),    // [nsh][smp]
                   nSamp * sizeof(float));

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


/* For TESTING.
 * For the GUI to display the omni EDCs by band (repeated up to ndir currently
 * for debugging) */
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
    
    if (pData->analysisStage >= EDC_DIR_DONE) { // ensure EDCs are rendered
        
        /* normalise to range [-1 1] for plotting */
        float maxVal, minVal, range, add, scale, sub;

        maxVal = edcIn[0][0];                       // intializse to first value of first channel
        
        int maxBndIdx = 0; // dbg vars
        for (int ib = 0; ib < nBand; ib++) {
            float val0 = edcIn[ib][0];              // first edc value in each channel
            if (val0 > maxVal) {
                maxVal = edcIn[ib][0];
                maxBndIdx = ib; // dbg vars
            }
        }
        
        // printf("\n>> max val found on omni bnd %d: %.2f\n\n", maxBndIdx, maxVal); // dbg
        
        minVal = maxVal - displayRange; // just display the uper displayRange in dB
        // check the first and last values of every channel
        for (int ib = 1; ib < nBand; ib++)
            if (edcIn[ib][0] > maxVal)
                maxVal = edcIn[ib][0];

        range = maxVal - minVal;
        add   = minVal * -1.f;
        scale = 2.0f/fabsf(range);
        sub   = 1.f;
        
        // printf("max %.1f, min %.1f, rng %.1f, add %.1f, scl %.1f, sub %.1f, ", maxVal, minVal, range, add, scale, sub); // dbg
        
        for(int i = 0; i < nDir; i++) {
            // int bndIdx = 0;                      // for just lowest band of all directions
            int bndIdx = i % nBand;                 // cycle through the bands of the chIdx

            utility_svsadd(&(edcIn[bndIdx][0]),
                           &add, nSamp,
                           &edcOut[i][0]);
            utility_svsmul(&edcOut[i][0],
                           &scale, nSamp,
                           &edcOut[i][0]);
            utility_svssub(&edcOut[i][0],
                           &sub, nSamp,
                           &edcOut[i][0]);
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

// dbg FUNCTION
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
    if((NORMALIZATION_TYPES)newType != FUMA_NORM
       || pData->analysisOrder==ANALYSIS_ORDER_FIRST) /* FUMA only supports 1st order */
    {
        pData->norm = (NORMALIZATION_TYPES)newType;
        hosirrlib_setInputNorm(hHS, newType); // quick fix, forward the new norm type to set the input normalization
    }
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
