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

#ifndef __HOSIRRLIB_H_INCLUDED__
#define __HOSIRRLIB_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */


//#ifndef __COMMON_H_INCLUDED__ // avoid conflict with saf if included in same project

/** Available static beamforming approaches */
typedef enum {
    CARDIOID_BEAM = 1,  /**< cardioid */
    HYPERCARDIOID_BEAM, /**< hyper-cardioid */
    MAX_EV_BEAM         /**< hyper-cardioid with max_rE weighting */

} BEAM_TYPES;

/**
 * Available Ambisonic channel ordering conventions
 *
 * @note CH_FUMA only supported for 1st order input.
 */
typedef enum _CH_ORDER {
    ACN_ORDER = 1, /**< Ambisonic Channel Numbering (ACN) */
    FUMA_ORDER     /**< (Obsolete) Furse-Malham/B-format (WXYZ) */

} CH_ORDERING;

/** Number of channel ordering options */
#define HOSIRR_NUM_CH_ORDERINGS ( 2 )

/**
 * Available Ambisonic normalisation conventions
 *
 * @note NORM_FUMA only supported for 1st order input and does NOT have the
 *       1/sqrt(2) scaling on the omni.
 */
typedef enum _NORM_TYPES {
    N3D_NORM = 1,   /**< orthonormalised (N3D) */
    SN3D_NORM,      /**< Schmidt semi-normalisation (SN3D) */
    FUMA_NORM       /**< (Obsolete) Same as NORM_SN3D for 1st order */

} NORMALIZATION_TYPES;

/** Number of normalisation options */
#define HOSIRR_NUM_NORM_TYPES ( 3 )

/**
 * Available loudspeaker array presets
 */
typedef enum _LS_ARRAY_PRESETS {
    LS_ARRAY_PRESET_DEFAULT = 1,
    LS_ARRAY_PRESET_15PX,
    LS_ARRAY_PRESET_16PX,
    LS_ARRAY_PRESET_5PX,
    LS_ARRAY_PRESET_7PX,
    LS_ARRAY_PRESET_8PX,
    LS_ARRAY_PRESET_9PX,
    LS_ARRAY_PRESET_10PX,
    LS_ARRAY_PRESET_11PX,
    LS_ARRAY_PRESET_11PX_7_4,
    LS_ARRAY_PRESET_13PX,
    LS_ARRAY_PRESET_22PX,
    LS_ARRAY_PRESET_AALTO_MCC,
    LS_ARRAY_PRESET_AALTO_MCC_SUBSET,
    LS_ARRAY_PRESET_AALTO_APAJA,
    LS_ARRAY_PRESET_AALTO_LR,
    LS_ARRAY_PRESET_DTU_AVIL,
    LS_ARRAY_PRESET_ZYLIA_LAB,
    LS_ARRAY_PRESET_T_DESIGN_4,
    LS_ARRAY_PRESET_T_DESIGN_12,
    LS_ARRAY_PRESET_T_DESIGN_24,
    LS_ARRAY_PRESET_T_DESIGN_36,
    LS_ARRAY_PRESET_T_DESIGN_48,
    LS_ARRAY_PRESET_T_DESIGN_60
    
} LS_ARRAY_PRESETS;

//#endif // __COMMON_H_INCLUDED__


/**
 * Available SPHERICAL DESIGN presets
 */
typedef enum _SPHDESIGN_PRESETS {
    SPHDESIGN_PRESET_DEFAULT = 1,
    SPHDESIGN_PRESET_15PX,
    SPHDESIGN_PRESET_16PX,
    SPHDESIGN_PRESET_5PX,
    SPHDESIGN_PRESET_7PX,
    SPHDESIGN_PRESET_8PX,
    SPHDESIGN_PRESET_9PX,
    SPHDESIGN_PRESET_10PX,
    SPHDESIGN_PRESET_11PX,
    SPHDESIGN_PRESET_11PX_7_4,
    SPHDESIGN_PRESET_13PX,
    SPHDESIGN_PRESET_22PX,
    SPHDESIGN_PRESET_AALTO_MCC,
    SPHDESIGN_PRESET_AALTO_MCC_SUBSET,
    SPHDESIGN_PRESET_AALTO_APAJA,
    SPHDESIGN_PRESET_AALTO_LR,
    SPHDESIGN_PRESET_DTU_AVIL,
    SPHDESIGN_PRESET_ZYLIA_LAB,
    SPHDESIGN_PRESET_T_DESIGN_4,
    SPHDESIGN_PRESET_T_DESIGN_12,
    SPHDESIGN_PRESET_T_DESIGN_24,
    SPHDESIGN_PRESET_T_DESIGN_36,
    SPHDESIGN_PRESET_T_DESIGN_48,
    SPHDESIGN_PRESET_T_DESIGN_60
    
} SPHDESIGN_PRESETS;

/**
 * Available analysis/rendering order options
 */
#define HOSIRR_MAX_SH_ORDER ( 7 )
typedef enum _ANALYSIS_ORDERS {
    ANALYSIS_ORDER_FIRST = 1, /**< 1st-order rendering (4 channel input) */
    ANALYSIS_ORDER_SECOND,    /**< 2nd-order rendering (9 channel input) */
    ANALYSIS_ORDER_THIRD,     /**< 3rd-order rendering (16 channel input) */
    ANALYSIS_ORDER_FOURTH,    /**< 4th-order rendering (25 channel input) */
    ANALYSIS_ORDER_FIFTH,     /**< 5th-order rendering (36 channel input) */
    ANALYSIS_ORDER_SIXTH,     /**< 6th-order rendering (49 channel input) */
    ANALYSIS_ORDER_SEVENTH    /**< 7th-order rendering (64 channel input) */
    
} ANALYSIS_ORDERS;


/** SRIR analysis stages */
typedef enum {
    RIR_NOT_LOADED = 1,     /**< No RIR is set/loaded */
    RIR_LOADED,             /**< RIR is read into member buffer */
    FILTERS_INTITIALIZED,   /**< Filterbank is initialized */
    ANALYSIS_BUFS_LOADED,   /**< Staging buffers for intermediate processing are loaded */
    BANDS_SPLIT,            /**< RIR is split into bands */
    DIRECT_ONSETS_FOUND,    /**< First arrival  position is found */
    DIFFUSENESS_ONSET_FOUND,/**< T60s are calculated for each beam/band */
    EDC_OMNI_DONE,          /**< EDCs are calculated for each band of the omni channel */
    T60_OMNI_DONE,          /**< T60s are calculated for each band of the omni channel */
    RDR_DONE,               /**< RDR is calculated for each band (omni) */
    BEAMFORMED,             /**< SHD signal has been beamformed into directional signals (by band)  */
    EDC_DIR_DONE,           /**< EDCs are calculated for each beam/band */
    T60_DIR_DONE,           /**< T60s are calculated for each beam/band */
//    RIR_PROCESSED,          /**< RIR has been fully processed and ready for FDN synthesis */
//    FDN_BANDS_SPLIT,        /**< FDN is split into bands */
//    FDN_EDC_DONE            /**< EDCs are calculated for each beam/band of the FDN channels */
    DIRGAIN_DONE            /**< Directional gain has been calculated between FDN and RIR channels/bands */
} ANALYSIS_STAGE;
    
/**
 * Status of the ambisonic RIR
 */
typedef enum _AMBI_RIR_STATUS {
    AMBI_RIR_STATUS_LOADED = 0,    /**< An Ambisonic RIR has been loaded */
    AMBI_RIR_STATUS_NOT_LOADED,    /**< An Ambisonic RIR has NOT been loaded */
    AMBI_RIR_STATUS_INVALID_FORMAT /**< A file that does not have (N+1)^2
                                    *   channels was loaded, and will not be
                                    *   used for rendering */
} AMBI_RIR_STATUS;
    
/**
 * Statis of the loudspeaker RIR
 */
typedef enum _LS_RIR_STATUS {
    LS_RIR_STATUS_RENDERED = 0,        /**< Loudspeaker RIR has been rendered */
    LS_RIR_STATUS_RENDEREDING_ONGOING, /**< Loudspeaker RIR is currently being
                                        *   rendered */
    LS_RIR_STATUS_NOT_RENDERED         /**< Loudspeaker RIR has not yet been
                                        *   rendered */
} LS_RIR_STATUS;

#define HOSIRR_MAX_NUM_OUTPUTS ( 64 )
#define HOSIRR_PROGRESSTEXT_CHAR_LENGTH ( 256 )


/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of hosirrlib
 *
 * @param[in] phHS (&) address of hosirrlib handle
 */
void hosirrlib_create(void** const phHS);

/**
 * Destroys an instance of hosirrlib
 *
 * @param[in] phHS (&) address of hosirrlib handle
 */
void hosirrlib_destroy(void** const phHS);

/**
 * Analyses the input Ambisonic RIR and synthesises a Loudspeaker RIR, based
 * on the configured settings, using the HO-SIRR algorithm [1]
 *
 * @param[in] hHS hosirrlib handle
 *
 * @see [1] McCormack, L., Politis, A., Scheuregger, O., and Pulkki, V. (2019).
 *          "Higher-order processing of spatial impulse responses". In
 *          Proceedings of the 23rd International Congress on Acoustics, 9--13
 *          September 2019 in Aachen, Germany.
 */
void hosirrlib_render(void* const hHS);

    
/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/* Create */

/* Set RIR, init resources */
int hosirrlib_setRIR(
                     void* const hHS,
                     const float** H,
                     int numChannels,
                     int numSamples,
                     int sampleRate);
void hosirrlib_initBandProcessing(
                                  void* const hHS,
                                  ANALYSIS_STAGE thisStage);
void hosirrlib_allocProcBufs(
                             void * const hHS,
                             ANALYSIS_STAGE thisStage);
void hosirrlib_setUninitialized(
                                void* const hHS);
void hosirrlib_processRIR(
                          void* const hHS,
                          ANALYSIS_STAGE endStage);
void hosirrlib_splitBands(
                          void* const hHS,
                          float** const inBuf,
                          float*** const bndBuf,
                          int removeFiltDelay,
                          ANALYSIS_STAGE thisStage);
void hosirrlib_setDirectOnsetIndices(
                                     void* const hHS,
                                     float* const brdbndBuf,
                                     float*** const bndBuf,
                                     const float thresh_dB,
                                     ANALYSIS_STAGE thisStage);
void hosirrlib_setDiffuseOnsetIndex_shd(
                                        void* const hHS,
                                        float ** const rirBuf_sh, // nsh x nsamp
                                        const float thresh_fac,
                                        const int nWin_smooth,
                                        ANALYSIS_STAGE thisStage);
void hosirrlib_setDiffuseOnsetIndex_mono(
                                         void* const hHS,
                                         float * const monoBuf,
                                         const float thresh_fac,
                                         const int winSize, // even
                                         const int hopSize,
                                         const int startSamp,
                                         ANALYSIS_STAGE thisStage
                                         );
void hosirrlib_calcRDR(
                       void* const hHS,
                       float*** const shInBuf,          // nband x nsh x nsamp
                       float* const rdrBuf_omn,         // nband x 1
                       const int nBand,
                       const int nSamp,
                       const int diffuseOnsetIdx,
                       int * const directOnsetIdx_bnd,  // bandwise direct onsets
                       float *  const t60Buf_omni,
                       const int srcDirectivityFlag,       // 0: omni, 1: "regular" loudspeaker
                       ANALYSIS_STAGE thisStage);
void hosirrlib_beamformRIR(
                           void* const hHS,
                           float*** const inBuf,
                           float*** const beamBuf,
                           ANALYSIS_STAGE thisStage);
void hosirrlib_calcEDC_beams(
                             void* const hHS,
                             float*** const beamBuf,
                             float*** const edcBuf,
                             const int nBand,
                             const int nDir,
                             const int nSamp,
                             ANALYSIS_STAGE thisStage);
void hosirrlib_calcEDC_omni(
                            void* const hHS,
                            float*** const shInBuf,
                            float** const edcBuf_omn,
                            const int nBand,
                            const int nSamp,
                            ANALYSIS_STAGE thisStage);
void hosirrlib_calcT60_beams(
                             void* const hHS,
                             float*** const edcBuf,
                             float** const t60Buf,
                             const int nBand,
                             const int nDir,
                             const int nSamp,
                             const float startDb,
                             const float spanDb,
                             const int beginIdx,
                             ANALYSIS_STAGE thisStage);
void hosirrlib_calcT60_omni(
                            void* const hHS,
                            float** const edcBuf_omn,
                            float* const t60Buf,
                            const int nBand,
                            const int nSamp,
                            const float startDb,
                            const float spanDb,
                            const int beginIdx,
                            ANALYSIS_STAGE thisStage);
void hosirrlib_calcDirectionalGainDB(
                                   void* const hHS,
                                   float** const dirGainBuf, // nBand x nChan
                                   const float start_db,
                                   const float span_db,
                                   const int beginIdx,
                                   ANALYSIS_STAGE thisStage);
/* Helpers */

int getDirectOnset_1ch(
                       const float * const chan,
                       float * const tmp, // buffer to hold copied channel data
                       const float thresh_dB,
                       const int nSamp);
void hosirrlib_calcEDC_1ch(
                           float* const dataBuf,
                           const int nSamp);
float hosirrlib_gainOffsetDB_1ch(
                              float* const srcEDC,
                              float* const targetEDC,
                              const int startIdx,
                               const int endIdx);
void hosirrlib_findDecayBounds(float* const edcBuf,
                               const int beginIdx,
                               const int bufLength,
                               const float start_db,
                               const float span_db,
                               int* const st_end_meas);
float hosirrlib_T60_lineFit(float* const edcBuf,
                            float* x_slopeBuf,
                            float* y_edc0mBuf,
                            float* stageBuf,
                            const int startIdx,
                            const int endIdx,
                            const float fs);
int hosirrlib_firstIndexLessThan(
                                 float* vec,
                                 int startIdx,
                                 int endIdx,
                                 float thresh);
int hosirrlib_firstIndexGreaterThan(
                                    float* vec,
                                    int startIdx,
                                    int endIdx,
                                    float thresh);
void hosirrlib_setSrcPosition(
                              void* const hHS,
                              const float x,
                              const float y,
                              const float z
                              );
void hosirrlib_setRecPosition(
                              void* const hHS,
                              const float x,
                              const float y,
                              const float z);
float hosirrlib_getSrcRecDistance(
                                  void* const hHS);
void hosirrlib_renderTMP(
                         void* const hHS);
void hosirrlib_copyNormalizedEDCs_dir(
                                      void* const hHS,
                                      float** edcCopy,
                                      float displayRange);
void hosirrlib_copyNormalizedEDCs_omni(
                                       void* const hHS,
                                       float** edcCopy,
                                       float displayRange);
void checkProperProcessingOrder(
                          void* const hHS,
                          ANALYSIS_STAGE currentStage,
                          const char *funcName
                          );

/* Getters */

int hosirrlib_getNumDirections(
                               void* const hHS);

/* Debug  func */
// void hosirrlib_inspectFilts(void* const hHS);

/**
 * Sets a flag, as to whether the renderer should isolate the first peak in the
 * Ambisonic RIR and process it based on broad-band analysis (0:disabled,
 * 1:enabled)
 *
 * This can help reduce timbral colouration in some cases.
 */
void hosirrlib_setBroadBandFirstPeakFLAG(void* const hHS, int newState);

/**
 * Sets the windowing length, in samples, used by the HOSIRR method
 */
void hosirrlib_setWindowLength(void* const hHS, int newValue);

/**
 * Sets the wet/dry balance; when 0: only dry (non-diffuse), 1: only wet
 * (diffuse)
 */
void hosirrlib_setWetDryBalance(void* const hHS, float newValue); 
    
/**
 * Load input Ambisonic (spherical harmonic) room impulse response (RIR) to be
 * rendered by hosirrlib.
 *
 * @param[in] hHS         hosirrlib handle
 * @param[in] H           The Ambisonic RIR; numChannels x numSamples
 * @param[in] numChannels Number of channels in H
 * @param[in] numSamples  Number of samples per channel in H
 * @param[in] sampleRate  Sample rate of the loaded H
 */
int hosirrlib_setAmbiRIR(void* const hHS,
                         const float** H,
                         int numChannels,
                         int numSamples,
                         int sampleRate);

/**
 * Sets the analysis/rendering order (see ANALYSIS_ORDERS enum)
 */
void hosirrlib_setAnalysisOrder(void* const hHS,  int newValue);
    
/**
 * Sets the azimuth of a specific loudspeaker
 */
void hosirrlib_setLoudspeakerAzi_deg(void* const hHS, int index, float newAzi_deg);

/**
 * Sets the elevation of a specific loudspeaker
 */
void hosirrlib_setLoudspeakerElev_deg(void* const hHS, int index, float newElev_deg);

/**
 * Sets the number of loudspeakers in the setup
 */
void hosirrlib_setNumLoudspeakers(void* const hHS, int new_nLoudspeakers);
 
/**
 * Sets the output loudspeaker preset
 *
 * For convenience, presets for several popular arrangements are included (see
 * LS_ARRAY_PRESETS enum).
 */
void hosirrlib_setOutputConfigPreset(void* const hHS, int newPresetID);

/**
 * Sets the Ambisonic channel ordering convention to use, which should match the
 * convention employed by the input RIR (see 'CH_ORDER' enum)
 */
void hosirrlib_setChOrder(void* const hHS, int newOrder);

/**
 * Sets the Ambisonic normalisation convention to use, which shoudl match the
 * convention employed by the input RIR (see 'NORM_TYPE' enum)
 */
void hosirrlib_setNormType(void* const hHS, int newType);

    
/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the order of the loaded Ambisonic RIR
 */
int hosirrlib_getAmbiRIRinputOrder(void* const hHS);

/**
 * Returns the length, in samples, of the loaded Ambisonic RIR
 */
int hosirrlib_getAmbiRIRlength_samples(void* const hHS);

/**
 * Returns the length, in seconds, of the loaded Ambisonic RIR
 */
float hosirrlib_getAmbiRIRlength_seconds(void* const hHS);

/**
 * Returns the sampling rate of the loaded Ambisonic RIR
 */
int hosirrlib_getAmbiRIRsampleRate(void* const hHS);
    
/**
 * Returns a flag, dictating whether the renderer should isolate the first peak
 * in the Ambisonic RIR and process it based on broad-band analysis (0:disabled,
 * 1:enabled)
 */
int hosirrlib_getBroadBandFirstPeakFLAG(void* const hHS);

/**
 * Returns the windowing length, in samples, used by the HOSIRR method
 */
int hosirrlib_getWindowLength (void* const hHS);

/**
 * Returns the wet/dry balance; when 0: only dry (non-diffuse), 1: only wet
 * (diffuse)
 */
float hosirrlib_getWetDryBalance(void* const hHS);

/**
 * Returns the status of the Ambisonic RIR (see AMBI_RIR_STATUS enum)
 */
int hosirrlib_getAmbiRIRstatus(void* const hHS);

/**
 * Returns the status of the Loudspeaker RIR (see LS_RIR_STATUS enum)
 */
int hosirrlib_getLsRIRstatus(void* const hHS);

/**
 * Returns the current rendering progress (0: just started, 0<x<1: ongoing,
 * 1: finished)
 */
float hosirrlib_getProgress0_1(void* const hHS);

/**
 * Returns a string of length 'HOSIRR_PROGRESSTEXT_CHAR_LENGTH', which describes
 * the progress of the current render
 */
void hosirrlib_getProgressText(void* const hHS, char* text);
    
void hosirrlib_getLsRIR(void* const hHS, float** lsRIR);
    
/**
 * Returns the master/maximum decoding order (see 'ANALYSIS_ORDERS' enum)
 */
int hosirrlib_getAnalysisOrder(void* const hHS);
    
/**
 * Returns the loudspeaker azimuth for a given index
 */
float hosirrlib_getLoudspeakerAzi_deg(void* const hHS, int index);

/**
 * Returns the loudspeaker elevation for a given index
 */
float hosirrlib_getLoudspeakerElev_deg(void* const hHS, int index);

/**
 * Returns the number of loudspeakers in the current layout
 */
int hosirrlib_getNumLoudspeakers(void* const hHS);

/**
 * Returns the maximum number of loudspeakers supported by hosirrlib
 */
int hosirrlib_getMaxNumLoudspeakers(void);
 
/**
 * Returns the Ambisonic channel ordering convention currently being, which
 * should match the convention employed by the input RIR (see 'CH_ORDER' enum)
 */
int hosirrlib_getChOrder(void* const hHS);

/**
 * Returns the Ambisonic normalisation convention currently being used, which
 * should match the convention employed by the input RIR (see 'NORM_TYPE' enum)
 */
int hosirrlib_getNormType(void* const hHS);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __HOSIRRLIB_H_INCLUDED__ */
