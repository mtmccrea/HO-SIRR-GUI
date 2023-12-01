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

void loadSphDesignPreset(SPHDESIGN_PRESETS preset,
    float dirs_deg[MAX_NUM_LOUDSPEAKERS_IN_PRESET]
                  [2], // TODO: SPHDESIGNs coincide with LS arrays for now
    int* newNCH)
{
    int ch, i, nCH;

    switch (preset) {
    default:
    case SPHDESIGN_PRESET_DEFAULT:
    case SPHDESIGN_PRESET_15PX:
        nCH = 15;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __SphCovering_15_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_16PX:
        nCH = 16;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __SphCovering_16_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_5PX:
        nCH = 5;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __5pX_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_7PX:
        nCH = 7;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __7pX_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_8PX:
        nCH = 8;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __8pX_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_9PX:
        nCH = 9;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __9pX_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_10PX:
        nCH = 10;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __10pX_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_11PX:
        nCH = 11;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __11pX_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_11PX_7_4:
        nCH = 11;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __11pX_7_4_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_13PX:
        nCH = 13;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __13pX_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_22PX:
        nCH = 22;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __22pX_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_AALTO_MCC:
        nCH = 45;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __Aalto_MCC_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_AALTO_MCC_SUBSET:
        nCH = 37;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __Aalto_MCCsubset_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_AALTO_APAJA:
        nCH = 29;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __Aalto_Apaja_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_AALTO_LR:
        nCH = 13;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __Aalto_LR_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_DTU_AVIL:
        nCH = 64;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __DTU_AVIL_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_ZYLIA_LAB:
        nCH = 22;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __Zylia_Lab_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_T_DESIGN_4:
        nCH = 4;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __Tdesign_degree_2_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_T_DESIGN_12:
        nCH = 12;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __Tdesign_degree_4_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_T_DESIGN_24:
        nCH = 24;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __Tdesign_degree_6_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_T_DESIGN_36:
        nCH = 36;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __Tdesign_degree_8_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_T_DESIGN_48:
        nCH = 48;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __Tdesign_degree_9_dirs_deg[ch][i];
        break;
    case SPHDESIGN_PRESET_T_DESIGN_60:
        nCH = 60;
        for (ch = 0; ch < nCH; ch++)
            for (i = 0; i < 2; i++)
                dirs_deg[ch][i] = __Tdesign_degree_10_dirs_deg[ch][i];
        break;
    }

    /* Fill remaining slots with default coords */
    for (; ch < MAX_NUM_LOUDSPEAKERS_IN_PRESET;
         ch++) // TODO: SPHDESIGNs coincide with LS arrays for now
        for (i = 0; i < 2; i++)
            dirs_deg[ch][i] = __default_LScoords64_rad[ch][i] * (180.0f / M_PI);

    /* specify new number of channels (for dynamically changing the number of TFT channels) */
    (*newNCH) = nCH;
}
