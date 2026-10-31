;Written by Athanasios (Thanos) Klidaras from 5/16/2024 to summer 2026 with input from graduate advisor Briony Horgan and contributions from Chandler Albright.Built with the assistance of ChatGPT and Gemini LLMs.
;
;Base code = heavily modified and improved from MASTER_paper1.pro (see Klidaras et al. 2025)

;PURPOSE: Runs through the following steps for a CRISM observation:
;         (1) Generate clay band center map, clay band depth map, clay detection mask
;         (2) Classify by continuum-removed band center into different alteration phases to generate clay phase map
;         (3) Generate parameters to quantify Al-unit diversity and isolate kaolin signatures (BD2165/BD2205 & IBD Asymmetry Ratio)
;         (4) Generate ferrous slope parameter (1200-1800 nm slope)
;         (5) Automatically generate ROIs by isolating the largest valid contiguous cluster that qualifies for various thresholds of spectral parameters, and extracts and saves a spectrum for each group in clay phase map
;         (6) Save Viviano et al. refined spectral summary parameters
;
;REQUIRED FUNCTIONS FOR EACH STEP:       
;         (1) Klidaras_band_center_depth (requires Klidaras_destripe_mask_filter)
;         (2) Klidaras_classify_alteration_phases
;         (3) Klidaras_Kaolinite_parameter (requires calc_band_depth, calc_integrated_band_depth, calc_IBD_asymmetry_ratio, classify_kaolin)
;         (4) calc_spectral_slope
;         (5) initial_ROI_pass, ROI_valid_largest_cluster (requires save_ROI_spectrum) 
;         (6) save_spectral_parameters
;
;TO RUN:
;If you want to run on a single CRISM cube, enter this into your console:
;
;         ENVI>test=master_paper2(if_1,sr_1,'1') 
;
;         where if_1 = TRDR .if cube, sr_1 = TRDR .sr cube, and '1' is the ID of that observation which will control the filenames of saved data
;         output variable 'test' will yield nothing (intended)
;         
;If you want to run on the entire spectral database, enter this into your console:
;         
;         ENVI>generate_database()
;         (required functions: generate_database, crism_load_data)
;        
;UPDATES: 
;         8/19/2025 iteration3 of spectral database adds: .csv of unmasked phyllosilicate array, band depth at band center array, the ROI arrays, and BD2290 & BD2210_2, big bug fix: changed refined spectral summary parameters to their correct array
;         9/22/2025 Modified master_paper2 function to output ROIs so I can experiment with identifying contiguous clusters
;         10/01/2025 Added iteration4 features (improved contiguous largest cluster ROIs, fixed threshold for masking phyllosilicate phase map, more spectral parameters)              
;         11/27/2025 Re-added iteration5 features (Contiguous ROIs & avg spectra based on spectral parameter + threshold, yet more spectral parameters) 
;         12/10/2025 Updated ROI logic so that invalid ROIs are skipped (files & spectra not written), moved plotting code into ROI_valid_largest_cluster to reduce clutter
;                    Added BD2500_2 (>0.5%) CINDEX2 (>0), SINDEX2 (>0), IBD Asymmetry Ratio (>2) ROIs 
;         12/12/2025 Completely rewrote continuum-removed band center & band depth at band center code to fix errors, make more efficient
;                    New function = Klidaras_band_center_depth, now uses shoulders 2.12-2.24 & 2.24-2.36 and merges based on which of the two has higher band depth
;                    Threshold for high band depth mask increased to >1%
;                    New function = Klidaras_classify_alteration_phases, classifies into 6 groups, shorter & more efficient
;         12/23/2025 Added new function (Klidaras_destripe_mask_filter) which filters out noisy columns from the high BD mask before it is outputted
;         1/7/2026 Added kaolin param ROIs, renamed GtMeanBD to HighBandDepth ROIs (name misleading, since threshold actually Gt mean + 1*std)
;         1/12/2026 Added high ferrous slope ROIs (threshold=gt 4e-5), tweaked CINDEX2 and SINDEX2 ROI thresholds, added BD1750_2/BD2100_2/BD2250/BD2265 ROIs
;         6/2/2026 Improved commenting of code to be more descriptive
;         6/8/2026 iteration8 of spectral database: added BD530_2, fixed bug with ferrous spectral slope. Simplified terminology ("merged band center wvl" -> clay band center map, "merged band depth" -> clay band depth map, "high band depth mask" -> clay detection mask, "masked alteration phase map" -> clay phase map)

FUNCTION master_paper2,if_qub,sr_qub,id

print,'Processing CRISM image with the CCS database ID: ',id
restore,'crism_wvl_array.sav'
bad_pixels=where(if_qub[*,*,0] gt 100.)

;Step (1): Generate clay band center map, clay band depth map, clay detection mask

output=Klidaras_band_center_depth(if_qub, crism_wvl_array, id)
clay_band_center_map=output.merged_center_wvl
clay_band_depth_map=output.merged_depth_map
clay_detection_mask=output.merged_high_BD_mask

;Step (2): Classify by continuum-removed band center into different alteration phases to generate clay phase map

clay_phase_map=Klidaras_classify_alteration_phases(clay_band_center_map, clay_detection_mask, id)

;Step(3): Generate parameters to quantify Al-unit diversity and isolate kaolin signatures (BD2165/BD2205 & IBD Asymmetry Ratio)
masked_kaolin_param=Klidaras_kaolinite_parameter(if_qub,crism_wvl_array,id,clay_phase_map)

;Step(4): Generate ferrous slope parameter (1200-1800 nm slope)
ferrous_slope=calc_spectral_slope(if_qub,1200,1800,7,7,'Ferrous_slope_Klidaras',id,wavelengths=crism_wvl_array, average_mode='median')

;Step(5): Automatically generate ROIs by isolating the largest valid contiguous cluster that qualifies for various thresholds of spectral parameters, and extracts and saves a spectrum for each group in clay phase map

;FERROUS SLOPE ROIs, threshold gt 3e-5
initial_ROIs=initial_ROI_pass(clay_phase_map,sr_qub,id,4e-5,clay_band_depth_map,masked_kaolin_param,ferrous_slope,ROI_type='FerrousSlope')
cluster_ROIs=ROI_valid_largest_cluster(initial_ROIs,id,if_qub,crism_wvl_array,[700,2500,0.05,0.3],ROI_type='FerrousSlope')

;KAOLINITE ROIs, threshold eq 1 (pixels meeting kaolinite criteria have value of 1)
initial_ROIs=initial_ROI_pass(clay_phase_map,sr_qub,id,1,clay_band_depth_map,masked_kaolin_param,ferrous_slope,ROI_type='Kaolin')
cluster_ROIs=ROI_valid_largest_cluster(initial_ROIs,id,if_qub,crism_wvl_array,[2100,2300,0.1,0.2],ROI_type='Kaolin')

;GREATER THAN MEAN + 1 STD BAND DEPTH ROIS
initial_ROIs=initial_ROI_pass(clay_phase_map,sr_qub,id,0,clay_band_depth_map,masked_kaolin_param,ferrous_slope,ROI_type='HighBandDepth')
cluster_ROIs=ROI_valid_largest_cluster(initial_ROIs,id,if_qub,crism_wvl_array,[2000,2500,0.1,0.2],ROI_type='HighBandDepth')

;HIGH BD860_2 ROIS, threshold = >0.01 (1%)
initial_ROIs=initial_ROI_pass(clay_phase_map,sr_qub,id,0.01,clay_band_depth_map,masked_kaolin_param,ferrous_slope,ROI_type='BD860_2')
cluster_ROIs=ROI_valid_largest_cluster(initial_ROIs,id,if_qub,crism_wvl_array,[500,1200,0.,0.3],ROI_type='BD860_2')

;HIGH BD2500_2 ROIS, threshold = >0.01 (1%)
initial_ROIs=initial_ROI_pass(clay_phase_map,sr_qub,id,0.01,clay_band_depth_map,masked_kaolin_param,ferrous_slope,ROI_type='BD2500_2')
cluster_ROIs=ROI_valid_largest_cluster(initial_ROIs,id,if_qub,crism_wvl_array,[2000,2600,0.1,0.2],ROI_type='BD2500_2')

;HIGH CINDEX2 ROIS, threshold = >0.15
initial_ROIs=initial_ROI_pass(clay_phase_map,sr_qub,id,0.15,clay_band_depth_map,masked_kaolin_param,ferrous_slope,ROI_type='CINDEX2')
cluster_ROIs=ROI_valid_largest_cluster(initial_ROIs,id,if_qub,crism_wvl_array,[2000,2600,0.1,0.2],ROI_type='CINDEX2')

;HIGH SINDEX2 ROIS, threshold = >0.05
initial_ROIs=initial_ROI_pass(clay_phase_map,sr_qub,id,0.05,clay_band_depth_map,masked_kaolin_param,ferrous_slope,ROI_type='SINDEX2')
cluster_ROIs=ROI_valid_largest_cluster(initial_ROIs,id,if_qub,crism_wvl_array,[1200,2400,0.,0.2],ROI_type='SINDEX2')

;HIGH BD1750_2 ROIS, threshold = >0.01 (1%)
initial_ROIs=initial_ROI_pass(clay_phase_map,sr_qub,id,0.01,clay_band_depth_map,masked_kaolin_param,ferrous_slope,ROI_type='BD1750_2')
cluster_ROIs=ROI_valid_largest_cluster(initial_ROIs,id,if_qub,crism_wvl_array,[1200,2400,0.,0.3],ROI_type='BD1750_2')

;HIGH BD2100_2 ROIS, threshold = >0.01 (1%)
initial_ROIs=initial_ROI_pass(clay_phase_map,sr_qub,id,0.01,clay_band_depth_map,masked_kaolin_param,ferrous_slope,ROI_type='BD2100_2')
cluster_ROIs=ROI_valid_largest_cluster(initial_ROIs,id,if_qub,crism_wvl_array,[1800,2400,0.,0.3],ROI_type='BD2100_2')

;HIGH BD2250 ROIS, threshold = >0.01 (1%)
initial_ROIs=initial_ROI_pass(clay_phase_map,sr_qub,id,0.01,clay_band_depth_map,masked_kaolin_param,ferrous_slope,ROI_type='BD2250')
cluster_ROIs=ROI_valid_largest_cluster(initial_ROIs,id,if_qub,crism_wvl_array,[1800,2400,0.,0.3],ROI_type='BD2250')

;HIGH BD2265 ROIS, threshold = >0.01 (1%)
initial_ROIs=initial_ROI_pass(clay_phase_map,sr_qub,id,0.01,clay_band_depth_map,masked_kaolin_param,ferrous_slope,ROI_type='BD2265')
cluster_ROIs=ROI_valid_largest_cluster(initial_ROIs,id,if_qub,crism_wvl_array,[1800,2400,0.,0.3],ROI_type='BD2265')

;Step(6): Save Viviano et al. refined spectral summary parameters
output=save_spectral_parameters(sr_qub,id,if_qub)

return,0

END

;+
; NAME:
;   KLIDARAS_BAND_CENTER_DEPTH
;
; PURPOSE:
;   Calculates spectral parameters for Al/OH (2120-2240 nm) and Fe/Mg (2240-2360 nm).
;   Merges them based on whichever feature is deeper.
;   Creates a SEPARATE mask for high band depth pixels (>1%) without modifying the data.
;
; RETURN:
;   result.MERGED_CENTER_WVL : Composite center (winner of Al vs FeMg) - ALL PIXELS
;   result.MERGED_DEPTH_MAP  : Composite depth (winner of Al vs FeMg) - ALL PIXELS
;   result.HIGH_BD_MASK      : Binary mask (1=High Depth, 0=Low Depth)
;-
FUNCTION Klidaras_band_center_depth, cube, wvl, id

  COMPILE_OPT IDL2

  ; --- 1. SETUP & ERROR CHECKING ---
  sz = SIZE(cube, /STRUCTURE)
  IF sz.n_dimensions NE 3 THEN BEGIN
    MESSAGE, 'Input must be a 3D spectral cube.', /CONTINUE
    RETURN, -1
  ENDIF

  nx = sz.dimensions[0]
  ny = sz.dimensions[1]
  nb = sz.dimensions[2]

  IF N_ELEMENTS(wvl) NE nb THEN BEGIN
    MESSAGE, 'Wavelength array size does not match cube bands.', /CONTINUE
    RETURN, -1
  ENDIF

  width = 3

  ; Initialize Output Arrays
  Al_center_wvl   = MAKE_ARRAY(nx, ny, /FLOAT, VALUE=0.0)
  Al_depth_map    = MAKE_ARRAY(nx, ny, /FLOAT, VALUE=0.0)
  FeMg_center_wvl = MAKE_ARRAY(nx, ny, /FLOAT, VALUE=0.0)
  FeMg_depth_map  = MAKE_ARRAY(nx, ny, /FLOAT, VALUE=0.0)

  ; ====================================================================
  ; BLOCK 1: AL/OH FEATURE (2120 nm - 2240 nm)
  ; ====================================================================
  PRINT, ' '
  PRINT, 'Processing Feature 1: Al/OH'

  target_L = 2120.0
  target_R = 2240.0
  diff_L = MIN(ABS(wvl - target_L), idx_left)
  diff_R = MIN(ABS(wvl - target_R), idx_right)

  s_L = (idx_left - (width / 2)) > 0
  e_L = (s_L + width - 1) < (nb - 1)
  s_R = (idx_right - (width / 2)) > 0
  e_R = (s_R + width - 1) < (nb - 1)

  x_L = FLOAT(wvl[idx_left])
  x_R = FLOAT(wvl[idx_right])
  val_L = TOTAL(cube[*,*,s_L:e_L], 3) / FLOAT(e_L - s_L + 1)
  val_R = TOTAL(cube[*,*,s_R:e_R], 3) / FLOAT(e_R - s_R + 1)

  dist_total = x_R - x_L
  IF dist_total EQ 0 THEN dist_total = 1.0
  slope_map = (val_R - val_L) / dist_total

  min_val_tracker = MAKE_ARRAY(nx, ny, /FLOAT, VALUE=999.0)

  FOR i = idx_left, idx_right DO BEGIN
    cur_wvl   = FLOAT(wvl[i])
    dist_curr = cur_wvl - x_L
    cont_i    = val_L + (slope_map * dist_curr)
    spec_i    = cube[*,*,i]
    ratio_i   = MAKE_ARRAY(nx, ny, /FLOAT, VALUE=1.0)

    good_mask = WHERE((cont_i NE 0.0) AND FINITE(cont_i), count)
    IF count GT 0 THEN ratio_i[good_mask] = spec_i[good_mask] / cont_i[good_mask]

    update_idx = WHERE(ratio_i LT min_val_tracker, u_count)
    IF u_count GT 0 THEN BEGIN
      min_val_tracker[update_idx] = ratio_i[update_idx]
      Al_center_wvl[update_idx]   = cur_wvl
    ENDIF
  ENDFOR

  Al_depth_map = 1.0 - min_val_tracker
  bad_mask = WHERE(min_val_tracker EQ 999.0, b_count)
  IF b_count GT 0 THEN BEGIN
    Al_depth_map[bad_mask]  = 0.0
    Al_center_wvl[bad_mask] = 0.0
  ENDIF

  ; ====================================================================
  ; BLOCK 2: Fe/Mg FEATURE (2240 nm - 2360 nm)
  ; ====================================================================
  PRINT, 'Processing Feature 2: Fe/Mg'

  target_L = 2240.0
  target_R = 2360.0
  diff_L = MIN(ABS(wvl - target_L), idx_left)
  diff_R = MIN(ABS(wvl - target_R), idx_right)

  s_L = (idx_left - (width / 2)) > 0
  e_L = (s_L + width - 1) < (nb - 1)
  s_R = (idx_right - (width / 2)) > 0
  e_R = (s_R + width - 1) < (nb - 1)

  x_L = FLOAT(wvl[idx_left])
  x_R = FLOAT(wvl[idx_right])
  val_L = TOTAL(cube[*,*,s_L:e_L], 3) / FLOAT(e_L - s_L + 1)
  val_R = TOTAL(cube[*,*,s_R:e_R], 3) / FLOAT(e_R - s_R + 1)

  dist_total = x_R - x_L
  IF dist_total EQ 0 THEN dist_total = 1.0
  slope_map = (val_R - val_L) / dist_total

  min_val_tracker = MAKE_ARRAY(nx, ny, /FLOAT, VALUE=999.0)

  FOR i = idx_left, idx_right DO BEGIN
    cur_wvl   = FLOAT(wvl[i])
    dist_curr = cur_wvl - x_L
    cont_i    = val_L + (slope_map * dist_curr)
    spec_i    = cube[*,*,i]
    ratio_i   = MAKE_ARRAY(nx, ny, /FLOAT, VALUE=1.0)

    good_mask = WHERE((cont_i NE 0.0) AND FINITE(cont_i), count)
    IF count GT 0 THEN ratio_i[good_mask] = spec_i[good_mask] / cont_i[good_mask]

    update_idx = WHERE(ratio_i LT min_val_tracker, u_count)
    IF u_count GT 0 THEN BEGIN
      min_val_tracker[update_idx] = ratio_i[update_idx]
      FeMg_center_wvl[update_idx] = cur_wvl
    ENDIF
  ENDFOR

  FeMg_depth_map = 1.0 - min_val_tracker
  bad_mask = WHERE(min_val_tracker EQ 999.0, b_count)
  IF b_count GT 0 THEN BEGIN
    FeMg_depth_map[bad_mask]  = 0.0
    FeMg_center_wvl[bad_mask] = 0.0
  ENDIF

  ; ====================================================================
  ; BLOCK 3: MERGING & MASKING (CORRECTED)
  ; ====================================================================

  ; Initialize merged arrays
  merged_band_center_wvl = MAKE_ARRAY(nx, ny, /FLOAT, VALUE=0.0)
  merged_depth_map       = MAKE_ARRAY(nx, ny, /FLOAT, VALUE=0.0)

  ; --- Step A: Compare Al vs FeMg ---
  ; 1. Where Al Depth >= FeMg Depth, Al wins
  idx_Al = WHERE(Al_depth_map GE FeMg_depth_map, count_Al)
  IF count_Al GT 0 THEN BEGIN
    merged_band_center_wvl[idx_Al] = Al_center_wvl[idx_Al]
    merged_depth_map[idx_Al]       = Al_depth_map[idx_Al]
  ENDIF

  ; 2. Where FeMg Depth > Al Depth, FeMg wins
  idx_Fe = WHERE(FeMg_depth_map GT Al_depth_map, count_Fe)
  IF count_Fe GT 0 THEN BEGIN
    merged_band_center_wvl[idx_Fe] = FeMg_center_wvl[idx_Fe]
    merged_depth_map[idx_Fe]       = FeMg_depth_map[idx_Fe]
  ENDIF

  ;Multiply merged_depth_map by 100 to convert from decimal to into percentage(%) band depth
  merged_depth_map=merged_depth_map*100

  ; --- Step B: Create Threshold Mask Only ---
  ; We DO NOT apply this to the data arrays. We just export it.
  merged_high_BD_mask = MAKE_ARRAY(nx, ny, /BYTE, VALUE=0)

  idx_good = WHERE(merged_depth_map GE 1, count_good)
  IF count_good GT 0 THEN merged_high_BD_mask[idx_good] = 1

  ;Filtering the high band depth mask to surpress noisy column stripe artifacts (generate clay detection mask)
  clay_detection_mask=Klidaras_destripe_mask_filter(merged_high_BD_mask)

  ;Outputting relevant files as .csvs
  filename='Saved data paper2/'+id+'_clay_band_center_map.csv'
  write_csv,filename,merged_band_center_wvl
  print,'File written: '+filename
  filename='Saved data paper2/'+id+'_clay_band_depth_map.csv'
  write_csv,filename,merged_depth_map
  print,'File written: '+filename
  filename='Saved data paper2/'+id+'_clay_detection_mask.csv'
  write_csv,filename,clay_detection_mask
  print,'File written: '+filename


  PRINT, 'Mask Created. Found ', count_good, ' pixels with Depth >= 1%.'

  ; ====================================================================
  ; FINAL OUTPUT
  ; ====================================================================

  result = { Al_center_wvl:      Al_center_wvl,      $
    Al_depth_map:       Al_depth_map,       $
    FeMg_center_wvl:    FeMg_center_wvl,    $
    FeMg_depth_map:     FeMg_depth_map,     $
    merged_center_wvl:  merged_band_center_wvl, $
    merged_depth_map:   merged_depth_map,   $
    merged_high_BD_mask:       merged_high_BD_mask }

  RETURN, result

END

;+
; NAME:
;   Klidaras_destripe_mask_filter
;
; PURPOSE:
;   Removes vertical streaks from a mask image using a two-step filtering process:
;   1. Morphological opening with a horizontal kernel (width=5).
;   2. Removal of all contiguous clusters smaller than 15 pixels.
;
; CALLING SEQUENCE:
;   Result = Klidaras_destripe_mask_filter(mask_image)
;
; INPUTS:
;   mask_image: A 2D array representing the mask (typically 0s and 1s).
;
; OUTPUTS:
;   Returns the filtered 2D mask array.
;-

FUNCTION Klidaras_destripe_mask_filter, mask_image
  COMPILE_OPT IDL2

  ; Step 1: Morphological Opening with Horizontal Kernel (w=5)
  ; Create a structuring element that is 5 pixels wide and 1 pixel tall.
  ; BYTARR(5, 1) creates the array, +1B sets all values to 1.
  struc = BYTARR(5, 1) + 1B

  ; Apply Morphological Opening (Erosion followed by Dilation).
  ; This removes vertical features narrower than the kernel width.
  mask_opened = MORPH_OPEN(mask_image, struc)

  ; Step 2: Remove Small Clusters (< 15 pixels)
  ; Label connected components.
  ; Default LABEL_REGION uses 4-neighbor connectivity, similar to standard
  ; Python ndimage.label defaults.
  regions = LABEL_REGION(mask_opened)

  ; Calculate the population (number of pixels) of each labeled region.
  ; HISTOGRAM creates an array where index i contains the count of pixels with label i.
  region_counts = HISTOGRAM(regions)

  ; Map the counts back to the pixel coordinates.
  ; region_counts[regions] creates an image where each pixel's value is
  ; the size of the cluster it belongs to.
  pixel_cluster_sizes = region_counts[regions]

  ; Create the final mask.
  ; We keep a pixel if:
  ;   1. It was in the opened mask (mask_opened NE 0)
  ;   2. Its cluster size is >= 15 (pixel_cluster_sizes GE 15)
  ; Note: Background (0) usually has a large count, so we must check mask_opened to
  ; ensure we don't turn the background into 1s.
  mask_final = mask_opened AND (pixel_cluster_sizes GE 15)

  RETURN, mask_final
END

;+
; NAME:
;   KLIDARAS_CLASSIFY_ALTERATION_PHASES
;
; PURPOSE:
;   Classifies spectral pixels into 6 alteration phases based on their
;   Band Center Wavelength.
;
;   INPUT:
;     merged_center_wvl : 2D array of band centers (can be in nm or microns).
;     high_BD_mask      : Binary mask to apply to the final output.
;     id                : String identifier for filename saving.
;
;   OUTPUT:
;     A BD-masked 2D Integer array (0-6) where:
;     0 = Unclassified / Background
;     1 = Alunite/Beidellite                                        (2155 - 2186 nm)
;     2 = Allophane/Imogolite/Mont/Kaol/Hall/Illite/Opal            (2186 - 2220 nm)
;     3 = Opal/Jarosite/Acid-leached nontronite                     (2220 - 2280 nm)
;     4 = High Fe/Mg dioct. smectites (e.g., Nontronite)            (2280 - 2290 nm)
;     5 = High Fe/Mg di/tri smectites                               (2290 - 2307 nm)
;     6 = Low Fe/Mg tri smectites (e.g., Saponite)                  (2307 - 2320 nm)
;-
FUNCTION Klidaras_classify_alteration_phases, merged_center_wvl, high_BD_mask, id

  COMPILE_OPT IDL2

  ; 1. Setup Output Array
  sz = SIZE(merged_center_wvl, /STRUCTURE)
  alteration_phase_int_map = MAKE_ARRAY(sz.dimensions[0], sz.dimensions[1], /INTEGER, VALUE=0)

  ; 2. Handle Units (Microns vs Nanometers)
  ;    We normalize everything to Nanometers for comparison.
  max_val = MAX(merged_center_wvl, /NAN)

  IF max_val LT 10.0 THEN BEGIN
    ; Input is likely in Microns (e.g., 2.20)
    wvl_ref = merged_center_wvl * 1000.0
    PRINT, 'Detected input in Microns. Converting to nm for classification.'
  ENDIF ELSE BEGIN
    ; Input is likely in Nanometers (e.g., 2200)
    wvl_ref = merged_center_wvl
  ENDELSE

  PRINT, 'Classifying alteration phases by band center...'

  ; 3. Perform Vectorized Classification

  ; --- Phase 1: 2155 - 2186 nm ---
  ; Alunite/Beidellite
  idx1 = WHERE((wvl_ref GE 2155.0) AND (wvl_ref LT 2186.0), c1)
  IF c1 GT 0 THEN alteration_phase_int_map[idx1] = 1

  ; --- Phase 2: 2186 - 2220 nm ---
  ; Allophane/Imogolite/Mont/Kaol/Hall/Illite/Opal
  idx2 = WHERE((wvl_ref GE 2186.0) AND (wvl_ref LT 2220.0), c2)
  IF c2 GT 0 THEN alteration_phase_int_map[idx2] = 2

  ; --- Phase 3: 2220 - 2280 nm ---
  ; Opal/Jarosite/Acid-leached nontronite
  idx3 = WHERE((wvl_ref GE 2220.0) AND (wvl_ref LT 2280.0), c3)
  IF c3 GT 0 THEN alteration_phase_int_map[idx3] = 3

  ; --- Phase 4: 2280 - 2290 nm ---
  ; High Fe/Mg dioct. smectites (e.g., Nontronite)
  idx4 = WHERE((wvl_ref GE 2280.0) AND (wvl_ref LT 2290.0), c4)
  IF c4 GT 0 THEN alteration_phase_int_map[idx4] = 4

  ; --- Phase 5: 2290 - 2307 nm ---
  ; High Fe/Mg di/tri smectites
  idx5 = WHERE((wvl_ref GE 2290.0) AND (wvl_ref LT 2307.0), c5)
  IF c5 GT 0 THEN alteration_phase_int_map[idx5] = 5

  ; --- Phase 6: 2307 - 2320 nm ---
  ; Low Fe/Mg tri smectites (e.g., Saponite)
  ; Note: Keeping LE for the upper bound to include exactly 2320.0
  idx6 = WHERE((wvl_ref GE 2307.0) AND (wvl_ref LE 2320.0), c6)
  IF c6 GT 0 THEN alteration_phase_int_map[idx6] = 6

  ; 4. Report Statistics (does not execute to reduce terminal spam)
  ;PRINT, '--- Classification Stats ---'
  ;PRINT, 'Phase 1 (Alunite/Beidellite): ', c1, ' pixels'
  ;PRINT, 'Phase 2 (Kaol/Mont/Illite):   ', c2, ' pixels'
  ;PRINT, 'Phase 3 (Opal/Jarosite):      ', c3, ' pixels'
  ;PRINT, 'Phase 4 (Nontronite):         ', c4, ' pixels'
  ;PRINT, 'Phase 5 (Hi-Fe/Mg Smectite):  ', c5, ' pixels'
  ;PRINT, 'Phase 6 (Saponite):           ', c6, ' pixels'
  ;PRINT, 'Unclassified:                 ', N_ELEMENTS(alteration_phase_int_map) - (c1+c2+c3+c4+c5+c6)

  ; 5. Save Output
  filename = 'Saved data paper2/' + id + '_alteration_phase_map.csv'
  write_csv, filename, alteration_phase_int_map
  PRINT, 'File written: ' + filename

  filename = 'Saved data paper2/' + id + '_alteration_phase_map_BD_masked.csv'
  write_csv, filename, alteration_phase_int_map * high_BD_mask
  PRINT, 'File written: ' + filename

  RETURN, alteration_phase_int_map * high_BD_mask

END

;+
; NAME:
;   KLIDARAS_CLASSIFY_ALTERATION_PHASES
;
; PURPOSE:
;   Classifies spectral pixels into 6 alteration phases based on their
;   Band Center Wavelength.
;
;   INPUT:
;     merged_center_wvl : 2D array of band centers (can be in nm or microns).
;     high_BD_mask      : Binary mask to apply to the final output.
;     id                : String identifier for filename saving.
;
;   OUTPUT:
;     A BD-masked 2D Integer array (0-6) where:
;     0 = Unclassified / Background
;     1 = Alunite/Beidellite                                        (2155 - 2186 nm)
;     2 = Allophane/Imogolite/Mont/Kaol/Hall/Illite/Opal            (2186 - 2220 nm)
;     3 = Opal/Jarosite/Acid-leached nontronite                     (2220 - 2280 nm)
;     4 = High Fe/Mg dioct. smectites (e.g., Nontronite)            (2280 - 2290 nm)
;     5 = High Fe/Mg di/tri smectites                               (2290 - 2307 nm)
;     6 = Low Fe/Mg tri smectites (e.g., Saponite)                  (2307 - 2320 nm)
;-
FUNCTION Klidaras_classify_alteration_phases, merged_center_wvl, high_BD_mask, id

  COMPILE_OPT IDL2

  ; 1. Setup Output Array
  sz = SIZE(merged_center_wvl, /STRUCTURE)
  alteration_phase_int_map = MAKE_ARRAY(sz.dimensions[0], sz.dimensions[1], /INTEGER, VALUE=0)

  ; 2. Handle Units (Microns vs Nanometers)
  ;    We normalize everything to Nanometers for comparison.
  max_val = MAX(merged_center_wvl, /NAN)

  IF max_val LT 10.0 THEN BEGIN
    ; Input is likely in Microns (e.g., 2.20)
    wvl_ref = merged_center_wvl * 1000.0
    PRINT, 'Detected input in Microns. Converting to nm for classification.'
  ENDIF ELSE BEGIN
    ; Input is likely in Nanometers (e.g., 2200)
    wvl_ref = merged_center_wvl
  ENDELSE

  PRINT, 'Classifying alteration phases by band center...'

  ; 3. Perform Vectorized Classification

  ; --- Phase 1: 2155 - 2186 nm ---
  ; Alunite/Beidellite
  idx1 = WHERE((wvl_ref GE 2155.0) AND (wvl_ref LT 2186.0), c1)
  IF c1 GT 0 THEN alteration_phase_int_map[idx1] = 1

  ; --- Phase 2: 2186 - 2220 nm ---
  ; Allophane/Imogolite/Mont/Kaol/Hall/Illite/Opal
  idx2 = WHERE((wvl_ref GE 2186.0) AND (wvl_ref LT 2220.0), c2)
  IF c2 GT 0 THEN alteration_phase_int_map[idx2] = 2

  ; --- Phase 3: 2220 - 2280 nm ---
  ; Opal/Jarosite/Acid-leached nontronite
  idx3 = WHERE((wvl_ref GE 2220.0) AND (wvl_ref LT 2280.0), c3)
  IF c3 GT 0 THEN alteration_phase_int_map[idx3] = 3

  ; --- Phase 4: 2280 - 2290 nm ---
  ; High Fe/Mg dioct. smectites (e.g., Nontronite)
  idx4 = WHERE((wvl_ref GE 2280.0) AND (wvl_ref LT 2290.0), c4)
  IF c4 GT 0 THEN alteration_phase_int_map[idx4] = 4

  ; --- Phase 5: 2290 - 2307 nm ---
  ; High Fe/Mg di/tri smectites
  idx5 = WHERE((wvl_ref GE 2290.0) AND (wvl_ref LT 2307.0), c5)
  IF c5 GT 0 THEN alteration_phase_int_map[idx5] = 5

  ; --- Phase 6: 2307 - 2320 nm ---
  ; Low Fe/Mg tri smectites (e.g., Saponite)
  ; Note: Keeping LE for the upper bound to include exactly 2320.0
  idx6 = WHERE((wvl_ref GE 2307.0) AND (wvl_ref LE 2320.0), c6)
  IF c6 GT 0 THEN alteration_phase_int_map[idx6] = 6

  ; 4. Report Statistics (does not execute to reduce terminal spam)
  ;PRINT, '--- Classification Stats ---'
  ;PRINT, 'Phase 1 (Alunite/Beidellite): ', c1, ' pixels'
  ;PRINT, 'Phase 2 (Kaol/Mont/Illite):   ', c2, ' pixels'
  ;PRINT, 'Phase 3 (Opal/Jarosite):      ', c3, ' pixels'
  ;PRINT, 'Phase 4 (Nontronite):         ', c4, ' pixels'
  ;PRINT, 'Phase 5 (Hi-Fe/Mg Smectite):  ', c5, ' pixels'
  ;PRINT, 'Phase 6 (Saponite):           ', c6, ' pixels'
  ;PRINT, 'Unclassified:                 ', N_ELEMENTS(alteration_phase_int_map) - (c1+c2+c3+c4+c5+c6)

  ; 5. Save Output
  filename = 'Saved data paper2/' + id + '_alteration_phase_map.csv'
  write_csv, filename, alteration_phase_int_map
  PRINT, 'File written: ' + filename

  filename = 'Saved data paper2/' + id + '_alteration_phase_map_BD_masked.csv'
  write_csv, filename, alteration_phase_int_map * high_BD_mask
  PRINT, 'File written: ' + filename

  RETURN, alteration_phase_int_map * high_BD_mask

END

;Usage:
;Kaolin_param = Klidaras_kaolinite_parameter(if_1,crism_wvl_array)
;Updates:
;12/3/2025: writes masked BD_2165/BD_2205 & masked IBD Assymetry Ratio out as .csv's
;1/7/2026: added kaolin param (mask where pixels which meet kaolinite criteria have value of 1)

FUNCTION Klidaras_Kaolinite_parameter, cube, wvl, id, masked_phyllosilicate_int_array

  ;Creating mask to exclude non-Al-unit pixels (which are noisy and thus have high variance of numbers)
  s=size(masked_phyllosilicate_int_array,/dim)
  Al_unit_mask=fltarr(s[0],s[1])
  ;Mask where band center phases 1 & 2 (Al-rich and Al/Si-rich) are present in the masked phyllosilicate phase map
  Al_unit_mask[where(masked_phyllosilicate_int_array eq 1)]=1
  Al_unit_mask[where(masked_phyllosilicate_int_array eq 2)]=1


  ;(1) Calculating the two requisite band depths (at 2165 nm and 2205 nm)
  BD_2165_Klidaras=calc_band_depth(cube,2120,2165,2240,3,1,3,id,Al_unit_mask,WAVELENGTHS=wvl)
  BD_2205_Klidaras=calc_band_depth(cube,2120,2205,2240,3,1,3,id,Al_unit_mask,WAVELENGTHS=wvl)
  ;Saving each to .csv files
  filename='Saved data paper2/'+id+'_Kaolin_param_BD2165_Klidaras.csv'
  write_csv,filename,BD_2165_Klidaras
  print,'File written: '+filename
  filename='Saved data paper2/'+id+'_Kaolin_param_BD2205_Klidaras.csv'
  write_csv,filename,BD_2205_Klidaras
  print,'File written: '+filename
  filename='Saved data paper2/'+id+'_Kaolin_param_masked_BD2165_Klidaras.csv'
  write_csv,filename,BD_2165_Klidaras*Al_unit_mask
  print,'File written: '+filename
  filename='Saved data paper2/'+id+'_Kaolin_param_masked_BD2205_Klidaras.csv'
  write_csv,filename,BD_2205_Klidaras*Al_unit_mask
  print,'File written: '+filename

  ;Changing bad pixels from value of 0. (results in divide by zero error) to 1
  bad_pixels=where(BD_2205_Klidaras eq 0.)
  BD_2165_Klidaras[bad_pixels]=1
  BD_2205_Klidaras[bad_pixels]=1

  ;(2) Ratio BD_2165/BD_2205 (kaolinite/halloysite = ~0.6)
  ratio=BD_2165_Klidaras/BD_2205_Klidaras
  ratio[bad_pixels]=0
  ;Saving...
  filename='Saved data paper2/'+id+'_Kaolin_param_BD2165_BD2205_ratio.csv'
  write_csv,filename,ratio
  print,'File written: '+filename
  filename='Saved data paper2/'+id+'_Kaolin_param_masked_BD2165_BD2205_ratio.csv'
  write_csv,filename,ratio*Al_unit_mask
  print,'File written: '+filename

  ;(3) Integrated Band Depth (IBD)
  IBD=calc_integrated_band_depth(cube,2120,2240,3,3,id,Al_unit_mask,wavelengths=wvl)
  ;Saving...
  filename='Saved data paper2/'+id+'_Kaolin_param_IBD.csv'
  write_csv,filename,IBD
  print,'File written: '+filename
  filename='Saved data paper2/'+id+'_Kaolin_param_masked_IBD.csv'
  write_csv,filename,IBD*Al_unit_mask
  print,'File written: '+filename


  ;(4) Calculating IBD Assymetry Ratio (IBD Left/IBD Right, kaolinite/halloysite = ~2-8)
  IBD_Asymmetry_Ratio=calc_IBD_asymmetry_ratio(cube,2120,2240,3,3,id,Al_unit_mask,wavelengths=wvl)
  ;Saving...
  filename='Saved data paper2/'+id+'_Kaolin_param_IBD_Asymmetry_Ratio.csv'
  write_csv,filename,IBD_Asymmetry_Ratio
  print,'File written: '+filename
  filename='Saved data paper2/'+id+'_Kaolin_param_masked_IBD_Asymmetry_Ratio.csv'
  write_csv,filename,IBD_Asymmetry_Ratio*Al_unit_mask
  print,'File written: '+filename

  ;(5) Classifying spectra dominated by kaolin
  Kaolin_param=classify_kaolin(ratio,IBD_Asymmetry_Ratio,Al_unit_mask)
  ;Saving...
  filename='Saved data paper2/'+id+'_Kaolin_param.csv'
  write_csv,filename,Kaolin_param
  print,'File written: '+filename
  filename='Saved data paper2/'+id+'_Kaolin_param_masked.csv'
  write_csv,filename,Kaolin_param*Al_unit_mask
  print,'File written: '+filename


  return,Kaolin_param*Al_unit_mask

END

;+
; NAME:
;   CALC_BAND_DEPTH
; Created by Thanos Klidaras on 11.24.2025
;
; PURPOSE:
;   Calculates the continuum-removed band depth for a specific absorption feature
;   within a spectral cube.
;
;   Formula: BD = 1 - (Reflectance_Center / Continuum_Center)
;
; ARGUMENTS:
;   cube            : Input 3D spectral cube [x, y, bands]
;   in_left         : Left shoulder position (Band Index OR Wavelength)
;   in_center       : Absorption center position (Band Index OR Wavelength)
;   in_right        : Right shoulder position (Band Index OR Wavelength)
;   in_left_width   : (Optional) Width of left shoulder in bands (Default: 1)
;   in_center_width : (Optional) Width of center band in bands (Default: 1)
;   in_right_width  : (Optional) Width of right shoulder in bands (Default: 1)
;                     * If width > 1, the code takes the MEAN of the bands
;                       centered on the target.
;
; KEYWORDS:
;   WAVELENGTHS : An array of wavelength values corresponding to the bands.
;                 IF SET: The inputs (in_left, etc.) are treated as Wavelengths.
;                         The code finds the closest spectral band to these values.
;                 IF NOT SET: The inputs are treated as integer Band Indices.
;
; RETURN:
;   A 2D array [x, y] containing the band depth values. Returns -1 on failure.
;-
FUNCTION CALC_BAND_DEPTH, cube, in_left, in_center, in_right, in_left_width, in_center_width, in_right_width, id, Al_unit_mask, WAVELENGTHS=wvl

  COMPILE_OPT IDL2

  ; 1. Get Dimensions
  sz = SIZE(cube, /STRUCTURE)
  IF sz.n_dimensions NE 3 THEN BEGIN
    MESSAGE, 'Input must be a 3D spectral cube.', /CONTINUE
    RETURN, -1
  ENDIF

  nx = sz.dimensions[0]
  ny = sz.dimensions[1]
  nb = sz.dimensions[2]

  ; 2. Process Width Arguments (Defaults to 1 if not set)
  IF N_ELEMENTS(in_left_width)   EQ 0 THEN w_L = 1 ELSE w_L = FIX(in_left_width)
  IF N_ELEMENTS(in_center_width) EQ 0 THEN w_C = 1 ELSE w_C = FIX(in_center_width)
  IF N_ELEMENTS(in_right_width)  EQ 0 THEN w_R = 1 ELSE w_R = FIX(in_right_width)

  ; 3. Determine Mode (Wavelengths vs Indices) and Find Indices
  IF N_ELEMENTS(wvl) GT 0 THEN BEGIN
    ; --- WAVELENGTH MODE ---
    IF N_ELEMENTS(wvl) NE nb THEN BEGIN
      MESSAGE, 'Wavelength array size does not match cube bands.', /CONTINUE
      RETURN, -1
    ENDIF

    ; Find closest band index for the requested wavelengths
    dummy = MIN(ABS(wvl - in_left), idx_left)
    dummy = MIN(ABS(wvl - in_right), idx_right)
    dummy = MIN(ABS(wvl - in_center), idx_center)

    ; Define X-axis values for interpolation using the ACTUAL wavelengths of the found bands
    x_L = float(wvl[idx_left])
    x_R = float(wvl[idx_right])
    x_C = float(wvl[idx_center])

    PRINT, '--- Auto-Selected Bands ---'

    ; --- LEFT PRINTOUT ---
    ; Logic: Start = Index - (Width/2). End = Start + Width - 1.
    s = (idx_left - (w_L / 2)) > 0
    e = (s + w_L - 1) < (nb - 1)
    ; Re-clamp start if end forced a clip? No, we prefer to clip range at edges.
    vals = wvl[s:e]
    v_str = STRJOIN(STRTRIM(vals, 2), ', ')
    PRINT, 'Left:   Req=', STRTRIM(in_left, 2), ', Width=', STRTRIM(w_L, 2), ', Actual=', v_str

    ; --- CENTER PRINTOUT ---
    s = (idx_center - (w_C / 2)) > 0
    e = (s + w_C - 1) < (nb - 1)
    vals = wvl[s:e]
    v_str = STRJOIN(STRTRIM(vals, 2), ', ')
    PRINT, 'Center: Req=', STRTRIM(in_center, 2), ', Width=', STRTRIM(w_C, 2), ', Actual=', v_str

    ; --- RIGHT PRINTOUT ---
    s = (idx_right - (w_R / 2)) > 0
    e = (s + w_R - 1) < (nb - 1)
    vals = wvl[s:e]
    v_str = STRJOIN(STRTRIM(vals, 2), ', ')
    PRINT, 'Right:  Req=', STRTRIM(in_right, 2), ', Width=', STRTRIM(w_R, 2), ', Actual=', v_str

  ENDIF ELSE BEGIN
    ; --- INDEX MODE ---
    idx_left   = FIX(in_left)
    idx_right  = FIX(in_right)
    idx_center = FIX(in_center)

    ; Define X-axis values (just indices)
    x_L = float(idx_left)
    x_R = float(idx_right)
    x_C = float(idx_center)
  ENDELSE

  ; 4. Error Checking: Bounds & Logic

  ; A. Check if indices are inside cube dimensions
  IF (idx_left LT 0) OR (idx_left GE nb) OR $
    (idx_right LT 0) OR (idx_right GE nb) OR $
    (idx_center LT 0) OR (idx_center GE nb) THEN BEGIN
    MESSAGE, 'Calculated or provided band indices are out of cube range.', /CONTINUE
    RETURN, -1
  ENDIF

  ; B. Check if Center is strictly between Left and Right
  IF (x_C LE x_L) OR (x_C GE x_R) THEN BEGIN
    MESSAGE, 'ERROR: Band Center must be strictly between Left and Right shoulders.', /CONTINUE
    PRINT, 'Left: ', x_L, ' Center: ', x_C, ' Right: ', x_R
    RETURN, -1
  ENDIF

  ; 5. Calculate Interpolation Weight
  weight = (x_C - x_L) / (x_R - x_L)

  ; 6. Extract 2D Slices (With Width Averaging)

  ; --- LEFT SHOULDER EXTRACTION ---
  s = (idx_left - (w_L / 2)) > 0
  e = (s + w_L - 1) < (nb - 1)

  IF s EQ e THEN val_L = cube[*,*,s] $
  ELSE val_L = TOTAL(cube[*,*,s:e], 3) / FLOAT(e - s + 1)

  ; --- RIGHT SHOULDER EXTRACTION ---
  s = (idx_right - (w_R / 2)) > 0
  e = (s + w_R - 1) < (nb - 1)

  IF s EQ e THEN val_R = cube[*,*,s] $
  ELSE val_R = TOTAL(cube[*,*,s:e], 3) / FLOAT(e - s + 1)

  ; --- CENTER EXTRACTION ---
  s = (idx_center - (w_C / 2)) > 0
  e = (s + w_C - 1) < (nb - 1)

  IF s EQ e THEN val_C = cube[*,*,s] $
  ELSE val_C = TOTAL(cube[*,*,s:e], 3) / FLOAT(e - s + 1)

  ; 7. Calculate Continuum at Center position
  ; Linear interpolation: y = y1 + (y2 - y1) * weight
  continuum_val = val_L + (val_R - val_L) * weight

  ; 8. Calculate Band Depth
  ; Standard definition: Depth = 1 - (Value / Continuum)
  depth_map = MAKE_ARRAY(nx, ny, /FLOAT, VALUE=0.0)

  ; Perform calculation handling division by zero
  valid_mask = WHERE(FINITE(continuum_val) AND (continuum_val NE 0.0), count)

  IF count GT 0 THEN BEGIN
    depth_map[valid_mask] = 1.0 - (val_C[valid_mask] / continuum_val[valid_mask])
  ENDIF

  RETURN, depth_map

END

;+
; NAME:
;   CALC_INTEGRATED_BAND_DEPTH
;
; PURPOSE:
;   Calculates the Integrated Band Depth (IBD) for an absorption feature.
;   Unlike standard Band Depth (which measures one point), IBD calculates the
;   sum of depths for all bands between the left and right shoulders.
;
;   Formula: IBD = SUM( 1 - (Reflectance_i / Continuum_i) )
;            for all bands i from Left_Shoulder to Right_Shoulder.
;
; ARGUMENTS:
;   cube            : Input 3D spectral cube [x, y, bands]
;   in_left         : Left shoulder position (Start of integration)
;   in_right        : Right shoulder position (End of integration)
;   in_left_width   : (Optional) Width of left shoulder anchor (averaged). Default: 1
;   in_right_width  : (Optional) Width of right shoulder anchor (averaged). Default: 1
;
; KEYWORDS:
;   WAVELENGTHS : An array of wavelength values.
;
; RETURN:
;   A 2D array [x, y] containing the integrated band depth area.
;-
FUNCTION CALC_INTEGRATED_BAND_DEPTH, cube, in_left, in_right, in_left_width, in_right_width, id, Al_unit_mask, WAVELENGTHS=wvl

  COMPILE_OPT IDL2

  ; 1. Get Dimensions
  sz = SIZE(cube, /STRUCTURE)
  IF sz.n_dimensions NE 3 THEN BEGIN
    MESSAGE, 'Input must be a 3D spectral cube.', /CONTINUE
    RETURN, -1
  ENDIF

  nx = sz.dimensions[0]
  ny = sz.dimensions[1]
  nb = sz.dimensions[2]

  ; 2. Process Width Arguments (Defaults to 1 if not set)
  IF N_ELEMENTS(in_left_width)   EQ 0 THEN w_L = 1 ELSE w_L = FIX(in_left_width)
  IF N_ELEMENTS(in_right_width)  EQ 0 THEN w_R = 1 ELSE w_R = FIX(in_right_width)

  ; 3. Determine Indices based on Wavelengths or Direct Input
  IF N_ELEMENTS(wvl) GT 0 THEN BEGIN
    ; --- WAVELENGTH MODE ---
    IF N_ELEMENTS(wvl) NE nb THEN BEGIN
      MESSAGE, 'Wavelength array size does not match cube bands.', /CONTINUE
      RETURN, -1
    ENDIF

    ; Find closest band indices
    diff_L = MIN(ABS(wvl - in_left), idx_left)
    diff_R = MIN(ABS(wvl - in_right), idx_right)

    ; Get actual wavelengths for reporting
    x_L = float(wvl[idx_left])
    x_R = float(wvl[idx_right])

    PRINT, '--- Auto-Selected Bands (Integration Range) ---'

    ; --- UNIT MISMATCH CHECK ---
    IF (diff_L GT 50.0) OR (diff_R GT 50.0) THEN BEGIN
      PRINT, '>>> WARNING: LARGE DISCREPANCY DETECTED <<<'
      PRINT, '    Requested wavelengths are significantly different from found bands.'
      PRINT, '    Likely Unit Mismatch (e.g., Micrometers vs Nanometers).'
      PRINT, '    Please check your inputs.'
    ENDIF

    ; Report Ranges
    s = (idx_left - (w_L / 2)) > 0
    e = (s + w_L - 1) < (nb - 1)
    v_str_L = STRJOIN(STRTRIM(wvl[s:e], 2), ', ')

    s = (idx_right - (w_R / 2)) > 0
    e = (s + w_R - 1) < (nb - 1)
    v_str_R = STRJOIN(STRTRIM(wvl[s:e], 2), ', ')

    PRINT, 'Start (Left): Req=', STRTRIM(in_left, 2), ', Anchor Wvl=', v_str_L
    PRINT, 'End (Right):  Req=', STRTRIM(in_right, 2), ', Anchor Wvl=', v_str_R
    PRINT, 'Integration spans from Band ', strtrim(idx_left,2), ' to ', strtrim(idx_right,2)

  ENDIF ELSE BEGIN
    ; --- INDEX MODE ---
    idx_left   = FIX(in_left)
    idx_right  = FIX(in_right)
    x_L = float(idx_left)
    x_R = float(idx_right)
  ENDELSE

  ; 4. Error Checking
  IF (idx_left LT 0) OR (idx_right GE nb) THEN BEGIN
    MESSAGE, 'Indices out of range.', /CONTINUE
    RETURN, -1
  ENDIF

  ; Ensure Left comes before Right
  IF idx_left GE idx_right THEN BEGIN
    MESSAGE, 'ERROR: Left shoulder must be at a lower index/wavelength than Right shoulder.', /CONTINUE
    RETURN, -1
  ENDIF

  ; 5. Calculate Anchor Values (Continuum Endpoints)
  ; --- LEFT ANCHOR ---
  s = (idx_left - (w_L / 2)) > 0
  e = (s + w_L - 1) < (nb - 1)
  IF s EQ e THEN val_L = cube[*,*,s] $
  ELSE val_L = TOTAL(cube[*,*,s:e], 3) / FLOAT(e - s + 1)

  ; --- RIGHT ANCHOR ---
  s = (idx_right - (w_R / 2)) > 0
  e = (s + w_R - 1) < (nb - 1)
  IF s EQ e THEN val_R = cube[*,*,s] $
  ELSE val_R = TOTAL(cube[*,*,s:e], 3) / FLOAT(e - s + 1)


  ; 6. Perform Integration Loop
  ; We iterate from Left Index to Right Index.

  ; Initialize Output Array
  ibd_map = MAKE_ARRAY(nx, ny, /FLOAT, VALUE=0.0)

  ; Pre-calculate Continuum Slope for every pixel
  ; Slope = (Right_Val - Left_Val) / Total_Distance
  dist_total = FLOAT(x_R - x_L)

  ; Avoid divide by zero if dist is 0 (unlikely given check above)
  IF dist_total EQ 0 THEN dist_total = 1.0

  slope_map = (val_R - val_L) / dist_total

  PRINT, 'Calculating Integrated Band Depth...'

  ; Loop through every band in the range
  FOR i = idx_left, idx_right DO BEGIN

    ; A. Get current X position (Wavelength or Index)
    IF N_ELEMENTS(wvl) GT 0 THEN cur_x = FLOAT(wvl[i]) ELSE cur_x = FLOAT(i)

    ; B. Calculate Continuum at this specific band i
    ; Cont = Start + (Slope * DistanceFromStart)
    dist_curr = cur_x - x_L
    cont_i    = val_L + (slope_map * dist_curr)

    ; C. Get actual Reflectance at band i
    spec_i    = cube[*,*,i]

    ; D. Calculate Depth: 1 - (Spectrum / Continuum)
    ; Handle division by zero carefully
    term = MAKE_ARRAY(nx, ny, /FLOAT, VALUE=0.0)

    ; Only calculate where continuum is valid
    good_mask = WHERE( (cont_i NE 0.0) AND FINITE(cont_i), count)

    IF count GT 0 THEN BEGIN
      term[good_mask] = 1.0 - (spec_i[good_mask] / cont_i[good_mask])
    ENDIF

    ; E. Add to Total (Integration)
    ; We ignore pixels where the calculation failed (remained 0)
    ibd_map = ibd_map + term

  ENDFOR

  ; 7. Finish

  RETURN, ibd_map

END

;+
; NAME:
;   CALC_IBD_ASYMMETRY_RATIO
;
; PURPOSE:
;   Calculates the Asymmetry Ratio (S) of an absorption feature using
;   Integrated Band Depth (IBD).
;
;   The function automatically finds the continuum-removed "Band Center"
;   (deepest point) for each pixel, then splits the integration into
;   Left and Right halves.
;
;   Formula: S = IBD_Left / IBD_Right
;     S = 1.0 : Symmetric
;     S > 1.0 : Left-Asymmetric (Tail to the left)
;     S < 1.0 : Right-Asymmetric (Tail to the right)
;
; ARGUMENTS:
;   cube            : Input 3D spectral cube [x, y, bands]
;   in_left         : Left shoulder position (Band Index OR Wavelength)
;   in_right        : Right shoulder position (Band Index OR Wavelength)
;   in_left_width   : (Optional) Width of left shoulder average (Default: 1)
;   in_right_width  : (Optional) Width of right shoulder average (Default: 1)
;
; KEYWORDS:
;   WAVELENGTHS : An array of wavelength values.
;                 IF SET: Inputs are treated as Wavelengths.
;
; RETURN:
;   A 2D array [x, y] containing the Asymmetry Ratio.
;   Returns -1 on failure.
;-
FUNCTION CALC_IBD_ASYMMETRY_RATIO, cube, in_left, in_right, in_left_width, in_right_width, id, Al_unit_mask, WAVELENGTHS=wvl

  COMPILE_OPT IDL2

  ; 1. Get Dimensions
  sz = SIZE(cube, /STRUCTURE)
  IF sz.n_dimensions NE 3 THEN BEGIN
    MESSAGE, 'Input must be a 3D spectral cube.', /CONTINUE
    RETURN, -1
  ENDIF

  nx = sz.dimensions[0]
  ny = sz.dimensions[1]
  nb = sz.dimensions[2]

  ; 2. Process Width Arguments
  IF N_ELEMENTS(in_left_width)   EQ 0 THEN w_L = 1 ELSE w_L = FIX(in_left_width)
  IF N_ELEMENTS(in_right_width)  EQ 0 THEN w_R = 1 ELSE w_R = FIX(in_right_width)

  ; 3. Determine Indices (Wavelengths vs Indices)
  IF N_ELEMENTS(wvl) GT 0 THEN BEGIN
    ; --- WAVELENGTH MODE ---
    IF N_ELEMENTS(wvl) NE nb THEN BEGIN
      MESSAGE, 'Wavelength array size does not match cube bands.', /CONTINUE
      RETURN, -1
    ENDIF

    ; Find closest band indices for shoulders
    dummy = MIN(ABS(wvl - in_left), idx_left)
    dummy = MIN(ABS(wvl - in_right), idx_right)

    ; Get actual wavelengths for calculation/reporting
    x_L = float(wvl[idx_left])
    x_R = float(wvl[idx_right])

    PRINT, '--- Asymmetry Calculation Setup ---'
    PRINT, 'Range: ', STRTRIM(x_L, 2), ' to ', STRTRIM(x_R, 2)
    PRINT, 'Indices: ', STRTRIM(idx_left, 2), ' to ', STRTRIM(idx_right, 2)

  ENDIF ELSE BEGIN
    ; --- INDEX MODE ---
    idx_left   = FIX(in_left)
    idx_right  = FIX(in_right)
    x_L = float(idx_left)
    x_R = float(idx_right)
  ENDELSE

  ; 4. Error Checking
  IF (idx_left LT 0) OR (idx_right GE nb) THEN BEGIN
    MESSAGE, 'Indices out of range.', /CONTINUE
    RETURN, -1
  ENDIF

  IF idx_left GE idx_right THEN BEGIN
    MESSAGE, 'ERROR: Left shoulder must be < Right shoulder.', /CONTINUE
    RETURN, -1
  ENDIF

  ; 5. Calculate Anchor Values (Continuum Endpoints)
  ; --- LEFT ANCHOR ---
  s = (idx_left - (w_L / 2)) > 0
  e = (s + w_L - 1) < (nb - 1)
  IF s EQ e THEN val_L = cube[*,*,s] $
  ELSE val_L = TOTAL(cube[*,*,s:e], 3) / FLOAT(e - s + 1)

  ; --- RIGHT ANCHOR ---
  s = (idx_right - (w_R / 2)) > 0
  e = (s + w_R - 1) < (nb - 1)
  IF s EQ e THEN val_R = cube[*,*,s] $
  ELSE val_R = TOTAL(cube[*,*,s:e], 3) / FLOAT(e - s + 1)

  ; 6. Prepare Continuum and Depth Cube
  ; We need to store the depth for every band in the range to finding the max depth.

  range_len = idx_right - idx_left + 1

  ; Create a temporary 3D array to hold Band Depths for the specific range
  ; Dimensions: [nx, ny, range_len]
  depth_cube = MAKE_ARRAY(nx, ny, range_len, /FLOAT, VALUE=0.0)

  ; Calculate Continuum Slope
  dist_total = FLOAT(x_R - x_L)
  IF dist_total EQ 0 THEN dist_total = 1.0
  slope_map = (val_R - val_L) / dist_total

  PRINT, 'Generating Continuum-Removed Cube...'

  ; Loop to fill depth_cube
  FOR i = 0, range_len - 1 DO BEGIN
    actual_band_idx = idx_left + i

    ; X-coordinate for this band
    IF N_ELEMENTS(wvl) GT 0 THEN cur_x = FLOAT(wvl[actual_band_idx]) $
    ELSE cur_x = FLOAT(actual_band_idx)

    ; Calculate Continuum
    dist_curr = cur_x - x_L
    cont_i    = val_L + (slope_map * dist_curr)

    ; Get Spectrum
    spec_i    = cube[*,*,actual_band_idx]

    ; Calculate Depth: 1 - (Refl/Cont)
    ; Avoid div by zero
    good_mask = WHERE((cont_i NE 0.0) AND FINITE(cont_i), count)
    IF count GT 0 THEN BEGIN
      depth_cube[*,*,i] = 1.0 - (spec_i / cont_i)
    ENDIF
  ENDFOR

  ; 7. Find Dynamic Center (Deepest Absorption)
  PRINT, 'Finding dynamic band centers...'

  ; MAX with DIMENSION=3 finds the max value along the Z-axis of our depth_cube.
  ; FIX: 'pos_idx' must be the second positional argument, not a keyword.
  max_val = MAX(depth_cube, pos_idx, DIMENSION=3)

  ; Convert Linear Index to Subscript (Z-index relative to range)
  ; In a [nx, ny, nz] array, Z = Linear / (nx*ny)
  plane_size = nx * ny
  center_offset_map = pos_idx / plane_size

  ; center_offset_map is now a 2D array [nx, ny] containing integers
  ; from 0 to (range_len-1), indicating where the deepest point is.

  ; 8. Calculate IBD_Left and IBD_Right
  PRINT, 'Integrating Split Areas...'

  ibd_left  = MAKE_ARRAY(nx, ny, /FLOAT, VALUE=0.0)
  ibd_right = MAKE_ARRAY(nx, ny, /FLOAT, VALUE=0.0)

  FOR i = 0, range_len - 1 DO BEGIN
    ; Extract depth slice for this step
    d_slice = depth_cube[*,*,i]

    ; --- LEFT SIDE INTEGRATION ---
    ; Include if current index (i) is LE center index
    mask_L = (i LE center_offset_map)
    ibd_left = ibd_left + (d_slice * mask_L)

    ; --- RIGHT SIDE INTEGRATION ---
    ; Include if current index (i) is GE center index
    ; Note: The center band is effectively counted in BOTH sums.
    ; This is standard for shape ratios to avoid dropping the peak data.
    mask_R = (i GE center_offset_map)
    ibd_right = ibd_right + (d_slice * mask_R)
  ENDFOR

  ; 9. Calculate Asymmetry Ratio (S = Left / Right)
  ratio_map = MAKE_ARRAY(nx, ny, /FLOAT, VALUE=-1.0) ; Default to -1 or NaN

  ; Handle Division by Zero (if right side has no area)
  valid_mask = WHERE((ibd_right NE 0.0) AND FINITE(ibd_right), count)

  IF count GT 0 THEN BEGIN
    ratio_map[valid_mask] = ibd_left[valid_mask] / ibd_right[valid_mask]
  ENDIF

  ; Optional: Clean up pixels where no absorption existed (max depth <= 0)
  ; If the feature is flat or convex, the ratio is meaningless.
  no_feat = WHERE(max_val LE 0.0, nf_count)
  IF nf_count GT 0 THEN ratio_map[no_feat] = 0.0

  PRINT, 'Calculation Complete.'

  RETURN, ratio_map

END

; NAME:
;   CLASSIFY_KAOLIN
;
; PURPOSE:
;   Classifies pixels based on specific thresholds for 'ratio' and
;   'IBD_Asymmetry_Ratio'. Returns a classification map with integer values
;   1 through 5, with 0 as the background.
;
; ARGUMENTS:
;   ratio               : Input array (any dimension) representing the ratio.
;   IBD_Asymmetry_Ratio : Input array (same dimension as ratio).
;
; RETURN:
;   An integer array of the same dimensions as the input.

FUNCTION CLASSIFY_KAOLIN, ratio, IBD_Asymmetry_Ratio,mask

  COMPILE_OPT IDL2

  s=size(ratio)


  ; 1. Initialize the result array with 0s.
  ; Multiplying by 0 and converting to FIX ensures we have an integer array
  ; with the exact same dimensions as the input 'ratio'.
  result = FIX(ratio * 0)

  ; 2. Define conditions using WHERE.
  ; The WHERE function finds the indices of the array that meet the criteria.
  ; We check if 'count' is greater than 0 to avoid errors if no pixels match.

  ; Kaolinite/halloysite
  ; Ratio: 0.4 to 0.9 | IBD Asymmetry Ratio: 2 to 7
  idx1 = WHERE((ratio GE 0.4) AND (ratio LE 0.9) AND $
    (IBD_Asymmetry_Ratio GE 2.0) AND (IBD_Asymmetry_Ratio LE 7.0), count1)
  IF (count1 GT 0) THEN result[idx1] = 1

  RETURN, result

END

; NAME:
;   CALC_SPECTRAL_SLOPE
;
; UPDATES:
;   6/8/2026: Added optional mode to switch kernel averaging mode from mean to median to fix bug caused by occasional noisy CRISM channel at ~1200 nm
;
; PURPOSE:
;   Calculates the spectral slope between two points in a spectral cube.
;
;   Formula: Slope = (Reflectance_Right - Reflectance_Left) / (X_Right - X_Left)
;
; ARGUMENTS:
;   cube             : Input 3D spectral cube [x, y, bands]
;   in_left          : Left shoulder position (Band Index OR Wavelength)
;   in_right         : Right shoulder position (Band Index OR Wavelength)
;   in_left_width    : (Optional) Width of left shoulder in bands (Default: 1)
;   in_right_width   : (Optional) Width of right shoulder in bands (Default: 1)
;                      * If width > 1, the code averages the bands centered
;                        on the target.
;
; KEYWORDS:
;   WAVELENGTHS  : An array of wavelength values corresponding to the bands.
;                  IF SET: The inputs (in_left, in_right) are treated as Wavelengths.
;                          The denominator (dX) is calculated using real wavelength units.
;                  IF NOT SET: The inputs are treated as integer Band Indices.
;                          The denominator (dX) is calculated using band indices.
;   AVERAGE_MODE : String specifying how to average across the kernel width.
;                  'mean'   - (Default) Uses the standard mean (TOTAL / n).
;                  'median' - Uses the MEDIAN to ignore extreme data spikes.
;
; RETURN:
;   A 2D array [x, y] containing the slope values. Returns -1 on failure.
;-
; EXAMPLE USAGE:
;   output=calc_spectral_slope(cube, 1200, 1800, 7, 7, name, id, WAVELENGTHS=wvl, AVERAGE_MODE='median')

FUNCTION CALC_SPECTRAL_SLOPE, cube, in_left, in_right, in_left_width, in_right_width, name, id, WAVELENGTHS=wvl, AVERAGE_MODE=average_mode

  COMPILE_OPT IDL2

  ; 1. Get Dimensions
  sz = SIZE(cube, /STRUCTURE)
  IF sz.n_dimensions NE 3 THEN BEGIN
    MESSAGE, 'Input must be a 3D spectral cube.', /CONTINUE
    RETURN, -1
  ENDIF

  nx = sz.dimensions[0]
  ny = sz.dimensions[1]
  nb = sz.dimensions[2]

  ; 2. Process Width Arguments (Defaults to 1 if not set)
  IF N_ELEMENTS(in_left_width)  EQ 0 THEN w_L = 1 ELSE w_L = FIX(in_left_width)
  IF N_ELEMENTS(in_right_width) EQ 0 THEN w_R = 1 ELSE w_R = FIX(in_right_width)

  ; Parse Average Mode (Default to Mean)
  use_median = 0
  IF N_ELEMENTS(average_mode) GT 0 THEN BEGIN
    IF STRLOWCASE(average_mode) EQ 'median' THEN use_median = 1
  ENDIF

  ; 3. Determine Mode (Wavelengths vs Indices) and Find Indices
  IF N_ELEMENTS(wvl) GT 0 THEN BEGIN
    ; --- WAVELENGTH MODE ---
    IF N_ELEMENTS(wvl) NE nb THEN BEGIN
      MESSAGE, 'Wavelength array size does not match cube bands.', /CONTINUE
      RETURN, -1
    ENDIF

    ; Find closest band index for the requested wavelengths
    dummy = MIN(ABS(wvl - in_left), idx_left)
    dummy = MIN(ABS(wvl - in_right), idx_right)

    ; Define X-axis values for slope calculation using ACTUAL wavelengths
    x_L = float(wvl[idx_left])
    x_R = float(wvl[idx_right])

    PRINT, '--- Auto-Selected Bands (Slope) ---'

    ; --- LEFT PRINTOUT ---
    s = (idx_left - (w_L / 2)) > 0
    e = (s + w_L - 1) < (nb - 1)
    vals = wvl[s:e]
    v_str = STRJOIN(STRTRIM(vals, 2), ', ')
    PRINT, 'Left:   Req=', STRTRIM(in_left, 2), ', Width=', STRTRIM(w_L, 2), ', Actual=', v_str

    ; --- RIGHT PRINTOUT ---
    s = (idx_right - (w_R / 2)) > 0
    e = (s + w_R - 1) < (nb - 1)
    vals = wvl[s:e]
    v_str = STRJOIN(STRTRIM(vals, 2), ', ')
    PRINT, 'Right:  Req=', STRTRIM(in_right, 2), ', Width=', STRTRIM(w_R, 2), ', Actual=', v_str

  ENDIF ELSE BEGIN
    ; --- INDEX MODE ---
    idx_left  = FIX(in_left)
    idx_right = FIX(in_right)

    ; Define X-axis values (just indices)
    x_L = float(idx_left)
    x_R = float(idx_right)
  ENDELSE

  ; 4. Error Checking: Bounds & Logic

  ; A. Check if indices are inside cube dimensions
  IF (idx_left LT 0) OR (idx_left GE nb) OR $
    (idx_right LT 0) OR (idx_right GE nb) THEN BEGIN
    MESSAGE, 'Calculated or provided band indices are out of cube range.', /CONTINUE
    RETURN, -1
  ENDIF

  ; B. Check if Left and Right are identical (Slope Undefined)
  IF (x_L EQ x_R) THEN BEGIN
    MESSAGE, 'ERROR: Left and Right positions cannot be identical.', /CONTINUE
    RETURN, -1
  ENDIF

  ; Note: Unlike Band Depth, Slope doesn't strictly require Left < Right.
  ; A negative dx is mathematically valid, though usually physically inverted.
  ; We proceed but dx will handle the sign.

  ; 5. Extract 2D Slices (With Width Averaging)

  ; --- LEFT SHOULDER EXTRACTION ---
  s = (idx_left - (w_L / 2)) > 0
  e = (s + w_L - 1) < (nb - 1)

  IF s EQ e THEN val_L = cube[*,*,s] $
  ELSE BEGIN
    IF use_median THEN val_L = MEDIAN(cube[*,*,s:e], DIMENSION=3) $
    ELSE val_L = TOTAL(cube[*,*,s:e], 3) / FLOAT(e - s + 1)
  ENDELSE

  ; --- RIGHT SHOULDER EXTRACTION ---
  s = (idx_right - (w_R / 2)) > 0
  e = (s + w_R - 1) < (nb - 1)

  IF s EQ e THEN val_R = cube[*,*,s] $
  ELSE BEGIN
    IF use_median THEN val_R = MEDIAN(cube[*,*,s:e], DIMENSION=3) $
    ELSE val_R = TOTAL(cube[*,*,s:e], 3) / FLOAT(e - s + 1)
  ENDELSE

  ; 6. Calculate Spectral Slope
  ; Formula: m = (y2 - y1) / (x2 - x1)

  dx = x_R - x_L

  ; Initialize Output
  slope_map = MAKE_ARRAY(nx, ny, /FLOAT, VALUE=0.0)

  ; Calculate
  ; Note: val_R and val_L are 2D arrays, dx is a scalar.
  slope_map = (val_R - val_L) / dx

  ; Save as .csv file
  filename='Saved data paper2/'+id+'_'+name+'.csv'
  write_csv,filename,slope_map
  print,'File written: '+filename

  RETURN, slope_map

END

FUNCTION initial_ROI_pass,phyllosilicate_int_array,sr_qub,id,threshold,local_band_depth,masked_kaolin_param,ferrous_slope,ROI_type=ROI_type

  ;Written by Thanos Klidaras on 10/31/2025 to allow different inputs (e.g., band depth, hematite parameter) for ROI generation
  ;Takes in:
  ;ROI type: 'gt_mean_BD','BD860_2','high_carbonate','high_sulfate'
  ;Outputs: variable originally called Greater_mean_BD_ROIs
  ;Call like: greater_mean_BD_ROIs=initial_ROI_pass(phyllosilicate_int_array,sr_1,id='1',ROI_type='GtMeanBD',0,fullrange_localBD)
  ;Updated 11/12/2025: changed variable names to be agnostic (e.g., variable BD860_2 -> index)
  ;Updated 1/12/2026: added ROIs for kaolin param, ferrous slope, BD1750_2, BD2100_2, BD2250, BD2265

  ;For each phyllosilicate phase, creates a mask of the pixels with the highest BD (An ROI for each to extract average spectra from)

  s=size(phyllosilicate_int_array[*,*],/dim)
  ROIs_for_each_phase=intarr(s[0],s[1],7)

  IF (n_elements(ROI_type) gt 0) AND (ROI_type EQ 'HighBandDepth') THEN BEGIN
    print,'Generating *high band depth* ROIs for each phyllosillicate phase...'
    print,'Threshold: local_band_depth gt mean + 1 std'

    phase1_mask=intarr(s[0],s[1])
    phase1_mask[where(phyllosilicate_int_array eq 1)]=1
    local_BD_phase1=local_band_depth*phase1_mask
    threshold=mean(local_BD_phase1[where(local_BD_phase1 ne 0.)]) + 1*stddev(local_BD_phase1[where(local_BD_phase1 ne 0.)])  ;threshold=mean+1*std
    phase1_ROI=intarr(s[0],s[1])
    phase1_ROI[where(local_BD_phase1 gt threshold)]=1
    ROIs_for_each_phase[*,*,1]=phase1_ROI

    phase2_mask=intarr(s[0],s[1])
    phase2_mask[where(phyllosilicate_int_array eq 2)]=1
    local_BD_phase2=local_band_depth*phase2_mask
    threshold=mean(local_BD_phase2[where(local_BD_phase2 ne 0.)]) + 1*stddev(local_BD_phase2[where(local_BD_phase2 ne 0.)])  ;threshold=mean+1*std
    phase2_ROI=intarr(s[0],s[1])
    phase2_ROI[where(local_BD_phase2 gt threshold)]=1
    ROIs_for_each_phase[*,*,2]=phase2_ROI

    phase3_mask=intarr(s[0],s[1])
    phase3_mask[where(phyllosilicate_int_array eq 3)]=1
    local_BD_phase3=local_band_depth*phase3_mask
    threshold=mean(local_BD_phase3[where(local_BD_phase3 ne 0.)]) + 1*stddev(local_BD_phase3[where(local_BD_phase3 ne 0.)])  ;threshold=mean+1*std
    phase3_ROI=intarr(s[0],s[1])
    phase3_ROI[where(local_BD_phase3 gt threshold)]=1
    ROIs_for_each_phase[*,*,3]=phase3_ROI

    phase4_mask=intarr(s[0],s[1])
    phase4_mask[where(phyllosilicate_int_array eq 4)]=1
    local_BD_phase4=local_band_depth*phase4_mask
    threshold=mean(local_BD_phase4[where(local_BD_phase4 ne 0.)]) + 1*stddev(local_BD_phase4[where(local_BD_phase4 ne 0.)])  ;threshold=mean+1*std
    phase4_ROI=intarr(s[0],s[1])
    phase4_ROI[where(local_BD_phase4 gt threshold)]=1
    ROIs_for_each_phase[*,*,4]=phase4_ROI

    phase5_mask=intarr(s[0],s[1])
    phase5_mask[where(phyllosilicate_int_array eq 5)]=1
    local_BD_phase5=local_band_depth*phase5_mask
    threshold=mean(local_BD_phase5[where(local_BD_phase5 ne 0.)]) + 1*stddev(local_BD_phase5[where(local_BD_phase5 ne 0.)])  ;threshold=mean+1*std
    phase5_ROI=intarr(s[0],s[1])
    phase5_ROI[where(local_BD_phase5 gt threshold)]=1
    ROIs_for_each_phase[*,*,5]=phase5_ROI

    phase6_mask=intarr(s[0],s[1])
    phase6_mask[where(phyllosilicate_int_array eq 6)]=1
    local_BD_phase6=local_band_depth*phase6_mask
    threshold=mean(local_BD_phase6[where(local_BD_phase6 ne 0.)]) + 1*stddev(local_BD_phase6[where(local_BD_phase6 ne 0.)])  ;threshold=mean+1*std
    phase6_ROI=intarr(s[0],s[1])
    phase6_ROI[where(local_BD_phase6 gt threshold)]=1
    ROIs_for_each_phase[*,*,6]=phase6_ROI

    all_phases_mask=intarr(s[0],s[1])
    all_phases_mask[where(phyllosilicate_int_array ne 0)]=1
    local_BD_all_phases=local_band_depth*all_phases_mask
    threshold=mean(local_BD_all_phases[where(local_BD_all_phases ne 0.)]) + 1*stddev(local_BD_all_phases[where(local_BD_all_phases ne 0.)])  ;threshold=mean+1*std
    all_phases_ROI=intarr(s[0],s[1])
    all_phases_ROI[where(local_BD_all_phases gt threshold)]=1
    ROIs_for_each_phase[*,*,0]=all_phases_ROI

  ENDIF

  IF (n_elements(ROI_type) gt 0) AND (ROI_type EQ 'BD860_2') THEN BEGIN
    print,'Generating *high BD_860 (hematite)* ROIs for each phyllosillicate phase...'
    print,'Threshold: BD860_2 gt ' +string(threshold)
    ;Isolating hematite refined spectral summary parameter
    index=sr_qub[*,*,6] ;BD860_2

    phase1_mask=intarr(s[0],s[1])
    phase1_mask[where(phyllosilicate_int_array eq 1)]=1
    high_index_phase1=index*phase1_mask
    phase1_ROI=intarr(s[0],s[1])
    phase1_ROI[where(high_index_phase1 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,1]=phase1_ROI

    phase2_mask=intarr(s[0],s[1])
    phase2_mask[where(phyllosilicate_int_array eq 2)]=1
    high_index_phase2=index*phase2_mask
    phase2_ROI=intarr(s[0],s[1])
    phase2_ROI[where(high_index_phase2 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,2]=phase2_ROI

    phase3_mask=intarr(s[0],s[1])
    phase3_mask[where(phyllosilicate_int_array eq 3)]=1
    high_index_phase3=index*phase3_mask
    phase3_ROI=intarr(s[0],s[1])
    phase3_ROI[where(high_index_phase3 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,3]=phase3_ROI

    phase4_mask=intarr(s[0],s[1])
    phase4_mask[where(phyllosilicate_int_array eq 4)]=1
    high_index_phase4=index*phase4_mask
    phase4_ROI=intarr(s[0],s[1])
    phase4_ROI[where(high_index_phase4 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,4]=phase4_ROI

    phase5_mask=intarr(s[0],s[1])
    phase5_mask[where(phyllosilicate_int_array eq 5)]=1
    high_index_phase5=index*phase5_mask
    phase5_ROI=intarr(s[0],s[1])
    phase5_ROI[where(high_index_phase5 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,5]=phase5_ROI

    phase6_mask=intarr(s[0],s[1])
    phase6_mask[where(phyllosilicate_int_array eq 6)]=1
    high_index_phase6=index*phase6_mask
    phase6_ROI=intarr(s[0],s[1])
    phase6_ROI[where(high_index_phase6 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,6]=phase6_ROI

    all_phases_mask=intarr(s[0],s[1])
    all_phases_mask[where(phyllosilicate_int_array ne 0)]=1
    all_phases_high_index=index*all_phases_mask
    all_phases_ROI=intarr(s[0],s[1])
    all_phases_ROI[where(all_phases_high_index gt threshold)]=1
    ROIs_for_each_phase[*,*,0]=all_phases_ROI

  ENDIF

  IF (n_elements(ROI_type) gt 0) AND (ROI_type EQ 'BD2500_2') THEN BEGIN
    print,'Generating *high BD2500_2 (carbonate)* ROIs for each phyllosillicate phase...'
    print,'Threshold: BD2500_2 gt ' +string(threshold)
    ;Isolating carbonate refined spectral summary parameter
    index=sr_qub[*,*,45] ;BD2500_2

    phase1_mask=intarr(s[0],s[1])
    phase1_mask[where(phyllosilicate_int_array eq 1)]=1
    high_index_phase1=index*phase1_mask
    phase1_ROI=intarr(s[0],s[1])
    phase1_ROI[where(high_index_phase1 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,1]=phase1_ROI

    phase2_mask=intarr(s[0],s[1])
    phase2_mask[where(phyllosilicate_int_array eq 2)]=1
    high_index_phase2=index*phase2_mask
    phase2_ROI=intarr(s[0],s[1])
    phase2_ROI[where(high_index_phase2 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,2]=phase2_ROI

    phase3_mask=intarr(s[0],s[1])
    phase3_mask[where(phyllosilicate_int_array eq 3)]=1
    high_index_phase3=index*phase3_mask
    phase3_ROI=intarr(s[0],s[1])
    phase3_ROI[where(high_index_phase3 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,3]=phase3_ROI

    phase4_mask=intarr(s[0],s[1])
    phase4_mask[where(phyllosilicate_int_array eq 4)]=1
    high_index_phase4=index*phase4_mask
    phase4_ROI=intarr(s[0],s[1])
    phase4_ROI[where(high_index_phase4 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,4]=phase4_ROI

    phase5_mask=intarr(s[0],s[1])
    phase5_mask[where(phyllosilicate_int_array eq 5)]=1
    high_index_phase5=index*phase5_mask
    phase5_ROI=intarr(s[0],s[1])
    phase5_ROI[where(high_index_phase5 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,5]=phase5_ROI

    phase6_mask=intarr(s[0],s[1])
    phase6_mask[where(phyllosilicate_int_array eq 6)]=1
    high_index_phase6=index*phase6_mask
    phase6_ROI=intarr(s[0],s[1])
    phase6_ROI[where(high_index_phase6 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,6]=phase6_ROI

    all_phases_mask=intarr(s[0],s[1])
    all_phases_mask[where(phyllosilicate_int_array ne 0)]=1
    all_phases_high_index=index*all_phases_mask
    all_phases_ROI=intarr(s[0],s[1])
    all_phases_ROI[where(all_phases_high_index gt threshold)]=1
    ROIs_for_each_phase[*,*,0]=all_phases_ROI

  ENDIF

  IF (n_elements(ROI_type) gt 0) AND (ROI_type EQ 'CINDEX2') THEN BEGIN
    print,'Generating *high CINDEX2 (carbonate)* ROIs for each phyllosillicate phase...'
    print,'Threshold: CINDEX2 gt ' +string(threshold)
    ;Isolating carbonate refined spectral summary parameter
    index=sr_qub[*,*,50] ;CINDEX2

    phase1_mask=intarr(s[0],s[1])
    phase1_mask[where(phyllosilicate_int_array eq 1)]=1
    high_index_phase1=index*phase1_mask
    phase1_ROI=intarr(s[0],s[1])
    phase1_ROI[where(high_index_phase1 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,1]=phase1_ROI

    phase2_mask=intarr(s[0],s[1])
    phase2_mask[where(phyllosilicate_int_array eq 2)]=1
    high_index_phase2=index*phase2_mask
    phase2_ROI=intarr(s[0],s[1])
    phase2_ROI[where(high_index_phase2 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,2]=phase2_ROI

    phase3_mask=intarr(s[0],s[1])
    phase3_mask[where(phyllosilicate_int_array eq 3)]=1
    high_index_phase3=index*phase3_mask
    phase3_ROI=intarr(s[0],s[1])
    phase3_ROI[where(high_index_phase3 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,3]=phase3_ROI

    phase4_mask=intarr(s[0],s[1])
    phase4_mask[where(phyllosilicate_int_array eq 4)]=1
    high_index_phase4=index*phase4_mask
    phase4_ROI=intarr(s[0],s[1])
    phase4_ROI[where(high_index_phase4 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,4]=phase4_ROI

    phase5_mask=intarr(s[0],s[1])
    phase5_mask[where(phyllosilicate_int_array eq 5)]=1
    high_index_phase5=index*phase5_mask
    phase5_ROI=intarr(s[0],s[1])
    phase5_ROI[where(high_index_phase5 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,5]=phase5_ROI

    phase6_mask=intarr(s[0],s[1])
    phase6_mask[where(phyllosilicate_int_array eq 6)]=1
    high_index_phase6=index*phase6_mask
    phase6_ROI=intarr(s[0],s[1])
    phase6_ROI[where(high_index_phase6 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,6]=phase6_ROI

    all_phases_mask=intarr(s[0],s[1])
    all_phases_mask[where(phyllosilicate_int_array ne 0)]=1
    all_phases_high_index=index*all_phases_mask
    all_phases_ROI=intarr(s[0],s[1])
    all_phases_ROI[where(all_phases_high_index gt threshold)]=1
    ROIs_for_each_phase[*,*,0]=all_phases_ROI

  ENDIF

  IF (n_elements(ROI_type) gt 0) AND (ROI_type EQ 'SINDEX2') THEN BEGIN
    print,'Generating *high SINDEX2 (sulfate)* ROIs for each phyllosillicate phase...'
    print,'Threshold: SINDEX2 gt ' +string(threshold)
    ;Isolating sulfate refined spectral summary parameter
    index=sr_qub[*,*,41] ;SINDEX2

    phase1_mask=intarr(s[0],s[1])
    phase1_mask[where(phyllosilicate_int_array eq 1)]=1
    high_index_phase1=index*phase1_mask
    phase1_ROI=intarr(s[0],s[1])
    phase1_ROI[where(high_index_phase1 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,1]=phase1_ROI

    phase2_mask=intarr(s[0],s[1])
    phase2_mask[where(phyllosilicate_int_array eq 2)]=1
    high_index_phase2=index*phase2_mask
    phase2_ROI=intarr(s[0],s[1])
    phase2_ROI[where(high_index_phase2 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,2]=phase2_ROI

    phase3_mask=intarr(s[0],s[1])
    phase3_mask[where(phyllosilicate_int_array eq 3)]=1
    high_index_phase3=index*phase3_mask
    phase3_ROI=intarr(s[0],s[1])
    phase3_ROI[where(high_index_phase3 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,3]=phase3_ROI

    phase4_mask=intarr(s[0],s[1])
    phase4_mask[where(phyllosilicate_int_array eq 4)]=1
    high_index_phase4=index*phase4_mask
    phase4_ROI=intarr(s[0],s[1])
    phase4_ROI[where(high_index_phase4 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,4]=phase4_ROI

    phase5_mask=intarr(s[0],s[1])
    phase5_mask[where(phyllosilicate_int_array eq 5)]=1
    high_index_phase5=index*phase5_mask
    phase5_ROI=intarr(s[0],s[1])
    phase5_ROI[where(high_index_phase5 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,5]=phase5_ROI

    phase6_mask=intarr(s[0],s[1])
    phase6_mask[where(phyllosilicate_int_array eq 6)]=1
    high_index_phase6=index*phase6_mask
    phase6_ROI=intarr(s[0],s[1])
    phase6_ROI[where(high_index_phase6 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,6]=phase6_ROI

    all_phases_mask=intarr(s[0],s[1])
    all_phases_mask[where(phyllosilicate_int_array ne 0)]=1
    all_phases_high_index=index*all_phases_mask
    all_phases_ROI=intarr(s[0],s[1])
    all_phases_ROI[where(all_phases_high_index gt threshold)]=1
    ROIs_for_each_phase[*,*,0]=all_phases_ROI

  ENDIF

  IF (n_elements(ROI_type) gt 0) AND (ROI_type EQ 'Kaolin') THEN BEGIN
    print,'Generating *kaolinite/halloysite* ROIs for each phyllosillicate phase...'
    print,'Threshold: kaolin_param eq 1'
    ;Isolating kaolin parameter
    index=masked_kaolin_param ;Kaolin parameter, masked to Al-unit only

    phase1_mask=intarr(s[0],s[1])
    phase1_mask[where(phyllosilicate_int_array eq 1)]=1
    high_index_phase1=index*phase1_mask
    phase1_ROI=intarr(s[0],s[1])
    phase1_ROI[where(high_index_phase1 eq threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,1]=phase1_ROI

    phase2_mask=intarr(s[0],s[1])
    phase2_mask[where(phyllosilicate_int_array eq 2)]=1
    high_index_phase2=index*phase2_mask
    phase2_ROI=intarr(s[0],s[1])
    phase2_ROI[where(high_index_phase2 eq threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,2]=phase2_ROI

    phase3_mask=intarr(s[0],s[1])
    phase3_mask[where(phyllosilicate_int_array eq 3)]=1
    high_index_phase3=index*phase3_mask
    phase3_ROI=intarr(s[0],s[1])
    phase3_ROI[where(high_index_phase3 eq threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,3]=phase3_ROI

    phase4_mask=intarr(s[0],s[1])
    phase4_mask[where(phyllosilicate_int_array eq 4)]=1
    high_index_phase4=index*phase4_mask
    phase4_ROI=intarr(s[0],s[1])
    phase4_ROI[where(high_index_phase4 eq threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,4]=phase4_ROI

    phase5_mask=intarr(s[0],s[1])
    phase5_mask[where(phyllosilicate_int_array eq 5)]=1
    high_index_phase5=index*phase5_mask
    phase5_ROI=intarr(s[0],s[1])
    phase5_ROI[where(high_index_phase5 eq threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,5]=phase5_ROI

    phase6_mask=intarr(s[0],s[1])
    phase6_mask[where(phyllosilicate_int_array eq 6)]=1
    high_index_phase6=index*phase6_mask
    phase6_ROI=intarr(s[0],s[1])
    phase6_ROI[where(high_index_phase6 eq threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,6]=phase6_ROI

    all_phases_mask=intarr(s[0],s[1])
    all_phases_mask[where(phyllosilicate_int_array ne 0)]=1
    all_phases_high_index=index*all_phases_mask
    all_phases_ROI=intarr(s[0],s[1])
    all_phases_ROI[where(all_phases_high_index eq threshold)]=1
    ROIs_for_each_phase[*,*,0]=all_phases_ROI

  ENDIF

  IF (n_elements(ROI_type) gt 0) AND (ROI_type EQ 'FerrousSlope') THEN BEGIN
    print,'Generating *high ferrous slope* ROIs for each phyllosillicate phase...'
    print,'Threshold: ferrous_slope_Klidaras gt ' +string(threshold)
    ;Isolating ferrous_slope_Klidaras spectral parameter
    index=ferrous_slope ;1.2-1.8 micron spectral slope

    phase1_mask=intarr(s[0],s[1])
    phase1_mask[where(phyllosilicate_int_array eq 1)]=1
    high_index_phase1=index*phase1_mask
    phase1_ROI=intarr(s[0],s[1])
    phase1_ROI[where(high_index_phase1 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,1]=phase1_ROI

    phase2_mask=intarr(s[0],s[1])
    phase2_mask[where(phyllosilicate_int_array eq 2)]=1
    high_index_phase2=index*phase2_mask
    phase2_ROI=intarr(s[0],s[1])
    phase2_ROI[where(high_index_phase2 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,2]=phase2_ROI

    phase3_mask=intarr(s[0],s[1])
    phase3_mask[where(phyllosilicate_int_array eq 3)]=1
    high_index_phase3=index*phase3_mask
    phase3_ROI=intarr(s[0],s[1])
    phase3_ROI[where(high_index_phase3 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,3]=phase3_ROI

    phase4_mask=intarr(s[0],s[1])
    phase4_mask[where(phyllosilicate_int_array eq 4)]=1
    high_index_phase4=index*phase4_mask
    phase4_ROI=intarr(s[0],s[1])
    phase4_ROI[where(high_index_phase4 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,4]=phase4_ROI

    phase5_mask=intarr(s[0],s[1])
    phase5_mask[where(phyllosilicate_int_array eq 5)]=1
    high_index_phase5=index*phase5_mask
    phase5_ROI=intarr(s[0],s[1])
    phase5_ROI[where(high_index_phase5 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,5]=phase5_ROI

    phase6_mask=intarr(s[0],s[1])
    phase6_mask[where(phyllosilicate_int_array eq 6)]=1
    high_index_phase6=index*phase6_mask
    phase6_ROI=intarr(s[0],s[1])
    phase6_ROI[where(high_index_phase6 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,6]=phase6_ROI

    all_phases_mask=intarr(s[0],s[1])
    all_phases_mask[where(phyllosilicate_int_array ne 0)]=1
    all_phases_high_index=index*all_phases_mask
    all_phases_ROI=intarr(s[0],s[1])
    all_phases_ROI[where(all_phases_high_index gt threshold)]=1
    ROIs_for_each_phase[*,*,0]=all_phases_ROI

  ENDIF

  IF (n_elements(ROI_type) gt 0) AND (ROI_type EQ 'BD2100_2') THEN BEGIN
    print,'Generating *high BD2100_2 (monohydrated sulfates, alunite)* ROIs for each phyllosillicate phase...'
    print,'Threshold: BD2100_2 gt ' +string(threshold)
    index=sr_qub[*,*,28] ;BD2100_2

    phase1_mask=intarr(s[0],s[1])
    phase1_mask[where(phyllosilicate_int_array eq 1)]=1
    high_index_phase1=index*phase1_mask
    phase1_ROI=intarr(s[0],s[1])
    phase1_ROI[where(high_index_phase1 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,1]=phase1_ROI

    phase2_mask=intarr(s[0],s[1])
    phase2_mask[where(phyllosilicate_int_array eq 2)]=1
    high_index_phase2=index*phase2_mask
    phase2_ROI=intarr(s[0],s[1])
    phase2_ROI[where(high_index_phase2 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,2]=phase2_ROI

    phase3_mask=intarr(s[0],s[1])
    phase3_mask[where(phyllosilicate_int_array eq 3)]=1
    high_index_phase3=index*phase3_mask
    phase3_ROI=intarr(s[0],s[1])
    phase3_ROI[where(high_index_phase3 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,3]=phase3_ROI

    phase4_mask=intarr(s[0],s[1])
    phase4_mask[where(phyllosilicate_int_array eq 4)]=1
    high_index_phase4=index*phase4_mask
    phase4_ROI=intarr(s[0],s[1])
    phase4_ROI[where(high_index_phase4 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,4]=phase4_ROI

    phase5_mask=intarr(s[0],s[1])
    phase5_mask[where(phyllosilicate_int_array eq 5)]=1
    high_index_phase5=index*phase5_mask
    phase5_ROI=intarr(s[0],s[1])
    phase5_ROI[where(high_index_phase5 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,5]=phase5_ROI

    phase6_mask=intarr(s[0],s[1])
    phase6_mask[where(phyllosilicate_int_array eq 6)]=1
    high_index_phase6=index*phase6_mask
    phase6_ROI=intarr(s[0],s[1])
    phase6_ROI[where(high_index_phase6 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,6]=phase6_ROI

    all_phases_mask=intarr(s[0],s[1])
    all_phases_mask[where(phyllosilicate_int_array ne 0)]=1
    all_phases_high_index=index*all_phases_mask
    all_phases_ROI=intarr(s[0],s[1])
    all_phases_ROI[where(all_phases_high_index gt threshold)]=1
    ROIs_for_each_phase[*,*,0]=all_phases_ROI

  ENDIF

  IF (n_elements(ROI_type) gt 0) AND (ROI_type EQ 'BD2250') THEN BEGIN
    print,'Generating *high BD2250 (Opal)* ROIs for each phyllosillicate phase...'
    print,'Threshold: BD2250 gt ' +string(threshold)
    index=sr_qub[*,*,35] ;BD2250

    phase1_mask=intarr(s[0],s[1])
    phase1_mask[where(phyllosilicate_int_array eq 1)]=1
    high_index_phase1=index*phase1_mask
    phase1_ROI=intarr(s[0],s[1])
    phase1_ROI[where(high_index_phase1 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,1]=phase1_ROI

    phase2_mask=intarr(s[0],s[1])
    phase2_mask[where(phyllosilicate_int_array eq 2)]=1
    high_index_phase2=index*phase2_mask
    phase2_ROI=intarr(s[0],s[1])
    phase2_ROI[where(high_index_phase2 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,2]=phase2_ROI

    phase3_mask=intarr(s[0],s[1])
    phase3_mask[where(phyllosilicate_int_array eq 3)]=1
    high_index_phase3=index*phase3_mask
    phase3_ROI=intarr(s[0],s[1])
    phase3_ROI[where(high_index_phase3 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,3]=phase3_ROI

    phase4_mask=intarr(s[0],s[1])
    phase4_mask[where(phyllosilicate_int_array eq 4)]=1
    high_index_phase4=index*phase4_mask
    phase4_ROI=intarr(s[0],s[1])
    phase4_ROI[where(high_index_phase4 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,4]=phase4_ROI

    phase5_mask=intarr(s[0],s[1])
    phase5_mask[where(phyllosilicate_int_array eq 5)]=1
    high_index_phase5=index*phase5_mask
    phase5_ROI=intarr(s[0],s[1])
    phase5_ROI[where(high_index_phase5 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,5]=phase5_ROI

    phase6_mask=intarr(s[0],s[1])
    phase6_mask[where(phyllosilicate_int_array eq 6)]=1
    high_index_phase6=index*phase6_mask
    phase6_ROI=intarr(s[0],s[1])
    phase6_ROI[where(high_index_phase6 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,6]=phase6_ROI

    all_phases_mask=intarr(s[0],s[1])
    all_phases_mask[where(phyllosilicate_int_array ne 0)]=1
    all_phases_high_index=index*all_phases_mask
    all_phases_ROI=intarr(s[0],s[1])
    all_phases_ROI[where(all_phases_high_index gt threshold)]=1
    ROIs_for_each_phase[*,*,0]=all_phases_ROI

  ENDIF

  IF (n_elements(ROI_type) gt 0) AND (ROI_type EQ 'BD2265') THEN BEGIN
    print,'Generating *high BD2265 (Jarosite, acid-leached nontronite)* ROIs for each phyllosillicate phase...'
    print,'Threshold: BD2265 gt ' +string(threshold)
    index=sr_qub[*,*,37] ;BD2265

    phase1_mask=intarr(s[0],s[1])
    phase1_mask[where(phyllosilicate_int_array eq 1)]=1
    high_index_phase1=index*phase1_mask
    phase1_ROI=intarr(s[0],s[1])
    phase1_ROI[where(high_index_phase1 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,1]=phase1_ROI

    phase2_mask=intarr(s[0],s[1])
    phase2_mask[where(phyllosilicate_int_array eq 2)]=1
    high_index_phase2=index*phase2_mask
    phase2_ROI=intarr(s[0],s[1])
    phase2_ROI[where(high_index_phase2 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,2]=phase2_ROI

    phase3_mask=intarr(s[0],s[1])
    phase3_mask[where(phyllosilicate_int_array eq 3)]=1
    high_index_phase3=index*phase3_mask
    phase3_ROI=intarr(s[0],s[1])
    phase3_ROI[where(high_index_phase3 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,3]=phase3_ROI

    phase4_mask=intarr(s[0],s[1])
    phase4_mask[where(phyllosilicate_int_array eq 4)]=1
    high_index_phase4=index*phase4_mask
    phase4_ROI=intarr(s[0],s[1])
    phase4_ROI[where(high_index_phase4 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,4]=phase4_ROI

    phase5_mask=intarr(s[0],s[1])
    phase5_mask[where(phyllosilicate_int_array eq 5)]=1
    high_index_phase5=index*phase5_mask
    phase5_ROI=intarr(s[0],s[1])
    phase5_ROI[where(high_index_phase5 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,5]=phase5_ROI

    phase6_mask=intarr(s[0],s[1])
    phase6_mask[where(phyllosilicate_int_array eq 6)]=1
    high_index_phase6=index*phase6_mask
    phase6_ROI=intarr(s[0],s[1])
    phase6_ROI[where(high_index_phase6 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,6]=phase6_ROI

    all_phases_mask=intarr(s[0],s[1])
    all_phases_mask[where(phyllosilicate_int_array ne 0)]=1
    all_phases_high_index=index*all_phases_mask
    all_phases_ROI=intarr(s[0],s[1])
    all_phases_ROI[where(all_phases_high_index gt threshold)]=1
    ROIs_for_each_phase[*,*,0]=all_phases_ROI

  ENDIF

  IF (n_elements(ROI_type) gt 0) AND (ROI_type EQ 'BD1750_2') THEN BEGIN
    print,'Generating *high BD1750_2 (Gypsum, alunite)* ROIs for each phyllosillicate phase...'
    print,'Threshold: BD1750_2 gt ' +string(threshold)
    index=sr_qub[*,*,24] ;BD1750_2

    phase1_mask=intarr(s[0],s[1])
    phase1_mask[where(phyllosilicate_int_array eq 1)]=1
    high_index_phase1=index*phase1_mask
    phase1_ROI=intarr(s[0],s[1])
    phase1_ROI[where(high_index_phase1 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,1]=phase1_ROI

    phase2_mask=intarr(s[0],s[1])
    phase2_mask[where(phyllosilicate_int_array eq 2)]=1
    high_index_phase2=index*phase2_mask
    phase2_ROI=intarr(s[0],s[1])
    phase2_ROI[where(high_index_phase2 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,2]=phase2_ROI

    phase3_mask=intarr(s[0],s[1])
    phase3_mask[where(phyllosilicate_int_array eq 3)]=1
    high_index_phase3=index*phase3_mask
    phase3_ROI=intarr(s[0],s[1])
    phase3_ROI[where(high_index_phase3 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,3]=phase3_ROI

    phase4_mask=intarr(s[0],s[1])
    phase4_mask[where(phyllosilicate_int_array eq 4)]=1
    high_index_phase4=index*phase4_mask
    phase4_ROI=intarr(s[0],s[1])
    phase4_ROI[where(high_index_phase4 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,4]=phase4_ROI

    phase5_mask=intarr(s[0],s[1])
    phase5_mask[where(phyllosilicate_int_array eq 5)]=1
    high_index_phase5=index*phase5_mask
    phase5_ROI=intarr(s[0],s[1])
    phase5_ROI[where(high_index_phase5 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,5]=phase5_ROI

    phase6_mask=intarr(s[0],s[1])
    phase6_mask[where(phyllosilicate_int_array eq 6)]=1
    high_index_phase6=index*phase6_mask
    phase6_ROI=intarr(s[0],s[1])
    phase6_ROI[where(high_index_phase6 gt threshold)]=1  ;Threshold = user-set
    ROIs_for_each_phase[*,*,6]=phase6_ROI

    all_phases_mask=intarr(s[0],s[1])
    all_phases_mask[where(phyllosilicate_int_array ne 0)]=1
    all_phases_high_index=index*all_phases_mask
    all_phases_ROI=intarr(s[0],s[1])
    all_phases_ROI[where(all_phases_high_index gt threshold)]=1
    ROIs_for_each_phase[*,*,0]=all_phases_ROI

  ENDIF

  return,ROIs_for_each_phase

END

FUNCTION ROI_valid_largest_cluster,ROIs,id,if_qub,crism_wvl_array,plot_limits,ROI_type=ROI_type
  ;PURPOSE: Takes in an array of ROIs, and isolates the largest contiguous cluster that is >15 pixels and has an aspect ratio of <2 (not tall)
  ;UPDATES:
  ;26/09/2025: added minimum ROI threshold (>15 pixels), added aspect ratio limits (<5)
  ;12/10/2025: Fixed crash on empty input groups, fixed writing empty CSVs, changed aspect ratio limit from (<5 to <3)
  ;            I now call save_ROI_spectrum here, to reduce clutter in MASTER_paper2
  ;12/23/2025: Changed aspect ratio limit from <3 to <2
  s=size(ROIs,/dim)
  cluster_ROIs=intarr(s[0],s[1],s[2])
  merged_ROIs=intarr(s[0],s[1])
  merged_cluster_ROIs=intarr(s[0],s[1])
  s2=size(ROIs[*,*,0],/dim)
  largest_cluster_mask=intarr(s2[0],s2[1])

  ;Identifying each array's largest valid cluster (>=15 pixels) with aspect ratio < 2
  for i = 0, 6 do begin
    largest_cluster_mask = INTARR(s2[0], s2[1])
    ROI      = ROIs[*,*,i]
    ROIlabel = LABEL_REGION(ROI)
    n        = MAX(ROIlabel)

    ; --- CRASH FIX: Check for empty group ---
    if n eq 0 then begin
      PRINT, 'No pixels found for phyllosilicate group ', i
      continue
    endif

    counts   = HISTOGRAM(ROIlabel, MIN=1, MAX=n, NBINS=n)

    ; Keep only clusters with size >= 15:
    allowed = WHERE(counts GE 15, n_allowed)

    if n_allowed gt 0 then begin
      counts_allowed = counts[allowed]

      ; Sort candidates by size (descending):
      ord_desc   = REVERSE(SORT(counts_allowed))
      cand_bins  = allowed[ord_desc]        ; 0-based bin indices
      max_ar     = 2.0
      found      = 0B

      dims = SIZE(ROIlabel, /DIM)           ; for bbox calcs once per slice

      for j = 0, N_ELEMENTS(cand_bins)-1 do begin
        label_val = cand_bins[j] + 1L       ; convert bin -> label (1..n)
        cluster_coords = WHERE(ROIlabel EQ label_val, ncoords)
        if ncoords le 0 then continue

        ;Aspect ratio from bounding box
        x = cluster_coords mod dims[0]
        y = cluster_coords /   dims[0]

        xmin = min(x, MAX=xmax)
        ymin = min(y, MAX=ymax)

        width  = (xmax - xmin + 1L)
        height = (ymax - ymin + 1L)

        ar = float(height) / float(width)   ; can be <1 (wide) or >1 (tall)

        if ar lt max_ar then begin
          largest_cluster_mask[cluster_coords] = 1
          PRINT, 'phyllosilicate group ', i, $
            ' chosen cluster label=', label_val, $
            ' size=', ncoords, ' aspect ratio=', ar
          found = 1B
          break
        endif else begin
          PRINT, 'phyllosilicate group ', i, $
            ' rejected cluster label=', label_val, $
            ' size=', ncoords, ' aspect ratio=', ar, ' (>=2)'
        endelse
      endfor

      if not found then $
        PRINT, 'No cluster >=15 px with aspect ratio < 2 for phyllosilicate group ', i

    endif else begin
      PRINT, 'No cluster >= 15 pixels for phyllosilicate group ', i
    endelse

    cluster_ROIs[*,*,i] = largest_cluster_mask
  endfor

  ;Creating merged ROI maps to enable easy evaluation of changes
  for i=1,6 do begin
    ; Check using count keyword to get true pixel count (avoids the -1/1 confusion)
    dummy = where(cluster_ROIs[*,*,i] eq 1, count)
    print,'Largest valid cluster ROI for phase ',i,' has:', count ,' pixels'

    merged_ROIs=merged_ROIs+(i*ROIs[*,*,i])
    merged_cluster_ROIs=merged_cluster_ROIs+(i*cluster_ROIs[*,*,i])
  endfor

  ;Saving largest cluster valid ROIs as .csv files
  ;Loop 1-6. Only write if mask contains valid pixels (MAX value is 1)

  for i=1, 6 do begin
    if MAX(cluster_ROIs[*,*,i]) EQ 1 then begin
      filename = 'Saved data paper2/' + id + '_Group' + STRTRIM(i, 2) + '_ROI_' + ROI_type + '.csv'
      write_csv, filename, cluster_ROIs[*,*,i]
      print, 'File written: ' + filename
    endif
  endfor

  ;Saving ROI spectra:

  x_short=plot_limits[0]
  x_long=plot_limits[1]
  y_low=plot_limits[2]
  y_high=plot_limits[3]

  ; Create a pointer array to store the return values (size 7 to allow indices 1-6)
  spectra_results = PTRARR(7)

  FOR i=1, 6 DO BEGIN
    ; Check if the current ROI slice has any valid pixels (MAX value will be 1)
    IF MAX(cluster_ROIs[*,*,i]) EQ 1 THEN BEGIN

      PRINT, 'Generating spectrum for Group ', i

      ; Call the function dynamically
      ; STRTRIM(i, 2) converts the integer i to a string ('1', '2', etc.)
      result = save_ROI_spectrum(if_qub, cluster_ROIs[*,*,i], crism_wvl_array, x_short, x_long, y_low, y_high, 1, id, STRTRIM(i, 2), ROI_type=ROI_type)

      ; Store the result in the pointer array
      spectra_results[i] = PTR_NEW(result)

    ENDIF ELSE BEGIN
      PRINT, 'Skipping spectrum generation for Group ', i, ' (Empty ROI)'
    ENDELSE

  ENDFOR

  return,cluster_ROIs

END

FUNCTION save_ROI_spectrum,cube,mask,wvl,lambda_short,lambda_long,y_short,y_long,do_oplot,id,Grp,ROI_type=ROI_type

  ;Written by Thanos Klidaras on 18/06/2024


  ;FUNCTION:
  ;          Plots the smoothed spectrum of a given pixel in a spectral cube, over a custom wavelength range
  ;          It can also plot the smoothed average spectrum across an ROI file that is passed to it
  ;          if do_oplot eq 0, the spectrum will be displayed in a window
  ;EXAMPLE OF USE:
  ;          plot_spectra(if_1,strict_al_mask,crism_wvl_array,1500,2400,0.12,0.16,0)



  ;spectra=smooth(reform(cube[i,j,*]),3,/nan,/edge_truncate) ;Functionality to extract spectrum from given pixel coordinates, not needed currently

  ;Counting number of pixels (value=1) in the ROI (mask array)
  pixelnum=n_elements(mask[where(mask eq 1)])

  ;Looping through every CRISM channel and masking pixels
  s=size(cube[*,*,*],/dim)
  average_spectrum=fltarr(s[2])
  standard_dev=fltarr(s[2])

  for i=0,s[2]-1 do begin
    masked_channel=reform(cube[*,*,i]*mask)
    average_spectrum[i]=mean(masked_channel[where((masked_channel ne 0) and (masked_channel lt 100))]) ;Calculates the average value of all the masked pixels for all 489 CRISM channels, skipping over pixels with no value or glitchy super high value,  and stores it as one spectrum
    standard_dev[i]=stddev(masked_channel[where((masked_channel ne 0) and (masked_channel lt 100))])
  endfor

  spectrum=smooth(average_spectrum,3,/nan,/edge_truncate)   ;Smooths the average spectrum
  unsmoothed_spectrum=average_spectrum
  ;Creating the spectral plot

  graph_wvl_range=[lambda_short,lambda_long]

  graph_channel_range=[closest(wvl,graph_wvl_range[0]),closest(wvl,graph_wvl_range[1])] ;The CRISM channels corresponding to the wavelength range of interest (e.g., 1.5 um - 2.4 um = 147 - 283)

  yrange=[min(spectrum[graph_channel_range[0]:graph_channel_range[1]]),max(spectrum[graph_channel_range[0]:graph_channel_range[1]])]

  ;yrange=[y_short,y_long]


  if do_oplot eq 0 then plot,wvl,spectrum,xrange=[graph_wvl_range[0],graph_wvl_range[1]],yrange=yrange,color=cgcolor("white"),xtitle='Wavelength (nm)',ytitle='I/F reflectance',xcharsize = 1.5, ycharsize = 1.5 else oplot,wvl,spectrum,color=cgcolor("white")

  ;Saving spectra as .csv files
  filename='Saved data paper2/'+id+'_Group'+Grp+'_spectrum_'+ROI_type+'.csv'
  write_csv,filename,wvl,unsmoothed_spectrum,standard_dev
  print,'File written: '+filename

  return,spectrum

END

FUNCTION save_spectral_parameters,sr_qub,id,if_qub
  ;Refined spectral summary parameters of interest:
  ;Updated 8/19/2025: Corrected code with right indexes, see Viviano_refined_spectral_summary_parameters_ORDER.xlsx spreadsheet
  ;Updated 10/1/2025: Added several new parameters for iteration4
  ;Updated 1/9/2026: Added BDI1000IR ferrous phase investigation
  ;Updated 1/12/2026: Added BD1750_2 to verify sulfates
  ;Updated 6/8/2026: Added BD530_2 to verify crystalline iron oxides, improved text of prints (iteration8)

  ;Saving spectral parameters of interest
  R770=sr_qub[*,*,0]
  BD860_2=sr_qub[*,*,6]
  BD1900_2=sr_qub[*,*,25]
  BD2265=sr_qub[*,*,37]
  BD2500_2=sr_qub[*,*,45]
  MIN2200=sr_qub[*,*,31]
  MIN2250=sr_qub[*,*,36]
  SINDEX2=sr_qub[*,*,41]
  CINDEX2=sr_qub[*,*,50]
  BD2210_2=sr_qub[*,*,32]
  BD2290=sr_qub[*,*,38]
  BD920_2=sr_qub[*,*,7]
  BD1300=sr_qub[*,*,15]
  OLINDEX3=sr_qub[*,*,13]
  LCPINDEX2=sr_qub[*,*,16]
  HCPINDEX2=sr_qub[*,*,17]
  BD2165=sr_qub[*,*,29]
  BD2190=sr_qub[*,*,30]
  D2300=sr_qub[*,*,39]
  BD1900r2=sr_qub[*,*,26]
  BD2250=sr_qub[*,*,35]
  BD2100_2=sr_qub[*,*,28]
  R443=if_qub[*,*,1]
  R527=if_qub[*,*,14]
  R599=if_qub[*,*,25]
  R801=if_qub[*,*,45]
  R912=if_qub[*,*,62]
  R1010=if_qub[*,*,77]
  BDI1000IR=sr_qub[*,*,12]
  BD1750_2=sr_qub[*,*,24]
  BD530_2=sr_qub[*,*,2]

  s=size(sr_qub,/dim)

  ;Writing each spectral parameter as a csv file
  print,'Saving relevant refined spectral summary parameters and I/F...'
  filename='Saved data paper2/'+id+'_R770.csv'
  write_csv,filename,R770
  print,'(1/31) File written: '+filename
  filename='Saved data paper2/'+id+'_BD860_2.csv'
  write_csv,filename,BD860_2
  print,'(2/31) File written: '+filename
  filename='Saved data paper2/'+id+'_BD1900_2.csv'
  write_csv,filename,BD1900_2
  print,'(3/31) File written: '+filename
  filename='Saved data paper2/'+id+'_BD2265.csv'
  write_csv,filename,BD2265
  print,'(4/31) File written: '+filename
  filename='Saved data paper2/'+id+'_BD2500_2.csv'
  write_csv,filename,BD2500_2
  print,'(5/31) File written: '+filename
  filename='Saved data paper2/'+id+'_MIN2200.csv'
  write_csv,filename,MIN2200
  print,'(6/31) File written: '+filename
  filename='Saved data paper2/'+id+'_MIN2250.csv'
  write_csv,filename,MIN2250
  print,'(7/31) File written: '+filename
  filename='Saved data paper2/'+id+'_SINDEX2.csv'
  write_csv,filename,SINDEX2
  print,'(8/31) File written: '+filename
  filename='Saved data paper2/'+id+'_CINDEX2.csv'
  write_csv,filename,CINDEX2
  print,'(9/31) File written: '+filename
  filename='Saved data paper2/'+id+'_BD2210_2.csv'
  write_csv,filename,BD2210_2
  print,'(10/31) File written: '+filename
  filename='Saved data paper2/'+id+'_BD2290.csv'
  write_csv,filename,BD2290
  print,'(11/31) File written: '+filename
  filename='Saved data paper2/'+id+'_BD920_2.csv'
  write_csv,filename,BD920_2
  print,'(12/31) File written: '+filename
  filename='Saved data paper2/'+id+'_BD1300.csv'
  write_csv,filename,BD1300
  print,'(13/31) File written: '+filename
  filename='Saved data paper2/'+id+'_OLINDEX3.csv'
  write_csv,filename,OLINDEX3
  print,'(14/31) File written: '+filename
  filename='Saved data paper2/'+id+'_LCPINDEX2.csv'
  write_csv,filename,LCPINDEX2
  print,'(15/31) File written: '+filename
  filename='Saved data paper2/'+id+'_HCPINDEX2.csv'
  write_csv,filename,HCPINDEX2
  print,'(16/31) File written: '+filename
  filename='Saved data paper2/'+id+'_BD2165.csv'
  write_csv,filename,BD2165
  print,'(17/31) File written: '+filename
  filename='Saved data paper2/'+id+'_BD2190.csv'
  write_csv,filename,BD2190
  print,'(18/31) File written: '+filename
  filename='Saved data paper2/'+id+'_D2300.csv'
  write_csv,filename,D2300
  print,'(19/31) File written: '+filename
  filename='Saved data paper2/'+id+'_BD1900r2.csv'
  write_csv,filename,BD1900r2
  print,'(20/31) File written: '+filename
  filename='Saved data paper2/'+id+'_BD2250.csv'
  write_csv,filename,BD2250
  print,'(21/31) File written: '+filename
  filename='Saved data paper2/'+id+'_BD2100_2.csv'
  write_csv,filename,BD2100_2
  print,'(22/31) File written: '+filename
  filename='Saved data paper2/'+id+'_R443.csv'
  write_csv,filename,R443
  print,'(23/31) File written: '+filename
  filename='Saved data paper2/'+id+'_R527.csv'
  write_csv,filename,R527
  print,'(24/31) File written: '+filename
  filename='Saved data paper2/'+id+'_R599.csv'
  write_csv,filename,R599
  print,'(25/31) File written: '+filename
  filename='Saved data paper2/'+id+'_R801.csv'
  write_csv,filename,R801
  print,'(26/31) File written: '+filename
  filename='Saved data paper2/'+id+'_R912.csv'
  write_csv,filename,R912
  print,'(27/31) File written: '+filename
  filename='Saved data paper2/'+id+'_R1010.csv'
  write_csv,filename,R1010
  print,'(28/31) File written: '+filename
  filename='Saved data paper2/'+id+'_BDI1000IR.csv'
  write_csv,filename,BDI1000IR
  print,'(29/31) File written: '+filename
  filename='Saved data paper2/'+id+'_BD1750_2.csv'
  write_csv,filename,BD1750_2
  print,'(30/31) File written: '+filename
  filename='Saved data paper2/'+id+'_BD530_2.csv'
  write_csv,filename,BD530_2
  print,'(31/31) File written: '+filename

  return,0.

END

FUNCTION generate_database
  ;Written by Chandler
  ;Updated by Thanos on Jul 1st 2025 to better comment code
  ;Updated by Thanos on Dec 23rd 2025 to update to handle the fewer cubes in iteration6 (only the 89/124 I validated had CCS)
  ;Updated by Thanos on Jun 8th 2026 to improve commenting, printing
  ;
  ;PURPOSE:
  ;
  ;By typing generate_database(), this code uses crism_load_data.pro to load the .if and .sr cubes
  ;for each CRISM image in sequence, and then runs master_paper2.pro on each to generate database
  ;Note: Running the entire database takes about 5 hours 10 mins!

  ;If you want to run on the entire spectral database, enter this into your console:
  ;
  ;         ENVI>generate_database()
  ;         (required functions: generate_database, crism_load_data)

  for i=0,88 do begin ;To run entire database, do: for i=0,88
    print,string(i+1),'/89'
    crism_load_data, rfl = rfl1, params=params1, loc = loc1,i

  endfor

END

PRO crism_load_data, rfl=rfl, params=params, loc=loc,value

  ;Chandler's code
  ;Thanos example: crism_load_data, rfl = rfl1, params=params1, loc = loc1,0
  ;Modified by Thanos Klidaras on Nov 20th 2024 so that file names are located automatically. The function generate_database() executes this program and loops through all 124 CRISM images :)
  ;Modified by Thanos Klidaras on Jul 1 2025:
  ;- Corrected file structure so that it navigates to the right folder
  ;- Implemented code so that the HRL cubes execute
  ;Modified by Thanos Klidaras on Dec 12 2025:
  ;- Updated filepath to new directory
  ;Modified by Thanos Klidaras on Dec 23 2025:
  ;- Added new filename_list (sixth_iteration_filename_list) reflecting only the cubes manually verified by me on 12.20.2025 to have CCS (72%)
  ;Modified by Thanos Klidaras on Dec 63 2025:
  ;- Variable sixth_iteration_filename_list was broken, made a fixed version that is also correctly numerically sorted
  ;How to save: save,filename='Saved cubes/1_Mawrth_Vallis_3bfb.sav',if_1,sr_1

  ;Read label file to determine array sizes (only lines and samples)

  ;Full filename list, all Carter et al. 2015 cubes validated by Ryan Navarre (n=124)
  filename_list=['100_None_1957c', '30_None_1011a', '101_None_8e5d', '31_Mawrth_Vallis_13e49', '102_None_168c0', '32_Circum_Hellas_9100', '103_None_7af6', '33_NW_Arabia_Terra_8438', '104_None_b278', '36_Mawrth_Vallis_129cf', '105_None_ca1b', '37_Mawrth_Vallis_12d3c', '106_None_bdd9', '38_None_137ad', '107_None_85d7', '39_Mawrth_Vallis_9326', '108_None_8ca3', '40_Terra_Cimmeria_Terra_Sirenum_62b6', '112_Terra_Cimmeria_Terra_Sirenum_c441', '41_NW_Arabia_Terra_a8ce', '114_Terra_Cimmeria_Terra_Sirenum_ac8b', '42_Circum_Hellas_147fb', '115_Terra_Cimmeria_Terra_Sirenum_11b46', '43_Circum_Hellas_96be', '116_Terra_Cimmeria_Terra_Sirenum_9bbc', '44_None_a865', '117_Terra_Cimmeria_Terra_Sirenum_96fe', '45_Circum_Hellas_16c3b', '118_Terra_Cimmeria_Terra_Sirenum_913a', '46_Mawrth_Vallis_bf57', '119_Terra_Cimmeria_Terra_Sirenum_9a7b', '48_None_1ec17', '11_Mawrth_Vallis_1176a', '50_NE_Syrtis_6408', '120_Terra_Cimmeria_Terra_Sirenum_951c', '53_Terra_Cimmeria_Terra_Sirenum_ada1', '121_Terra_Cimmeria_Terra_Sirenum_8a1c', '55_Terra_Cimmeria_Terra_Sirenum_92f7', '123_Terra_Cimmeria_Terra_Sirenum_97c5', '56_Mawrth_Vallis_1ec38', '124_Terra_Cimmeria_Terra_Sirenum_91eb', '57_Mawrth_Vallis_c467', '125_Terra_Cimmeria_Terra_Sirenum_10817', '58_None_d23e', '126_Terra_Cimmeria_Terra_Sirenum_18002', '59_None_18fcf', '127_Terra_Cimmeria_Terra_Sirenum_8377', '5_Mawrth_Vallis_89f7', '128_Terra_Cimmeria_Terra_Sirenum_167fd', '60_Mawrth_Vallis_bb21', '129_Terra_Cimmeria_Terra_Sirenum_7f7e', '61_Circum_Hellas_92b4', '130_Terra_Cimmeria_Terra_Sirenum_a106', '62_Circum_Hellas_8936', '131_Terra_Cimmeria_Terra_Sirenum_94c2', '63_Circum_Hellas_9665', '132_Terra_Cimmeria_Terra_Sirenum_1c624', '64_Circum_Hellas_a5a1', '133_Terra_Cimmeria_Terra_Sirenum_17199', '65_Circum_Hellas_a2b8', '134_Terra_Cimmeria_Terra_Sirenum_5120', '66_Circum_Hellas_136f4', '135_Circum_Hellas_9318', '68_Circum_Hellas_12455', '136_Circum_Hellas_98e8', '69_Circum_Hellas_14838', '137_Mawrth_Vallis_c066', '6_Mawrth_Vallis_b643', '138_Circum_Hellas_59df', '70_Circum_Hellas_920a', '139_Mawrth_Vallis_8f1a', '71_Circum_Hellas_14241', '13_NE_Syrtis_9d44', '73_Circum_Hellas_989c', '140_Mawrth_Vallis_a955', '74_Mawrth_Vallis_90af', '141_Terra_Cimmeria_Terra_Sirenum_9d7e', '75_Mawrth_Vallis_9cae', '143_Mawrth_Vallis_ac83', '76_Mawrth_Vallis_d4ec', '144_Mawrth_Vallis_a12a', '77_Mawrth_Vallis_18627', '145_Mawrth_Vallis_b7ca', '78_Mawrth_Vallis_c596', '146_NW_Arabia_Terra_a8f3', '79_Mawrth_Vallis_c18e', '147_Terra_Cimmeria_Terra_Sirenum_ab32', '7_Mawrth_Vallis_a600', '148_Mawrth_Vallis_4eca', '80_Mawrth_Vallis_17c51', '149_NE_Syrtis_b404', '82_Mawrth_Vallis_b2a2', '14_Mawrth_Vallis_aa7d', '83_Mawrth_Vallis_648c', '150_Terra_Cimmeria_Terra_Sirenum_9e4c', '84_Mawrth_Vallis_8838', '15_Mawrth_Vallis_109e0', '85_Mawrth_Vallis_1654f', '17_Mawrth_Vallis_bb59', '86_NE_Syrtis_50f2', '18_Mawrth_Vallis_94f6', '87_NE_Syrtis_a8ec', '19_None_a819', '88_NE_Syrtis_a053', '1_Mawrth_Vallis_3bfb', '89_NE_Syrtis_9971', '20_None_1b8ae', '8_Mawrth_Vallis_b141', '21_Mawrth_Vallis_119dc', '90_NE_Syrtis_21c5a', '22_Circum_Hellas_13d3a', '92_NW_Arabia_Terra_8690', '23_Mawrth_Vallis_1903b', '93_NW_Arabia_Terra_8351', '24_Mawrth_Vallis_16931', '96_None_16f51', '25_Mawrth_Vallis_a425', '97_None_b092', '26_None_8ebf', '98_None_11e2e', '27_NE_Syrtis_821f', '9_NE_Syrtis_48b2', '29_Circum_Hellas_91b1', '2_Mawrth_Vallis_7c34']
  sorted_filename_list = ['1_Mawrth_Vallis_3bfb', '2_Mawrth_Vallis_7c34', '5_Mawrth_Vallis_89f7','6_Mawrth_Vallis_b643', '7_Mawrth_Vallis_a600', '8_Mawrth_Vallis_b141','9_NE_Syrtis_48b2', '11_Mawrth_Vallis_1176a', '13_NE_Syrtis_9d44','14_Mawrth_Vallis_aa7d', '15_Mawrth_Vallis_109e0', '17_Mawrth_Vallis_bb59','18_Mawrth_Vallis_94f6', '19_None_a819', '20_None_1b8ae','21_Mawrth_Vallis_119dc', '22_Circum_Hellas_13d3a', '23_Mawrth_Vallis_1903b','24_Mawrth_Vallis_16931', '25_Mawrth_Vallis_a425', '26_None_8ebf','27_NE_Syrtis_821f', '29_Circum_Hellas_91b1', '30_None_1011a','31_Mawrth_Vallis_13e49', '32_Circum_Hellas_9100', '33_NW_Arabia_Terra_8438','36_Mawrth_Vallis_129cf', '37_Mawrth_Vallis_12d3c', '38_None_137ad','39_Mawrth_Vallis_9326', '40_Terra_Cimmeria_Terra_Sirenum_62b6','41_NW_Arabia_Terra_a8ce', '42_Circum_Hellas_147fb', '43_Circum_Hellas_96be','44_None_a865', '45_Circum_Hellas_16c3b', '46_Mawrth_Vallis_bf57','48_None_1ec17', '50_NE_Syrtis_6408', '53_Terra_Cimmeria_Terra_Sirenum_ada1','55_Terra_Cimmeria_Terra_Sirenum_92f7', '56_Mawrth_Vallis_1ec38','57_Mawrth_Vallis_c467', '58_None_d23e', '59_None_18fcf','60_Mawrth_Vallis_bb21', '61_Circum_Hellas_92b4', '62_Circum_Hellas_8936','63_Circum_Hellas_9665', '64_Circum_Hellas_a5a1', '65_Circum_Hellas_a2b8','66_Circum_Hellas_136f4', '68_Circum_Hellas_12455', '69_Circum_Hellas_14838','70_Circum_Hellas_920a', '71_Circum_Hellas_14241', '73_Circum_Hellas_989c','74_Mawrth_Vallis_90af', '75_Mawrth_Vallis_9cae', '76_Mawrth_Vallis_d4ec','77_Mawrth_Vallis_18627', '78_Mawrth_Vallis_c596', '79_Mawrth_Vallis_c18e','80_Mawrth_Vallis_17c51', '82_Mawrth_Vallis_b2a2', '83_Mawrth_Vallis_648c','84_Mawrth_Vallis_8838', '85_Mawrth_Vallis_1654f', '86_NE_Syrtis_50f2','87_NE_Syrtis_a8ec', '88_NE_Syrtis_a053', '89_NE_Syrtis_9971','90_NE_Syrtis_21c5a', '92_NW_Arabia_Terra_8690', '93_NW_Arabia_Terra_8351','96_None_16f51', '97_None_b092', '98_None_11e2e', '100_None_1957c','101_None_8e5d', '102_None_168c0', '103_None_7af6', '104_None_b278','105_None_ca1b', '106_None_bdd9', '107_None_85d7', '108_None_8ca3','112_Terra_Cimmeria_Terra_Sirenum_c441', '114_Terra_Cimmeria_Terra_Sirenum_ac8b','115_Terra_Cimmeria_Terra_Sirenum_11b46', '116_Terra_Cimmeria_Terra_Sirenum_9bbc','117_Terra_Cimmeria_Terra_Sirenum_96fe', '118_Terra_Cimmeria_Terra_Sirenum_913a','119_Terra_Cimmeria_Terra_Sirenum_9a7b', '120_Terra_Cimmeria_Terra_Sirenum_951c','121_Terra_Cimmeria_Terra_Sirenum_8a1c', '123_Terra_Cimmeria_Terra_Sirenum_97c5','124_Terra_Cimmeria_Terra_Sirenum_91eb', '125_Terra_Cimmeria_Terra_Sirenum_10817','126_Terra_Cimmeria_Terra_Sirenum_18002', '127_Terra_Cimmeria_Terra_Sirenum_8377','128_Terra_Cimmeria_Terra_Sirenum_167fd', '129_Terra_Cimmeria_Terra_Sirenum_7f7e','130_Terra_Cimmeria_Terra_Sirenum_a106', '131_Terra_Cimmeria_Terra_Sirenum_94c2','132_Terra_Cimmeria_Terra_Sirenum_1c624', '133_Terra_Cimmeria_Terra_Sirenum_17199','134_Terra_Cimmeria_Terra_Sirenum_5120', '135_Circum_Hellas_9318','136_Circum_Hellas_98e8', '137_Mawrth_Vallis_c066', '138_Circum_Hellas_59df','139_Mawrth_Vallis_8f1a', '140_Mawrth_Vallis_a955','141_Terra_Cimmeria_Terra_Sirenum_9d7e', '143_Mawrth_Vallis_ac83','144_Mawrth_Vallis_a12a', '145_Mawrth_Vallis_b7ca', '146_NW_Arabia_Terra_a8f3','147_Terra_Cimmeria_Terra_Sirenum_ab32', '148_Mawrth_Vallis_4eca','149_NE_Syrtis_b404', '150_Terra_Cimmeria_Terra_Sirenum_9e4c']

  ;Sixth iteration filename list (n=89), only including the cubes validated by Thanos Klidaras to have CCS using spectral parameter color composites on 12.20.2025
  sixth_iteration_filename_list = ['1_Mawrth_Vallis_3bfb', '5_Mawrth_Vallis_89f7','6_Mawrth_Vallis_b643', '7_Mawrth_Vallis_a600', '8_Mawrth_Vallis_b141', '11_Mawrth_Vallis_1176a', '13_NE_Syrtis_9d44','14_Mawrth_Vallis_aa7d', '17_Mawrth_Vallis_bb59','18_Mawrth_Vallis_94f6', '20_None_1b8ae','21_Mawrth_Vallis_119dc', '22_Circum_Hellas_13d3a', '23_Mawrth_Vallis_1903b','24_Mawrth_Vallis_16931', '25_Mawrth_Vallis_a425', '33_NW_Arabia_Terra_8438','36_Mawrth_Vallis_129cf', '37_Mawrth_Vallis_12d3c','39_Mawrth_Vallis_9326', '40_Terra_Cimmeria_Terra_Sirenum_62b6','41_NW_Arabia_Terra_a8ce', '46_Mawrth_Vallis_bf57','48_None_1ec17', '50_NE_Syrtis_6408', '53_Terra_Cimmeria_Terra_Sirenum_ada1','55_Terra_Cimmeria_Terra_Sirenum_92f7', '56_Mawrth_Vallis_1ec38','57_Mawrth_Vallis_c467', '58_None_d23e', '62_Circum_Hellas_8936','63_Circum_Hellas_9665', '64_Circum_Hellas_a5a1', '65_Circum_Hellas_a2b8','66_Circum_Hellas_136f4', '68_Circum_Hellas_12455', '73_Circum_Hellas_989c','74_Mawrth_Vallis_90af', '76_Mawrth_Vallis_d4ec','77_Mawrth_Vallis_18627', '78_Mawrth_Vallis_c596', '79_Mawrth_Vallis_c18e','80_Mawrth_Vallis_17c51', '82_Mawrth_Vallis_b2a2', '83_Mawrth_Vallis_648c','84_Mawrth_Vallis_8838', '85_Mawrth_Vallis_1654f', '86_NE_Syrtis_50f2','87_NE_Syrtis_a8ec', '88_NE_Syrtis_a053', '89_NE_Syrtis_9971','90_NE_Syrtis_21c5a', '92_NW_Arabia_Terra_8690', '93_NW_Arabia_Terra_8351','96_None_16f51', '100_None_1957c','101_None_8e5d', '103_None_7af6', '104_None_b278', '107_None_85d7', '108_None_8ca3','112_Terra_Cimmeria_Terra_Sirenum_c441', '114_Terra_Cimmeria_Terra_Sirenum_ac8b', '116_Terra_Cimmeria_Terra_Sirenum_9bbc','117_Terra_Cimmeria_Terra_Sirenum_96fe', '118_Terra_Cimmeria_Terra_Sirenum_913a','119_Terra_Cimmeria_Terra_Sirenum_9a7b', '120_Terra_Cimmeria_Terra_Sirenum_951c','121_Terra_Cimmeria_Terra_Sirenum_8a1c','124_Terra_Cimmeria_Terra_Sirenum_91eb', '125_Terra_Cimmeria_Terra_Sirenum_10817','126_Terra_Cimmeria_Terra_Sirenum_18002', '127_Terra_Cimmeria_Terra_Sirenum_8377','130_Terra_Cimmeria_Terra_Sirenum_a106', '131_Terra_Cimmeria_Terra_Sirenum_94c2','132_Terra_Cimmeria_Terra_Sirenum_1c624', '133_Terra_Cimmeria_Terra_Sirenum_17199','134_Terra_Cimmeria_Terra_Sirenum_5120', '135_Circum_Hellas_9318','136_Circum_Hellas_98e8', '137_Mawrth_Vallis_c066', '140_Mawrth_Vallis_a955','141_Terra_Cimmeria_Terra_Sirenum_9d7e', '143_Mawrth_Vallis_ac83', '145_Mawrth_Vallis_b7ca', '146_NW_Arabia_Terra_a8f3', '148_Mawrth_Vallis_4eca','149_NE_Syrtis_b404', '150_Terra_Cimmeria_Terra_Sirenum_9e4c']
  chosen_file=sixth_iteration_filename_list[value]
  dir='C:/Users/aklidara/OneDrive - purdue.edu/Paper2 data backup 11.24.2025/paper2 CRISM MRDR images/'+sixth_iteration_filename_list[value]+'/'



  split_filename=strsplit(chosen_file,/EXTRACT,'_')
  id=split_filename[0]
  end_code=split_filename[-1]
  if strlen(end_code) eq 5 then file="frt000"+end_code+"_07"
  if strlen(end_code) eq 4 then file="frt0000"+end_code+"_07"
  ;Exceptions for our non-FRT misfits
  if id eq '40' then file="hrl000062b6_07"
  if id eq '20' then file="hrl0001b8ae_07"
  if id eq '11' then file="hrs0001176a_07"
  if id eq '50' then file="hrl00006408_07"
  if id eq '68' then file="hrl00012455_07"
  if id eq '80' then file="hrl00017c51_07"
  if id eq '83' then file="hrl0000648c_07"
  if id eq '87' then file="hrl0000a8ec_07"
  if id eq '126' then file="hrl00018002_07"
  if id eq '149' then file="hrl0000b404_07"




  print,'Now calculating for CRISM image '+chosen_file+' (filename='+file+')'

  label_files = FILE_SEARCH(dir + file + '_if*_mtr3.lbl')
  ;print,label_files

  IF N_ELEMENTS(label_files) EQ 0 THEN BEGIN
    MESSAGE, 'Label file not found: ' + dir + file + '_if*_mtr3.lbl'
    RETURN
  ENDIF
  label_file = label_files[0]

  OPENR, unit, label_file, /GET_LUN
  temp = ''
  samples = 0
  lines = 0
  WHILE NOT EOF(unit) DO BEGIN
    READF, unit, temp
    ; Remove whitespace
    temp = STRCOMPRESS(temp, /REMOVE_ALL)
    ; Extract number of samples
    IF STRPOS(temp, 'LINE_SAMPLES=') GE 0 THEN BEGIN
      samples_line = STRSPLIT(temp, '=', /EXTRACT)
      samples = FIX(samples_line[1])
    ENDIF
    ; Extract number of lines
    IF STRPOS(temp, 'LINES=') GE 0 THEN BEGIN
      lines_line = STRSPLIT(temp, '=', /EXTRACT)
      lines = FIX(lines_line[1])
    ENDIF
  ENDWHILE
  CLOSE, unit
  FREE_LUN, unit

  ; Check if necessary parameters were found
  IF samples EQ 0 OR lines EQ 0 THEN BEGIN
    MESSAGE, 'Failed to read array dimensions from label file.'
    RETURN
  ENDIF

  data_type = 4  ; 4 corresponds to floating-point data type in IDL

  ; Initialize data arrays
  rfl = FLTARR(samples, lines, 489)
  params = FLTARR(samples, lines, 60)
  loc = FLTARR(samples, lines, 14)

  ; Read in location data
  loc_files = FILE_SEARCH(dir + file + '_de*l_mtr1.img')
  IF N_ELEMENTS(loc_files) GT 0 THEN BEGIN
    loc_file = loc_files[0]
    OPENR, unit, loc_file, /GET_LUN
    READU, unit, loc
    CLOSE, unit
    FREE_LUN, unit
    loc = TRANSPOSE(loc, [0, 1, 2])
  ENDIF ELSE BEGIN
    loc = 0
    PRINT, 'Location file not found: ', dir + file + '_de*l_mtr1.img'
  ENDELSE

  ; Read in parameter data
  params_files = FILE_SEARCH(dir + file + '_sr*j_mtr3.img')
  IF N_ELEMENTS(params_files) GT 0 THEN BEGIN
    params_file = params_files[0]
    OPENR, unit, params_file, /GET_LUN
    READU, unit, params
    CLOSE, unit
    FREE_LUN, unit
    params = TRANSPOSE(params, [0, 1, 2])
  ENDIF ELSE BEGIN
    params = 0
    PRINT, 'Parameter file not found: ', dir + file + '_sr*j_mtr3.img'
  ENDELSE

  ; Load reflectance data
  data_files = FILE_SEARCH(dir + file + '_if*j_mtr3.img')
  IF N_ELEMENTS(data_files) GT 0 THEN BEGIN
    data_file = data_files[0]
    OPENR, unit, data_file, /GET_LUN
    READU, unit, rfl
    CLOSE, unit
    FREE_LUN, unit
    rfl = TRANSPOSE(rfl, [0, 1, 2])
  ENDIF ELSE BEGIN
    MESSAGE, 'Data file not found: ' + dir + file + '_if*j_mtr3.img'
    RETURN
  ENDELSE

  test=master_paper2(rfl,params,id)


END
