# Chapter 1 Introduction - Purpose and Theme Points

## Intended use

These notes are for the short opening directly below `\chapter{Introduction}`. The opening should introduce the broad problem, the purpose of the thesis, and the connection between SIF enhancement and winter-wheat yield prediction. Detailed motivation, research questions, methods, and numerical results should remain in their later sections.

## Core points suitable for the introduction

1. Photosynthesis drives vegetation growth, terrestrial carbon uptake, and the formation of crop biomass. Observing photosynthetic activity over large areas is therefore relevant to both carbon-cycle research and agricultural monitoring.

2. Conventional reflectance-based vegetation indices mainly describe canopy greenness and structure. Solar-induced chlorophyll fluorescence (SIF) is emitted by chlorophyll during photosynthetic energy conversion and provides a more functional indicator of current photosynthetic activity.

3. Satellite SIF has enabled regional and global monitoring of photosynthesis, vegetation stress, GPP, and crop productivity. Its agricultural value is particularly relevant when heat or drought reduces photosynthetic activity before strong changes become visible in canopy greenness.

4. Existing satellite observations involve an important spatial-sampling trade-off. OCO-2 provides comparatively small SIF footprints, but records them only along narrow, widely separated tracks. TROPOMI provides much more complete and frequent coverage, but at a coarser spatial resolution.

5. SIF is a weak signal relative to reflected sunlight and individual retrievals contain substantial noise. Spatial or temporal averaging improves signal stability but reduces the spatial detail needed in heterogeneous agricultural landscapes.

6. Machine learning provides a way to combine sparse or coarse SIF observations with spatially complete auxiliary data. Reflectance, vegetation indices, radiation, environmental variables, seasonal information, and land-cover or crop information can be used to reconstruct missing SIF and infer how the signal varies within a coarse observational support.

7. Previous products such as SIFoco2_005 demonstrate that sparse OCO-2 observations can be reconstructed into spatially continuous SIF fields. SIFnet further shows that a spatial model can combine a parent SIF observation with finer auxiliary rasters to infer sub-grid spatial patterns.

8. Nevertheless, most established products remain at approximately 500 m or coarser resolution. This remains limiting for European agricultural regions where individual fields and crop mixtures vary over much shorter distances.

9. The central theme of this thesis is therefore scale bridging: learning from sparse footprint- or grid-supported OCO-2 observations while using high-resolution Sentinel-2 and ancillary variables to produce spatially detailed SIF estimates.

10. The thesis also investigates whether the choice of input representation and spatial support matters. Tabular Random Forest and XGBoost models, an area-to-point neural network, and spatial CNN/U-Net models are compared, while the effect of aggregating noisy OCO-2 supervision to 4 km cells is examined separately.

11. The predicted maps are produced on a 20 m grid, but they are learned from aggregate or footprint-level SIF observations. The introduction should call them high-resolution model estimates rather than direct 20 m SIF observations.

12. The agricultural purpose is to isolate SIF over winter-wheat pixels and test whether this crop-pure enhanced signal provides more useful information for Bavarian winter-wheat yield prediction than raw OCO-2 footprints containing mixed crops and other land cover.

13. The overall purpose can be stated as developing and evaluating an OCO-2-based SIF enhancement workflow for Germany and assessing whether its additional spatial detail has practical value for crop-specific yield modelling.

## Evidence and ideas from the four main source papers

### SIFnet - Gensheimer et al. (2022)

- GPP is an important part of the global carbon cycle, but it cannot be measured directly from present satellite systems; remotely observed quantities are therefore used as proxies.
- Vegetation indices describe photosynthetic capacity or canopy greenness, whereas SIF more directly indicates photosynthetic activity.
- Processes affecting GPP occur at spatial scales finer than most satellite SIF products can resolve.
- SIFnet demonstrates that coarse SIF can be combined with finer auxiliary rasters in a CNN to estimate a finer SIF field.
- The parent SIF observation was the most important input to SIFnet, while NIRv was the most important auxiliary variable. This supports combining actual SIF measurements with high-resolution optical predictors rather than treating reflectance as a complete replacement for SIF.
- Independent OCO-2/OCO-3 evaluation illustrates the importance of validating enhanced predictions at the spatial support of an independent observation.

### SIFoco2_005 - Yu et al. (2019)

- OCO-2 offers relatively fine SIF footprints along its orbit, but the narrow tracks leave large spatial gaps and prevent spatially continuous monitoring.
- The SIF signal is weak and noisy; aggregation reduces random error but creates a trade-off between signal stability and spatial resolution.
- Coarse SIF is especially limiting in heterogeneous ecosystems because mixed land-cover signals complicate comparison with local GPP and vegetation conditions.
- Spatially complete MODIS observations can be used to learn relationships along OCO-2 tracks and extend them into the gaps.
- The relationship between SIF and reflectance varies with biome, phenology, and stress. Reflectance mainly describes structural change, while SIF may respond to physiological stress before visible structural change occurs.
- Spatially continuous SIF can support drought monitoring, agricultural planning, and yield estimation, but its inferred fine-scale variation still depends on the predictors and modelling assumptions.

### Global GOSAT, OCO-2, and OCO-3 SIF datasets - Doughty et al. (2022)

- SIF is produced when a small part of the energy absorbed by chlorophyll is re-emitted during photosynthesis.
- Satellite SIF retrieval is comparatively recent and requires instruments with high spectral resolution and signal-to-noise ratio because fluorescence is weak relative to background radiance.
- OCO-2 and OCO-3 provide some of the smallest spaceborne SIF footprints, but their narrow swaths create sparse spatial coverage.
- SIF interpretation must account for retrieval uncertainty, observation geometry, cloud filtering, and the indirect relationship between observed canopy SIF and photosynthesis.
- These points justify careful target preparation, uncertainty filtering, footprint-aware spatial matching, and cautious interpretation of the final 20 m estimates.

### Machine Learning for Satellite SIF - Verrelst et al. (2026)

- Current satellite missions differ greatly in spatial resolution, temporal frequency, and sampling strategy; together they provide useful but incomplete observations of vegetation photosynthesis.
- Scientific use is limited by low signal-to-noise ratios, atmospheric and directional effects, coarse resolution, irregular sampling, and scale mismatch with heterogeneous vegetation.
- Machine learning is increasingly used after the original SIF retrieval to reconstruct, gap-fill, and downscale SIF by combining it with multi-source auxiliary information.
- ML can produce higher-level products useful for GPP and vegetation-stress applications, but robust validation, uncertainty reporting, and testing under geographic or temporal domain shift remain essential.
- This supports the thesis focus on spatially separated evaluation, comparison of data representations, and explicit distinction between output-grid resolution and validation support.

## Thesis-specific purpose and contribution points

- **Study domain:** selected Sentinel-2 MGRS tiles covering diverse agricultural and agro-climatic regions of Germany.
- **SIF target:** quality-controlled OCO-2 observations from 2019-2024, including a noise-reduced target formed from the 757 and 771 nm retrievals.
- **High-resolution information:** Sentinel-2 spectral indices on a 20 m grid, combined with FAPAR, PAR/APAR, seasonal variables, and annual crop-composition channels.
- **Methodological comparison:** tabular footprint means, area-to-point learning, single- and multi-footprint CNNs, and a spatially aggregated 4 km U-Net representation.
- **Aggregation question:** determine whether averaging nearby OCO-2 observations produces more stable supervision and improves prediction, while acknowledging that it narrows the target distribution.
- **Generalization question:** evaluate the selected model on OCO-2 footprints from previously unseen MGRS tiles.
- **Agricultural application:** predict SIF over winter-wheat pixels, summarize it by month and Bavarian NUTS-3 region, and compare its yield-prediction value with raw mixed-footprint OCO-2 SIF.
- **Interpretation boundary:** a 20 m output grid provides spatial detail guided by high-resolution predictors; validation after aggregation to OCO-2 support does not prove that each individual 20 m value is correct.

## Suggested two-paragraph sequence

### Paragraph 1 - Broad problem and need

1. Begin with the importance of photosynthesis for vegetation productivity, the carbon cycle, and crop production.
2. Introduce SIF as a satellite-observable signal connected to photosynthetic activity and distinguish it briefly from structural greenness indices.
3. State the main observational limitation: available SIF is either spatially sparse or comparatively coarse and noisy.
4. Explain why this is particularly problematic for heterogeneous agricultural landscapes and crop-specific monitoring.

### Paragraph 2 - Purpose and theme of this thesis

1. Introduce machine learning as the bridge between OCO-2 observations and spatially complete high-resolution predictors.
2. State that the thesis combines OCO-2 SIF with Sentinel-2, radiation, seasonal, and crop-composition information over Germany.
3. Mention the comparison of tabular and spatial models and the investigation of aggregation as a way to stabilize noisy supervision.
4. End with the downstream purpose: assessing whether crop-pure enhanced SIF improves winter-wheat yield prediction relative to raw mixed OCO-2 footprints.

## Papers cited by the main sources that are suitable for direct citation

The opening should cite original studies for its central scientific claims. The following is a deliberately small selection from the reference lists of the four main papers.

### Core SIF and photosynthesis background

1. **Porcar-Castell et al. (2014), _Linking chlorophyll a fluorescence to photosynthesis for remote sensing applications: Mechanisms and challenges_.** Useful for explaining the physical connection between fluorescence and photosynthesis and the limits of interpreting SIF as a direct measurement of GPP. DOI: `10.1093/jxb/eru191`. Not currently found in `manual.bib`.

2. **Mohammed et al. (2019), _Remote sensing of solar-induced chlorophyll fluorescence (SIF) in vegetation: 50 years of progress_.** Broad background source for the development and significance of SIF remote sensing. DOI: `10.1016/j.rse.2019.04.030`. Not currently found in `manual.bib`.

3. **Frankenberg et al. (2011), _New global observations of the terrestrial carbon cycle from GOSAT: Patterns of plant fluorescence with gross primary productivity_.** Foundational evidence for the relationship between spaceborne SIF and GPP at large scales. DOI: `10.1029/2011GL048738`. Not currently found in `manual.bib`.

4. **Guanter et al. (2014), _Global and time-resolved monitoring of crop photosynthesis with chlorophyll fluorescence_.** Particularly relevant for introducing agricultural monitoring with satellite SIF. DOI: `10.1073/pnas.1320008111`. Not currently found in `manual.bib`.

### OCO-2 observations and their limitations

5. **Sun et al. (2017), _OCO-2 advances photosynthesis observation from space via solar-induced chlorophyll fluorescence_.** Strong original source for the scientific value of OCO-2 SIF and its relatively fine footprint observations. DOI: `10.1126/science.aam5747`. Not currently found in `manual.bib`.

6. **Sun et al. (2018), _Overview of Solar-Induced Chlorophyll Fluorescence (SIF) from the Orbiting Carbon Observatory-2: Retrieval, cross-mission comparison, and global monitoring for GPP_.** Useful for OCO-2 retrieval characteristics, wavelength combination, quality control, and GPP relevance. DOI: `10.1016/j.rse.2018.02.016`. Not currently found in `manual.bib`.

7. **Doughty et al. (2022), _Global GOSAT, OCO-2, and OCO-3 solar-induced chlorophyll fluorescence datasets_.** Appropriate source for comparing satellite sampling, footprint size, retrieval uncertainty, and observation geometry. DOI: `10.5194/essd-14-1513-2022`. Not currently found in `manual.bib`.

### Reconstruction and high-resolution enhancement

8. **Yu et al. (2019), _High-resolution global contiguous SIF of OCO-2_.** Direct precedent for reconstructing sparse OCO-2 SIF and for discussing the need to represent biome and temporal variation. Existing key: `\cite{yu2019high}`.

9. **Gensheimer et al. (2022), _A convolutional neural network for spatial downscaling of satellite-based solar-induced chlorophyll fluorescence (SIFnet)_.** Direct precedent for CNN-based SIF downscaling, use of the parent SIF observation, and independent footprint validation. Existing key: `\cite{gensheimer2022convolutional}`.

10. **Dechant et al. (2020), _Canopy structure explains the relationship between photosynthesis and sun-induced chlorophyll fluorescence in crops_.** Useful for explaining why crop composition and canopy structure matter when interpreting or decomposing mixed SIF signals. DOI: `10.1016/j.rse.2020.111733`. Not currently found in `manual.bib`.

### Agricultural application

11. **Peng et al. (2020), _Assessing the benefit of satellite-based solar-induced chlorophyll fluorescence in crop yield prediction_.** Supports the connection between spatially enhanced SIF and maize/soybean yield prediction. Existing key: `\cite{peng2020assessing}`.

12. **Song et al. (2018), _Satellite sun-induced chlorophyll fluorescence detects early response of winter wheat to heat stress in the Indian Indo-Gangetic Plains_.** Supports the value of SIF for detecting winter-wheat physiological stress before greenness indices. Existing key: `\cite{song2018satellite}`.

13. **Qiu et al. (2022), _Monitoring drought impacts on crop productivity of the U.S. Midwest with solar-induced fluorescence_.** Supports the sensitivity of reconstructed SIF to drought-related changes in crop productivity and yield. Existing key: `\cite{qiu2022monitoring}`.

## Recommended citation restraint for the opening

- Use approximately four to six citations across the two opening paragraphs.
- Prioritize one foundational SIF/photosynthesis paper, one OCO-2 data paper, one enhancement paper, and one agricultural application paper.
- A suitable compact set would be Porcar-Castell et al. (2014), Sun et al. (2017 or 2018), Yu et al. (2019), Gensheimer et al. (2022), and either Peng et al. (2020) or Song et al. (2018).
- Reserve detailed product comparisons, performance values, and validation limitations for Chapters 3, 5, and 6.

## Statements to avoid or qualify

- Do not describe satellite SIF as a direct measurement of GPP; it is an emitted radiative signal related to photosynthetic processes and commonly used as a proxy for GPP.
- Do not claim that SIF is always a better yield predictor than every vegetation index; its advantage varies with crop, scale, stress condition, and validation design.
- Do not describe the generated 20 m values as observed or independently validated 20 m SIF.
- Do not imply that the model establishes causal effects of individual predictors.
- Do not present the current yield experiment as an operational yield model; its role is to test the relative value of crop-pure enhanced SIF.
