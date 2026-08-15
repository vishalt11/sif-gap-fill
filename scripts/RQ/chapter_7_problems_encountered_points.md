# Chapter 7: Problems Encountered — Writing Points

The section can be written in the following order. For each problem, briefly state the reason for the methodological decision, the action taken in this thesis, and the limitation that remains.

## 1. Construction of the 4 km spatial target

- OCO-2 observations occur along narrow and unevenly sampled tracks, so the number and spatial arrangement of observations differ strongly among locations and dates.
- Density-based clustering was considered for grouping nearby observations. However, this could produce clusters with irregular shapes and different spatial extents, making the resulting SIF targets less comparable and difficult to represent using fixed-size U-Net inputs.
- A fixed 4 km grid was therefore used to provide a consistent spatial support and a regular $200\times200$ input at 20 m resolution.
- Each retained grid cell was required to contain at least five OCO-2 observations from the same date, MGRS tile, and measurement mode. Averaging these observations reduced individual-retrieval noise and provided a more stable training target.
- This decision also introduced limitations. Sparse grid cells were discarded, footprints close to a grid boundary could be separated into different cells, and the retained dataset became biased towards locations with denser OCO-2 sampling.
- Spatial averaging also reduced the observed SIF range and removed many extreme values. The resulting model was therefore better at representing the central SIF range than unusually low or high observations.

## 2. Why observations were not grouped across time

- Temporal aggregation over 8 or 16 days was not used because OCO-2 does not provide spatially complete observations during these periods.
- Keeping the spatial support fixed at 4 km would often still provide too few repeat observations. Obtaining enough observations by expanding the grouping would instead combine separate OCO-2 track segments that may be located far apart.
- The resulting target would no longer represent one clearly defined 4 km area, and its effective spatial support could vary greatly between samples.
- This differs from wide-swath TROPOMI data, for which multi-day composites are more likely to contain repeated or near-complete observations of the same location.
- Combining several OCO-2 dates could also mix different weather conditions and short-term vegetation responses. Same-date spatial grouping was therefore retained even though it produced fewer training samples.

## 3. Temporal mismatch between SIF and the predictors

- The input datasets do not describe vegetation conditions over exactly the same period. OCO-2 SIF is attached to an individual acquisition date, VIIRS PAR is a daily mean, GLASS FAPAR is an eight-day composite, and Sentinel-2 WASP reflectance represents a monthly composite.
- The predictors were matched as closely as their temporal definitions allowed: PAR by exact date, FAPAR by the eight-day interval containing the SIF date, and Sentinel-2 by its monthly product interval.
- Nevertheless, the model output is not a fully instantaneous estimate of the vegetation state at the exact OCO-2 overpass time. It combines information representing daily, eight-day, and monthly conditions.
- Native OCO-2 validation can consequently be noisy because the observed SIF may contain short-term physiological responses that are not represented in the temporally averaged predictors.
- This limitation should be considered when interpreting the weaker performance against individual OCO-2 footprints. A future workflow could use optical and environmental predictors acquired closer to the SIF overpass time, provided that cloud-free high-resolution observations are available.

## 4. Predictor resolutions and alignment

The following table can be converted into a LaTeX table in the thesis.

| Dataset or variable | Variables used | Native spatial resolution | Temporal coverage or support | Aggregation/downscaling and temporal matching |
|---|---|---:|---|---|
| OCO-2 SIF target | Combined 757 and 771 nm SIF | Irregular footprints, approximately 1 km$^2$ after processing | Individual dates, February–July 2019–2024 | At least five same-date and same-mode observations averaged within each fixed 4 km cell; individual footprint masks retained for external validation |
| Sentinel-2 L3A WASP | NDVI, NDMI, NDRE, NIRv and EVI derived from reflectance bands | Bands at 10 or 20 m | Monthly composites, February–July 2019–2024 | 10 m bands area-averaged to 20 m; native 20 m bands retained; matched using the monthly product interval |
| GLASS FAPAR | FAPAR | Approximately 231.7 m | Eight-day composites, 2019–2024 | Bilinear interpolation to the 20 m grid; matched to the composite interval containing the OCO-2 date |
| VIIRS VNP18A2 | PAR | Approximately 927 m | Daily mean, 2019–2024 | Quality filtering on the native grid followed by bilinear interpolation to 20 m; matched by exact date |
| Derived APAR | FAPAR $\times$ PAR | 20 m model grid after alignment | Combined eight-day FAPAR and daily PAR support | Calculated after FAPAR and PAR were aligned to the same 20 m grid |
| German crop-type maps | Crop fractions, active-crop fraction and non-crop fraction | 10 m | Annual maps, 2019–2024 in this study | Binary crop masks area-averaged from 10 to 20 m; matched by year, with crop activity defined by month |
| Seasonal encoding | Sine and cosine of month | Spatially constant channels | Observation month | Repeated over every 20 m pixel in the input window |

- Resampling FAPAR and PAR to 20 m only aligns their grids; it does not create genuine 20 m information from the original 231.7 m and 927 m products.
- Similarly, the model produces a 20 m SIF grid, but the observational supervision remains at OCO-2 footprint or 4 km aggregate support. Accuracy at the coarse validation support does not by itself prove the accuracy of each 20 m pixel.

## 5. Sentinel-2 data volume, processing time, and geographical scope

- Many previous SIF-reconstruction studies used MODIS predictors because they provide regular global composites at moderate resolution and are easier to store and process. Sentinel-2 provides much finer field detail, but substantially increases data volume and processing complexity.
- Moving from 500 m to 20 m increases the number of pixels covering the same area by a factor of 625. A complete Sentinel-2 MGRS tile contains approximately 30 million 20 m pixels per layer, before adding dates, predictor channels, masks, and intermediate products.
- The Sentinel-2 archive for the selected 11 tiles approached approximately 500 GB before all model-ready processing was completed. Reprojection, masking, band alignment, index calculation, chip construction, and map inference therefore required considerable storage and computing time.
- These requirements limited the geographical and temporal scope of the study. The selected tiles represent several important German agricultural regions, but they do not provide complete wall-to-wall coverage of all German landscapes and conditions.
- Increasing the geographical scope could add more examples from the low and high ends of the SIF distribution and could reduce the model's tendency to predict values near the centre of the training range. However, this would require substantially more data acquisition, storage, preprocessing, and model-inference time.

## 6. Uncertainty inherited from the crop-type rasters

- The annual crop-type rasters were themselves produced using Random Forest classification and spatial-temporal filtering of Sentinel-1 and Sentinel-2 observations. They are therefore model outputs rather than error-free field records.
- The seven-year German product reported annual overall accuracies of approximately 0.81–0.83. Winter wheat was among a group of crops with multi-year average F1-scores between approximately 0.79 and 0.84 \cite{gessner_crop_2025}.
- The earlier 2018 Germany-wide map reported a winter-wheat F1-score of 0.82 and an overall map accuracy of 75.5\% \cite{asam_crop_2022}.
- These are classification results, so $R^2$ and RMSE were not reported and are not the appropriate metrics for this dataset.
- Classification errors can propagate into the SIF model because crop fractions are used as predictors. They can also affect the yield analysis because the enhanced SIF is extracted only from pixels classified as winter wheat.
- Confusion among spectrally and phenologically similar winter cereals is particularly relevant. False winter-wheat classifications can introduce other crops into the supposedly crop-pure SIF mean, while omitted wheat pixels reduce the available crop area.
- Requiring all four underlying 10 m pixels to be classified as winter wheat before accepting one 20 m pixel reduced mixed boundary pixels, but it could not remove classification errors occurring within fields.

## 7. Raman-scattering bias in the OCO-2 SIF input

- The V11r OCO-2 SIF product used in this thesis does not include the new explicit correction for rotational Raman scattering described by Sanghavi et al. \cite{sanghavi2025impact}.
- Raman scattering can fill solar Fraunhofer lines in a way that resembles the smooth fluorescence signal. The study showed that this can create artificial positive SIF values as high as approximately $0.35\,\mathrm{mW\,m^{-2}\,sr^{-1}\,nm^{-1}}$ over bright non-vegetated surfaces.
- Existing zero-level corrections remove part of this effect, but residual biases of approximately $\pm0.15\,\mathrm{mW\,m^{-2}\,sr^{-1}\,nm^{-1}}$ can remain and vary with surface albedo, surface pressure or elevation, solar zenith angle, and viewing geometry.
- These biases originate in the satellite retrieval and may therefore be inherited by a machine-learning enhancement model. Because reflectance and seasonal variables are also used as predictors, the model could partly learn retrieval artefacts associated with surface brightness or observation geometry.
- Sanghavi et al. proposed a lookup-table correction based on high-resolution physical simulations. The lookup table estimates the inelastic signal from surface albedo, surface pressure, solar zenith angle, and viewing angle so that it can be subtracted from the retrieved SIF.
- The paper states that this correction is being incorporated into Build 12 of the operational OCO-2/3 Level-2 products. Repeating the workflow with a corrected future product would help reduce this source of spatial and temporal bias.

## 8. Limitations of the winter-wheat yield experiment

- Producing the 20 m SIF map for one MGRS tile required approximately 45 minutes in the implemented workflow. Extending inference over more tiles, months, and years therefore became computationally expensive.
- The final yield dataset contained 197 region-year training records and only 48 NUTS-3 regions in the independent 2024 test set. The small test set produced wide bootstrap intervals, and the observed improvement from crop-pure SIF was not statistically conclusive.
- Although the SIF maps were generated at 20 m, the winter-wheat pixels were averaged to one monthly value for each NUTS-3 region. Much of the field-scale spatial information was therefore averaged out before yield prediction.
- Each MGRS tile covers roughly 10–12 NUTS-3 regions, of which approximately 4–6 may be independent cities (\emph{kreisfreie St\"adte}). These urban districts generally contain little agricultural land and often do not report crop yield, so processing another complete tile produces only a limited number of additional usable yield records.
- More NUTS-3 regions and years are needed, but parcel-level yield observations would be more suitable for testing whether 20 m crop-pure SIF provides a genuine field-scale advantage.
- Brandt et al. used measured parcel-level yields to develop yield models for German crops \cite{brandt2024ensemble}. Those farm-level reference records were not publicly downloadable because of farmer confidentiality, so this thesis had to rely on publicly available NUTS-3 statistics.
- The official response category was \emph{Winterweizen (einschl. Dinkel und Einkorn)}, meaning that the reported yield combines winter wheat with spelt and einkorn, whereas the crop mask was intended to identify winter wheat. This creates a remaining mismatch between the satellite predictor and the statistical response.
- The yield experiment was limited to winter wheat and omitted detailed weather, soil, cultivar, and crop-management variables. It was designed primarily to compare crop-pure enhanced SIF with raw mixed-footprint SIF, not to construct a complete operational yield model.
- The results therefore represent predictive associations within the selected Bavarian regions and years. They should not be interpreted as causal effects or assumed to transfer automatically to other crops, years, or geographical areas.

## Suggested paragraph order in Chapter 7

1. Explain the sparse OCO-2 sampling problem, the fixed 4 km grid, and the minimum-five-observation rule.
2. Explain why same-date spatial grouping was preferred to 8- or 16-day temporal grouping.
3. Discuss the different temporal and spatial supports of OCO-2, Sentinel-2, FAPAR, and PAR; place the predictor table here.
4. Discuss Sentinel-2 data volume, processing time, limited geographical coverage, and the restricted SIF range.
5. Discuss uncertainty inherited from the machine-learning crop maps.
6. Discuss the uncorrected Raman-scattering contribution in V11r and the planned Build 12 correction.
7. End with the limitations of the NUTS-3 winter-wheat yield experiment and the need for additional regions, years, predictors, and confidential parcel-level yield observations.

# Outlook points

## 1. Use of FLEX SIF observations

- ESA's FLEX mission is scheduled for launch in September 2026 and will use the FLORIS imaging spectrometer to provide SIF observations at approximately 300 m resolution.
- FLEX will provide a much denser and finer base SIF dataset than the sparse OCO-2 tracks used in this thesis. However, 300 m remains too coarse for many individual agricultural fields.
- A future study could combine FLEX SIF with Sentinel-2 predictors to generate field-scale estimates while constraining the mean fine-resolution prediction to reproduce the original FLEX observation.
- The reviewed outlook literature also identifies FLEX downscaling, temporal reconstruction, multi-sensor fusion, physically constrained modelling, and uncertainty propagation as important future applications of machine learning.
- FLEX's tandem operation with Sentinel-3 could additionally provide near-simultaneous information about surface reflectance, temperature, atmospheric conditions, and vegetation state.
- Relevant FLEX reference: Drusch et al., *The fluorescence explorer mission concept—ESA's Earth Explorer 8*.

## 2. Reprocess the model using corrected OCO-2/3 products

- When the corrected Build 12 OCO-2/3 Level-2 SIF products become available, the workflow should be repeated using these observations instead of V11r.
- Build 12 is expected to include the lookup-table correction for the false SIF-like signal caused by rotational Raman scattering \cite{sanghavi2025impact}.
- This would reduce retrieval biases associated with surface albedo, surface pressure, solar zenith angle, and viewing geometry before the SIF values are used as machine-learning targets.

## 3. Use TROPOMI as an observation-preserving target

- A future model could use spatially complete TROPOMI SIF as the parent observation and Sentinel-2 as the high-resolution spatial predictor.
- The fine-resolution predictions could be constrained so that their average within each TROPOMI pixel reproduces the observed coarse SIF value. This would prevent the enhanced map from replacing the measured signal with an unconstrained predictor-based estimate.
- This approach would provide stronger observation preservation than the present OCO-2 reconstruction. Its main limitation is that a suitable TROPOMI observation must exist for the location and date being downscaled.

## 4. Improve temporal matching of Sentinel-2 and SIF

- Instead of monthly WASP composites, future work could search for Sentinel-2 Level-2A scenes acquired closest to each OCO-2 or other SIF observation date.
- Closer temporal matching could make the predictors more representative of vegetation conditions at the SIF overpass and reduce noise during native-footprint validation.
- This would require explicit cloud and shadow filtering, selection or compositing of multiple candidate scenes, and rules for the maximum acceptable time difference.
- Some locations may have no valid Sentinel-2 observation for several days or weeks because of cloud cover. Consequently, closer matching would considerably increase data acquisition and processing time and could also reduce the number of usable SIF samples.

## 5. Acquire parcel-level crop-yield observations

- Future yield evaluation should seek confidential parcel-level yield measurements from the relevant German agricultural or statistical authorities, similar to the dataset used by Brandt et al. \cite{brandt2024ensemble}.
- Access would probably require a formal data-use agreement, secure processing, and publication of only anonymized or aggregated results to protect farmer confidentiality.
- Parcel-level observations would allow the 20 m crop-pure SIF estimates to be evaluated near their intended spatial scale instead of first averaging them to NUTS-3 regions.
- A larger parcel dataset should also include weather, soil, cultivar, and management variables so that enhanced SIF can be assessed as one component of a more complete yield model.

## Overall future direction

- Future development should combine denser SIF missions, temporally closer predictors, observation-preserving constraints, corrected retrieval products, and field-scale reference data.
- Uncertainty should be carried from the original SIF retrieval through reconstruction or downscaling and into the final crop-yield estimates, rather than reporting model error only at the last processing stage.
