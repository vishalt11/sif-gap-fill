# Literature Review

## Section 1 - SIF enhancement papers

Abstract-level screening criterion: a paper is treated as relevant when it predicts, reconstructs,
gap-fills, spatially or temporally downscales, forecasts, or directly improves the retrieval precision
of SIF. Papers that only use or normalize existing SIF without producing an enhanced SIF estimate or
product are treated as non-relevant to Section 3.1.

After removal of the duplicate files and relocation of the non-relevant paper, the Section 3.1
folder contains 21 unique, relevant papers.

### Relevant

1. **A Global, 0.05-Degree Product of Solar-Induced Chlorophyll Fluorescence Derived from OCO-2, MODIS, and Reanalysis Data** - predicts spatially continuous, 0.05-degree, 8-day GOSIF from sparse OCO-2 observations.
2. **A Spatial Downscaling Method for Solar-Induced Chlorophyll Fluorescence Product Using Random Forest Regression and Drought Monitoring in Henan Province** - random-forest spatial downscaling of GOSIF for regional drought monitoring.
3. **A Spatially Downscaled Sun-Induced Fluorescence Global Product for Enhanced Monitoring of Vegetation Productivity** - downscales GOME-2 SIF to 0.05-degree and 8-day resolution.
4. **A Spatially Downscaled TROPOMI SIF Product at 0.005-Degree Resolution With Bias Correction** - reconstructs bias-corrected, finer-resolution TROPOMI SIF.
5. **Downscaling Solar-Induced Chlorophyll Fluorescence Based on Convolutional Neural Network Method to Monitor Agricultural Drought** - CNN-based downscaling from 0.05 to 0.008 degrees.
6. **Downscaling Solar-Induced Chlorophyll Fluorescence to a 0.05-Degree Monthly Product Using AVHRR Data in East Asia (1995-2003)** - random-forest kriging reconstruction of pre-2000 fine-resolution SIF.
7. **Forecasting of Fine-Grained SIF of OCO-2 Using Multi-source Data and AI-Based Techniques** - forecasts monthly and seasonal OCO-2 SIF using EVI, data fusion, and machine learning.
8. **Generation of a Global Spatially Continuous TanSat Solar-Induced Chlorophyll Fluorescence Product by Considering the Impact of the Solar Radiation Intensity** - predicts spatially continuous TanSat SIF with random forest and MODIS/environmental predictors.
9. **DOSIF: Long-Term Daily SIF from OCO-3 with Global Contiguous Coverage** - reconstructs a long-term, daily, spatially continuous OCO-3-based SIF product.
10. **A More Precise Retrieval of Sun-Induced Chlorophyll Fluorescence from Satellite Data Using Artificial Neural Networks** - uses an ANN to reduce retrieval noise while retaining finer effective spatial and temporal resolution. This is retrieval enhancement rather than gap-filling or downscaling, but it is still relevant to the broader enhancement-method discussion.
11. **RTSIF: A Long-Term Reconstructed TROPOMI Solar-Induced Fluorescence Dataset Using Machine Learning Algorithms** - reconstructs TROPOMI-like SIF for 2001-2020 at 0.05-degree, 8-day resolution.
12. **Tracking Drought in Dryland Vegetation through the Photosynthetic Afternoon Depression Index of Sun-Induced Chlorophyll Fluorescence** - reconstructs high-spatial-resolution (500 m) SIF from GOCI and OCO-3 before deriving a drought indicator.
13. **Hybrid Machine Learning and Geostatistical Methods for Gap Filling and Predicting Solar-Induced Fluorescence Values** - explicitly compares and combines machine-learning and geostatistical SIF gap-filling/prediction methods.
14. **A Spatiotemporal Constrained Machine Learning Method for OCO-2 Solar-Induced Chlorophyll Fluorescence (SIF) Reconstruction** - reconstructs contiguous OCO-2 SIF using spatiotemporally constrained LightGBM.
15. **A Convolutional Neural Network for Spatial Downscaling of Satellite-Based Solar-Induced Chlorophyll Fluorescence (SIFnet)** - CNN-based TROPOMI SIF downscaling to 500 m.
16. **An Improved Spatially Downscaled Solar-Induced Chlorophyll Fluorescence Dataset from the TROPOMI Product (TroDSIF)** - redistributes coarse TROPOMI SIF to a 500 m, 16-day product.
17. **High-Resolution Global Contiguous SIF of OCO-2 (MODIS-SIF)** - machine-learning reconstruction of globally contiguous OCO-2 SIF at 0.05-degree, 16-day resolution.
18. **Generating High-Resolution Total Canopy SIF Emission from TROPOMI Data: Algorithm and Application (HSIF)** - derives 1 km observed and total-canopy SIF from coarse TROPOMI observations.
19. **A Global Spatially Contiguous Solar-Induced Fluorescence (CSIF) Dataset Using Neural Networks** - predicts globally contiguous OCO-2-based SIF from MODIS reflectance using a neural network.
20. **CNSIF: A Reconstructed Monthly 500-Meter Spatial Resolution Solar-Induced Chlorophyll Fluorescence Dataset in China** - deep-learning reconstruction using Landsat/Sentinel-2 reflectance and thermal data.
21. **Spatial Statistical Prediction of Solar-Induced Chlorophyll Fluorescence (SIF) from Multivariate OCO-2 Data** - produces a denoised and gap-filled Level-3 OCO-2 SIF product using cokriging.

### Detailed extraction - papers 1-21

The performance evidence below distinguishes three forms of evaluation. **Internal evaluation**
uses held-out observations from the same SIF source that supplied the model target. **External SIF
validation** uses SIF observations from a different instrument or ground system. **Proxy/application
validation** compares the enhanced SIF with variables such as GPP, yield, or drought indicators;
these comparisons demonstrate usefulness but are not direct measurements of SIF prediction error.

#### 1. GOSIF

- **Reference/title:** *A Global, 0.05-Degree Product of Solar-Induced Chlorophyll Fluorescence Derived from OCO-2, MODIS, and Reanalysis Data*.
- **Product/acronym:** GOSIF (global OCO-2 SIF).
- **Enhancement type:** Data-driven reconstruction and gap filling without an aggregation constraint. The model transfers the relationship learned at sampled OCO-2 locations to a global, spatially continuous grid and retrospectively extends the record before the OCO-2 period.
- **Supervisory SIF source:** Quality-filtered OCO-2 SIF at 757 nm, aggregated to 0.05-degree, 8-day cells; cells required more than five soundings. Nadir-mode observations were used for model development.
- **Predictors:** MODIS EVI plus MERRA-2 photosynthetically active radiation (PAR), vapor-pressure deficit (VPD), and air temperature. Land-cover type was tested but omitted from the final model because it added little predictive benefit and could introduce classification uncertainty.
- **Model family:** Cubist regression trees (rule-based model trees with multivariate linear regressions at terminal rules).
- **Output resolution/coverage:** Global; 0.05 degree; 8-day; 2000-2017.
- **Training/testing design and performance:** Approximately 1.3 million samples from each of 2015 and 2016 formed the training pool, and approximately 1.1 million samples from 2017 formed a temporally held-out test set. The selected model used half of the available training samples and obtained test R-squared = 0.79 and RMSE = 0.07 W m^-2 micrometre^-1 sr^-1. Performance was above R-squared = 0.75 for several productive biomes, including croplands, but weaker for evergreen broadleaf forest (R-squared = 0.43) and open shrubland (R-squared = 0.46).
- **External validation:** No direct independent SIF validation was reported. As proxy validation, GOSIF was compared with FLUXNET2015 eddy-covariance GPP from 91 sites and obtained R-squared = 0.73 (p < 0.001). It was also compared with MODIS GPP and coarse OCO-2 aggregates for spatial and seasonal consistency.
- **Thesis-use note:** GOSIF is a close conceptual baseline for the thesis because it predicts sparse OCO-2 SIF from spatially complete optical and meteorological covariates. The important contrast is scale: GOSIF learns and predicts on a 0.05-degree grid, whereas the thesis uses Sentinel-2 information and explicitly studies aggregation at much finer scales. The biome-specific degradation also motivates agro-climatic and land-cover-aware evaluation rather than relying only on a global score.

#### 2. Random-forest downscaling of GOSIF in Henan Province

- **Reference/title:** *A Spatial Downscaling Method for Solar-Induced Chlorophyll Fluorescence Product Using Random Forest Regression and Drought Monitoring in Henan Province*.
- **Product/acronym:** No official SIF-product acronym stated; described here as RF-downscaled GOSIF. The SIF anomaly index is a downstream drought indicator, not the product name.
- **Enhancement type:** Regional statistical downscaling of GOSIF from 0.05 degree to 1 km.
- **Supervisory SIF source:** Monthly GOSIF, which is itself an OCO-2-based reconstructed product.
- **Predictors:** MODIS monthly NDVI (MOD13A3) and daytime LST (MOD11A2). The 1 km predictors were first resampled to the 0.05-degree GOSIF grid to train the coarse-scale relationship, then supplied at 1 km for prediction.
- **Model family:** Random-forest regression.
- **Output resolution/coverage:** Henan Province, China; 1 km; monthly; 2001-2020. May 2018 was used for the paper's principal downscaling comparison.
- **Training/testing design and performance:** NA - the paper does not report a conventional held-out train/test split or RMSE/R-squared against unseen SIF. After the 1 km predictions were aggregated back to 0.05 degree, their correlation with GOSIF was R = 0.74; this is an aggregation-consistency check against the supervisory product, not an independent accuracy test.
- **External validation:** No independent SIF observations were used. Proxy/application validation showed a May 2018 correlation of R = 0.74 between downscaled SIF and MODIS GPP, compared with R = 0.68 for original GOSIF. The SIF anomaly index correlated with wheat yield at R = 0.93 and maize yield at R = 0.89, while its correlation with drought-affected area was R = -0.58.
- **Thesis-use note:** This is directly relevant to the thesis's agricultural objective and demonstrates the value of higher-resolution spatial detail. However, it learns from an already reconstructed coarse product and lacks a genuinely held-out or independent SIF evaluation. Its very high yield correlations should therefore be presented as application evidence, not proof that the 1 km SIF values are physically accurate.

#### 3. Global downscaled GOME-2 SIF (Duveiller et al.)

- **Reference/title:** *A Spatially Downscaled Sun-Induced Fluorescence Global Product for Enhanced Monitoring of Vegetation Productivity*.
- **Product/acronym:** No single formal acronym is introduced in the paper; the product is a downscaled GOME-2 SIF dataset with separate Joiner et al. (JJ) and Koehler et al. (PK) retrieval variants.
- **Enhancement type:** Aggregation-constrained spatial downscaling. Fine-grid SIF is spatially redistributed while preserving the original coarse GOME-2 signal.
- **Supervisory SIF source:** Two GOME-2 retrievals (JJ and PK), gridded to 0.5 degree and 16-day composites.
- **Predictors:** The selected semi-empirical configuration uses MODIS NIRv as the vegetation/absorbed-light proxy, NDWI as the water-stress proxy, and afternoon LST as the temperature-stress proxy. Earlier configurations using NDVI, evapotranspiration, and LST were also evaluated.
- **Model family:** A locally calibrated, semi-empirical light-use-efficiency model rather than a generic machine-learning regressor. It retains the coarse SIF signal and uses fine-resolution remote-sensing drivers to allocate it among subpixels.
- **Output resolution/coverage:** Global; 0.05 degree; 8-day; 2007-2018.
- **Training/testing design and performance:** NA - there is no standard random train/test split against GOME-2. Independent OCO-2 observations from 2015-2017 were used to select and benchmark configurations. For the new PK configuration, the reported comparison was absolute bias = 0.0615, r = 0.8134, and agreement index lambda = 0.7208; the new JJ configuration obtained absolute bias = 0.063, r = 0.8245, and lambda = 0.7325.
- **External validation:** OCO-2 is an independent SIF validation source here, although it was also used for configuration selection. The product was additionally compared with TROPOMI during their 2018 overlap. Pixel-wise bias correction reduced the root-mean-square deviation of median latitudinal profiles from 0.188 to 0.0328 mW m^-2 sr^-1 nm^-1.
- **Thesis-use note:** This paper provides a useful contrast between unconstrained prediction and mass/aggregation-preserving downscaling. Its constraint prevents fine-grid estimates from changing the parent-cell mean, an idea closely related to the thesis's aggregation experiments. The paper also shows that cross-sensor bias correction may be necessary even when spatial and temporal patterns agree.

#### 4. BCSIF

- **Reference/title:** *A Spatially Downscaled TROPOMI SIF Product at 0.005-Degree Resolution With Bias Correction*.
- **Product/acronym:** BCSIF (bias-corrected downscaled SIF); BCSIFas denotes the all-sky version.
- **Enhancement type:** TROPOMI SIF downscaling with an explicit coarse-scale prediction-bias correction. The correction is designed to preserve information in the original TROPOMI observations that a universal statistical model fails to reproduce.
- **Supervisory SIF source:** TROPOMI L2B clear-sky SIF, aggregated to 0.05 degree and 16-day composites for model development.
- **Predictors:** MODIS reflectance bands 1-4, NDWI, cosine of solar zenith angle, and ERA5 air temperature, VPD, and soil moisture. PAR variables were considered in data preparation, while the stated final RF input set contained the nine variables listed above.
- **Model family:** Random-forest regression followed by spatial bias correction derived from the difference between observed and RF-predicted SIF at the coarse scale.
- **Output resolution/coverage:** China; 0.005 degree (approximately 500 m); 16-day; 2019-2020.
- **Training/testing design and performance:** Three repeated 70/30 train/test experiments were used for each 16-day RF model. Mean test R-squared was 0.887 and mean RMSE was 0.0716 mW m^-2 nm^-1 sr^-1. When the high-resolution products were reaggregated to 0.05 degree, bias correction improved agreement with TROPOMI from R-squared = 0.853 to 0.882 and reduced RMSE from 0.0898 to 0.0810 mW m^-2 nm^-1 sr^-1.
- **External validation:** Direct independent SIF validation used ChinaSpec tower observations at the DaMan maize and Arou alpine-meadow sites. Mean R-squared increased from 0.590 for 0.05-degree TROPOMI SIF to 0.798 for 0.005-degree BCSIF, while mean RMSE decreased from 0.201 to 0.141 mW m^-2 nm^-1 sr^-1. As proxy validation, mean biome-level spatial-temporal R-squared with FLUXCOM GPP increased from 0.472 for the original TROPOMI product to 0.877 for BCSIF; GLASS GPP was used for spatial-pattern comparison.
- **Thesis-use note:** Among the first five papers, BCSIF supplies the strongest direct evidence that reducing the spatial mismatch can improve agreement with ground SIF. Its coarse-scale residual/bias-preservation step is particularly relevant to the thesis: a high-capacity model can generate realistic spatial detail while still missing information present in the supervising coarse observation.

#### 5. CNN downscaling of GOSIF for agricultural drought

- **Reference/title:** *Downscaling Solar-Induced Chlorophyll Fluorescence Based on Convolutional Neural Network Method to Monitor Agricultural Drought*.
- **Product/acronym:** No official SIF-product acronym stated. TFDI means temperature fluorescence dryness index and is the downstream drought index derived from the downscaled SIF, not the SIF product itself.
- **Enhancement type:** Regional CNN-based spatial downscaling of GOSIF from 0.05 degree to 0.008 degree.
- **Supervisory SIF source:** GOSIF, an OCO-2-based reconstructed product, aggregated to monthly values for this application.
- **Predictors:** MODIS NDVI, EVI, daytime LST, and nighttime LST. The 1 km predictor fields were resampled to 0.05 degree to learn the relationship with GOSIF and then used at their native fine scale for downscaling. Other datasets such as precipitation, soil moisture, elevation, and land use supported drought analysis rather than the core SIF prediction model.
- **Model family:** CNN with three 3 x 3 convolutional layers, intended to learn nonlinear relationships and spatial dependencies among adjacent pixels.
- **Output resolution/coverage:** Henan Province, China; 0.008 degree (approximately 1 km); monthly; June-October 2013-2017.
- **Training/testing design and performance:** NA - no explicit held-out split, R-squared, RMSE, or MAE for unseen SIF is reported. When the downscaled SIF was aggregated back to the 0.05-degree GOSIF grid, the maximum reported correlation with GOSIF was R = 0.9174. This measures coarse-scale consistency rather than independent fine-scale accuracy.
- **External validation:** No independent SIF dataset was used. Proxy validation against MOD17A2H GPP produced a maximum monthly correlation of R = 0.8346 in June 2017. For the downstream drought application, annual mean TFDI correlated with Henan summer-corn yield at R = -0.84.
- **Thesis-use note:** This paper supports the use of spatial-context models for heterogeneous agricultural landscapes, but it does not establish that its 1 km estimates are correct at 1 km. That limitation is important for positioning the thesis's footprint-level validation and aggregation analysis. The paper's use of GOSIF as supervision also means that OCO-2 information has already passed through an earlier reconstruction model before the CNN downscaling stage.

#### 6. DSIFRFK\*_EA0.05

- **Reference/title:** *Downscaling Solar-Induced Chlorophyll Fluorescence to a 0.05-Degree Monthly Product Using AVHRR Data in East Asia (1995-2003)*.
- **Product/acronym:** DSIFRFK\*_EA0.05.
- **Enhancement type:** Regional statistical downscaling and historical reconstruction. A random-forest trend is predicted at fine resolution and coarse-grid residuals are redistributed by kriging.
- **Supervisory SIF source:** Monthly Level-3 GOME/ERS-2 SIF from channel 4 (734-758 nm, referenced at 740 nm), modeled on a 2-degree grid. The original GOME footprint was approximately 40 x 320 km.
- **Predictors:** AVHRR and ERA5 variables. From 15 initial variables, the selected six were skin temperature, downward shortwave radiation, NIRv, NDVI, near-infrared reflectance, and total precipitation.
- **Model family:** Random forest kriging (RFK): random-forest regression plus ordinary kriging of residuals. RF alone and RF with inverse-distance-weighted residual interpolation were comparison methods. The global RF configuration was selected instead of separate monthly models.
- **Output resolution/coverage:** East Asia (72 degrees E-150 degrees E and 16 degrees N-54 degrees N); 0.05 degree; monthly; July 1995-June 2003.
- **Training/testing design and performance:** RF modeling was repeated with random 80/20 training/test splits. The selected downscaled product had the highest reported mean R-squared of 0.83 and the lowest RMSE of 0.08 mW m^-2 nm^-1 sr^-1 when compared with original GOME SIF.
- **External validation:** Comparisons with other SIF products produced correlations of r = 0.73 with 1-degree SCIAMACHY SIF, r = 0.88 with 0.05-degree downscaled SCIAMACHY SIF, and r = 0.89 with the 0.05-degree LT_SIFc* product. These are cross-product consistency checks rather than fine-scale ground truth. Proxy validation against GPP from eight flux sites gave R-squared = 0.73, with site-specific R-squared values from 0.61 to 0.99.
- **Thesis-use note:** The combination of a learned trend and spatially redistributed residuals is directly relevant to aggregation-preserving SIF enhancement. Importantly, the authors acknowledge that ordinary kriging does not account for the change in spatial support between 2-degree and 0.05-degree data, explicitly identifying the modifiable areal unit problem. This supports the thesis's focus on how aggregation and observation support affect model evaluation.

#### 7. Fine-grained OCO-2 SIF forecasting

- **Reference/title:** *Forecasting of Fine-Grained SIF of OCO-2 Using Multi-source Data and AI-Based Techniques*.
- **Product/acronym:** NA - no formal product acronym is stated.
- **Enhancement type:** Two-stage gap filling and temporal forecasting. Random forest is first used to construct spatially complete SIF fields needed as lagged inputs; machine-learning models then forecast future monthly and seasonal SIF.
- **Supervisory SIF source:** OCO-2 SIF observations with an approximately 1.29 x 2.25 km footprint, using data from 2017-2019.
- **Predictors:** Two preceding SIF values (t-1 and t-2), MODIS EVI, and spatial-temporal features. EVI was selected after comparison with NDVI, NIRv, blue and red reflectance, and land-surface temperature.
- **Model family:** Random forest for the initial gap-filled SIF series. Lasso, Ridge, random forest, LightGBM, CNN, and LSTM were compared for forecasting; Lasso was best for the monthly experiment and LightGBM for the seasonal experiment.
- **Output resolution/coverage:** Karnataka and Maharashtra, India; approximately 1 x 1 km; monthly and seasonal forecasts; 2019 forecast period.
- **Training/testing design and performance:** Models were trained using 2017-2018 and tested by forecasting 2019. The headline monthly Lasso result was RMSE = 0.0355 W m^-2 nm^-1 sr^-1 and MAPE = 16.9093%. The seasonal LightGBM result was RMSE = 0.0389 W m^-2 nm^-1 sr^-1 and MAPE = 17.4895%. A state-level EVI/LightGBM experiment reported RMSE = 0.0561 and MAPE = 18.3890% for Karnataka, and RMSE = 0.0594 and MAPE = 19.2179% for Maharashtra.
- **External validation:** NA - no independent SIF instrument or ground SIF observations were used. Forecasts were compared with an OCO-2-derived 0.05-degree SIF product and NDVI, so these provide same-source/product and proxy comparisons rather than independent SIF validation.
- **Thesis-use note:** Unlike most Section 3.1 papers, this work predicts future SIF rather than only reconstructing SIF for the same date. Its use of lagged, initially gap-filled SIF introduces serial information but also means first-stage prediction error can propagate into the forecast. Its nearest-pixel data fusion also contrasts with the thesis's explicit averaging of Sentinel-2 pixels within OCO-2 footprints.

#### 8. Spatially continuous TanSat SIF

- **Reference/title:** *Generation of a Global Spatially Continuous TanSat Solar-Induced Chlorophyll Fluorescence Product by Considering the Impact of the Solar Radiation Intensity*.
- **Product/acronym:** No unambiguous formal acronym is introduced; the product is described as continuous TanSat SIF and is also called TanSIF in the paper.
- **Enhancement type:** Global statistical reconstruction and gap filling without an explicit aggregation-preservation constraint.
- **Supervisory SIF source:** TanSat nadir SIF at 758 nm, with an approximately 2 x 2 km footprint, from February 2017 to August 2019.
- **Predictors:** MODIS BRDF-corrected red, near-infrared, blue, and green reflectance; NDVI; cosine of the noon solar zenith angle; and air temperature.
- **Model family:** Random-forest regression applied to normalized apparent SIF yield. Including solar-radiation geometry was a central methodological contribution.
- **Output resolution/coverage:** Global; 0.05 degree; 4-day; 2017-2019.
- **Training/testing design and performance:** For each year, observations were randomly divided into 70% training and 30% validation data. Validation R-squared values for 2017, 2018, and 2019 were 0.75, 0.73, and 0.81, with corresponding RMSE values of 0.32, 0.30, and 0.32 mW m^-2 nm^-1 sr^-1. For 2018, adding the cosine of solar zenith angle improved R-squared from 0.65 to approximately 0.72-0.73 and reduced RMSE from 0.34 to 0.30 mW m^-2 nm^-1 sr^-1. The abstract and main text differ slightly on whether the latter R-squared is 0.72 or 0.73.
- **External validation:** Independent cross-sensor comparison against TROPOMI SIF produced R-squared = 0.73 after spatial matching to a 0.2-degree grid. No ground SIF validation was reported.
- **Thesis-use note:** This study demonstrates that illumination geometry can materially improve a reflectance-based SIF reconstruction and should be considered when observations span different sun-sensor conditions. Its random observation split may nevertheless overstate transfer performance relative to spatially or temporally separated tests.

#### 9. DOSIF

- **Reference/title:** *DOSIF: Long-Term Daily SIF from OCO-3 with Global Contiguous Coverage*.
- **Product/acronym:** DOSIF (daily OCO-3 SIF).
- **Enhancement type:** Long-term daily reconstruction and gap filling using day-of-year-specific and biome-specific sampling.
- **Supervisory SIF source:** Quality-controlled OCO-3 SIF from 2019-2023, with a native footprint of approximately 1.3 x 2.25 km. The processing combines the 757 and 771 nm bands and applies five-nearest-neighbor smoothing to reduce observation noise.
- **Predictors:** Seven daily MODIS MCD43A4 BRDF-corrected surface-reflectance bands. MCD12C1 land cover and location information define vegetated masks and spatial sub-biomes rather than serving only as ordinary continuous predictors.
- **Model family:** Moving Spatial-Temporal Window Sampling (MSTWS) with CatBoost. ANN, random forest, and CatBoost were benchmarked; CatBoost was selected for stronger validation robustness. For each target day of year, samples came from a 16-day moving window across available OCO-3 years and were divided into 19 sub-biome models.
- **Output resolution/coverage:** Global vegetated land; 0.05 degree; daily; 2001-present. The generated dataset described in the paper extends through 2024.
- **Training/testing design and performance:** Samples were split 70/30 for training and validation, with five-fold cross-validation for hyperparameter tuning. Overall CatBoost performance was training R-squared = 0.92, RMSE = 0.05 W m^-2 micrometre^-1 sr^-1, slope = 0.94; validation R-squared = 0.81, RMSE = 0.07 W m^-2 micrometre^-1 sr^-1, slope = 0.91. On the same reserved comparison set, MSTWS obtained R-squared = 0.85, RMSE = 0.07, and slope = 0.93, compared with R-squared = 0.81, RMSE = 0.08, and slope = 0.87 for a universal model. The conclusion separately states MSTWS R-squared = 0.92 versus 0.85 and an independent-validation R-squared of 0.88, but these values do not match the detailed Results section; the detailed Results values should therefore be preferred when quoting performance.
- **External validation:** Independent airborne CFIS observations from 2016 were aggregated to the DOSIF 0.05-degree grid and divided by 1.12 to account for the 755 versus 757 nm wavelength difference. The comparison produced R-squared = 0.75, RMSE = 0.11 W m^-2 micrometre^-1 sr^-1, and slope = 0.96.
- **Thesis-use note:** DOSIF provides strong evidence for context-specific rather than universal SIF-reflectance relationships: day-of-year and sub-biome sampling improved prediction over a single global model. It also provides genuine independent SIF validation, although the CFIS wavelength correction, temporal separation from the OCO-3 training period, and aggregation to 0.05 degree should be stated. The result does not demonstrate accuracy at individual fine-resolution pixels.

#### 10. ANN-based TROPOMI SIF retrieval

- **Reference/title:** *A More Precise Retrieval of Sun-Induced Chlorophyll Fluorescence from Satellite Data Using Artificial Neural Networks*.
- **Product/acronym:** NA - no formal product acronym is stated; described as ANN-based TROPOMI SIF.
- **Enhancement type:** Retrieval-precision enhancement rather than spatial downscaling or gap filling. The method estimates SIF directly from TROPOMI radiance while reducing the noise of individual retrievals, thereby avoiding the spatial or temporal aggregation normally used for noise reduction.
- **Supervisory SIF source:** OCO-2/3 SIF at 757 nm, converted to 740 nm by a factor of 1.5. At least five time- and geometry-matched OCO-2/3 observations inside each TROPOMI footprint were averaged to form a low-noise reference SIF target.
- **Predictors:** TROPOMI Level-1B top-of-atmosphere radiance across 743-758 nm. The models were trained separately for detector columns because their spectral responses differ.
- **Model family:** Feed-forward ANN with one hidden layer and eight neurons; separate models were fitted for TROPOMI columns 26-422.
- **Output resolution/coverage:** Global; daily; original TROPOMI footprints of approximately 7 x 3.5 km before August 2019 and 5.5 x 3.5 km thereafter; May 2018-December 2024.
- **Training/testing design and performance:** Optimized matched samples were randomly divided 50/50 into calibration and validation sets; the calibration half was further divided 60/20/20 for training, validation, and testing during model fitting. Against the held-out matched OCO-2/3 reference, combined-column performance was R-squared = 0.85 and RMSE = 0.217 mW m^-2 nm^-1 sr^-1. SCOPE-MODTRAN simulations produced R-squared = 0.96 and RMSE = 0.304 mW m^-2 nm^-1 sr^-1. The assumed noise of a single ANN retrieval was approximately half the reported noise of SVD-based TROPOMI products.
- **External validation:** NA for direct independent validation of the ANN product. The paper notes that OCO-2 SIF has previously been validated using airborne measurements, but those measurements were not a separate direct test of the final ANN-TROPOMI product in this study. Proxy validation used daily GPP from 96 AmeriFlux/FLUXNET sites: vegetation-type R-squared ranged from 0.45 to 0.79 for ANN SIF, compared with 0.12-0.66 for ESA SVD SIF and 0.09-0.65 for Caltech SVD SIF. Near-zero SIF and lower dispersion at ten non-vegetated pseudo-invariant calibration sites provided an additional retrieval-noise check.
- **Thesis-use note:** This paper belongs in a distinct retrieval-enhancement category. It improves the SIF observation supplied to later reconstruction models but does not fill spatial gaps. Its strategy of aggregating multiple accurate OCO-2/3 retrievals within a coarser TROPOMI footprint to create low-noise labels is conceptually relevant to the thesis's footprint-level aggregation, while also illustrating the trade-off between label precision and spatial support.

#### 11. RTSIF

- **Reference/title:** *A Long-Term Reconstructed TROPOMI Solar-Induced Fluorescence Dataset Using Machine Learning Algorithms*.
- **Product/acronym:** RTSIF (reconstructed TROPOMI SIF).
- **Enhancement type:** Global reconstruction, gap filling, and temporal extension of the short TROPOMI record to 2001.
- **Supervisory SIF source:** Clear-sky Caltech TROPOMI SIF at 740 nm from March 2018 to December 2020, aggregated to 0.05-degree, 8-day cells containing more than four valid retrievals.
- **Predictors:** Seven MODIS BRDF-corrected NBAR reflectance bands, MODIS land-surface temperature, MODIS land cover, CERES PAR, and ISLSCP II C3/C4 vegetation fractions.
- **Model family:** XGBoost gradient-boosted decision trees.
- **Output resolution/coverage:** Global clear-sky land; 0.05 degree; 8-day; 2001-2020.
- **Training/testing design and performance:** Data were divided into 80% training and 20% testing, with a 10-fold cross-validated grid search for hyperparameters. Training performance was R-squared = 0.916 and RMSE = 0.059 mW m^-2 nm^-1 sr^-1. Testing performance was R-squared = 0.907, RMSE = 0.062 mW m^-2 nm^-1 sr^-1, and regression slope = 1.001.
- **External validation:** Direct tower-SIF validation used PhotoSpec observations at Southern Old Black Spruce and Niwot Ridge. Agreement was R-squared = 0.754 and 0.84, respectively. Cross-sensor comparisons with OCO-2 SIF (2015-2020) and GOME-2 SIF (2007-2019), aggregated to 1 degree and 8 days, generally produced regional R-squared values above 0.7, although no single global score was stated. Proxy evaluation used FLUXNET2015 GPP from 76 homogeneous sites; a relationship was reported at 8-day and annual scales, but a numerical performance value was not clearly stated.
- **Thesis-use note:** RTSIF shows that a high-quality, short-duration SIF record can supervise a much longer reconstruction when the predictor archive extends farther back. Its independent tower-SIF evaluation is valuable, but the two tower sites are insufficient to establish universal fine-grid accuracy. The model is trained and evaluated at the same 0.05-degree support, unlike the thesis's explicit use of much finer Sentinel-2 pixels aggregated within OCO-2 footprints.

#### 12. ROSIF and the photosynthetic afternoon-depression study

- **Reference/title:** *Tracking Drought in Dryland Vegetation through the Photosynthetic Afternoon Depression Index of Sun-Induced Chlorophyll Fluorescence*.
- **Product/acronym:** ROSIF (reconstructed OCO-3 SIF).
- **Enhancement type:** Regional, sub-daily reconstruction and spatial gap filling of OCO-3 SIF, followed by derivation of an afternoon-depression drought indicator.
- **Supervisory SIF source:** Quality-filtered OCO-3 version 10 nadir SIF from 2019-2020, separated into observations around 10:00 and 14:00 local time.
- **Predictors:** The eight GOCI visible-to-near-infrared reflectance bands at the corresponding observation time. ERA5-Land PAR, VPD, precipitation, land-surface temperature, and soil moisture supported SIF-yield and drought analysis but were not inputs to the stated core SIF reconstruction equation.
- **Model family:** XGBoost and random forest were compared, with XGBoost selected for the final ROSIF reconstruction. Bayesian optimization was used for hyperparameters.
- **Output resolution/coverage:** Liaoning Province, China; 500 m; two sub-daily snapshots at approximately 10:00 and 14:00 local time; 2019-2020.
- **Training/testing design and performance:** Ten-fold cross-validation was used. At 10:00, XGBoost obtained R-squared = 0.684, RMSE = 0.388, and MAE = 0.290 mW m^-2 sr^-1 nm^-1; random forest obtained R-squared = 0.683, RMSE = 0.389, and MAE = 0.291. At 14:00, XGBoost obtained R-squared = 0.598, RMSE = 0.348, and MAE = 0.259; random forest obtained R-squared = 0.592, RMSE = 0.351, and MAE = 0.262 in the same units.
- **External validation:** NA for a quantitative independent fine-resolution SIF accuracy test. ROSIF was compared qualitatively with the coarser reconstructed GOSIF and RTSIF products. Ground SIF measurements for corn and soybean from 2017 were used to examine seasonal SIF-yield behavior, but they were not reported as a direct collocated validation of the 2019-2020 ROSIF estimates. As application validation, the derived afternoon-depression SIF-yield index correlated with SPEI at r = 0.71 (p < 0.01) and with the soil-moisture anomaly indicator SMZ at r = 0.53 (p < 0.05).
- **Thesis-use note:** The study is useful because it demonstrates sub-daily reconstruction and shows how a derived SIF signal can reveal drought information hidden in a single-time observation. For the thesis, it also illustrates that environmental variables used for downstream interpretation are not necessarily model predictors. The reported cross-validation only establishes agreement with held-out OCO-3 samples, not the accuracy of every 500 m pixel.

#### 13. Hybrid kriging with external drift

- **Reference/title:** *Hybrid Machine Learning and Geostatistical Methods for Gap Filling and Predicting Solar-Induced Fluorescence Values*.
- **Product/acronym:** NA - no formal product acronym is stated. The proposed method is kriging with external drift (KED).
- **Enhancement type:** Global SIF gap filling that combines data-driven SIF estimates with the spatial covariance structure of observed SIF.
- **Supervisory SIF source:** Clear-sky OCO-2 Level-2 SIF at 757 nm for 2019, mean-aggregated to 0.05-degree cells; cells with fewer than five retrievals were excluded.
- **Predictors:** The external-drift covariate is the existing 0.05-degree, 4-day CSIF neural-network product. CSIF itself was generated from the first four MODIS NBAR reflectance bands, but this paper uses the resulting CSIF estimate rather than retraining the original neural network. Nearby OCO-2 observations and their modeled covariance provide the geostatistical information.
- **Model family:** Moving-window ordinary kriging, CSIF machine-learning prediction, and their hybrid through kriging with external drift.
- **Output resolution/coverage:** Global; 0.05 degree; 4-day support inherited from CSIF; 2019 method-evaluation period.
- **Training/testing design and performance:** Leave-one-out cross-validation was performed on more than 400,000 OCO-2 aggregate samples. The hybrid method obtained MAE = 0.1183, RMSE = 0.1556, R-squared = 0.8523, and bias = -0.0002 mW m^-2 sr^-1 nm^-1. Ordinary kriging obtained MAE = 0.1318, RMSE = 0.1752, R-squared = 0.8129, and bias = -0.0003. The standalone ML/CSIF estimate obtained MAE = 0.1399, RMSE = 0.1823, R-squared = 0.8004, and bias = 0.0103.
- **External validation:** NA - all three approaches were evaluated against withheld aggregates from the same OCO-2 dataset. No independent satellite, airborne, or tower SIF observations were used.
- **Thesis-use note:** This is the clearest paper in the batch for comparing two complementary information sources: ancillary-variable relationships and spatial covariance among SIF observations. Its improvement over both parent methods suggests a useful future extension for the thesis. However, leave-one-out validation may be easier than predicting large continuous gaps because nearby OCO-2 observations remain available for most held-out points.

#### 14. ST-LGBM OCO-2 SIF reconstruction

- **Reference/title:** *A Spatiotemporal Constrained Machine Learning Method for OCO-2 Solar-Induced Chlorophyll Fluorescence (SIF) Reconstruction*.
- **Product/acronym:** ST-LGBM (spatiotemporally constrained LightGBM) is the method acronym; no separate product acronym is clearly stated.
- **Enhancement type:** Global reconstruction and gap filling designed specifically to improve extrapolation into regions between sparse OCO-2 swaths.
- **Supervisory SIF source:** Quality-filtered OCO-2 version 10r daily-corrected nadir SIF at 757 nm. Observations were aggregated to 0.05-degree, 8-day cells when more than five footprints were available. Data from 2015-2018 were used in model development.
- **Predictors:** MODIS NIRv, MERRA-2 VPD, air temperature and PAR, MODIS land cover, plus two derived constraints: a spatial SIF factor based on geographically nearby, NIRv-similar OCO-2 observations and a temporal factor based on the same phenological period in other years.
- **Model family:** LightGBM with spatial and temporal SIF constraints (ST-LGBM). Ordinary LightGBM was the main baseline; DBN, ANN, and Cubist models were also used in sensitivity experiments.
- **Output resolution/coverage:** Global; 0.05 degree; 8-day; September 2014-December 2019.
- **Training/testing design and performance:** Two validation designs were reported. Random ten-fold validation used 90% training and 10% testing per fold: ST-LGBM obtained training R-squared = 0.84 and RMSE = 0.06 W m^-2 micrometre^-1 sr^-1, and test R-squared = 0.83 and RMSE = 0.07. Ordinary LightGBM obtained training/test R-squared = 0.82/0.81 and test RMSE = 0.07. A more demanding swath-held-out experiment simulated regions without training observations: ST-LGBM achieved test R-squared = 0.79, compared with 0.74 for ordinary LightGBM. The ST-LGBM decrease from random to spatial-gap testing was 0.04, versus 0.07 for LightGBM.
- **External validation:** Independent cross-sensor comparison used 2018 TROPOMI SIF at 740 nm. For annual averages, R-squared with TROPOMI increased from 0.89 for LightGBM to 0.91 for ST-LGBM. Across 8-day scenes, ST-LGBM improved R-squared relative to LightGBM by 0.04 on average and relative to original OCO-2 observations by 0.13 on average. The comparison is affected by the 740 versus 757 nm wavelength difference and independent retrieval errors.
- **Thesis-use note:** This paper provides especially strong methodological support for spatially separated validation. It demonstrates that random sample splits can overestimate performance in genuinely unobserved regions, which closely matches the thesis's need to evaluate gap filling beyond sampled OCO-2 footprints. The improvement from spatiotemporal SIF neighbors also suggests that Sentinel-2-only prediction could later be extended with structured spatial and temporal context.

#### 15. SIFnet

- **Reference/title:** *A Convolutional Neural Network for Spatial Downscaling of Satellite-Based Solar-Induced Chlorophyll Fluorescence (SIFnet)*.
- **Product/acronym:** SIFnet.
- **Enhancement type:** CNN-based spatial downscaling of TROPOMI SIF by a factor of ten. Coarse SIF is retained as a model input, but the stated MSE-plus-structural-similarity loss does not explicitly force fine-pixel means to equal the parent coarse observation.
- **Supervisory SIF source:** TROPOMI SIF at 740 nm. During scale-factor training, 0.05-degree TROPOMI SIF was the target and SIF coarsened to 0.5 degree was an input. The optimized network was then supplied with 0.05-degree TROPOMI SIF to estimate 0.005-degree SIF.
- **Predictors:** Coarse TROPOMI SIF; seven MODIS reflectance bands; NIRv, kNDVI, NDVI, and EVI; cosine of solar zenith angle; ERA5-Land temperature and precipitation at the current and previous 16-day period; SMAP surface and subsurface soil moisture; elevation; fractional land-cover classes; forest share; and forest-edge share. Coarse SIF was the most important feature, followed by NIRv and cosine of solar zenith angle.
- **Model family:** Residual convolutional neural network with convolutional and ReLU layers, optimized using a combined mean-squared-error and structural-dissimilarity loss.
- **Output resolution/coverage:** Conterminous United States; 0.005 degree (approximately 500 m); 16-day; April 2018-March 2021.
- **Training/testing design and performance:** Five non-US geographic folds across Asia, Europe, southern Africa, and South America were used for training, while North America was held out for validation and hyperparameter selection. At the emulated 0.05-degree target scale, SIFnet obtained R-squared = 0.92, RMSE = 0.17 mW m^-2 sr^-1 nm^-1, and structural similarity index = 0.87 against TROPOMI. This is a geographically separated same-sensor test and demonstrates continental transfer, but it evaluates the learned tenfold scaling relationship at 0.05 degree rather than directly observing truth at 0.005 degree.
- **External validation:** Independent OCO-2 observations from April 2018-March 2021 and OCO-3 observations from July 2019-March 2021 were used over CONUS. Fine SIFnet cells were averaged within each OCO footprint. Against OCO-2, SIFnet obtained Pearson r = 0.77, R-squared = 0.57, and RMSE = 0.21 mW m^-2 sr^-1 nm^-1. Against OCO-3, it obtained r = 0.79, R-squared = 0.62, and RMSE = 0.19 mW m^-2 sr^-1 nm^-1. Against the combined OCO-2 and OCO-3 observations, it obtained r = 0.78, R-squared = 0.59, and RMSE = 0.20 mW m^-2 sr^-1 nm^-1.
- **Thesis-use note:** SIFnet is a close deep-learning analogue for learning spatial detail from coarse SIF plus fine auxiliary data. Its geographic holdout and OCO-2/3 footprint validation are major strengths. The method's averaging of fine predictions within OCO footprints is directly comparable to the thesis evaluation design, while its lack of an explicit coarse-mean conservation constraint leaves room for the thesis's aggregation-focused analysis.

#### 16. TroDSIF

- **Reference/title:** *An Improved Spatially Downscaled Solar-Induced Chlorophyll Fluorescence Dataset from the TROPOMI Product*.
- **Product/acronym:** TroDSIF.
- **Enhancement type:** Global aggregation-preserving spatial downscaling and gap filling. Random-forest-predicted fine-resolution SIF is used as a spatial weight to redistribute the original TROPOMI signal rather than replacing it directly.
- **Supervisory SIF source:** Caltech daily-corrected far-red TROPOMI SIF at 740 nm, represented as 0.05-degree, 16-day original SIF (OSIF).
- **Predictors:** MODIS red, near-infrared, blue, and green BRDF-corrected reflectance; NDVI; cosine of solar zenith angle; and ERA5 air temperature.
- **Model family:** Random forest for the intermediate predicted SIF (PSIF), followed by weighted redistribution of OSIF using a two-dimensional Gaussian spatial weighting scheme. A Gaussian-weighted fallback is also used where predictor or original-SIF data are missing.
- **Output resolution/coverage:** Global; 500 m; 16-day; April 2018-July 2021.
- **Training/testing design and performance:** Each year's samples were randomly divided into 70% training and 30% validation data. For the reported 2019 model, training R-squared = 0.908 and RMSE = 0.059 mW m^-2 nm^-1 sr^-1; validation R-squared = 0.893 and RMSE = 0.064 mW m^-2 nm^-1 sr^-1. After reaggregation to 0.05 degree, TroDSIF agreed with OSIF at R-squared = 0.948 and RMSE = 0.057 on day of year 14, and R-squared = 0.934 and RMSE = 0.067 on day of year 206 in 2019.
- **External validation:** Direct validation used 16-day tower SIF from six ChinaSpec sites. Site-level TroDSIF RMSE ranged from 0.104 to 0.223 and MAE from 0.077 to 0.163 mW m^-2 nm^-1 sr^-1; a single combined tower R-squared was not stated. Proxy validation against GPP from 67 AmeriFlux sites in 2019 produced an overall R-squared of 0.542 for TroDSIF, compared with 0.483 for original TROPOMI SIF.
- **Thesis-use note:** TroDSIF is one of the closest aggregation-constrained comparators for the thesis. Its fine spatial pattern comes from the learned auxiliary-data relationship, while the original coarse observation remains the controlling measured signal. The approach nevertheless assumes a linear proportional relationship between OSIF and the predicted SIF used for redistribution, and tower comparisons still involve spatial-support differences.

#### 17. SIFoco2_005

- **Reference/title:** *High-Resolution Global Contiguous SIF of OCO-2*.
- **Product/acronym:** SIFoco2_005; the native along-orbit observations are denoted SIFoco2_orbit.
- **Enhancement type:** Global reconstruction and gap filling of sparse OCO-2 observations, with separate prediction models by biome and 16-day time step.
- **Supervisory SIF source:** Quality-controlled OCO-2 nadir observations. SIF at 757 and 771 nm was combined as (SIF757 + 1.5 x SIF771) / 2, converted to a daily mean, and the five nearest valid footprints of the same land-cover class were aggregated for model samples.
- **Predictors:** MODIS MCD43A4/MCD43C4 BRDF-corrected seven-band surface reflectance. MODIS land cover was used to stratify the models by biome and geographic sub-biome.
- **Model family:** Feed-forward multilayer-perceptron ANN. Hidden-layer and neuron counts were optimized separately for each biome-time model through five-fold cross-validation.
- **Output resolution/coverage:** Global; 0.05 degree; 16-day; September 2014-August 2018.
- **Training/testing design and performance:** Five-fold cross-validation withheld 20% of the observations within every sub-biome and time step. Pooled validation performance was R-squared = 0.83, RMSE = 0.065 W m^-2 micrometre^-1 sr^-1, and regression slope = 0.91. In the August 2015 example, croplands, deciduous broadleaf forests, and savannas had R-squared above 0.75 and RMSE at or below 0.075 in the same units.
- **External validation:** Independent 2016 airborne CFIS measurements were scaled by 1/1.12 for the 755 versus 757 nm wavelength difference and aggregated to 0.05 degree. In areas without OCO-2 overpasses, SIFoco2_005 versus CFIS obtained R-squared = 0.72 and slope = 0.96; RMSE was not clearly stated. As a baseline sensor check, native OCO-2 observations under coordinated CFIS flights obtained R-squared = 0.81 and slope = 1.01.
- **Thesis-use note:** The biome- and time-specific design addresses the fact that reflectance-SIF relationships change with vegetation physiology and water stress. Its validation against airborne SIF specifically outside OCO-2 swaths is a major strength. However, the reconstructed output and its independent test are still evaluated at 0.05-degree support, rather than at individual sub-kilometre pixels.

#### 18. HSIF

- **Reference/title:** *Generating High-Resolution Total Canopy SIF Emission from TROPOMI Data: Algorithm and Application*.
- **Product/acronym:** HSIFtotal for high-resolution total canopy SIF emission; HSIFobs for the corresponding high-resolution observed SIF.
- **Enhancement type:** Physically constrained downscaling of TROPOMI observed SIF to 1 km, together with conversion from observed directional SIF to total canopy fluorescence emission.
- **Supervisory SIF source:** Clear-sky TROPOMI SIF retrieved over 743-758 nm and normalized to 740 nm. At least six TROPOMI retrievals were averaged daily within each 0.2-degree parent grid before redistribution.
- **Predictors:** Sentinel-3 OLCI chlorophyll FPAR; fluorescence efficiency derived at coarse scale and assigned by vegetation type and moving window; fluorescence escape probability derived from MODIS MCD19A3 BRDF parameters, NIRv, canopy interception, and leaf/soil reflectance; ERA5 shortwave radiation/PAR; and MODIS land cover.
- **Model family:** Semi-empirical, energy-conservation and light-use-efficiency framework rather than a machine-learning model. It constrains the mean fine-grid emission to the parent TROPOMI observation while allocating spatial variation through FPAR, fluorescence efficiency, and escape probability.
- **Output resolution/coverage:** Global; 1 km; source-day estimates, with 8-day averages used for tower comparisons; 2018-2020 study period.
- **Training/testing design and performance:** NA - the method is physically formulated and has no conventional training/test split. Its primary accuracy assessment is the independent OCO-2 comparison reported below.
- **External validation:** HSIFtotal and OCO-2 total SIF were both aggregated to 0.05 degree, using cells with at least six OCO-2 observations. The comparison obtained R-squared = 0.78 and RMSE = 1.33 mW m^-2 nm^-1. As proxy validation, HSIFtotal averaged within 1 km of 135 flux towers correlated with 8-day tower GPP at R-squared = 0.70, compared with R-squared = 0.64 for low-resolution TROPOMI SIF averaged within 10 km.
- **Thesis-use note:** HSIF provides a physically motivated alternative to purely data-driven downscaling and explicitly conserves energy at the parent scale. Its finding that SIF-GPP agreement improves when spatial support approaches the flux footprint strongly reinforces the thesis's focus on footprint-aware aggregation. It cannot reconstruct dates without an original TROPOMI observation, so it is downscaling rather than full temporal gap filling.

#### 19. CSIF

- **Reference/title:** *A Global Spatially Contiguous Solar-Induced Fluorescence (CSIF) Dataset Using Neural Networks*.
- **Product/acronym:** CSIF; specifically CSIFclear-inst for clear-sky instantaneous SIF and CSIFall-daily for all-sky daily-average SIF.
- **Enhancement type:** Global reconstruction, spatial gap filling, and temporal extension of sparse OCO-2 SIF using continuous MODIS reflectance.
- **Supervisory SIF source:** Cloud-free OCO-2 nadir SIF at 757 nm, aggregated to 0.05-degree cells containing more than five valid soundings. Data from 2015-2016 were used for training and 2014 plus 2017 for temporal validation.
- **Predictors:** The first four MODIS MCD43C4 BRDF-corrected NBAR bands: red, near-infrared, blue, and green reflectance. Tests with all seven bands and temperature/VPD improved R-squared by less than 0.01 and were therefore not retained.
- **Model family:** Feed-forward neural network with one hidden layer and five neurons.
- **Output resolution/coverage:** Global; 0.05 degree; 4-day. CSIFclear-inst covers 2000-2017, while CSIFall-daily covers 2000-2016.
- **Training/testing design and performance:** The temporally separated validation used 2014 and 2017 observations after training on 2015-2016. Training performance was R-squared = 0.796 and RMSE = 0.182 mW m^-2 nm^-1 sr^-1; validation performance was R-squared = 0.786 and RMSE = 0.177 mW m^-2 nm^-1 sr^-1.
- **External validation:** NA for direct independent SIF validation. Cross-product comparisons showed consistency with GOME-2 SIF and GOME-2-based reconstructions, but no single external headline score was stated. Proxy validation against GPP at 40 FLUXNET2015 sites gave site-level R-squared values from 0.01 to 0.93, with median R-squared = 0.64 for CSIFall-daily and mean GPP RMSE = 1.67 g C m^-2 day^-1.
- **Thesis-use note:** CSIF is a foundational reflectance-only OCO-2 reconstruction and a useful simple neural-network baseline. The temporally held-out years are stronger than a random split, but direct fine-resolution SIF validation is absent. The product primarily captures vegetation structure and absorbed-light variation, while environmental down-regulation can remain in the residual between observed OCO-2 SIF and CSIF.

#### 20. CNSIF

- **Reference/title:** *CNSIF: A Reconstructed Monthly 500-Meter Spatial Resolution Solar-Induced Chlorophyll Fluorescence Dataset in China*.
- **Product/acronym:** CNSIF.
- **Enhancement type:** National long-term SIF reconstruction and spatial downscaling of GOSIF to 500 m through deep learning and transfer learning.
- **Supervisory SIF source:** Monthly 0.05-degree GOSIF, an OCO-2-based reconstructed SIF product.
- **Predictors:** NIRv, EVI, near-infrared reflectance, and land-surface temperature derived from Landsat 7 ETM+ for January 2003-March 2017 and from Sentinel-2 MSI reflectance plus Landsat 8 TIRS temperature thereafter. GOSIF is included as the coarse SIF input/target within the CNN framework, and MODIS land cover supports classification and analysis.
- **Model family:** Two-dimensional CNN with four 3 x 3 convolutional layers and one fully connected layer. A 9 x 9 fine-pixel input window was selected, and transfer learning moved the relationship learned at coarse support to the 500 m output grid.
- **Output resolution/coverage:** China; 500 m; monthly; 2003-2022.
- **Training/testing design and performance:** Monthly data from 2003-2018 and 2020-2022 were used for model development, while 2019 GOSIF was excluded for temporally independent validation. The selected 9 x 9 model obtained R-squared = 0.856 and RMSE = 0.089 mW m^-2 sr^-1 nm^-1 against 2019 GOSIF. This is a same-product consistency/generalization test rather than absolute fine-scale truth.
- **External validation:** Monthly ChinaSpec tower SIF from nine sites, scaled to the CNSIF pixel using NIRv, produced pooled R-squared = 0.581 and RMSE = 0.152 mW m^-2 sr^-1 nm^-1; site-level R-squared ranged from 0.324 to 0.947. At the homogeneous DaMan maize site, R-squared was 0.760 for CNSIF versus 0.591 for GOSIF. Proxy validation against monthly tower GPP from ten flux sites produced R-squared = 0.550.
- **Thesis-use note:** CNSIF is particularly relevant because it combines a coarse reconstructed SIF product with Landsat/Sentinel-2 information and explicitly targets agricultural heterogeneity. However, errors can propagate from GOSIF, and the tower observations require NIRv-based spatial scaling. Its 2019 validation demonstrates temporal generalization to its supervisory product, while the tower comparison provides the more meaningful independent accuracy evidence.

#### 21. coSIF

- **Reference/title:** *Spatial Statistical Prediction of Solar-Induced Chlorophyll Fluorescence (SIF) from Multivariate OCO-2 Data*.
- **Product/acronym:** coSIF.
- **Enhancement type:** De-noising and spatial gap filling of OCO-2 SIF with formal pixel-level uncertainty quantification.
- **Supervisory SIF source:** Quality-controlled OCO-2 version 10r SIF at 740 nm, monthly averaged to 0.05 degree. OCO-2 XCO2 from the following month is the secondary observed process used for cokriging.
- **Predictors:** Spatial covariance among SIF observations; cross-covariance with one-month-lagged OCO-2 XCO2; spatial bisquare basis functions for large-scale trends; and MODIS land cover only as a land/water prediction mask.
- **Model family:** Multivariate geostatistical cokriging with a full bivariate Matern covariance model, non-stationary measurement errors, and coherent root-mean-squared prediction errors for every output cell. SIF-only kriging and a trend-surface-only model are baselines.
- **Output resolution/coverage:** North America from 22-58 degrees N and 125-65 degrees W; 0.05 degree; monthly. Demonstrations are provided for February, April, July, and October 2021 using XCO2 from the following months.
- **Training/testing design and performance:** Rather than random cross-validation, the study withheld two spatially contiguous 5 x 5 degree blocks in July 2021: an Iowa Corn Belt block and a Colorado/Kansas cropland block. For cokriging, Corn Belt BIAS = -0.06 and RASPE = 0.58; cropland BIAS = 0.01 and RASPE = 0.56. Kriging produced RASPE = 0.59 and 0.56, respectively, while the trend-only model produced 0.62 and 0.60. coSIF prediction uncertainties were on average about four times smaller than the gridded OCO-2 measurement-error standard deviations.
- **External validation:** NA - validation used spatially withheld observations from the same OCO-2 SIF dataset; no independent satellite, airborne, or tower SIF data were used.
- **Thesis-use note:** coSIF is methodologically distinct because it exploits the spatial covariance of SIF and returns uncertainty for each gap-filled value. Its contiguous-block validation is well aligned with the practical task of predicting genuinely unobserved areas. The method improves uncertainty characterization more clearly than point accuracy, and its 0.05-degree monthly implementation does not itself produce Sentinel-2-scale spatial detail.


## Section 2 - Crop Yield Prediction
