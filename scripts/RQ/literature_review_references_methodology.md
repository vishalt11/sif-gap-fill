# Literature Review - Methodology Notes

## Purpose and use

This file is a compact writing aid for Thesis Sections 3.1--3.3. It supplements, rather than repeats, `literature_review_references_table.md`. Each entry explains the methodological workflow in brief and records only a small number of additional points that may help synthesize the literature. Numerical performance results and detailed validation values remain in the reference table.

The Section 3.1 classification is based on the paper's primary contribution. Some studies span more than one category.

## 1. GOSIF

- **Primary Section 3.1 group:** Reconstruction and gap filling.
- **Brief methodology:** Quality-controlled OCO-2 soundings are aggregated to the modelling grid and paired with MODIS EVI and meteorological drivers. A Cubist model learns the relationship between SIF, vegetation state, incoming radiation, atmospheric demand, and temperature at observed locations. The fitted relationship is then applied to spatially complete predictor grids, filling orbital gaps and extending OCO-2-like SIF back to the start of the MODIS record.
- **Useful synthesis point:** This is reconstruction rather than observation-preserving downscaling: the predicted value is the final SIF estimate, so no constraint forces it to reproduce a contemporaneous OCO-2 parent-cell value. The authors found that land-cover class added little once the selected continuous predictors were present, but performance still varied among biomes.
- **Section 3.3 comparison:** Yes - representative tabular machine-learning reconstruction baseline.

## 2. Random-forest downscaling of GOSIF in Henan Province

- **Primary Section 3.1 group:** Spatial downscaling.
- **Brief methodology:** GOSIF, NDVI, and daytime land-surface temperature are first placed on the same coarse grid, where a random forest learns the relationship between GOSIF and the two fine-resolution covariates. The trained relationship is then evaluated with NDVI and temperature at 1 km to create the downscaled field. The enhanced SIF is subsequently combined into an anomaly-based indicator for agricultural drought monitoring.
- **Useful synthesis point:** The workflow assumes that a relationship calibrated at coarse spatial support remains valid at 1 km. Its reaggregation check and GPP/yield application results demonstrate consistency and usefulness, but not direct accuracy of individual 1 km SIF pixels. The supervisory signal is already reconstructed GOSIF, so upstream errors may be inherited.
- **Section 3.3 comparison:** No - useful regional application example, but weaker validation than the selected representative products.

## 3. Global downscaled GOME-2 SIF (Duveiller et al.)

- **Primary Section 3.1 group:** Spatial downscaling.
- **Brief methodology:** A semi-empirical light-use-efficiency formulation relates coarse GOME-2 SIF to fine-resolution vegetation and stress proxies. NIRv represents absorbed-light and canopy effects, while NDWI and afternoon land-surface temperature represent water and temperature limitations. These fine-grid drivers spatially redistribute the parent SIF signal, and independent OCO-2 observations are used to select the explanatory-variable configuration; a later pixel-wise correction harmonizes the product with TROPOMI.
- **Useful synthesis point:** Unlike unconstrained prediction, this method retains the coarse measured signal while allocating its spatial distribution. It is therefore an early example of aggregation-constrained downscaling and a useful conceptual precedent for conservation across spatial supports. The authors also present the framework as transferable to future finer-resolution sensors, including combinations involving Sentinel-2 and Landsat.
- **Section 3.3 comparison:** No - methodologically important and worth discussing in Section 3.1, but TroDSIF and HSIF provide more direct representatives of constrained downscaling in the focused comparison.

## 4. BCSIF

- **Primary Section 3.1 group:** Spatial downscaling.
- **Brief methodology:** A random forest is trained at the coarse TROPOMI support using optical, illumination, and environmental variables, then applied to the same variables at 0.005 degrees. The authors calculate the discrepancy between observed and RF-predicted SIF at coarse scale, spatially distribute that bias to the fine grid, and add it to the direct fine-resolution prediction. This correction is intended to retain physiological information present in TROPOMI but absent from a universal predictor-SIF relationship.
- **Useful synthesis point:** BCSIF separates fine-scale pattern generation from preservation of coarse-scale observation information. This is especially relevant when predictors reproduce canopy structure well but miss physiological variability. Its tower comparison suggests that reduced spatial mismatch can improve validation, although evidence is limited to a small number of SIF sites.
- **Section 3.3 comparison:** Yes - representative residual/bias-corrected downscaling method.

## 5. CNN downscaling of GOSIF for agricultural drought

- **Primary Section 3.1 group:** Spatial downscaling.
- **Brief methodology:** Monthly GOSIF and fine-resolution NDVI, EVI, daytime temperature, and nighttime temperature are aligned at the coarse scale to train a CNN. Convolutional filters use neighbouring pixels to learn spatially structured relationships rather than treating every grid cell independently. The trained network is then supplied with the fine-grid variables to estimate SIF at approximately 1 km, after which downscaled SIF and temperature are combined in the temperature fluorescence dryness index.
- **Useful synthesis point:** The paper motivates CNNs through their ability to represent adjacent-pixel context in fragmented agricultural landscapes. However, its evidence for the fine SIF field comes mainly from reaggregated agreement with GOSIF and correlation with MODIS GPP, so the claimed 1 km accuracy remains indirect.
- **Section 3.3 comparison:** No - useful for motivating spatial input representations, but SIFnet and CNSIF provide stronger CNN comparators.

## 6. DSIFRFK_EA0.05

- **Primary Section 3.1 group:** Spatial downscaling, with historical reconstruction.
- **Brief methodology:** Coarse GOME SIF is modelled as a large-scale trend plus a spatial residual. A random forest estimates the trend from AVHRR vegetation variables and ERA5 climate data, and coarse residuals are interpolated to the target grid using ordinary kriging before being added back to the fine-resolution trend. Alternative global/month-specific models, predictor subsets, random forest alone, and inverse-distance residual interpolation are compared before selecting the final RF-kriging product.
- **Useful synthesis point:** The method demonstrates how residual information can supplement a predictor-driven downscaling model and extend fine-resolution SIF into the pre-MODIS period. However, interpolating residuals from 2-degree support to 0.05 degrees does not remove the change-of-support problem; apparent detail can arise from model and interpolation assumptions rather than direct fine-scale SIF evidence.
- **Section 3.3 comparison:** No - valuable historical and geostatistical example, but coSIF and BCSIF cover the selected comparison dimensions more directly.

## 7. Fine-grained OCO-2 SIF forecasting

- **Primary Section 3.1 group:** Other and hybrid enhancement approaches - forecasting.
- **Brief methodology:** The workflow first uses a random forest and optical/temperature variables to construct a spatially complete SIF series because sparse OCO-2 observations cannot directly provide consistent lagged inputs. The gap-filled SIF from the two preceding periods, EVI, and spatial-temporal variables are then fused with OCO-2 samples and supplied to linear, tree-based, CNN, and LSTM forecasting models. Models trained on 2017-2018 predict monthly or seasonal SIF in 2019.
- **Useful synthesis point:** This paper addresses a different target from same-date reconstruction: future SIF. It also exposes a two-stage dependency, because errors and smoothing introduced during initial gap filling become predictors in the forecasting stage. Nearest-pixel matching is simpler than footprint-based aggregation and may introduce additional scale mismatch.
- **Section 3.3 comparison:** No - retain as the principal forecasting example in Section 3.1.

## 8. Spatially continuous TanSat SIF

- **Primary Section 3.1 group:** Reconstruction and gap filling.
- **Brief methodology:** The study first analyses field and satellite evidence to identify drivers of canopy SIF, emphasizing illumination as well as vegetation state. A random forest is then trained on TanSat SIF using MODIS visible-to-near-infrared reflectance, NDVI, noon solar-zenith-angle cosine, and temperature. Applying that model to globally complete predictor grids creates continuous TanSat-like SIF between the sparse observation tracks.
- **Useful synthesis point:** Its main methodological contribution is showing that solar illumination geometry is not merely a nuisance correction but an important predictor of SIF magnitude. Omitting the solar-zenith proxy measurably degraded reconstruction performance, which supports including observation geometry or radiation information when samples span changing illumination conditions.
- **Section 3.3 comparison:** No - useful illumination-focused reconstruction example, but GOSIF and ST-LGBM are more central to the focused comparison.

## 9. DOSIF

- **Primary Section 3.1 group:** Reconstruction and gap filling.
- **Brief methodology:** OCO-3 observations and daily MODIS reflectance are linked through Moving Spatial-Temporal Window Sampling. For each target day of year, training samples are drawn from a surrounding temporal window across available years and divided into geographic sub-biomes, allowing separate models to represent local phenology and vegetation behaviour. Several algorithms are benchmarked and CatBoost is selected to generate a global daily SIF record throughout the MODIS period.
- **Useful synthesis point:** DOSIF replaces one universal mapping with day-specific and sub-biome-specific mappings, trading computational complexity for better contextual transfer. This is useful evidence that SIF-reflectance relationships are conditional on phenology and region. Its independent airborne CFIS validation is valuable, but the detailed Results and Conclusion report inconsistent headline values; the detailed Results values in the reference table should be used when writing.
- **Section 3.3 comparison:** No - a strong reconstruction paper, but omitted to keep the comparative subset focused and manageable.

## 10. ANN-based TROPOMI SIF retrieval

- **Primary Section 3.1 group:** Other and hybrid enhancement approaches - retrieval-precision enhancement.
- **Brief methodology:** OCO-2/3 observations are matched to TROPOMI footprints under restrictions on overpass time, solar geometry, viewing geometry, and reference noise. Multiple OCO observations within a TROPOMI footprint are averaged to create a low-noise SIF label, while the corresponding high-spectral-resolution TROPOMI radiance between 743 and 758 nm forms the input. Separate feed-forward ANNs for detector columns learn to retrieve SIF directly from the Fraunhofer-line spectral structure.
- **Useful synthesis point:** This is not reflectance-based reconstruction: the network extracts fluorescence information from hyperspectral radiance itself and improves single-retrieval precision, reducing the need for spatial or temporal averaging. The distinction matters because retrieval enhancement improves the observation supplied to later gap-filling models rather than estimating SIF where no radiance observation exists. Remaining uncertainty includes cross-sensor geometry/time mismatch in the training labels.
- **Section 3.3 comparison:** No - retain as the principal retrieval-enhancement example in Section 3.1.

## 11. RTSIF

- **Primary Section 3.1 group:** Reconstruction and gap filling.
- **Brief methodology:** Clear-sky TROPOMI SIF is aggregated and used to supervise XGBoost with MODIS reflectance, land-surface temperature and land cover, CERES PAR, and C3/C4 vegetation fractions. The variable selection is motivated by a light-use-efficiency decomposition: reflectance represents absorbed radiation and canopy state, PAR represents illumination, and temperature and vegetation type help represent variation in fluorescence efficiency. The fitted model is applied to the longer ancillary-data archive to reconstruct TROPOMI-like SIF back to 2001.
- **Useful synthesis point:** RTSIF exemplifies temporal extension from a short but information-rich sensor record. The model can reproduce a TROPOMI-consistent historical series only under the assumption that the learned predictor-SIF relationship is temporally transferable; historical physiological responses not captured by the covariates cannot be recovered directly.
- **Section 3.3 comparison:** No - useful long-term reconstruction example, but GOSIF and ST-LGBM better match the focused comparison questions.

## 12. ROSIF and the photosynthetic afternoon-depression study

- **Primary Section 3.1 group:** Reconstruction and gap filling, with a downstream physiological application.
- **Brief methodology:** Separate random-forest and XGBoost mappings are developed between the eight GOCI reflectance bands and OCO-3 SIF near 10:00 and 14:00 local time; XGBoost is selected to reconstruct continuous 500 m morning and afternoon fields. SIF is then normalized by PAR and NIRv to approximate fluorescence yield. The difference between morning and afternoon yield forms an afternoon-depression indicator that is compared with independent drought indices.
- **Useful synthesis point:** The study shows why sub-daily timing can matter: drought information may appear in the change between morning and afternoon rather than in a single SIF snapshot. The physiological interpretation depends on simplifying assumptions in the NIRv-based yield decomposition, especially in sparse canopies, and the ground SIF record does not constitute a direct validation of the reconstructed study-period product.
- **Section 3.3 comparison:** No - retain as the sub-daily reconstruction and drought-response example in Section 3.1.

## 13. Hybrid kriging with external drift

- **Primary Section 3.1 group:** Other and hybrid enhancement approaches - machine learning plus geostatistics.
- **Brief methodology:** The existing CSIF machine-learning estimate provides a spatially complete external drift, while observed OCO-2 aggregates provide the local SIF values and empirical covariance structure. Within moving windows, kriging with external drift adjusts the broad ancillary-variable estimate according to neighbouring OCO-2 residual behaviour. Ordinary kriging and standalone CSIF are evaluated as the two parent methods.
- **Useful synthesis point:** The hybrid combines complementary information: machine learning supplies spatially complete environmental structure, whereas kriging retains local information from actual SIF observations. Its improvement over both parents supports hybrid extensions to the thesis, although leave-one-out validation with neighbouring observations available is less demanding than filling a broad region with no nearby OCO-2 swath.
- **Section 3.3 comparison:** No - methodologically relevant, but coSIF is retained as the principal geostatistical comparator and the focused set already includes several hybrids.

## 14. ST-LGBM OCO-2 SIF reconstruction

- **Primary Section 3.1 group:** Reconstruction and gap filling.
- **Brief methodology:** A LightGBM baseline relates OCO-2 SIF to NIRv, PAR, vapour-pressure deficit, temperature, and land cover. ST-LGBM adds two observation-derived covariates: a spatial factor formed from geographically nearby pixels with similar vegetation state and a temporal factor formed from phenologically corresponding observations in other years. These constraints inject information from the sparse SIF swaths themselves into predictions between swaths.
- **Useful synthesis point:** The paper directly demonstrates that random validation is optimistic for spatial gap filling: performance declines when complete regions rather than scattered samples are withheld. The added constraints reduce that decline and remain helpful across alternative model families and input combinations. This validation logic is especially relevant to evaluating transfer beyond observed OCO-2 footprints.
- **Section 3.3 comparison:** Yes - representative spatiotemporally constrained reconstruction and spatial-holdout study.

## 15. SIFnet

- **Primary Section 3.1 group:** Spatial downscaling.
- **Brief methodology:** SIFnet learns a tenfold scale transformation by first coarsening TROPOMI SIF and training a residual CNN to reproduce the original resolution from coarse SIF plus fine auxiliary rasters. Once optimized, the same learned scale factor is applied with native TROPOMI and finer covariates to estimate 0.005-degree SIF. Convolutional neighbourhoods encode spatial context, and the loss combines numerical error with structural dissimilarity so that both values and spatial patterns influence training.
- **Useful synthesis point:** The continental holdout and OCO-2/OCO-3 footprint comparison are stronger than random same-region validation. Feature interpretation shows that the parent SIF remains the dominant input and NIRv is the leading auxiliary feature. Nevertheless, the loss does not explicitly require the mean of the fine cells to equal the coarse TROPOMI observation, and performance remains heterogeneous in low-signal drylands and complex urban landscapes.
- **Section 3.3 comparison:** Yes - principal CNN downscaling comparator and the closest published analogue to footprint-averaged validation in the thesis.

## 16. TroDSIF

- **Primary Section 3.1 group:** Spatial downscaling, with observation-preserving redistribution.
- **Brief methodology:** A random forest first predicts fine-resolution SIF from reflectance, NDVI, illumination geometry, and temperature. These predictions are not used as the final product; instead, their relative spatial pattern supplies weights for redistributing the original TROPOMI SIF through Gaussian neighbourhoods. Special rules handle negative values and gaps, while reaggregation is used to verify preservation of the parent signal.
- **Useful synthesis point:** TroDSIF explicitly distinguishes measured coarse information from model-generated fine structure, which makes it an important aggregation-aware comparator. Its central limitation is the assumed linear proportionality between the original SIF and the machine-learning prediction used as a weight; violations of that relationship can distort the subpixel allocation even when the coarse signal is preserved.
- **Section 3.3 comparison:** Yes - representative aggregation-preserving redistribution method.

## 17. SIFoco2_005

- **Primary Section 3.1 group:** Reconstruction and gap filling.
- **Brief methodology:** Noise-reduced OCO-2 samples are paired with MODIS seven-band reflectance, and separate feed-forward ANN models are fitted for each biome and 16-day time step. Geographic sub-biomes balance the training distribution, while automated cross-validation selects the ANN size. The fitted local-in-ecology and local-in-time mappings are applied to spatially continuous MODIS fields in OCO-2 orbital gaps.
- **Useful synthesis point:** The paper's “physiological constraint” is implemented mainly through biome and time stratification rather than an explicit process equation. Removing those strata caused seasonal under- or overestimation and weakened drought sensitivity, showing that a universal reflectance-SIF relationship can smooth physiologically important variation. Validation with airborne CFIS specifically away from OCO-2 tracks is a notable strength.
- **Section 3.3 comparison:** No - important reconstruction evidence, but CSIF and ST-LGBM cover the selected ANN and constrained-reconstruction contrasts.

## 18. HSIF

- **Primary Section 3.1 group:** Spatial downscaling, with physical and aggregation constraints.
- **Brief methodology:** The method begins from the requirement that the mean of the 1 km subpixels reproduce the low-resolution TROPOMI observation. A light-use-efficiency and radiative-transfer formulation separates total canopy fluorescence emission from directional observed SIF using chlorophyll FPAR, fluorescence efficiency, and escape probability. Fine-resolution optical and canopy variables then allocate total and observed SIF within each coarse support while maintaining energy conservation.
- **Useful synthesis point:** HSIF is an important alternative to empirical scale transfer because its fine values are governed by an explicit physical decomposition and coarse conservation rule. It also demonstrates that matching the SIF averaging radius to a flux-tower footprint strengthens the SIF-GPP relationship. The product cannot fill dates without a parent TROPOMI observation and remains sensitive to uncertainty in BRDF, FPAR, and escape-probability inputs.
- **Section 3.3 comparison:** Optional - include if space permits to represent physically constrained downscaling; otherwise discuss it in Section 3.1 beside TroDSIF.

## 19. CSIF

- **Primary Section 3.1 group:** Reconstruction and gap filling.
- **Brief methodology:** OCO-2 SIF is aggregated and paired with the first four MODIS BRDF-corrected reflectance bands. A small feed-forward neural network learns the clear-sky instantaneous reflectance-SIF relationship, which is then applied to gap-filled MODIS reflectance over the full record. Radiation scaling converts the reconstructed clear-sky instantaneous signal into clear-sky and all-sky daily products.
- **Useful synthesis point:** CSIF deliberately uses a parsimonious reflectance-only mapping because additional bands and meteorology brought little average validation improvement. The resulting product can be interpreted as the SIF expected from canopy optical state and illumination; deviations of observed OCO-2 SIF from CSIF may therefore retain information on short-term environmental down-regulation. That interpretation is useful, but it also highlights what the reconstruction itself may fail to reproduce under stress.
- **Section 3.3 comparison:** Yes - representative neural-network reconstruction and reflectance-only baseline.

## 20. CNSIF

- **Primary Section 3.1 group:** Spatial downscaling combined with long-term reconstruction.
- **Brief methodology:** A 2D CNN receives GOSIF together with fine-resolution NIRv, EVI, near-infrared reflectance, and land-surface temperature. Neighbourhood windows represent the internal spatial structure hidden inside each coarse SIF cell; the selected network is pretrained at coarse support and transferred to 500 m using Landsat/Sentinel-2 reflectance and Landsat thermal information. Separate yearly reconstruction models generate the long China-wide record, with one year withheld for temporal consistency testing.
- **Useful synthesis point:** CNSIF is highly relevant to heterogeneous agriculture because it explicitly uses subpixel image structure and shows improved delineation of winter-wheat fields and fragmented landscapes. At the same time, GOSIF is already a modelled product, so CNSIF can inherit upstream assumptions, and the 2019 GOSIF comparison tests temporal consistency rather than absolute fine-pixel truth. The independent tower analysis is therefore the more important validation evidence.
- **Section 3.3 comparison:** Yes - representative 500 m deep-learning method using Landsat/Sentinel-2 information and an agricultural-landscape use case.

## 21. coSIF

- **Primary Section 3.1 group:** Other and hybrid enhancement approaches - multivariate geostatistical gap filling.
- **Brief methodology:** Monthly OCO-2 SIF and the following month's OCO-2 XCO2 are represented as noisy observations of two latent spatial processes. Large-scale trends are fitted with spatial basis functions, while a full bivariate Matern model describes residual covariance within each variable and cross-covariance between them. Cokriging predicts SIF on the complete grid and returns a coherent prediction-error estimate for every cell; SIF-only kriging and trend-only prediction form nested baselines.
- **Useful synthesis point:** coSIF shifts emphasis from generating visually detailed predictions to quantifying what is known at unobserved locations. Spatial-block validation and explicit measurement errors make it especially rigorous for real gaps. XCO2 is helpful only when meaningful lagged cross-dependence exists; during weak-dependence periods, simpler SIF-only kriging can perform comparably and more efficiently.
- **Section 3.3 comparison:** Yes - principal uncertainty-aware geostatistical comparator.


# Section 3.2 - Crop Yield and GPP Estimation

The entries below remain in the numbered folder order. They are not assigned to literature groups yet. These notes emphasize how each study links its input variables to yield, productivity, or GPP and record only the most apparent limitations.

## 1. Saxony regional yield uncertainty

- **Brief methodology:** Random-forest models are fitted separately for each crop, forecast date, and geographical setup. The paper contrasts one Saxony-wide model with models trained within agriculturally comparable regions (ACRs), and predicts a year kept outside the training period.
- **Variables used:** Meteorological variables accumulated to different dates and historical regional crop-yield records.
- **Geographical scope:** Saxony, Germany; several crops, including winter wheat, over 2015--2022.
- **Important issue:** Agricultural landscapes and crop responses are too heterogeneous for one model to work equally well everywhere. Strong training fits sometimes collapsed in the withheld year, whereas local ACR models were often more stable.

## 2. SIF benefit for crop-yield prediction

- **Brief methodology:** Comparable county-yield models are fitted with coarse GOME-2 SIF, downscaled SIF, SIF-derived productivity variables, NDVI, or MODIS GPP. Their out-of-sample in-season forecasts are compared at progressively later dates in the growing season.
- **Variables used:** GOME-2 SIF, statistically downscaled SIF, SIF productivity metrics, MODIS NDVI, MODIS GPP, and county effects.
- **Geographical scope:** Major agricultural regions of the United States at county scale.
- **Important issue:** A physiological variable is not automatically a better operational predictor. Downscaling introduced uncertainty, and NDVI was at least as useful as SIF in several tests and more useful under drought.

## 3. Spaceborne SIF crop productivity

- **Brief methodology:** The study converts GOME-2 SIF into crop carbon uptake and then into harvested production. Instead of using one empirical SIF--yield line, it accounts for differences in C3/C4 photosynthesis, crop carbon composition, and respiratory carbon losses before comparing results with county statistics.
- **Variables used:** GOME-2 SIF, crop type/photosynthetic pathway, parameters for crop stoichiometry and respiration, temperature information, and reported crop yields.
- **Geographical scope:** U.S. croplands and county statistics during 2007--2012.
- **Important issue:** The framework has a clearer physiological interpretation than a direct regression, but yield depends on additional processes such as respiration and biomass allocation. These conversion assumptions can introduce uncertainty even when SIF represents photosynthesis well.

## 4. Multi-year field wheat yield

- **Brief methodology:** Sentinel-2 observations and climate variables are summarized for individual wheat fields and used in stepwise linear-regression and random-forest models. The paper tests raw bands, spectral indices, and water-balance variables across several growing seasons.
- **Variables used:** Sentinel-2 bands, NDRE, REIP, NDWI, precipitation, reference and crop evapotranspiration, and crop water requirement.
- **Geographical scope:** Three agricultural regions in southern Germany during 2016--2018, using directly measured field yields.
- **Important issue:** Raw bands plus water/weather information performed better than relying on vegetation indices alone. The limited number of regions and seasons may make the fitted relationships sensitive to local climate and management.

## 5. NCP SIF wheat-yield machine learning

- **Brief methodology:** Municipal yield records are paired with growing-season SIF, vegetation, climate, and soil variables. Five algorithms are evaluated using repeated train/test splits, after which the models predict the independent year 2020.
- **Variables used:** SIF; NIRv, NDVI, EVI, and EVI2; precipitation, soil moisture, minimum temperature, vapour-pressure deficit, and radiation; and soil texture, pH, cation-exchange capacity, and organic matter.
- **Geographical scope:** Fifty-eight municipalities in five provinces of the North China Plain, 2007--2020; pixels with more than 60% winter-wheat cover were retained.
- **Important issue:** XGBoost achieved the strongest individual test result, whereas random forest was more stable. SIF was usually more useful than NIRv, and adding more predictor groups did not always improve performance because information overlapped.

## 6. Lightweight Australian wheat model

- **Brief methodology:** A compact process-guided model uses SIF and environmental stress variables to estimate the fraction of open PSII reaction centres ($q_L$). It uses this regulation term to estimate GPP and then converts accumulated productivity into wheat yield and regional production.
- **Variables used:** Satellite SIF, air and dew-point temperature, soil-water information, vapour-pressure deficit, harvested area, and crop conversion parameters.
- **Geographical scope:** Australian wheat regions and states during 2019--2022.
- **Important issue:** The method needs fewer inputs than a full crop-simulation model, but coarse SIF and uncertain harvested-area statistics affect production estimates. The fitted $q_L$ response may also require recalibration outside Australian conditions.

## 7. Scalable MLR-SIF yield framework

- **Brief methodology:** A mechanistic light-reaction (MLR) model corrects observed SIF for canopy escape effects and uses absorbed radiation to estimate GPP. GPP is converted to NPP and then to crop yield using an NPP:GPP ratio, harvest index, and harvested area. Calibrated and uncalibrated versions are compared with ANN and random forest.
- **Variables used:** OCO-2 SIF, PAR, MODIS NDVI/NIR/fPAR, an estimated fluorescence escape probability, MODIS NPP:GPP ratios, harvested area, and harvest index.
- **Geographical scope:** Corn in 210 U.S. Corn Belt counties and wheat in 55 districts of the Indian Indo-Gangetic Plain.
- **Important issue:** Machine learning was competitive where abundant, crop-pure observations were available, while the process-guided model transferred better to data-limited India. Coarse mixed OCO-2 footprints and cloud cover reduced the reliability of crop-specific SIF there.

## 8. Ground-to-space crop productivity

- **Brief methodology:** The study first links field-measured SIF with eddy-covariance GPP and uses SCOPE simulations to interpret canopy and C3/C4 effects. It then aggregates TROPOMI SIF over cropland and compares seasonal SIF with tower productivity and county crop statistics.
- **Variables used:** Field PhotoSpec SIF, tower GPP, SCOPE simulations, TROPOMI SIF, planted area, C3/C4 crop fractions, county production/yield, and NPP estimates.
- **Geographical scope:** Corn and soybean systems from field sites to Iowa/U.S. Corn Belt county scales.
- **Important issue:** SIF--productivity ratios differ between C3 and C4 crops, so mixed pixels require crop-composition correction. Cloud sampling, high-light behaviour, and the fraction of fluorescence escaping the canopy also affect scale comparisons.

## 9. Belgium NDVI wheat yield

- **Brief methodology:** Seasonal NDVI is summarized as its maximum and time integral for each farm field and supplied to random forest together with monthly weather. The study checks whether the satellite variables improve farm-yield prediction across years.
- **Variables used:** Maximum NDVI, integrated NDVI (aNDVI), monthly precipitation, and recorded winter-wheat yield.
- **Geographical scope:** Northern Belgium, with 666 fields in 2016, 609 in 2017, and 210 in 2018.
- **Important issue:** NDVI summaries were less informative than rainfall near tillering and anthesis. Greenness can remain high while water or heat stress reduces photosynthesis and final grain formation.

## 10. Practical SIF-$q_L$ crop model

- **Brief methodology:** Leaf measurements are used to build a simple model of $q_L$, the fraction of open PSII reaction centres. Satellite SIF is combined with $q_L$ and phenology to estimate daily crop GPP, which is accumulated and converted to county yield.
- **Variables used:** Top-of-canopy satellite SIF, air temperature, vapour-pressure deficit, soil moisture, leaf $q_L$ measurements, crop phenology, harvested area, and conversion parameters.
- **Geographical scope:** Corn and soybean across more than 700 counties in 12 U.S. Midwestern states during 2018--2023, with two AmeriFlux evaluation sites.
- **Important issue:** Adding $q_L$ helps separate changes in photochemical regulation from changes in SIF magnitude. The parameters are crop-specific, and coarse satellite observations still mix crops and management conditions.

## 11. Regional GPP from downscaled SIF

- **Brief methodology:** Moving-window regressions relate coarse GOME-2 SIF to MODIS variables, and the window with the best local fit is used to create finer SIF. A statistical SIF--GPP relationship then estimates regional GPP, which is compared with eddy-covariance observations.
- **Variables used:** GOME-2 SIF, MODIS NDVI, fPAR, soil-moisture index, land-surface temperature, and flux-tower GPP.
- **Geographical scope:** Regional demonstrations evaluated at eddy-covariance sites; the exact overall domain is less central than the comparison of local and fixed regression windows.
- **Important issue:** Local windows can represent spatially varying SIF--predictor relationships better than one global model. However, fine-scale patterns partly come from MODIS proxies rather than direct fine-resolution SIF measurements.

## 12. German ensemble crop yield

- **Brief methodology:** Six regression algorithms generate parcel-level predictions from Earth-observation and environmental features. Stacking and majority-voting ensembles combine the base predictions, which are then aggregated to districts and a regular grid for comparison with agricultural statistics.
- **Variables used:** Multi-source satellite observations/vegetation features, meteorological time series, soil variables, parcel boundaries, and official crop-yield statistics.
- **Geographical scope:** Two German federal states, approximately 140,000--155,000 parcels per year during 2019--2022; winter wheat, winter barley, and winter rapeseed.
- **Important issue:** The ensemble gives scalable and stable estimates across several crops, but requires a large parcel/yield training system. District aggregation improves reported accuracy while masking some parcel-level errors.

## 13. SIF-EC transfer-learning GPP

- **Brief methodology:** A neural model is first pretrained on spatially extensive GPP estimates derived from SIF. Its learned representation is then fine-tuned with eddy-covariance GPP so that tower observations correct the SIF-based relationship without discarding its broad spatial information.
- **Variables used:** Three alternative SIF-based GPP datasets, eddy-covariance GPP, and the environmental/remote-sensing inputs used by the common GPP model; the complete predictor list was not extracted for these brief notes.
- **Geographical scope:** Long-term global terrestrial GPP with evaluation across the available flux-tower network.
- **Important issue:** The design addresses the central trade-off between wide satellite coverage and sparse but more direct ground GPP. Bias remains possible where tower sites do not represent the climate, vegetation, or management conditions of unsampled regions.

## 14. Sentinel-2 wheat-yield mapping

- **Brief methodology:** Combine-harvester yield points are spatially matched with Sentinel-2 and environmental layers. A random forest learns this point-to-predictor relationship and applies it across wheat fields to produce 10 m yield maps.
- **Variables used:** Sentinel-2 bands/vegetation information, meteorological variables, topography, soil moisture, field boundaries, and combine-yield measurements.
- **Geographical scope:** Thirty-nine wheat fields in the United Kingdom, with more than 8,000 yield observations.
- **Important issue:** The approach produces useful within-field detail, but that detail is supervised by local combine data. The paper's limited seasonal and geographic sample leaves model transfer without new yield measurements uncertain.

## 15. Reconstructed SIF versus GPP

- **Brief methodology:** Four reconstructed SIF datasets are sampled at flux-tower locations and compared with GPP at 4-day, 8-day, and monthly scales. Their performance is also compared with NDVI, EVI, and NIRv, including a separate analysis of the 2018 European drought.
- **Variables used:** CSIF, GOSIF, LUE-SIF, HSIF, FLUXNET2015/ICOS GPP, NDVI, EVI, NIRv, land-cover class, and drought conditions.
- **Geographical scope:** Global flux-tower sites from 2007--2018 plus European sites affected by the 2018 drought.
- **Important issue:** This paper evaluates products rather than predicting crop yield. Reconstructed SIF generally tracks GPP better than greenness indices, but relationships weaken at a substantial subset of drought-affected sites and differ with time scale and product.

## 16. GOSIF drought crop productivity

- **Brief methodology:** Seasonal anomalies in GOSIF, coarse GOME-2 SIF, and MODIS vegetation indices are compared with tower GPP and crop-yield anomalies. Regression tests determine which satellite variable best represents normal interannual changes and the severe 2012 drought.
- **Variables used:** GOSIF, GOME-2 SIF, NDVI, EVI, NIRv, eddy-covariance GPP, corn/soybean yield statistics, and drought anomalies.
- **Geographical scope:** U.S. Midwest croplands during 2008--2018.
- **Important issue:** GOSIF responded much more strongly to the drought than greenness indices and closely followed observed yield loss. Because GOSIF is itself reconstructed using environmental inputs, its advantage cannot be attributed to the original OCO-2 observations alone.

## 17. Sentinel-2 crop-model wheat yield

- **Brief methodology:** Time-series Sentinel-2 indices describe canopy amount and chlorophyll status, while a crop-simulation model supplies a water-stress index. Statistical yield models test the satellite variables alone and in combination with the simulated stress variable, followed by cross-validation and an independent field test.
- **Variables used:** OSAVI, NDVI, chlorophyll index (CI), NDRE, a crop-model stress index (SI), and measured field yield.
- **Geographical scope:** One hundred and three dryland wheat fields in northeastern Australia during 2016--2017.
- **Important issue:** Structural, chlorophyll, and water-stress information are complementary, giving a large improvement when combined. Prediction remains difficult for severely stressed, very low-yield fields, where small absolute errors become large percentages.

## 18. Northeast Germany sensor suitability

- **Brief methodology:** Dense combine-harvester yield maps are matched with images from six optical satellite sensors. The authors calculate 15 indices and examine field-wise correlations with yield as a function of sensor resolution, red-edge availability, acquisition date, soil, terrain, and field heterogeneity.
- **Variables used:** Multispectral reflectance and 15 spectral indices from Sentinel-2, RapidEye, Landsat and other sensors; yield-monitor data; soil and relief variables; acquisition timing and field characteristics.
- **Geographical scope:** Northeast Germany, 13 years, 947 field-yield datasets, and 755 satellite images for cereals and canola.
- **Important issue:** Image timing relative to phenology and spatial resolution were more important than any universally best index. The study reports correlations rather than independently validated yield forecasts, so isolated high coefficients are not directly comparable with model test scores.

## 19. CONUS SIF winter-wheat yield

- **Brief methodology:** County yield is modelled with XGBoost and random forest using separate and combined SIF, vegetation-index, climate, and soil feature groups. Comparing these configurations measures the additional contribution of SIF beyond more conventional predictors.
- **Variables used:** SIF, NIRv and two other vegetation indices, 13 climate variables, four soil variables, and county winter-wheat yield.
- **Geographical scope:** Conterminous United States at county scale over 14 years.
- **Important issue:** SIF improves the vegetation-index model, but the improvement is moderate because predictor groups contain overlapping seasonal information. National models must also absorb large differences in climate, soils, management, and wheat phenology.

## 20. German CNN wheat yield

- **Brief methodology:** Weekly environmental sequences are aligned with crop phenology and supplied to a one-dimensional CNN that learns short- and long-term patterns across the growing season. Entire recent years are withheld and compared with eight statistical and machine-learning baselines.
- **Variables used:** Weekly minimum/maximum temperature, radiation, precipitation, relative humidity, wind, soil properties, and sowing, flowering, and harvest dates.
- **Geographical scope:** Two hundred and seventy-one German counties, 1999--2019, with 2017--2019 used as separate temporal tests.
- **Important issue:** The CNN gains accuracy by learning temporal patterns without hand-designed seasonal summaries. It uses no satellite observations and omits important cultivar and management variables, which may contribute to the higher errors in eastern and northeastern Germany.

## 21. Spectral-index wheat yield

- **Brief methodology:** Canopy hyperspectral measurements are collected repeatedly under eight nitrogen treatments. Four indices are evaluated across phenological stages; the indices with useful variation are regressed against measured grain and biomass yield to identify the best observation stage.
- **Variables used:** NDVI, NDWI, SAVI, NDNI, nitrogen treatment, phenological stage, and measured grain/biomass yield.
- **Geographical scope:** An irrigated wheat experiment in central India over three winter seasons from 2014 to 2017.
- **Important issue:** Heading-stage NDVI and NDWI have very strong within-experiment relationships with yield. The model is based on canopy spectroscopy under controlled nitrogen treatments and has no independent regional validation, so it should not be treated as satellite-scale operational performance.

## 22. SIF wheat heat-stress response

- **Brief methodology:** Satellite SIF and greenness indices are tracked before, during, and after the 2010 heatwave. SIF is normalized using absorbed-light/illumination information to distinguish fluorescence-efficiency changes from canopy amount, and satellite anomalies are compared with temperature and reported wheat-yield loss.
- **Variables used:** SIF, fluorescence/SIF yield, NDVI, EVI, fPAR or absorbed-light information, air temperature, crop calendar, and regional wheat-yield statistics.
- **Geographical scope:** Winter wheat in the northwestern Indian Indo-Gangetic Plains during the 2010 heat-stress event.
- **Important issue:** SIF detects reduced photosynthetic function earlier than NDVI/EVI, which respond mainly after visible canopy changes. However, the SIF-based loss is larger than the reported yield loss, and coarse mixed pixels plus a single extreme event limit generalization.
