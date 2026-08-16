# Chapter 2 Background: Proposed Content

This file is a sentence-level plan for approximately four paragraphs. Each numbered point can become one sentence, although closely related points may be combined while writing the LaTeX text.

## Proposed subsection structure

```latex
\section{Background\label{sec:background}}

\subsection{Satellite SIF Retrieval and Photosynthetic Interpretation}

\subsection{Spectral Predictors and Spatial Reference Framework}
```

The first subsection should contain Paragraphs 1 and 2 below. The second subsection should contain Paragraphs 3 and 4.

## Subsection 1: Satellite SIF Retrieval and Photosynthetic Interpretation

### Paragraph 1 - SIF signal, satellite retrieval, and quality control

1. Solar-induced chlorophyll fluorescence (SIF) is the weak radiation emitted by chlorophyll, mainly between approximately 650 and 800 nm, after vegetation absorbs photosynthetically active radiation (PAR); it is produced alongside photochemistry and non-photochemical energy dissipation \cite{doughty2022global,porcar2014linking}.
2. Satellite instruments do not measure photosynthesis directly; they retrieve SIF because the emitted fluorescence partly fills narrow solar Fraunhofer absorption lines in the measured radiance spectrum \cite{doughty2022global}.
3. The OCO-2/OCO-3 retrievals use fitting windows near 757 and 771 nm, and the retrieved fractional fluorescence contribution is converted to absolute SIF in units of \(\mathrm{W\,m^{-2}\,sr^{-1}\,\mu\mathrm{m}^{-1}}\) \cite{doughty2022global}.
4. Retrieval quality depends on spectral-fit quality, continuum radiance, gas-band ratios, solar zenith angle, cloud contamination, and the land fraction inside the footprint; quality flag 0 denotes the best observations and flag 1 denotes good observations recommended for most analyses \cite{doughty2022global}.
5. Small negative retrieved values are physically possible as measurement estimates even though actual fluorescence emission cannot be negative, because single-sounding retrieval noise can be large relative to a weak true signal and its absolute uncertainty increases with radiance \cite{doughty2022global,kohler2018global}.
6. Negative values should therefore not be removed automatically: Doughty et al. suggest accepting a value when \(SIF+2\sigma\geq0\), treating it as questionable when only \(SIF+3\sigma\geq0\), and rejecting it when \(SIF+3\sigma<0\) \cite{doughty2022global}.
7. The thesis follows this logic for its combined 757/771 nm target, propagates the daily-corrected uncertainty, retains only the uncertainty-based accepted class in the principal datasets, and applies a broad range filter of \(-0.5\leq SIF<2.0\).
8. Source quality flags 0 and 1 and the nadir/glint measurement mode are retained as separate variables so that their behaviour can be compared rather than concealed within the uncertainty filter; cite the actual products as \cite{oco2_sif_v11r,oco2_sif_v112r}.
9. The instantaneous retrieval is also multiplied by the supplied daily correction factor, which scales the overpass-time observation according to the ratio between daily-mean clear-sky illumination and illumination at the satellite overpass \cite{doughty2022global}.

**Useful detail if more explanation is needed:** the thesis target is \(S=(S_{757}+1.5S_{771})/2\), with \(\sigma_S=\tfrac{1}{2}\sqrt{\sigma_{757}^{2}+(1.5\sigma_{771})^{2}}\). This definition belongs mainly in the methods chapter, so the Background only needs it if the quality-control sentence otherwise becomes unclear.

### Paragraph 2 - Illumination geometry and the physical link between SIF and GPP

1. Gross primary production (GPP) can be expressed with the light-use-efficiency relation \(GPP=LUE\times FAPAR\times PAR=LUE\times APAR\) \cite{monteith1972solar,monteith1977climate,liu2020estimating}.
2. PAR is incoming radiation in the photosynthetically active wavelength range of approximately 400--700 nm, FAPAR is the fraction of that radiation absorbed by vegetation, and APAR is the absorbed radiation calculated as \(FAPAR\times PAR\).
3. Light-use efficiency (LUE) describes how efficiently the absorbed radiation is converted into fixed carbon, so it connects the available light energy to GPP.
4. Canopy SIF can be written conceptually as \(SIF=f_{esc}\times APAR\times\Phi_F\), where \(\Phi_F\) is fluorescence yield at the photosystem level and \(f_{esc}\) is the fraction of emitted fluorescence that escapes the canopy and reaches the sensor \cite{berry2013new,porcar2014linking,liu2020estimating}.
5. SIF and GPP are often correlated because both contain APAR, but they are not identical: GPP additionally depends on LUE, whereas measured SIF depends on fluorescence yield, stress-related energy partitioning, canopy structure, and escape probability \cite{doughty2022global,liu2020estimating}.
6. Illumination geometry must also be considered because the cosine of the solar zenith angle (SZA) is a proxy for PAR under cloud-free conditions, although it is not SIF itself \cite{chen2020upscaling}.
7. Phase angle describes the relative alignment of the Sun, the observed surface, and the sensor; near-zero angles can produce a hotspot-like increase in observed SIF, while the reported effect is small once the absolute phase angle is above roughly \(20^{\circ}\) \cite{doughty2022global,kohler2018global}.
8. In `sif_data_builder.R`, relative azimuth was calculated as the wrapped difference between solar and viewing azimuth, and SZA, viewing zenith angle, and relative azimuth were then combined with the spherical cosine relation to obtain phase angle; its sign only records which side of the relative-azimuth direction the observation occupies.
9. The thesis diagnostic across nine MGRS tiles did not contain the low-phase-angle observations associated with the strongest hotspot effect and found no systematic increase of SIF with phase angle; `RQ/other_figures/sif_phase_angle_and_value_satellite_2022-05-07.pdf` can be cited as an example visual comparison.

## Subsection 2: Spectral Predictors and Spatial Reference Framework

### Paragraph 3 - Sentinel-2 reflectance bands and vegetation indices

1. Surface reflectance is the fraction of incoming radiation reflected by the land surface in a defined wavelength interval, and differences among visible, near-infrared, red-edge, and short-wave-infrared bands describe different canopy properties.
2. The final five spectral-index channels use six Sentinel-2 bands: B2 blue (about 490 nm, 10 m), B4 red (about 665 nm, 10 m), B5 red edge (about 705 nm, 20 m), B8 broad NIR (about 833 nm, 10 m), B8A narrow NIR (about 865 nm, 20 m), and B11 SWIR (about 1610 nm, 20 m).
3. B4 lies in a strong chlorophyll-absorption region, B5 samples the rapid red-edge transition associated with chlorophyll condition, B8/B8A respond strongly to leaf and canopy scattering, B11 is sensitive to canopy water absorption, and B2 supplies the blue correction term used by EVI.
4. NDVI, \((B8-B4)/(B8+B4)\), represents the amount and condition of photosynthetically active green vegetation, although it can become less sensitive in dense canopies \cite{tucker1979red}.
5. NDRE, \((B8A-B5)/(B8A+B5)\), combines the narrow NIR and red-edge bands and is useful for canopy chlorophyll or nitrogen status, particularly after a crop canopy becomes dense \cite{barnes2000coincident}.
6. NIRv, \(B8\times NDVI\), estimates the NIR reflectance attributable to vegetation rather than soil or other background components and is therefore related to intercepted radiation and GPP \cite{badgley2017canopy}.
7. NDMI, \((B8-B11)/(B8+B11)\), uses the contrast between NIR and SWIR reflectance to represent vegetation moisture; it follows the water-sensitive logic of Gao's NDWI, although this thesis uses Sentinel-2 B11 near 1.61 \(\mu\mathrm{m}\) rather than Gao's original 1.24 \(\mu\mathrm{m}\) band \cite{gao1996ndwi}.
8. EVI, \(2.5(B8-B4)/(B8+6B4-7.5B2+1)\), represents vegetation greenness while the blue and coefficient terms reduce atmospheric and canopy-background effects and improve sensitivity in high-biomass conditions \cite{huete2002overview}.
9. For the model inputs, the 10 m bands were area-averaged to the common 20 m grid before index calculation, while B5, B8A, and B11 remained at their native 20 m spacing; the full WASP product also supplied B3, B6, B7, and B12, but these bands were not used directly by the final five spectral indices \cite{hagolle_wasp_2018}.

**Band-to-index check:** B2 is used by EVI; B4 by NDVI, NIRv, and EVI; B5 by NDRE; B8 by NDVI, NDMI, NIRv, and EVI; B8A by NDRE; and B11 by NDMI.

### Paragraph 4 - CRS, MGRS tiles, and the meaning of spatial resolution

1. A coordinate reference system (CRS) defines how numerical coordinates correspond to positions on Earth; a geographic CRS records angular latitude and longitude, whereas a projected CRS gives planar coordinates that allow distances, areas, and pixel sizes to be expressed in metres.
2. OCO-2 footprint corners are initially represented in WGS 84 geographic coordinates (EPSG:4326), but footprint areas and centroids are calculated after transformation to the European equal-area CRS ETRS89-LAEA (EPSG:3035).
3. Sentinel-2 rasters are processed in the UTM CRS associated with each tile so that the common output grid has a consistent 20 m pixel spacing and aligns with the native optical imagery.
4. The Military Grid Reference System (MGRS) is not a separate image CRS; it is a grid-reference system built mainly on UTM zones and 100 km grid squares, and Sentinel-2 uses these tile identifiers to divide and distribute its imagery.
5. A resolution written in degrees describes an angular geographic grid rather than a fixed ground distance: \(0.5^{\circ}\) and \(0.05^{\circ}\) cells become narrower east--west toward higher latitudes.
6. Near \(50^{\circ}\) N in Germany, \(0.05^{\circ}\) is approximately 5.6 km north--south by 3.6 km east--west, while \(0.5^{\circ}\) is approximately 55.6 km by 35.7 km; these are approximate dimensions rather than universal cell sizes.
7. Spatial resolution must also be distinguished from spatial support: a raster may be resampled to 20 m pixels, but its independent information remains limited by the native resolution and footprint of the source observation.
8. In this thesis, the 20 m grid provides the output locations and Sentinel-2 spatial pattern, whereas observed OCO-2 SIF constrains footprint or 4 km aggregate means; the resulting maps should therefore be described as 20 m downscaled estimates rather than directly observed 20 m SIF.

## Citation and writing notes

- The two main sources for Subsection 1 are `doughty2022global` and `liu2020estimating`.
- The citation keys requested for the spectral indices are present in `manual.bib`: `badgley2017canopy`, `barnes2000coincident`, `tucker1979red`, `gao1996ndwi`, and `huete2002overview`.
- The illumination statement uses `chen2020upscaling`; the SIF/GPP equations use `monteith1972solar`, `monteith1977climate`, `berry2013new`, and `porcar2014linking`.
- The Sentinel-2 band wavelengths and native resolutions were checked against the official Sentinel-2 documentation. The CRS/MGRS statements were checked against official UTM/MGRS descriptions. These basic technical statements currently have no dedicated CRS/MGRS entry in `manual.bib`.
- Avoid stating that cosine SZA is SIF. It is only an illumination proxy under clear-sky conditions.
- Avoid describing negative SIF as necessarily invalid. The uncertainty test determines whether a negative retrieval is statistically compatible with a true value near zero.
- Avoid implying that the daily correction removes all weather and geometry effects. It normalizes the clear-sky diurnal illumination component but does not remove every atmospheric, directional, or physiological influence.

## External papers cited within the two main sources (titles only)

- New global observations of the terrestrial carbon cycle from GOSAT: Patterns of plant fluorescence with gross primary productivity
- Global retrievals of solar-induced chlorophyll fluorescence with TROPOMI: First results and intersensor comparison to OCO-2
- Overview of Solar-Induced Chlorophyll Fluorescence (SIF) from the Orbiting Carbon Observatory-2: Retrieval, cross-mission comparison, and global monitoring for GPP
- TROPOMI reveals dry-season increase of solar-induced chlorophyll fluorescence in the Amazon forest
- Disentangling changes in the spectral shape of chlorophyll fluorescence: Implications for remote sensing of photosynthesis
- Sun-induced chlorophyll fluorescence is more strongly related to absorbed light than to photosynthesis at half-hourly resolution in a rice paddy
- On the covariation of chlorophyll fluorescence and photosynthesis across scales
- Chlorophyll fluorescence---a practical guide
- Canopy structure explains the relationship between photosynthesis and sun-induced chlorophyll fluorescence in crops
- First observations of global and seasonal terrestrial chlorophyll fluorescence from space
- Systematic Orbital Geometry-Dependent Variations in Satellite Solar-Induced Fluorescence (SIF) Retrievals
- Solar radiation and productivity in tropical ecosystems
- Climate and the efficiency of crop production in Britain
- New methods for measurements of photosynthesis from space
- Linking chlorophyll a fluorescence to photosynthesis for remote sensing applications: mechanisms and challenges
