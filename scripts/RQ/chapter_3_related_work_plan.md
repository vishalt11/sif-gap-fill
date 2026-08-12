# Chapter 3 Related Work - Proposed Structure

## Proposed chapter structure

- **3.1 SIF enhancement techniques** (target: approximately 2,000-2,500 words)
- **3.2 Crop-yield estimation from remote-sensing data**
- **3.3 Comparative analysis and research gap**

Section 3.1 can organize the literature into three primary-purpose groups. The categories are not completely exclusive because some products combine reconstruction, downscaling, gap filling, and de-noising. Hybrid papers should be assigned to the group that best represents their main contribution, with the overlap stated explicitly.

## Suggested Section 3.1 structure

### Opening taxonomy: reconstruction versus downscaling (approximately 250-350 words)

Define the distinction before reviewing individual studies:

- **Reconstruction and gap filling** estimate SIF where observations are spatially or temporally missing, usually on a fixed output grid. They learn relationships between observed SIF and spatially complete predictors or use the spatial/temporal covariance of SIF. Reconstruction may also extend a short satellite record backward in time.
- **Spatial downscaling** estimates the sub-pixel distribution of an observed or reconstructed coarse SIF signal at a finer grid. Some methods directly predict fine-grid values, whereas aggregation-constrained methods redistribute the coarse observation and preserve its parent-cell mean or total.
- **Other and hybrid enhancement approaches** include retrieval-noise reduction, temporal forecasting, and combinations of machine learning with geostatistics. These approaches address limitations other than, or in addition to, spatial resolution and missing coverage.

The introduction should stress that spatially detailed estimates are not automatically observations at that fine scale. Their reliability depends on the assumptions relating fine-resolution predictors to SIF and on the spatial support used for validation.

### Reconstruction and gap filling (approximately 650-800 words)

Representative papers include GOSIF, spatially continuous TanSat SIF, DOSIF, RTSIF, ROSIF, ST-LGBM, SIFoco2_005, and CSIF. The subsection should synthesize families of approaches instead of presenting one isolated paragraph per paper.

Discuss:

- the supervisory SIF sensor and why reconstruction is necessary;
- optical, meteorological, land-cover, spatial, and temporal predictors;
- model families such as Cubist, random forest, neural networks, boosted trees, and context-specific models;
- output resolution, temporal coverage, and retrospective extension;
- whether validation is random, temporally separated, spatially separated, or independent;
- limitations such as weak transfer to unobserved regions, biome-dependent performance, and the inability of reflectance-only predictors to represent rapid physiological stress.

Important thesis links are the use of sparse OCO-2 supervision, comparison of tree-based and neural-network models, the value of spatially separated evaluation, and variation among vegetation or agro-climatic regions.

### Spatial downscaling (approximately 800-1,000 words)

Representative papers include the Henan GOSIF downscaling study, the Duveiller et al. GOME-2 product, BCSIF, CNN downscaling of GOSIF, DSIFRFK_EA0.05, SIFnet, TroDSIF, HSIF, and CNSIF.

Organize the discussion around methodological contrasts:

- **Unconstrained statistical downscaling:** a model trained at coarse support is applied to fine-resolution predictors.
- **Aggregation-constrained redistribution:** fine-grid patterns are estimated while preserving the original coarse observation.
- **Bias- or residual-corrected methods:** information missing from the statistical prediction is reintroduced at the coarse scale.
- **Deep-learning approaches:** CNNs use spatial neighbourhoods and auxiliary imagery to represent within-pixel heterogeneity.
- **Physically guided approaches:** energy conservation, light-use efficiency, fluorescence escape probability, or physiological variables guide the allocation.

For each family, briefly cover the supervisory data, predictors, model, output scale, and most informative validation result. Give particular attention to whether reported performance demonstrates only coarse-scale consistency or genuine fine-scale accuracy. External validation against OCO-2/OCO-3 footprints, airborne CFIS, or tower SIF is more informative than comparison with the supervisory product or with GPP alone.

This subsection should connect directly to the thesis's Sentinel-2-to-OCO-2 footprint design, aggregation experiments, CNN-versus-tabular comparison, and crop-pure fine-resolution predictions.

### Other and hybrid enhancement approaches (approximately 300-450 words)

Use this subsection for studies whose main contribution does not fit cleanly into the first two groups:

- fine-grained OCO-2 SIF forecasting;
- ANN-based TROPOMI retrieval-precision enhancement;
- hybrid machine learning and kriging with external drift;
- coSIF cokriging with uncertainty quantification.

Focus on what additional problem each method solves: future-value prediction, reduction of retrieval noise, use of nearby SIF covariance, or formal uncertainty estimation. Explain that these methods can complement reconstruction or downscaling rather than replace them.

### Short synthesis and transition (approximately 150-250 words)

End Section 3.1 by identifying recurring gaps:

- random train/test splits may overestimate performance in genuinely unobserved regions;
- agreement after reaggregation does not prove fine-pixel accuracy;
- reconstructed products used as supervision can propagate upstream errors;
- external SIF validation is less common than same-source or GPP validation;
- few approaches combine very high-resolution Sentinel-2 information, explicit OCO-2 footprint support, agricultural masks, and aggregation-aware evaluation.

These gaps should lead naturally to the thesis methodology without presenting the thesis results in the literature-review chapter.

## What to report for each paper or paper family

The methodology should be summarized briefly, but the review should also evaluate the evidence. Include only the details needed for comparison:

1. enhancement objective and category;
2. supervisory SIF source;
3. main predictors and model family;
4. spatial/temporal output and study domain;
5. validation design and one or two representative metrics;
6. whether validation is internal SIF, external SIF, or proxy/application validation;
7. central assumption or limitation;
8. relevance to the thesis research questions.

Avoid listing every metric from `literature_review_references_table.md` in the prose. Use the table for complete factual coverage and use the narrative to explain methodological patterns, strengths, limitations, and research gaps.

## Proposed Section 3.3: Comparative analysis and research gap

An acronym is useful for readable comparison, but it should not be the only selection rule. Select a smaller set that represents different methodological choices and is especially relevant to the thesis. A strong candidate set is:

- GOSIF - tabular machine-learning reconstruction;
- CSIF - neural-network reconstruction from reflectance;
- ST-LGBM - spatially and temporally constrained reconstruction;
- BCSIF - bias-corrected downscaling;
- SIFnet - CNN downscaling and OCO-2/OCO-3 footprint validation;
- TroDSIF - aggregation-preserving redistribution;
- CNSIF - 500 m deep-learning reconstruction using Landsat/Sentinel-2 information;
- coSIF - geostatistical prediction with uncertainty quantification.

Optionally replace or supplement one of these with HSIF to include a physically constrained method.

Compare the selected studies using common dimensions: enhancement task, supervisory sensor, predictors, model family, output resolution, coarse-scale conservation, validation split, external validation, uncertainty reporting, and relevance to agricultural fine-scale SIF. Do not rank papers solely by R-squared or RMSE because their sensors, wavelengths, units, domains, periods, spatial supports, and validation datasets differ.

Section 3.3 should finish by integrating Sections 3.1 and 3.2: explain why a high-resolution, crop-aware, footprint-validated SIF product may be useful for yield estimation, and state the specific methodological gap addressed by the thesis.
