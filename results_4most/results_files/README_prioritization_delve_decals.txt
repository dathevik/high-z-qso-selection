# High redshift quasars candidates list selected for 4MOST ChANGES using SED fitting (Mkrtchyan et al., in prep.) 
# 2024-08-05: Last update TM; TM will be on leave between 2024-08 and 2025-02 included, and will not reply emails. For any question refer to Chiara Mazzucchelli (chiara.mazzucchelli@mail.udp)


This dataset contains High redshift quasar candidates sample from Viktoria's sources, providing various measurements and fitting results for each source.

Columns Description:
- ls_id: ID number of the source from Viktoria's sample.
- ra: Right ascension (degrees).
- dec: Declination (degrees).
- mag_auto_g, mag_auto_r, mag_auto_i, mag_auto_z (AB): Magnitudes from DECALS DR10 or DELVE DR2 (refer to mag_type column for the reference).
- snr_g, snr_r, snr_i, snr_z: Signal to Noise ratios from DECALS DR10 or DELVE DR2 (refer to mag_type column for the reference).
- magerr_auto_g, magerr_auto_r, magerr_auto_i, magerr_auto_z: Magnitude errors from DECALS DR10 or DELVE DR2 (refer to mag_type column for the reference).
- mag_type: Indicates which griz magnitudes set was used for the source: DECALS DR10 or DELVE DR2.
- ypetromag, jpetromag, hpetromag, kspetromag (AB): YJHKs magnitudes from VHS.
- ypetromagerr, jpetromagerr, hpetromagerr, kspetromagerr: YJHKs magnitude errors from VHS.
- w1mpro, w2mpro (AB): Magnitudes from ALLWISE.
- w1sigmpro, w2sigmpro: Magnitude errors from ALLWISE.
- BD_chi2_min: Minimum chi-squared value fitted with brown dwarf templates.
- BD_chi2_template: Corresponding brown dwarf template of the minimum chi-squared value.
- QSO_chi2_min: Minimum chi-squared value fitted with quasar templates.
- QSO_z: Corresponding quasar photometric redshift.
- R_chi2_best: Ratio of QSO_chi2_min/ BD_chi2_min. R_chi2_best < 1 indicates a better fit with QSO templates.
- F_test_value: A statistical assessment where the test statistic follows an F-distribution under the null hypothesis.
- BIC_value: Bayesian Information Criterion evaluates a model's predictive likelihood.
- QSO_EBV: EBV of the QSO template (0, 0, 0.1, 0.25).
- QSO_EMline: Emission lines of QSO template (weak = -2, normal = 0, strong = 2).
- Data_Points: Number of photometric detections out of 10 photometric bands.
- QSO_Par: Type of QSO template (Par = 1 is from Selsing+2016, Par = 3 is from Temple+2021).
- ebv: Galactic E(B-V) in the direction of the source

This file contains detailed photometric and fitting data useful for various astronomical analyses. Ensure to refer to the appropriate columns for specific details.


