# Dynamic Treatment Regimes - Switch Project

This repository contains the code to replicate the target trial emulation for switching ventilation modes. 

The analysis evaluates dynamic treatment regimes comparing $PaO_2/FiO_2$ thresholds of **>150 vs >200 vs >250 mmHg**, stratified by whether PEEP is **<8 cmH₂O or ≥8 cmH₂O**.

### Primary Outcome
* **Restricted Mean Time Lost (RMTL)** of successful extubation within 28 days (672 hours).

---

## Required Data Format
The input dataset (`df_survival`) must be a longitudinal, **long-format** panel dataset grouped by `stay_id`.

### 1. Key Structuring & Core Columns
* **`stay_id`** *(integer/character)*: Unique identifier for each patient's ICU stay. This variable acts as the panel grouping variable used to compute lags, trajectories, and cluster-robust summaries.
* **`hour`** *(numeric)*: Time counter in hours, indicating the time elapsed since meeting the study's baseline eligibility criteria. 
* **`time`** *(POSIXct/datetime)*: The actual calendar timestamp of the observation, required to align absolute times across events.

### 2. Intervention & Censoring Anchors
* **`first_switch_time`** *(POSIXct/datetime)*: The absolute timestamp when the transition from controlled to assisted ventilation first occurred for a patient. 
* **`switched`** *(binary: 0 or 1)*: Indicator of whether the patient has been switched to assisted ventilation at that specific time point.
* **`nmb`** *(binary: 0 or 1)*: Neuromuscular blockade use, used to monitor adherence to protocols.

### 3. Competing Risk Event Indicators
The outcomes are evaluated using a multi-state competing risks framework via the **Aalen-Johansen estimator**. The target trial emulation requires these mutual hazards to be tracked at every interval step: 
* **`successfully_extubated`** *(binary: 0 or 1)*: The primary successful event indicator.
* **`died`** *(binary: 0 or 1)*: The competing risk event indicator.
* **`ltfu_discharged`** *(binary: 0 or 1)*: Lost to follow-up (e.g., early ICU discharge before extubation).

### 4. Dynamic Regulating Variables (Confounders/Triggers)
The emulation logic dynamically evaluates compatibility with the assigned treatment regimes based on physiological parameters recorded at each hour step:
* **`peep`** *(numeric)*: Positive End-Expiratory Pressure, used for regime stratification (<8 vs. ≥8 cmH₂O).
* **`PF_ratio`** *(numeric)*: The $PaO_2/FiO_2$ ratio, evaluated against the dynamic threshold rules.

### 5. Propensity Score Time-Varying Covariates
To accurately compute Inverse Probability of Compatibility Weights (IPCW) and mitigate time-dependent confounding, the dataset must provide sequential physiological and clinical state variables over time:
* **Respiratory Metrics**: `ph`, `pco2`, `rr_set` (set respiratory rate), `vt` (tidal volume), and `driving` (driving pressure).
* **Clinical Severity**: `cv_sofa` (Cardiovascular SOFA score, indicating hemodynamic stability).
* **Neurological/Sedation Depth**: `sedation` (Factor levels: `'Awake'`, `'Moderate sedation'`, `'Deep sedation'`).
* **Fixed Baseline Demographics**: `age`.
