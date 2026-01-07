"""
PCA Score Statistical Analysis with General Linear Models
=========================================================
PCA Score Statistical Analysis Tool - Supports complex designs with multiple covariates

This tool performs statistical analysis on pseudobulk or sample-level data,
with PC score as the response variable.

Main features:
- OLS (Ordinary Least Squares) with robust standard errors
- LMM (Linear mixed model) with random effects
- Freedman-Lane permutation testing
- Estimated Marginal Means (EMM)

Application examples:
- scRNA-seq pseudobulk data: PC score and sex/cell type association
- Batch effect or age etc. covariate adjusted analysis
- Hierarchical model considering donor variation
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
from statsmodels.regression.mixed_linear_model import MixedLM
from patsy import dmatrix
import warnings
import zipfile
import io
from datetime import datetime
warnings.filterwarnings('ignore')

# Import helper functions for sample name processing
sys.path.insert(0, '/home/ichiro/streamlit')
from helper_func import remove_after_space, remove_sample_num

# Page configuration
st.set_page_config(page_title="PCA Statistical Analysis", page_icon="📊", layout="wide")

# Initialize session state for temp directory
if "pca_temp_dir" not in st.session_state:
    temp_dir = f"temp/pca_{int(time.time())}"
    os.makedirs(temp_dir, exist_ok=True)
    st.session_state.pca_temp_dir = temp_dir

st.title("📊 PCA Score Statistical Analysis")
st.markdown("""
### Generalized Linear Model Analysis with PC Score as Response Variable

This app statistically tests the effects of **sex, cell subtype, batch, age**, etc.
on PC scores obtained from PCA (Principal Component Analysis).

---
""")

with st.expander("📚 Feature Details (Click to Expand)", expanded=False):
    st.markdown("""
#### 📌 Main use cases:
- **scRNA-seq pseudobulk analysis**: PC score and phenotype association for sample-unit aggregated data
- **Unbalanced design supported**: Supports different group sizes and missing data
- **Hierarchical structure**: Consider hierarchy such as donor and technical batch
- **Covariate adjustment**: Adjust for age, depth, quality metrics, etc.

#### 🔬 Implemented Statistical Methods:

**1. OLS (Ordinary Least Squares)**
- Fixed effect model (sex, subtype, batch, age, etc.)
- **HC3 robust SE**: Standard error robust to heterogeneous variance (recommended)
- **Cluster-robust SE**: Considering within-cluster correlation (e.g., donor unit)
- Type II ANOVA: Significance test for each variable

**2. LMM (Linear Mixed Model)**
- **Random effect**: Hierarchical structure support like `(1|donor)`, `(1|batch)`
- **REML Estimation**: Unbiased variance estimation (recommended)
- **ML Estimation**: For model comparison (AIC/BIC)
- **ICC calculation**: Intraclass correlation coefficient (random effect contribution rate)
- Warning displayed when donor count is low (<5)

**3. Permutation Test**
- **Freedman-Lane method**: Test strictly controlling covariates (recommended)
- **Simple permutation**: Simple label permutation
- Effective for small samples when distribution assumptions are questionable
- Iteration count: 1,000 to 50,000 (customizable)

**4. EMM (Estimated Marginal Means)**
- Covariate-adjusted group mean estimation
- Visualization of subtype sex differences, etc.
- With confidence intervals for easier interpretation

#### ⚙️ Data Quality Checks:
- **Complete collinearity detection**: Warning for impossible estimation effects in advance
- **VIF calculation**: Multicollinearity diagnosis for numeric covariates
- **Sample size confirmation**: Sample count per group and cross-tabulation
- **Missing value handling**: Automatic removal and reporting

#### 📊 Visualization:
- **Forest plot**: Coefficient estimation value and 95% confidence interval
- **Diagnostic plot**: Residual, Q-Q plot, Scale-Location
- **EMM plot**: Estimated mean value by group
- **Permutation histogram**: Null distribution and observed statistic

#### 💾 Results Download:
- Coefficient table (CSV)
- ANOVA table (CSV)
- Variance component (LMM)
- Permutation distribution (CSV)
- All figures (can be saved from Streamlit)
""")

# =============================================================================
# Practical Guide
# =============================================================================
with st.expander("📖 Practical Guide: Sex Difference Axis Exploration in Cell Type Pseudobulk Data", expanded=False):
    st.markdown("""
### 🎯 Use Case: Finding "Sex Axis" in cell.type 4 types × sex 8 samples

#### **Input Data Prerequisite**

**1 Row = 1 Sample (pseudobulk)**

Column example:
```
sample_id    sex       cell.type    PC1      PC2      PC3      PC4      ...
sample_001   Female    art         -2.34     1.52    -0.87     0.45
sample_002   Male      art          1.87    -1.23     0.92    -0.35
sample_003   Female    cap1         0.23     2.11     1.34    -1.12
sample_004   Male      cap1        -1.45    -0.98     0.67     1.23
sample_005   Female    cap2         1.12     0.34    -1.56     0.89
sample_006   Male      cap2        -0.89     1.67     1.12    -0.67
sample_007   Female    vein         0.56    -1.34     0.23     1.78
sample_008   Male      vein        -1.23     0.45    -0.89     0.34
```

---

#### **App Settings (Column Mapping)**

1. **Main variable 1 (main_var)** = `sex`
2. **Main variable 2 / Blocking variable (blocking_var)** = `cell.type`
3. **Analysis PC** = Select PC1 to PCk sequentially (or create "sex composite axis" at once as described later)

---

#### **StepA：PCper （one dimension at a time）"Main effectscreening**

**Purpose:** which PCis sex best describeswhether, cell.type as covariatewhile adjusting determine。

**Steps:**

1. **PCselect one**（e.g.：PC3）

2. **in analysis tab OLSselect**and build Model formula as follows：
   ```
   PC3 ~ C(sex) + C(cell.type)
   ```
   - Interactionnot included（Main effectonly）
   - Standard Error is  **HC3 Robust** recommended（Small samplerobust to）

3. **Output Confirm**：
   - **sex Coefficient**（Female vs Male）and its  HC3robust SEof **p-value** confirm
   - Coefficientof**signnumber**and**Effect size**（coefficient absolute value）record

4. **Permutation Test on**for more accurate p-valueobtain：
   - **Freedman-Lanemethod** select（recommended）
   - strata (= cell.typemaintaining  sex label permutation
   - Small sample（n=8）also robust TestisPossible

5. **EMM（Estimated marginal means）plot**for consistent significance confirm：
   - cell.type for each  Female vs Male difference **whether signs are consistent**Confirm
   - Allin all cell typesSamedirectiondirectionofdifferenceare present、sexofMain effectis evidence of consistency

6. **PC1to PCk repeat**comprehensive evaluation as follows：
   - ✅ **Effect size**（sexcoefficient magnitude）
   - ✅ **consistent significance **（subtypealignment of signs per：EMMplot to confirm）
   - ✅ **p-value**（HC3 and Permutation ofbothdirection）

   → thisrathetotalgoand「**sextosaimocontributegivedoPC**」theselectbu

---

#### **📝 Statistical interpretation**

thisOLSModelis 、**supportedpresenttTest（Δ = Female − Male ofMean）** andthisqualityallyetcvalueis：

- **Degrees of freedom**: df = cell.typecount − 1 = 3
- **Testtarget**: cell.typethelayerandandcontroldoneUpperwithof sex ofMain effect
- **smallmarkthisalsorobust**: HC3 Robust SE and Permutation thecombineusedothatwith、n=8alsotrustrelycanruinfermeasureisPossible

**Notepoint:**
- cell.typeper toSamplecountisetcshii（baransudesign）case、saimostatisticaldetectionpowerisHigh
- Interaction `C(sex) * C(cell.type)` theenterrerucase、eachcell（sex × cell.type）tosailow3SampleRequired
- thisUse casewithis eachcellto1Samplezutsushikanotedme、Interactionis Estimationnotpossible（completepartseparate）

---

#### **🔍 Results readdirectionexample**

| PC  | sexCoefficient | HC3 p-value | Perm p-value | consistent significance  | totalgoevaluation |
|-----|---------|---------|----------|--------|----------|
| PC1 | 0.45    | 0.234   | 0.251    | ✅ All+ | weaki    |
| PC2 | -1.87   | 0.032   | 0.041    | ✅ All- | **strong** |
| PC3 | 0.23    | 0.678   | 0.712    | ± mixexist  | noneshi    |
| PC4 | 1.12    | 0.098   | 0.112    | ✅ All+ | middleapproximately  |

→ koexamplewithis  **PC2 is sex axis** andandsaimostrong（negativeofdirectiondirectiontopartseparate）

---

#### **💡 nextofSteps（responduse）**

1. **multiplecountPCtheunifygodone"sex score"ofmakebecome**
   - Effect sizeisLargePC（e.g.: PC2, PC4）thelinearresultgo
   - `sex_score = -1.87 × PC2 + 1.12 × PC4` ofyoutoweightattachkeMean
   - thiscomposite axisthenewedarespondanswervariablecountandandagainanalysis

2. **IndividualGenereberutoofdorirudaun**
   - sexaxis（e.g.: PC2）of loading isLargeGenetheextract
   - significance differencetocontributegivedoGenegroupandandinterpretation

3. **multipleTestcorrection**
   - multiplecountPC（e.g.: PC1to PC10）theexplorationdonecase、Benjamini-Hochberg FDRcorrectiontheApply
   - 10pieceofPCtheTestdocase、p < 0.05 → FDR < 0.05 tocorrection

---

#### **📚 thisUse casewithuseumachineability**

- ✅ OLS with HC3 Robust SE
- ✅ Freedman-Lane Permutation Test
- ✅ EMM（Estimated marginal means）plot
- ✅ diagnosisplot（Residualnormalisignificance ,etcVariancesignificance ）
- ✅ Coefficientte-buruofDownload（multiplecountPCofresulttheunifygo）

---
""")

# =============================================================================
# Practical guide 2: timesystemColumn × GeneType
# =============================================================================
with st.expander("📖 Practical guide2：WT vs KO × timesystemColumn3pointinInteractionanalysis", expanded=False):
    st.markdown("""
### 🎯 Use case：WT/KO × 3timepoint（eachn=3）for "genotypeeffectistimeBetweendependsonka」theTest

#### **InputDataprerequisite**

**1Row = 1Sample（independentSamplethinkset）**

- **design**: 2GeneType (WT/KO) × 3timepoint (t1/t2/t3) × eachcelln=3 = **total18Sample**
- **assumption**: eachtimepointwithDifferentpiecebody（independentSample）
- **ifsamepiecebodyofrepeated measuresif** → Lowerki「repeated measurescase of」thereferreference

Columnexample：
```
sample_id    genotype    time    donor_id    PC1      PC2      PC3      PC4      ...
WT_t1_1      WT          t1      D01        -1.23     0.45    -0.67     1.12
WT_t1_2      WT          t1      D02         0.87    -1.34     0.92    -0.45
WT_t1_3      WT          t1      D03        -0.45     1.23    -1.01     0.78
WT_t2_1      WT          t2      D04         1.56    -0.23     1.34    -0.89
...（followingLower、KO_t1, KO_t2, KO_t3mosamemanner）
```

---

#### **statisticalModelofbattleomit**

**important**: first**InteractionthemiddlehearttoTest** → Requiredifsimpleeffecttodigru

1. **InteractionisSignificant** → genotypeeffectis timepointwithDifferent（timeBetweendependson）
   - eachtimepointinWT vs KOtheIndividualtoTest（simpleeffect）
   - eachgenotypeInsideintimesystemColumnvariableizationtheTest

2. **Interactionisnon-Significant** → genotypeeffectis timepointtodependsonsezuoneset
   - InteractionthedropandandMain effectonlyofModelwithinterpretation
   - genotypeMain effectandtimeMain effecttheIndividualtoevaluation

---

#### **apuro-chiA：timethecategoryandandhandle（generalunderstand,recommended）**

**Modelformula:**
```
PC ~ C(genotype) * C(time)
```

**specialsign:**
- eachtimepointtheindependentdonecategoryandandhandle
- non-linearatimeBetweenvariableizationtomosupported
- 3timepointonlycase of、thisissaimosoftsoftwithsafeall

---

##### **🔧 apuriinSettingsStep（categoryversion）**

**1. DataPreparation**
- `genotype` Column: "WT", "KO"
- `time` Column: "t1", "t2", "t3" （characterColumnandand）
- `PC1`, `PC2`, ... Column

**2. ColumnMapping**
- **Main variable 1**: `genotype`
- **Main variable 2 / Blocking variable**: `time`
- **analysisPC**: PC1to PCk fromSelect
- **Donor/Subject ID**: （independentSampleifNot required）

**3. ModelSettings**
- **Modeltype**: OLS
- **Standard Error**: **HC3 Robust** （recommended、Small samplerobust to）
- **Interactionterm**: `genotype:time` theAdd ← **important！**

**4. TestOption**
- **ANOVA Type**: Type II（recommended）
- **Permutation Test**: ON（Freedman-Lanemethod、recommended）
  - n=18is Small sampleaofwithPermutation testwithsupplementstrong

**5. Visualization**
- **EMM（Estimated marginal means）**: ON
  - genotype and time ofbothdirectionselect
  - plotwithtimesystemColumntorendoofgroupBetweendifferenceconfirm
- **diagnosisplot**: ON（Residualnormalisignificance ,etcVariancesignificance Confirm）

---

##### **📊 Results readdirection（categoryversion）**

**Step 1: ANOVAtablewithInteractionconfirm**

| term | Sum Sq | df | F value | p value |
|----|--------|----|---------| --------|
| C(genotype) | 15.3 | 1 | 8.23 | 0.012 |
| C(time) | 45.7 | 2 | 12.34 | 0.001 |
| **C(genotype):C(time)** | **23.4** | **2** | **6.32** | **0.011** ← important |
| Residual | 22.3 | 12 | - | - |

**judgment:**
- **Interaction p=0.011 < 0.05** → **Significant**
- → **genotypeeffectis timepointwithDifferent**（timeBetweendependsonofKOeffect）

**Step 2: EMMtablewitheachtimepointofgroupBetweendifferenceconfirm**

| genotype | time | mean | SE | 95% CI lower | 95% CI upper |
|----------|------|------|----|--------------| -------------|
| WT | t1 | -0.45 | 0.23 | -0.95 | 0.05 |
| KO | t1 | -0.52 | 0.23 | -1.02 | -0.02 |
| WT | t2 | 0.87 | 0.23 | 0.37 | 1.37 |
| KO | t2 | 1.89 | 0.23 | 1.39 | 2.39 |
| WT | t3 | 1.34 | 0.23 | 0.84 | 1.84 |
| KO | t3 | 2.78 | 0.23 | 2.28 | 3.28 |

**interpretation:**
- **t1**: KO - WT = -0.07（hobodifferencenone）
- **t2**: KO - WT = 1.02（KOisincreaseadd）
- **t3**: KO - WT = 1.44（differenceissaratoexpandlarge）

→ **KOeffectis timeBetweenandandmotoincreaselarge**（t1withis differencenone → t3withlargekiadifference）

**Step 3: EMMplotwithVisualization**

timeBetweenthehorizontalaxis、PC scoretheverticalaxistoand：
- WT（blueline）is looseorkatoUpperrise
- KO（redline）is urgentintensetoUpperrise
- 2oflineisflatRowwithnot = Interactionpresent

---

##### **📝 simpleeffectofTest（InteractionisSignificantacase）**

InteractionisSignificantif、**eachtimepointinWT vs KOtheIndividualtoTest**：

**directionmethod1: sabusetoanalysis（handaction）**
1. Datathetimepointper topartdivide（t1only、t2only、t3only）
2. eachsabusetowith `PC ~ C(genotype)` theTest
3. **Bonferronicorrection**: p-valueofThresholdthe 0.05/3 = 0.0167 toSettings
   or  **BH-FDRcorrection**theApply

**directionmethod2: EMM contrast（recommended、generalcomerealinstalladvanceset）**
- selfactionallyeachtimepointinpeawaizucomparisonthecalculation
- multipleTestcorrectiontheselfactionApply

**presentexistofapuriinrealapplydirectionmethod:**
1. originofFiletheExcelwithtimepointper topartdivide（t1.tsv, t2.tsv, t3.tsv）
2. eachFiletheIndividualtouploadand `PC ~ C(genotype)` withTest
3. 3ofp-valuemanuactionwithBHcorrection

---

#### **apuro-chiB：timethecountValuetorendo（linear）andandhandle**

**Modelformula:**
```
PC ~ C(genotype) * time_numeric
```

**specialsign:**
- timeBetweentheContinuous variable（0, 1, 2）andandhandle
- **InteractionCoefficient = Slopeofdifference **（KOandWTwithtimeBetweenhookdistributeisdifferuka）
- 1ofCoefficientwithtimeBetweendependsonsignificance theneedabout canru
- **Note**: lineartheassumption（3pointshikanotofwith2nextis Overfittingrisuku）

---

##### **🔧 apuriinSettingsStep（countValueversion）**

**1. DataPreparation**
- `time_numeric` ColumntheAdd: t1→0, t2→1, t3→2
- or  `time` Columnthestraightconnect 0, 1, 2 tovariablefurther

**2. ColumnMapping**
- **Main variable 1**: `genotype`（category）
- **AddofCovariate**: `time_numeric`（countvalue and andadmitrecognizedoneru）
- **Interactionterm**: `genotype:time_numeric` theAdd

**3. ModelSettings**
- otheris apuro-chiAandSame

---

##### **📊 Results readdirection（countValueversion）**

**Coefficienttable:**

| term | Coef | SE | t | p value |
|----|------|----|---|---------|
| Intercept | -0.45 | 0.18 | -2.50 | 0.027 |
| C(genotype)[T.KO] | -0.07 | 0.25 | -0.28 | 0.783 |
| time_numeric | 0.90 | 0.15 | 6.00 | <0.001 |
| **C(genotype)[T.KO]:time_numeric** | **0.75** | **0.21** | **3.57** | **0.003** ← important |

**interpretation:**
- **time_numericCoefficient = 0.90**: WTgroupwithis timeBetween1unitaedriPC scoreis0.90increaseadd
- **InteractionCoefficient = 0.75**: KOgroupwithis timeBetween1unitaedri**sarato0.75manyku**increaseadd
  - tsumariKOgroupofSlope = 0.90 + 0.75 = **1.65**
  - WTgroupofSlope = 0.90
- **Interaction p=0.003** → **KOandWTwithtimeBetweentorendoisSignificanttoDifferent**

**benefitpoint:**
- 1ofTestfor "Slopeofdifference」theevaluation
- Effect sizeisstraightfeelteki（Slopeofdifferencepart）

**lackpoint:**
- linearassumption（3pointwithis inspectprovetroubledifficult）
- non-linearavariableizationtheseeescapesuPossiblesignificance 

---

#### **apuro-chiC：repeated measures（samepiecebodythechasetrace）case of**

**Beforeprovide:** SameDonor（D01, D02, D03）the3timepointwithmeasure

**Modelformula（LMM）:**
```
PC ~ C(genotype) * C(time) + (1|donor)
```

or timeBetweenSlopemorandamuto：
```
PC ~ C(genotype) * time_numeric + (1 + time_numeric|donor)
```

**Note:**
- DonorcountisFew（<5）case、randamuSlope `(1+time|donor)` is unstable
- firstis  `(1|donor)` fromstartmeru

---

##### **🔧 apuriinSettingsStep（LMMversion）**

**1. DataPreparation**
- `donor` Columnisrequired
- sameDonorIDismultiplecountRow（3Row）topresentreru

**2. ColumnMapping**
- **Main variable 1**: `genotype`
- **Main variable 2**: `time`
- **Donor/Subject ID**: `donor` ← **important！**
- **Interactionterm**: `genotype:time`

**3. ModelSettings**
- **Modeltype**: **LMM**
- **Estimationmethod**: REML（recommended）
- otherofSettingsis OLSandsamemanner

---

##### **📊 Results readdirection（LMMversion）**

**Fixed effecttable:**

| term | Coef | SE | z | p value |
|----|------|----|---|---------|
| Intercept | -0.45 | 0.28 | -1.61 | 0.108 |
| C(genotype)[T.KO] | -0.07 | 0.40 | -0.18 | 0.860 |
| C(time)[T.t2] | 1.32 | 0.20 | 6.60 | <0.001 |
| C(time)[T.t3] | 1.79 | 0.20 | 8.95 | <0.001 |
| C(genotype)[T.KO]:C(time)[T.t2] | 1.09 | 0.28 | 3.89 | <0.001 |
| C(genotype)[T.KO]:C(time)[T.t3] | 1.51 | 0.28 | 5.39 | <0.001 |

**Random effect:**
- `Var(donor)` = 0.45（DonorBetweenofbaratsuki）
- `Var(Residual)` = 0.23（measureInsideofbaratsuki）
- **ICC = 0.45 / (0.45 + 0.23) = 0.66**
  - 66%ofvariableactionisDonorBetweenofpiecebodydifferencewithDescriptiondoneru
  - repeated measuresofCorrelationisstrong → LMMissuitablecut

**interpretation:**
- InteractionisSignificant → genotypeeffectis timeBetweendependson
- Donorvariation betweenthesuitablecuttothinkconsiderdoneUpperinresulttheory

---

#### **🔍 reportmakemethod（recommended）**

**1. Interactionofresulttheorythesaiexcellentahead**
```
「genotype × time ofInteractionisSignificantwithaed（F(2,12)=6.32, p=0.011, HC3 robust SE;
Freedman-Lane permutation p=0.015）。thisis 、KOeffectistimeBetweendependsonwithexistthattheshowsu。」
```

**2. EMMfigurecreate**
- timeBetweenthehorizontalaxis、PC scoretheverticalaxis
- WTandKOofMean ± 95%CI thefoldrelinewith
- Possibleif **Δ(t) = KO - WT** ofdifferencepartplottheAdd

**3. simpleeffect（InteractionisSignificantacase）**
```
「matterAfterTestandandeachtimepointinWT vs KOthecomparisondoneresult：
- t1: Δ=-0.07 (95%CI: -0.68, 0.54), p=0.78（Significantdifferencenone）
- t2: Δ=1.02 (95%CI: 0.41, 1.63), p=0.004（Significant）
- t3: Δ=1.44 (95%CI: 0.83, 2.05), p<0.001（Significant）
BonferronicorrectionAftermo t2, t3 is Significantthefiberhold。」
```

**4. Effect sizethereporttell**
- eachtimepointofdifferenceestimationvalue and 95%CI
- Overallof genotype effect（categoryversionifη², countValueversionifSlopeofdifference）

**5. Permutation pthecombineki**
- HC3andarrangegosurebarobustsignificance ofprovebasis

---

#### **⚠️ Notematterterm**

**1. Samplesize**
- eachcelln=3is saismalllimit（statisticaldetectionpowerisLow）
- Possibleif n≥5 recommended
- **Permutation testwithsupplementstrong**dothatthestrongkurecommended

**2. multipleTestcorrection**
- simpleeffect（3comparison）tois  **Bonferroni** or  **BH-FDR** theApply
- multiplecountPCtheTestdocasemosamemanner

**3. linear vs category**
- **3pointonlyiflinearassumptionis dangerdanger** → categoryversionrecommended
- linearversionis  exploratory tousei、categoryversionwithConfirm

**4. completepartseparate（Perfect Separation）**
- Interactiontermtheenterreruand 2×3=6cell × n=3 = 18Sample
- df = 18 - 6 - 1 = 11（remainderaffluentis Few）
- lackdamageorOutlierisexistandEstimationunstable → diagnosisplot to confirm

**5. repeated measurescase of**
- DonorcountisFew（<5）and `(1+time|donor)` is unstable
- firstis  `(1|donor)` withTest、ConvergencesurebaSlopemotrysu

---

#### **📚 thisUse casewithuseumachineability**

- ✅ OLS with HC3 Robust SE（independentSample）
- ✅ LMM with (1|donor)（repeated measures）
- ✅ Type II ANOVA（InteractionofTest）
- ✅ Freedman-Lane Permutation Test（Small sampleofsupplementstrong）
- ✅ EMM（Estimated marginal means）plot（timesystemColumntorendoofVisualization）
- ✅ diagnosisplot（assumptionofConfirm）
- ⚠️ simpleeffectofselfactioncalculation（generalcomerealinstalladvanceset、presentexistis handactionwithsabusetoanalysis）

---

#### **💡 realpracticee.g.ofmaandme**

| statesituation | recommendedModel | Interactionterm | Testmethod | EMM |
|------|-----------|----------|--------|-----|
| independentSample、category | `PC ~ C(genotype) * C(time)` | present | OLS + HC3 + Perm | genotype × time |
| independentSample、linear | `PC ~ C(genotype) * time_numeric` | present | OLS + HC3 + Perm | Slopeofdifference |
| repeated measures | `PC ~ C(genotype) * C(time) + (1\|donor)` | present | LMM (REML) | genotype × time |

---

**hiandthatwith:**
1. **firstInteraction `genotype × time` theTest**
2. **Significantifeachtimepointofsimpleeffectto**（presentexistis handactionsabusetoanalysis）
3. **non-SignificantifMain effectModelwithinterpretation**
4. **Small sample（n=3/cell）aofwith Permutation withsupplementstrongrequired**

thisflowrewith、WT/KO × 3timepoint × n=3 andiupresentrealtekiarulemodelalso、**timeBetweendependsonofgenotypeeffect**thestatisticalallyevaluationcaning。

---
""")

st.markdown("---")

# =============================================================================
# Helper Functions
# =============================================================================

def remove_common_suffix(strings):
    """endtailoftogetherthroughneedelementexcludingremove

    DESeq2-LRT.pyandSameAlgorithm
    """
    if not strings or len(strings) == 0:
        return []
    # saimoShortcharacterColumn longsaobtain
    min_length = min(len(s) for s in strings)
    # togetherthroughofendtailPartialoflongsafind 
    suffix_length = 0
    for i in range(1, min_length + 1):
        suffix = strings[0][-i:]
        if all(s.endswith(suffix) for s in strings):
            suffix_length = i
        else:
            break
    # togetherthroughofendtailPartialisseetsufromnotcaseis originofrisutothereturnsu
    if suffix_length == 0:
        return strings
    # togetherthroughofendtailPartialtheDeleteandNewrisutocreate
    return [s[:-suffix_length] for s in strings]

def detect_collinearity(df, formula_terms):
    """completecollinearity detection

    Categorical variableBetweenwithcompletetotogetherroledoing（onedirectionisdecidemarebaotherdirectionmodecidemaru）
    casethedetectiondo。thisyouacase、Modelis Estimationcannot。

    e.g.: existsubtypeisAllFemalecase of、sexeffectandsubtypeeffectis partseparatecannot。
    """
    issues = []

    # Categorical variableonlythecheck
    cat_vars = [col for col in formula_terms if col in df.columns and df[col].dtype == 'object']

    for i, var1 in enumerate(cat_vars):
        for var2 in cat_vars[i+1:]:
            cross_tab = pd.crosstab(df[var1], df[var2])
            # eachRowor eachColumnto1tsushikanon-zeroneedelementisnotcaseis completetogetherrole
            if (cross_tab > 0).sum(axis=0).min() == 1 or (cross_tab > 0).sum(axis=1).min() == 1:
                issues.append(f"⚠️ **completetogetherrole (Perfect confounding)**: `{var1}` and `{var2}` iscompletetoonereachdoinging")

    return issues

def calculate_vif(df, numeric_cols):
    """VIF (Variance Inflation Factor) calculation

    VIF > 10: strongmultiplecollinearity present
    VIF > 5: middleapproximatelyofmultiplecollinearity present
    VIF < 5: questiontopicnone
    """
    from statsmodels.stats.outliers_influence import variance_inflation_factor

    vif_data = pd.DataFrame()
    vif_data["Variable"] = numeric_cols
    vif_data["VIF"] = [variance_inflation_factor(df[numeric_cols].values, i)
                       for i in range(len(numeric_cols))]
    return vif_data

def freedman_lane_permutation(y, X_full_df, interest_var_name, n_perm=10000, random_state=42,
                               stratify_by=None, one_sided=None, design_info=None):
    """Freedman-LanemethodtoyoruPermutation test（firmrobustversion）

    Covariatethecontrolshiaisra、specificofvariablecount（e.g.: sex）ofeffectthetested。

    Step:
    1. fullModel（Allofvariablecount）andnullModel（TesttargetfollowingOutside）thefuito
    2. nullModelofResidualpermute（stratifiedOptionpresent）
    3. placechangedoneedResidualfromquasisimilarrespondanswervariablecountthestructurebuild
    4. fullModeltheagainfuitoandStatisticthecalculation
    5. Iterationandnulldistributioncreate

    Parameters:
    -----------
    y : array-like
        respondanswervariablecount (PC scores)
    X_full_df : DataFrame
        fulldesignRowColumn（ColumnNamewith）
    interest_var_name : str
        TesttargetofvariablecountName（e.g.: 'sex'）。thisvariablecountincluding AllofColumnexcludingOutsideandnullModelcreate
    n_perm : int
        placechangetimescount
    random_state : int
        Random seed
    stratify_by : array-like, optional
        stratifiedvariablecount（e.g.: subtype）。pointsetdoand、eachlayerInsidewithonlyResidualpermute
    one_sided : str, optional
        one-sidedTestofdirectiondirection。'greater' (obs > null) or 'less' (obs < null)
        Nonecase ofis two-sidedTest
    design_info : patsy.DesignInfo, optional
        designRowColumn inforeport（firmrobustaColumnspecificfor）

    Returns:
    --------
    dict: observed (observedStatistic), null_distribution (nulldistribution), p_value,
          method_info (Reproducibilityforofinforeport)
    """
    import streamlit as st
    np.random.seed(random_state)

    # fullModelthefuito
    model_full = sm.OLS(y, X_full_df).fit()

    # Testtargetofvariablecountincluding Columnthespecific（design_infousewithfirmrobustization）
    if design_info is not None:
        # patsy design_infotheuseedfirmrobustaColumnspecific
        interest_cols = []
        for term in design_info.terms:
            term_name = term.name()
            # Testtargetvariablecountincluding ta-mu（Main effect + Interaction）thespecific
            if term_name != 'Intercept' and interest_var_name in term_name:
                # thista-musupporteddoColumnobtain
                slice_obj = design_info.term_name_slices[term_name]
                cols = X_full_df.columns[slice_obj]
                interest_cols.extend(cols)
    else:
        # fuo-rubaku: characterColumnMatching（Afterdirectionmutualchangesignificance ）
        interest_cols = [col for col in X_full_df.columns
                         if f'C({interest_var_name})' in col and col != 'Intercept']

    if not interest_cols:
        raise ValueError(f"Testtargetvariablecount '{interest_var_name}' supporteddoColumnisseetsukarinot")

    # observedStatistic（Main effectofColumn cutagainsttValueuse - firmrobustization）
    main_effect_col = [col for col in interest_cols if ':' not in col]
    if not main_effect_col:
        raise ValueError(f"Main effectofColumnisseetsukarinot: {interest_cols}")

    # absolute valuetStatisticuse（signnumberdependsontheavoidkeru）
    obs_t = float(model_full.tvalues[main_effect_col[0]])
    obs_stat = np.abs(obs_t)

    # nullModel（TesttargetvariablecountofAllofColumnexcludingOutside）
    X_reduced = X_full_df.drop(columns=interest_cols)

    if X_reduced.shape[1] == 1:  # Interceptonly
        residuals = y - np.mean(y)
        fitted_reduced = np.full_like(y, np.mean(y))
    else:
        model_reduced = sm.OLS(y, X_reduced).fit()
        residuals = np.array(model_reduced.resid)  # Convert to numpy array
        fitted_reduced = np.array(model_reduced.predict())  # Convert to numpy array


    # nullModelandagainststandModelofformulatheSave（reportuse）
    reduced_formula = "y ~ " + " + ".join([col for col in X_reduced.columns if col != 'Intercept'])
    if reduced_formula == "y ~ ":
        reduced_formula = "y ~ 1"  # Interceptonly
    full_formula = "y ~ " + " + ".join([col for col in X_full_df.columns if col != 'Intercept'])

    # placechange
    null_dist = []
    null_dist_abs = []  # absolute valuemoSave
    for _ in range(n_perm):
        # Residualpermute（stratifiedOption）
        if stratify_by is not None:
            # stratified permutation：eachlayerInsidewithonlyResidualshufflefull
            perm_residuals = np.zeros_like(residuals)
            for stratum in np.unique(stratify_by):
                stratum_mask = (stratify_by == stratum)
                stratum_indices = np.where(stratum_mask)[0]
                perm_indices = np.random.permutation(stratum_indices)
                perm_residuals[stratum_indices] = residuals[perm_indices]
        else:
            # non-stratified permutation：Overallshufflefull
            perm_idx = np.random.permutation(len(residuals))
            perm_residuals = residuals[perm_idx]

        y_perm = fitted_reduced + perm_residuals

        # fullModeltheagainfuito
        model_perm = sm.OLS(y_perm, X_full_df).fit()
        perm_t = float(model_perm.tvalues[main_effect_col[0]])
        null_dist.append(perm_t)
        null_dist_abs.append(np.abs(perm_t))

    null_dist = np.array(null_dist)
    null_dist_abs = np.array(null_dist_abs)

    # p-valuecalculation（one-sided/two-sided）+ continuouscontinuesignificance correction
    if one_sided == 'greater':
        # one-sided（Largedirectiondirection）: obs_tuse
        p_value = (np.sum(null_dist >= obs_t) + 1) / (n_perm + 1)
    elif one_sided == 'less':
        # one-sided（Smalldirectiondirection）: obs_tuse
        p_value = (np.sum(null_dist <= obs_t) + 1) / (n_perm + 1)
    else:
        # two-sidedTest: absolute valuewithcomparison
        p_value = (np.sum(null_dist_abs >= obs_stat) + 1) / (n_perm + 1)

    return {
        'observed': obs_t,  # originoftValue（signnumberwith）thereturnsu
        'observed_abs': obs_stat,  # absolute valuemoreturnsu
        'null_distribution': null_dist,  # signnumberwithtValueofdistribution
        'null_distribution_abs': null_dist_abs,  # absolute valueofdistribution
        'p_value': p_value,
        'stratified': stratify_by is not None,
        'one_sided': one_sided,
        'interest_cols': interest_cols,
        'main_effect_col': main_effect_col[0],
        'method_info': {  # Reproducibilityforofinforeport
            'n_permutations': n_perm,
            'random_seed': random_state,
            'test_statistic': '|t|' if one_sided is None else 't',
            'continuity_correction': True,
            'reduced_model': reduced_formula,
            'full_model': full_formula,
            'stratification_var': 'None' if stratify_by is None else 'Provided'
        }
    }

def calculate_emm(model, df, factors, numeric_vars_at_mean=True):
    """Estimated marginal means (Estimated Marginal Means) calculation

    CovariatethespecificofValue（normalis Mean）fixed toand、
    mainneedacausechildoflevelfor each PredictionMeanValuethecalculationdo。

    e.g.: agetheMeanValuefixed toand、significance by×sabutypeofgroupmigowasefor each 
        PC scoreofPredictionMeanthecalculation

    Parameters:
    -----------
    model : fitted statsmodels model
    df : DataFrame
        originData
    factors : list
        EMM thecalculationdocausechildofrisuto
    numeric_vars_at_mean : bool
        countValuevariablecounttheMeanValuefixed todoka

    Returns:
    --------
    DataFrame: eachgroupmigowaseofPredictionMean、SE、95%CI
    """
    from itertools import product

    # referreferenceguridoofmakebecome
    ref_grid = {}

    for col in df.columns:
        if col in factors:
            ref_grid[col] = df[col].unique()
        elif pd.api.types.is_numeric_dtype(df[col]):
            if numeric_vars_at_mean:
                ref_grid[col] = [df[col].mean()]
            else:
                ref_grid[col] = df[col].unique()

    # allgroupmigowasegenerate
    grid_combos = list(product(*[ref_grid[col] for col in df.columns if col in ref_grid]))

    # PredictionuseDatafure-mucreate
    pred_df = pd.DataFrame(grid_combos, columns=[col for col in df.columns if col in ref_grid])

    # Prediction
    predictions = model.get_prediction(pred_df)
    pred_summary = predictions.summary_frame(alpha=0.05)

    result_df = pred_df.copy()
    result_df['mean'] = pred_summary['mean']
    result_df['se'] = pred_summary['mean_se']
    result_df['ci_lower'] = pred_summary['mean_ci_lower']
    result_df['ci_upper'] = pred_summary['mean_ci_upper']

    return result_df

def plot_forest(coef_df, title="Forest Plot"):
    """CoefficientofForest plotmakebecome"""
    fig, ax = plt.subplots(figsize=(8, max(6, len(coef_df) * 0.5)))

    y_pos = np.arange(len(coef_df))
    ax.errorbar(coef_df['Coef'], y_pos,
                xerr=[coef_df['Coef'] - coef_df['CI_lower'],
                      coef_df['CI_upper'] - coef_df['Coef']],
                fmt='o', markersize=8, capsize=5, capthick=2, color='steelblue')

    ax.axvline(0, color='red', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(coef_df.index)
    ax.set_xlabel('Coefficient Estimate (95% CI)', fontsize=11)
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    plt.tight_layout()

    return fig

def plot_diagnostic(model, title="Diagnostic Plots"):
    """Modeldiagnosisplot

    4ofdiagnosisplotcreate:
    1. Residuals vs Fitted: non-linearsignificance 、etcVariancesignificance ofcheck
    2. Q-Q plot: Residualnormalisignificance check
    3. Scale-Location: etcVariancesignificance ofcheck（marklevelizationResidual）
    4. Residual histogram: Residualofdistribution
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    fitted = model.fittedvalues
    resid = model.resid

    # Residuals vs Fitted
    axes[0, 0].scatter(fitted, resid, alpha=0.6, s=50, edgecolors='k', linewidths=0.5)
    axes[0, 0].axhline(0, color='red', linestyle='--', linewidth=1.5)
    axes[0, 0].set_xlabel('Fitted values', fontsize=11)
    axes[0, 0].set_ylabel('Residuals', fontsize=11)
    axes[0, 0].set_title('Residuals vs Fitted', fontsize=12, fontweight='bold')
    axes[0, 0].grid(alpha=0.3)

    # Q-Q plot
    stats.probplot(resid, dist="norm", plot=axes[0, 1])
    axes[0, 1].set_title('Normal Q-Q Plot', fontsize=12, fontweight='bold')
    axes[0, 1].grid(alpha=0.3)

    # Scale-Location
    standardized_resid = resid / np.std(resid)
    axes[1, 0].scatter(fitted, np.sqrt(np.abs(standardized_resid)),
                      alpha=0.6, s=50, edgecolors='k', linewidths=0.5)
    axes[1, 0].set_xlabel('Fitted values', fontsize=11)
    axes[1, 0].set_ylabel('√|Standardized residuals|', fontsize=11)
    axes[1, 0].set_title('Scale-Location', fontsize=12, fontweight='bold')
    axes[1, 0].grid(alpha=0.3)

    # Residuals histogram
    axes[1, 1].hist(resid, bins=min(30, len(resid)//3),
                   edgecolor='black', alpha=0.7, color='steelblue')
    axes[1, 1].set_xlabel('Residuals', fontsize=11)
    axes[1, 1].set_ylabel('Frequency', fontsize=11)
    axes[1, 1].set_title('Residual Distribution', fontsize=12, fontweight='bold')
    axes[1, 1].grid(alpha=0.3, axis='y')

    plt.suptitle(title, fontsize=14, fontweight='bold', y=1.00)
    plt.tight_layout()
    return fig

def create_results_zip(pc_col):
    """
    analysisresultandgurafuthezipFiletomaandmeru

    Parameters:
    -----------
    pc_col : str
        PCkaramuName（FileNametouse）

    Returns:
    --------
    bytes : zipFileofbainariData
    """
    zip_buffer = io.BytesIO()

    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base_name = f"pca_stats_{pc_col}_{timestamp}"

        # 0. WelchTypetTest results
        if 'welch_result' in st.session_state:
            welch_res = st.session_state.welch_result
            welch_df = pd.DataFrame({
                'test': ['Welch\'s t-test'],
                'group1': [welch_res['group1']],
                'group2': [welch_res['group2']],
                'group1_mean': [welch_res['group1_mean']],
                'group2_mean': [welch_res['group2_mean']],
                't_value': [welch_res['t']],
                'p_value': [welch_res['p']],
                'cohens_d': [welch_res['cohens_d']],
                'hedges_g': [welch_res['hedges_g']]
            })
            zip_file.writestr(f"{base_name}/tables/welch_ttest.csv", welch_df.to_csv(index=False))

        # 1. OLSCoefficienttable
        if 'ols_coef' in st.session_state:
            csv_data = st.session_state.ols_coef.to_csv()
            zip_file.writestr(f"{base_name}/tables/ols_coefficients.csv", csv_data)

        # 2. LMMCoefficienttable
        if 'lmm_coef' in st.session_state:
            csv_data = st.session_state.lmm_coef.to_csv()
            zip_file.writestr(f"{base_name}/tables/lmm_coefficients.csv", csv_data)

        # 3. Permutation testresult
        if 'perm_result' in st.session_state:
            perm_result = st.session_state.perm_result

            # Statisticsummary
            perm_summary = pd.DataFrame({
                'observed_stat': [perm_result['observed']],
                'p_value': [perm_result['p_value']],
                'n_permutations': [len(perm_result['null_distribution'])],
                'method': [perm_result.get('method', 'unknown')],
                'stratified': [perm_result.get('stratified', False)]
            })
            zip_file.writestr(f"{base_name}/tables/permutation_summary.csv",
                            perm_summary.to_csv(index=False))

            # nulldistribution
            perm_df = pd.DataFrame({
                'permutation_stat': perm_result['null_distribution']
            })
            zip_file.writestr(f"{base_name}/tables/permutation_null_distribution.csv",
                            perm_df.to_csv(index=False))

        # 4. ANOVAtable（OLS）
        if 'ols_anova' in st.session_state:
            csv_data = st.session_state.ols_anova.to_csv()
            zip_file.writestr(f"{base_name}/tables/ols_anova.csv", csv_data)

        # 5. EMMresult
        if 'ols_emm' in st.session_state:
            csv_data = st.session_state.ols_emm.to_csv(index=False)
            zip_file.writestr(f"{base_name}/tables/estimated_marginal_means.csv", csv_data)

        # 6. gurafuofSave
        # Forest plot (OLS)
        if 'fig_forest_ols' in st.session_state:
            img_buffer = io.BytesIO()
            st.session_state.fig_forest_ols.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/ols_forest_plot.png", img_buffer.getvalue())

            # PDFversionmoSave
            pdf_buffer = io.BytesIO()
            st.session_state.fig_forest_ols.savefig(pdf_buffer, format='pdf', bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/ols_forest_plot.pdf", pdf_buffer.getvalue())

        # Diagnostic plot
        if 'fig_diagnostic' in st.session_state:
            img_buffer = io.BytesIO()
            st.session_state.fig_diagnostic.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/ols_diagnostic_plots.png", img_buffer.getvalue())

            pdf_buffer = io.BytesIO()
            st.session_state.fig_diagnostic.savefig(pdf_buffer, format='pdf', bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/ols_diagnostic_plots.pdf", pdf_buffer.getvalue())

        # EMM plot
        if 'fig_emm' in st.session_state:
            img_buffer = io.BytesIO()
            st.session_state.fig_emm.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/estimated_marginal_means.png", img_buffer.getvalue())

            pdf_buffer = io.BytesIO()
            st.session_state.fig_emm.savefig(pdf_buffer, format='pdf', bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/estimated_marginal_means.pdf", pdf_buffer.getvalue())

        # Forest plot (LMM)
        if 'fig_forest_lmm' in st.session_state:
            img_buffer = io.BytesIO()
            st.session_state.fig_forest_lmm.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/lmm_forest_plot.png", img_buffer.getvalue())

            pdf_buffer = io.BytesIO()
            st.session_state.fig_forest_lmm.savefig(pdf_buffer, format='pdf', bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/lmm_forest_plot.pdf", pdf_buffer.getvalue())

        # Null distribution plot (Permutation test)
        if 'fig_null_dist' in st.session_state:
            img_buffer = io.BytesIO()
            st.session_state.fig_null_dist.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/permutation_null_distribution.png", img_buffer.getvalue())

            pdf_buffer = io.BytesIO()
            st.session_state.fig_null_dist.savefig(pdf_buffer, format='pdf', bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/permutation_null_distribution.pdf", pdf_buffer.getvalue())

        # READMEmakebecome
        readme_content = f"""PCA Statistics Analysis Results
================================

Analysis Date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
PC Column: {pc_col}

Directory Structure:
-------------------
tables/          - CSV files with statistical results
  - ols_coefficients.csv           : OLS regression coefficients
  - lmm_coefficients.csv           : LMM fixed effects coefficients
  - ols_anova.csv                  : ANOVA table
  - estimated_marginal_means.csv   : EMM estimates
  - permutation_summary.csv        : Permutation test summary
  - permutation_null_distribution.csv : Null distribution values

figures/         - Publication-quality figures (PNG 300dpi + PDF)
  - ols_forest_plot.*              : OLS coefficient forest plot
  - lmm_forest_plot.*              : LMM coefficient forest plot
  - ols_diagnostic_plots.*         : Model diagnostic plots
  - estimated_marginal_means.*     : EMM visualization
  - permutation_null_distribution.*: Permutation test null distribution

Generated by: PCA Statistics Streamlit App
"""
        zip_file.writestr(f"{base_name}/README.txt", readme_content)

    zip_buffer.seek(0)
    return zip_buffer.getvalue()

# =============================================================================
# Sidebar Options
# =============================================================================

with st.sidebar:
    st.title("⚙️ analysisOption")

    st.markdown("### Statistical methodsofSelect")
    analysis_method = st.radio(
        "analysishandmethod:",
        ["OLS (Ordinary Least Squares)",
         "LMM (Linear mixed model)",
         "bothdirection (OLS + LMM)"],
        index=0,
        help="OLS: Fixed effectonly。LMM: Random effect（donor, batchetc）including floorlayerModel"
    )

    if "OLS" in analysis_method:
        st.markdown("#### OLS Option")
        se_type = st.radio(
            "Standard erroroftype:",
            ["Classical (normal)", "HC3 (robustEstimation,recommended)", "Cluster-robust (cluster-robust)"],
            index=1,
            help="HC3: heterogeneous Variancerobust to。Cluster-robust: within-cluster Correlationconsidering "
        )

        if se_type == "Cluster-robust (cluster-robust)":
            st.info("📌 cluster-variablecountis DatauploadAftertoSelectdo")

        # ANOVA Type selection
        anova_type = st.radio(
            "ANOVA Type:",
            ["Type II (recommended)", "Type III"],
            index=0,
            help="Type II: each variable Main effecttheother variablescountwithadjustmentandTest（Interactionnonecase oftorecommended）\nType III: AllofInteractiontheincludemeteTest（Interactionpresentcase oftorecommended）"
        )

        # WLS option
        use_wls = st.checkbox(
            "WLS (addweightsaismalltworidemethod) use",
            value=False,
            help="Varianceisheterogeneous acaseto、eachobservedValuetoweightthetsuketeEstimationdo。weightis ResidualofinverseVariancewithselfactioncalculationdoneing"
        )

        if use_wls:
            with st.expander("⚠️ WLSuseUpperofimportantaNote（necessaryread）"):
                st.markdown("""
                ### WLS (Weighted Least Squares) ofsuitablecutause

                **⚠️ important**: WLSis "weight"ofdecidemedirectionnextordinalwithEstimationisunstabletoariorsuku、**Small samplewithis recommendeddonenot**

                ---

                ### 🚫 WLStheavoidkerushouldstatesituation（important）

                #### **Small samplesettotalwithis basethisallyNG**

                e.g.: subtype=4group、eachgroupn=1-2approximatelyofextremesmallsettotal
                - ❌ **rootbasisexistweightismakeritokui**（VarianceEstimationisunstable）
                - ❌ **Estimationisbureorsui**（weightofnotconfirmrealsignificance isHigh）
                - ❌ **oversurplussuitablegoofrisuku**（Degrees of freedomisFew）

                **substitutewaritouseushoulddirectionmethod**:
                - ✅ **WelchType + Type II ANOVA** (mainanalysis)
                - ✅ **OLS + HC3robustSE** (Coefficient,95%CI)
                - ✅ **Freedman-LanePermutation test** (Conditionwithnulldistribution)

                thisrais **weighttheassumptionsezu**differVariancerobust towith、smallmarkthisalsooverdegreetoburetokuiis。

                ---

                ### ✅ WLSisValidastatesituation（limitsetteki）

                **tenpartaIterationispresent、groupBetweenVariancedifferenceisclearcleartoEstimationcanrucaseonly**

                #### RequiredCondition（AllsatisfysuRequiredpresent）:

                1. **each group Samplesizeistenpart**
                   - sailowalsoeachgroupn ≥ 10-15
                   - VariancethesafesetallyEstimationcanru

                2. **groupBetweenwithVariancedifferenceisclearcleartosonexist**
                   - LeveneTestorBartlettTestwithdetectionPossible
                   - diagnosisplotwithviewfeelallyclearraka

                3. **weightofrootbasisisclearconfirm**
                   - groupInsideinverseVariance（σ²ᵢ）ismatterBeforeinforeportfromalreadyknow
                   - skilltechniquetekimultiplemakecountetcOutsidepartinforeportbased on

                4. **Modelassumptionissatisfydoneru**
                   - linearsignificance 、independentsignificance 
                   - weightwithAfterofResidualiscorrectruledistribution

                #### kocase ofofrecommended:

                - **mainanalysisweathersupplement**: WLS（inverseVarianceweight、weightexplicitly）
                - **HC3is normalNot required**（WLSwithetcVarianceizationcanteireba）
                - Requiredif**"reference"**andandcombineki（Sensitivityanalysis）

                ---

                ### 📊 useipartkeofpointneedle

                #### aaedofsettotalisSmall sample（4group×eachn=1-2etc）case of:

                ```
                【mainanalysis】
                ✅ WelchType + Type II ANOVA
                   PC_k ~ C(sex) + C(subtype)
                   → differVariancerobust、smallmarkthisalsosafeset

                【combineki】
                ✅ OLS + HC3robustSE
                   → CoefficientEstimation,95%CI
                   → p-valueis robust

                ✅ Freedman-LanePermutation test
                   → Conditionwithnulldistribution
                   → distributionofassumptionnone

                【Sensitivityanalysisonly】
                ⚠️ WLS
                   → weight=groupInsideinverseVarianceetcclearshow
                   → mainanalysistois usewanot
                ```

                #### tenpartaIterationisexistcase（eachgroupn ≥ 15etc）:

                ```
                【mainanalysisweathersupplement】
                ✅ WLS（inverseVarianceweight）
                   → effectratetekiEstimation
                   → weightexplicitly

                【Not required/reference】
                - HC3is normalNot required
                  （WLSwithetcVarianceizationcompletemi）
                - thoughtforSensitivityanalysisandandcombinekipossible
                ```

                ---

                ### ⚠️ WLSofdangerdangersignificance 

                #### Small samplewithWLStheuseuand:

                1. **weightestimationErrorisLarge**
                   - trueofVariancethecorrectconfirmtoEstimationcannot
                   - erroredweightwithEstimationisdistortmu

                2. **oversurplussuitablego（Overfitting）**
                   - observedValuetooffuiteinguisstrongsugiru
                   - generalizedsignificance abilityislowLower

                3. **extremeedgeaweight**
                   - Outliertoextremeedgeaweightisattachku
                   - countValuetekiunstablesignificance 

                4. **interpretationoftroubledifficultsa**
                   - weightwithMeaneffectofmeaningisnotclearconfirm
                   - ReproducibilityisLow

                ---

                ### 🔬 This appinrealinstall

                **directionmethod**: Residualbe-suofinverseVarianceweight

                ```python
                1. normalofOLSwithResidualthecalculation
                2. weight = 1 / (Residual²)
                3. weighttheNormalization（Mean=1）
                4. WLSwithagainEstimation
                ```

                **⚠️ thisdirectionmethodofquestiontopic**:
                - Residualbe-suaofwithtrueofVarianceandis Different
                - Small samplewithis Residualisunstable
                - **recommended**: matterBeforeinforeport（skilltechniquetekimultiplemakecount、alreadyknowofVariance）fromweightthedecideset

                ---

                ### 💡 recommendedmatterterm（maandme）

                #### Small sample（n < 15/group）:
                1. ❌ **WLSis usewanot**
                2. ✅ **WelchType + HC3themainanalysisto**
                3. ✅ **Permutation testwithImputation**

                #### largeSample（n ≥ 15/group）:
                1. ✅ **Variancedifferencevisualization,Test**
                2. ✅ **weightexplicitlyandWLS**
                3. ⚠️ **HC3and comparison（Sensitivityanalysis）**

                #### always:
                - 📊 **diagnosisplotrequired**（Residualpattern-nConfirm）
                - 📝 **weightofrootbasistheclearki**（theorytextreporttelltime）
                - 🔄 **multiplecounthandmethodwithResults robustsignificance Confirm**

                ---

                ### 📚 reference: differVariancetoofagainstprocessmethodofcomparison

                | directionmethod | Samplesizeneedcase | robustsignificance  | recommendeddegree（smalln） | recommendeddegree（largen） |
                |------|------------------|--------|--------------|--------------|
                | **Welch tTest** | smallalsoOK | High | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ |
                | **HC3robustSE** | smallalsoOK | High | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ |
                | **Permutation test** | smallalsoOK | non-alwaysHigh | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
                | **WLS** | largeisrequired | Low（smallntime） | ⭐ | ⭐⭐⭐⭐⭐ |
                | **normalofOLS** | - | none | ❌ | ⭐ |

                **resulttheory**: aaedofsettotalisSmall sampleif、WLSis avoidketeWelch+HC3+Permutation testtheuseimashiyou。
                """)


    if "LMM" in analysis_method:
        st.markdown("#### LMM Option")
        reml_method = st.radio(
            "Estimationdirectionmethod:",
            ["REML (recommended)", "ML"],
            index=0,
            help="REML: unbiased VarianceEstimation。ML: Modelfor comparison（AIC/BIC）"
        )

    st.markdown("---")
    st.markdown("### WelchTypetTest")
    run_welch = st.checkbox("WelchTypetTesttheRun", value=False,
                           help="**Main variableis2groupcase ofonlyValid**\n\n"
                                "etcVariancetheassumptionshinot2groupcomparison（Covariateadjustmentnone）\n"
                                "simplea2groupcomparisonandandreferenceValuethecalculateout")

    if run_welch:
        st.info("ℹ️ **WelchTypetTestabout:**\n\n"
                "- Main variableof2groupBetweenwithstraightconnectcomparison（Covariateadjustmentnone）\n"
                "- etcVariancesignificance theassumptionshinotrobustaTest\n"
                "- Main variableis3groupfollowingUppercase ofis Rundonenot")

    st.markdown("---")
    st.markdown("### Permutation test (Permutation Test)")
    run_permutation = st.checkbox("Permutation testtheRun", value=False,
                                  help="**Testtarget: Main variable1only**\n\n"
                                       "Small sampleordistributionassumptionisquestionablecasetorecommended\n"
                                       "Main variable2orsoofotherofCovariateis controlvariablecountandandhandlewareing")

    if run_permutation:
        st.info("ℹ️ **Permutation testoftargetvariablecount:** Main variable1onlyisTestdoneing。\n\n"
                "Main variable2orsoofother variablescountis Covariateandandcontroldoneing。")

        n_permutations = st.slider("placechangetimescount:",
                                    min_value=1000, max_value=50000,
                                    value=10000, step=1000,
                                    help="ManyaboutcorrectconfirmdaistimeBetweeniskakariing")
        perm_method = st.radio(
            "placechangebattleomit:",
            ["Freedman-Lanemethod (recommended)", "Simple label permutation"],
            index=0,
            help="Freedman-Lane: Covariatestrictly controlling Test\n"
                 "Main variable1ofeffecttheTestshi、other variablescountis controldoneing"
        )

        # one-sidedTestOption
        perm_sided = st.radio(
            "Testofdirectiondirection:",
            ["two-sidedTest (two-sided)", "one-sidedTest: Main variable1 > Criterion (greater)", "one-sidedTest: Main variable1 < Criterion (less)"],
            index=0,
            help="**two-sidedTest（recommended）**: effectofdirectiondirectionthequestionwazu、differenceisexistkatheTest\n\n"
                 "**one-sidedTest (greater)**: matterBeforeto「Main variable1isCriterionfromLarge」andadvancethinkdoingcaseonlyuse\n"
                 "e.g.: Female > Male theHypothesisanddocase\n\n"
                 "**one-sidedTest (less)**: matterBeforeto「Main variable1isCriterionfromSmall」andadvancethinkdoingcaseonlyuse\n\n"
                 "⚠️ one-sidedTestis matterBeforeHypothesisisexistcaseonlyuseplease do"
        )

    st.markdown("---")
    st.markdown("### VisualizationOption")
    show_diagnostics = st.checkbox("diagnosisplotdisplayed", value=True,
                                   help="Residualplot、Q-Qplotetc")
    show_emm = st.checkbox("Estimated marginal means (EMM) thecalculation", value=True,
                          help="groupby AdjustedMeanValuethecalculateout")

# =============================================================================
# File Upload
# =============================================================================

st.markdown("## 1️⃣ data upload")

# DataInputdirectionmethodofSelect
input_method = st.radio(
    "DataInputdirectionmethodselect:",
    ("PCA scoresFile + metaDatamanuactionInput", "PCA scoresUpload file（metaDatawith）"),
    help="pca.pyfromDownloaddonePCA scoresFileusedocase、metaDatamanuactionwithInputcaning。"
)

if input_method == "PCA scoresFile + metaDatamanuactionInput":
    # metaDatamanuactionInput
    st.markdown("""
    **PCA scoresFileandmetaDatatheby々toInput**

    1. pca.pyfromDownloaddonePCA scoresFile（SampleNameandPC1, PC2, ...including ）theupload
    2. Loweroffuo-muwithmetaDatatheInput
    """)

    uploaded_file = st.file_uploader("PCA scoresFileselect", type=["tsv", "csv", "txt"],
                                     key="pca_scores_file")

else:  # metaDatawithUpload file
    st.markdown("""
    **RequiredaDatashapeformula**: TSV or  CSV File

    - **1Row = 1Sample** (pseudobulk reberu)
    - **PC score Column**: PC1, PC2, PC3, ... etc
    - **Categorical variable**: sex, subtype, condition etc
    - **arbitrary: Donor/Subject ID** (Iterationmeasureisexistcase)
    - **arbitrary: countValueCovariate**: age, batch, depth etc
    """)

    uploaded_file = st.file_uploader("Fileselect", type=["tsv", "csv", "txt"])

if uploaded_file is not None:
    # Filereadloading
    try:
        df = pd.read_csv(uploaded_file, sep="\t", index_col=0)
    except:
        uploaded_file.seek(0)
        df = pd.read_csv(uploaded_file, sep=",", index_col=0)

    # indextheResetandnormalofColumntodo
    df = df.reset_index()

    # indexColumnNamethesuitablecuttoSettings
    if df.columns[0] == 'index':
        df.rename(columns={'index': 'sample_id'}, inplace=True)

    st.success(f"✅ FilereadloadingComplete: {df.shape[0]} Row × {df.shape[1]} Column")

    with st.expander("📋 Datapurebiyu- (first10Row)", expanded=True):
        st.dataframe(df.head(10))

    # session_stateofInitialization
    if 'metadata_submitted' not in st.session_state:
        st.session_state.metadata_submitted = False
    if 'metadata_df' not in st.session_state:
        st.session_state.metadata_df = None
    if 'combined_df' not in st.session_state:
        st.session_state.combined_df = None

    # metaDatahandactionInputcase of
    if input_method == "PCA scoresFile + metaDatamanuactionInput":
        st.markdown("## 2️⃣ metadata Input")
        st.markdown("Sampleper tometaDatathePlease enter")

        # SampleNameobtain（firstColumnor RowName）
        if df.index.name or (df.index[0] != 0):  # indexisSampleNamecase of
            sample_names = df.index.tolist()
        else:  # firstColumnisSampleNamecase of
            sample_names = df.iloc[:, 0].tolist()

        # metaDataColumnNameofInput
        st.markdown("### metaDataColumn Settings")
        col1, col2 = st.columns(2)

        with col1:
            meta_cols_input = st.text_area(
                "metaDataColumnName（kanmaareacutrior reformRowareacutri）:",
                value="sex, subtype",
                help="e.g.: sex, subtype, donor, age, batch"
            )

        # metaDataColumnNameofparse
        if ',' in meta_cols_input:
            meta_col_names = [x.strip() for x in meta_cols_input.split(',') if x.strip()]
        else:
            meta_col_names = [x.strip() for x in meta_cols_input.split('\n') if x.strip()]

        with col2:
            st.markdown("**InputdoColumn:**")
            for col in meta_col_names:
                st.write(f"- {col}")

        # metaDataInputfuo-mu
        st.markdown("### metadata Input")

        with st.form("metadata_input"):
            # eachSampleofmetaDatatheDataedeitawithInput
            meta_df = pd.DataFrame(index=sample_names, columns=meta_col_names)

            # DESeq2-LRTandSameAlgorithmwithdefaultValuetheSettings
            # 1. endtailoftogetherthroughneedelementexcludingremove
            sample_names_str = [str(s) for s in sample_names]
            group_base = remove_common_suffix(sample_names_str)
            # 2. endtailofcountcharacterexcludingremoveandgroup-puNametheEstimation
            group_names = [remove_sample_num(s) for s in group_base]

            # eachmetaDataColumn defaultValuetheSettings
            for col in meta_col_names:
                col_lower = col.lower()

                if col_lower in ['sex', 'gender']:
                    # group-puNamefromsignificance bytheinfermeasure（male/femaleisincludemarerucase）
                    def infer_sex(name):
                        name_lower = name.lower()
                        if 'male' in name_lower and 'female' not in name_lower:
                            return 'Male'
                        elif 'female' in name_lower or 'fem' in name_lower:
                            return 'Female'
                        else:
                            return 'Male'  # default
                    meta_df[col] = [infer_sex(g) for g in group_names]

                elif col_lower in ['donor', 'subject', 'patient', 'individual']:
                    # SampleNamefromDonorIDtheinfermeasure
                    # pattern-n1: sample_D01_A → D01
                    # pattern-n2: D01_typeA → D01
                    # pattern-n3: continuousnumberwithattachgive
                    donor_ids = []
                    for s in sample_names_str:
                        parts = s.split('_')
                        # D01ofyouapattern-nthesearchsu
                        donor_found = False
                        for part in parts:
                            if part.startswith('D') and len(part) >= 2 and part[1:].isdigit():
                                donor_ids.append(part)
                                donor_found = True
                                break
                        if not donor_found:
                            # group-puNameuse
                            donor_ids.append(group_names[len(donor_ids)])
                    meta_df[col] = donor_ids

                elif col_lower in ['subtype', 'celltype', 'cell_type', 'type', 'tissue']:
                    # group-puNamethesoofmamause（thisissabutypethetablesuthatisMany）
                    meta_df[col] = group_names

                elif col_lower in ['batch', 'library', 'plate', 'run']:
                    # batchNumbertheinfermeasure
                    # pattern-n1: sample_B1_typeA → B1
                    # pattern-n2: allSamplewithSamecaseis  Batch1
                    batch_ids = []
                    for s in sample_names_str:
                        parts = s.split('_')
                        batch_found = False
                        for part in parts:
                            if part.startswith('B') and len(part) >= 2 and part[1:].isdigit():
                                batch_ids.append(part)
                                batch_found = True
                                break
                            elif 'batch' in part.lower():
                                batch_ids.append(part)
                                batch_found = True
                                break
                        if not batch_found:
                            batch_ids.append('Batch1')
                    meta_df[col] = batch_ids

                elif col_lower in ['age', 'depth', 'weight', 'height', 'score']:
                    # countValueTypecase ofis 0withInitialization
                    meta_df[col] = 0

                elif col_lower in ['condition', 'treatment', 'group', 'status']:
                    # group-puNameuse
                    meta_df[col] = group_names

                else:
                    # soofotheris group-puNameuse
                    meta_df[col] = group_names

            st.write("eachSampleofmetaDatatheEditplease do:")
            st.info("💡 SampleNamefromselfactioninfermeasuredoneeddefaultValueisInputdoneteiing。RequiredtorespondjiteEditplease do。")
            edited_meta_df = st.data_editor(meta_df, use_container_width=True)

            submitted = st.form_submit_button("metaDatatheconfirmset", type="primary")

        if submitted:
            # session_statetoSave
            st.session_state.metadata_submitted = True
            st.session_state.metadata_df = edited_meta_df.copy()
            st.session_state.meta_col_names = meta_col_names

            # PCA scoresandmetaDatatheresultgo
            if df.index.name or (df.index[0] != 0):
                df_combined = pd.concat([df, edited_meta_df], axis=1)
            else:
                # firstColumnisSampleNamecase of
                df_temp = df.set_index(df.columns[0])
                df_combined = pd.concat([df_temp, edited_meta_df], axis=1)
                df_combined.reset_index(inplace=True)
                df_combined.rename(columns={'index': df.columns[0]}, inplace=True)

            st.session_state.combined_df = df_combined.copy()

        # metaDataisconfirmsetdoneteirucase、resultgodoneDatafure-muuse
        if st.session_state.metadata_submitted and st.session_state.combined_df is not None:
            df = st.session_state.combined_df.copy()
            st.success("✅ metaDataInputComplete")

            st.write("metaDataConfirm:")
            for col in st.session_state.meta_col_names:
                st.write(f"**{col}**: " + ', '.join(df[col].unique().astype(str)))

    # =============================================================================
    # Column Mapping
    # =============================================================================

    step_number = "3️⃣" if input_method == "PCA scoresFile + metaDatamanuactionInput" else "2️⃣"
    st.markdown(f"## {step_number} Column Mapping")
    st.markdown("eachvariablecountsupporteddoColumnthePlease select")

    with st.form("column_mapping"):
        col1, col2, col3 = st.columns(3)

        with col1:
            st.markdown("### 🔹 requiredColumn")

            # Sample ID
            sample_col = st.selectbox("SampleIDColumn:", df.columns,
                                       index=0 if 'sample' in df.columns[0].lower() else 0)

            # PC score columns
            numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
            pc_cols = [col for col in numeric_cols if col.upper().startswith('PC')]

            if len(pc_cols) == 0:
                st.warning("'PC'withstartmaruColumnisseetsukarinot。AllofcountValueColumndisplayeddo。")
                pc_cols = numeric_cols

            pc_col = st.selectbox("analysisdoPC scoreColumn:", pc_cols,
                                  index=0,
                                  help="e.g.: PC3, PC4 etc。Multiple PCtheanalysisdocaseis 、1tsuzutsuPlease run")

        with col2:
            st.markdown("### 🔹 Main variable")

            # categoryColumnobtain
            cat_cols = df.select_dtypes(include=['object', 'category']).columns.tolist()

            # SampleIDColumnexcludingOutside（yuni-kuValueismanysugiruColumnis excludeOutside）
            cat_cols_filtered = []
            for col in cat_cols:
                n_unique = df[col].nunique()
                n_samples = len(df)
                unique_ratio = n_unique / n_samples

                # FilteringCriterion:
                # 1. SampleID（Allonemean: ratio >= 0.9）excludingOutside
                # 2. yuni-kuValueismanysugiru（>50）ColumnexcludingOutside
                # 3. however、yuni-kuValueis2followingUpperexistthat（setcountColumnexcludingOutside）
                if unique_ratio < 0.9 and n_unique >= 2 and n_unique <= 50:
                    cat_cols_filtered.append(col)

            if not cat_cols_filtered:
                st.warning("⚠️ suitablecutaCategorical variableisseetsukarinot。AllofcategoryColumndisplayeddo。")
                cat_cols_filtered = cat_cols

            # Main variable (e.g.: sex)
            main_var = st.selectbox("🔴 Main variable 1 (required):",
                                    cat_cols_filtered,
                                    index=0 if len(cat_cols_filtered) > 0 else None,
                                    help="Testdoneimainneedavariablecount（e.g.: sex, treatment etc）。\n"
                                         "thisvariablecountis ModelofMain effectandandnecessaryzuincludemareing。\n"
                                         "SampleIDofyouaonemeanofColumnis selfactionallyexcludeOutsidedoneing。")

            # Blocking variable (e.g.: subtype)
            other_cats = [col for col in cat_cols_filtered if col != main_var]
            blocking_var = st.selectbox("🔴 Main variable 2 / Blocking variable (arbitrary):",
                                        ["none"] + other_cats,
                                        index=1 if len(other_cats) > 0 else 0,
                                        help="**important:** thisvariablecountmoModelofMain effectandandincludemareing。\n\n"
                                             "e.g.: sex theMain variable1、subtype theMain variable2andandSelectdoand、\n"
                                             "Modelformulais  `PC ~ C(sex) + C(subtype)` andariing。\n\n"
                                             "**usee.g.:**\n"
                                             "- 2ofMain variablethesametimetoTestdoneicase（sex and subtype etc）\n"
                                             "- burokingudesign（donor, batch etc）withStratificationizationdoneicase\n\n"
                                             "InteractiontermtheAdddoand、2ofvariablecountofmutualmakeusemoTestcaning。")
            blocking_var = None if blocking_var == "none" else blocking_var

        with col3:
            st.markdown("### 🔹 Hierarchical structure")

            # Donor ID for random effects
            donor_col = st.selectbox("Donor/Subject ID (Random effectuse):",
                                     ["none", "SampleIDandSame"] + df.columns.tolist(),
                                     index=0,
                                     help="Iterationmeasureisexistcasetopointset。LMMwith (1|donor) ofRandom effectandanduse")

            if donor_col == "none":
                donor_col = None
            elif donor_col == "SampleIDandSame":
                donor_col = sample_col

            # Additional random effect
            additional_re = st.selectbox("AddofRandom effect (e.g.: batch, library):",
                                         ["none"] + df.columns.tolist(),
                                         index=0,
                                         help="skilltechniquetekiabatch effectetc Random effectandandhandlecasetouse")
            additional_re = None if additional_re == "none" else additional_re

        # Additional covariates
        st.markdown("### 🔹 AddofCovariate（Fixed effect）")

        # alreadytousedoneColumnexcludingOutside
        used_cols = [sample_col, pc_col, main_var, blocking_var, donor_col, additional_re]
        remaining_cols = [col for col in df.columns if col not in used_cols and col is not None]

        additional_covars = st.multiselect(
            "AddofCovariateselect (countValueor category):",
            remaining_cols,
            help="age、batch、depthetc. adjustmentdoneivariablecount。thisrais Fixed effectandandModeltoincludemareing"
        )

        # Interaction terms
        st.markdown("### 🔹 Interactionterm")
        st.markdown("_Interactionis 、2ofvariablecountofeffectisindependentwithnotcasetousedo（e.g.: significance differenceissabutypebyDifferent）_")

        potential_interactions = []
        if blocking_var:
            potential_interactions.append(f"{main_var} : {blocking_var}")

        for covar in additional_covars:
            if df[covar].dtype in ['object', 'category']:
                potential_interactions.append(f"{main_var} : {covar}")

        interactions = st.multiselect(
            "Interactiontermselect (arbitrary):",
            potential_interactions,
            help="e.g.: sex:subtype is 、significance differenceissabutypeBetweenwithDifferentkathetested。Dataistenparttoexistcaseonlyuserecommended"
        )

        # fuo-muofsendtrustbutton
        mapping_submitted = st.form_submit_button("Column Mappingtheconfirmset", type="primary")

    # =============================================================================
    # Run Analysis Button
    # =============================================================================

    # Column Mappingisconfirmsetdoneedcaseonlyanalysisbuttondisplayed
    if mapping_submitted or 'mapping_confirmed' in st.session_state:
        if mapping_submitted:
            # ColumnMappingisvariablefurtherdoneedcase、followingBeforeofanalysisresulttheallclear
            st.session_state.mapping_confirmed = True
            st.session_state.sample_col = sample_col
            st.session_state.pc_col = pc_col
            st.session_state.main_var = main_var
            st.session_state.blocking_var = blocking_var
            st.session_state.donor_col = donor_col
            st.session_state.additional_re = additional_re
            st.session_state.additional_covars = additional_covars
            st.session_state.interactions = interactions

            # followingBeforeofanalysisresulttheclear
            results_to_clear = [
                'ols_model', 'ols_coef', 'ols_anova', 'ols_emm',
                'lmm_result', 'lmm_coef',
                'perm_result',
                'fig_forest_ols', 'fig_diagnostic', 'fig_emm',
                'fig_forest_lmm', 'fig_null_dist'
            ]
            for key in results_to_clear:
                if key in st.session_state:
                    del st.session_state[key]

        st.success("✅ Column MappingComplete")

        # session_statefromValueobtain
        sample_col = st.session_state.sample_col
        pc_col = st.session_state.pc_col
        main_var = st.session_state.main_var
        blocking_var = st.session_state.blocking_var
        donor_col = st.session_state.donor_col
        additional_re = st.session_state.additional_re
        additional_covars = st.session_state.additional_covars
        interactions = st.session_state.interactions

        if st.button("🚀 analysistheRun", type="primary"):

            step_number_check = "4️⃣" if input_method == "PCA scoresFile + metaDatamanuactionInput" else "3️⃣"
            st.markdown(f"## {step_number_check} DataQuality checks")

            # analysisuseDatafure-muofPreparation
            analysis_cols = [sample_col, pc_col, main_var]
            if blocking_var:
                analysis_cols.append(blocking_var)
            if donor_col and donor_col != sample_col:
                analysis_cols.append(donor_col)
            if additional_re:
                analysis_cols.append(additional_re)

            # additional_covarsisemptywithnotcaseonlyAdd
            if additional_covars:
                # risutowithexistthatconfirmandfromAdd
                if isinstance(additional_covars, list):
                    analysis_cols.extend(additional_covars)
                else:
                    analysis_cols.append(additional_covars)

            # NoneexcludingOutsideshi、weightmultiplemoexcluderemove
            analysis_cols = [col for col in analysis_cols if col is not None]
            analysis_cols = list(dict.fromkeys(analysis_cols))  # weightmultiplethekeepholdorderwithexcluderemove

            # debaguinforeport
            st.write("**SelectdoneedColumn:**", analysis_cols)

            analysis_df = df[analysis_cols].copy()

            # Missing valueofexcluderemove
            n_before = len(analysis_df)
            analysis_df = analysis_df.dropna()
            n_after = len(analysis_df)

            if n_before > n_after:
                st.warning(f"⚠️ Missing valueincluding  {n_before - n_after} RowtheDeletediddone")

            # Categorical variableofmarklevelization
            for col in analysis_df.columns:
                # colischaracterColumnwithexistthatconfirm
                if isinstance(col, str):
                    try:
                        if analysis_df[col].dtype == 'object':
                            analysis_df[col] = analysis_df[col].astype('category')
                    except (KeyError, AttributeError):
                        # Columnissonexistshinot、or Typeisnotcaseis sukipu
                        continue

            # SamplesizeofConfirm（lackmeasureinforeportincluding ）
            st.markdown("### 📊 Samplesizeofneedabout ")

            # lackmeasurecheck
            original_n = len(df)  # originofDatacount
            analysis_n = len(analysis_df)  # analysistousedoDatacount
            excluded_n = original_n - analysis_n

            if excluded_n > 0:
                st.warning(f"⚠️ **lackmeasureValueexcludeOutside**: {excluded_n} caseofSampleislackmeasureValuetofromexcludeOutsidedonemadone（originData: {original_n} → analysisuse: {analysis_n}）")

            col1, col2 = st.columns(2)

            with col1:
                st.metric("analysistousedoneSamplecount", len(analysis_df))

                # Main variableofdistribution
                st.write(f"**{main_var} ofdistribution:**")
                main_var_counts = analysis_df[main_var].value_counts()
                st.dataframe(main_var_counts)

                # notaveragebalanceofWarning
                if len(main_var_counts) > 1:
                    max_count = main_var_counts.max()
                    min_count = main_var_counts.min()
                    if max_count / min_count > 3:
                        st.warning(f"⚠️ **groupofnotaveragebalance**: sailargegroup({max_count})andsaismallgroup({min_count})with3timesfollowingUpperofdifferenceispresenting")

            with col2:
                if blocking_var:
                    st.write(f"**cross tabulationtable: {main_var} × {blocking_var}**")
                    cross_tab = pd.crosstab(analysis_df[main_var], analysis_df[blocking_var])
                    st.dataframe(cross_tab)

                    # Small samplesizeofWarning
                    min_cell_size = cross_tab.min().min()
                    if min_cell_size < 3:
                        st.warning(f"⚠️ onepartofcellofSamplesizeisnon-alwaysSmall（saismall: {min_cell_size}）is。EstimationisunstabletoaruPossiblesignificance ispresenting。")

                if donor_col and donor_col != sample_col:
                    n_donors = analysis_df[donor_col].nunique()
                    st.metric("Donorcount", n_donors)

                    if n_donors < 5:
                        st.warning("⚠️ Donorcountis5less thanis。Random effectofVarianceEstimationisunstabletoaruPossiblesignificance ispresenting。")

            # collinearity check
            st.markdown("### 🔍 collinearity diagnosis")

            formula_terms = [main_var]
            if blocking_var:
                formula_terms.append(blocking_var)
            formula_terms.extend(additional_covars)

            collinearity_issues = detect_collinearity(analysis_df, formula_terms)

            if collinearity_issues:
                st.error("**completecollinearity isdetectiondonemadone:**")
                for issue in collinearity_issues:
                    st.markdown(issue)
                st.warning("""
                ⚠️ **completetotogetherroledoingvariablecountis 、effectthepartseparateandEstimationcannot。**

                e.g.: existsubtypeisAllFemalecase of、sexeffectandsubtypeeffectis recognizebynotabilityis。

                **againstprocessmethod:**
                - anykaofvariablecounttheDelete
                - fromaveragebalanceoftakeedDatathecollectgather
                - Interactiontermusesezu、kidescribetekiacomparisontostaymeru
                """)
            else:
                st.success("✅ completecollinearity is detectiondonenotwithdone")

            # VIF for numeric variables
            numeric_covars = [col for col in additional_covars
                             if pd.api.types.is_numeric_dtype(analysis_df[col])]

            if len(numeric_covars) > 1:
                st.markdown("**countValueCovariateofmultiplecollinearity  (VIF):**")
                st.caption("VIF > 10: strongmultiplecollinearity present。VIF > 5: middleapproximatelyofmultiplecollinearity present。")

                try:
                    vif_df = calculate_vif(analysis_df, numeric_covars)
                    vif_df_display = vif_df.copy(); vif_df_display['VIF'] = vif_df_display['VIF'].round(2); st.dataframe(vif_df_display, use_container_width=True)

                    if (vif_df['VIF'] > 10).any():
                        st.warning("⚠️ VIF > 10 ofvariablecountispresenting。strongmultiplecollinearity isshowsuggestdoneing。")
                    elif (vif_df['VIF'] > 5).any():
                        st.info("ℹ️ VIF > 5 ofvariablecountispresenting。middleapproximatelyofmultiplecollinearity isshowsuggestdoneing。")
                except Exception as e:
                    st.info(f"VIFcalculationthesukipudiddone: {str(e)}")

            # =============================================================================
            # Build Formula
            # =============================================================================

            st.markdown("## 4️⃣ Modelformulaofstructurebuild")

            # debaguinforeport
            with st.expander("🔍 Modelstructurebuildinforeport（SelectdoneedvariablecountofConfirm）", expanded=True):
                st.write(f"**Main variable 1 (main_var):** {main_var}")
                st.write(f"**Main variable 2 / Blocking variable (blocking_var):** {blocking_var}")
                st.write(f"**AddCovariate (additional_covars):** {additional_covars}")
                st.write(f"**Interactionterm (interactions):** {interactions}")
                st.info("💡 Modelformulatoincludemeedivariablecountis 、Upperof 'Main variable 1' or  'Main variable 2' withPlease select。\n\n"
                        "'AddofCovariate' is adjustmentvariablecountandandusedone、mainaTesttargetwithis presentnot。")

            # Fixed effectofformula
            fixed_terms = [f"C({main_var})"]
            if blocking_var:
                fixed_terms.append(f"C({blocking_var})")

            for col in additional_covars:
                if analysis_df[col].dtype in ['object', 'category']:
                    fixed_terms.append(f"C({col})")
                else:
                    fixed_terms.append(col)

            # InteractiontermofAdd
            if interactions:
                for interaction in interactions:
                    vars_in_interaction = [v.strip() for v in interaction.split(":")]
                    interaction_term = " * ".join([
                        f"C({v})" if analysis_df[v].dtype in ['object', 'category'] else v
                        for v in vars_in_interaction
                    ])
                    # alreadytoincludemareteinotcaseonlyAdd
                    if interaction_term not in fixed_terms:
                        fixed_terms.append(interaction_term)

            formula = f"{pc_col} ~ " + " + ".join(fixed_terms)

            # Modelcontained invariablecountexplicitlyallyDisplay
            st.markdown("**📊 Modelcontained invariablecount:**")
            col_a, col_b = st.columns(2)
            with col_a:
                st.success(f"✅ Main variable 1: **{main_var}**")
                if blocking_var:
                    st.success(f"✅ Main variable 2: **{blocking_var}**")
                else:
                    st.warning("⚠️ Main variable 2: Selectnone")
            with col_b:
                if additional_covars:
                    st.info(f"📝 adjustmentvariablecount: {', '.join(additional_covars)}")
                if interactions:
                    st.info(f"🔗 Interaction: {', '.join(interactions)}")

            st.markdown("**Modelformula:**")
            st.code(formula, language="r")

            if donor_col:
                st.markdown(f"**Random effect (LMMonly):** `(1|{donor_col})`")

            # =============================================================================
            # Run OLS
            # =============================================================================

            if "OLS" in analysis_method:
                st.markdown("---")
                st.markdown("## 5️⃣ OLS analysisresult")

                try:
                    # WLS vs OLS
                    if use_wls:
                        # First fit OLS to get residuals for weights
                        ols_temp = smf.ols(formula, data=analysis_df).fit()
                        # Calculate weights as inverse of squared residuals
                        weights = 1.0 / (ols_temp.resid ** 2)
                        # Normalize weights
                        weights = weights / weights.mean()
                        # Fit WLS
                        ols_model = smf.wls(formula, data=analysis_df, weights=weights).fit()
                        st.info(f"ℹ️ WLS（addweightsaismalltworidemethod）usedoinging。weightis ResidualofinverseVariancewithcalculationdonemadone")
                    else:
                        # Modeloffuito (OLS)
                        ols_model = smf.ols(formula, data=analysis_df).fit()

                    # suitablecutaCovarianceRowColumnobtain
                    if se_type == "HC3 (robustEstimation,recommended)":
                        ols_model = ols_model.get_robustcov_results(cov_type='HC3')
                        st.info("ℹ️ HC3robustStandard errorusedoinging（heterogeneous Variancesupported）")
                    elif se_type == "Cluster-robust (cluster-robust)":
                        cluster_var = donor_col if donor_col else additional_re
                        if cluster_var:
                            ols_model = ols_model.get_robustcov_results(
                                cov_type='cluster',
                                groups=analysis_df[cluster_var]
                            )
                            st.info(f"ℹ️ cluster-robustStandard errorusedoinging（cluster-variablecount: {cluster_var}）")
                        else:
                            st.warning("cluster-variablecountispointsetdoneteinot。normalofHC3usedo。")
                            ols_model = ols_model.get_robustcov_results(cov_type='HC3')

                    # Modelformulaandreferreferencelevelofclearshow
                    st.markdown("### 📋 Coefficienttable")
                    st.caption(f"**Modelformula**: `{formula}`")

                    # referreferencelevelofextractandDisplay
                    reference_levels = {}
                    for col in analysis_df.columns:
                        if analysis_df[col].dtype == 'object' or analysis_df[col].dtype.name == 'category':
                            if col in formula:
                                ref_level = analysis_df[col].iloc[0] if len(analysis_df) > 0 else "N/A"
                                # Patsyis arufuabetoorderofsaifirsttheCriteriontodo
                                unique_vals = sorted(analysis_df[col].unique())
                                if len(unique_vals) > 0:
                                    reference_levels[col] = unique_vals[0]

                    if reference_levels:
                        ref_info = ", ".join([f"`{var}` ofCriterion={ref}" for var, ref in reference_levels.items()])
                        st.caption(f"**referreferencelevel**: {ref_info}")
                    conf_int = ols_model.conf_int()

                    # conf_intisDataFramekandarraykathejudgeset
                    if isinstance(conf_int, pd.DataFrame):
                        ci_lower = conf_int.iloc[:, 0]
                        ci_upper = conf_int.iloc[:, 1]
                    else:
                        # numpydistributeColumn case
                        ci_lower = conf_int[:, 0]
                        ci_upper = conf_int[:, 1]

                    coef_table = pd.DataFrame({
                        'Coef': ols_model.params,
                        'Std_Error': ols_model.bse,
                        't_value': ols_model.tvalues,
                        'p_value': ols_model.pvalues,
                        'CI_lower': ci_lower,
                        'CI_upper': ci_upper
                    })

                    # Significantsignificance ofkinumbertheAdd
                    def sig_symbol(p):
                        if p < 0.001:
                            return '***'
                        elif p < 0.01:
                            return '**'
                        elif p < 0.05:
                            return '*'
                        elif p < 0.1:
                            return '.'
                        else:
                            return ''

                    coef_table['Sig'] = coef_table['p_value'].apply(sig_symbol)

                    # countValueroundingDisplayuseofDataFramecreate
                    display_table = coef_table.copy()
                    display_table['Coef'] = display_table['Coef'].round(4)
                    display_table['Std_Error'] = display_table['Std_Error'].round(4)
                    display_table['t_value'] = display_table['t_value'].round(3)
                    display_table['p_value'] = display_table['p_value'].apply(lambda x: f"{x:.4g}")
                    display_table['CI_lower'] = display_table['CI_lower'].round(4)
                    display_table['CI_upper'] = display_table['CI_upper'].round(4)

                    st.dataframe(display_table, use_container_width=True)

                    st.caption("Significantlevel: *** p<0.001, ** p<0.01, * p<0.05, . p<0.1")

                    # Coefficienttableofinterpretationguide
                    with st.expander("📘 Coefficienttableofseedirection（p_valueColumn meaning）"):
                        st.markdown("""
                        ### p_valueColumn meaning

                        **eachCoefficientofSignificantsignificance theshowdo（other variablescounttheadjustmentdoneUpperwith）**

                        - ✅ **Covariateadjustmentpresent**: Allofvariablecountthesametimetothinkconsiderdoneresult
                        - koof p-valueis 、other variablescountofshadoweffectthetakeriexcludeied「purepureaeffect」theTestdoinging

                        ---

                        ### tablekiofmeaning
                        """)

                        # realoccasionofModelcontained invariablecountfromDescriptiongenerate
                        example_vars = []
                        for idx in display_table.index:
                            idx_str = str(idx)  # Ensure idx is a string
                            if idx_str != 'Intercept' and 'C(' in idx_str:
                                example_vars.append(idx_str)

                        if example_vars:
                            st.markdown("**thisModelexample:**")
                            # first1-2pieceofvariablecountthee.g.andanduse
                            for var_name in example_vars[:2]:
                                # C(sex)[T.M] ofyouashapeformulafrompartunderstand
                                if '[T.' in var_name:
                                    var_part = var_name.split('[')[0]  # C(sex)
                                    level_part = var_name.split('[T.')[1].rstrip(']')  # M
                                    base_var = var_part.replace('C(', '').replace(')', '')  # sex

                                    st.markdown(f"""
                                    - **`{var_name}`**
                                      - `C({base_var})`: {base_var}thecategorykaruvariablecountandandhandle
                                      - `[T.{level_part}]`: Treatment coding with {level_part}levelthetablesu
                                      - **meaning**: 「{level_part}is Criterioncategoryandratiobetewhichonlydifferuka」
                                      - **p-valueofinterpretation**: {level_part}andCriterioncategorytoSignificantadifferenceisexistka（**other variablescounttheAdjusted**）
                                    """)

                        st.markdown("""
                        ---

                        ### generaltekiatableki

                        | tableki | meaning |
                        |-----|------|
                        | `C(variablecountName)` | categorykaruvariablecountandandhandle |
                        | `[T.levelName]` | Treatment codingwithpointsetlevelthetablesu |
                        | `C(variablecountName)[T.levelName]` | 「thislevelis Criterioncategoryandratiobetewhichonlydifferuka」 |

                        ### Criterioncategory（Reference level）

                        - default: **arufuabetoorderwithfirstValue**
                        - CriterioncategoryofRowis CoefficienttabletoDisplaydonenot（Coefficient=0andandhandlewareru）
                        - otheroflevelis 「Criterionand difference」andandinterpretationdoneing

                        ### Welch tTestand differi

                        | Testdirectionmethod | Covariateadjustment | Displayplaceplace |
                        |---------|----------|---------|
                        | Welch tTest | ❌ none（simple2groupcomparison） | inforeportbokusu |
                        | Coefficienttableofp_value | ✅ present（othervariablecountconsidering ） | thistable |

                        **recommended**: correctconfirmastatisticalanalysistois 、thisCoefficienttableofp-valueuseplease do。
                        """)

                    # ModelGoodness of fit
                    st.markdown("### 📈 ModelGoodness of fit")
                    col1, col2, col3, col4 = st.columns(4)
                    with col1:
                        st.metric("R²", f"{ols_model.rsquared:.4f}")
                    with col2:
                        st.metric("Adjusted R²", f"{ols_model.rsquared_adj:.4f}")
                    with col3:
                        st.metric("AIC", f"{ols_model.aic:.2f}")
                    with col4:
                        st.metric("BIC", f"{ols_model.bic:.2f}")

                    # ANOVA
                    try:
                        if anova_type == "Type II (recommended)":
                            st.markdown("### 📊 Type II ANOVA")
                            st.caption("each variable OveralltekiaSignificantsignificance theTest（other variablescountincluding Modelincontributegive）")
                            anova_table = anova_lm(ols_model, typ=2)
                        else:  # Type III
                            st.markdown("### 📊 Type III ANOVA")
                            st.caption("AllofInteractiontheincludemeedeach variable Significantsignificance theTest")
                            anova_table = anova_lm(ols_model, typ=3)

                        # importantaNotemattertermdisplayed（SE typetorespondjite）
                        if se_type == "HC3 (robustEstimation,recommended)":
                            st.warning("⚠️ **ANOVA FTestis etcVariancesignificance andcorrectrulesignificance theassumptiondoinging**。UpperkiofCoefficienttable（HC3robustSE）andis p-valueofmeaningisdifferariing。heterogeneous Varianceisconcerndonerucaseis Coefficienttableofp-valuetheexcellentaheadplease do。")
                        elif se_type == "Cluster-robust (cluster-robust)":
                            st.warning("⚠️ **ANOVA FTestis etcVariancesignificance andcorrectrulesignificance theassumptiondoinging**。UpperkiofCoefficienttable（ClusterrobustSE）andis p-valueofmeaningisdifferariing。within-cluster Correlationorheterogeneous Varianceisconcerndonerucaseis Coefficienttableofp-valuetheexcellentaheadplease do。")
                        elif se_type == "Classical (normal)":
                            st.info("ℹ️ **ANOVA FTestandCoefficienttableis Sameassumption**（etcVariancesignificance ,correctrulesignificance ）usedoinging。heterogeneous Varianceisconcerndonerucaseis 、HC3robustSEofusetheinspectdiscussplease do。")

                        # countValueroundingDisplay
                        anova_display = anova_table.copy()
                        if 'sum_sq' in anova_display.columns:
                            anova_display['sum_sq'] = anova_display['sum_sq'].round(4)
                        if 'F' in anova_display.columns:
                            anova_display['F'] = anova_display['F'].round(3)
                        if 'PR(>F)' in anova_display.columns:
                            anova_display['PR(>F)'] = anova_display['PR(>F)'].apply(lambda x: f"{x:.4g}" if pd.notna(x) else x)

                        st.dataframe(anova_display, use_container_width=True)
                        st.session_state.ols_anova = anova_table  # zipDownloadfor useSave

                        # ANOVAtableofinterpretationguide
                        with st.expander("📘 ANOVAtableofseedirection（PR(>F)Column meaning）"):
                            st.markdown("""
                            ### PR(>F)Column meaning

                            **eachvariablecountOverallofSignificantsignificance theshowdo（other variablescounttheadjustmentdoneUpperwith）**

                            - ✅ **Covariateadjustmentpresent**: Allofvariablecountthesametimetothinkconsiderdoneresult
                            - koof p-valueis 、other variablescountofshadoweffectthetakeriexcludeiedUpperwith、soofvariablecountOverallisSignificanttocontributegivedoingkatheTestdoinging

                            ---

                            ### Coefficienttableand differi

                            | table | Testtarget | Covariateadjustment |
                            |----|---------|----------|
                            | **Coefficienttable** | eachlevelper （e.g.: Male vs Female） | ✅ present |
                            | **ANOVAtable** | variablecountOverall（e.g.: sexOverallofeffect） | ✅ present |
                            """)

                            # realoccasionofModelcontained invariablecountfromDescriptiongenerate
                            anova_vars = []
                            for idx in anova_display.index:
                                idx_str = str(idx)  # Ensure idx is a string
                                if idx_str != 'Residual' and 'C(' in idx_str:
                                    anova_vars.append(idx_str)

                            if anova_vars:
                                st.markdown("**thisModelexample:**")
                                for var_name in anova_vars[:2]:
                                    # C(sex) ofyouashapeformulafromvariablecountNametheextract
                                    base_var = var_name.replace('C(', '').replace(')', '')

                                    st.markdown(f"""
                                    - **`{var_name}`** of PR(>F)
                                      - **meaning**: {base_var}andiuvariablecountOverallis、PCValuetoSignificantashadoweffectthegiveeteiruka
                                      - **TestInsidecontent**: Allof{base_var}levelofCoefficientissametimetozerowhether
                                      - **adjustment**: other variablescount（Main variable2、Covariateetc）ofshadoweffectexcludingiedUpperinTest
                                    """)

                            st.markdown("""
                            ---

                            ### toolbodye.g.withreasonunderstanddo

                            **Model**: `PC4 ~ C(sex) + C(subtype)`

                            | variablecount | PR(>F) | meaning |
                            |------|--------|------|
                            | C(sex) | 0.0234 | subtypeadjustmentAftermo、sexis SignificanttoPC4toshadoweffectdo |
                            | C(subtype) | 0.0567 | sexadjustmentAftermo、subtypeis SignificanttoPC4toshadoweffectdo |

                            ### Welch tTestand differi

                            | Testdirectionmethod | Covariateadjustment | Testtarget | Displayplaceplace |
                            |---------|----------|---------|---------|
                            | Welch tTest | ❌ none | Main variable1only（2groupcomparison） | inforeportbokusu |
                            | ANOVAtable PR(>F) | ✅ present | eachvariablecountOverall | thistable |

                            ### Notematterterm

                            ⚠️ **ANOVA FTestis etcVariancesignificance andcorrectrulesignificance theassumptiondoinging**

                            heterogeneous Varianceisconcerndonerucaseis 、Coefficienttableofp-value（HC3/ClusterrobustSEusetime）theexcellentaheadplease do。

                            ### Type II vs Type III

                            - **Type II ANOVA（recommended）**: InteractiontheincludemanotModelwitheach variable Main effecttheTest
                            - **Type III ANOVA**: AllofInteractiontheincludemeedeach variable effecttheTest

                            **recommended**: correctconfirmavariablecountOverallofSignificantsignificance tois 、thisANOVAtableofPR(>F)useplease do。
                            """)

                    except Exception as e:
                        st.warning(f"ANOVAtablethecalculationcannotwithdone: {str(e)}")

                    # Forest plot
                    st.markdown("### 🌲 Forest Plot (CoefficientofVisualization)")
                    st.caption("eachCoefficientestimationvalue and 95%Confidence interval。zeroandweightifnotcaseis Significant")

                    # InterceptexcludingOutside
                    coef_for_plot = coef_table.iloc[1:].copy()
                    if len(coef_for_plot) > 0:
                        fig_forest = plot_forest(coef_for_plot, title=f"OLS Coefficient Estimates ({se_type})")
                        st.pyplot(fig_forest)
                        st.session_state.fig_forest_ols = fig_forest  # zipDownloadfor useSave
                    else:
                        st.info("plotPossibleaCoefficientispresentnot（Interceptonly）")

                    # diagnosisplot
                    if show_diagnostics:
                        st.markdown("### 🔬 diagnosisplot")
                        st.caption("""
                        - **Residuals vs Fitted**: pattern-nisakerebalinearsignificance andetcVariancesignificance issatisfydoneteiru
                        - **Normal Q-Q**: straightlineUppertoriderebacorrectrulesignificance issatisfydoneteiru
                        - **Scale-Location**: waterflatabeltstatewiththatbaetcVariancesignificance issatisfydoneteiru
                        - **Residual Distribution**: correctruledistributiontonearikaconfirm
                        """)
                        fig_diag = plot_diagnostic(ols_model, title="OLS Model Diagnostics")
                        st.pyplot(fig_diag)
                        st.session_state.fig_diagnostic = fig_diag  # zipDownloadfor useSave

                    # EMM
                    if show_emm and blocking_var:
                        st.markdown("### 📊 Estimated marginal means (Estimated Marginal Means)")
                        st.caption("covariate adjusted Afterof、each group PredictionMeanValue")

                        try:
                            emm_df = calculate_emm(ols_model, analysis_df,
                                                  [main_var, blocking_var] if blocking_var else [main_var])

                            # countValueroundingDisplay
                            emm_display = emm_df.copy()
                            for col in ['mean', 'se', 'ci_lower', 'ci_upper']:
                                if col in emm_display.columns:
                                    emm_display[col] = emm_display[col].round(4)

                            st.dataframe(emm_display, use_container_width=True)

                            # EMM plot
                            fig, ax = plt.subplots(figsize=(10, 6))

                            unique_blocks = emm_df[blocking_var].unique()
                            colors = plt.cm.Set2(np.linspace(0, 1, len(unique_blocks)))

                            for idx, block in enumerate(unique_blocks):
                                block_data = emm_df[emm_df[blocking_var] == block]
                                x_labels = block_data[main_var].values
                                x_pos = np.arange(len(block_data)) + idx * 0.2

                                ax.errorbar(x_pos, block_data['mean'],
                                           yerr=[block_data['mean'] - block_data['ci_lower'],
                                                 block_data['ci_upper'] - block_data['mean']],
                                           label=str(block), fmt='o-', capsize=5,
                                           markersize=8, linewidth=2, color=colors[idx])

                            ax.set_xticks(np.arange(len(x_labels)) + (len(unique_blocks)-1) * 0.1)
                            ax.set_xticklabels(x_labels)
                            ax.set_xlabel(main_var, fontsize=12)
                            ax.set_ylabel(f'Estimated Mean {pc_col}', fontsize=12)
                            ax.set_title('Estimated Marginal Means (95% CI)', fontsize=13, fontweight='bold')
                            ax.legend(title=blocking_var, fontsize=10)
                            ax.grid(True, alpha=0.3)
                            plt.tight_layout()
                            st.pyplot(fig)
                            st.session_state.fig_emm = fig  # zipDownloadfor useSave
                            st.session_state.ols_emm = emm_df  # DatamoSave

                        except Exception as e:
                            st.warning(f"EMMofCalculatingtoErrorisoccurlivediddone: {str(e)}")

                    # resulttheSave
                    st.session_state.ols_model = ols_model
                    st.session_state.ols_coef = coef_table

                except Exception as e:
                    st.error(f"❌ OLSanalysistoFaileddiddone")
                    st.exception(e)

            # =============================================================================
            # Run LMM
            # =============================================================================

            if "LMM" in analysis_method and donor_col:
                st.markdown("---")
                st.markdown("## 6️⃣ Linear mixed model (LMM) analysisresult")

                try:
                    # Random effectofformulathePreparation
                    re_formula = "1"  # randamuIntercept
                    groups = analysis_df[donor_col]

                    # LMMoffuito
                    lmm_model = smf.mixedlm(formula, data=analysis_df,
                                            groups=groups,
                                            re_formula=re_formula)

                    lmm_result = lmm_model.fit(reml=(reml_method == "REML (recommended)"))

                    st.info(f"ℹ️ Estimationdirectionmethod: {reml_method}")

                    # Coefficienttable
                    st.markdown("### 📋 Fixed effectofCoefficienttable")
                    coef_table_lmm = pd.DataFrame({
                        'Coef': lmm_result.params,
                        'Std_Error': lmm_result.bse,
                        'z_value': lmm_result.tvalues,
                        'p_value': lmm_result.pvalues,
                        'CI_lower': lmm_result.conf_int()[0],
                        'CI_upper': lmm_result.conf_int()[1]
                    })

                    coef_table_lmm['Sig'] = coef_table_lmm['p_value'].apply(
                        lambda p: '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else '.' if p < 0.1 else ''
                    )

                    # countValueroundingDisplay
                    lmm_display = coef_table_lmm.copy()
                    lmm_display['Coef'] = lmm_display['Coef'].round(4)
                    lmm_display['Std_Error'] = lmm_display['Std_Error'].round(4)
                    lmm_display['z_value'] = lmm_display['z_value'].round(3)
                    lmm_display['p_value'] = lmm_display['p_value'].apply(lambda x: f"{x:.4g}")
                    lmm_display['CI_lower'] = lmm_display['CI_lower'].round(4)
                    lmm_display['CI_upper'] = lmm_display['CI_upper'].round(4)

                    st.dataframe(lmm_display, use_container_width=True)

                    st.caption("Significantlevel: *** p<0.001, ** p<0.01, * p<0.05, . p<0.1")

                    # Variancecomponent
                    st.markdown("### 📊 Variancecomponent (Variance Components)")
                    st.caption("data variableactioniswherefromcometeirukatheshowdo")

                    var_random = lmm_result.cov_re.values[0][0]
                    var_residual = lmm_result.scale
                    var_total = var_random + var_residual

                    var_comp = pd.DataFrame({
                        'Component': [f'Random effect ({donor_col})', 'Residual (piecebodyInsidevariableaction)'],
                        'Variance': [var_random, var_residual],
                        'Std_Dev': [np.sqrt(var_random), np.sqrt(var_residual)],
                        'Proportion': [var_random / var_total, var_residual / var_total]
                    })

                    # countValueroundingDisplay
                    var_comp_display = var_comp.copy()
                    var_comp_display['Variance'] = var_comp_display['Variance'].round(4)
                    var_comp_display['Std_Dev'] = var_comp_display['Std_Dev'].round(4)
                    var_comp_display['Proportion'] = (var_comp_display['Proportion'] * 100).round(2).astype(str) + '%'

                    st.dataframe(var_comp_display, use_container_width=True)

                    # ICC
                    icc = var_random / var_total
                    st.metric("intraclass correlation coefficient (ICC)", f"{icc:.4f}",
                             help="sameDonorInsideofSampleBetweenofCorrelation。0toneariandRandom effectisSmall、1toneariandDonorBetweenofvariableactionisLarge")

                    if icc < 0.05:
                        st.info("ℹ️ ICC < 0.05: Random effectisnon-alwaysSmallis。OLSalsotenpartkaifrenot。")
                    elif icc > 0.5:
                        st.success("✅ ICC > 0.5: Random effectisLargeis。LMMofuseissuitablecutis。")

                    # ModelGoodness of fit
                    st.markdown("### 📈 ModelGoodness of fit")
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Log-likelihood", f"{lmm_result.llf:.2f}")
                    with col2:
                        st.metric("AIC", f"{lmm_result.aic:.2f}")
                    with col3:
                        st.metric("BIC", f"{lmm_result.bic:.2f}")

                    # Forest plot
                    st.markdown("### 🌲 Forest Plot (Fixed effect)")
                    coef_for_plot_lmm = coef_table_lmm.iloc[1:].copy()
                    if len(coef_for_plot_lmm) > 0:
                        fig_forest_lmm = plot_forest(coef_for_plot_lmm,
                                                     title="LMM Fixed Effects Estimates")
                        st.pyplot(fig_forest_lmm)
                        st.session_state.fig_forest_lmm = fig_forest_lmm  # zipDownloadfor useSave

                    # resulttheSave
                    st.session_state.lmm_result = lmm_result
                    st.session_state.lmm_coef = coef_table_lmm

                except Exception as e:
                    st.error(f"❌ LMManalysistoFaileddiddone")
                    st.exception(e)
                    st.info("hinto: When donor count is lowor、data variableactionisSmallcase、LMMofConvergencetoFaileddothatispresenting。OLSthetryplease do。")

            elif "LMM" in analysis_method and not donor_col:
                st.warning("⚠️ LMMtheRundotois 、Donor/Subject ID thePlease specify。")

            # =============================================================================
            # Welch's t-test
            # =============================================================================

            if run_welch:
                st.markdown("---")
                st.markdown("## 7️⃣ WelchTypetTest results")
                st.caption("etcVariancetheassumptionshinot2groupcomparison（Covariateadjustmentnone）")
                st.warning(f"🎯 **Testtargetvariablecount:** Main variable1 = **{main_var}**")

                # Main variableoflevelcountconfirm
                unique_main_levels = analysis_df[main_var].nunique()
                if unique_main_levels == 2:
                    from scipy import stats
                    groups = analysis_df[main_var].unique()
                    group1_data = analysis_df[analysis_df[main_var] == groups[0]][pc_col]
                    group2_data = analysis_df[analysis_df[main_var] == groups[1]][pc_col]
                    welch_t, welch_p = stats.ttest_ind(group1_data, group2_data, equal_var=False)

                    # kidescribestatistical
                    group1_mean = group1_data.mean()
                    group1_std = group1_data.std()
                    group1_n = len(group1_data)
                    group2_mean = group2_data.mean()
                    group2_std = group2_data.std()
                    group2_n = len(group2_data)

                    # Effect size (Cohen's d, Hedges' g)
                    pooled_std = np.sqrt(((group1_n - 1) * group1_std**2 + (group2_n - 1) * group2_std**2) / (group1_n + group2_n - 2))
                    cohens_d = (group1_mean - group2_mean) / pooled_std
                    correction_factor = 1 - (3 / (4 * (group1_n + group2_n) - 9))
                    hedges_g = cohens_d * correction_factor

                    # Results Display
                    st.success("✅ WelchTypetTesttheRundiddone")

                    col1, col2 = st.columns(2)

                    with col1:
                        st.markdown("### 📊 kidescribestatistical")
                        desc_df = pd.DataFrame({
                            'group': [groups[0], groups[1]],
                            'n': [group1_n, group2_n],
                            'Mean': [f"{group1_mean:.4f}", f"{group2_mean:.4f}"],
                            'marklevelbiasdifference': [f"{group1_std:.4f}", f"{group2_std:.4f}"]
                        })
                        st.dataframe(desc_df, use_container_width=True)

                    with col2:
                        st.markdown("### 📈 Test results")
                        st.metric("tValue", f"{welch_t:.4f}")
                        st.metric("p-value", f"{welch_p:.4g}")
                        if welch_p < 0.001:
                            st.success("***  p < 0.001")
                        elif welch_p < 0.01:
                            st.success("**  p < 0.01")
                        elif welch_p < 0.05:
                            st.success("*  p < 0.05")
                        elif welch_p < 0.1:
                            st.info(".  p < 0.1")
                        else:
                            st.info("n.s.  p ≥ 0.1")

                    st.markdown("### 📏 Effect size")
                    effect_col1, effect_col2 = st.columns(2)
                    with effect_col1:
                        st.metric("Cohen's d", f"{cohens_d:.4f}")
                    with effect_col2:
                        st.metric("Hedges' g (correctionversion)", f"{hedges_g:.4f}")

                    st.caption("Effect sizeofeyesafe: |d| < 0.2 (small)、0.2-0.5 (middle)、0.5-0.8 (large)、≥ 0.8 (non-alwayslarge)")

                    # importantaNotematterterm
                    st.warning(
                        f"⚠️ **importantaNote**\n\n"
                        f"thisTestis  **simple2groupcomparison** is：\n"
                        f"- other variablescount（Main variable2、AddCovariateetc）ofshadoweffectis  **adjustmentdoneteinot**\n"
                        f"- CovariateadjustmentisRequiredacaseis 、OLS/LMM of「Coefficienttable」or 「ANOVAtable」thegoConfirmkudasai\n\n"
                        f"💡 **useipartke**：\n"
                        f"- WelchTypetTest: Main variableonlyofsimplecomparison（referenceValue）\n"
                        f"- OLSCoefficienttable: CovariateadjustmentAfterofeffectEstimation（mainanalysis）\n"
                        f"- Permutation test: distributionassumptionnoneofrobustaTest（supplementstrong）"
                    )

                    # session_statetoSave
                    st.session_state.welch_result = {
                        't': welch_t,
                        'p': welch_p,
                        'cohens_d': cohens_d,
                        'hedges_g': hedges_g,
                        'group1': groups[0],
                        'group2': groups[1],
                        'group1_mean': group1_mean,
                        'group2_mean': group2_mean
                    }

                else:
                    st.error(f"❌ WelchTypetTestis 2groupcomparisononlysupporteddoinging。Main variable1 ({main_var}) is  {unique_main_levels} levelis。")
                    st.info("💡 3groupfollowingUppercase ofis 、OLSofANOVAtableor Permutation testthegousekudasai。")

            # =============================================================================
            # Permutation Test
            # =============================================================================

            if run_permutation:
                st.markdown("---")
                section_number = "8️⃣" if run_welch else "7️⃣"
                st.markdown(f"## {section_number} Permutation test (Permutation Test) result")
                st.caption("parametorikuaassumptiontodependsonshinot、robustaTestdirectionmethod")
                st.warning(f"🎯 **Testtargetvariablecount:** Main variable1 = **{main_var}**")

                with st.spinner(f"{n_permutations:,} timesofplacechangetheRunmiddle..."):
                    try:
                        # designRowColumn Preparation（design_infotakegain）
                        y = analysis_df[pc_col].values
                        X_full = dmatrix(formula.split('~')[1], data=analysis_df, return_type='dataframe')
                        design_info = X_full.design_info  # patsy design_infoobtain

                        # stratifiedvariablecountofPreparation
                        stratify_var = None
                        if blocking_var is not None:
                            stratify_var = analysis_df[blocking_var].values
                            st.warning(f"🔹 **stratified permutationuse**: {blocking_var} InsidewithonlyResidualpermutedo\n\n"
                                      f"thistofrom、{blocking_var}structurebuildthekeepholddonemama{main_var}effectofnulldistributiongeneratedo")
                        else:
                            st.info("ℹ️ Main variable2ispointsetdoneteinotedme、non-stratified permutationusedo")

                        # one-sidedTestofSettings
                        one_sided_param = None
                        if perm_sided == "one-sidedTest: Main variable1 > Criterion (greater)":
                            one_sided_param = 'greater'
                            st.info(f"📊 **one-sidedTest (greater)**: {main_var}ofordinal1level > Criterion（ordinal0level）theTest")
                        elif perm_sided == "one-sidedTest: Main variable1 < Criterion (less)":
                            one_sided_param = 'less'
                            st.info(f"📊 **one-sidedTest (less)**: {main_var}ofordinal1level < Criterion（ordinal0level）theTest")

                        # Permutation testofRun
                        if perm_method == "Freedman-Lanemethod (recommended)":
                            perm_result = freedman_lane_permutation(
                                y, X_full, main_var,
                                n_perm=n_permutations,
                                stratify_by=stratify_var,
                                one_sided=one_sided_param,
                                design_info=design_info  # firmrobustaColumnspecificforAdd
                            )

                            info_msg = f"ℹ️ **Freedman-Lanemethod:** Main variable1（**{main_var}**）ofeffecttheTest\n\n"
                            info_msg += f"other variablescount（Main variable2、AddCovariateetc）is controlvariablecountandandhandlewareing\n\n"
                            if perm_result.get('stratified'):
                                info_msg += f"✅ stratified permutation: {blocking_var}structurebuildthekeephold\n\n"
                            if perm_result.get('one_sided'):
                                info_msg += f"📊 Testdirectiondirection: {perm_result['one_sided']}\n\n"
                            st.success(info_msg)
                        else:
                            # Simple label permutation（yValueplacechangeversion：highfast）
                            np.random.seed(42)  # Freedman-LaneandSameshi-douse

                            # observedStatistic
                            obs_model = sm.OLS(y, X_full).fit()
                            interest_cols = [col for col in X_full.columns
                                           if f'C({main_var})' in col and col != 'Intercept']
                            main_effect_col = [col for col in interest_cols if ':' not in col][0]
                            obs_stat = float(obs_model.tvalues[main_effect_col])

                            null_dist = []
                            for _ in range(n_permutations):
                                # yValuepermute（stratifiedsupported）
                                if stratify_var is not None:
                                    # stratified permutation：eachlayerInsidewithonlyyValueshufflefull
                                    y_perm = np.zeros_like(y)
                                    for stratum in np.unique(stratify_var):
                                        stratum_mask = (stratify_var == stratum)
                                        stratum_indices = np.where(stratum_mask)[0]
                                        perm_indices = np.random.permutation(stratum_indices)
                                        y_perm[stratum_indices] = y[perm_indices]
                                else:
                                    # non-stratified permutation：Overallshufflefull
                                    y_perm = np.random.permutation(y)

                                # fuito（settotalRowColumnis solidset）
                                perm_model = sm.OLS(y_perm, X_full).fit()
                                null_dist.append(float(perm_model.tvalues[main_effect_col]))

                            null_dist = np.array(null_dist)

                            # p-valuecalculation（continuouscontinuesignificance correction）
                            if one_sided_param == 'greater':
                                p_value = (np.sum(null_dist >= obs_stat) + 1) / (n_permutations + 1)
                            elif one_sided_param == 'less':
                                p_value = (np.sum(null_dist <= obs_stat) + 1) / (n_permutations + 1)
                            else:
                                p_value = (np.sum(np.abs(null_dist) >= np.abs(obs_stat)) + 1) / (n_permutations + 1)

                            perm_result = {
                                'observed': obs_stat,
                                'null_distribution': null_dist,
                                'p_value': p_value,
                                'stratified': stratify_var is not None,
                                'one_sided': one_sided_param
                            }

                            info_msg = f"ℹ️ **Simple permutation:** yValuepermute（highfastversion）\n\n"
                            info_msg += "📊 **nullHypothesisofdifferi:**\n"
                            info_msg += "- **Freedman-Lane**: CovariatewithConditiondukeednullHypothesis（Residualpermute）\n"
                            info_msg += "- **Simple**: CovariatestructurebuildthenoneviewdonenullHypothesis（yValuepermute）\n\n"
                            info_msg += "⚠️ **p-valueofdifferi**: nullHypothesisisDifferentedme、p-valuemoDifferentofiscorrectnormalis。\n"
                            if stratify_var is not None:
                                info_msg += f"✅ stratified permutation: {blocking_var}InsidewithonlyyValuetheenterrereplaceeing\n\n"
                            else:
                                info_msg += "allSamplewithyValuetheenterrereplaceeing\n\n"
                            info_msg += "💡 **recommended**: Covariatethecontroldoneicaseis **Freedman-Lanemethod**useplease do。"
                            st.warning(info_msg)

                        # Results Display
                        col1, col2, col3, col4 = st.columns(4)
                        with col1:
                            st.metric("observedStatistic (tValue)", f"{perm_result['observed']:.4f}")
                        with col2:
                            test_type = "one-sided" if perm_result.get('one_sided') else "two-sided"
                            st.metric("Testoftypes", test_type)
                        with col3:
                            st.metric("Permutation p-value", f"{perm_result['p_value']:.4g}")
                        with col4:
                            sig_level = "***" if perm_result['p_value'] < 0.001 else \
                                        "**" if perm_result['p_value'] < 0.01 else \
                                        "*" if perm_result['p_value'] < 0.05 else "n.s."
                            st.metric("Significantsignificance ", sig_level)

                        # stratifiedinforeportofDisplay
                        if perm_result.get('stratified'):
                            st.success(f"✅ **stratified permutation**: {blocking_var}ofeachlevelInsidewithonlyResidualpermutediddone\n\n"
                                      f"each{blocking_var}Insidewith{main_var}effecttheInvalidizationdonenulldistributiongenerate")
                        else:
                            st.info("ℹ️ **non-stratified permutation**: allSampletowaedteResidualpermutediddone")

                        # ReproducibilityinforeportofDisplay（method_info）
                        if 'method_info' in perm_result:
                            with st.expander("🔬 analysisofDetails（Reproducibilityforofinforeport）", expanded=False):
                                info = perm_result['method_info']
                                st.markdown(f"""
**Modelformula:**
- **nullModel**: `{info['reduced_model']}`
- **againststandModel**: `{info['full_model']}`

**TestSettings:**
- **placechangetimescount**: {info['n_permutations']:,} times
- **Random seed**: {info['random_seed']}
- **Statistic**: {info['test_statistic']}
- **continuouscontinuesignificance correction**: {'present (+1correction)' if info['continuity_correction'] else 'none'}
- **stratified**: {info['stratification_var']}

**Testtarget**: `{main_var}` ofeffect（other variablescountis control）
""")

                        # histogram
                        st.markdown("### 📊 nulldistribution (Null Distribution)")
                        st.caption("ashcolor: nullHypothesisofLowerinStatisticofdistribution。redline: observeddoneedStatistic")

                        fig, ax = plt.subplots(figsize=(10, 6))
                        ax.hist(perm_result['null_distribution'], bins=50,
                               alpha=0.7, edgecolor='black', density=True, color='lightgray',
                               label='Null distribution')
                        ax.axvline(perm_result['observed'], color='red',
                                  linestyle='--', linewidth=2.5, label=f'Observed (t={perm_result["observed"]:.3f})')
                        ax.axvline(-perm_result['observed'], color='red',
                                  linestyle='--', linewidth=2.5)

                        # p-valueofterritoryareathepaintritsubushi
                        null_dist = perm_result['null_distribution']
                        obs = perm_result['observed']
                        extreme_vals = null_dist[np.abs(null_dist) >= np.abs(obs)]
                        if len(extreme_vals) > 0:
                            ax.hist(extreme_vals, bins=20, alpha=0.5, color='red',
                                   density=True, label=f'p={perm_result["p_value"]:.4g}')

                        ax.set_xlabel('Test Statistic (t-value)', fontsize=12)
                        ax.set_ylabel('Density', fontsize=12)
                        ax.set_title(f'Permutation Null Distribution ({n_permutations:,} permutations)',
                                    fontsize=13, fontweight='bold')
                        ax.legend(fontsize=10)
                        ax.grid(alpha=0.3)
                        plt.tight_layout()
                        st.pyplot(fig)
                        st.session_state.fig_null_dist = fig  # zipDownloadfor useSave

                        # resulttheSave
                        st.session_state.perm_result = perm_result

                    except Exception as e:
                        st.error(f"❌ Permutation testtoFaileddiddone")
                        st.exception(e)

            # =============================================================================
            # Download Results
            # =============================================================================

            st.markdown("---")
            section_number_download = "9️⃣" if run_welch else "8️⃣"
            st.markdown(f"## {section_number_download} Results Download")
            st.caption("analysisresultandgurafuthemaandmeteZIPFilewithDownloadcaning")

            # ZIPDownloadbuttonthemiddlecentertodistributeplace
            col1, col2, col3 = st.columns([1, 2, 1])
            with col2:
                # resultis1tsufollowingUpperexistcaseonlyDownloadbuttondisplayed
                has_results = any([
                    'ols_coef' in st.session_state,
                    'lmm_coef' in st.session_state,
                    'perm_result' in st.session_state,
                    'welch_result' in st.session_state
                ])

                if has_results:
                    # FileNamegenerate
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    zip_filename = f"pca_stats_{pc_col}_{timestamp}.zip"

                    # ZIPmakebecome
                    zip_data = create_results_zip(pc_col)

                    st.download_button(
                        label="📦 allofresulttheDownload (ZIP)",
                        data=zip_data,
                        file_name=zip_filename,
                        mime="application/zip",
                        help="statisticalte-buru（CSV）andgurafu（PNG 300dpi + PDF）theZIPFiletomaandmeteDownload",
                        type="primary",
                        use_container_width=True
                    )

                    # includemareruInsidecontentdisplayed
                    st.markdown("---")
                    st.markdown("#### 📋 ZIPFileofInsidecontent:")

                    contents = []
                    if 'welch_result' in st.session_state:
                        contents.append("✅ WelchTypetTest results (CSV)")
                    if 'ols_coef' in st.session_state:
                        contents.append("✅ OLSCoefficienttable (CSV)")
                    if 'ols_anova' in st.session_state:
                        contents.append("✅ ANOVAtable (CSV)")
                    if 'ols_emm' in st.session_state:
                        contents.append("✅ Estimated marginal means (CSV)")
                    if 'lmm_coef' in st.session_state:
                        contents.append("✅ LMMCoefficienttable (CSV)")
                    if 'perm_result' in st.session_state:
                        contents.append("✅ Permutation testresult (CSV)")
                    if 'fig_forest_ols' in st.session_state:
                        contents.append("✅ OLS Forest plot (PNG + PDF)")
                    if 'fig_diagnostic' in st.session_state:
                        contents.append("✅ diagnosisplot (PNG + PDF)")
                    if 'fig_emm' in st.session_state:
                        contents.append("✅ EMM plot (PNG + PDF)")
                    if 'fig_forest_lmm' in st.session_state:
                        contents.append("✅ LMM Forest plot (PNG + PDF)")
                    if 'fig_null_dist' in st.session_state:
                        contents.append("✅ Null distribution plot (PNG + PDF)")

                    # 2ColumnwithDisplay
                    col_a, col_b = st.columns(2)
                    half = len(contents) // 2 + len(contents) % 2
                    with col_a:
                        for item in contents[:half]:
                            st.markdown(f"- {item}")
                    with col_b:
                        for item in contents[half:]:
                            st.markdown(f"- {item}")

                else:
                    st.info("💡 analysistheRundoand、resulttheDownloadcanruyoutoariing")

else:
    st.info("👆 DataUpload fileandanalysistheStartplease do")

    st.markdown("---")
    st.markdown("### 📋 Datashapeformulaexample")
    st.markdown("""
    Datais  **TSV** or  **CSV** shapeformulawith、followingLowerofyouastructurebuildthethinksetdoinging:

    - **1Row = 1Sample** (pseudobulk reberu)
    - **PC score Column**: PC1, PC2, PC3, PC4 etc
    - **Categorical variable**: sex, subtype, condition etc
    - **arbitrary: Donor/Subject ID** (Iterationmeasureisexistcase)
    - **arbitrary: countValueCovariate**: age, batch, depth etc

    **SampleDatastructurebuild:**
    ```
    sample_id    sex       subtype    donor_id    age    batch    PC1      PC2      PC3      PC4
    sample_001   Female    TypeA      donor_01    45     1        -2.3     1.2      0.5      -0.8
    sample_002   Male      TypeA      donor_02    38     1         1.5    -0.8     -1.2      0.3
    sample_003   Female    TypeB      donor_03    52     1         0.2     2.1      1.8     -1.5
    sample_004   Male      TypeB      donor_04    41     2        -1.1    -1.5      0.9      1.2
    ...
    ```

    testDatais  `/home/ichiro/streamlit/temp/pca_test_data.tsv` topresenting。
    """)

st.markdown("---")

with st.expander("📚 Statistical methodsofDetails", expanded=False):
    st.markdown("""
#### **OLS (Ordinary Least Squares - Ordinary Least Squares)**
Fixed effectModel。AllofSampleisindependentwithexistandassumptiondo。

- **Classical SE**: markleveltekiaStandard error（etcVariancesignificance theassumption）
- **HC3 Robust SE**: heterogeneous Variancerobust to（SamplesizeisSmallcasetorecommended）
- **Cluster-robust SE**: within-cluster Correlationconsidering （e.g.: sameDonorofSampleBetweenofCorrelation）
- **WelchType**: etcVariancetheassumptionshinot（2groupcomparisononly）
- **WLS (addweightsaismalltworidemethod)**: heterogeneous Variancesupported

**Applye.g.**: eachDonorfrom1Sampleonly、or DonorBetweenofvariableactionisSmallcase

---

#### **LMM (Linear Mixed Model - Linear mixed model)**
Hierarchical structuretheholdtsuDatasupported。Random effectwithpiecebodyvariation betweenconsidering do。

- **randamuIntercept `(1|donor)`**: eachDonorwithInterceptisDifferentthatthepermitcontent
- **REML**: unbiased VarianceEstimation（recommended）
- **ML**: Modelfor comparison（AIC/BICofcomparison）
- **ICC (intraclass correlation coefficient)**: Donorvariation betweenofdividegotheshowsu

**Applye.g.**: eachDonorfrommultiplecountSample、or skilltechniquetekirepeatrireturnshiisexistcase

**Note**: Donorcountis5less thancase of、VarianceEstimationisunstabletoariing。

---

#### **Permutation Test (Permutation test)**
distributionassumptiontodependsonshinot、nonparametorikuaTestdirectionmethod。

**🎯 Testtarget:** Main variable1onlyisTestdoneing。Main variable2orsoofotherofCovariateis controlvariablecountandandhandlewareing。

- **Freedman-Lanemethod**: Covariatestrictly controlling Permutation test（recommended）
  - Main variable1ofeffecttheTest
  - nullModelandagainststandModelthefuito
  - nullModelofResidualpermuteandquasisimilarDatagenerate
  - againststandModeltheagainfuitoandStatisticthecalculation

- **Simple permutation**: raberuthesimpletoenterrereplaceerudirectionmethod
  - Main variable1ofraberutheenterrereplacee
  - Covariateisnot、or Fewcasetouse

**Applye.g.**:
- SamplesizeisSmall（<30）
- Residualnormalisignificance isquestionable
- Outlierissonexistdo

---

#### **ANOVA Type II vs Type III**
- **Type II** (recommended): each variable Main effecttheother variablescountwithadjustmentandTest（Interactionnonecase of）
- **Type III**: AllofInteractiontheincludemeteTest（Interactionpresentcase of）

---

#### **EMM (Estimated Marginal Means - Estimated marginal means)**
covariate adjusted Afterof、each group PredictionMeanValue。

- countValueCovariatetheMeanValuefixed to
- Categorical variableofeachlevelinPredictionValuethecalculation
- Confidence intervaltheattachgiveandcomparison

**Applye.g.**:
- subtype  significance differencevisualization
- ageorbatchtheadjustmentdonegroupBetweencomparison

---

### ⚠️ importantaNotematterterm

#### **1. completecollinearity  (Perfect Collinearity)**
2ofvariablecountiscompletetoonereachdoingcase、effectthepartseparatecannot。

**e.g.**: exist subtype isAll Female case of、sex effectand subtype effectis recognizebynotability

**againstprocessmethod**:
- anykaofvariablecounttheDelete
- fromaveragebalanceoftakeedDatathecollectgather
- kidescribetekiacomparisontostaymeru

#### **2. Samplesize**
- **LMMofDonorcount**: 5followingUpperrecommended（FewandVarianceEstimationisunstable）
- **eachcellofsaiSmall samplecount**: 3followingUpperrecommended
- **Permutation test**: SamplesizeisSmallcase（<30）tospecialtoexistuse

#### **3. Interactionterm**
Dataistenparttoexistcaseonlyuserecommended。eachcelltotenpartaSample（≥5）isRequired。

#### **4. multipleTest**
Multiple PCtheanalysisdocaseis 、**Benjamini-Hochberg FDRcorrection**theApplyplease do。

#### **5. countValueCovariateofmarklevelization**
ageetc、suke-ruislargekikuDifferentcountValueCovariateis 、matterBeforetomarklevelizationdothatrecommendeddo。

---

### 📖 referencetextdedicate

1. **Robust standard errors**:
   - Long, J. S., & Ervin, L. H. (2000). Using heteroscedasticity consistent standard errors in the linear regression model. *The American Statistician*, 54(3), 217-224.

2. **Linear Mixed Models**:
   - Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting linear mixed-effects models using lme4. *Journal of Statistical Software*, 67(1), 1-48.

3. **Freedman-Lane permutation**:
   - Freedman, D., & Lane, D. (1983). A nonstochastic interpretation of reported significance levels. *Journal of Business & Economic Statistics*, 1(4), 292-298.
   - Winkler, A. M., et al. (2014). Permutation inference for the general linear model. *NeuroImage*, 92, 381-397.

4. **Estimated Marginal Means**:
   - Searle, S. R., Speed, F. M., & Milliken, G. A. (1980). Population marginal means in the linear model: An alternative to least squares means. *The American Statistician*, 34(4), 216-221.

---

**openoccurperson**: Claude Code
**version**: 1.0
**FinalUpdate**: 2025-01
""")
