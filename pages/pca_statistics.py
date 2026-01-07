"""
PCA Score Statistical Analysis with General Linear Models
=========================================================
PCAスコアの統計解析ツール - 複数の共変量を含む複雑なデザインに対応

このツールは、pseudobulkまたはサンプルレベルのデータに対して、
PC scoreを応答変数とした統計解析を行います。

主な機能:
- OLS (通常最小二乗法) with robust standard errors
- LMM (線形混合モデル) with random effects
- Freedman-Lane permutation testing
- Estimated Marginal Means (EMM)

適用例:
- scRNA-seq pseudobulkデータのPC scoreと性別・細胞タイプの関連
- バッチ効果や年齢などの共変量を調整した解析
- ドナー間変動を考慮した階層モデル
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
### PC scoreを応答変数とした一般化線形モデル解析

このアプリは、PCA（主成分分析）で得られたPC scoreに対して、
**性別・細胞サブタイプ・バッチ・年齢**などの効果を統計的に検定します。

---
""")

with st.expander("📚 Feature Details (Click to expand)", expanded=False):
    st.markdown("""
#### 📌 主な用途:
- **scRNA-seq pseudobulk解析**: サンプル単位に集計したデータのPC scoreと表現型の関連
- **不均衡デザイン対応**: 群サイズが異なる、欠損がある場合にも対応
- **階層構造**: ドナーや技術的バッチなどの階層を考慮
- **共変量調整**: 年齢、深度、品質指標などを調整

#### 🔬 実装されている統計手法:

**1. OLS (Ordinary Least Squares - 通常最小二乗法)**
- 固定効果モデル（sex, subtype, batch, age等）
- **HC3 robust SE**: 不均一分散に頑健な標準誤差（推奨）
- **Cluster-robust SE**: クラスター内相関を考慮（例: donor単位）
- Type II ANOVA: 各変数の有意性検定

**2. LMM (Linear Mixed Model - 線形混合モデル)**
- **ランダム効果**: `(1|donor)`, `(1|batch)` など階層構造に対応
- **REML推定**: 不偏な分散推定（推奨）
- **ML推定**: モデル比較用（AIC/BIC）
- **ICC計算**: 級内相関係数（ランダム効果の寄与率）
- ドナー数が少ない場合（<5）は警告を表示

**3. Permutation Test (置換検定)**
- **Freedman-Lane法**: 共変量を制御した厳密な検定（推奨）
- **Simple permutation**: シンプルなラベル入替
- 小サンプルや分布仮定が疑わしい場合に有効
- 反復数: 1,000〜50,000（カスタマイズ可能）

**4. EMM (Estimated Marginal Means - 推定周辺平均)**
- 共変量を調整した群平均の推定
- サブタイプ別の性差などを可視化
- 信頼区間付きで解釈しやすい

#### ⚙️ データ品質チェック:
- **完全共線性の検出**: 推定不可能な効果を事前に警告
- **VIF計算**: 数値共変量の多重共線性診断
- **サンプルサイズ確認**: 各群のサンプル数とクロス集計表
- **欠損値処理**: 自動除去と報告

#### 📊 可視化:
- **Forest plot**: 係数推定値と95%信頼区間
- **診断プロット**: 残差、Q-Q plot、Scale-Location
- **EMM plot**: 群別の推定平均値
- **Permutation histogram**: 帰無分布と観測統計量

#### 💾 結果のダウンロード:
- 係数表（CSV）
- ANOVA表（CSV）
- 分散成分（LMM）
- Permutation分布（CSV）
- 全ての図（Streamlitから保存可能）
""")

# =============================================================================
# Practical Guide
# =============================================================================
with st.expander("📖 Practical Guide：Exploring sex difference axis in cell-type pseudobulk data", expanded=False):
    st.markdown("""
### 🎯 ユースケース：cell.type が4種類 × sex の8サンプルで「sex軸」を見つける

#### **入力データの前提**

**1行 = 1サンプル（pseudobulk）**

列の例：
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

#### **アプリ側の設定（列マッピング）**

1. **主要変数 1 (main_var)** = `sex`
2. **主要変数 2 / ブロッキング変数 (blocking_var)** = `cell.type`
3. **解析PC** = PC1〜PCk を順に選択（または後述の"sex 合成軸"を一発で作る）

---

#### **手順A：PCごと（1次元ずつ）の"sex主効果"スクリーニング**

**目的:** どのPCが sex を最もよく説明するかを、cell.type を共変量として調整しながら判定。

**ステップ:**

1. **PCを1つ選ぶ**（例：PC3）

2. **解析タブでOLSを選択**し、モデル式を以下のように構築：
   ```
   PC3 ~ C(sex) + C(cell.type)
   ```
   - 交互作用は入れない（主効果のみ）
   - Standard Error は **HC3 Robust** を推奨（小サンプルに頑健）

3. **出力の確認**：
   - **sex 係数**（Female vs Male）とその HC3ロバストSEの **p値** を確認
   - 係数の**符号**と**効果量**（coefficient の絶対値）を記録

4. **Permutation Test をON**にして、より正確なp値を取得：
   - **Freedman-Lane法** を選択（推奨）
   - 層（= cell.type）を維持したまま sex ラベルを入替
   - 小サンプル（n=8）でも頑健な検定が可能

5. **EMM（推定周辺平均）プロット**で一貫性を確認：
   - cell.type ごとの Female vs Male の差の**符号が揃っているか**確認
   - 全てのcell.typeで同じ方向の差があれば、sexの主効果が一貫している証拠

6. **PC1〜PCk で繰り返し**、以下を総合評価：
   - ✅ **効果量**（sex係数の大きさ）
   - ✅ **一貫性**（subtypeごとの符号の揃い方：EMMプロットで確認）
   - ✅ **p値**（HC3 と Permutation の両方）

   → これらを総合して「**sexに最も寄与するPC**」を選ぶ

---

#### **📝 統計的解釈**

このOLSモデルは、**対応ありt検定（Δ = Female − Male の平均）** と本質的に等価です：

- **自由度**: df = cell.type数 − 1 = 3
- **検定対象**: cell.typeを層として制御した上での sex の主効果
- **小標本でも頑健**: HC3 Robust SE と Permutation を併用することで、n=8でも信頼できる推測が可能

**注意点:**
- cell.typeごとにサンプル数が等しい（バランスデザイン）場合、最も統計検出力が高い
- 交互作用 `C(sex) * C(cell.type)` を入れる場合、各セル（sex × cell.type）に最低3サンプル必要
- 本ユースケースでは各セルに1サンプルずつしかないため、交互作用は推定不可（完全分離）

---

#### **🔍 結果の読み方の例**

| PC  | sex係数 | HC3 p値 | Perm p値 | 一貫性 | 総合評価 |
|-----|---------|---------|----------|--------|----------|
| PC1 | 0.45    | 0.234   | 0.251    | ✅ 全て+ | 弱い    |
| PC2 | -1.87   | 0.032   | 0.041    | ✅ 全て- | **強い** |
| PC3 | 0.23    | 0.678   | 0.712    | ± 混在  | 無し    |
| PC4 | 1.12    | 0.098   | 0.112    | ✅ 全て+ | 中程度  |

→ この例では **PC2 が sex 軸** として最も強い（負の方向に分離）

---

#### **💡 次のステップ（応用）**

1. **複数PCを統合した"sex score"の作成**
   - 効果量が大きいPC（例: PC2, PC4）を線形結合
   - `sex_score = -1.87 × PC2 + 1.12 × PC4` のように重み付け平均
   - この合成軸を新たな応答変数として再解析

2. **個別遺伝子レベルへのドリルダウン**
   - sex軸（例: PC2）の loading が大きい遺伝子を抽出
   - 性差に寄与する遺伝子群として解釈

3. **多重検定補正**
   - 複数PC（例: PC1〜PC10）を探索した場合、Benjamini-Hochberg FDR補正を適用
   - 10個のPCを検定する場合、p < 0.05 → FDR < 0.05 に補正

---

#### **📚 このユースケースで使う機能**

- ✅ OLS with HC3 Robust SE
- ✅ Freedman-Lane Permutation Test
- ✅ EMM（推定周辺平均）プロット
- ✅ 診断プロット（残差の正規性・等分散性）
- ✅ 係数テーブルのダウンロード（複数PCの結果を統合）

---
""")

# =============================================================================
# Practical Guide 2: 時系列 × 遺伝子型
# =============================================================================
with st.expander("📖 Practical Guide2：WT vs KO × 時系列3点での交互作用解析", expanded=False):
    st.markdown("""
### 🎯 ユースケース：WT/KO × 3時点（各n=3）で「genotype効果が時間依存か」を検定

#### **入力データの前提**

**1行 = 1サンプル（独立サンプル想定）**

- **デザイン**: 2遺伝子型 (WT/KO) × 3時点 (t1/t2/t3) × 各セルn=3 = **計18サンプル**
- **仮定**: 各時点で異なる個体（独立サンプル）
- **もし同一個体の繰返し測定なら** → 下記「繰返し測定の場合」を参照

列の例：
```
sample_id    genotype    time    donor_id    PC1      PC2      PC3      PC4      ...
WT_t1_1      WT          t1      D01        -1.23     0.45    -0.67     1.12
WT_t1_2      WT          t1      D02         0.87    -1.34     0.92    -0.45
WT_t1_3      WT          t1      D03        -0.45     1.23    -1.01     0.78
WT_t2_1      WT          t2      D04         1.56    -0.23     1.34    -0.89
...（以下、KO_t1, KO_t2, KO_t3も同様）
```

---

#### **統計モデルの戦略**

**重要**: まず**交互作用を中心に検定** → 必要なら単純効果へ掘る

1. **交互作用が有意** → genotype効果は時点で異なる（時間依存）
   - 各時点でのWT vs KOを個別に検定（単純効果）
   - 各genotype内での時系列変化を検定

2. **交互作用が非有意** → genotype効果は時点に依存せず一定
   - 交互作用を落として主効果のみのモデルで解釈
   - genotype主効果とtime主効果を個別に評価

---

#### **アプローチA：timeをカテゴリとして扱う（一般解・推奨）**

**モデル式:**
```
PC ~ C(genotype) * C(time)
```

**特徴:**
- 各時点を独立したカテゴリとして扱う
- 非線形な時間変化にも対応
- 3時点のみの場合、これが最も柔軟で安全

---

##### **🔧 アプリでの設定手順（カテゴリ版）**

**1. データ準備**
- `genotype` 列: "WT", "KO"
- `time` 列: "t1", "t2", "t3" （文字列として）
- `PC1`, `PC2`, ... 列

**2. 列マッピング**
- **主要変数 1**: `genotype`
- **主要変数 2 / ブロッキング変数**: `time`
- **解析PC**: PC1〜PCk から選択
- **Donor/Subject ID**: （独立サンプルなら不要）

**3. モデル設定**
- **モデルタイプ**: OLS
- **Standard Error**: **HC3 Robust** （推奨、小サンプルに頑健）
- **交互作用項**: `genotype:time` を追加 ← **重要！**

**4. 検定オプション**
- **ANOVA Type**: Type II（推奨）
- **Permutation Test**: ON（Freedman-Lane法、推奨）
  - n=18は小サンプルなので置換検定で補強

**5. 可視化**
- **EMM（推定周辺平均）**: ON
  - genotype と time の両方を選択
  - プロットで時系列トレンドの群間差を確認
- **診断プロット**: ON（残差の正規性・等分散性確認）

---

##### **📊 結果の読み方（カテゴリ版）**

**Step 1: ANOVA表で交互作用を確認**

| 項 | Sum Sq | df | F value | p value |
|----|--------|----|---------| --------|
| C(genotype) | 15.3 | 1 | 8.23 | 0.012 |
| C(time) | 45.7 | 2 | 12.34 | 0.001 |
| **C(genotype):C(time)** | **23.4** | **2** | **6.32** | **0.011** ← 重要 |
| Residual | 22.3 | 12 | - | - |

**判断:**
- **交互作用 p=0.011 < 0.05** → **有意**
- → **genotype効果は時点で異なる**（時間依存のKO効果）

**Step 2: EMM表で各時点の群間差を確認**

| genotype | time | mean | SE | 95% CI lower | 95% CI upper |
|----------|------|------|----|--------------| -------------|
| WT | t1 | -0.45 | 0.23 | -0.95 | 0.05 |
| KO | t1 | -0.52 | 0.23 | -1.02 | -0.02 |
| WT | t2 | 0.87 | 0.23 | 0.37 | 1.37 |
| KO | t2 | 1.89 | 0.23 | 1.39 | 2.39 |
| WT | t3 | 1.34 | 0.23 | 0.84 | 1.84 |
| KO | t3 | 2.78 | 0.23 | 2.28 | 3.28 |

**解釈:**
- **t1**: KO - WT = -0.07（ほぼ差なし）
- **t2**: KO - WT = 1.02（KOが増加）
- **t3**: KO - WT = 1.44（差がさらに拡大）

→ **KO効果は時間とともに増大**（t1では差なし → t3で大きな差）

**Step 3: EMMプロットで可視化**

時間を横軸、PC scoreを縦軸にして：
- WT（青線）は緩やかに上昇
- KO（赤線）は急激に上昇
- 2つの線が平行でない = 交互作用あり

---

##### **📝 単純効果の検定（交互作用が有意な場合）**

交互作用が有意なら、**各時点でのWT vs KOを個別に検定**：

**方法1: サブセット解析（手動）**
1. データを時点ごとに分割（t1のみ、t2のみ、t3のみ）
2. 各サブセットで `PC ~ C(genotype)` を検定
3. **Bonferroni補正**: p値の閾値を 0.05/3 = 0.0167 に設定
   または **BH-FDR補正**を適用

**方法2: EMM contrast（推奨、将来実装予定）**
- 自動的に各時点でのペアワイズ比較を計算
- 多重検定補正を自動適用

**現在のアプリでの実施方法:**
1. 元のファイルをExcelで時点ごとに分割（t1.tsv, t2.tsv, t3.tsv）
2. 各ファイルを個別にアップロードして `PC ~ C(genotype)` で検定
3. 3つのp値を手動でBH補正

---

#### **アプローチB：timeを数値トレンド（線形）として扱う**

**モデル式:**
```
PC ~ C(genotype) * time_numeric
```

**特徴:**
- 時間を連続変数（0, 1, 2）として扱う
- **交互作用係数 = 傾きの差**（KOとWTで時間勾配が違うか）
- 1つの係数で時間依存性を要約できる
- **注意**: 線形を仮定（3点しかないので2次は過学習リスク）

---

##### **🔧 アプリでの設定手順（数値版）**

**1. データ準備**
- `time_numeric` 列を追加: t1→0, t2→1, t3→2
- または `time` 列を直接 0, 1, 2 に変更

**2. 列マッピング**
- **主要変数 1**: `genotype`（カテゴリ）
- **追加の共変量**: `time_numeric`（数値として認識される）
- **交互作用項**: `genotype:time_numeric` を追加

**3. モデル設定**
- 他はアプローチAと同じ

---

##### **📊 結果の読み方（数値版）**

**係数表:**

| 項 | Coef | SE | t | p value |
|----|------|----|---|---------|
| Intercept | -0.45 | 0.18 | -2.50 | 0.027 |
| C(genotype)[T.KO] | -0.07 | 0.25 | -0.28 | 0.783 |
| time_numeric | 0.90 | 0.15 | 6.00 | <0.001 |
| **C(genotype)[T.KO]:time_numeric** | **0.75** | **0.21** | **3.57** | **0.003** ← 重要 |

**解釈:**
- **time_numeric係数 = 0.90**: WT群では時間1単位あたりPC scoreが0.90増加
- **交互作用係数 = 0.75**: KO群では時間1単位あたり**さらに0.75多く**増加
  - つまりKO群の傾き = 0.90 + 0.75 = **1.65**
  - WT群の傾き = 0.90
- **交互作用 p=0.003** → **KOとWTで時間トレンドが有意に異なる**

**利点:**
- 1つの検定で「傾きの差」を評価
- 効果量が直感的（傾きの差分）

**欠点:**
- 線形仮定（3点では検証困難）
- 非線形な変化を見逃す可能性

---

#### **アプローチC：繰返し測定（同一個体を追跡）の場合**

**前提:** 同じドナー（D01, D02, D03）を3時点で測定

**モデル式（LMM）:**
```
PC ~ C(genotype) * C(time) + (1|donor)
```

または時間傾きもランダムに：
```
PC ~ C(genotype) * time_numeric + (1 + time_numeric|donor)
```

**注意:**
- ドナー数が少ない（<5）場合、ランダム傾き `(1+time|donor)` は不安定
- まずは `(1|donor)` から始める

---

##### **🔧 アプリでの設定手順（LMM版）**

**1. データ準備**
- `donor` 列が必須
- 同一ドナーIDが複数行（3行）に現れる

**2. 列マッピング**
- **主要変数 1**: `genotype`
- **主要変数 2**: `time`
- **Donor/Subject ID**: `donor` ← **重要！**
- **交互作用項**: `genotype:time`

**3. モデル設定**
- **モデルタイプ**: **LMM**
- **推定法**: REML（推奨）
- 他の設定はOLSと同様

---

##### **📊 結果の読み方（LMM版）**

**固定効果表:**

| 項 | Coef | SE | z | p value |
|----|------|----|---|---------|
| Intercept | -0.45 | 0.28 | -1.61 | 0.108 |
| C(genotype)[T.KO] | -0.07 | 0.40 | -0.18 | 0.860 |
| C(time)[T.t2] | 1.32 | 0.20 | 6.60 | <0.001 |
| C(time)[T.t3] | 1.79 | 0.20 | 8.95 | <0.001 |
| C(genotype)[T.KO]:C(time)[T.t2] | 1.09 | 0.28 | 3.89 | <0.001 |
| C(genotype)[T.KO]:C(time)[T.t3] | 1.51 | 0.28 | 5.39 | <0.001 |

**ランダム効果:**
- `Var(donor)` = 0.45（ドナー間のばらつき）
- `Var(Residual)` = 0.23（測定内のばらつき）
- **ICC = 0.45 / (0.45 + 0.23) = 0.66**
  - 66%の変動がドナー間の個体差で説明される
  - 繰返し測定の相関が強い → LMMが適切

**解釈:**
- 交互作用が有意 → genotype効果は時間依存
- ドナー間変動を適切に考慮した上での結論

---

#### **🔍 レポート作法（推奨）**

**1. 交互作用の結論を最優先**
```
「genotype × time の交互作用が有意であった（F(2,12)=6.32, p=0.011, HC3 robust SE;
Freedman-Lane permutation p=0.015）。これは、KO効果が時間依存であることを示す。」
```

**2. EMM図を作成**
- 時間を横軸、PC scoreを縦軸
- WTとKOの平均 ± 95%CI を折れ線で
- 可能なら **Δ(t) = KO - WT** の差分プロットを追加

**3. 単純効果（交互作用が有意な場合）**
```
「事後検定として各時点でのWT vs KOを比較した結果：
- t1: Δ=-0.07 (95%CI: -0.68, 0.54), p=0.78（有意差なし）
- t2: Δ=1.02 (95%CI: 0.41, 1.63), p=0.004（有意）
- t3: Δ=1.44 (95%CI: 0.83, 2.05), p<0.001（有意）
Bonferroni補正後も t2, t3 は有意を維持。」
```

**4. 効果量を報告**
- 各時点の差の推定値と95%CI
- 全体の genotype 効果（カテゴリ版ならη², 数値版なら傾きの差）

**5. Permutation pを併記**
- HC3と整合すれば頑健性の証拠

---

#### **⚠️ 注意事項**

**1. サンプルサイズ**
- 各セルn=3は最小限（統計検出力が低い）
- 可能なら n≥5 を推奨
- **Permutation testで補強**することを強く推奨

**2. 多重検定補正**
- 単純効果（3比較）には **Bonferroni** または **BH-FDR** を適用
- 複数PCを検定する場合も同様

**3. 線形 vs カテゴリ**
- **3点のみなら線形仮定は危険** → カテゴリ版を推奨
- 線形版は exploratory に使い、カテゴリ版で確認

**4. 完全分離（Perfect Separation）**
- 交互作用項を入れると 2×3=6セル × n=3 = 18サンプル
- df = 18 - 6 - 1 = 11（余裕は少ない）
- 欠損や外れ値があると推定不安定 → 診断プロットで確認

**5. 繰返し測定の場合**
- ドナー数が少ない（<5）と `(1+time|donor)` は不安定
- まずは `(1|donor)` で検定、収束すれば傾きも試す

---

#### **📚 このユースケースで使う機能**

- ✅ OLS with HC3 Robust SE（独立サンプル）
- ✅ LMM with (1|donor)（繰返し測定）
- ✅ Type II ANOVA（交互作用の検定）
- ✅ Freedman-Lane Permutation Test（小サンプルの補強）
- ✅ EMM（推定周辺平均）プロット（時系列トレンドの可視化）
- ✅ 診断プロット（仮定の確認）
- ⚠️ 単純効果の自動計算（将来実装予定、現在は手動でサブセット解析）

---

#### **💡 実践例のまとめ**

| 状況 | 推奨モデル | 交互作用項 | 検定法 | EMM |
|------|-----------|----------|--------|-----|
| 独立サンプル、カテゴリ | `PC ~ C(genotype) * C(time)` | あり | OLS + HC3 + Perm | genotype × time |
| 独立サンプル、線形 | `PC ~ C(genotype) * time_numeric` | あり | OLS + HC3 + Perm | 傾きの差 |
| 繰返し測定 | `PC ~ C(genotype) * C(time) + (1\|donor)` | あり | LMM (REML) | genotype × time |

---

**ひとことで:**
1. **まず交互作用 `genotype × time` を検定**
2. **有意なら各時点の単純効果へ**（現在は手動サブセット解析）
3. **非有意なら主効果モデルで解釈**
4. **小サンプル（n=3/セル）なので Permutation で補強必須**

この流れで、WT/KO × 3時点 × n=3 という現実的な規模でも、**時間依存のgenotype効果**を統計的に評価できます。

---
""")

st.markdown("---")

# =============================================================================
# Helper Functions
# =============================================================================

def remove_common_suffix(strings):
    """末尾の共通要素を除去

    DESeq2-LRT.pyと同じアルゴリズム
    """
    if not strings or len(strings) == 0:
        return []
    # 最も短い文字列の長さを取得
    min_length = min(len(s) for s in strings)
    # 共通の末尾部分の長さを見つける
    suffix_length = 0
    for i in range(1, min_length + 1):
        suffix = strings[0][-i:]
        if all(s.endswith(suffix) for s in strings):
            suffix_length = i
        else:
            break
    # 共通の末尾部分が見つからない場合は元のリストを返す
    if suffix_length == 0:
        return strings
    # 共通の末尾部分を削除して新しいリストを作成
    return [s[:-suffix_length] for s in strings]

def detect_collinearity(df, formula_terms):
    """完全共線性の検出

    カテゴリ変数間で完全に共役している（一方が決まれば他方も決まる）
    場合を検出します。このような場合、モデルは推定できません。

    例: あるsubtypeが全てFemaleの場合、sex効果とsubtype効果は分離できません。
    """
    issues = []

    # カテゴリ変数のみをチェック
    cat_vars = [col for col in formula_terms if col in df.columns and df[col].dtype == 'object']

    for i, var1 in enumerate(cat_vars):
        for var2 in cat_vars[i+1:]:
            cross_tab = pd.crosstab(df[var1], df[var2])
            # 各行または各列に1つしか非ゼロ要素がない場合は完全共役
            if (cross_tab > 0).sum(axis=0).min() == 1 or (cross_tab > 0).sum(axis=1).min() == 1:
                issues.append(f"⚠️ **完全共役 (Perfect confounding)**: `{var1}` と `{var2}` が完全に一致しています")

    return issues

def calculate_vif(df, numeric_cols):
    """VIF (Variance Inflation Factor) の計算

    VIF > 10: 強い多重共線性あり
    VIF > 5: 中程度の多重共線性あり
    VIF < 5: 問題なし
    """
    from statsmodels.stats.outliers_influence import variance_inflation_factor

    vif_data = pd.DataFrame()
    vif_data["Variable"] = numeric_cols
    vif_data["VIF"] = [variance_inflation_factor(df[numeric_cols].values, i)
                       for i in range(len(numeric_cols))]
    return vif_data

def freedman_lane_permutation(y, X_full_df, interest_var_name, n_perm=10000, random_state=42,
                               stratify_by=None, one_sided=None, design_info=None):
    """Freedman-Lane法による置換検定（堅牢版）

    共変量を制御しながら、特定の変数（例: sex）の効果を検定します。

    手順:
    1. フルモデル（全ての変数）と帰無モデル（検定対象以外）をフィット
    2. 帰無モデルの残差を置換（層化オプションあり）
    3. 置換された残差から擬似応答変数を構築
    4. フルモデルを再フィットして統計量を計算
    5. 反復して帰無分布を作成

    Parameters:
    -----------
    y : array-like
        応答変数 (PC scores)
    X_full_df : DataFrame
        フルデザイン行列（列名付き）
    interest_var_name : str
        検定対象の変数名（例: 'sex'）。この変数を含む全ての列を除外して帰無モデルを作成
    n_perm : int
        置換回数
    random_state : int
        乱数シード
    stratify_by : array-like, optional
        層化変数（例: subtype）。指定すると、各層内でのみ残差を置換
    one_sided : str, optional
        片側検定の方向。'greater' (obs > null) or 'less' (obs < null)
        Noneの場合は両側検定
    design_info : patsy.DesignInfo, optional
        デザイン行列の情報（堅牢な列特定のため）

    Returns:
    --------
    dict: observed (観測統計量), null_distribution (帰無分布), p_value,
          method_info (再現性のための情報)
    """
    import streamlit as st
    np.random.seed(random_state)

    # フルモデルをフィット
    model_full = sm.OLS(y, X_full_df).fit()

    # 検定対象の変数を含む列を特定（design_info使用で堅牢化）
    if design_info is not None:
        # patsy design_infoを使った堅牢な列特定
        interest_cols = []
        for term in design_info.terms:
            term_name = term.name()
            # 検定対象変数を含むターム（主効果 + 交互作用）を特定
            if term_name != 'Intercept' and interest_var_name in term_name:
                # このタームに対応する列を取得
                slice_obj = design_info.term_name_slices[term_name]
                cols = X_full_df.columns[slice_obj]
                interest_cols.extend(cols)
    else:
        # フォールバック: 文字列マッチング（後方互換性）
        interest_cols = [col for col in X_full_df.columns
                         if f'C({interest_var_name})' in col and col != 'Intercept']

    if not interest_cols:
        raise ValueError(f"検定対象変数 '{interest_var_name}' に対応する列が見つかりません")

    # 観測統計量（主効果の列の絶対t値を使用 - 堅牢化）
    main_effect_col = [col for col in interest_cols if ':' not in col]
    if not main_effect_col:
        raise ValueError(f"主効果の列が見つかりません: {interest_cols}")

    # 絶対値t統計量を使用（符号依存を避ける）
    obs_t = float(model_full.tvalues[main_effect_col[0]])
    obs_stat = np.abs(obs_t)

    # 帰無モデル（検定対象変数の全ての列を除外）
    X_reduced = X_full_df.drop(columns=interest_cols)

    if X_reduced.shape[1] == 1:  # 切片のみ
        residuals = y - np.mean(y)
        fitted_reduced = np.full_like(y, np.mean(y))
    else:
        model_reduced = sm.OLS(y, X_reduced).fit()
        residuals = np.array(model_reduced.resid)  # Convert to numpy array
        fitted_reduced = np.array(model_reduced.predict())  # Convert to numpy array


    # 帰無モデルと対立モデルの式を保存（レポート用）
    reduced_formula = "y ~ " + " + ".join([col for col in X_reduced.columns if col != 'Intercept'])
    if reduced_formula == "y ~ ":
        reduced_formula = "y ~ 1"  # 切片のみ
    full_formula = "y ~ " + " + ".join([col for col in X_full_df.columns if col != 'Intercept'])

    # 置換
    null_dist = []
    null_dist_abs = []  # 絶対値も保存
    for _ in range(n_perm):
        # 残差を置換（層化オプション）
        if stratify_by is not None:
            # 層化置換：各層内でのみ残差をシャッフル
            perm_residuals = np.zeros_like(residuals)
            for stratum in np.unique(stratify_by):
                stratum_mask = (stratify_by == stratum)
                stratum_indices = np.where(stratum_mask)[0]
                perm_indices = np.random.permutation(stratum_indices)
                perm_residuals[stratum_indices] = residuals[perm_indices]
        else:
            # 非層化置換：全体をシャッフル
            perm_idx = np.random.permutation(len(residuals))
            perm_residuals = residuals[perm_idx]

        y_perm = fitted_reduced + perm_residuals

        # フルモデルを再フィット
        model_perm = sm.OLS(y_perm, X_full_df).fit()
        perm_t = float(model_perm.tvalues[main_effect_col[0]])
        null_dist.append(perm_t)
        null_dist_abs.append(np.abs(perm_t))

    null_dist = np.array(null_dist)
    null_dist_abs = np.array(null_dist_abs)

    # p値の計算（片側/両側）+ 連続性補正
    if one_sided == 'greater':
        # 片側（大きい方向）: obs_tを使用
        p_value = (np.sum(null_dist >= obs_t) + 1) / (n_perm + 1)
    elif one_sided == 'less':
        # 片側（小さい方向）: obs_tを使用
        p_value = (np.sum(null_dist <= obs_t) + 1) / (n_perm + 1)
    else:
        # 両側検定: 絶対値で比較
        p_value = (np.sum(null_dist_abs >= obs_stat) + 1) / (n_perm + 1)

    return {
        'observed': obs_t,  # 元のt値（符号付き）を返す
        'observed_abs': obs_stat,  # 絶対値も返す
        'null_distribution': null_dist,  # 符号付きt値の分布
        'null_distribution_abs': null_dist_abs,  # 絶対値の分布
        'p_value': p_value,
        'stratified': stratify_by is not None,
        'one_sided': one_sided,
        'interest_cols': interest_cols,
        'main_effect_col': main_effect_col[0],
        'method_info': {  # 再現性のための情報
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
    """推定周辺平均 (Estimated Marginal Means) の計算

    共変量を特定の値（通常は平均）に固定して、
    主要な因子の水準ごとの予測平均値を計算します。

    例: 年齢を平均値に固定して、性別×サブタイプの組み合わせごとの
        PC scoreの予測平均を計算

    Parameters:
    -----------
    model : fitted statsmodels model
    df : DataFrame
        元データ
    factors : list
        EMM を計算する因子のリスト
    numeric_vars_at_mean : bool
        数値変数を平均値に固定するか

    Returns:
    --------
    DataFrame: 各組み合わせの予測平均、SE、95%CI
    """
    from itertools import product

    # 参照グリッドの作成
    ref_grid = {}

    for col in df.columns:
        if col in factors:
            ref_grid[col] = df[col].unique()
        elif pd.api.types.is_numeric_dtype(df[col]):
            if numeric_vars_at_mean:
                ref_grid[col] = [df[col].mean()]
            else:
                ref_grid[col] = df[col].unique()

    # 全組み合わせを生成
    grid_combos = list(product(*[ref_grid[col] for col in df.columns if col in ref_grid]))

    # 予測用データフレームを作成
    pred_df = pd.DataFrame(grid_combos, columns=[col for col in df.columns if col in ref_grid])

    # 予測
    predictions = model.get_prediction(pred_df)
    pred_summary = predictions.summary_frame(alpha=0.05)

    result_df = pred_df.copy()
    result_df['mean'] = pred_summary['mean']
    result_df['se'] = pred_summary['mean_se']
    result_df['ci_lower'] = pred_summary['mean_ci_lower']
    result_df['ci_upper'] = pred_summary['mean_ci_upper']

    return result_df

def plot_forest(coef_df, title="Forest Plot"):
    """係数のForest plot作成"""
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
    """モデル診断プロット

    4つの診断プロットを作成:
    1. Residuals vs Fitted: 非線形性、等分散性のチェック
    2. Q-Q plot: 残差の正規性チェック
    3. Scale-Location: 等分散性のチェック（標準化残差）
    4. Residual histogram: 残差の分布
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
    解析結果とグラフをzipファイルにまとめる

    Parameters:
    -----------
    pc_col : str
        PCカラム名（ファイル名に使用）

    Returns:
    --------
    bytes : zipファイルのバイナリデータ
    """
    zip_buffer = io.BytesIO()

    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base_name = f"pca_stats_{pc_col}_{timestamp}"

        # 0. Welch型t検定結果
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

        # 1. OLS係数表
        if 'ols_coef' in st.session_state:
            csv_data = st.session_state.ols_coef.to_csv()
            zip_file.writestr(f"{base_name}/tables/ols_coefficients.csv", csv_data)

        # 2. LMM係数表
        if 'lmm_coef' in st.session_state:
            csv_data = st.session_state.lmm_coef.to_csv()
            zip_file.writestr(f"{base_name}/tables/lmm_coefficients.csv", csv_data)

        # 3. Permutation test結果
        if 'perm_result' in st.session_state:
            perm_result = st.session_state.perm_result

            # 統計量サマリー
            perm_summary = pd.DataFrame({
                'observed_stat': [perm_result['observed']],
                'p_value': [perm_result['p_value']],
                'n_permutations': [len(perm_result['null_distribution'])],
                'method': [perm_result.get('method', 'unknown')],
                'stratified': [perm_result.get('stratified', False)]
            })
            zip_file.writestr(f"{base_name}/tables/permutation_summary.csv",
                            perm_summary.to_csv(index=False))

            # 帰無分布
            perm_df = pd.DataFrame({
                'permutation_stat': perm_result['null_distribution']
            })
            zip_file.writestr(f"{base_name}/tables/permutation_null_distribution.csv",
                            perm_df.to_csv(index=False))

        # 4. ANOVA表（OLS）
        if 'ols_anova' in st.session_state:
            csv_data = st.session_state.ols_anova.to_csv()
            zip_file.writestr(f"{base_name}/tables/ols_anova.csv", csv_data)

        # 5. EMM結果
        if 'ols_emm' in st.session_state:
            csv_data = st.session_state.ols_emm.to_csv(index=False)
            zip_file.writestr(f"{base_name}/tables/estimated_marginal_means.csv", csv_data)

        # 6. グラフの保存
        # Forest plot (OLS)
        if 'fig_forest_ols' in st.session_state:
            img_buffer = io.BytesIO()
            st.session_state.fig_forest_ols.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/ols_forest_plot.png", img_buffer.getvalue())

            # PDF版も保存
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

        # README作成
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
    st.title("⚙️ Analysis Options")

    st.markdown("### Select Statistical Method")
    analysis_method = st.radio(
        "解析手法:",
        ["OLS (通常最小二乗法)",
         "LMM (線形混合モデル)",
         "両方 (OLS + LMM)"],
        index=0,
        help="OLS: 固定効果のみ。LMM: ランダム効果（donor, batch等）を含む階層モデル"
    )

    if "OLS" in analysis_method:
        st.markdown("#### OLS オプション")
        se_type = st.radio(
            "標準誤差のタイプ:",
            ["Classical (通常)", "HC3 (頑健推定・推奨)", "Cluster-robust (クラスター頑健)"],
            index=1,
            help="HC3: 不均一分散に頑健。Cluster-robust: クラスター内相関を考慮"
        )

        if se_type == "Cluster-robust (クラスター頑健)":
            st.info("📌 クラスター変数はデータアップロード後に選択します")

        # ANOVA Type selection
        anova_type = st.radio(
            "ANOVA Type:",
            ["Type II (推奨)", "Type III"],
            index=0,
            help="Type II: 各変数の主効果を他の変数で調整して検定（交互作用なしの場合に推奨）\nType III: 全ての交互作用を含めて検定（交互作用ありの場合に推奨）"
        )

        # WLS option
        use_wls = st.checkbox(
            "WLS (加重最小二乗法) を使用",
            value=False,
            help="分散が不均一な場合に、各観測値に重みをつけて推定します。重みは残差の逆分散で自動計算されます"
        )

        if use_wls:
            with st.expander("⚠️ WLS使用上の重要な注意（必読）"):
                st.markdown("""
                ### WLS (Weighted Least Squares) の適切な使用

                **⚠️ 重要**: WLSは"重み"の決め方次第で推定が不安定になりやすく、**小サンプルでは推奨されません**

                ---

                ### 🚫 WLSを避けるべき状況（重要）

                #### **小サンプル設計では基本的にNG**

                例: subtype=4群、各群n=1-2程度の極小設計
                - ❌ **根拠ある重みが作りにくい**（分散推定が不安定）
                - ❌ **推定がブレやすい**（重みの不確実性が高い）
                - ❌ **過剰適合のリスク**（自由度が少ない）

                **代わりに使うべき方法**:
                - ✅ **Welch型 + Type II ANOVA** (主解析)
                - ✅ **OLS + HC3頑健SE** (係数・95%CI)
                - ✅ **Freedman-Lane置換検定** (条件付き帰無分布)

                これらは**重みを仮定せず**異分散に頑健で、小標本でも過度にブレにくいです。

                ---

                ### ✅ WLSが有効な状況（限定的）

                **十分な反復があり、群間分散差が明瞭に推定できる場合のみ**

                #### 必要条件（全て満たす必要あり）:

                1. **各群のサンプルサイズが十分**
                   - 最低でも各群n ≥ 10-15
                   - 分散を安定的に推定できる

                2. **群間で分散差が明瞭に存在**
                   - Levene検定やBartlett検定で検出可能
                   - 診断プロットで視覚的に明らか

                3. **重みの根拠が明確**
                   - 群内逆分散（σ²ᵢ）が事前情報から既知
                   - 技術的複製数など外部情報に基づく

                4. **モデル仮定が満たされる**
                   - 線形性、独立性
                   - 重み付き後の残差が正規分布

                #### この場合の推奨:

                - **主解析候補**: WLS（逆分散重み、重みを明示）
                - **HC3は通常不要**（WLSで等分散化できていれば）
                - 必要なら**"参考"**として併記（感度解析）

                ---

                ### 📊 使い分けの指針

                #### あなたの設計が小サンプル（4群×各n=1-2など）の場合:

                ```
                【主解析】
                ✅ Welch型 + Type II ANOVA
                   PC_k ~ C(sex) + C(subtype)
                   → 異分散頑健、小標本でも安定

                【併記】
                ✅ OLS + HC3頑健SE
                   → 係数推定・95%CI
                   → p値は頑健

                ✅ Freedman-Lane置換検定
                   → 条件付き帰無分布
                   → 分布の仮定なし

                【感度解析のみ】
                ⚠️ WLS
                   → 重み=群内逆分散など明示
                   → 主解析には使わない
                ```

                #### 十分な反復がある場合（各群n ≥ 15など）:

                ```
                【主解析候補】
                ✅ WLS（逆分散重み）
                   → 効率的推定
                   → 重みを明示

                【不要/参考】
                - HC3は通常不要
                  （WLSで等分散化済み）
                - 念のため感度解析として併記可
                ```

                ---

                ### ⚠️ WLSの危険性

                #### 小サンプルでWLSを使うと:

                1. **重みの推定誤差が大きい**
                   - 真の分散を正確に推定できない
                   - 誤った重みで推定が歪む

                2. **過剰適合（Overfitting）**
                   - 観測値へのフィッティングが強すぎる
                   - 一般化性能が低下

                3. **極端な重み**
                   - 外れ値に極端な重みが付く
                   - 数値的不安定性

                4. **解釈の困難さ**
                   - 重み付き平均効果の意味が不明確
                   - 再現性が低い

                ---

                ### 🔬 このアプリでの実装

                **方法**: 残差ベースの逆分散重み

                ```python
                1. 通常のOLSで残差を計算
                2. 重み = 1 / (残差²)
                3. 重みを正規化（平均=1）
                4. WLSで再推定
                ```

                **⚠️ この方法の問題**:
                - 残差ベースなので真の分散とは異なる
                - 小サンプルでは残差が不安定
                - **推奨**: 事前情報（技術的複製数、既知の分散）から重みを決定

                ---

                ### 💡 推奨事項（まとめ）

                #### 小サンプル（n < 15/群）:
                1. ❌ **WLSは使わない**
                2. ✅ **Welch型 + HC3を主解析に**
                3. ✅ **置換検定で補完**

                #### 大サンプル（n ≥ 15/群）:
                1. ✅ **分散差を可視化・検定**
                2. ✅ **重みを明示してWLS**
                3. ⚠️ **HC3との比較（感度解析）**

                #### 常に:
                - 📊 **診断プロット必須**（残差パターン確認）
                - 📝 **重みの根拠を明記**（論文報告時）
                - 🔄 **複数手法で結果の頑健性確認**

                ---

                ### 📚 参考: 異分散への対処法の比較

                | 方法 | サンプルサイズ要件 | 頑健性 | 推奨度（小n） | 推奨度（大n） |
                |------|------------------|--------|--------------|--------------|
                | **Welch t検定** | 小でもOK | 高い | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ |
                | **HC3頑健SE** | 小でもOK | 高い | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ |
                | **置換検定** | 小でもOK | 非常に高い | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
                | **WLS** | 大が必須 | 低い（小n時） | ⭐ | ⭐⭐⭐⭐⭐ |
                | **通常のOLS** | - | なし | ❌ | ⭐ |

                **結論**: あなたの設計が小サンプルなら、WLSは避けてWelch+HC3+置換検定を使いましょう。
                """)


    if "LMM" in analysis_method:
        st.markdown("#### LMM オプション")
        reml_method = st.radio(
            "推定方法:",
            ["REML (推奨)", "ML"],
            index=0,
            help="REML: 不偏な分散推定。ML: モデル比較用（AIC/BIC）"
        )

    st.markdown("---")
    st.markdown("### Welch型t検定")
    run_welch = st.checkbox("Welch型t検定を実行", value=False,
                           help="**主要変数が2群の場合のみ有効**\n\n"
                                "等分散を仮定しない2群比較（共変量調整なし）\n"
                                "単純な2群比較として参考値を算出")

    if run_welch:
        st.info("ℹ️ **Welch型t検定について:**\n\n"
                "- 主要変数の2群間で直接比較（共変量調整なし）\n"
                "- 等分散性を仮定しない頑健な検定\n"
                "- 主要変数が3群以上の場合は実行されません")

    st.markdown("---")
    st.markdown("### 置換検定 (Permutation Test)")
    run_permutation = st.checkbox("置換検定を実行", value=False,
                                  help="**検定対象: 主要変数1のみ**\n\n"
                                       "小サンプルや分布仮定が疑わしい場合に推奨\n"
                                       "主要変数2やその他の共変量は制御変数として扱われます")

    if run_permutation:
        st.info("ℹ️ **Permutation testの対象変数:** 主要変数1のみが検定されます。\n\n"
                "主要変数2やその他の変数は共変量として制御されます。")

        n_permutations = st.slider("置換回数:",
                                    min_value=1000, max_value=50000,
                                    value=10000, step=1000,
                                    help="多いほど正確だが時間がかかります")
        perm_method = st.radio(
            "置換戦略:",
            ["Freedman-Lane法 (推奨)", "Simple label permutation"],
            index=0,
            help="Freedman-Lane: 共変量を制御した厳密な検定\n"
                 "主要変数1の効果を検定し、他の変数は制御されます"
        )

        # 片側検定オプション
        perm_sided = st.radio(
            "検定の方向:",
            ["両側検定 (two-sided)", "片側検定: 主要変数1 > 基準 (greater)", "片側検定: 主要変数1 < 基準 (less)"],
            index=0,
            help="**両側検定（推奨）**: 効果の方向を問わず、差があるかを検定\n\n"
                 "**片側検定 (greater)**: 事前に「主要変数1が基準より大きい」と予想している場合のみ使用\n"
                 "例: Female > Male を仮説とする場合\n\n"
                 "**片側検定 (less)**: 事前に「主要変数1が基準より小さい」と予想している場合のみ使用\n\n"
                 "⚠️ 片側検定は事前仮説がある場合のみ使用してください"
        )

    st.markdown("---")
    st.markdown("### 可視化オプション")
    show_diagnostics = st.checkbox("診断プロットを表示", value=True,
                                   help="残差プロット、Q-Qプロット等")
    show_emm = st.checkbox("推定周辺平均 (EMM) を計算", value=True,
                          help="群別の調整済み平均値を算出")

# =============================================================================
# File Upload
# =============================================================================

st.markdown("## 1️⃣ データのアップロード")

# データ入力方法の選択
input_method = st.radio(
    "データ入力方法を選択:",
    ("PCA scoresファイル + メタデータを手動入力", "PCA scoresファイルをアップロード（メタデータ付き）"),
    help="pca.pyからダウンロードしたPCA scoresファイルを使用する場合、メタデータを手動で入力できます。"
)

if input_method == "PCA scoresファイル + メタデータを手動入力":
    # メタデータを手動入力
    st.markdown("""
    **PCA scoresファイルとメタデータを別々に入力**

    1. pca.pyからダウンロードしたPCA scoresファイル（サンプル名とPC1, PC2, ...を含む）をアップロード
    2. 下のフォームでメタデータを入力
    """)

    uploaded_file = st.file_uploader("PCA scoresファイルを選択", type=["tsv", "csv", "txt"],
                                     key="pca_scores_file")

else:  # メタデータ付きファイルをアップロード
    st.markdown("""
    **必要なデータ形式**: TSV または CSV ファイル

    - **1行 = 1サンプル** (pseudobulk レベル)
    - **PC score 列**: PC1, PC2, PC3, ... など
    - **カテゴリ変数**: sex, subtype, condition など
    - **任意: Donor/Subject ID** (反復測定がある場合)
    - **任意: 数値共変量**: age, batch, depth など
    """)

    uploaded_file = st.file_uploader("ファイルを選択", type=["tsv", "csv", "txt"])

if uploaded_file is not None:
    # ファイル読み込み
    try:
        df = pd.read_csv(uploaded_file, sep="\t", index_col=0)
    except:
        uploaded_file.seek(0)
        df = pd.read_csv(uploaded_file, sep=",", index_col=0)

    # インデックスをリセットして通常の列にする
    df = df.reset_index()

    # インデックス列名を適切に設定
    if df.columns[0] == 'index':
        df.rename(columns={'index': 'sample_id'}, inplace=True)

    st.success(f"✅ ファイル読み込み完了: {df.shape[0]} 行 × {df.shape[1]} 列")

    with st.expander("📋 データプレビュー (最初の10行)", expanded=True):
        st.dataframe(df.head(10))

    # session_stateの初期化
    if 'metadata_submitted' not in st.session_state:
        st.session_state.metadata_submitted = False
    if 'metadata_df' not in st.session_state:
        st.session_state.metadata_df = None
    if 'combined_df' not in st.session_state:
        st.session_state.combined_df = None

    # メタデータ手動入力の場合
    if input_method == "PCA scoresファイル + メタデータを手動入力":
        st.markdown("## 2️⃣ メタデータの入力")
        st.markdown("サンプルごとにメタデータを入力してください")

        # サンプル名を取得（最初の列または行名）
        if df.index.name or (df.index[0] != 0):  # インデックスがサンプル名の場合
            sample_names = df.index.tolist()
        else:  # 最初の列がサンプル名の場合
            sample_names = df.iloc[:, 0].tolist()

        # メタデータ列名の入力
        st.markdown("### メタデータ列の設定")
        col1, col2 = st.columns(2)

        with col1:
            meta_cols_input = st.text_area(
                "メタデータ列名（カンマ区切りまたは改行区切り）:",
                value="sex, subtype",
                help="例: sex, subtype, donor, age, batch"
            )

        # メタデータ列名のパース
        if ',' in meta_cols_input:
            meta_col_names = [x.strip() for x in meta_cols_input.split(',') if x.strip()]
        else:
            meta_col_names = [x.strip() for x in meta_cols_input.split('\n') if x.strip()]

        with col2:
            st.markdown("**入力する列:**")
            for col in meta_col_names:
                st.write(f"- {col}")

        # メタデータ入力フォーム
        st.markdown("### メタデータの入力")

        with st.form("metadata_input"):
            # 各サンプルのメタデータをデータエディタで入力
            meta_df = pd.DataFrame(index=sample_names, columns=meta_col_names)

            # DESeq2-LRTと同じアルゴリズムでデフォルト値を設定
            # 1. 末尾の共通要素を除去
            sample_names_str = [str(s) for s in sample_names]
            group_base = remove_common_suffix(sample_names_str)
            # 2. 末尾の数字を除去してグループ名を推定
            group_names = [remove_sample_num(s) for s in group_base]

            # 各メタデータ列のデフォルト値を設定
            for col in meta_col_names:
                col_lower = col.lower()

                if col_lower in ['sex', 'gender']:
                    # グループ名から性別を推測（male/femaleが含まれる場合）
                    def infer_sex(name):
                        name_lower = name.lower()
                        if 'male' in name_lower and 'female' not in name_lower:
                            return 'Male'
                        elif 'female' in name_lower or 'fem' in name_lower:
                            return 'Female'
                        else:
                            return 'Male'  # デフォルト
                    meta_df[col] = [infer_sex(g) for g in group_names]

                elif col_lower in ['donor', 'subject', 'patient', 'individual']:
                    # サンプル名からドナーIDを推測
                    # パターン1: sample_D01_A → D01
                    # パターン2: D01_typeA → D01
                    # パターン3: 連番で付与
                    donor_ids = []
                    for s in sample_names_str:
                        parts = s.split('_')
                        # D01のようなパターンを探す
                        donor_found = False
                        for part in parts:
                            if part.startswith('D') and len(part) >= 2 and part[1:].isdigit():
                                donor_ids.append(part)
                                donor_found = True
                                break
                        if not donor_found:
                            # グループ名を使用
                            donor_ids.append(group_names[len(donor_ids)])
                    meta_df[col] = donor_ids

                elif col_lower in ['subtype', 'celltype', 'cell_type', 'type', 'tissue']:
                    # グループ名をそのまま使用（これがサブタイプを表すことが多い）
                    meta_df[col] = group_names

                elif col_lower in ['batch', 'library', 'plate', 'run']:
                    # バッチ番号を推測
                    # パターン1: sample_B1_typeA → B1
                    # パターン2: 全サンプルで同じ場合は Batch1
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
                    # 数値型の場合は0で初期化
                    meta_df[col] = 0

                elif col_lower in ['condition', 'treatment', 'group', 'status']:
                    # グループ名を使用
                    meta_df[col] = group_names

                else:
                    # その他はグループ名を使用
                    meta_df[col] = group_names

            st.write("各サンプルのメタデータを編集してください:")
            st.info("💡 サンプル名から自動推測されたデフォルト値が入力されています。必要に応じて編集してください。")
            edited_meta_df = st.data_editor(meta_df, use_container_width=True)

            submitted = st.form_submit_button("メタデータを確定", type="primary")

        if submitted:
            # session_stateに保存
            st.session_state.metadata_submitted = True
            st.session_state.metadata_df = edited_meta_df.copy()
            st.session_state.meta_col_names = meta_col_names

            # PCA scoresとメタデータを結合
            if df.index.name or (df.index[0] != 0):
                df_combined = pd.concat([df, edited_meta_df], axis=1)
            else:
                # 最初の列がサンプル名の場合
                df_temp = df.set_index(df.columns[0])
                df_combined = pd.concat([df_temp, edited_meta_df], axis=1)
                df_combined.reset_index(inplace=True)
                df_combined.rename(columns={'index': df.columns[0]}, inplace=True)

            st.session_state.combined_df = df_combined.copy()

        # メタデータが確定されている場合、結合したデータフレームを使用
        if st.session_state.metadata_submitted and st.session_state.combined_df is not None:
            df = st.session_state.combined_df.copy()
            st.success("✅ メタデータ入力完了")

            st.write("メタデータ確認:")
            for col in st.session_state.meta_col_names:
                st.write(f"**{col}**: " + ', '.join(df[col].unique().astype(str)))

    # =============================================================================
    # Column Mapping
    # =============================================================================

    step_number = "3️⃣" if input_method == "PCA scoresファイル + メタデータを手動入力" else "2️⃣"
    st.markdown(f"## {step_number} 列のマッピング")
    st.markdown("各変数に対応する列を選択してください")

    with st.form("column_mapping"):
        col1, col2, col3 = st.columns(3)

        with col1:
            st.markdown("### 🔹 必須列")

            # Sample ID
            sample_col = st.selectbox("サンプルID列:", df.columns,
                                       index=0 if 'sample' in df.columns[0].lower() else 0)

            # PC score columns
            numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
            pc_cols = [col for col in numeric_cols if col.upper().startswith('PC')]

            if len(pc_cols) == 0:
                st.warning("'PC'で始まる列が見つかりません。全ての数値列を表示します。")
                pc_cols = numeric_cols

            pc_col = st.selectbox("解析するPC score列:", pc_cols,
                                  index=0,
                                  help="例: PC3, PC4 など。複数のPCを解析する場合は、1つずつ実行してください")

        with col2:
            st.markdown("### 🔹 主要変数")

            # カテゴリ列を取得
            cat_cols = df.select_dtypes(include=['object', 'category']).columns.tolist()

            # サンプルID列を除外（ユニーク値が多すぎる列は除外）
            cat_cols_filtered = []
            for col in cat_cols:
                n_unique = df[col].nunique()
                n_samples = len(df)
                unique_ratio = n_unique / n_samples

                # フィルタリング基準:
                # 1. サンプルID（全て一意: ratio >= 0.9）を除外
                # 2. ユニーク値が多すぎる（>50）列を除外
                # 3. ただし、ユニーク値が2以上あること（定数列を除外）
                if unique_ratio < 0.9 and n_unique >= 2 and n_unique <= 50:
                    cat_cols_filtered.append(col)

            if not cat_cols_filtered:
                st.warning("⚠️ 適切なカテゴリ変数が見つかりません。全てのカテゴリ列を表示します。")
                cat_cols_filtered = cat_cols

            # 主要変数 (例: sex)
            main_var = st.selectbox("🔴 主要変数 1 (必須):",
                                    cat_cols_filtered,
                                    index=0 if len(cat_cols_filtered) > 0 else None,
                                    help="検定したい主要な変数（例: sex, treatment など）。\n"
                                         "この変数はモデルの主効果として必ず含まれます。\n"
                                         "サンプルIDのような一意の列は自動的に除外されます。")

            # ブロッキング変数 (例: subtype)
            other_cats = [col for col in cat_cols_filtered if col != main_var]
            blocking_var = st.selectbox("🔴 主要変数 2 / ブロッキング変数 (任意):",
                                        ["なし"] + other_cats,
                                        index=1 if len(other_cats) > 0 else 0,
                                        help="**重要:** この変数もモデルの主効果として含まれます。\n\n"
                                             "例: sex を主要変数1、subtype を主要変数2として選択すると、\n"
                                             "モデル式は `PC ~ C(sex) + C(subtype)` となります。\n\n"
                                             "**使用例:**\n"
                                             "- 2つの主要変数を同時に検定したい場合（sex と subtype など）\n"
                                             "- ブロッキングデザイン（donor, batch など）で層別化したい場合\n\n"
                                             "交互作用項を追加すると、2つの変数の相互作用も検定できます。")
            blocking_var = None if blocking_var == "なし" else blocking_var

        with col3:
            st.markdown("### 🔹 階層構造")

            # Donor ID for random effects
            donor_col = st.selectbox("Donor/Subject ID (ランダム効果用):",
                                     ["なし", "サンプルIDと同じ"] + df.columns.tolist(),
                                     index=0,
                                     help="反復測定がある場合に指定。LMMで (1|donor) のランダム効果として使用")

            if donor_col == "なし":
                donor_col = None
            elif donor_col == "サンプルIDと同じ":
                donor_col = sample_col

            # Additional random effect
            additional_re = st.selectbox("追加のランダム効果 (例: batch, library):",
                                         ["なし"] + df.columns.tolist(),
                                         index=0,
                                         help="技術的なバッチ効果などをランダム効果として扱う場合に使用")
            additional_re = None if additional_re == "なし" else additional_re

        # Additional covariates
        st.markdown("### 🔹 追加の共変量（固定効果）")

        # 既に使用した列を除外
        used_cols = [sample_col, pc_col, main_var, blocking_var, donor_col, additional_re]
        remaining_cols = [col for col in df.columns if col not in used_cols and col is not None]

        additional_covars = st.multiselect(
            "追加の共変量を選択 (数値またはカテゴリ):",
            remaining_cols,
            help="年齢、バッチ、深度などの調整したい変数。これらは固定効果としてモデルに含まれます"
        )

        # Interaction terms
        st.markdown("### 🔹 交互作用項")
        st.markdown("_交互作用は、2つの変数の効果が独立でない場合に使用します（例: 性差がサブタイプによって異なる）_")

        potential_interactions = []
        if blocking_var:
            potential_interactions.append(f"{main_var} : {blocking_var}")

        for covar in additional_covars:
            if df[covar].dtype in ['object', 'category']:
                potential_interactions.append(f"{main_var} : {covar}")

        interactions = st.multiselect(
            "交互作用項を選択 (任意):",
            potential_interactions,
            help="例: sex:subtype は、性差がサブタイプ間で異なるかを検定します。データが十分にある場合のみ使用を推奨"
        )

        # フォームの送信ボタン
        mapping_submitted = st.form_submit_button("列のマッピングを確定", type="primary")

    # =============================================================================
    # Run Analysis Button
    # =============================================================================

    # 列のマッピングが確定された場合のみ解析ボタンを表示
    if mapping_submitted or 'mapping_confirmed' in st.session_state:
        if mapping_submitted:
            # 列マッピングが変更された場合、以前の解析結果をすべてクリア
            st.session_state.mapping_confirmed = True
            st.session_state.sample_col = sample_col
            st.session_state.pc_col = pc_col
            st.session_state.main_var = main_var
            st.session_state.blocking_var = blocking_var
            st.session_state.donor_col = donor_col
            st.session_state.additional_re = additional_re
            st.session_state.additional_covars = additional_covars
            st.session_state.interactions = interactions

            # 以前の解析結果をクリア
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

        st.success("✅ 列のマッピング完了")

        # session_stateから値を取得
        sample_col = st.session_state.sample_col
        pc_col = st.session_state.pc_col
        main_var = st.session_state.main_var
        blocking_var = st.session_state.blocking_var
        donor_col = st.session_state.donor_col
        additional_re = st.session_state.additional_re
        additional_covars = st.session_state.additional_covars
        interactions = st.session_state.interactions

        if st.button("🚀 解析を実行", type="primary"):

            step_number_check = "4️⃣" if input_method == "PCA scoresファイル + メタデータを手動入力" else "3️⃣"
            st.markdown(f"## {step_number_check} データ品質チェック")

            # 解析用データフレームの準備
            analysis_cols = [sample_col, pc_col, main_var]
            if blocking_var:
                analysis_cols.append(blocking_var)
            if donor_col and donor_col != sample_col:
                analysis_cols.append(donor_col)
            if additional_re:
                analysis_cols.append(additional_re)

            # additional_covarsが空でない場合のみ追加
            if additional_covars:
                # リストであることを確認してから追加
                if isinstance(additional_covars, list):
                    analysis_cols.extend(additional_covars)
                else:
                    analysis_cols.append(additional_covars)

            # Noneを除外し、重複も除去
            analysis_cols = [col for col in analysis_cols if col is not None]
            analysis_cols = list(dict.fromkeys(analysis_cols))  # 重複を保持順で除去

            # デバッグ情報
            st.write("**選択された列:**", analysis_cols)

            analysis_df = df[analysis_cols].copy()

            # 欠損値の除去
            n_before = len(analysis_df)
            analysis_df = analysis_df.dropna()
            n_after = len(analysis_df)

            if n_before > n_after:
                st.warning(f"⚠️ 欠損値を含む {n_before - n_after} 行を削除しました")

            # カテゴリ変数の標準化
            for col in analysis_df.columns:
                # colが文字列であることを確認
                if isinstance(col, str):
                    try:
                        if analysis_df[col].dtype == 'object':
                            analysis_df[col] = analysis_df[col].astype('category')
                    except (KeyError, AttributeError):
                        # 列が存在しない、または型がない場合はスキップ
                        continue

            # サンプルサイズの確認（欠測情報を含む）
            st.markdown("### 📊 サンプルサイズの要約")

            # 欠測チェック
            original_n = len(df)  # 元のデータ数
            analysis_n = len(analysis_df)  # 解析に使用するデータ数
            excluded_n = original_n - analysis_n

            if excluded_n > 0:
                st.warning(f"⚠️ **欠測値除外**: {excluded_n} 件のサンプルが欠測値により除外されました（元データ: {original_n} → 解析使用: {analysis_n}）")

            col1, col2 = st.columns(2)

            with col1:
                st.metric("解析に使用したサンプル数", len(analysis_df))

                # 主要変数の分布
                st.write(f"**{main_var} の分布:**")
                main_var_counts = analysis_df[main_var].value_counts()
                st.dataframe(main_var_counts)

                # 不均衡の警告
                if len(main_var_counts) > 1:
                    max_count = main_var_counts.max()
                    min_count = main_var_counts.min()
                    if max_count / min_count > 3:
                        st.warning(f"⚠️ **群の不均衡**: 最大群({max_count})と最小群({min_count})で3倍以上の差があります")

            with col2:
                if blocking_var:
                    st.write(f"**クロス集計表: {main_var} × {blocking_var}**")
                    cross_tab = pd.crosstab(analysis_df[main_var], analysis_df[blocking_var])
                    st.dataframe(cross_tab)

                    # 小サンプルサイズの警告
                    min_cell_size = cross_tab.min().min()
                    if min_cell_size < 3:
                        st.warning(f"⚠️ 一部のセルのサンプルサイズが非常に小さい（最小: {min_cell_size}）です。推定が不安定になる可能性があります。")

                if donor_col and donor_col != sample_col:
                    n_donors = analysis_df[donor_col].nunique()
                    st.metric("ドナー数", n_donors)

                    if n_donors < 5:
                        st.warning("⚠️ ドナー数が5未満です。ランダム効果の分散推定が不安定になる可能性があります。")

            # 共線性チェック
            st.markdown("### 🔍 共線性診断")

            formula_terms = [main_var]
            if blocking_var:
                formula_terms.append(blocking_var)
            formula_terms.extend(additional_covars)

            collinearity_issues = detect_collinearity(analysis_df, formula_terms)

            if collinearity_issues:
                st.error("**完全共線性が検出されました:**")
                for issue in collinearity_issues:
                    st.markdown(issue)
                st.warning("""
                ⚠️ **完全に共役している変数は、効果を分離して推定できません。**

                例: あるsubtypeが全てFemaleの場合、sex効果とsubtype効果は識別不能です。

                **対処法:**
                - いずれかの変数を削除
                - より均衡の取れたデータを収集
                - 交互作用項を使用せず、記述的な比較に留める
                """)
            else:
                st.success("✅ 完全共線性は検出されませんでした")

            # VIF for numeric variables
            numeric_covars = [col for col in additional_covars
                             if pd.api.types.is_numeric_dtype(analysis_df[col])]

            if len(numeric_covars) > 1:
                st.markdown("**数値共変量の多重共線性 (VIF):**")
                st.caption("VIF > 10: 強い多重共線性あり。VIF > 5: 中程度の多重共線性あり。")

                try:
                    vif_df = calculate_vif(analysis_df, numeric_covars)
                    vif_df_display = vif_df.copy(); vif_df_display['VIF'] = vif_df_display['VIF'].round(2); st.dataframe(vif_df_display, use_container_width=True)

                    if (vif_df['VIF'] > 10).any():
                        st.warning("⚠️ VIF > 10 の変数があります。強い多重共線性が示唆されます。")
                    elif (vif_df['VIF'] > 5).any():
                        st.info("ℹ️ VIF > 5 の変数があります。中程度の多重共線性が示唆されます。")
                except Exception as e:
                    st.info(f"VIF計算をスキップしました: {str(e)}")

            # =============================================================================
            # Build Formula
            # =============================================================================

            st.markdown("## 4️⃣ モデル式の構築")

            # デバッグ情報
            with st.expander("🔍 モデル構築情報（選択された変数の確認）", expanded=True):
                st.write(f"**主要変数 1 (main_var):** {main_var}")
                st.write(f"**主要変数 2 / ブロッキング変数 (blocking_var):** {blocking_var}")
                st.write(f"**追加共変量 (additional_covars):** {additional_covars}")
                st.write(f"**交互作用項 (interactions):** {interactions}")
                st.info("💡 モデル式に含めたい変数は、上の '主要変数 1' または '主要変数 2' で選択してください。\n\n"
                        "'追加の共変量' は調整変数として使用され、主な検定対象ではありません。")

            # 固定効果の式
            fixed_terms = [f"C({main_var})"]
            if blocking_var:
                fixed_terms.append(f"C({blocking_var})")

            for col in additional_covars:
                if analysis_df[col].dtype in ['object', 'category']:
                    fixed_terms.append(f"C({col})")
                else:
                    fixed_terms.append(col)

            # 交互作用項の追加
            if interactions:
                for interaction in interactions:
                    vars_in_interaction = [v.strip() for v in interaction.split(":")]
                    interaction_term = " * ".join([
                        f"C({v})" if analysis_df[v].dtype in ['object', 'category'] else v
                        for v in vars_in_interaction
                    ])
                    # 既に含まれていない場合のみ追加
                    if interaction_term not in fixed_terms:
                        fixed_terms.append(interaction_term)

            formula = f"{pc_col} ~ " + " + ".join(fixed_terms)

            # モデルに含まれる変数を明示的に表示
            st.markdown("**📊 モデルに含まれる変数:**")
            col_a, col_b = st.columns(2)
            with col_a:
                st.success(f"✅ 主要変数 1: **{main_var}**")
                if blocking_var:
                    st.success(f"✅ 主要変数 2: **{blocking_var}**")
                else:
                    st.warning("⚠️ 主要変数 2: 選択なし")
            with col_b:
                if additional_covars:
                    st.info(f"📝 調整変数: {', '.join(additional_covars)}")
                if interactions:
                    st.info(f"🔗 交互作用: {', '.join(interactions)}")

            st.markdown("**モデル式:**")
            st.code(formula, language="r")

            if donor_col:
                st.markdown(f"**ランダム効果 (LMMのみ):** `(1|{donor_col})`")

            # =============================================================================
            # Run OLS
            # =============================================================================

            if "OLS" in analysis_method:
                st.markdown("---")
                st.markdown("## 5️⃣ OLS 解析結果")

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
                        st.info(f"ℹ️ WLS（加重最小二乗法）を使用しています。重みは残差の逆分散で計算されました")
                    else:
                        # モデルのフィット (OLS)
                        ols_model = smf.ols(formula, data=analysis_df).fit()

                    # 適切な共分散行列を取得
                    if se_type == "HC3 (頑健推定・推奨)":
                        ols_model = ols_model.get_robustcov_results(cov_type='HC3')
                        st.info("ℹ️ HC3頑健標準誤差を使用しています（不均一分散に対応）")
                    elif se_type == "Cluster-robust (クラスター頑健)":
                        cluster_var = donor_col if donor_col else additional_re
                        if cluster_var:
                            ols_model = ols_model.get_robustcov_results(
                                cov_type='cluster',
                                groups=analysis_df[cluster_var]
                            )
                            st.info(f"ℹ️ クラスター頑健標準誤差を使用しています（クラスター変数: {cluster_var}）")
                        else:
                            st.warning("クラスター変数が指定されていません。通常のHC3を使用します。")
                            ols_model = ols_model.get_robustcov_results(cov_type='HC3')

                    # モデル式と参照水準の明示
                    st.markdown("### 📋 係数表")
                    st.caption(f"**モデル式**: `{formula}`")

                    # 参照水準の抽出と表示
                    reference_levels = {}
                    for col in analysis_df.columns:
                        if analysis_df[col].dtype == 'object' or analysis_df[col].dtype.name == 'category':
                            if col in formula:
                                ref_level = analysis_df[col].iloc[0] if len(analysis_df) > 0 else "N/A"
                                # Patsyはアルファベット順の最初を基準にする
                                unique_vals = sorted(analysis_df[col].unique())
                                if len(unique_vals) > 0:
                                    reference_levels[col] = unique_vals[0]

                    if reference_levels:
                        ref_info = ", ".join([f"`{var}` の基準={ref}" for var, ref in reference_levels.items()])
                        st.caption(f"**参照水準**: {ref_info}")
                    conf_int = ols_model.conf_int()

                    # conf_intがDataFrameかndarrayかを判定
                    if isinstance(conf_int, pd.DataFrame):
                        ci_lower = conf_int.iloc[:, 0]
                        ci_upper = conf_int.iloc[:, 1]
                    else:
                        # numpy配列の場合
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

                    # 有意性の記号を追加
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

                    # 数値を丸めて表示用のDataFrameを作成
                    display_table = coef_table.copy()
                    display_table['Coef'] = display_table['Coef'].round(4)
                    display_table['Std_Error'] = display_table['Std_Error'].round(4)
                    display_table['t_value'] = display_table['t_value'].round(3)
                    display_table['p_value'] = display_table['p_value'].apply(lambda x: f"{x:.4g}")
                    display_table['CI_lower'] = display_table['CI_lower'].round(4)
                    display_table['CI_upper'] = display_table['CI_upper'].round(4)

                    st.dataframe(display_table, use_container_width=True)

                    st.caption("有意水準: *** p<0.001, ** p<0.01, * p<0.05, . p<0.1")

                    # 係数表の解釈ガイド
                    with st.expander("📘 係数表の見方（p_value列の意味）"):
                        st.markdown("""
                        ### p_value列の意味

                        **各係数の有意性を示します（他の変数を調整した上で）**

                        - ✅ **共変量調整あり**: 全ての変数を同時に考慮した結果
                        - この p値は、他の変数の影響を取り除いた「純粋な効果」を検定しています

                        ---

                        ### 表記の意味
                        """)

                        # 実際のモデルに含まれる変数から説明を生成
                        example_vars = []
                        for idx in display_table.index:
                            idx_str = str(idx)  # Ensure idx is a string
                            if idx_str != 'Intercept' and 'C(' in idx_str:
                                example_vars.append(idx_str)

                        if example_vars:
                            st.markdown("**このモデルの例:**")
                            # 最初の1-2個の変数を例として使用
                            for var_name in example_vars[:2]:
                                # C(sex)[T.M] のような形式から分解
                                if '[T.' in var_name:
                                    var_part = var_name.split('[')[0]  # C(sex)
                                    level_part = var_name.split('[T.')[1].rstrip(']')  # M
                                    base_var = var_part.replace('C(', '').replace(')', '')  # sex

                                    st.markdown(f"""
                                    - **`{var_name}`**
                                      - `C({base_var})`: {base_var}をカテゴリカル変数として扱う
                                      - `[T.{level_part}]`: Treatment coding で {level_part}水準を表す
                                      - **意味**: 「{level_part}は基準カテゴリと比べてどれだけ違うか」
                                      - **p値の解釈**: {level_part}と基準カテゴリに有意な差があるか（**他の変数を調整済み**）
                                    """)

                        st.markdown("""
                        ---

                        ### 一般的な表記

                        | 表記 | 意味 |
                        |-----|------|
                        | `C(変数名)` | カテゴリカル変数として扱う |
                        | `[T.水準名]` | Treatment codingで指定水準を表す |
                        | `C(変数名)[T.水準名]` | 「この水準は基準カテゴリと比べてどれだけ違うか」 |

                        ### 基準カテゴリ（Reference level）

                        - デフォルト: **アルファベット順で最初の値**
                        - 基準カテゴリの行は係数表に表示されません（係数=0として扱われる）
                        - 他の水準は「基準との差」として解釈されます

                        ### Welch t検定との違い

                        | 検定方法 | 共変量調整 | 表示場所 |
                        |---------|----------|---------|
                        | Welch t検定 | ❌ なし（単純2群比較） | 情報ボックス |
                        | 係数表のp_value | ✅ あり（他変数を考慮） | この表 |

                        **推奨**: 正確な統計解析には、この係数表のp値を使用してください。
                        """)

                    # モデル適合度
                    st.markdown("### 📈 モデル適合度")
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
                        if anova_type == "Type II (推奨)":
                            st.markdown("### 📊 Type II ANOVA")
                            st.caption("各変数の全体的な有意性を検定（他の変数を含むモデルでの寄与）")
                            anova_table = anova_lm(ols_model, typ=2)
                        else:  # Type III
                            st.markdown("### 📊 Type III ANOVA")
                            st.caption("全ての交互作用を含めた各変数の有意性を検定")
                            anova_table = anova_lm(ols_model, typ=3)

                        # 重要な注意事項を表示（SE typeに応じて）
                        if se_type == "HC3 (頑健推定・推奨)":
                            st.warning("⚠️ **ANOVA F検定は等分散性と正規性を仮定しています**。上記の係数表（HC3頑健SE）とはp値の意味が異なります。不均一分散が懸念される場合は係数表のp値を優先してください。")
                        elif se_type == "Cluster-robust (クラスター頑健)":
                            st.warning("⚠️ **ANOVA F検定は等分散性と正規性を仮定しています**。上記の係数表（Cluster頑健SE）とはp値の意味が異なります。クラスター内相関や不均一分散が懸念される場合は係数表のp値を優先してください。")
                        elif se_type == "Classical (通常)":
                            st.info("ℹ️ **ANOVA F検定と係数表は同じ仮定**（等分散性・正規性）を使用しています。不均一分散が懸念される場合は、HC3頑健SEの使用を検討してください。")

                        # 数値を丸めて表示
                        anova_display = anova_table.copy()
                        if 'sum_sq' in anova_display.columns:
                            anova_display['sum_sq'] = anova_display['sum_sq'].round(4)
                        if 'F' in anova_display.columns:
                            anova_display['F'] = anova_display['F'].round(3)
                        if 'PR(>F)' in anova_display.columns:
                            anova_display['PR(>F)'] = anova_display['PR(>F)'].apply(lambda x: f"{x:.4g}" if pd.notna(x) else x)

                        st.dataframe(anova_display, use_container_width=True)
                        st.session_state.ols_anova = anova_table  # zipダウンロード用に保存

                        # ANOVA表の解釈ガイド
                        with st.expander("📘 ANOVA表の見方（PR(>F)列の意味）"):
                            st.markdown("""
                            ### PR(>F)列の意味

                            **各変数全体の有意性を示します（他の変数を調整した上で）**

                            - ✅ **共変量調整あり**: 全ての変数を同時に考慮した結果
                            - この p値は、他の変数の影響を取り除いた上で、その変数全体が有意に寄与しているかを検定しています

                            ---

                            ### 係数表との違い

                            | 表 | 検定対象 | 共変量調整 |
                            |----|---------|----------|
                            | **係数表** | 各水準ごと（例: Male vs Female） | ✅ あり |
                            | **ANOVA表** | 変数全体（例: sex全体の効果） | ✅ あり |
                            """)

                            # 実際のモデルに含まれる変数から説明を生成
                            anova_vars = []
                            for idx in anova_display.index:
                                idx_str = str(idx)  # Ensure idx is a string
                                if idx_str != 'Residual' and 'C(' in idx_str:
                                    anova_vars.append(idx_str)

                            if anova_vars:
                                st.markdown("**このモデルの例:**")
                                for var_name in anova_vars[:2]:
                                    # C(sex) のような形式から変数名を抽出
                                    base_var = var_name.replace('C(', '').replace(')', '')

                                    st.markdown(f"""
                                    - **`{var_name}`** の PR(>F)
                                      - **意味**: {base_var}という変数全体が、PC値に有意な影響を与えているか
                                      - **検定内容**: 全ての{base_var}水準の係数が同時にゼロかどうか
                                      - **調整**: 他の変数（主要変数2、共変量など）の影響を除いた上での検定
                                    """)

                            st.markdown("""
                            ---

                            ### 具体例で理解する

                            **モデル**: `PC4 ~ C(sex) + C(subtype)`

                            | 変数 | PR(>F) | 意味 |
                            |------|--------|------|
                            | C(sex) | 0.0234 | subtype調整後も、sexは有意にPC4に影響する |
                            | C(subtype) | 0.0567 | sex調整後も、subtypeは有意にPC4に影響する |

                            ### Welch t検定との違い

                            | 検定方法 | 共変量調整 | 検定対象 | 表示場所 |
                            |---------|----------|---------|---------|
                            | Welch t検定 | ❌ なし | 主要変数1のみ（2群比較） | 情報ボックス |
                            | ANOVA表 PR(>F) | ✅ あり | 各変数全体 | この表 |

                            ### 注意事項

                            ⚠️ **ANOVA F検定は等分散性と正規性を仮定しています**

                            不均一分散が懸念される場合は、係数表のp値（HC3/Cluster頑健SE使用時）を優先してください。

                            ### Type II vs Type III

                            - **Type II ANOVA（推奨）**: 交互作用を含まないモデルで各変数の主効果を検定
                            - **Type III ANOVA**: 全ての交互作用を含めた各変数の効果を検定

                            **推奨**: 正確な変数全体の有意性には、このANOVA表のPR(>F)を使用してください。
                            """)

                    except Exception as e:
                        st.warning(f"ANOVA表を計算できませんでした: {str(e)}")

                    # Forest plot
                    st.markdown("### 🌲 Forest Plot (係数の可視化)")
                    st.caption("各係数の推定値と95%信頼区間。ゼロと重ならない場合は有意")

                    # 切片を除外
                    coef_for_plot = coef_table.iloc[1:].copy()
                    if len(coef_for_plot) > 0:
                        fig_forest = plot_forest(coef_for_plot, title=f"OLS Coefficient Estimates ({se_type})")
                        st.pyplot(fig_forest)
                        st.session_state.fig_forest_ols = fig_forest  # zipダウンロード用に保存
                    else:
                        st.info("プロット可能な係数がありません（切片のみ）")

                    # 診断プロット
                    if show_diagnostics:
                        st.markdown("### 🔬 診断プロット")
                        st.caption("""
                        - **Residuals vs Fitted**: パターンがなければ線形性と等分散性が満たされている
                        - **Normal Q-Q**: 直線上に乗れば正規性が満たされている
                        - **Scale-Location**: 水平な帯状であれば等分散性が満たされている
                        - **Residual Distribution**: 正規分布に近いかを確認
                        """)
                        fig_diag = plot_diagnostic(ols_model, title="OLS Model Diagnostics")
                        st.pyplot(fig_diag)
                        st.session_state.fig_diagnostic = fig_diag  # zipダウンロード用に保存

                    # EMM
                    if show_emm and blocking_var:
                        st.markdown("### 📊 推定周辺平均 (Estimated Marginal Means)")
                        st.caption("共変量を調整した後の、各群の予測平均値")

                        try:
                            emm_df = calculate_emm(ols_model, analysis_df,
                                                  [main_var, blocking_var] if blocking_var else [main_var])

                            # 数値を丸めて表示
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
                            st.session_state.fig_emm = fig  # zipダウンロード用に保存
                            st.session_state.ols_emm = emm_df  # データも保存

                        except Exception as e:
                            st.warning(f"EMMの計算中にエラーが発生しました: {str(e)}")

                    # 結果を保存
                    st.session_state.ols_model = ols_model
                    st.session_state.ols_coef = coef_table

                except Exception as e:
                    st.error(f"❌ OLS解析に失敗しました")
                    st.exception(e)

            # =============================================================================
            # Run LMM
            # =============================================================================

            if "LMM" in analysis_method and donor_col:
                st.markdown("---")
                st.markdown("## 6️⃣ 線形混合モデル (LMM) 解析結果")

                try:
                    # ランダム効果の式を準備
                    re_formula = "1"  # ランダム切片
                    groups = analysis_df[donor_col]

                    # LMMのフィット
                    lmm_model = smf.mixedlm(formula, data=analysis_df,
                                            groups=groups,
                                            re_formula=re_formula)

                    lmm_result = lmm_model.fit(reml=(reml_method == "REML (推奨)"))

                    st.info(f"ℹ️ 推定方法: {reml_method}")

                    # 係数表
                    st.markdown("### 📋 固定効果の係数表")
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

                    # 数値を丸めて表示
                    lmm_display = coef_table_lmm.copy()
                    lmm_display['Coef'] = lmm_display['Coef'].round(4)
                    lmm_display['Std_Error'] = lmm_display['Std_Error'].round(4)
                    lmm_display['z_value'] = lmm_display['z_value'].round(3)
                    lmm_display['p_value'] = lmm_display['p_value'].apply(lambda x: f"{x:.4g}")
                    lmm_display['CI_lower'] = lmm_display['CI_lower'].round(4)
                    lmm_display['CI_upper'] = lmm_display['CI_upper'].round(4)

                    st.dataframe(lmm_display, use_container_width=True)

                    st.caption("有意水準: *** p<0.001, ** p<0.01, * p<0.05, . p<0.1")

                    # 分散成分
                    st.markdown("### 📊 分散成分 (Variance Components)")
                    st.caption("データの変動がどこから来ているかを示します")

                    var_random = lmm_result.cov_re.values[0][0]
                    var_residual = lmm_result.scale
                    var_total = var_random + var_residual

                    var_comp = pd.DataFrame({
                        'Component': [f'ランダム効果 ({donor_col})', '残差 (個体内変動)'],
                        'Variance': [var_random, var_residual],
                        'Std_Dev': [np.sqrt(var_random), np.sqrt(var_residual)],
                        'Proportion': [var_random / var_total, var_residual / var_total]
                    })

                    # 数値を丸めて表示
                    var_comp_display = var_comp.copy()
                    var_comp_display['Variance'] = var_comp_display['Variance'].round(4)
                    var_comp_display['Std_Dev'] = var_comp_display['Std_Dev'].round(4)
                    var_comp_display['Proportion'] = (var_comp_display['Proportion'] * 100).round(2).astype(str) + '%'

                    st.dataframe(var_comp_display, use_container_width=True)

                    # ICC
                    icc = var_random / var_total
                    st.metric("級内相関係数 (ICC)", f"{icc:.4f}",
                             help="同一ドナー内のサンプル間の相関。0に近いとランダム効果が小さい、1に近いとドナー間の変動が大きい")

                    if icc < 0.05:
                        st.info("ℹ️ ICC < 0.05: ランダム効果が非常に小さいです。OLSでも十分かもしれません。")
                    elif icc > 0.5:
                        st.success("✅ ICC > 0.5: ランダム効果が大きいです。LMMの使用が適切です。")

                    # モデル適合度
                    st.markdown("### 📈 モデル適合度")
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("対数尤度", f"{lmm_result.llf:.2f}")
                    with col2:
                        st.metric("AIC", f"{lmm_result.aic:.2f}")
                    with col3:
                        st.metric("BIC", f"{lmm_result.bic:.2f}")

                    # Forest plot
                    st.markdown("### 🌲 Forest Plot (固定効果)")
                    coef_for_plot_lmm = coef_table_lmm.iloc[1:].copy()
                    if len(coef_for_plot_lmm) > 0:
                        fig_forest_lmm = plot_forest(coef_for_plot_lmm,
                                                     title="LMM Fixed Effects Estimates")
                        st.pyplot(fig_forest_lmm)
                        st.session_state.fig_forest_lmm = fig_forest_lmm  # zipダウンロード用に保存

                    # 結果を保存
                    st.session_state.lmm_result = lmm_result
                    st.session_state.lmm_coef = coef_table_lmm

                except Exception as e:
                    st.error(f"❌ LMM解析に失敗しました")
                    st.exception(e)
                    st.info("ヒント: ドナー数が少ない場合や、データの変動が小さい場合、LMMの収束に失敗することがあります。OLSを試してください。")

            elif "LMM" in analysis_method and not donor_col:
                st.warning("⚠️ LMMを実行するには、Donor/Subject ID を指定してください。")

            # =============================================================================
            # Welch's t-test
            # =============================================================================

            if run_welch:
                st.markdown("---")
                st.markdown("## 7️⃣ Welch型t検定結果")
                st.caption("等分散を仮定しない2群比較（共変量調整なし）")
                st.warning(f"🎯 **検定対象変数:** 主要変数1 = **{main_var}**")

                # 主要変数の水準数を確認
                unique_main_levels = analysis_df[main_var].nunique()
                if unique_main_levels == 2:
                    from scipy import stats
                    groups = analysis_df[main_var].unique()
                    group1_data = analysis_df[analysis_df[main_var] == groups[0]][pc_col]
                    group2_data = analysis_df[analysis_df[main_var] == groups[1]][pc_col]
                    welch_t, welch_p = stats.ttest_ind(group1_data, group2_data, equal_var=False)

                    # 記述統計
                    group1_mean = group1_data.mean()
                    group1_std = group1_data.std()
                    group1_n = len(group1_data)
                    group2_mean = group2_data.mean()
                    group2_std = group2_data.std()
                    group2_n = len(group2_data)

                    # 効果量 (Cohen's d, Hedges' g)
                    pooled_std = np.sqrt(((group1_n - 1) * group1_std**2 + (group2_n - 1) * group2_std**2) / (group1_n + group2_n - 2))
                    cohens_d = (group1_mean - group2_mean) / pooled_std
                    correction_factor = 1 - (3 / (4 * (group1_n + group2_n) - 9))
                    hedges_g = cohens_d * correction_factor

                    # 結果の表示
                    st.success("✅ Welch型t検定を実行しました")

                    col1, col2 = st.columns(2)

                    with col1:
                        st.markdown("### 📊 記述統計")
                        desc_df = pd.DataFrame({
                            '群': [groups[0], groups[1]],
                            'n': [group1_n, group2_n],
                            '平均': [f"{group1_mean:.4f}", f"{group2_mean:.4f}"],
                            '標準偏差': [f"{group1_std:.4f}", f"{group2_std:.4f}"]
                        })
                        st.dataframe(desc_df, use_container_width=True)

                    with col2:
                        st.markdown("### 📈 検定結果")
                        st.metric("t値", f"{welch_t:.4f}")
                        st.metric("p値", f"{welch_p:.4g}")
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

                    st.markdown("### 📏 効果量")
                    effect_col1, effect_col2 = st.columns(2)
                    with effect_col1:
                        st.metric("Cohen's d", f"{cohens_d:.4f}")
                    with effect_col2:
                        st.metric("Hedges' g (補正版)", f"{hedges_g:.4f}")

                    st.caption("効果量の目安: |d| < 0.2 (小)、0.2-0.5 (中)、0.5-0.8 (大)、≥ 0.8 (非常に大)")

                    # 重要な注意事項
                    st.warning(
                        f"⚠️ **重要な注意**\n\n"
                        f"この検定は **単純2群比較** です：\n"
                        f"- 他の変数（主要変数2、追加共変量など）の影響は **調整されていません**\n"
                        f"- 共変量調整が必要な場合は、OLS/LMM の「係数表」または「ANOVA表」をご確認ください\n\n"
                        f"💡 **使い分け**：\n"
                        f"- Welch型t検定: 主要変数のみの単純比較（参考値）\n"
                        f"- OLS係数表: 共変量調整後の効果推定（主解析）\n"
                        f"- 置換検定: 分布仮定なしの頑健な検定（補強）"
                    )

                    # session_stateに保存
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
                    st.error(f"❌ Welch型t検定は2群比較のみ対応しています。主要変数1 ({main_var}) は {unique_main_levels} 水準です。")
                    st.info("💡 3群以上の場合は、OLSのANOVA表または置換検定をご使用ください。")

            # =============================================================================
            # Permutation Test
            # =============================================================================

            if run_permutation:
                st.markdown("---")
                section_number = "8️⃣" if run_welch else "7️⃣"
                st.markdown(f"## {section_number} 置換検定 (Permutation Test) 結果")
                st.caption("パラメトリックな仮定に依存しない、頑健な検定方法")
                st.warning(f"🎯 **検定対象変数:** 主要変数1 = **{main_var}**")

                with st.spinner(f"{n_permutations:,} 回の置換を実行中..."):
                    try:
                        # デザイン行列の準備（design_info取得）
                        y = analysis_df[pc_col].values
                        X_full = dmatrix(formula.split('~')[1], data=analysis_df, return_type='dataframe')
                        design_info = X_full.design_info  # patsy design_infoを取得

                        # 層化変数の準備
                        stratify_var = None
                        if blocking_var is not None:
                            stratify_var = analysis_df[blocking_var].values
                            st.warning(f"🔹 **層化置換を使用**: {blocking_var} 内でのみ残差を置換します\n\n"
                                      f"これにより、{blocking_var}構造を保持したまま{main_var}効果の帰無分布を生成します")
                        else:
                            st.info("ℹ️ 主要変数2が指定されていないため、非層化置換を使用します")

                        # 片側検定の設定
                        one_sided_param = None
                        if perm_sided == "片側検定: 主要変数1 > 基準 (greater)":
                            one_sided_param = 'greater'
                            st.info(f"📊 **片側検定 (greater)**: {main_var}の第1水準 > 基準（第0水準）を検定")
                        elif perm_sided == "片側検定: 主要変数1 < 基準 (less)":
                            one_sided_param = 'less'
                            st.info(f"📊 **片側検定 (less)**: {main_var}の第1水準 < 基準（第0水準）を検定")

                        # 置換検定の実行
                        if perm_method == "Freedman-Lane法 (推奨)":
                            perm_result = freedman_lane_permutation(
                                y, X_full, main_var,
                                n_perm=n_permutations,
                                stratify_by=stratify_var,
                                one_sided=one_sided_param,
                                design_info=design_info  # 堅牢な列特定のため追加
                            )

                            info_msg = f"ℹ️ **Freedman-Lane法:** 主要変数1（**{main_var}**）の効果を検定\n\n"
                            info_msg += f"他の変数（主要変数2、追加共変量など）は制御変数として扱われます\n\n"
                            if perm_result.get('stratified'):
                                info_msg += f"✅ 層化置換: {blocking_var}構造を保持\n\n"
                            if perm_result.get('one_sided'):
                                info_msg += f"📊 検定方向: {perm_result['one_sided']}\n\n"
                            st.success(info_msg)
                        else:
                            # Simple label permutation（y値置換版：高速）
                            np.random.seed(42)  # Freedman-Laneと同じシードを使用

                            # 観測統計量
                            obs_model = sm.OLS(y, X_full).fit()
                            interest_cols = [col for col in X_full.columns
                                           if f'C({main_var})' in col and col != 'Intercept']
                            main_effect_col = [col for col in interest_cols if ':' not in col][0]
                            obs_stat = float(obs_model.tvalues[main_effect_col])

                            null_dist = []
                            for _ in range(n_permutations):
                                # y値を置換（層化対応）
                                if stratify_var is not None:
                                    # 層化置換：各層内でのみy値をシャッフル
                                    y_perm = np.zeros_like(y)
                                    for stratum in np.unique(stratify_var):
                                        stratum_mask = (stratify_var == stratum)
                                        stratum_indices = np.where(stratum_mask)[0]
                                        perm_indices = np.random.permutation(stratum_indices)
                                        y_perm[stratum_indices] = y[perm_indices]
                                else:
                                    # 非層化置換：全体をシャッフル
                                    y_perm = np.random.permutation(y)

                                # フィット（設計行列は固定）
                                perm_model = sm.OLS(y_perm, X_full).fit()
                                null_dist.append(float(perm_model.tvalues[main_effect_col]))

                            null_dist = np.array(null_dist)

                            # p値計算（連続性補正）
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

                            info_msg = f"ℹ️ **Simple permutation:** y値を置換（高速版）\n\n"
                            info_msg += "📊 **帰無仮説の違い:**\n"
                            info_msg += "- **Freedman-Lane**: 共変量で条件づけた帰無仮説（残差を置換）\n"
                            info_msg += "- **Simple**: 共変量構造を無視した帰無仮説（y値を置換）\n\n"
                            info_msg += "⚠️ **p値の違い**: 帰無仮説が異なるため、p値も異なるのが正常です。\n"
                            if stratify_var is not None:
                                info_msg += f"✅ 層化置換: {blocking_var}内でのみy値を入れ替えます\n\n"
                            else:
                                info_msg += "全サンプルでy値を入れ替えます\n\n"
                            info_msg += "💡 **推奨**: 共変量を制御したい場合は**Freedman-Lane法**を使用してください。"
                            st.warning(info_msg)

                        # 結果の表示
                        col1, col2, col3, col4 = st.columns(4)
                        with col1:
                            st.metric("観測統計量 (t値)", f"{perm_result['observed']:.4f}")
                        with col2:
                            test_type = "片側" if perm_result.get('one_sided') else "両側"
                            st.metric("検定の種類", test_type)
                        with col3:
                            st.metric("Permutation p値", f"{perm_result['p_value']:.4g}")
                        with col4:
                            sig_level = "***" if perm_result['p_value'] < 0.001 else \
                                        "**" if perm_result['p_value'] < 0.01 else \
                                        "*" if perm_result['p_value'] < 0.05 else "n.s."
                            st.metric("有意性", sig_level)

                        # 層化情報の表示
                        if perm_result.get('stratified'):
                            st.success(f"✅ **層化置換**: {blocking_var}の各水準内でのみ残差を置換しました\n\n"
                                      f"各{blocking_var}内で{main_var}効果を無効化した帰無分布を生成")
                        else:
                            st.info("ℹ️ **非層化置換**: 全サンプルにわたって残差を置換しました")

                        # 再現性情報の表示（method_info）
                        if 'method_info' in perm_result:
                            with st.expander("🔬 解析の詳細（再現性のための情報）", expanded=False):
                                info = perm_result['method_info']
                                st.markdown(f"""
**モデル式:**
- **帰無モデル**: `{info['reduced_model']}`
- **対立モデル**: `{info['full_model']}`

**検定設定:**
- **置換回数**: {info['n_permutations']:,} 回
- **乱数シード**: {info['random_seed']}
- **統計量**: {info['test_statistic']}
- **連続性補正**: {'あり (+1補正)' if info['continuity_correction'] else 'なし'}
- **層化**: {info['stratification_var']}

**検定対象**: `{main_var}` の効果（他の変数は制御）
""")

                        # ヒストグラム
                        st.markdown("### 📊 帰無分布 (Null Distribution)")
                        st.caption("灰色: 帰無仮説の下での統計量の分布。赤線: 観測された統計量")

                        fig, ax = plt.subplots(figsize=(10, 6))
                        ax.hist(perm_result['null_distribution'], bins=50,
                               alpha=0.7, edgecolor='black', density=True, color='lightgray',
                               label='Null distribution')
                        ax.axvline(perm_result['observed'], color='red',
                                  linestyle='--', linewidth=2.5, label=f'Observed (t={perm_result["observed"]:.3f})')
                        ax.axvline(-perm_result['observed'], color='red',
                                  linestyle='--', linewidth=2.5)

                        # p値の領域を塗りつぶし
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
                        st.session_state.fig_null_dist = fig  # zipダウンロード用に保存

                        # 結果を保存
                        st.session_state.perm_result = perm_result

                    except Exception as e:
                        st.error(f"❌ 置換検定に失敗しました")
                        st.exception(e)

            # =============================================================================
            # Download Results
            # =============================================================================

            st.markdown("---")
            section_number_download = "9️⃣" if run_welch else "8️⃣"
            st.markdown(f"## {section_number_download} 結果のダウンロード")
            st.caption("解析結果とグラフをまとめてZIPファイルでダウンロードできます")

            # ZIPダウンロードボタンを中央に配置
            col1, col2, col3 = st.columns([1, 2, 1])
            with col2:
                # 結果が1つ以上ある場合のみダウンロードボタンを表示
                has_results = any([
                    'ols_coef' in st.session_state,
                    'lmm_coef' in st.session_state,
                    'perm_result' in st.session_state,
                    'welch_result' in st.session_state
                ])

                if has_results:
                    # ファイル名を生成
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    zip_filename = f"pca_stats_{pc_col}_{timestamp}.zip"

                    # ZIP作成
                    zip_data = create_results_zip(pc_col)

                    st.download_button(
                        label="📦 すべての結果をダウンロード (ZIP)",
                        data=zip_data,
                        file_name=zip_filename,
                        mime="application/zip",
                        help="統計テーブル（CSV）とグラフ（PNG 300dpi + PDF）をZIPファイルにまとめてダウンロード",
                        type="primary",
                        use_container_width=True
                    )

                    # 含まれる内容を表示
                    st.markdown("---")
                    st.markdown("#### 📋 ZIPファイルの内容:")

                    contents = []
                    if 'welch_result' in st.session_state:
                        contents.append("✅ Welch型t検定結果 (CSV)")
                    if 'ols_coef' in st.session_state:
                        contents.append("✅ OLS係数表 (CSV)")
                    if 'ols_anova' in st.session_state:
                        contents.append("✅ ANOVA表 (CSV)")
                    if 'ols_emm' in st.session_state:
                        contents.append("✅ 推定周辺平均 (CSV)")
                    if 'lmm_coef' in st.session_state:
                        contents.append("✅ LMM係数表 (CSV)")
                    if 'perm_result' in st.session_state:
                        contents.append("✅ Permutation test結果 (CSV)")
                    if 'fig_forest_ols' in st.session_state:
                        contents.append("✅ OLS Forest plot (PNG + PDF)")
                    if 'fig_diagnostic' in st.session_state:
                        contents.append("✅ 診断プロット (PNG + PDF)")
                    if 'fig_emm' in st.session_state:
                        contents.append("✅ EMM plot (PNG + PDF)")
                    if 'fig_forest_lmm' in st.session_state:
                        contents.append("✅ LMM Forest plot (PNG + PDF)")
                    if 'fig_null_dist' in st.session_state:
                        contents.append("✅ Null distribution plot (PNG + PDF)")

                    # 2列で表示
                    col_a, col_b = st.columns(2)
                    half = len(contents) // 2 + len(contents) % 2
                    with col_a:
                        for item in contents[:half]:
                            st.markdown(f"- {item}")
                    with col_b:
                        for item in contents[half:]:
                            st.markdown(f"- {item}")

                else:
                    st.info("💡 解析を実行すると、結果をダウンロードできるようになります")

else:
    st.info("👆 データファイルをアップロードして解析を開始してください")

    st.markdown("---")
    st.markdown("### 📋 データ形式の例")
    st.markdown("""
    データは **TSV** または **CSV** 形式で、以下のような構造を想定しています:

    - **1行 = 1サンプル** (pseudobulk レベル)
    - **PC score 列**: PC1, PC2, PC3, PC4 など
    - **カテゴリ変数**: sex, subtype, condition など
    - **任意: Donor/Subject ID** (反復測定がある場合)
    - **任意: 数値共変量**: age, batch, depth など

    **サンプルデータ構造:**
    ```
    sample_id    sex       subtype    donor_id    age    batch    PC1      PC2      PC3      PC4
    sample_001   Female    TypeA      donor_01    45     1        -2.3     1.2      0.5      -0.8
    sample_002   Male      TypeA      donor_02    38     1         1.5    -0.8     -1.2      0.3
    sample_003   Female    TypeB      donor_03    52     1         0.2     2.1      1.8     -1.5
    sample_004   Male      TypeB      donor_04    41     2        -1.1    -1.5      0.9      1.2
    ...
    ```

    テストデータは `/home/ichiro/streamlit/temp/pca_test_data.tsv` にあります。
    """)

st.markdown("---")

with st.expander("📚 統計手法の詳細", expanded=False):
    st.markdown("""
#### **OLS (Ordinary Least Squares - 通常最小二乗法)**
固定効果モデル。全てのサンプルが独立であると仮定します。

- **Classical SE**: 標準的な標準誤差（等分散性を仮定）
- **HC3 Robust SE**: 不均一分散に頑健（サンプルサイズが小さい場合に推奨）
- **Cluster-robust SE**: クラスター内相関を考慮（例: 同一ドナーのサンプル間の相関）
- **Welch型**: 等分散を仮定しない（2群比較のみ）
- **WLS (加重最小二乗法)**: 不均一分散に対応

**適用例**: 各ドナーから1サンプルのみ、またはドナー間の変動が小さい場合

---

#### **LMM (Linear Mixed Model - 線形混合モデル)**
階層構造を持つデータに対応。ランダム効果で個体間変動を考慮します。

- **ランダム切片 `(1|donor)`**: 各ドナーで切片が異なることを許容
- **REML**: 不偏な分散推定（推奨）
- **ML**: モデル比較用（AIC/BICの比較）
- **ICC (級内相関係数)**: ドナー間変動の割合を示す

**適用例**: 各ドナーから複数サンプル、または技術的繰り返しがある場合

**注意**: ドナー数が5未満の場合、分散推定が不安定になります。

---

#### **Permutation Test (置換検定)**
分布仮定に依存しない、ノンパラメトリックな検定方法。

**🎯 検定対象:** 主要変数1のみが検定されます。主要変数2やその他の共変量は制御変数として扱われます。

- **Freedman-Lane法**: 共変量を制御した厳密な置換検定（推奨）
  - 主要変数1の効果を検定
  - 帰無モデルと対立モデルをフィット
  - 帰無モデルの残差を置換して擬似データを生成
  - 対立モデルを再フィットして統計量を計算

- **Simple permutation**: ラベルを単純に入れ替える方法
  - 主要変数1のラベルを入れ替え
  - 共変量がない、または少ない場合に使用

**適用例**:
- サンプルサイズが小さい（<30）
- 残差の正規性が疑わしい
- 外れ値が存在する

---

#### **ANOVA Type II vs Type III**
- **Type II** (推奨): 各変数の主効果を他の変数で調整して検定（交互作用なしの場合）
- **Type III**: 全ての交互作用を含めて検定（交互作用ありの場合）

---

#### **EMM (Estimated Marginal Means - 推定周辺平均)**
共変量を調整した後の、各群の予測平均値。

- 数値共変量を平均値に固定
- カテゴリ変数の各水準での予測値を計算
- 信頼区間を付与して比較

**適用例**:
- サブタイプ別の性差を可視化
- 年齢やバッチを調整した群間比較

---

### ⚠️ 重要な注意事項

#### **1. 完全共線性 (Perfect Collinearity)**
2つの変数が完全に一致している場合、効果を分離できません。

**例**: ある subtype が全て Female の場合、sex 効果と subtype 効果は識別不能

**対処法**:
- いずれかの変数を削除
- より均衡の取れたデータを収集
- 記述的な比較に留める

#### **2. サンプルサイズ**
- **LMMのドナー数**: 5以上を推奨（少ないと分散推定が不安定）
- **各セルの最小サンプル数**: 3以上を推奨
- **置換検定**: サンプルサイズが小さい場合（<30）に特に有用

#### **3. 交互作用項**
データが十分にある場合のみ使用を推奨。各セルに十分なサンプル（≥5）が必要。

#### **4. 多重検定**
複数のPCを解析する場合は、**Benjamini-Hochberg FDR補正**を適用してください。

#### **5. 数値共変量の標準化**
年齢など、スケールが大きく異なる数値共変量は、事前に標準化することを推奨します。

---

### 📖 参考文献

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

**開発者**: Claude Code
**バージョン**: 1.0
**最終更新**: 2025-01
""")
