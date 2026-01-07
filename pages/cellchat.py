import scanpy as sc
import numpy as np
import streamlit as st
import pandas as pd
#import fireducks.pandas as pd
import pickle
import re
import os
import time
import matplotlib
matplotlib.use("cairo") # PDFoftransparentclearizepairplan interactivewithisumakurowkanot?
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from matplotlib.colors import LinearSegmentedColormap
#from scipy import sparse
import scipy.sparse
#from scipy.stats import zscore
from statsmodels.stats.multitest import multipletests
import traceback
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter
import rpy2.robjects.packages as rpackages
from rpy2.robjects.packages import importr
from joblib import Parallel, delayed
from tqdm import tqdm
#print("rpy2 version:", rpy2.__version__)
import scipy
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.offline import plot
from numba import njit, prange
# import scipy.stats as stats
from typing import Union, Optional, List
from streamlit_sortables import sort_items
from matplotlib.patches import FancyArrowPatch
from matplotlib.path import Path
import math
import logging
import matplotlib.cm as cm
import shutil
import io
import zipfile
from datetime import datetime, timedelta
logger = logging.getLogger("streamlit.runtime.scriptrunner_utils.script_run_context")
logger.disabled = True

from pages.cellchat_vis import (
    netVisual_circle, netVisual_circle_individual,
    netVisual_aggregate, netVisual_heatmap,
    netAnalysis_signalingRole_network, plotGeneExpression)

def get_r_permutation(nC, nboot, seed=1):
    """
    Rsaywordofsample.intrelatedcounttheuseandplacechangeDatathegenerateshi、
    Errorprocesstheincludemetafirmrobustnaimplement
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr
        
        # Rofbasemainpake-jitheinpo-to
        base = importr('base')
        
        # Rwithrandomcountshi-dotheSettings
        ro.r('set.seed({})'.format(seed))
        
        # RsidewithpermutationDatathegenerate（toraburushiyu-teinguforofkomentoattachki）
        r_code = """
        nC <- {}
        nboot <- {}
        cat("Generating permutation in R with nC =", nC, "and nboot =", nboot, "\\n")
        permutation <- replicate(nboot, sample.int(nC, size = nC))
        cat("R permutation dimensions:", dim(permutation), "\\n")
        permutation
        """.format(nC, nboot)
        
        # Rko-dotheRunandResulttheget
        r_permutation = ro.r(r_code)
        
        # RmatrixtheNumPydistributecolumntoconvert
        permutation = np.array(r_permutation).astype(int)
        
        # Ris1be-sunaofwith、Pythonof0be-sutoconvert
        permutation = permutation - 1
        
        # distributecolumnofshapestatetheConfirmshi、requiredtorespondjiterotateplace
        if permutation.shape[0] == nboot and permutation.shape[1] == nC:
            permutation = permutation.T  # CellChatofperiodwaitdoshapeformattomatchwaseru
        
        print(f"Successfully generated R permutation with shape: {permutation.shape}")
        print(f"Permutation range: min={np.min(permutation)}, max={np.max(permutation)}")
        print(f"Expected range: 0 to {nC-1}")
        return permutation
        
    except Exception as e:
        print(f"Error generating R permutation: {str(e)}")
        print("Falling back to numpy permutation")
        
        # fo-rubaku: NumPyofrandomcountgeneratetheuse
        np.random.seed(seed)
        permutation = np.zeros((nC, nboot), dtype=int)
        for i in range(nboot):
            permutation[:, i] = np.random.permutation(nC)
        return permutation

@st.cache_data
def create_cell_color_mapping(cell_list, palette_name):
    """
    Cellnameandcolorofonepierceshitamapinguthecreatedorelatedcount

    Parameters
    ----------
    cell_list : list
        Cellnameoflist
    palette_name : str
        usedoseparatescattercolor-paretoname

    Returns
    -------
    dict
        Cellnamethekey、colorthebariyu-anddowordwrite

    Note
    ----
    pointsetsaretaseparatescatterparetoofcase、mazudeforutoofcolorcount（base_palette oflength）thegetshi、
    Cellcountissoofcountorbelownarasoofmamausefor、supererucaseis sns.color_palette() of n_colors to
    Cellcountthepointsetandsupplementbetweencolorthegenerateshimasu。
    """
    n_cells = len(cell_list)
    # deforutoofseparatescatterparetotheget（Example: "Set1"nara9color）
    base_palette = sns.color_palette(palette_name)
    base_n = len(base_palette)
    
    if n_cells <= base_n:
        # basemainparetoofcolorcountorinnara、aheadheadfromrequirednacounttheusefor
        colors = base_palette[:n_cells]
    else:
        # Cellcountisbasemainparetoofcolorcountthesupererucase、requirednacountofcolorthegenerate（lineshapesupplementbetween）
        colors = sns.color_palette(palette_name, n_colors=n_cells)
        
    # Cellnameandcolorthepairrespondattachketawordwritethereturnsu
    return {cell: color for cell, color in zip(cell_list, colors)}

def should_clear_session_state(uploaded_file):
    """trueofnewFilejudgeset（sizedependexisttheexcluderemove）"""
    if not uploaded_file:
        return False
    
    # firsttimeUpload
    if 'last_file_name' not in st.session_state:
        return True
    
    # Filenameischangefurthersareta
    if uploaded_file.name != st.session_state.get('last_file_name', ''):
        return True
    
    # AnalysisResultisexistatshinot（Errorafteretc）
    if st.session_state.get('cellchat_res') is None:
        return True
    
    return False

def clear_cellchat_session_state():
    """newFiletimeoffitcutnasetionstatestatekuria"""
    # newFiletimeisSettingsmokuria（FileischangewaruandCell typemochangewarufor）
    keys_to_clear = ['cellchat_res', 'cellchat_temp_dir', 'sorted_order', 
                     'cell_color_map', 'current_cmap']
    for key in keys_to_clear:
        if key in st.session_state:
            del st.session_state[key]
    st.cache_data.clear()

@st.cache_data
def sanitize_filename(filename, max_length=20):
    """Filenamethesafealltoandmaximumlengththecontrollimitdo"""
    # specialspecialtextcharthedeleteremovealsoisplacechange
    filename = re.sub(r'[\\/*?:"<>|]', "_", filename)
    # lengththecontrollimit
    if len(filename) > max_length:
        base, ext = os.path.splitext(filename)
        filename = base[:max_length] + ext
    return filename

def save_cellchat_result(result, uploaded_filename, selected_types, output_dir="temp"):
    """CellChatofResultthepickleFileandandSave"""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Filenamethesafealltoandlengththecontrollimit
    safe_filename = sanitize_filename(uploaded_filename, 20)
    
    # signaltaiputheshortishapeformatto
    signal_types = "_".join([sig_type[:3] for sig_type in selected_types])
    
    # SaveFilenamethecreate
    save_filename = f"cellchat_{safe_filename}_{signal_types}.pkl"
    save_path = os.path.join(output_dir, save_filename)
    
    # ResulttheSave
    with open(save_path, 'wb') as f:
        pickle.dump(result, f)
    
    return save_path


def identify_overexpressed_genes(
    adata,
    group_by: str = None,
    idents_use: list = None,
    invert: bool = False,
    features_name: str = "features",
    only_pos: bool = True,
    features: list = None,
    return_object: bool = True,
    thresh_pct: float = 0,
    thresh_fc: float = 0,
    thresh_p: float = 0.05,
    do_de: bool = True,
    do_fast: bool = True,
    min_cells: int = 10,
    min_cells_expr: int = 10
):
    if not hasattr(pd.DataFrame, "iteritems"): #iteritemsisnotcaseofpairrespond
        pd.DataFrame.iteritems = pd.DataFrame.items

    # selfmoveconvertthevalidize
    numpy2ri.activate()
    pandas2ri.activate()
    presto = rpackages.importr("presto")
    # Datamatrixoflevelprepare（herewithis genes x cells thebeforeprovideanddo）
    if scipy.sparse.issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = np.asarray(adata.X)
    # ※ AnnData withis X ispassnormal cells×genes naofwith、requiredtorespondjiterotateplaceandplease
    # Example: X = X.T  # korewith X ofrowisGene、columnisCellandbecome

    # serulabelofget
    labels = adata.obs[group_by] if group_by else adata.obs.iloc[:, 0]
    labels = pd.Categorical(labels)
    
    # Genelistofget
    features_use = list(adata.var_names) if features is None else list(set(features) & set(adata.var_names))
    
    # oversurplusExpressionAnalysisforofmatrixthecreate（herewithis X ofrowisGene、columnisCellandtempset）
   # data_use = X[np.array([adata.var_names.get_loc(gene) for gene in features_use]), :]
    X_gene = X.T 
    data_use = X_gene[np.array([adata.var_names.get_loc(gene) for gene in features_use]), :]
    
    # eachkurasutaoflabel
    level_use = list(labels.categories)
    
    if do_de:
        if do_fast:
            # R of presto::wilcoxauc thecallbioutsu
            # data_use: rowisGene、columnisCellofDatamatrix
            # labels: CellgoandofGroup（lengthisCellcount）
            with localconverter(ro.default_converter + pandas2ri.converter):
                # data_use_df andand pandas.DataFrame toconvert
                data_use_df = pd.DataFrame(data_use, index=features_use, columns=adata.obs_names)
                # R withisrownameis features、columnnameisSamplename
            
            # Dataofnextsourcechieku
            print(f"data_use_df shape: {data_use_df.shape}, labelcount: {len(labels)}")
            
            # nextsourceisonecauseshinotcase、Filepassreasonwithfo-rubaku
            if data_use_df.shape[1] != len(labels):
                print("Warning: Dataofnextsourcenotonecause。FilepassreasonwithRtoprocessthedependrelyshimasu")
                genes_de = run_presto_through_files(data_use_df, labels, level_use)
            else:
                try:
                    # passnormalofrpy2passreasonofprocessthetrymiru
                    # labels the R of factor toconvert
                    r_labels = ro.FactorVector(labels.astype(str).tolist())
                    # groups_use Parametertois、pairobjectanddoGroupofbekutoruthetransfersu
                    r_groups = ro.StrVector(level_use)
                    
                    # presto of wilcoxauc callbioutshi
                    res = presto.wilcoxauc(data_use_df, r_labels, groups_use=r_groups)
                    
                    # res is R of DataFrame tobecomeofwith、pandas toconvert
                    with localconverter(ro.default_converter + pandas2ri.converter):
                        genes_de = pd.DataFrame(ro.conversion.rpy2py(res))
                except Exception as e:
                    print(f"rpy2passreasonofRunError: {str(e)}")
                    print("FilepassreasonwithRtoprocessthedependrelyshimasu")
                    genes_de = run_presto_through_files(data_use_df, labels, level_use)
            
            # Rversionandsamemannerto、columnnameofchangefurther
            genes_de.rename(columns={
                "group": "clusters",
                "feature": "features",
                "pval": "pvalues",
                "logFC": "logFC",
                "pct_in": "pct.1",
                "pct_out": "pct.2",
                "padj": "pvalues.adj"
            }, inplace=True)
            genes_de.loc[:, "logFC_abs"]  = genes_de["logFC"].abs().copy()
            # pct.max is pct.1 and pct.2 ofmaximumvalue
            genes_de.loc[:, "pct.max"] = genes_de[["pct.1", "pct.2"]].max(axis=1)
            
            # Filtering
            markers_all = genes_de[(genes_de["pvalues"] < thresh_p) &
                                    (genes_de["logFC_abs"] >= thresh_fc) &
                                    (genes_de["pct.max"] > thresh_pct * 100)].copy()
            markers_all.sort_values("pvalues", inplace=True)
            
        else:
            # followcomeof Python implement（scipy of Mann–Whitney U checksettheuse）
            markers = []
            for cluster in level_use:
                cluster_mask = (labels == cluster)
                cluster_data = data_use[:, cluster_mask]
                other_data = data_use[:, ~cluster_mask]
                for i, gene in enumerate(features_use):
                    pct1 = np.mean(cluster_data[i, :] > 0)
                    pct2 = np.mean(other_data[i, :] > 0)
                    stat, p_val = scipy.stats.mannwhitneyu(cluster_data[i, :], other_data[i, :], alternative='two-sided')
                    log_fc = np.log(np.mean(cluster_data[i, :] + 1e-10) / np.mean(other_data[i, :] + 1e-10))
                    markers.append({
                        "clusters": cluster,
                        "features": gene,
                        "pvalues": p_val,
                        "logFC": log_fc,
                        "pct.1": pct1,
                        "pct.2": pct2
                    })
            markers_all = pd.DataFrame(markers)
            markers_all["logFC_abs"] = markers_all["logFC"].abs()
            markers_all["pct.max"] = markers_all[["pct.1", "pct.2"]].max(axis=1)
            markers_all = markers_all[(markers_all["pvalues"] < thresh_p) &
                                      (markers_all["logFC_abs"] >= thresh_fc) &
                                      (markers_all["pct.max"] > thresh_pct)].copy()
            markers_all.sort_values("pvalues", inplace=True)
        
        features_sig = markers_all["features"].unique()
    else:
        st.warning("DEG filter is off.")
        # do.DE=False ofcase：eachGeneis min_cells_expr orupExpressionandexistkawithselectsep
        nCells = (data_use > 0).sum(axis=1)
        markers_all = pd.DataFrame({
            "features": features_use,
            "nCells": nCells
        })
        markers_all = markers_all[markers_all["nCells"] >= min_cells_expr]
        features_sig = markers_all["features"].unique()
    
    # Resultthewordwriteandandreturnsu
    print(f"{len(features_sig)} genes passed filtering")
    return {"features": features_sig, "markers_all": markers_all}


def run_presto_through_files(data_df, labels, groups_use):
    """
    FilepassreasonwithR prestotheRundo
    
    Parameters:
    -----------
    data_df : pandas.DataFrame
        Gene×Cellofmatrix
    labels : pandas.Series
        CellofGrouplabel
    groups_use : list
        usedoGrouplist
        
    Returns:
    --------
    pandas.DataFrame
        wilcoxaucResult
    """
    import os
    import tempfile
    import subprocess
    import uuid
    
    # onetimedeirekutorithecreate
    temp_dir = tempfile.mkdtemp(prefix="presto_")
    
    try:
        # uni-kunasetionID
        session_id = str(uuid.uuid4())
        
        # Filepasu
        data_path = os.path.join(temp_dir, f"data_{session_id}.csv")
        labels_path = os.path.join(temp_dir, f"labels_{session_id}.csv")
        groups_path = os.path.join(temp_dir, f"groups_{session_id}.csv")
        output_path = os.path.join(temp_dir, f"result_{session_id}.csv")
        
        # DatatheCSVandandSave
        data_df.to_csv(data_path)
        pd.DataFrame({"label": labels.astype(str)}).to_csv(labels_path, index=False)
        pd.DataFrame({"group": groups_use}).to_csv(groups_path, index=False)
        
        # Rsukuriputothecreate
        r_script_path = os.path.join(temp_dir, f"script_{session_id}.R")
        r_code = f"""
        # requirednaraiburarithero-do
        library(presto)
        
        # DataofLoading
        data <- read.csv("{data_path}", row.names=1)
        labels <- read.csv("{labels_path}")$label
        groups <- read.csv("{groups_path}")$group
        
        # DatatheConfirm
        print(dim(data))
        print(length(labels))
        print(groups)
        
        # DataofPreprocessing
        data_matrix <- as.matrix(data)
        
        # wilcoxauc ofRun
        result <- presto::wilcoxauc(data_matrix, labels, groups_use=groups)
        
        # ResulttheSave
        write.csv(result, "{output_path}", row.names=FALSE)
        
        print("processComplete")
        """
        
        with open(r_script_path, "w") as f:
            f.write(r_code)
        
        # RsukuriputotheRun
        try:
            subprocess.run(["Rscript", r_script_path], check=True)
        except subprocess.CalledProcessError as e:
            print(f"RsukuriputoRunError: {e}")
            raise
        
        # Resultthereadmiintomu
        if os.path.exists(output_path):
            result_df = pd.read_csv(output_path)
            return result_df
        else:
            raise FileNotFoundError(f"ResultFileisviewtsukarimasen: {output_path}")
            
    finally:
        # onetimeFileofkuri-napu（Option）
        # shutil.rmtree(temp_dir)
        print(f"onetimeFileis {temp_dir} toSavesareteimasu")
    
    return None


@njit
def compute_hill_outer(dataL, dataR, k, n):
    numCluster = dataL.shape[0]
    result = np.empty((numCluster, numCluster), dtype=np.float64)
    kn = k ** n
    for i in prange(numCluster):
        for j in range(numCluster):
            x = dataL[i] * dataR[j]
            xn = x ** n
            result[i, j] = xn / (kn + xn)
    return result

# PythonimplementofCellChatDBprocessthedebagudoko-do
# get_cellchatdb_from_rrelatedcountofaftertoorbelowofrelatedcounttheadd

def debug_cellchatdb(db):
    """CellChatDBofincontentthedetailstoOutputandverify"""
    print("==== CellChatDBverify ====")
    
    if 'interaction' in db:
        interaction = db['interaction']
        print(f"Ligand-Receptorpaircount: {len(interaction)}")
        if len(interaction) > 0:
            print("firstof5pair:")
            cols_to_show = ['ligand', 'receptor']
            if 'interaction_name' in interaction.columns:
                cols_to_show.insert(0, 'interaction_name')
            if 'pathway_name' in interaction.columns:
                cols_to_show.append('pathway_name')
            print(interaction[cols_to_show].head(5))
            
            # Ligand,ReceptorofkindtypetheConfirm
            ligands = interaction['ligand'].unique()
            receptors = interaction['receptor'].unique()
            
            print(f"\nuni-kunaLigandcount: {len(ligands)}")
            print(f"LigandofExample (firstof10): {', '.join(map(str, ligands[:10]))}")
            
            print(f"\nuni-kunaReceptorcount: {len(receptors)}")
            print(f"ReceptorofExample (firstof10): {', '.join(map(str, receptors[:10]))}")
            
            # Ligand,ReceptoroftypetheConfirm
            print(f"\nLigandofDatatype: {type(interaction['ligand'].iloc[0])}")
            print(f"ReceptorofDatatype: {type(interaction['receptor'].iloc[0])}")
            
            # complexofLigand,ReceptorofConfirm
            complex_ligands = [l for l in ligands if isinstance(l, str) and '_' in l]
            complex_receptors = [r for r in receptors if isinstance(r, str) and '_' in r]
            
            print(f"\ncomplexLigandcount: {len(complex_ligands)}")
            if complex_ligands:
                print(f"complexLigandofExample: {', '.join(complex_ligands[:5])}")
                
            print(f"complexReceptorcount: {len(complex_receptors)}")
            if complex_receptors:
                print(f"complexReceptorofExample: {', '.join(complex_receptors[:5])}")

            # Databe-sugetafter、structmakethedetailstoConfirm
            print("DBstructmake:")
            for key in db:
                print(f"{key}: {type(db[key])}, rowcount: {len(db[key])}")

            print("Ligandoffirstof5tsuofExample:")
            print(interaction['ligand'].head(5).tolist())
            print("LigandofDatatype:", interaction['ligand'].dtype)
    else:
        print("CellChatDBtointeractionisincludedmasen")
    
    # complex, cofactorofConfirm
    if 'complex' in db and not db['complex'].empty:
        print(f"\ncomplexinfocount: {len(db['complex'])}")
        print("complexofExample (firstof3):")
        print(db['complex'].head(3))
    else:
        print("\ncomplexinfoisincludednotkaemptywithsu")
        
    if 'cofactor' in db and not db['cofactor'].empty:
        print(f"\ncofactorinfocount: {len(db['cofactor'])}")
        print("cofactorofExample (firstof3):")
        print(db['cofactor'].head(3))
    else:
        print("\ncofactorinfoisincludednotkaemptywithsu")
        
    print("==== CellChatDBverifyComplete ====\n")
    
    return db  # sourceofDBthesoofmamareturnandsourceofprocesstheinheritcontinuepossibleto

    
def hill_function(x, k, n):
    """
    bekutoruizesaretahirurelatedcountofCalculation
    countformat: y = x^n / (k^n + x^n)
    x: Inputdistributecolumn
    k: hirusetcount
    n: hirurelatecount
    """
    x_n = np.power(x, n)
    return x_n / (np.power(k, n) + x_n)

def compute_hill_outer_vectorized(dataL, dataR, k, n):
    """
    dataLanddataRofoutaccumulatetopairandhirurelatedcountthefitfordo（bekutoruizeversion）
    dataL: 1nextsourcedistributecolumn（LigandExpressionvalue）
    dataR: 1nextsourcedistributecolumn（ReceptorExpressionvalue）
    k, n: hirurelatedcountofParameter
    returnrivalue: Calculationsareta2nextsourcedistributecolumn
    """
    outer = np.outer(dataL, dataR)
    return hill_function(outer, k, n)

@st.cache_data
def get_cellchatdb_from_r(species="human"):
    """
    RofCellChatpake-jifromCellChatDBthegetdo（rpy2 3.5.1topairrespond）
    
    Parameters
    ----------
    species : str
        'human' alsois 'mouse' thepointsetandDatabe-sutheSelect
        
    Returns
    -------
    dict
        'interaction', 'complex', 'cofactor', 'geneInfo' theincludemuCellChatDBofwordwrite（eachvalueis Pandas DataFrame）
    """
    try:
        # selfmoveconvertofvalidize
        pandas2ri.activate()
        # basepake-jiis importr theusetereadmiintomu
        base = importr("base")
        # CellChatpake-jiofro-do
        ro.r("library(CellChat)")
        
        # pairobjectofDatasetothero-do
        if species.lower() == "human":
            r_command = """
            data("CellChatDB.human")
            CellChatDB.human <- CellChatDB.human
            """
            db_name = "CellChatDB.human"
        elif species.lower() == "mouse":
            r_command = """
            data("CellChatDB.mouse")
            CellChatDB.mouse <- CellChatDB.mouse
            """
            db_name = "CellChatDB.mouse"
        else:
            raise ValueError("species is 'human' alsois 'mouse' thepointsetandplease")
        
        # RkomandoofRun
        ro.r(r_command)
        
        db = {}
        components = ['interaction', 'complex', 'cofactor', 'geneInfo']
        for comp in components:
            # eachkonpo-nentotheget（selfmoveconverttothan Pandas DataFrame tonateexistiszu）
            r_data = ro.r(f"{db_name}${comp}")
            # moshi r_data is Pandas DataFrame narasoofmamauseu
            if isinstance(r_data, pd.DataFrame):
                py_data = r_data
            else:
                # converttherowu（passnormalisheretoiscomenotiszuwithsu）
                with ro.default_converter + pandas2ri.converter as cv:
                    py_data = ro.conversion.rpy2py(r_data)
            db[comp] = py_data


        # get_cellchatdb_from_rrelatedcountinwithofmodifycorrect
        # Datafure-mugetafter、typetheclearshowaltoconvert
        if 'ligand' in db['interaction'].columns:
            db['interaction'].loc[:, 'ligand'] = db['interaction']['ligand'].astype(str)
        if 'receptor' in db['interaction'].columns:
            db['interaction'].loc[:, 'receptor'] = db['interaction']['receptor'].astype(str)



        return db
    except Exception as e:
        import traceback
        print("CellChatDBofgetintoErrorisoccurred: " + str(e))
        traceback.print_exc()
        return None


def debug_cellchatdb(db):
    """CellChatDBofincontentthedetailstoOutputandverify"""
    print("==== CellChatDBverify ====")
    
    if 'interaction' in db:
        interaction = db['interaction']
        print(f"Ligand-Receptorpaircount: {len(interaction)}")
        if len(interaction) > 0:
            print("firstof5pair:")
            cols_to_show = ['ligand', 'receptor']
            if 'interaction_name' in interaction.columns:
                cols_to_show.insert(0, 'interaction_name')
            if 'pathway_name' in interaction.columns:
                cols_to_show.append('pathway_name')
            print(interaction[cols_to_show].head(5))
            
            # Ligand,ReceptorofkindtypetheConfirm
            ligands = interaction['ligand'].unique()
            receptors = interaction['receptor'].unique()
            
            print(f"\nuni-kunaLigandcount: {len(ligands)}")
            print(f"LigandofExample (firstof10): {', '.join(map(str, ligands[:10]))}")
            
            print(f"\nuni-kunaReceptorcount: {len(receptors)}")
            print(f"ReceptorofExample (firstof10): {', '.join(map(str, receptors[:10]))}")
            
            # Ligand,ReceptoroftypetheConfirm
            print(f"\nLigandofDatatype: {type(interaction['ligand'].iloc[0])}")
            print(f"ReceptorofDatatype: {type(interaction['receptor'].iloc[0])}")
            
            # complexofLigand,ReceptorofConfirm
            complex_ligands = [l for l in ligands if isinstance(l, str) and '_' in l]
            complex_receptors = [r for r in receptors if isinstance(r, str) and '_' in r]
            
            print(f"\ncomplexLigandcount: {len(complex_ligands)}")
            if complex_ligands:
                print(f"complexLigandofExample: {', '.join(complex_ligands[:5])}")
                
            print(f"complexReceptorcount: {len(complex_receptors)}")
            if complex_receptors:
                print(f"complexReceptorofExample: {', '.join(complex_receptors[:5])}")
    else:
        print("CellChatDBtointeractionisincludedmasen")
    
    # complex, cofactorofConfirm
    if 'complex' in db and not db['complex'].empty:
        print(f"\ncomplexinfocount: {len(db['complex'])}")
        print("complexofExample (firstof3):")
        print(db['complex'].head(3))
    else:
        print("\ncomplexinfoisincludednotkaemptywithsu")
        
    if 'cofactor' in db and not db['cofactor'].empty:
        print(f"\ncofactorinfocount: {len(db['cofactor'])}")
        print("cofactorofExample (firstof3):")
        print(db['cofactor'].head(3))
    else:
        print("\ncofactorinfoisincludednotkaemptywithsu")
        
    print("==== CellChatDBverifyComplete ====\n")
    
    return db  # sourceofDBthesoofmamareturnandsourceofprocesstheinheritcontinuepossibleto

# deirekutoriprocessforherupa-relatedcount
@st.cache_data
def clear_old_directories(path):
    
    now = datetime.now()
    for dir_name in os.listdir(path):
        dir_path = os.path.join(path, dir_name)
        if os.path.isdir(dir_path):
            try:
                # deirekutorinamefromtaimusutanputheget
                timestamp = float(dir_name)
                dir_time = datetime.fromtimestamp(timestamp)
                # 24timebetweenorupbeforeofdeirekutorithedeleteremove
                if now - dir_time > timedelta(hours=24):
                    shutil.rmtree(dir_path)
            except:
                pass

@st.cache_data
def clear_old_files(path):
    
    now = datetime.now()
    for file_name in os.listdir(path):
        file_path = os.path.join(path, file_name)
        if os.path.isfile(file_path):
            file_time = datetime.fromtimestamp(os.path.getmtime(file_path))
            # 24timebetweenorupbeforeofFilethedeleteremove
            if now - file_time > timedelta(hours=24):
                os.remove(file_path)

@st.cache_data
def check_species_index(gene_list):
    """Geneshinboruofbigtextcharoutpresentfrequentdegreefromhitokamausukatheinfermeasuredo"""
    if not gene_list:
        return 0  # emptyoflistofcaseisdeforutowithhito
    
    # sanpuringu（Genelistisbigkiicase）
    sample_genes = gene_list[:500] if len(gene_list) > 500 else gene_list
    
    # bigtextcharofoutpresentfrequentdegreetheCalculation
    uppercase_ratios = []
    for gene in sample_genes:
        if not gene or not isinstance(gene, str):
            continue
        uppercase_count = sum(1 for char in gene if char.isupper())
        ratio = uppercase_count / len(gene) if len(gene) > 0 else 0
        uppercase_ratios.append(ratio)
    
    # meanvaluetheCalculation
    avg_uppercase_ratio = sum(uppercase_ratios) / len(uppercase_ratios) if uppercase_ratios else 0
    
    # ResulttherogutoOutput
    print(f"Geneshinboruofbigtextcharoutpresentfrequentdegree: {avg_uppercase_ratio:.2f}")
    
    # hito：bigtextcharismanyi（BRCA1）、mausu：smalltextcharismanyi（Brca1）
    return 1 if avg_uppercase_ratio > 0.5 else 0  # 0.5not yetfullnaramausu、soreorupnarahito

def find_first_index_or_zero(lst, elements):
    for element in elements:
        try:
            return lst.index(element)
        except ValueError:
            continue
    return 0

@st.cache_data
def read_h5ad(file):
    adata = sc.read_h5ad(file)
    return adata

# mostfitizesaretamatrixprocessrelatedcount
def optimize_matrix_operations(X, adata):
    """
    matrixoperatemakethemostfitizedorelatedcount
    """
    logger.info(f"Xshapestate: {X.shape}, adata.obs_nameslength: {len(adata.obs_names)}, adata.var_nameslength: {len(adata.var_names)}")
    
    # supa-sumatrixofconvert（onedegreetoconvertandmemorieffectratethedirectionup）
    if scipy.sparse.issparse(X):
        logger.info("supa-sumatrixthedensematrixtoconvertandimasu")
        
        # memorieffectrateoffortopartdivaltoconvert
        chunk_size = 5000  # fitcutnachiyankusizetheSelect
        
        if X.shape[0] > chunk_size:
            logger.info(f"bigkinamatrixthe{chunk_size}rowzutsuconvertshimasu")
            result = []
            for i in range(0, X.shape[0], chunk_size):
                end = min(i + chunk_size, X.shape[0])
                result.append(X[i:end].toarray())
            X = np.vstack(result)
        else:
            X = X.toarray()
    
    return X


def calculate_mean_expression_optimized(data_use, cell_labels, cell_types, min_cells, FunMean):
    """
    Calculate mean expression for ALL cell types (matching R behavior)
    """
    data_use_avg_dict = {}
    cell_counts = {}
    
    # Get cell indices for all cell types
    cell_type_indices = {}
    for cell_type in cell_types:
        indices = np.where(cell_labels == cell_type)[0]
        cell_counts[cell_type] = len(indices)
        cell_type_indices[cell_type] = indices
    
    # Calculate for ALL cell types (don't filter yet)
    for cell_type in cell_types:
        indices = cell_type_indices[cell_type]
        # Perform calculation even for cell types with few cells
        gene_chunk_size = 1000
        n_genes = data_use.shape[1]
        avg_expr = np.zeros(n_genes)
        
        for i in range(0, n_genes, gene_chunk_size):
            end = min(i + gene_chunk_size, n_genes)
            data_subset = data_use[indices, i:end]
            avg_expr[i:end] = np.apply_along_axis(FunMean, 0, data_subset)
        
        data_use_avg_dict[cell_type] = avg_expr
    
    
    return data_use_avg_dict, cell_counts



# logger = logging.getLogger("CellChat")


def apply_mean_function(data_subset, fun_type='triMean', trim=0.1):
    """meanCalculationrelatedcount（Numbanashi）"""
    if fun_type == 'triMean':
        # fourdivrankcountofCalculationandmean（Numbanot yetuse）
        q1 = np.quantile(data_subset, 0.25, axis=0, method='linear')
        q2 = np.quantile(data_subset, 0.5, axis=0, method='linear')  # median
        q3 = np.quantile(data_subset, 0.75, axis=0, method='linear')
        return (q1 + 2*q2 + q3) / 4
    elif fun_type == 'truncatedMean':
      # correctshiitorimumeanofimplement
      from scipy.stats import trim_mean
      if data_subset.ndim == 1:
          return trim_mean(data_subset, proportiontocut=trim)

      else:
          # 2nextsourcedistributecolumnofcase、eachcolumn（Gene）totsuiteCalculation
          return np.array([trim_mean(data_subset[:, i], proportiontocut=trim)
                             for i in range(data_subset.shape[1])])
    else:
        # deforutoissimplepuremean
        return np.mean(data_subset, axis=0)

def process_single_permutation(data_use, cluster_indices, cell_types, fun_type='triMean', trim=0.1):
    """simpleonepermutationofallkurasutaandallGeneofmeanExpressiontheCalculation（Numbanashi）"""
    n_genes = data_use.shape[1]
    numCluster = len(cell_types)
    result = np.zeros((n_genes, numCluster), dtype=np.float32)
    
    # eachClustertotsuiteprocess
    for ct_idx in range(numCluster):
        cells = cluster_indices[ct_idx]
        if len(cells) > 0:
            data_subset = data_use[cells]
            result[:, ct_idx] = apply_mean_function(data_subset, fun_type, trim)
    
    return result

def process_permutation_batch(batch_indices, data_use, cell_labels, permutation, cell_types, 
                             fun_type='triMean', trim=0.1):
    """permutationofbachiprocessrelatedcount"""
    n_genes = data_use.shape[1]
    numCluster = len(cell_types)
    results = np.zeros((n_genes, numCluster, len(batch_indices)), dtype=np.float32)
    
    for idx, j in enumerate(batch_indices):
        # jnumbereyeofPermutationafterofCell type
        group_boot = cell_labels.values[permutation[:, j]]
        
        # eachClustertobelongdoCellofindextheget
        cluster_indices = [np.where(group_boot == ct)[0] for ct in cell_types]
        
        # allofGeneandCell typeofmeanExpressiontheCalculation
        results[:, :, idx] = process_single_permutation(
            data_use, cluster_indices, cell_types, fun_type, trim
        )
    
    return results

def precompute_gene_expressions(data_use, cell_labels, permutation, cell_types, FunMean, nboot, n_jobs=32,  fun_type="triMean", trim=0.1):
    """
    thingbeforetoallGene×allCluster×allpermutationofmeanExpressionthemaandmeteCalculation
    （Numbanashiofba-jiyon）
    
    Parameters
    ----------
    data_use : np.ndarray
        NormalizationsaretaExpressionmatrix (cells × genes)
    cell_labels : pd.Series
        Cell typelabel
    permutation : np.ndarray
        placechangetesutoforofindexdistributecolumn (cells × nboot)
    cell_types : list
        validnaCell typeoflist
    FunMean : function
        meanCalculationrelatedcount
    nboot : int
        placechangetesutotimecount
    n_jobs : int
        aligncolumnprocesswithusedokoacount
        
    Returns
    -------
    np.ndarray
        thingbeforeCalculationsaretameanExpressiondistributecolumn (genes × clusters × nboot)
    list
        simplepureizesaretaplacechangekurasumeanakusesuforoflist
    """
    logger.info("allGene×allCluster×allplacechangeofmeanExpressionthethingbeforeCalculationin...")
    
    # FunMeanfromCalculationtaiputhedecideset（koofpartdivissimpleabbreviateize）
  #  fun_type = 'triMean'  # deforuto
  #  trim = 0.1  # deforutovalue
    
    numCluster = len(cell_types)
    n_genes = data_use.shape[1]
    
    # memoriuseamountthecontrolcontroldoforofbachiprocess
    batch_size = min(10, nboot)  # onedegreetoprocessdopermutationofcounttheadjustarrange
    
    # allGene×allCluster×allpermutationofmeanExpressiontheformatstoredodistributecolumn
    all_gene_expr = np.zeros((n_genes, numCluster, nboot), dtype=np.float32)
    
    # bachitodivratioandaligncolumnprocess
    all_batches = [list(range(i, min(i+batch_size, nboot))) for i in range(0, nboot, batch_size)]
    
    # koacountthecontrollimit
    n_jobs_effective = min(n_jobs, len(all_batches), 32)  # maximum8koatocontrollimit
    
    if n_jobs_effective > 1:
        # aligncolumnprocesswithCalculation
        logger.info(f"{n_jobs_effective}koawithaligncolumnprocesstheRunin...")
        results = Parallel(n_jobs=n_jobs_effective)(
            delayed(process_permutation_batch)(
                batch_indices, data_use, cell_labels, permutation, cell_types, fun_type, trim
            ) for batch_indices in all_batches
        )
        
        # Resultofunifymatch
        for batch_idx, batch_indices in enumerate(all_batches):
            for idx, j in enumerate(batch_indices):
                if j < nboot:  # rangesurroundchieku
                    all_gene_expr[:, :, j] = results[batch_idx][:, :, idx]
    else:
        # shingurukoawithprocess
        logger.info("shingurukoawithprocesstheRunin...")
        for j in range(nboot):
            if j % 10 == 0:
                logger.info(f"Permutation {j+1}/{nboot} theProcessing...")
            
            # jnumbereyeofPermutationafterofCell type
            group_boot = cell_labels.values[permutation[:, j]]
            
            # eachClustertobelongdoCellofindextheget
            cluster_indices = [np.where(group_boot == ct)[0] for ct in cell_types]
            
            # allofGeneandCell typeofmeanExpressiontheCalculation
            all_gene_expr[:, :, j] = process_single_permutation(
                data_use, cluster_indices, cell_types, fun_type, trim
            )
    
    logger.info("All mean gene expression is calclulated.")
    
    return all_gene_expr


def precompute_complex_mapping(complex_input, gene_to_index):
    """
    complexandsoofstructbecomeGeneofindexmapinguthethingbeforeCalculation
    """
    complex_mapping = {}
    
    if complex_input.empty:
        return complex_mapping
    
    for complex_name in complex_input.index:
        # complexofsubunitotheget
        subunits_cols = [col for col in complex_input.columns if 'subunit' in col]
        subunits = complex_input.loc[complex_name, subunits_cols].dropna().astype(str)
        subunits = [s for s in subunits if s != "" and s in gene_to_index]
        
        if subunits:
            # subunitoofGeneindextheSave
            complex_mapping[complex_name] = [gene_to_index[s] for s in subunits]
    
    return complex_mapping

import pandas as pd

def check_gene_symbol(gene_set, gene_info):
    """
    gene_info ofindextoincludednotGeneisareba、WarningtheOutputdo（modifycorrectisrowwanot）
    """
    missing = [gene for gene in gene_set if gene not in gene_info.index]
    if missing:
        print("Warning: The following genes are not in geneInfo:", missing)

def extract_gene_subset(genes, complex_input, gene_info):
    """
    complexData（complex_input）toexistatdocase、eachGeneofrowfromsubunitoinfotheextractoutdo。
    existatshinakerebasourceofGenenamethesoofmamareturnsu。
    """
    extracted = []
    for gene in genes:
        if gene in complex_input.index:
            row = complex_input.loc[gene]
            # NaN yaemptytextcharcolumntheremoveremoveandsubunitotheextractout
            subunits = [str(x).strip() for x in row if pd.notna(x) and str(x).strip() != ""]
            if subunits:
                extracted.extend(subunits)
            else:
                extracted.append(gene)
        else:
            extracted.append(gene)
    return list(set(extracted))

def extractGene(db):
    """
    CellChatDB fromrelatedgiveandexistGeneshinborutheextractoutdorelatedcount
    （Rversion extractGene ofraisemovetomatchwaseru）

    Parameters
    ----------
    db : dict
        orbelowofkeytheincludemuwordwrite:
          - "interaction": pandas DataFrame（mustmust）→ column "ligand", "receptor", "agonist", "antagonist", "co_A_receptor", "co_I_receptor theincludemu
          - "complex": pandas DataFrame（rownametocomplexname、eachcolumntosubunitoinfo）
          - "cofactor": pandas DataFrame（rownameiskofuakuta-Gene、"cofactor" withstartmarucolumnissubunitoinfo）
          - "geneInfo": pandas DataFrame（publicformatnaGeneshinborutheindexandandkeephold）
          
    Returns
    -------
    list
        weightmultitheremoveita、relatedgiveandexistuni-kunaGeneshinboruoflist
    """
    # getwaymethodisyu-za-sidewithorbelowofliketorowteimasu
    # resource = db.get('interaction', pd.DataFrame())
    # complex_input = db.get('complex', pd.DataFrame())
    # cofactor_input = db.get('cofactor', pd.DataFrame())
    # gene_info = db.get('geneInfo', pd.DataFrame())
    resource = db.get('interaction', pd.DataFrame())
    complex_input = db.get('complex', pd.DataFrame())
    cofactor_input = db.get('cofactor', pd.DataFrame())
    gene_info = db.get('geneInfo', pd.DataFrame())

    # complexandkofuakuta-inofallneedelementtotsuitepublicformatGeneshinborukathechieku（Warningtheoutsu）
    if not complex_input.empty:
        complex_genes = complex_input.values.flatten()
        check_gene_symbol(list(complex_genes), gene_info)
    if not cofactor_input.empty:
        cofactor_genes = cofactor_input.values.flatten()
        check_gene_symbol(list(cofactor_genes), gene_info)

    # resource from ligand and receptor ofuni-kunaGenetheget
    geneL = resource["ligand"].dropna().astype(str).unique().tolist() if "ligand" in resource.columns else []
    geneR = resource["receptor"].dropna().astype(str).unique().tolist() if "receptor" in resource.columns else []
    geneLR = geneL + geneR

    # resource inwithcomplexDatatoincludednotGenetotsuitechieku（Warningofmi）
    if not complex_input.empty:
        genes_not_in_complex = [gene for gene in geneLR if gene not in complex_input.index]
    else:
        genes_not_in_complex = geneLR
    check_gene_symbol(genes_not_in_complex, gene_info)

    # complexinfotobaseduite、ligand, receptor ofeachlistthefurthernew
    if not complex_input.empty:
        geneL = extract_gene_subset(geneL, complex_input, gene_info)
        geneR = extract_gene_subset(geneR, complex_input, gene_info)
    geneLR = geneL + geneR

    # resource from agonist, antagonist, co_A_receptor, co_I_receptor columnthegetshi、kofuakuta-oflistthecreate
    cofactor_list = []
    for col in ["agonist", "antagonist", "co_A_receptor", "co_I_receptor"]:
        if col in resource.columns:
            cofactor_list.extend(resource[col].dropna().astype(str).tolist())
    cofactor_list = list(set([g for g in cofactor_list if g != ""]))

    # cofactor_input from、rownameis cofactor_list toonecausedorowof "cofactor" withstartmarucolumntheextractout
    gene_cofactor = []
    if not cofactor_input.empty:
        cofactor_cols = [col for col in cofactor_input.columns if col.startswith("cofactor")]
        matched_rows = cofactor_input.loc[cofactor_input.index.intersection(cofactor_list)]
        if not matched_rows.empty and cofactor_cols:
            cofactorsubunits = matched_rows[cofactor_cols]
            cofactorsubunits_v = cofactorsubunits.values.flatten()
            gene_cofactor = list(set([str(x).strip() for x in cofactorsubunits_v if pd.notna(x) and str(x).strip() != ""]))

    # finalalto、interaction reasoncomeofGeneand cofactor reasoncomeofGenetheunifymatchanduni-kunalistthereturnsu
    gene_use = list(set(geneLR + gene_cofactor))
    return gene_use, resource, complex_input, cofactor_input, gene_info

def identify_overexpressed_interactions(features_sig, gene_use, resource, complex_input=None):
    """
    RversionofidentifyOverExpressedInteractionstocompletealltoonecausesaseruimplement
    """
    # typeofNormalization
    features_sig = [str(f) for f in features_sig]
    gene_use = [str(g) for g in gene_use]
    
    # debaguOutput
    print(f"DEG genes count: {len(features_sig)}")
    
    # 1. complexofFiltering - Rversionandcompletealltosamerojiku
    expressed_complex_names = []
    complexes_with_subunits = {}
    
    if complex_input is not None and not complex_input.empty:
        for complex_name in complex_input.index:
            # Rversionandsamesubunitogetrojiku
            subunit_cols = [col for col in complex_input.columns if "subunit" in col]
            subunits = [str(s) for s in complex_input.loc[complex_name, subunit_cols].values if isinstance(s, str) and s != ""]
            
            # Rversionofrojiku：subunitoofizurekaisDEGwith、allteisAnalysispairobjectGenetoincludemareru
            has_deg_subunit = any(s in features_sig for s in subunits)
            all_in_gene_use = all(s in gene_use for s in subunits)
            
            if has_deg_subunit and all_in_gene_use:
                expressed_complex_names.append(complex_name)
                complexes_with_subunits[complex_name] = subunits
    
    print(f"Complex subunits with DEGs count: {len(expressed_complex_names)}")
    
    # 2. LRpairofFiltering
    valid_elements = set(features_sig).union(set(expressed_complex_names))
    valid_pairs = []
    
    for idx, row in resource.iterrows():
        ligand = str(row['ligand'])
        receptor = str(row['receptor'])
        
        # Rversionof all(unlist(pairLR[x,]) %in% c(features.sig, rownames(complexSubunits.sig)))
        if ligand in valid_elements and receptor in valid_elements:
            valid_pairs.append(idx)
    
    # FilteringsaretaLRpair
    resource_filtered = resource.loc[valid_pairs].copy() if valid_pairs else resource.iloc[0:0].copy()
    print(f"Filtered LR pairs: {len(resource_filtered)}")
    
    # ResulttheCSVtoekusupo-to（RResultandComparisonfor）
    pd.Series(features_sig).to_csv("py_features_sig.csv")
    pd.Series(expressed_complex_names).to_csv("py_complex_names.csv")
    resource_filtered.to_csv("py_filtered_pairs.csv")
    
    # 3. relatedconnectGeneofcollectgather
    lr_related_genes = set(features_sig)
    
    for _, row in resource_filtered.iterrows():
        ligand = str(row['ligand'])
        receptor = str(row['receptor'])
        
        if ligand in gene_use:
            lr_related_genes.add(ligand)
        if receptor in gene_use:
            lr_related_genes.add(receptor)
        
        if ligand in complexes_with_subunits:
            lr_related_genes.update(complexes_with_subunits[ligand])
        if receptor in complexes_with_subunits:
            lr_related_genes.update(complexes_with_subunits[receptor])
    
    return resource_filtered, list(lr_related_genes)

def debug_cellchatdb(db):
    """CellChatDBofincontentthedetailstoOutputandverify"""
    print("==== CellChatDBverify ====")
    
    if 'interaction' in db:
        interaction = db['interaction']
        print(f"Ligand-Receptorpaircount: {len(interaction)}")
        if len(interaction) > 0:
            print("firstof5pair:")
            cols_to_show = ['ligand', 'receptor']
            if 'interaction_name' in interaction.columns:
                cols_to_show.insert(0, 'interaction_name')
            if 'pathway_name' in interaction.columns:
                cols_to_show.append('pathway_name')
            print(interaction[cols_to_show].head(5))
            
            # Ligand,ReceptorofkindtypetheConfirm
            ligands = interaction['ligand'].unique()
            receptors = interaction['receptor'].unique()
            
            print(f"\nuni-kunaLigandcount: {len(ligands)}")
            print(f"LigandofExample (firstof10): {', '.join(map(str, ligands[:10]))}")
            
            print(f"\nuni-kunaReceptorcount: {len(receptors)}")
            print(f"ReceptorofExample (firstof10): {', '.join(map(str, receptors[:10]))}")
            
            # Ligand,ReceptoroftypetheConfirm
            print(f"\nLigandofDatatype: {type(interaction['ligand'].iloc[0])}")
            print(f"ReceptorofDatatype: {type(interaction['receptor'].iloc[0])}")
            
            # complexofLigand,ReceptorofConfirm
            complex_ligands = [l for l in ligands if isinstance(l, str) and '_' in l]
            complex_receptors = [r for r in receptors if isinstance(r, str) and '_' in r]
            
            print(f"\ncomplexLigandcount: {len(complex_ligands)}")
            if complex_ligands:
                print(f"complexLigandofExample: {', '.join(complex_ligands[:5])}")
                
            print(f"complexReceptorcount: {len(complex_receptors)}")
            if complex_receptors:
                print(f"complexReceptorofExample: {', '.join(complex_receptors[:5])}")

            # Databe-sugetafter、structmakethedetailstoConfirm
            print("DBstructmake:")
            for key in db:
                print(f"{key}: {type(db[key])}, rowcount: {len(db[key])}")

            print("Ligandoffirstof5tsuofExample:")
            print(interaction['ligand'].head(5).tolist())
            print("LigandofDatatype:", interaction['ligand'].dtype)
    else:
        print("CellChatDBtointeractionisincludedmasen")
    
    # complex, cofactorofConfirm
    if 'complex' in db and not db['complex'].empty:
        print(f"\ncomplexinfocount: {len(db['complex'])}")
        print("complexofExample (firstof3):")
        print(db['complex'].head(3))
    else:
        print("\ncomplexinfoisincludednotkaemptywithsu")
        
    if 'cofactor' in db and not db['cofactor'].empty:
        print(f"\ncofactorinfocount: {len(db['cofactor'])}")
        print("cofactorofExample (firstof3):")
        print(db['cofactor'].head(3))
    else:
        print("\ncofactorinfoisincludednotkaemptywithsu")
        
    print("==== CellChatDBverifyComplete ====\n")
    
    return db  # sourceofDBthesoofmamareturnandsourceofprocesstheinheritcontinuepossibleto



def cellchat_analysis(
    adata,
    groupby,
    #db,
    gene_use, #koreisLR genesofinfo
    complex_input,
    cofactor_input,
    resource,
    use_layer=None,
    min_cells=10,
    min_cells_expr=10,
    expr_prop=0.1,
    pseudocount=1.0,
    trim_threshold=0.05,
    k=0.5,
    n=1,
    fun_type="triMean",
    raw_use=True,
    population_size=False,
    nboot=100,
    seed=1,
    n_jobs=8,
    key_added="cellchat_res",
   # debug_mode=False,
    trim=0.1,
    apply_pval_filter=True,  # deforutotheTruetochangefurther
    features=None,  # add：specificofGenelisttheusedoOption
    r_patcher=False,
    do_de=True
):
    """
    mostfitizesaretaCellChatAlgorithm
    
    Parameters
    ----------
    adata : AnnData
        AnnDataobujiekuto
    groupby : str
        Cell type/Clustertheincludemu.obsofcolormuname
    db : dict
        CellChatDBwordwrite（get_cellchatdb_from_r()withgetshitamoof）
        'interaction', 'complex', 'cofactor', 'geneInfo'theincludemu
    
    [soofotherofParameterisomitabbreviate]
    apply_pval_filter : bool, optional
        SignificantwaterleveltofulltanotInteractionthedeleteremovedowhether。Default: True（Rofimplementandonecause）

    features : list, optional
        specificofGenelisttheusedocasetopointset。Noneofcaseisoverexpressed_genesfromspecific。
    
    Returns
    -------
    dict
        CellbetweencommunicationAnalysisResult
    """
    try:
        logger.info("CellChatAlgorithmtheRunin...")

        progress_bar = st.progress(0)
        status_area = st.empty()
        status_area.text("Preparing data...")

        # shi-dotheSettings
        np.random.seed(seed)


        # Preprocessingsuteputheadd
        if features is not None:
            st.write("Using union/intersection genes...")
            # provideprovidesaretaGenelisttheuse
            logger.info(f"Using provided gene list: {len(features)} genes")
            # existatdoGeneofmithekeephold

            # newshikusetmeaningshitarelatedcounttheuseandLRpairandrelatedconnectGenetheFiltering

            adata_filtered, resource_filtered = preprocess_data(adata, groupby, complex_input, gene_use=gene_use, min_cells=min_cells, min_cells_expr=min_cells_expr,
                thresh_pct=expr_prop, resource=resource, features =features, do_de=do_de)
            
            if not resource_filtered.empty:
                st.write(resource_filtered.head(3))
        else:
            st.write("Features is None")
            # followcomepassrioverexpressed_genesthespecific
            adata_filtered, resource_filtered = preprocess_data(adata, groupby, complex_input, gene_use=gene_use, min_cells=min_cells, min_cells_expr=min_cells_expr,
                thresh_pct=expr_prop, resource=resource, features =None, do_de=do_de)
        
        
        #after、reource_filterdtheresourceanddo
        resource=resource_filtered.copy()
        progress_bar.progress(0.05)
        
       
        # Datashapeformatofverify
        if not hasattr(adata_filtered, 'obs') or not hasattr(adata_filtered, 'var'):
            raise ValueError("InputisAnnDataobujiekutowithisnot exist")
            
        if groupby not in adata_filtered.obs.columns:
            raise ValueError(f"pointsetsaretaGroupcolumn '{groupby}' isadata.obstoexistatshimasen")
        
        # ExpressionDataofget
        if use_layer is not None and use_layer in adata_filtered.layers:
            logger.info(f"reiya- '{use_layer}' theuseshimasu")
            X = adata_filtered.layers[use_layer]
        else:
            logger.info("deforutoXmatorikusutheuseshimasu")
            X = adata_filtered.X
        
        # supa-sumatrixofconvert（mostfitize）
        X = optimize_matrix_operations(X, adata_filtered)

        progress_bar.progress(0.1)
        
        # ExpressionmatrixtheDataFrametoconvert
        expr_df = pd.DataFrame(X, index=adata_filtered.obs_names, columns=adata_filtered.var_names)
        
        # Cell typelabelofgetandverify
        cell_labels = adata_filtered.obs[groupby].copy()
        cell_types = np.array(sorted(cell_labels.unique()))
        
        logger.info(f"Cell typecount: {len(cell_types)}")
        if len(cell_types) < 2:
            raise ValueError(f"Cell typeis1tsushikanot exist。fewnakuandmo2tsurequiredwithsu。")
        
        # RofimplementandsamemannertomaximumvaluewithNormalization（R: data.use <- data/max(data)）
        data_use = X / np.max(X).astype(np.float64)
        
        # debagufor：DataofStatisticalinfotheOutput
        print(f"Data statistics after max normalization - min: {np.min(data_use)}, max: {np.max(data_use)}, mean: {np.mean(data_use)}")

        nC = data_use.shape[0]
        progress_bar.progress(0.2)
        status_area.text("Filtering data...")
        # FunMeanofSelect - Rofimplementandonecausesaseru
        if fun_type == "triMean":
            def FunMean(x):
                # NANtheclearshowaltoremoveout（na.rm=TRUEtopairrespond）
                x_no_nan = x[~np.isnan(x)]
                if len(x_no_nan) == 0:
                    return np.nan
                # Randcompletealltoonecausesaserufor、Rof quantile relatedcountthedirectconnectuse
                try:
                    import rpy2.robjects as ro
                    from rpy2.robjects import numpy2ri
                    numpy2ri.activate()
                    r_quantile = ro.r['quantile']
                    r_mean = ro.r['mean']
                    
                    # RwithtriMeantheCalculation
                    r_vector = ro.FloatVector(x_no_nan)
                    r_result = r_mean(r_quantile(r_vector, ro.FloatVector([0.25, 0.5, 0.5, 0.75])))
                    return float(r_result[0])
                except:
                    # fo-rubaku: Pythonofimplement
                    q1 = np.percentile(x_no_nan, 25, interpolation='linear')
                    q2 = np.percentile(x_no_nan, 50, interpolation='linear')
                    q3 = np.percentile(x_no_nan, 75, interpolation='linear')
                    return (q1 + 2*q2 + q3) / 4
        elif fun_type == "truncatedMean":
            from scipy.stats import trim_mean
            def FunMean(x):
                return trim_mean(x, proportiontocut=trim)
        elif fun_type == "median":
            def FunMean(x):
                return np.median(x)
        else:
            def FunMean(x):
                return np.mean(x)
        
        # eachGeneofeachCell typetookerumeanExpressionamounttheCalculation（mostfitize）
        data_use_avg_dict, cell_counts = calculate_mean_expression_optimized(
            data_use, cell_labels, cell_types, min_cells, FunMean
        )
        
        logger.info(f"Cell typegoandofCellcount: {cell_counts}")
        
        # validnaCell typetheConfirm
      #  if len(cell_types) < 2:
      #      raise ValueError(f"minimumCellcount {min_cells} thefulltasuCell typeis2tsunot yetfullwithsu。min_cellsthereduceraandplease。")
        progress_bar.progress(0.3)
        status_area.text("Calculating LR expression levels...")
        # meanExpressionamounttheDataFrametoconvert
        data_use_avg_df = pd.DataFrame(data_use_avg_dict, index=adata_filtered.var_names)
        
        # Genetheindextomapingu
        gene_to_index = {gene: i for i, gene in enumerate(adata_filtered.var_names)}
        
        # LigandandReceptorofmeanExpressionamounttheCalculation
        logger.info("LigandandReceptorofExpressionamounttheCalculationin...")
        dataLavg = computeExpr_LR(resource['ligand'].values, data_use_avg_df, complex_input)
        dataRavg = computeExpr_LR(resource['receptor'].values, data_use_avg_df, complex_input)
        
        # debagufor：Expressionvaluetheekusupo-to
        pd.DataFrame(dataLavg, columns=cell_types, index=resource.index).to_csv("py_dataLavg.csv")
        pd.DataFrame(dataRavg, columns=cell_types, index=resource.index).to_csv("py_dataRavg.csv")
        
        # togetheractivenatureizeandtogetherblockharmReceptorofeffectresulttheconsiderconsider
        dataRavg_co_A_receptor = computeExpr_coreceptor(cofactor_input, data_use_avg_df, resource, "A")
        dataRavg_co_I_receptor = computeExpr_coreceptor(cofactor_input, data_use_avg_df, resource, "I")
        dataRavg = dataRavg * dataRavg_co_A_receptor / dataRavg_co_I_receptor
        
        # Cellcountofeffectresulttheconsiderconsider
        if population_size:
            # eachCell typeofCellcountofratiomatchtheCalculation
            cell_proportions = np.array([np.sum(cell_labels == ct) for ct in cell_types]) / nC
            # eachLigand,Receptorpairtopairandsameratiomatchtheuse
            dataLavg2 = np.tile(cell_proportions, (len(resource), 1))
            dataRavg2 = dataLavg2
        else:
            dataLavg2 = np.ones((len(resource), len(cell_types)))
            dataRavg2 = np.ones((len(resource), len(cell_types)))
        
        # agonisutoandantagonisutoofindexthespecific
        index_agonist = np.where(resource['agonist'].notna() & (resource['agonist'] != ""))[0] if 'agonist' in resource.columns else []
        index_antagonist = np.where(resource['antagonist'].notna() & (resource['antagonist'] != ""))[0] if 'antagonist' in resource.columns else []
        
        progress_bar.progress(0.4)
        status_area.text("Preparing permutation data...")
        # placechangechecksetofforofDatathelevelprepare
        permutation = np.zeros((nC, nboot), dtype=int)

        # RandofcompleteallonecauseoffortoRofrandomcountgeneratetheuse
        try:
            permutation = get_r_permutation(nC, nboot, seed=seed)
            print("Using R permutation for exact compatibility")
        except Exception as e:
            print(f"R permutation failed, using Python fallback: {e}")
            for i in range(nboot):
                permutation[:, i] = np.random.permutation(nC)
        progress_bar.progress(0.5)
        
        # allGeneofmeanExpressionthethingbeforeCalculation（mostfitize）
        #all_gene_expr, data_use_avg_boot = precompute_gene_expressions(
        all_gene_expr = precompute_gene_expressions(
            data_use, cell_labels, permutation, cell_types, FunMean, nboot, fun_type=fun_type, trim=trim
        )
        
        # komiyunike-tionprobabilityandSignificantnaturetheCalculation
        numCluster = len(cell_types)
        nLR = len(resource)
        Prob = np.zeros((numCluster, numCluster, nLR))
        Pval = np.zeros((numCluster, numCluster, nLR))
        
        logger.info(f"Ligand-ReceptorpairofAnalysisthestart: {nLR}pair")
        progress_bar.progress(0.6)
        status_area.text("Starting LR pair analysis...")
        
        # complexandsoofstructbecomeGeneofmapinguthethingbeforetoCalculation
        complex_mapping = precompute_complex_mapping(complex_input, gene_to_index)
        
        # Ligand,ReceptorofGeneindextheget
        ligand_indices = []
        receptor_indices = []
        
        nLR = len(resource)
        for i in range(nLR):
            ligand = resource['ligand'].iloc[i]
            receptor = resource['receptor'].iloc[i]
            
            # simpleoneGenekacomplexkathejudgeset
            if isinstance(ligand, str) and ligand in gene_to_index:
                # simpleoneGene
                ligand_indices.append((i, [gene_to_index[ligand]], False))
            elif isinstance(ligand, str) and ligand in complex_mapping:
                # complex
                ligand_indices.append((i, complex_mapping[ligand], True))
            else:
                # not yetknowofGene/complex
                ligand_indices.append((i, [], None))
            
            if isinstance(receptor, str) and receptor in gene_to_index:
                receptor_indices.append((i, [gene_to_index[receptor]], False))
            elif isinstance(receptor, str) and receptor in complex_mapping:
                receptor_indices.append((i, complex_mapping[receptor], True))
            else:
                receptor_indices.append((i, [], None))
        
        # LRpairgoandofru-pu
        st.write("Long loop processes...")
        progress_bar_lr = st.progress(0)
        for i in range(nLR):
            dataLR = np.outer(dataLavg[i, :], dataRavg[i, :])
            
            # hirurelatedcounttheuse（Rofko-doandonecause）
            P1 = compute_hill_outer_vectorized(dataLavg[i, :], dataRavg[i, :], k, n)
            
            # agonisutoeffectresulttheCalculation
            P2 = np.ones((numCluster, numCluster))
            if i in index_agonist:
                data_agonist = computeExpr_agonist(data_use_avg_df, resource, cofactor_input, i, k, n)
                P2 = np.outer(data_agonist, data_agonist)
            
            # antagonisutoeffectresulttheCalculation
            P3 = np.ones((numCluster, numCluster))
            if i in index_antagonist:
                data_antagonist = computeExpr_antagonist(data_use_avg_df, resource, cofactor_input, i, k, n)
                P3 = np.outer(data_antagonist, data_antagonist)
            
            # CellcounteffectresulttheCalculation
            P4 = np.ones((numCluster, numCluster))
            if population_size:
                P4 = np.outer(dataLavg2[i, :], dataRavg2[i, :])
            
            # finalalnaprobability
            Pnull = P1 * P2 * P3 * P4
            Prob[:, :, i] = Pnull
            
            # placechangechecksettoyoruPvalueCalculation
            if np.sum(Pnull) == 0:
                # Interactionisnotcase
                Pval[:, :, i] = 1
                continue
            
            Pnull_vec = Pnull.flatten()
            
            # koofLRpairofLigandandReceptorofindexinfotheget
            ligand_info = ligand_indices[i]
            receptor_info = receptor_indices[i]
            
            # Expressionisgetwithkinotcaseissukipu
            if ligand_info[2] is None or receptor_info[2] is None:
                Pval[:, :, i] = 1
                continue
            
            # eachplacechangewithofbu-tosutorapuprobabilitytheCalculation
            Pboot = np.zeros((numCluster * numCluster, nboot))
            
            # bachiprocesstheguideenterandmemoriuseamountthemanagereason
            batch_size = min(20, nboot)
            
            # aligncolumnprocessofSettings
            n_jobs_to_use = min(n_jobs, os.cpu_count() or 1)
            
            # aligncolumnprocessforrelatedcount

            def compute_permutation_batch_vectorized(batch_indices):
                batch_results = np.zeros((numCluster * numCluster, len(batch_indices)))
                for idx, j in enumerate(batch_indices):
                    lr_i, ligand_gene_indices, is_ligand_complex = ligand_info
                    # LigandofExpressionvalueget
                    if not is_ligand_complex:
                        if ligand_gene_indices:
                            ligand_idx = ligand_gene_indices[0]
                            dataLavgB = all_gene_expr[ligand_idx, :, j].reshape(1, -1)
                        else:
                            dataLavgB = np.zeros((1, numCluster))
                    else:
                        expr_values = np.array([all_gene_expr[l_idx, :, j] for l_idx in ligand_gene_indices])
                        if expr_values.size > 0:
                            log_values = np.log(expr_values + 1e-10)
                            dataLavgB = np.exp(np.mean(log_values, axis=0)).reshape(1, -1)
                        else:
                            dataLavgB = np.zeros((1, numCluster))
                    # ReceptorofExpressionvalueget
                    # ReceptorofExpressionget
                    lr_i, receptor_gene_indices, is_receptor_complex = receptor_info
                    if not is_receptor_complex:
                        if receptor_gene_indices:
                            receptor_idx = receptor_gene_indices[0]
                            dataRavgB = all_gene_expr[receptor_idx, :, j].reshape(1, -1)
                        else:
                            dataRavgB = np.zeros((1, numCluster))
                    else:
                        expr_values = np.array([all_gene_expr[r_idx, :, j] for r_idx in receptor_gene_indices])
                        if expr_values.size > 0:
                            log_values = np.log(expr_values + 1e-10)
                            dataRavgB = np.exp(np.mean(log_values, axis=0)).reshape(1, -1)
                        else:
                            dataRavgB = np.zeros((1, numCluster))
                    # outaccumulatetheCalculationshi、bekutoruizehirurelatedcountthefitfor
                    dataLRB = np.outer(dataLavgB[0, :], dataRavgB[0, :])
                    P1_boot = hill_function(dataLRB, k, n)
                    batch_results[:, idx] = P1_boot.flatten()
                return batch_results
            # bachiprocesswithaligncolumnCalculation
            for b_start in range(0, nboot, batch_size):
                b_end = min(b_start + batch_size, nboot)
                batch_indices = list(range(b_start, b_end))
                
                if n_jobs_to_use > 1 and len(batch_indices) > 1:
                    # aligncolumnprocess
                    batch_results_list = Parallel(n_jobs=n_jobs_to_use, backend="loky")(
                        delayed(compute_permutation_batch_vectorized)([j]) for j in batch_indices
                    )
                    # Resultofresultmatch
                    for j_idx, j in enumerate(batch_indices):
                        Pboot[:, j] = batch_results_list[j_idx][:, 0]
                else:
                    # simpleonesuredoprocess
                    batch_results = compute_permutation_batch_vectorized(batch_indices)
                    for j_idx, j in enumerate(batch_indices):
                        Pboot[:, j] = batch_results[:, j_idx]
            
            # pvalueofCalculation
            #nReject = np.sum(Pboot > np.expand_dims(Pnull_vec, 1), axis=1)
            # RofimplementtomatchwasetediffdivwithComparison
            if r_patcher:
                # floatmovesmallcountpointprecisedegreeofquestiontopicthetimeavoid
                nReject = np.sum((Pboot - np.expand_dims(Pnull_vec, 1)) > 1e-10, axis=1)
            else:
                # Randsame: rowSums(Pboot - Pnull > 0)
                nReject = np.sum((Pboot - np.expand_dims(Pnull_vec, 1)) > 0, axis=1)
            p = nReject / nboot
            Pval[:, :, i] = p.reshape(numCluster, numCluster)
            progress_bar_lr.progress((i + 1) / nLR)
        

        progress_bar_lr.empty()
        progress_bar.progress(0.8)

        # purobuis0ofcaseofpvaluethe1toSettings (koreisnormaltorowu)
        Pval[Prob == 0] = 1

        # PvalueFilteringbeforeofstatestatetheverify
        print(f"CellChat_analysis Before pval filtering - Prob > 0 count: {np.sum(Prob > 0)}")
        
        # apply_pval_filterParametertobaseduitepvaluewithofFilteringthefitfor
        if apply_pval_filter:
            if r_patcher:
                #Prob[(Pval >= trim_threshold) > 1.49e-8] = 0
                Prob[Pval >= (trim_threshold - 1.49e-8)] = 0
            else:
                Prob[Pval >= trim_threshold] = 0
        # PvalueFilteringafterofstatestatetheverify
        print(f"After pval filtering - Prob > 0 count: {np.sum(Prob > 0)}")


        # Filter out communication for cell types with few cells
        cell_counts = {ct: np.sum(cell_labels == ct) for ct in cell_types}
        cell_excludes = [ct for ct, count in cell_counts.items() if count <= min_cells]

        if cell_excludes:
            print(f"The cell-cell communication related with the following cell groups are excluded due to the few number of cells: {cell_excludes}")
            
            # Get the indices of excluded cell types
            exclude_indices = [i for i, ct in enumerate(cell_types) if ct in cell_excludes]
            
            # Set entire rows and columns to zero for excluded cell types
            for idx in exclude_indices:
                Prob[idx, :, :] = 0  # Sender is excluded
                Prob[:, idx, :] = 0  # Receiver is excluded
            
            # Update p-values
            Pval[Prob == 0] = 1

        # ResulttonamebeforetheSettings
        dimnames = [list(cell_types), list(cell_types), list(resource.index)]
        
        # debagufor：InteractioncountofdetailstheOutput
        print(f"\n=== Interaction Count Summary ===")
        print(f"Total LR pairs analyzed: {nLR}")
        print(f"LR pairs with Prob > 0 (before p-value filter): {np.sum(np.any(Prob > 0, axis=(0,1)))}")
        print(f"Total interactions with Prob > 0: {np.sum(Prob > 0)}")
        print(f"Total interactions with p-value < {trim_threshold}: {np.sum(Pval < trim_threshold)}")
        print(f"Total significant interactions (Prob > 0 AND p < {trim_threshold}): {np.sum((Prob > 0) & (Pval < trim_threshold))}")
        
        # Cell typepairgoandofInteractioncount
        interaction_counts = np.zeros((len(cell_types), len(cell_types)))
        for i in range(len(cell_types)):
            for j in range(len(cell_types)):
                if apply_pval_filter:
                    interaction_counts[i, j] = np.sum((Prob[i, j, :] > 0) & (Pval[i, j, :] < trim_threshold))
                else:
                    interaction_counts[i, j] = np.sum(Prob[i, j, :] > 0)
        
        print(f"\nInteraction counts per cell type pair:")
        interaction_df = pd.DataFrame(interaction_counts, index=cell_types, columns=cell_types)
        print(interaction_df)
        print(f"================================\n")
        
        # signalPathwayreberuwithofcommunicationprobabilitytheCalculation
        netP = computeCommunProbPathway({"prob": Prob, "pval": Pval}, resource,
            thresh=trim_threshold, apply_pval_filter=apply_pval_filter, r_patcher=r_patcher)
        
        # incenternaturepointmarktheCalculation
        if netP["pathways"] is not None and len(netP["pathways"]) > 0:
            netP["centr"] = {}
            for p_idx in range(len(netP["pathways"])):
                # eachpathwayofprobabilitymatrixtheget
                pathway_prob = netP["prob"][:, :, p_idx]
                
                # simpleoneofpathwaytopairandincenternaturetheCalculation
                # 3nextsourcedistributecolumnthecreateandCalculationrelatedcounttotransfersu
                pathway_prob_3d = np.expand_dims(pathway_prob, axis=2)
                netP["centr"][p_idx] = netAnalysis_computeCentrality(pathway_prob_3d)[0]
        
        # aggregateNetworkofincenternaturemoCalculation
        prob_sum = np.sum(Prob, axis=2)
        prob_sum_3d = np.expand_dims(prob_sum, axis=2)
        net_centr = netAnalysis_computeCentrality(prob_sum_3d)[0]
        
        # aggregateNetworktheCalculation
        net_summary = aggregateCell_Cell_Communication({"prob": Prob, "pval": Pval},
            cell_types, pval_threshold=trim_threshold, apply_pval_filter=apply_pval_filter,
            r_patcher=r_patcher)
        progress_bar.empty()
        logger.info("CellChatAnalysisisCompleteshimashita")
        
        # ResulttheDatafure-mutoconvert
        results_data = {
            'source': [],
            'target': [],
            'interaction_name': [],
            'ligand': [],
            'receptor': [],
            'prob': [],
            'pval': []
        }
        
        for i in range(len(cell_types)):
            for j in range(len(cell_types)):
                for k in range(nLR):
                    if Prob[i, j, k] > 0:
                        results_data['source'].append(cell_types[i])
                        results_data['target'].append(cell_types[j])
                        results_data['interaction_name'].append(resource.index[k])
                        results_data['ligand'].append(resource['ligand'].iloc[k])
                        results_data['receptor'].append(resource['receptor'].iloc[k])
                        results_data['prob'].append(Prob[i, j, k])
                        results_data['pval'].append(Pval[i, j, k])
        
        results_df = pd.DataFrame(results_data)
        
        return {
            'adata': adata_filtered,
            'results': results_df,
            'net': {
                "prob": Prob,
                "pval": Pval,
                "dimnames": dimnames,
                "centr": net_centr
            },
            'netP': netP,
            'network': net_summary,
            'groupby': groupby 
        }
    
    except Exception as e:
        logger.error(f"Errorisoccurred: {str(e)}")
        logger.error(traceback.format_exc())
        return {'error': str(e), 'traceback': traceback.format_exc()}


def geometricMean(expr_values):
    """
    RofimplementtomatchwasetahowwhatmeanCalculation
    """
    if expr_values.ndim == 1:
        # 1nextsourcedistributecolumnofcase
        # Randsamemovemake: log(0) = -Inf, mean with -Inf = -Inf, exp(-Inf) = 0
        with np.errstate(divide='ignore'):
            log_values = np.log(expr_values)
        # -Infisincludemarerucase、meanmo-Inftobecome（Rofmovemaketheagainpresent）
        if np.any(np.isneginf(log_values)):
            return 0.0
        else:
            return np.exp(np.mean(log_values))
    else:
        # 2nextsourcedistributecolumnofcase（columngoandtoCalculation）
        result = np.zeros(expr_values.shape[1])
        for i in range(expr_values.shape[1]):
            with np.errstate(divide='ignore'):
                log_values = np.log(expr_values[:, i])
            if np.any(np.isneginf(log_values)):
                result[i] = 0.0
            else:
                result[i] = np.exp(np.mean(log_values))
        return result

@st.cache_data
def computeExpr_coreceptor(cofactor_input, data_use, pairLRsig, type_coreceptor):
    """
    Ligand-ReceptorofInteractiontookerutogetherReceptoreffectresultthemoderuize
    """
    if cofactor_input.empty or pairLRsig.empty:
        return np.ones((len(pairLRsig), data_use.shape[1]))
    
    coreceptor_col = 'co_A_receptor' if type_coreceptor == "A" else 'co_I_receptor'
    
    if coreceptor_col not in pairLRsig.columns:
        return np.ones((len(pairLRsig), data_use.shape[1]))
    
    coreceptor_all = pairLRsig[coreceptor_col].values
    index_coreceptor = np.where((pd.notna(coreceptor_all)) & (coreceptor_all != ""))[0]
    numCluster = data_use.shape[1]
    data_coreceptor = np.ones((len(coreceptor_all), numCluster))
    
    if len(index_coreceptor) > 0:
        for idx in index_coreceptor:
            coreceptor = coreceptor_all[idx]
            if coreceptor in cofactor_input.index:
                # cofactortheget
                cofactor_cols = [col for col in cofactor_input.columns if 'cofactor' in col]
                cofactors = cofactor_input.loc[coreceptor, cofactor_cols].dropna().astype(str)
                # cofactors = [c for c in cofactors if c != "" and c in data_use.index]
                cofactors = [str(c) for c in cofactors if c != "" and str(c) in data_use.index]
                
                if len(cofactors) == 1:
                    data_coreceptor[idx] = 1 + data_use.loc[cofactors[0]].values
                elif len(cofactors) > 1:
                    # 1 + eachtogetherReceptorofExpressionofaccumulate
                    prod = np.ones(numCluster)
                    for c in cofactors:
                        prod *= (1 + data_use.loc[c].values)
                    data_coreceptor[idx] = prod
    
    return data_coreceptor

@st.cache_data
def computeExpr_agonist(data_use, pairLRsig, cofactor_input, index_agonist, Kh, n):
    """
    agonisutoisLigand-ReceptorInteractiontogiveerueffectresultthemoderuize
    """
    if cofactor_input.empty or pairLRsig.empty or 'agonist' not in pairLRsig.columns:
        return np.ones(data_use.shape[1])
    
    agonist = pairLRsig['agonist'].iloc[index_agonist]
    if pd.isna(agonist) or agonist == "" or agonist not in cofactor_input.index:
        return np.ones(data_use.shape[1])
    
    # agonisutoGenetheget
    cofactor_cols = [col for col in cofactor_input.columns if 'cofactor' in col]
    agonist_genes = cofactor_input.loc[agonist, cofactor_cols].dropna().astype(str)
    agonist_genes = [g for g in agonist_genes if g != "" and g in data_use.index]
    
    if len(agonist_genes) == 1:
        data_avg = data_use.loc[agonist_genes[0]].values
        data_agonist = 1 + data_avg**n / (Kh**n + data_avg**n)
    elif len(agonist_genes) > 1:
        # eachagonisutoGeneofeffectresultofaccumulate
        data_agonist = np.ones(data_use.shape[1])
        for g in agonist_genes:
            data_avg = data_use.loc[g].values
            data_agonist *= (1 + data_avg**n / (Kh**n + data_avg**n))
    else:
        data_agonist = np.ones(data_use.shape[1])
    
    return data_agonist

@st.cache_data
def computeExpr_antagonist(data_use, pairLRsig, cofactor_input, index_antagonist, Kh, n):
    """
    antagonisutoisLigand-ReceptorInteractiontogiveerueffectresultthemoderuize
    """
    if cofactor_input.empty or pairLRsig.empty or 'antagonist' not in pairLRsig.columns:
        return np.ones(data_use.shape[1])
    
    antagonist = pairLRsig['antagonist'].iloc[index_antagonist]
    if pd.isna(antagonist) or antagonist == "" or antagonist not in cofactor_input.index:
        return np.ones(data_use.shape[1])
    
    # antagonisutoGenetheget
    cofactor_cols = [col for col in cofactor_input.columns if 'cofactor' in col]
    antagonist_genes = cofactor_input.loc[antagonist, cofactor_cols].dropna().astype(str)
    antagonist_genes = [g for g in antagonist_genes if g != "" and g in data_use.index]
    
    if len(antagonist_genes) == 1:
        data_avg = data_use.loc[antagonist_genes[0]].values
        data_antagonist = Kh**n / (Kh**n + data_avg**n)
    elif len(antagonist_genes) > 1:
        # eachantagonisutoGeneofeffectresultofaccumulate
        data_antagonist = np.ones(data_use.shape[1])
        for g in antagonist_genes:
            data_avg = data_use.loc[g].values
            data_antagonist *= Kh**n / (Kh**n + data_avg**n)
    else:
        data_antagonist = np.ones(data_use.shape[1])
    
    return data_antagonist

@st.cache_data
def computeCommunProbPathway(net, pairLR_use, thresh=0.05, apply_pval_filter=True, r_patcher=False):
    """
    signalPathwayreberuwithofcommunicationprobabilitytheCalculation
    
    Parameters
    ----------
    net : dict
        communicationprobabilityandSignificantnaturetheincludemuwordwrite
    pairLR_use : pd.DataFrame
        Ligand-Receptorpairinfo
    thresh : float, optional
        SignificantnaturejudgesetofpvalueThreshold
    apply_pval_filter : bool, optional
        pvalueFilteringthefitfordowhether。deforutoisTrue
        
    Returns
    -------
    dict
        pathwayreberuwithofcommunicationprobabilityinfo
    """
    prob = net["prob"].copy()
    pval = net["pval"]
    
    # PvalueFilteringbeforeofstatestatetheverify
    print(f"computeCommunProbPathway - Before filtering - prob > 0 count: {np.sum(prob > 0)}")

    # apply_pval_filterParametertobaseduitepvaluewithofFilteringthefitfor
    if apply_pval_filter:
        if r_patcher:
           # prob[(pval - thresh) >= 1.49e-8] = 0
           prob[pval >= (thresh-1.49e-8)] = 0
        else:
            prob[pval >= thresh] = 0
    # PvalueFilteringafterofstatestatetheverify
    print(f"computeCommunProbPathway - After filtering - prob > 0 count: {np.sum(prob > 0)}")
    # pathwayinfoisnotcaseisprocessthesukipu
    if 'pathway_name' not in pairLR_use.columns:
        return {
            "pathways": [],
            "prob": np.zeros((prob.shape[0], prob.shape[1], 0))
        }
    
    # signalPathwaygoandtoaggregate
    pathways = pairLR_use['pathway_name'].dropna().unique()
    prob_pathways = np.zeros((prob.shape[0], prob.shape[1], len(pathways)))
    
    for i, pathway in enumerate(pathways):
        idx = np.where(pairLR_use['pathway_name'] == pathway)[0]
        prob_pathways[:, :, i] = np.sum(prob[:, :, idx], axis=2)
    
    # totalcommunicationstrengthtobaseduitealignbereplacee
    pathway_sums = np.sum(prob_pathways, axis=(0, 1))
    significant_idx = np.where(pathway_sums > 0)[0]
    
    if len(significant_idx) == 0:
        return {
            "pathways": [],
            "prob": np.zeros((prob.shape[0], prob.shape[1], 0))
        }
    
    sort_idx = significant_idx[np.argsort(-pathway_sums[significant_idx])]
    
    pathways_sig = pathways[sort_idx]
    prob_pathways_sig = prob_pathways[:, :, sort_idx]
    
    return {
        "pathways": pathways_sig,
        "prob": prob_pathways_sig
    }


@st.cache_data
def computeExpr_complex(complex_input, data_use, complex_genes):
    """complexofExpressiontheCalculationdorelatedcount（Rofimplementtopairrespond）"""
    result = np.zeros((len(complex_genes), data_use.shape[1]))
    
    for i, complex_gene in enumerate(complex_genes):
        if complex_gene in complex_input.index:
            # complexofsubunitotheget
            subunits_cols = [col for col in complex_input.columns if 'subunit' in col]
            subunits = complex_input.loc[complex_gene, subunits_cols].dropna().astype(str)
            subunits = [s for s in subunits if s != "" and s in data_use.index]
            
            if len(subunits) > 0:
                # Expressionvalueofget
                expr_values = data_use.loc[subunits].values
                # howwhatmeanofCalculation
                result[i] = geometricMean(expr_values)
    
    return result

@st.cache_data
def computeExpr_LR(geneLR, data_use, complex_input):
    """LigandalsoisReceptorofExpressiontheCalculationdorelatedcount（Rofimplementtopairrespond）"""
    geneLR = [str(gene) for gene in geneLR]  # clearshowaltotextcharcolumntoconvert
    print(f"First 5 geneLR values: {geneLR[:5]}")
    print(f"Data_use.index type: {type(data_use.index)}")
    print(f"First 5 gene names in data_use: {list(data_use.index)[:5]}")
    print(f"Gene name types in data_use: {[type(g) for g in data_use.index[:5]]}")

    nLR = len(geneLR)
    numCluster = data_use.shape[1]
    
    # simpleoneGeneofprocess（Randsamemannerofwaymethod）
    index_singleL = [i for i, gene in enumerate(geneLR) if gene in data_use.index]
    dataLavg = np.zeros((nLR, numCluster))
    
    if index_singleL:
        # simpleoneGeneofExpressionthemaandmeteget
        gene_indices = [geneLR[i] for i in index_singleL if geneLR[i] in data_use.index]
        if gene_indices:
            dataL1avg = data_use.loc[gene_indices].values
            for idx, gene_idx in enumerate(index_singleL):
                if idx < len(dataL1avg):
                    dataLavg[gene_idx] = dataL1avg[idx]
    
    # complexofprocess（Randsamemannerofwaymethod）
    index_complexL = [i for i in range(nLR) if i not in index_singleL]
    if index_complexL and not complex_input.empty:
        complex_genes = [geneLR[i] for i in index_complexL]
        data_complex = computeExpr_complex(complex_input, data_use, complex_genes)
        for idx, complex_idx in enumerate(index_complexL):
            dataLavg[complex_idx] = data_complex[idx]
    
    return dataLavg





@st.cache_data
def netAnalysis_computeCentrality(prob):
    """
    NetworkincenternaturepointmarktheCalculation (reformgoodversion)
    
    Parameters
    ----------
    prob : numpy.ndarray
        communicationprobabilitymatrix (shapestate: [cell_types, cell_types, interactions])
    debug_mode : bool, optional
        debagumo-do
        
    Returns
    -------
    centrality : dict
        eachintarakution（Network）goandofincenternaturepointmarkofdeikushiyonari
    """
    
    # Inputis2nextsourceofcaseis3nextsourcetoexpandstretch（axisErrortimeavoid）
    if prob.ndim == 2:
        prob = np.expand_dims(prob, axis=2)
    
    # rpy2ofSettingsandrequirednaRpake-jiofinpo-to
    r_available = False
    try:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr
        from rpy2.robjects import numpy2ri
        
        # selfmoveconvertofvalidize
        numpy2ri.activate()
        
        # Rraiburariofinpo-to
        base = importr('base')
        sna = importr('sna')
        
        # igraphpake-jiofLoading
        try:
            ro.r('library(igraph, quietly=TRUE)')
            r_available = True
            print("RraiburariofLoadingSuccess: igraph, sna")
        except Exception as e:
            print(f"igraphLoadingError: {str(e)}")
            r_available = False
            
    except Exception as e:
        print(f"RraiburariofLoadingError: {str(e)}")
        print(traceback.format_exc())
        print("NetworkXwithfo-rubakushimasu。")
    
    # NetworkXtheuseshitaCalculation
  #  centrality_nx = netAnalysis_computeCentrality_nx(prob)
    

    try:
        # Rofrelatedcountthesetmeaning
        ro.r('''
        computeCentralityLocal <- function(net) {
          centr <- list()
          G <- igraph::graph_from_adjacency_matrix(net, mode = "directed", weighted = TRUE)
          centr$outdeg_unweighted <- rowSums(net > 0)
          centr$indeg_unweighted <- colSums(net > 0)
          centr$outdeg <- igraph::strength(G, mode="out")
          centr$indeg <- igraph::strength(G, mode="in")
          centr$hub <- igraph::hub_score(G)$vector
          centr$authority <- igraph::authority_score(G)$vector
          centr$eigen <- igraph::eigen_centrality(G)$vector
          centr$page_rank <- igraph::page_rank(G)$vector
          igraph::E(G)$weight <- 1/igraph::E(G)$weight
          centr$betweenness <- igraph::betweenness(G)
          centr$flowbet <- tryCatch({
            sna::flowbet(net)
          }, error = function(e) {
            rep(0, nrow(net))
          })
          centr$info <- tryCatch({
            sna::infocent(net, diag = TRUE, rescale = TRUE, cmode = "lower")
          }, error = function(e) {
            rep(0, nrow(net))
          })
          return(centr)
        }
        ''')

        centrality_r = {}
        for i in range(prob.shape[2]):
            try:
                net_mat = prob[:, :, i]
                
                # NumPymatrixtheRofmatrixtoconvert
                r_mat = ro.r.matrix(net_mat, nrow=net_mat.shape[0], ncol=net_mat.shape[1])
                
                # RrelatedcounttheRun
                r_centr = ro.r['computeCentralityLocal'](r_mat)
                
                # RofResultthePythonofwordwritetoconvert
                centrality_r[i] = {}
                for metric in ["outdeg", "indeg", "outdeg_unweighted", "indeg_unweighted", 
                              "hub", "authority", "eigen", "page_rank", "betweenness", "flowbet", "info"]:
                    try:
                        if metric in r_centr.names:
                            centrality_r[i][metric] = np.array(r_centr.rx2(metric))
                    except Exception as e:
                        print(f"pointmark {metric} ofgetintoError: {str(e)}")
                        centrality_r[i][metric] = np.zeros(net_mat.shape[0])  # deforutovaluetheSettings
            
            except Exception as e:
                print(f"RCalculationwithError (interaction {i}): {str(e)}")
                centrality_r[i] = centrality_nx[i]  # NetworkXofResulttheuse
        
        # debagumo-doofcase、ResulttheComparison
  #      if debug_mode:
  #          compare_centrality_results(centrality_r, centrality_nx)
        
        return centrality_r
    except Exception as e:
        print(f"RwithofincenternatureCalculationtolosefailshimashita: {str(e)}")
        print(traceback.format_exc())
        print("NetworkXofResulttheuseshimasu。")

@st.cache_data
def aggregateCell_Cell_Communication(net, cell_types, pval_threshold=0.05, apply_pval_filter=True, r_patcher=False):
    """
    CellbetweencommunicationNetworktheaggregate
    
    Parameters
    ----------
    net : dict
        communicationprobabilityandSignificantnaturetheincludemuwordwrite
    cell_types : list
        validnaCell typeoflist
    pval_threshold : float, optional
        SignificantnaturejudgesetofpvalueThreshold。deforutois0.05
    apply_pval_filter : bool, optional
        pvalueFilteringthefitfordowhether。deforutoisTrue
        
    Returns
    -------
    dict
        aggregatesaretaNetworkinfo
    """
    prob = net["prob"]
    pval = net["pval"]
    
    # PythonimplementwithFilteringbeforeafterofvaluetheComparison
    print(f"aggregateCell_Cell_Communication Before filtering - sig_prob > 0 count: {np.sum(prob > 0)}")
    print('pval_threshold')
    print(pval_threshold)
    # apply_pval_filterParametertobaseduitepvaluewithofFilteringthefitfor
    sig_prob = prob.copy()
    if apply_pval_filter:
        if r_patcher:
            sig_prob[(pval - pval_threshold) >= 1.49e-8] = 0
        else:
            sig_prob[pval >= pval_threshold] = 0
    print(f"aggregateCell_Cell_Communication After filtering - sig_prob > 0 count: {np.sum(sig_prob > 0)}")
    # diffbecomeThresholdwithoftryshi
    test_thresholds = [0, 0.0001, 0.001, 0.005, 0.006, 0.007, 0.008, 0.009, 0.0095, 0.0099, 0.01, 0.02, 0.05, 0.1]
    for thresh in test_thresholds:
        count = np.sum(sig_prob > thresh)
        print(f"Count with threshold {thresh}: {count}")


    # strength_matrix and count_matrix ofCalculationitemplacetoorbelowtheadd
    print(f"sig_prob shape: {sig_prob.shape}")
    print(f"sig_prob min: {np.min(sig_prob)}, max: {np.max(sig_prob)}, sum: {np.sum(sig_prob)}")
    print(f"sig_prob > 0 count: {np.sum(sig_prob > 0)}")
    print(f"sig_prob > 1e-10 count: {np.sum(sig_prob > 1e-10)}")
    print(f"sig_prob > 0.01 count: {np.sum(sig_prob > 0.01)}")
    
    # eachCell typepairbetweenoftotalintarakutionstrength
    strength_matrix = np.sum(sig_prob, axis=2)
    
    # intarakutioncount
    print("countmatrix >0")
    count_matrix = np.sum(sig_prob > 0, axis=2) #koredaandRtoratiobetecountismanykubecome
    print(count_matrix)
    print("countmatrix >1.49e-8")
    count_matrix = np.sum(sig_prob > 1.49e-8, axis=2)
    print(count_matrix)
    print("countmatrix >0.001") #r_patcherisexistandkiiskorewithcleanup
    count_matrix = np.sum(sig_prob > 0.001, axis=2)
    print(count_matrix)
    if not r_patcher:
        count_matrix = np.sum(sig_prob > 0, axis=2)

    # eachLigand-Receptorpairofcontribution
    lr_contribution = np.sum(np.sum(sig_prob, axis=0), axis=0)
    
    # eachCell typeofsendreceivetotalamount
    outgoing = np.sum(strength_matrix, axis=1)
    incoming = np.sum(strength_matrix, axis=0)
    
    # matrixtheDataFrametoconvert
    strength_df = pd.DataFrame(strength_matrix, index=cell_types, columns=cell_types)
    count_df = pd.DataFrame(count_matrix, index=cell_types, columns=cell_types)
    
    # Networkpointmark
    try:
        network_centrality = calculate_network_centrality({"strength_matrix": strength_df})
    except Exception as e:
        print(f"NetworkincenternatureCalculationError: {str(e)}")
        import traceback
        print(traceback.format_exc())
        network_centrality = pd.DataFrame()

    
    # matrixofStatisticalinfo
    strength_stats = {
        'min': float(np.min(strength_matrix)),
        'max': float(np.max(strength_matrix)),
        'mean': float(np.mean(strength_matrix)),
        'std': float(np.std(strength_matrix)),
        'non_zero': int(np.sum(strength_matrix > 0))
    }
    
    return {
        'strength_matrix': strength_df,
        'count_matrix': count_df,
        'lr_contribution': lr_contribution,
        'outgoing': pd.Series(outgoing, index=cell_types),
        'incoming': pd.Series(incoming, index=cell_types),
        'network_centrality': network_centrality,
        'strength_stats': strength_stats
    }


def debug_interaction_matrix(matrix, title="Interaction Matrix Debug"):
    """
    InteractionmatrixofdebaguinfotheOutput
    
    Parameters
    ----------
    matrix : pd.DataFrame or numpy.ndarray
        Interactionmatrix
    title : str
        title
    """
    print(f"===== {title} =====")
    
    # Check if matrix is empty - handle both DataFrame and ndarray
    if hasattr(matrix, 'empty'):
        # For pandas DataFrame
        if matrix.empty:
            print("emptyofmatrixwithsu")
            return
        values = matrix.values
        print(f"Shape: {matrix.shape}")
        print(f"Min value: {np.min(values)}")
        print(f"Max value: {np.max(values)}")
        print(f"Mean value: {np.mean(values)}")
        print(f"Std value: {np.std(values)}")
        print(f"Non-zero entries: {np.count_nonzero(values)} / {matrix.size}")
        print(f"Sample values:\n{matrix.iloc[:min(3, matrix.shape[0]), :min(3, matrix.shape[1])]}")
    else:
        # For numpy ndarray
        if matrix.size == 0:
            print("emptyofmatrixwithsu")
            return
        print(f"Shape: {matrix.shape}")
        print(f"Min value: {np.min(matrix)}")
        print(f"Max value: {np.max(matrix)}")
        print(f"Mean value: {np.mean(matrix)}")
        print(f"Std value: {np.std(matrix)}")
        print(f"Non-zero entries: {np.count_nonzero(matrix)} / {matrix.size}")
        print(f"Sample values:\n{matrix[:min(3, matrix.shape[0]), :min(3, matrix.shape[1])]}")
    
    print("================\n")


def plot_circle_communication(network_summary, title="Cell-Cell Communication Network", figsize=(10, 10)):
    """
    sa-kiyura-plotwithCellbetweencommunicationthedraw
    
    Parameters
    ----------
    network_summary : dict
        NetworkaggregateResult
    title : str
        plotoftitle
    figsize : tuple
        figureofsize
        
    Returns
    -------
    matplotlib.figure.Figure
        plotofFigureobujiekuto
    """
    try:
        fig, ax = plt.subplots(figsize=figsize)
        
        # strengthmatrixtheget
        if 'strength_matrix' not in network_summary or isinstance(network_summary['strength_matrix'], pd.DataFrame) and network_summary['strength_matrix'].empty:
            ax.text(0.5, 0.5, "Dataisnot exist", ha='center', va='center')
            ax.axis('off')
            plt.title(title)
            return fig
            
        matrix = network_summary['strength_matrix']
        
        # matrixtheNormalization
        if isinstance(matrix, pd.DataFrame):
            if matrix.sum().sum() > 0:
                norm_matrix = matrix / matrix.max().max()
            else:
                norm_matrix = matrix
        else:
            if np.sum(matrix) > 0:
                norm_matrix = matrix / np.max(matrix)
            else:
                norm_matrix = matrix
        
        # color-maptheSettings
        colors = plt.cm.Set3(np.linspace(0, 1, len(matrix)))
        
        # sa-kiyura-plotofdraw（simpleeasyversion）
        # realoccasionofCellChatwithiscirclizeetcofpake-jitheusetethanwashpracticesaretaplotthecreate
        plt.title(title)
        plt.axis('equal')
        plt.axis('off')
        
        # simpleeasyalnasa-kiyura-plotofdami-implement
        theta = np.linspace(0, 2*np.pi, len(matrix)+1)[:-1]
        
        # no-dothedistributeplace
        x = np.cos(theta)
        y = np.sin(theta)
        
        # no-doofdraw
        if isinstance(matrix, pd.DataFrame):
            labels = matrix.index
        else:
            labels = [f"Cluster{i}" for i in range(len(matrix))]
            
        for i, (xi, yi, label) in enumerate(zip(x, y, labels)):
            ax.scatter(xi, yi, s=300, color=colors[i], edgecolor='black', zorder=10)
            ax.text(xi*1.15, yi*1.15, label, ha='center', va='center', fontsize=12, fontweight='bold')
        
        # ejiofdraw
        for i, sender in enumerate(range(len(matrix))):
            for j, receiver in enumerate(range(len(matrix))):
                # DataFrameandNumpyofbothwaytopairrespond
                if isinstance(matrix, pd.DataFrame):
                    sender_label = matrix.index[i]
                    receiver_label = matrix.columns[j]
                    value = matrix.loc[sender_label, receiver_label]
                    if i != j and value > 0:
                        strength = norm_matrix.loc[sender_label, receiver_label]
                else:
                    if i != j and matrix[i, j] > 0:
                        strength = norm_matrix.iloc[i, j]
                    else:
                        continue
                
                # bejiecurvelinewitharcthedraw
                xi, yi = x[i], y[i]
                xj, yj = x[j], y[j]
                # inbetweenpointthefewshizurasu
                xm = (xi + xj) / 2
                ym = (yi + yj) / 2
                # incenterfromseparatesu
                dx = xm
                dy = ym
                d = np.sqrt(dx**2 + dy**2)
                xm += dx / d * 0.3
                ym += dy / d * 0.3
                
                # thicksaisstrengthtoratioExample
                width = 1 + 4 * strength
                ax.plot([xi, xm, xj], [yi, ym, yj], 'gray', linewidth=width, alpha=0.6)
                
                # arrowmarkthedraw
                ax.arrow(xm, ym, (xj-xm)*0.2, (yj-ym)*0.2, width=0.01*strength, 
                        head_width=0.05*width, head_length=0.1*width, 
                        fc='black', ec='black', zorder=5)
        
        return fig
    except Exception as e:
        st.error(f"sa-kiyura-plotcreateError: {str(e)}")
        st.error(traceback.format_exc())
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, f"sa-kiyura-plotcreateError: {str(e)}", ha='center', va='center', wrap=True)
        ax.axis('off')
        plt.title(title)
        return fig

def plot_dot_lr_network(results_df, source_cells, target_cells, top_n=20, pval_threshold=0.05, figsize=(12, 10)):
    """
    dotoplotwithCellbetweenofLRInteractionthedraw
    
    Parameters
    ----------
    results_df : pd.DataFrame
        CellChatofCalculationResult
    source_cells : list
        sendsideofCell typelist
    target_cells : list
        receivesideofCell typelist
    top_n : int
        Displaydouprankpairofcount
    pval_threshold : float
        SignificantandminasuPvalueofThreshold
    figsize : tuple
        figureofsize
        
    Returns
    -------
    matplotlib.figure.Figure
        plotofFigureobujiekuto
    """
    try:
        if results_df.empty:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, "validnaInteractionDataisnot exist", ha='center', va='center')
            ax.axis('off')
            plt.title('Ligand-Receptor Interactions')
            return fig
            
        # SelectshitaCell typetheFiltering
        filtered_df = results_df[
            (results_df['source'].isin(source_cells)) & 
            (results_df['target'].isin(target_cells)) &
            (results_df['pval'] <= pval_threshold)
        ]
        
        if len(filtered_df) == 0:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, f"pointsetshitaConditionwithvalidnaInteractionisnot exist。\nSelectshitaCell typebetweenwithisSignificantnaInteractionischeckoutsaremasenwithshita。\nPvalueThresholdtheupgeruka、sepofCell typetheSelectandplease。", 
                   ha='center', va='center', wrap=True)
            ax.axis('off')
            plt.title('Ligand-Receptor Interactions')
            return fig
        
        # LRpairgoandtoaggregate
        if 'interaction_prob_normalized' not in filtered_df.columns:
            filtered_df['interaction_prob_normalized'] = filtered_df['prob']
            
        lr_summary = filtered_df.groupby(['ligand', 'receptor'])['interaction_prob_normalized'].sum().reset_index()
        lr_summary = lr_summary.sort_values('interaction_prob_normalized', ascending=False).head(top_n)
        
        # Ligand-Receptorpairtheresultmatch
        lr_summary['interaction'] = lr_summary['ligand'] + '-' + lr_summary['receptor']
        
        # eachpairofCell typepairgoandofstrengththeCalculation
        dot_data = []
        
        for lr_pair in lr_summary[['ligand', 'receptor']].itertuples(index=False):
            ligand, receptor = lr_pair
            
            for source in source_cells:
                for target in target_cells:
                    # koofLRpairandkoofCell typepairtopairresponddorowtheextractout
                    subset = filtered_df[
                        (filtered_df['ligand'] == ligand) & 
                        (filtered_df['receptor'] == receptor) &
                        (filtered_df['source'] == source) & 
                        (filtered_df['target'] == target)
                    ]
                    
                    if len(subset) > 0:
                        row = subset.iloc[0]
                        dot_data.append({
                            'interaction': f"{ligand}-{receptor}",
                            'source': source,
                            'target': target,
                            'strength': row['interaction_prob_normalized'],
                            'pvalue': row['pval']
                        })
        
        dot_df = pd.DataFrame(dot_data)
        
        if len(dot_df) == 0:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, "dotoplotforofDataisnot exist", ha='center', va='center')
            ax.axis('off')
            plt.title('Ligand-Receptor Interactions')
            return fig
        
        # sort
        unique_interactions = lr_summary['interaction'].tolist()
        
        # plotcreate
        fig, ax = plt.subplots(figsize=figsize)
        
        # dotoplotofDatalevelprepare
        plot_data = {}
        
        for source in source_cells:
            for target in target_cells:
                cell_pair = f"{source}->{target}"
                plot_data[cell_pair] = {}
                
                source_target_df = dot_df[(dot_df['source'] == source) & (dot_df['target'] == target)]
                
                for interaction in unique_interactions:
                    interaction_df = source_target_df[source_target_df['interaction'] == interaction]
                    
                    if len(interaction_df) > 0:
                        plot_data[cell_pair][interaction] = {
                            'strength': interaction_df['strength'].values[0],
                            'pvalue': interaction_df['pvalue'].values[0]
                        }
                    else:
                        plot_data[cell_pair][interaction] = {
                            'strength': 0,
                            'pvalue': 1
                        }
        
        # eachserupairtopairandplot
        cell_pairs = list(plot_data.keys())
        n_pairs = len(cell_pairs)
        
        if n_pairs == 0:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, "dotoplotforofDataisnot exist", ha='center', va='center')
            ax.axis('off')
            plt.title('Ligand-Receptor Interactions')
            return fig
            
        # verticaltosubplotthecreate
        fig, axes = plt.subplots(n_pairs, 1, figsize=(figsize[0], figsize[1] * n_pairs / 3), sharex=True)
        
        if n_pairs == 1:
            axes = [axes]
        
        for i, (cell_pair, ax) in enumerate(zip(cell_pairs, axes)):
            strengths = []
            sizes = []
            interactions = []
            
            for interaction in unique_interactions:
                data = plot_data[cell_pair][interaction]
                strengths.append(data['strength'])
                # sizeis-log10(pvalue)toratioExample
                sizes.append(max(20, -np.log10(data['pvalue'] + 1e-10) * 20))
                interactions.append(interaction)
            
            # zerowithnotstrengthofmiplot
            non_zero = [i for i, s in enumerate(strengths) if s > 0]
            
            if non_zero:
                # color-mapofSettings
                max_strength = max(np.array(strengths)[non_zero]) if non_zero else 1
                colors = plt.cm.YlOrRd(np.array(strengths)[non_zero] / max_strength)
                
                # scatter plotplot
                scatter = ax.scatter(
                    [interactions[i] for i in non_zero],
                    [0] * len(non_zero),
                    s=[sizes[i] for i in non_zero],
                    c=colors,
                    alpha=0.7,
                    edgecolors='gray'
                )
                
                # Yaxislabeltheserupairnameto
                ax.set_ylabel(cell_pair, fontsize=12)
                
                # YaxiseyeflourishrithenonDisplay
                ax.set_yticks([])
            else:
                ax.text(0.5, 0, f"SignificantnaInteractionnashi", ha='center', va='center')
                ax.set_ylabel(cell_pair, fontsize=12)
                ax.set_yticks([])
            
            # mostafterofsubplotofmiXaxislabelthetimerotate
            if i == len(axes) - 1:
                plt.xticks(rotation=90, ha='right')
                
        plt.suptitle('Ligand-Receptor Interactions Between Cell Types', fontsize=14)
        plt.tight_layout()
        plt.subplots_adjust(top=0.95)
        
        return fig
    except Exception as e:
        st.error(f"dotoplotcreateError: {str(e)}")
        st.error(traceback.format_exc())
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, f"dotoplotcreateError: {str(e)}", ha='center', va='center', wrap=True)
        ax.axis('off')
        plt.title('Ligand-Receptor Interactions')
        return fig

def plot_aggregated_network(network_summary, threshold=0.01, figsize=(12, 10)):
    """
    aggregatesaretaNetworkthegraphshapeformatwithdraw
    
    Parameters
    ----------
    network_summary : dict
        NetworkaggregateResult
    threshold : float
        DisplaydominimumInteractionstrength
    figsize : tuple
        figureofsize
        
    Returns
    -------
    matplotlib.figure.Figure
        plotofFigureobujiekuto
    """
    try:
        # strengthmatrixtheget
        if 'strength_matrix' not in network_summary or isinstance(network_summary['strength_matrix'], pd.DataFrame) and network_summary['strength_matrix'].empty:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, "validnaInteractionDataisnot exist", ha='center', va='center')
            ax.axis('off')
            plt.title('Cell-Cell Communication Network')
            return fig
            
        matrix = network_summary['strength_matrix']
        
        # shikiivalueorbelowtheFiltering
        if isinstance(matrix, pd.DataFrame):
            filtered_matrix = matrix.copy()
            filtered_matrix[filtered_matrix < threshold] = 0
        else:
            filtered_matrix = matrix.copy()
            filtered_matrix[filtered_matrix < threshold] = 0
        
        # graphthecreate
        if isinstance(filtered_matrix, pd.DataFrame):
            G = nx.from_pandas_adjacency(filtered_matrix, create_using=nx.DiGraph())
        else:
            # safesetoffor、matrixthe numpy distributecolumntoconvertandfromgraphgenerate
            adj_mat = np.array(filtered_matrix)
            G = nx.from_numpy_array(adj_mat, create_using=nx.DiGraph())
        # alonestandno-dothedeleteremove
        G.remove_nodes_from(list(nx.isolates(G)))
        
        if len(G) == 0:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, f"Threshold {threshold} orupofInteractionisnot exist。\nThresholdthebelowgetemiteplease。", 
                   ha='center', va='center')
            ax.axis('off')
            plt.title('Cell-Cell Communication Network')
            return fig
        
        # ejiofweighttheseto
        if isinstance(filtered_matrix, pd.DataFrame):
            for u, v, d in G.edges(data=True):
                d['weight'] = filtered_matrix.loc[u, v]
        else:
            for u, v, d in G.edges(data=True):
                d['weight'] = filtered_matrix[u, v]
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # no-doofcolortheSettings
        node_color = plt.cm.Set3(np.linspace(0, 1, len(G)))[:len(G)]
        
        # no-doofrankplacethedecideset
        pos = nx.spring_layout(G, seed=42)
        
        # ejiofthicksatheweighttoratioExamplesaseru
        if isinstance(filtered_matrix, pd.DataFrame):
            max_weight = filtered_matrix.max().max()
        else:
            max_weight = np.max(filtered_matrix)
            
        if max_weight > 0:
            edge_weights = [d['weight'] * 5 / max_weight for u, v, d in G.edges(data=True)]
        else:
            edge_weights = [1.0 for u, v, d in G.edges(data=True)]
        
        # no-dothedraw
        nx.draw_networkx_nodes(G, pos, node_size=700, node_color=node_color, alpha=0.8, ax=ax)
        
        # ejithedraw
        nx.draw_networkx_edges(
            G, pos, width=edge_weights, alpha=0.6, 
            edge_color='gray', arrows=True, 
            arrowstyle='->', arrowsize=20, ax=ax
        )
        
        # labelthedraw
        nx.draw_networkx_labels(G, pos, font_size=12, font_family='sans-serif', ax=ax)

        # titleandadjustarrange
        plt.title('Cell-Cell Communication Network', fontsize=14)
        plt.axis('off')
        plt.tight_layout()
        
        return fig
    except Exception as e:
        st.error(f"NetworkfigurecreateError: {str(e)}")
        st.error(traceback.format_exc())
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, f"NetworkfigurecreateError: {str(e)}", ha='center', va='center', wrap=True)
        ax.axis('off')
        plt.title('Cell-Cell Communication Network')
        return fig


def preprocess_data(
    adata, 
    group_by,
    complex_input,
    gene_use=None,
    min_cells=10, 
    thresh_pct=0,
    thresh_p=0.05,
    resource=None,
    features=None,
    do_de=True,
    min_cells_expr=10
):
    """
    CellChatofDataPreprocessingtherowurelatedcount（Rversionofprocessordertomatchwasetamodifycorrectversion）
    """
    st.write(f"original data: {adata.shape[0]} cells x {adata.shape[1]} genes")
    
    # 1. CellGroupofCellcountchieku koreisRwithisfirsttoisrowwanot
#    cell_counts = adata.obs[group_by].value_counts()
#    valid_groups = cell_counts[cell_counts >= min_cells].index.tolist()
#    if len(valid_groups) < len(cell_counts):
#        st.warning(f"{len(cell_counts) - len(valid_groups)}pieceofCellGroupis{min_cells}Cellnot yetfullofforremoveoutsaremashita")
    
    # validnaGroupofCelldakethekeephold
#    adata_filtered = adata[adata.obs[group_by].isin(valid_groups)].copy()
    
    valid_signaling_genes = [g for g in gene_use if g in adata.var_names]
    st.write(f"LR genes found in data: {len(valid_signaling_genes)}")
    # signalingGeneofmithekeephold（RversionofsubsetDatamutualcurrentofprocess）
    adata_filtered = adata[:, valid_signaling_genes].copy()
    st.write(f"LR gene data: {adata_filtered.shape[0]} cells x {adata_filtered.shape[1]} genes")
    
    if features is None:
        # 3. highExpressionGeneofsameset（RversionofidentifyOverExpressedGenesmutualcurrent）
        overexpressed_genes_result = identify_overexpressed_genes(
            adata_filtered,  # suwithtosignalingGenetoFilteringdonemi
            group_by=group_by,
            thresh_pct=thresh_pct,
            thresh_p=thresh_p,
            only_pos=True,
            do_de=do_de,
            do_fast=True,
            min_cells=min_cells,
            min_cells_expr=min_cells_expr
        )
        
        # highExpressionGenetheextractout
        features_sig = overexpressed_genes_result['features']
    else:
        st.write("Using features...")
        features_sig = features.copy()
    
    # 4. LRpairFiltering（RversionofidentifyOverExpressedInteractionsmutualcurrent）
    # tapuruthereceiveketakeru
    result = identify_overexpressed_interactions(
        features_sig=features_sig, 
        gene_use=valid_signaling_genes,
        resource=resource,
        complex_input=complex_input
    )
    
    # tapurutheexpandopen
    resource_filtered, lr_genes_from_function = result
    
    st.write(f"Filtered LR pairs: {len(resource_filtered)}")
    
    # riso-suisemptywithnotcaseofmiDisplay
    if not resource_filtered.empty:
        st.write(resource_filtered.head(3))
        
    # 5. ExpressiondiffofexistGene＋FilteringsaretaLRpairtorelatedconnectdoGeneofmithekeephold
    lr_related_genes = set(features_sig)  # DEGGenethebasemainandandkeephold
    
    for _, row in resource_filtered.iterrows():  # herethemodifycorrect
        ligand, receptor = row['ligand'], row['receptor']
        # simplepurenaGeneofcaseissoofmamaadd
        if ligand in valid_signaling_genes:
            lr_related_genes.add(ligand)
        if receptor in valid_signaling_genes:
            lr_related_genes.add(receptor)
        
        # complexofcaseissubunitotheadd
        if complex_input is not None:
            if ligand in complex_input.index:
                subunit_cols = [col for col in complex_input.columns if 'subunit' in col]
                subunits = complex_input.loc[ligand, subunit_cols].dropna().astype(str)
                subunits = [s for s in subunits if s != "" and s in valid_signaling_genes]
                lr_related_genes.update(subunits)
            
            if receptor in complex_input.index:
                subunit_cols = [col for col in complex_input.columns if 'subunit' in col]
                subunits = complex_input.loc[receptor, subunit_cols].dropna().astype(str)
                subunits = [s for s in subunits if s != "" and s in valid_signaling_genes]
                lr_related_genes.update(subunits)
    
    # relatedconnectdoGeneofmithekeephold
    lr_related_genes = list(lr_related_genes)
    adata_filtered = adata_filtered[:, lr_related_genes].copy()
    
    st.write(f"Preprocessingafter: {adata_filtered.shape[0]}Cell x {adata_filtered.shape[1]}Gene")
    
    return adata_filtered, resource_filtered


# plotforofyu-teiriteirelatedcount
def scPalette(n):
    """
    Generate colors from a customed color palette
    
    Parameters:
    ----------
    n : int
        Number of colors
        
    Returns:
    -------
    colors : list
        A color palette for plotting
    """
    colorSpace = ['#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC',
                 '#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A','#E3BE00','#FB9A99',
                 '#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B',
                 '#00441B','#DEDC00','#DCF0B9','#8DD3C7','#999999']
    if n <= len(colorSpace):
        colors = colorSpace[:n]
    else:
        # colorRampPalettetomutualcurrentdorelatedcount
        from matplotlib.colors import LinearSegmentedColormap
        cmap = LinearSegmentedColormap.from_list('custom_cmap', colorSpace)
        colors = [cmap(i / n) for i in range(n)]
    return colors

def ggPalette(n):
    """
    Generate ggplot2 colors
    
    Parameters:
    ----------
    n : int
        Number of colors to generate
        
    Returns:
    -------
    colors : list
        A color palette for plotting
    """
    from matplotlib.colors import hsv_to_rgb
    hues = np.linspace(15, 375, n + 1) / 360
    saturations = np.ones(n + 1) * 0.65
    values = np.ones(n + 1) * 0.65
    hsv_colors = np.column_stack((hues, saturations, values))
    rgb_colors = [hsv_to_rgb(hsv) for hsv in hsv_colors]
    return rgb_colors[:n]

# ko-dodaiaguramuVisualizationofforofrelatedcount
def plot_chord_cell(net, color_use=None, group=None, cell_order=None, 
                   sources_use=None, targets_use=None, lab_cex=0.8,
                   small_gap=1, big_gap=10, remove_isolate=False, 
                   transparency=0.4, title_name=None):
    """
    Chord diagram for visualizing cell-cell communication
    
    Parameters
    ----------
    net : numpy.ndarray or pandas.DataFrame
        A weighted matrix or a dataframe with three columns defining the cell-cell communication network
    color_use : dict or list, optional
        Colors for the cell groups
    group : dict, optional
        A named group labels for making multiple-group Chord diagrams
    cell_order : list, optional
        A vector defining the cell type orders
    sources_use : list, optional
        A vector giving the index or the name of source cell groups
    targets_use : list, optional
        A vector giving the index or the name of target cell groups
    lab_cex : float, optional
        Font size for the text
    small_gap : int, optional
        Small gap between sectors
    big_gap : int, optional
        Gap between the different sets of sectors
    remove_isolate : bool, optional
        Whether remove the isolate nodes in the communication network
    transparency : float, optional
        Transparency of link colors
    title_name : str, optional
        Title of the plot
        
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        print("Please install plotly: pip install plotly")
        return None
    
    # InputisDataFrameofcase、matrixtoconvert
    if isinstance(net, pd.DataFrame):
        if all(c in net.columns for c in ["source", "target", "prob"]):
            cell_levels = sorted(list(set(net['source'].unique()) | set(net['target'].unique())))
            net_matrix = np.zeros((len(cell_levels), len(cell_levels)))
            for _, row in net.iterrows():
                src_idx = cell_levels.index(row['source'])
                tgt_idx = cell_levels.index(row['target'])
                net_matrix[src_idx, tgt_idx] = row['prob']
            net = net_matrix
            
    # so-su,ta-getoofFiltering
    if sources_use is not None or targets_use is not None:
        cells_level = np.arange(net.shape[0])
        df_net = pd.DataFrame()
        for i in range(net.shape[0]):
            for j in range(net.shape[1]):
                if net[i, j] > 0:
                    df_net = df_net.append({
                        'source': i,
                        'target': j,
                        'value': net[i, j]
                    }, ignore_index=True)
        
        if sources_use is not None:
            if isinstance(sources_use[0], int):
                df_net = df_net[df_net['source'].isin(sources_use)]
            else:
                indices = [cell_levels.index(s) for s in sources_use if s in cell_levels]
                df_net = df_net[df_net['source'].isin(indices)]
        
        if targets_use is not None:
            if isinstance(targets_use[0], int):
                df_net = df_net[df_net['target'].isin(targets_use)]
            else:
                indices = [cell_levels.index(t) for t in targets_use if t in cell_levels]
                df_net = df_net[df_net['target'].isin(indices)]
        
        # matrixtoconvertshidirectsu
        if not df_net.empty:
            net = np.zeros((len(cell_levels), len(cell_levels)))
            for _, row in df_net.iterrows():
                net[int(row['source']), int(row['target'])] = row['value']
    
    # alonestandno-doofdeleteremove（Option）
    if remove_isolate:
        idx1 = np.where(np.sum(net, axis=0) == 0)[0]
        idx2 = np.where(np.sum(net, axis=1) == 0)[0]
        idx = np.intersect1d(idx1, idx2)
        if len(idx) > 0:
            net = np.delete(np.delete(net, idx, axis=0), idx, axis=1)
    
    # color-Settings
    if color_use is None:
        color_use = scPalette(net.shape[0])
    
    # labelSettings
    if cell_order is None:
        labels = [f"Cell{i+1}" for i in range(net.shape[0])]
    else:
        labels = cell_order
    
    # Plotlytheuseshitako-dodaiaguramuofcreate
    source_list = []
    target_list = []
    value_list = []
    color_list = []
    
    for i in range(net.shape[0]):
        for j in range(net.shape[1]):
            if net.iloc[i, j] > 0:
                source_list.append(i)
                target_list.append(j)
                value_list.append(net.iloc[i, j])
                # colorisso-suofcolortheuse
                if isinstance(color_use, list):
                    color_list.append(color_use[i % len(color_use)])
                else:
                    color_list.append(color_use.get(i, '#377EB8'))
    
    # color-the16advanceshapeformattoconvert（requirednacase）
    for i in range(len(color_list)):
        if isinstance(color_list[i], tuple) and len(color_list[i]) in [3, 4]:
            if max(color_list[i]) <= 1:
                color_list[i] = f'rgba({int(color_list[i][0]*255)},{int(color_list[i][1]*255)},{int(color_list[i][2]*255)},{transparency if len(color_list[i]) == 3 else color_list[i][3]})'
            else:
                color_list[i] = f'rgba({int(color_list[i][0])},{int(color_list[i][1])},{int(color_list[i][2])},{transparency if len(color_list[i]) == 3 else color_list[i][3]/255})'
    
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=20,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels
        ),
        link=dict(
            source=source_list,
            target=target_list,
            value=value_list,
            color=color_list
        )
    )])
    
    fig.update_layout(
        title=title_name,
        font_size=lab_cex * 10
    )
    
    return fig




@st.cache_data
def calculate_network_centrality(net):
    """
    NetworkincenternaturetheCalculationdo - netAnalysis_signalingRole_networkofCalculationmethodtheuse
    
    Parameters
    ----------
    net : dict
        NetworkData（strength_matrixtheincludemu）
        
    Returns
    -------
    pd.DataFrame
        incenternatureCalculationResult
    """
    import networkx as nx
    import pandas as pd
    import numpy as np
    
    # strength_matrixtheget
    if 'strength_matrix' in net:
        if isinstance(net['strength_matrix'], pd.DataFrame):
            adj = net['strength_matrix'].values
            cell_types = net['strength_matrix'].index.tolist()
        else:
            adj = np.array(net['strength_matrix'])
            # cell_typestheget
            if 'net' in net and 'dimnames' in net['net'] and len(net['net']['dimnames']) > 0:
                cell_types = net['net']['dimnames'][0]
            else:
                cell_types = [f"Cell{i+1}" for i in range(adj.shape[0])]
    elif 'net' in net and 'prob' in net['net']:
        # net['net']['prob']theuse
        adj = np.array(net['net']['prob'])
        if adj.ndim == 3:  # 3nextsourcedistributecolumnnara2nextsourcetogatherabout
            adj = np.sum(adj, axis=2)
        
        if 'dimnames' in net['net'] and len(net['net']['dimnames']) > 0:
            cell_types = net['net']['dimnames'][0]
        else:
            cell_types = [f"Cell{i+1}" for i in range(adj.shape[0])]
    else:
        print("validnaInteractionDataisviewtsukarimasen")
        return pd.DataFrame()
    
    # graphthecreate
    G = nx.from_numpy_array(adj, create_using=nx.DiGraph)
    mapping = {i: cell_types[i] for i in range(len(cell_types))}
    G = nx.relabel_nodes(G, mapping)
    
    # incenternaturepointmarkofCalculation
    # 1. Sender & Receiver: sourceofgraphwithCalculation
    sender = dict(G.out_degree(weight='weight'))
    receiver = dict(G.in_degree(weight='weight'))
    
    # 2. Betweenness: weightreversecountconvertshitagraphwithCalculation
    G_betweenness = G.copy()
    for u, v, d in G_betweenness.edges(data=True):
        if 'weight' in d and d['weight'] > 0:
            d['weight'] = 1.0 / d['weight']
    
    try:
        betweenness = nx.betweenness_centrality(G_betweenness, weight='weight')
    except:
        print("mediummediateincenternatureofCalculationtolosefailshimashita")
        betweenness = {node: 0.0 for node in G.nodes()}
    
    # 3. Mediator & Influencer: specialspecialrelatedcountwithCalculation
    try:
        from pages.cellchat_vis import compute_flow_betweenness, compute_information_centrality
        mediator = compute_flow_betweenness(G)
        influencer = compute_information_centrality(G)
    except Exception as e:
        print(f"specialspecialincenternatureCalculationError: {str(e)}")
        # fo-rubaku: marklevelalnaincenternaturepointmarkwithsubstitutefor
        mediator = {node: 0.0 for node in G.nodes()}
        influencer = {node: 0.0 for node in G.nodes()}
    
    # Datafure-muthecreate
    cell_type_list = list(G.nodes())
    centrality_df = pd.DataFrame({
        'cell_type': cell_type_list,
        'degree': [sender[node] + receiver[node] for node in cell_type_list],
        'in_degree': [receiver[node] for node in cell_type_list],
        'out_degree': [sender[node] for node in cell_type_list],
        'betweenness': [betweenness[node] for node in cell_type_list],
        'flowbet': [mediator[node] for node in cell_type_list],  # Mediator
        'info': [influencer[node] for node in cell_type_list]    # Influencer
    })
    
    return centrality_df


# meinprocess
if __name__ == "__main__":
    st.set_page_config(page_title="CellChat-Python", page_icon="💬")
    st.sidebar.title("Options")
    if "cellchat_res" not in st.session_state:
        st.session_state.cellchat_res = None

    if "cellchat_temp_dir" not in st.session_state:
        st.session_state.cellchat_temp_dir = True
        cellchat_temp_dir = "temp/" + str(round(time.time()))
        if not os.path.exists('temp'):
            os.mkdir('temp')
        else:
            clear_old_directories("temp")
            clear_old_files("temp")
        os.mkdir(cellchat_temp_dir)
        st.session_state.cellchat_temp_dir = cellchat_temp_dir
    else:
        cellchat_temp_dir = st.session_state.cellchat_temp_dir

    st.markdown("### CellChat Python implemetation")
 #   st.markdown("###### weight(strength)isRandhoboonecause。interaction countismutualrelateddoissweetme。")

    uploaded_file = st.file_uploader("H5AD for calculation or cellchat result pkl file for visualization", type=['h5ad', 'pkl'], help="H5ADFilewithnewruleAnalysis、PKLFilewithSavedonemiResulttheLoadingmasu")

    if uploaded_file is not None:
        file_type = uploaded_file.name.split('.')[-1].lower()
        
        # PKLFileofcaseisSavedonemiCellChatResultthereadmiintomu
        if file_type == 'pkl':
            try:
                # SavedonemiFilethereadmiintomu
                bytes_data = uploaded_file.read()
                result = pickle.loads(bytes_data)
                st.session_state.cellchat_res = result
                cell_list = []
                if 'net' in result and 'dimnames' in result['net'] and len(result['net']['dimnames']) > 0:
                    cell_list = sorted(result['net']['dimnames'][0])

                if "sorted_order" not in st.session_state:
                    st.session_state.sorted_order = cell_list

                st.markdown("##### Cells:")
                st.write(", ".join(cell_list))
                if 'groupby' in result:
                    # guro-baruchangecountofgroupbythefurthernew
                    groupby = result['groupby']
                st.success("CellChatResultthecorrectnormaltoLoadingmashita！Visualizationthestartwithkimasu")
                
                # UploadsaretaFilenametheSave（Savebuttonfor）
                if 'uploaded_filename' not in st.session_state:
                    st.session_state.uploaded_filename = uploaded_file.name
                
                # herefromisVisualizationtabutodirectconnectadvancemu
                if st.session_state.get('cellchat_res') is not None:
                    # cell_listofextractout - multicountofwaymethodthetrymiru
                    pass
                    
            except Exception as e:
                st.error(f"PKLFileofLoadingtolosefailshimashita: {str(e)}")
                st.error(traceback.format_exc())
        
        # H5ADFileofcaseispassnormalofAnalysisfuro-
        elif file_type == 'h5ad':
            # newshiiFilejudgesetrojiku（sizedependexisttheexcluderemove）
            if should_clear_session_state(uploaded_file):
                clear_cellchat_session_state()
                st.session_state.last_file_name = uploaded_file.name
                st.success("newshiiFileischeckoutsaremashita。kiyashiyuthekuriashimashita。")

            # UploadsaretaFilenametheSave（Savebuttonfor）
            st.session_state.uploaded_filename = uploaded_file.name

            adata = read_h5ad(uploaded_file)

            st.write("Uploaded data:")
            st.write(adata)
            temp_df = pd.DataFrame(
            adata.X[:5,:8].toarray() if scipy.sparse.issparse(adata.X) else adata.X[:5,:8],
            index=adata.obs_names[:5],
            columns=adata.var_names[:8]
            )
            st.write("adata.X (default data layer)")
            st.dataframe(temp_df) 


            meta = adata.obs.columns.to_list()

            for i in ['nFeature_RNA', 'nCount_RNA', 'percent.mt', 'Cell_id']:
                try:
                    meta.remove(i)
                except:
                    pass

            # useforpossiblenareiya-theDisplay
            available_layers = list(adata.layers.keys()) if hasattr(adata, 'layers') else []
            st.write(f"Additional data layers: {', '.join(available_layers) if available_layers else 'NA'}")

            groupby = st.selectbox("Cell type:", meta, 
                                 index=find_first_index_or_zero(meta, ["cell.ident", "seurat_clusters", "cell_type", "louvain"]))
            
            cell_list = sorted(adata.obs[groupby].cat.categories.to_list() if hasattr(adata.obs[groupby], 'cat') else sorted(adata.obs[groupby].unique()))
            st.write(", ".join(cell_list))

            # sorted_orderofinitializetheherewithrowu
            if "sorted_order" not in st.session_state:
                st.session_state.sorted_order = cell_list

            with st.form("Basic settings:"):
                col1, col2 = st.columns(2)
                
                with col1:
                    species = st.radio("Species:", ('mouse', 'human'), 
                                      index=check_species_index(adata.var.index.to_list()[:50]))
                    
                    data_layer = st.selectbox("Using layer (should be normalized data):", 
                                             ['X (default)'] + available_layers,
                                             index=0)
                    
                    selected_types = st.multiselect(
                        "Signaling types to analyze (can choose multiple):",
                        ["Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact", "Non-protein Signaling"],
                        default=["Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"]
                    )

                    st.markdown("---")
                    st.markdown("##### Do not change the following parameters unless necessary.")
                    n_cpus = st.slider("Num CPUs:", 
                                  min_value=1, max_value=os.cpu_count(), step=1, value=os.cpu_count()-2)
                   # r_patcher = st.checkbox("Fine-tune to match R calc results", value=False,
                   #     help="RandofCalculationResultofdifferitheexistdegreedegreeabsorbcollectdo。weighttotsuiteischiekushinakutemonearsimilardo。countischiekushinotandmanykubecomeis、leandirectionisonecausedo。countofthresholdofbaselevelisyayaarbitrarymeaningal。")
                    r_patcher = False

                    n_perms = st.slider("Permutation number:", 
                                min_value=20, max_value=500, step=20, value=100)     


                    min_cells = st.number_input('Min cells in a cell-type:', 
                                               min_value=0, max_value=100, step=1, value=10,
                                               help="SCALAofcellchatwithisFilteringisoff。")

                    expr_prop = st.number_input('Min fraction of expressing cells:', 
                                               min_value=0.0, max_value=0.9, step=0.01, value=0.0)


                with col2:

                    do_de = st.checkbox("Filter pathways by differential expression in at least one cell type.",
                                        value=True, help="dorekaonetsuofClusterwithDEGofmooftheselectsep。offofcase、ExpressionCellcountwithFiltering。")

                    thresh_p = st.number_input('Threshold p value for overexpressed genes:', 
                                               min_value=0.0, max_value=0.9, step=0.01, value=0.05,
                                               help = "WilcoxontoyoruchangeizeGeneSelectofthreshold")

                    # nesutoshitacolormuwithindento
                    subcol1, subcol2 = st.columns([0.1, 0.9])
                    with subcol2:
                        st.write("If DEG filtering is off:")
                        min_cells_expr = st.number_input('Gene filtering for min expressing cells:', 
                                                   min_value=0, max_value=100, step=1, value=10,
                                                   help="DEGtoyoruFilteringisoffofandkitoallCellofinwithExpressionandexistCellcountwithFiltering。kanarilowkushinotand、mushirostrictshikubecome。")
                   
                    st.markdown("---")
                    fun_type = st.radio("Method to calculate average expression per cell group", ['triMean','truncatedMean'], index = 0, help="EstimationsareruLigand-Receptorpairofcountis、CellGroupgoandofmeanGeneExpressionamounttheCalculationdowaymethodtodependexistdo。trimeanisotherofwaymethodthanmofewnotInteractionthegenerate。CellChatisthanstrongiInteractiontheadvancemeasuredonatureableisexcellentreteexistkoandiswakateori、realtestalverifyoffortoInteractionthesqueezeriintomuupwithnonnormaltohavefor。trimeanis25%cutcutmeantonearsimilarandori、koreis1tsuofGroupwithExpressionandexistCellofratiomatchis25%not yetfullofcase、meanGeneExpressionamountiszerotobecomekoandthemeaningtastedo。truncatedMeanwithis5%and10%cutcutmeanequaltheSettingswithkiru")

                    st.markdown("*To reduce stringency for pathway selection, use truncatedMean.*")
                    trim = st.number_input("Trimming for truncated mean:", min_value=0.00, max_value=0.25, value=0.10)

                    st.markdown("---")                    
                    apply_pval_filter = st.checkbox("Filter interaction pathways by P value", value=True)

                    trim_threshold = st.number_input('Filtering P threshold:', 
                                               min_value=0.01, max_value=0.20, step=0.01, value=0.05)
                    st.markdown('---')
                    population_size = st.checkbox("population.size", value=False, help="CellcountismanyiClustersamepersonofcommunicationhodostrongku（probabilityishighmeto）evaluatepricedocase。")
      
                
                submitted_basic = st.form_submit_button("Apply settings")

            if submitted_basic:
                st.session_state.cellchat_res = None

            st.markdown("### Click [Apply settings] before run!")

            if st.button('Run calc') or (st.session_state.get('cellchat_res') is not None):
                if st.session_state.get('cellchat_res') is None:
                    # CellChatDBtheget
                    try:
                        start_cpu = time.time()
                        with st.spinner('Getting CellChatDB from R...'):
                            cellchatdb = get_cellchatdb_from_r(species=species)
                            # debaguOutputtheadd
                        cellchatdb = debug_cellchatdb(cellchatdb)
                        st.success(f"{species} CellChatDB is succcesfully obtained")
                    except Exception as e:
                        st.error(f"CellChatDBofgetintoErrorisoccurred: {str(e)}")
                        st.stop()
                    interaction = cellchatdb['interaction']
                    interaction_filtered = interaction[interaction['annotation'].isin(selected_types)]
                    cellchatdb['interaction'] = interaction_filtered
                    
                    st.write(f"LR pair number: {len(cellchatdb['interaction'])}")
                    st.write(cellchatdb['interaction'].head(3))
                    
                    # Datareiya-ofSettings
                    use_layer = None if data_layer == 'X (default)' else data_layer

                    #herewithdbtheexpandopenandoku
                    gene_use, resource, complex_input, cofactor_input, gene_info = extractGene(cellchatdb)
                    # herewith、resource isFilteringdonemiof interaction（tsumari、selected_types tothiscurrentdopair）ofmitheincludemuiszuwithsu。
                    st.write(f"Resource LR pair numbers:{len(resource)}")
                    print(f'gene_use {len(gene_use)}')
                    print(f'complex_input {len(complex_input)}')
                    print(f'cofactor_input {len(cofactor_input)}')
                    print(f'gene_info {len(gene_info)}')


                    # Databe-suofverify
                    if resource.empty:
                        raise ValueError("DBofInteractioninfo(interaction)isemptywithsu。validnaCellChatDBkaConfirmandplease。")
                        
                    # requirednacolumnofexistattheConfirm
                    for required_col in ['ligand', 'receptor']:
                        if required_col not in resource.columns:
                            raise ValueError(f"riso-suDatafure-mutois '{required_col}' columnisrequiredwithsu")

                    st.write(f"LR genes in db: {len(gene_use)}")

                    with st.spinner('CellChat calculation...'):
                        result = cellchat_analysis(
                            adata,
                            groupby=groupby,
                          #  db=cellchatdb,
                            gene_use=gene_use,
                            complex_input=complex_input,
                            cofactor_input=cofactor_input,
                            resource=resource,
                            use_layer=use_layer,
                            min_cells=min_cells,
                            min_cells_expr=min_cells_expr,
                            expr_prop=expr_prop,
                         #   log_scale=log_transformed,
                            pseudocount=1.0,
                            k=0.5,
                            trim_threshold=trim_threshold,  
                            nboot=n_perms,
                            seed=1,
                            n_jobs=n_cpus,
                            fun_type=fun_type,
                            trim=trim,
                            apply_pval_filter=apply_pval_filter,
                            r_patcher=r_patcher,
                            population_size=population_size,
                            do_de=do_de
                        )
                        
                        st.session_state.cellchat_res = result
                    
                    if 'error' in result:
                        st.error(f"AnalysisintoErrorisoccurred: {result['error']}")
                        st.code(result['traceback'])
                    else:
                        st.success('CellChatAnalysisisCompleteshimashita！')
                        time_cpu = time.time() - start_cpu
                        st.write(f"Computing time: {round(time_cpu)}")
                        print(f"Computing time: {round(time_cpu)}")
                
            # orbelowisResultofDisplaypartdiv
            
        result = st.session_state.get('cellchat_res')

        if result is not None:
            
            if 'error' in result:
                st.error(f"AnalysisintoErrorisoccurred: {result['error']}")
                if 'traceback' in result:
                    st.code(result['traceback'])
            else:
                # ResultDisplay
                st.subheader("1. CellChat results")
                
                st.write("Interactionte-buru (uprank10item):")
                if 'results' in result and len(result['results']) > 0:
                    st.write(result['results'].sort_values('prob', ascending=False).head(10))
                    

                    # mazu、ResultofCSV（Interactionte-buru）thetextcharcolumntoconvert
                    csv_str = result['results'].to_csv(index=False, sep = '\t')

                    # nextto、strength_matrixandcount_matrixtheTSVtoconvert
                    # DataFrameofcase
                    if isinstance(result['network']['strength_matrix'], pd.DataFrame):
                        strength_tsv = result['network']['strength_matrix'].to_csv(sep='\t', index=True)
                    else:
                        # numpydistributecolumnofcase、DataFrameizeandfromTSVize（herewithistemptodimnamesisexistandtempset）
                        strength_df = pd.DataFrame(result['network']['strength_matrix'], 
                                                   index=result['network']['dimnames'][0], 
                                                   columns=result['network']['dimnames'][1])
                        strength_tsv = strength_df.to_csv(sep='\t', index=True)

                    if isinstance(result['network']['count_matrix'], pd.DataFrame):
                        count_tsv = result['network']['count_matrix'].to_csv(sep='\t', index=True)
                    else:
                        count_df = pd.DataFrame(result['network']['count_matrix'], 
                                                index=result['network']['dimnames'][0], 
                                                columns=result['network']['dimnames'][1])
                        count_tsv = count_df.to_csv(sep='\t', index=True)

                    # memoriupwithZIPFilethecreate
                    zip_buffer = io.BytesIO()
                    with zipfile.ZipFile(zip_buffer, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
                        zf.writestr("cellchat_results.tsv", csv_str)
                        zf.writestr("strength_matrix.tsv", strength_tsv)
                        zf.writestr("count_matrix.tsv", count_tsv)

                    # bafuaofpointatheaheadheadtoreturnsu
                    zip_buffer.seek(0)

                    # st.download_buttonwithZIPFiletheprovideprovide
                    st.download_button(
                        label="ResulttheZIPandandDownload",
                        data=zip_buffer,
                        file_name="cellchat_results.zip",
                        mime="application/zip"
                    )


                else:
                    st.warning("validnaInteractionResultisnot exist。Parametertheadjustarrangeandagaintryrowandplease。")
                    
                # debaguinfotheDisplay
                with st.expander("debaguinfo"):
                    if 'network' in result and 'strength_matrix' in result['network']:
                        debug_interaction_matrix(result['network']['strength_matrix'], "Interactionstrengthmatrix")
                    if 'network' in result and 'count_matrix' in result['network']:
                        debug_interaction_matrix(result['network']['count_matrix'], "Interactioncountmatrix")
                
                # Visualization
                st.subheader("2. Visualization of communication network")
                st.markdown("###### Visualization options are displayed at the bottom of the left side panel")
                plt.style.use('default')
                with st.sidebar:
                    colormap_options = [
                        "YlOrRd", "OrRd", "YlOrBr", "Oranges", "Reds", "RdPu", "Purples", 
                        "PuRd", "Blues", "Greens", "YlGn", "YlGnBu", "GnBu", "Greys", "binary",
                        "viridis", "plasma", "inferno", "magma", "cividis", 
                        "viridis_r", "plasma_r", "inferno_r", "magma_r", "cividis_r"
                    ]
                    color_heatmap = st.selectbox(
                        "Heatmap color:",
                        colormap_options,
                        index=2
                    )
                    cmap_name = st.sidebar.selectbox(
                        "Select colormap for other plots:",
                        ["tab10","Set1", "Set2", "Set3", "tab20", "Paired", "Dark2",
                        "tab20b", "tab20c","Pastel1",
                        "Pastel2",  "Accent", "viridis", "plasma", "inferno", "magma"],
                        index=0
                    )

                    # sorted_order ofeternalcontinueizeimplement（sourceofcellchat.pyofFunctiontherecoversource）
                    if "sorted_order" not in st.session_state:
                        st.session_state.sorted_order = cell_list
                        sorted_order = None
                    else:
                        sorted_order = st.session_state.get('sorted_order')

                    # color-mapofeternalcontinueizeimplement
                    # color-mapischangefurthersaretacasealsoisfirstmeteofcase
                    if 'cell_color_map' not in st.session_state or st.session_state.get('current_cmap', '') != cmap_name:
                        if sorted_order:
                            st.session_state.cell_color_map = create_cell_color_mapping(sorted_order, cmap_name)
                        else:
                            st.session_state.cell_color_map = create_cell_color_mapping(cell_list, cmap_name)
                        st.session_state.current_cmap = cmap_name

                    sort_cell = st.checkbox("Change cell order?")
                    # session_statetoenterreru
                    if sort_cell:
                        with st.form("Sorter"):
                            #sorted_order = sort_items(sorted(adata.obs[groupby].unique()))
                            sorted_order = sort_items(st.session_state.get('sorted_order', cell_list).copy())
                            submitted_sort = st.form_submit_button("Done sorting")
                        st.session_state.sorted_order = sorted_order
                        # sorted_orderchangefurthertimeiscurrentofcolor-paretowithcolormapthefurthernew
                        current_cmap = st.session_state.get('current_cmap', cmap_name)
                        st.session_state.cell_color_map = create_cell_color_mapping(sorted_order, current_cmap)
                     #   st.write(f"Color map changed{st.session_state.cell_color_map }")
                    else:
                        sorted_order = None
                    
                    # Add source and target filtering options
                    st.markdown("#### Cell Type Filtering")
                    filter_cells = st.checkbox("Filter source/target cells?")
                    if filter_cells:
                        sources_use = st.multiselect(
                            "Select source cell types:",
                            options=cell_list,
                            default=cell_list,
                            help="Select which cell types to use as sources (signal senders)"
                        )
                        targets_use = st.multiselect(
                            "Select target cell types:",
                            options=cell_list,
                            default=cell_list,
                            help="Select which cell types to use as targets (signal receivers)"
                        )
                    else:
                        sources_use = None
                        targets_use = None
                
                tabs = st.tabs([
                    "Heatmap", 
                    "Chord", 
                    "LR contrib", 
                    "Dot", 
                    "Centrality", 
                    "Role sctr",
                    "Sig contrib",
                    "Circle",   # New tab for circle plots
                    "Sig role",   # New tab for signaling role analysis 
                    "Expression"     # New tab for gene expression analysis
                ])

                # Keep your existing tab implementations
                with tabs[0]:
                    st.markdown("#### Cell interaction heatmap")
                    
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        
                        heatmap_type = st.radio(
                            "Heatmap type:",
                            ["Interaction strength", "Interaction number"],
                            horizontal=True
                        )
                        heatmap_annot = st.checkbox("Show values?", value=True)
                        show_color_bar = st.checkbox("Show cell-type color bars?", value=False)

                    with col2:
                        heatmap_title = st.slider("Title font size:", min_value=0, max_value=30, value=20)
                        heatmap_font = st.slider("Font size:", min_value=0, max_value=30, value=14)
                        heatmap_x = st.slider("Fig width:", min_value=1.0, max_value=20.0, value=10.0, step=0.2)
                        heatmap_y = st.slider("Fig height:", min_value=1.0, max_value=20.0, value=8.0, step=0.2)


                    if st.button("Generate heatmap"):
                        
                        
                        if 'network' in result:
                            if heatmap_type == "Interaction strength":
                                matrix = result['network']['strength_matrix']
                                title = "Cell-Cell Interaction Strength"
                            else:
                                matrix = result['network']['count_matrix']
                                title = "Number of Significant Interactions"

                            # Use the enhanced heatmap visualization
                            measure = "weight" if heatmap_type == "Interaction strength" else "count"
                            fig_heatmap = netVisual_heatmap(
                                result, 
                                measure=measure, 
                                color_heatmap=color_heatmap, 
                                font_size=heatmap_font, 
                                font_size_title=heatmap_title,
                                sorted_order=sorted_order,
                                annot=heatmap_annot,
                                color_use=st.session_state.get('cell_color_map', None),  # colormapingutheadd
                                show_color_bar=show_color_bar,
                                figx=heatmap_x,
                                figy=heatmap_y,
                                sources_use=sources_use,
                                targets_use=targets_use
                            ) 
                            st.pyplot(fig_heatmap)
                            
                            # PDF saving and download unchanged
                            pdf_path = f"{cellchat_temp_dir}/heatmap.pdf"
                            fig_heatmap.savefig(pdf_path, bbox_inches='tight')
                            
                            with open(pdf_path, "rb") as pdf_file:
                                pdf_bytes = pdf_file.read()
                                st.download_button(
                                    label="Download heatmap pdf",
                                    data=pdf_bytes,
                                    file_name=f'cellchat_{heatmap_type}_heatmap.pdf',
                                    mime='application/octet-stream'
                                )
                        else:
                            st.warning("Networkinfoisuseforwithkimasen")

                with tabs[1]:  # ko-dodaiaguramu tab
                    st.markdown("#### Chord diagram")
                    
                    if 'netP' in result and 'pathways' in result['netP']:

                        col1, col2 = st.columns(2)

                        with col1:
                            # Selectlimb: "Aggregate" alsois pieceseppathwayofmulticountSelect
                            option_type = st.radio(
                                "Display type:",
                                ["Aggregate", "Specific pathways"],
                                horizontal=True
                            )
                            if option_type == "Aggregate":
                                selected_pathway = "Aggregate"
                                # Interaction measure theSelectdorajiobutton
                                diagram_measure = st.radio(
                                    "Select value:",
                                    ["Interaction weight", "Interaction number"],
                                    horizontal=True
                                )
                                measure = "weight" if diagram_measure == "Interaction weight" else "count"
                            else:
                                # multicountofpathwaytheSelectpossibletodo
                                pathways = list(result['netP']['pathways'])
                                selected_pathways = st.multiselect(
                                    "Select pathways (multiple selection allowed):",
                                    sorted(pathways),
                                    default=[pathways[0]] if pathways else []
                                )
                                selected_pathway = selected_pathways  # listandandtransfersu
                                measure = "weight"  # specificpathwaywithisnormaltoweight
                        
                        with col2:
                            # DisplaykasutamaizuOption
                            sector_space = st.number_input("Space between sectors (degree):", min_value=0, max_value=20, value=5)
                            alpha_edge = st.slider("Transparency of edge:", min_value=0.1, max_value=1.0, value=0.4, step=0.1)
                            edge_border_width = st.slider("Edge border width", min_value=0.0, max_value=1.0, value=0.4, step=0.1)
                           # show_edge_border = st.checkbox("Show edge border", value=True)
                            chord_x = st.slider("Fig X size:", min_value=4.0, max_value=20.0, value=10.0, step=0.5)
                            chord_y = st.slider("Fig Y size:", min_value=4.0, max_value=20.0, value=10.0, step=0.5)
                        
                        if st.button("Generate chord diagram"):
                            if option_type == "Specific pathways" and not selected_pathways:
                                st.warning("Please select at least one pathway.")
                            else:
                                with st.spinner("Making chord diagram..."):
                                    from pages.cellchat_vis import create_chord_diagram_pycirclize
                                    try:

                                        fig_chord = create_chord_diagram_pycirclize(
                                            net=result,
                                            pathway_name=selected_pathway,
                                            color_use=st.session_state.get('cell_color_map', None),
                                            sorted_order=sorted_order,
                                            measure=measure,
                                            figsize=(chord_x, chord_y),
                                            space=sector_space,               # sekuta-betweenofemptywhite 5degree
                                            alpha_edge=alpha_edge,        # rinkuoftransparentcleardegree 0.6
                                        #    show_edge_border=show_edge_border, # rinkuofblackiwheeloutlinethenonDisplay
                                            edge_border_width=edge_border_width,
                                            sources_use=sources_use,
                                            targets_use=targets_use
                                        )
                                        st.pyplot(fig_chord)
                                        
                                        # PDF andandSave,Download
                                        import io
                                        pdf_buffer = io.BytesIO()
                                        fig_chord.savefig(pdf_buffer, format="pdf", bbox_inches="tight")
                                        pdf_data = pdf_buffer.getvalue()
                                        
                                        if option_type == "Aggregate":
                                            filename = 'cellchat_chord_aggregate.pdf'
                                        else:
                                            # multicountpathwayofcaseispathwaycounttheDisplay
                                            if len(selected_pathway) <= 2:
                                                pathway_str = '_'.join(selected_pathway)
                                            else:
                                                pathway_str = f'{len(selected_pathway)}_pathways'
                                            filename = f'cellchat_chord_{pathway_str}.pdf'
                                        
                                        st.download_button(
                                            label="Download chord diagram",
                                            data=pdf_data,
                                            file_name=filename,
                                            mime='application/octet-stream'
                                        )
                                    except Exception as e:
                                        st.error(f"ko-dodaiaguramugenerateError: {str(e)}")
                                        import traceback
                                        st.error(traceback.format_exc())
                    else:
                        st.warning("Resulttopathwayinfoisnot exist。")
                
                
                with tabs[2]:
                    st.markdown("#### Ligand-Receptorpairofcontribution")
                    
                    st.info("""
                    koofgraphiseachLigand-ReceptorpairisallbodyofsignalingNetworktodoofdegreedegreecontributegiveandexistkatheshowandimasu。
                    contributionishighipairhodo、Cell-cell communicationtooiteweightneednaroleratiotheresulttaandimasu。
                    possiblenacaseis、eachpairisbelongdosignalingpathwaymoDisplaysaremasu。
                    """)
                    
                    if 'network' in result and 'lr_contribution' in result['network']:
                        lr_contribution = result['network']['lr_contribution']
                        
                        if isinstance(lr_contribution, np.ndarray) and len(lr_contribution) > 0:
                            # DataFrameofinfofromDataoflengththeget
                            lr_contrib_len = np.sum(lr_contribution > 0)
                            
                            # DisplaydotopuLRpaircountofSelect（maximumvaluetheadjustarrange）
                            top_n = st.slider("DisplaydotopuLRpaircount:", 
                                             min_value=1, 
                                             max_value=min(30, int(lr_contrib_len)) if lr_contrib_len > 0 else 10, 
                                             value=min(10, int(lr_contrib_len)) if lr_contrib_len > 0 else 10)
                            
                            # color-sukeymuSelect
                            color_scheme = st.selectbox(
                                "color-sukeymu:", 
                                ["viridis", "plasma", "inferno", "magma", "cividis", "Blues", "Greens", "Oranges", "Reds"],
                                index=0
                            )

                            if st.button("Generate LR contribution graph"):
                            
                                # plotofhighsathemovealtoadjustarrange（paircounttoratioExample）
                                plot_height = 6 + 0.3 * top_n  # paircounttorespondjitehighsatheadjustarrange
                                from pages.cellchat_vis import  plot_lr_contribution

                                # Networksamari-toLRnameinfotheadd
                                network_data = result['network'].copy()
                                
                                # interaction_namesofget
                                interaction_names = result['net']['dimnames'][2]
                                
                                # LRpairinfoofwordwritethecreate
                                lr_pairs = {}
                                for i, interaction_name in enumerate(interaction_names):
                                    # resultsDatafure-mufromLRinfothechecksearch
                                    matching_rows = result['results'][result['results']['interaction_name'] == interaction_name]
                                    if not matching_rows.empty:
                                        row = matching_rows.iloc[0]
                                        lr_pairs[interaction_name] = f"{row['ligand']}-{row['receptor']}"
                                
                                # network_datatoLRinfotheadd
                                network_data['lr_names'] = lr_pairs
                                network_data['interaction_names'] = interaction_names
                                
                                # LRcontributionofgraphthegenerate
                                fig_lr = plot_lr_contribution(
                                    network_data, 
                                    top_n=top_n, 
                                    figsize=(12, plot_height)
                                )

                                
                                # figureofsutairuthefurthernew
                                for ax in fig_lr.get_axes():
                                    # color-mapthechangefurther
                                    if hasattr(ax, 'patches') and len(ax.patches) > 0:
                                        cmap = plt.cm.get_cmap(color_scheme)
                                        colors = cmap(np.linspace(0.2, 0.9, len(ax.patches)))
                                        for i, patch in enumerate(ax.patches):
                                            patch.set_color(colors[i])
                                    
                                    # backscenecolorofSettings
                                 #   ax.set_facecolor('#f8f9fa')  # thinigure-
                                    
                                    # guridolineofSettings
                                 #   ax.grid(axis='x', linestyle='--', alpha=0.3)
                                    
                                    # framelineofSettings
                                    for spine in ax.spines.values():
                                        spine.set_color('#cccccc')
                                
                                st.pyplot(fig_lr)
                                
                                # PDFSaveandDownload
                                pdf_path = f"{cellchat_temp_dir}/lr_contribution.pdf"
                                fig_lr.savefig(pdf_path, bbox_inches='tight')
                                
                                with open(pdf_path, "rb") as pdf_file:
                                    pdf_bytes = pdf_file.read()
                                    st.download_button(
                                        label="LRcontributionPDFtheDownload",
                                        data=pdf_bytes,
                                        file_name='cellchat_lr_contribution.pdf',
                                        mime='application/octet-stream'
                                    )
                                
                                # contributionuprankofLRpairinfothete-buruwithmoDisplay（Option）
                                if st.checkbox("contributionuprankofLRpairdetailstheDisplay"):
                                    try:
                                        # contributionDataofget,arrangeshape
                                        lr_data = []
                                        for i, val in enumerate(lr_contribution):
                                            if val > 0:
                                                # LRpairinfoofgettryrow
                                                lr_info = {"index": i, "contribution": val}
                                                
                                                # ResultDatafure-mufromofinfoget
                                                if 'results' in result and not result['results'].empty:
                                                    # intarakutionnameofgettryrow
                                                    if 'net' in result and 'dimnames' in result['net'] and len(result['net']['dimnames']) > 2:
                                                        interaction_names = result['net']['dimnames'][2]
                                                        if i < len(interaction_names):
                                                            matching_rows = result['results'][
                                                                result['results']['interaction_name'] == interaction_names[i]
                                                            ]
                                                            if not matching_rows.empty:
                                                                row = matching_rows.iloc[0]
                                                                lr_info.update({
                                                                    "interaction_name": row.get('interaction_name', f"LR_{i}"),
                                                                    "ligand": row.get('ligand', ''),
                                                                    "receptor": row.get('receptor', ''),
                                                                    "pathway": row.get('pathway_name', '')
                                                                })
                                                
                                                # infoisgetwithkinakatacaseoffo-rubaku
                                                if 'interaction_name' not in lr_info:
                                                    lr_info.update({
                                                        "interaction_name": f"LR_{i}",
                                                        "ligand": "notclear",
                                                        "receptor": "notclear",
                                                        "pathway": "notclear"
                                                    })
                                                
                                                lr_data.append(lr_info)
                                        
                                        # contributionwithdescendingsort
                                        lr_df = pd.DataFrame(lr_data).sort_values('contribution', ascending=False).head(top_n)
                                        
                                        # columnnamethedaymainwordtochangefurther
                                        lr_df.columns = [
                                            'index' if col == 'index' else
                                            'contribution' if col == 'contribution' else
                                            'intarakutionname' if col == 'interaction_name' else
                                            'Ligand' if col == 'ligand' else
                                            'Receptor' if col == 'receptor' else
                                            'pathway' if col == 'pathway' else col
                                            for col in lr_df.columns
                                        ]
                                        
                                        # DisplaycolumnoforderSettings
                                        display_cols = ['intarakutionname', 'Ligand', 'Receptor', 'pathway', 'contribution']
                                        display_cols = [col for col in display_cols if col in lr_df.columns]
                                        
                                        # te-buruDisplay
                                        st.dataframe(lr_df[display_cols])
                                    
                                    except Exception as e:
                                        st.error(f"LRpairdetailsinfoofgetError: {str(e)}")
                        else:
                            st.warning("LRcontributionDataisnot exist")
                    else:
                        st.warning("LRcontributioninfoisuseforwithkimasen")
                
                with tabs[3]:
                    st.markdown("#### LRInteractiondotoplot")
                    
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        source_cells = st.multiselect(
                            'LigandExpressionCellofSelect:',
                            cell_list,
                            default=[cell_list[0]] if len(cell_list) > 0 else []
                        )
                    
                    with col2:
                        target_cells = st.multiselect(
                            'ReceptorExpressionCellofSelect:',
                            cell_list,
                            default=[cell_list[0]] if len(cell_list) > 0 else []
                        )
                    
                    if 'results' in result:

                    
                        pval_col = 'pval' if 'pval' in result['results'].columns else 'pvalue'
                        if pval_col in result['results'].columns:
                            sig_lr_count = len(result['results'][result['results'][pval_col] <= 0.05])
                        else:
                            sig_lr_count = len(result['results'])
                            
                        top_n_dot = st.slider("DisplaydotopuLRpaircount (dotoplot):", 
                                             min_value=1, 
                                             max_value=min(30, sig_lr_count) if sig_lr_count > 0 else 15, 
                                             value=min(15, sig_lr_count) if sig_lr_count > 0 else 15)
                        pval_threshold = st.slider("PvalueThreshold:", min_value=0.001, max_value=0.1, value=0.05, step=0.001)

                        if st.button("Generate dotplot"):
                            
                            if source_cells and target_cells:
                                fig_dot = plot_dot_lr_network(
                                    result['results'],
                                    source_cells,
                                    target_cells,
                                    top_n=top_n_dot,
                                    pval_threshold=pval_threshold
                                )
                                
                                st.pyplot(fig_dot)
                                
                                # PDFSaveandDownload
                                pdf_path = f"{cellchat_temp_dir}/dot_plot.pdf"
                                fig_dot.savefig(pdf_path, bbox_inches='tight')
                                
                                with open(pdf_path, "rb") as pdf_file:
                                    pdf_bytes = pdf_file.read()
                                    st.download_button(
                                        label="dotoplotPDFtheDownload",
                                        data=pdf_bytes,
                                        file_name='cellchat_dot_plot.pdf',
                                        mime='application/octet-stream'
                                    )
                            else:
                                st.warning("sendsideandreceivesideofCell typetheSelectandplease。")
                    else:
                        st.warning("Resultinfoisuseforwithkimasen")
                
                with tabs[4]:
                    st.markdown("#### Networkincenternature")
                    
                    if 'network' in result and 'network_centrality' in result['network']:
                        from pages.cellchat_vis import plot_network_centrality

                        st.subheader("In/Out activity")
                        
                        if 'outgoing' in result['network'] and 'incoming' in result['network']:
                            outgoing = result['network']['outgoing']
                            incoming = result['network']['incoming']
                            
                            # sendreceiveDatathema-ji
                            if isinstance(outgoing, pd.Series) and not outgoing.empty:

                                if st.button("Generate activity plot"):
                                    io_df = pd.DataFrame({
                                        'cell_type': outgoing.index,
                                        'outgoing': outgoing.values,
                                        'incoming': incoming.values
                                    })
                                    
                                    st.write(io_df)
                                    
                                    fig, ax = plt.subplots(figsize=(10, 6))
                                    
                                    x = np.arange(len(io_df))
                                    width = 0.35
                                    
                                    ax.bar(x - width/2, io_df['outgoing'], width, label='Outgoing')
                                    ax.bar(x + width/2, io_df['incoming'], width, label='Incoming')
                                    
                                    ax.set_xticks(x)
                                    ax.set_xticklabels(io_df['cell_type'], rotation=45, ha='right')
                                    ax.legend()
                                    
                                    ax.set_title('Outgoing vs Incoming Communication Activity')
                                    ax.set_ylabel('Interaction Strength')
                                    
                                    plt.tight_layout()
                                    
                                    st.pyplot(fig)


                                    # PDFSaveandDownload
                                    pdf_path = f"{cellchat_temp_dir}/centrality.pdf"
                                    fig.savefig(pdf_path, bbox_inches='tight')
                                    
                                    with open(pdf_path, "rb") as pdf_file:
                                        pdf_bytes = pdf_file.read()
                                        st.download_button(
                                            label="Download in/out PDF",
                                            data=pdf_bytes,
                                            file_name='cellchat_in_out_activity.pdf',
                                            mime='application/octet-stream'
                                        )

                            else:
                                st.warning("sendreceiveactivenatureDataisemptywithsu")
                        else:
                            st.warning("sendreceiveactivenatureinfoisuseforwithkimasen")
                    else:
                        st.warning("Networkincenternatureinfoisuseforwithkimasen")
                
                with tabs[5]:
                    st.markdown("#### signalroleratioanalysis")
                    
                    st.info("""
                    signalroleratioanalysisis、communicationNetworktookerueachCell typeof
                    sendpersonalsoisreceivepersonandandofcontributededicatedegreeandallbodyalnashadoweffectpowertheshowshimasu。
                    """)
                    
                    # UIthe2colormutodivratioanddistributeplace
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        # pathwaySelect（Aggregate alsoisspecificpathway）
                        if 'netP' in result and 'pathways' in result['netP']:
                            pathways = result['netP']['pathways']
                            selected_path = st.selectbox(
                                "roleratioanalysisofpathwaytheSelect:",
                                options=["Aggregate"] + sorted(list(pathways)),
                                index=0
                            )
                        else:
                            st.warning("Resulttopathwayinfoisnot exist。")
                            selected_path = "Aggregate"
                        
                        # XaxisandYaxisofincenternaturepointmarktheSelect
                        measure_options = ["outdeg", "indeg", "outdeg_unweighted", "indeg_unweighted"]
                        if 'netP' in result and 'centr' in result['netP'] and len(result['netP']['centr']) > 0:
                            first_key = 0 if isinstance(result['netP']['centr'], list) else list(result['netP']['centr'].keys())[0]
                            if isinstance(result['netP']['centr'][first_key], dict):
                                measure_options = list(result['netP']['centr'][first_key].keys())
                        
                        x_measure = st.selectbox("X axis:", options=measure_options, 
                                                index=measure_options.index("outdeg") if "outdeg" in measure_options else 0)
                        y_measure = st.selectbox("Y axis:", options=measure_options, 
                                                index=measure_options.index("indeg") if "indeg" in measure_options else 0)
                    
                    with col2:
                        # DisplaykasutamaizuOption
                        do_label = st.checkbox("Cell typelabeltheDisplay", value=True)
                        dot_alpha = st.slider("pointoftransparentcleardegree:", min_value=0.1, max_value=1.0, value=0.7, step=0.1)
                        use_count = st.checkbox("Dot size reflects link number", value=True)
                        cell_groups = None
                    # detailsSettingsforofekusupanda-
                    with st.expander("Options"):
                        label_size = st.slider("labelsize:", min_value=4, max_value=16, value=8)
                        dot_size_min = st.slider("minimumpointsize:", min_value=1, max_value=10, value=2)
                        dot_size_max = st.slider("maximumpointsize:", min_value=dot_size_min+1, max_value=20, value=6)
                        
                        # Weight Min/MaxSettings
                        weight_minmax_input = st.text_input(
                            "Weight MinMax (minimumvalue,maximumvalue):",
                            value=""
                        )
                        weight_MinMax = None
                        if weight_minmax_input.strip():
                            try:
                                min_val, max_val = map(float, weight_minmax_input.split(","))
                                weight_MinMax = (min_val, max_val)
                            except:
                                st.warning("Weight MinMaxofshapeformatiscorrectshikunot exist。'minimumvalue,maximumvalue'ofshapeformatwithInputandplease。")
                        
                        # titleandaxislabel
                        title = st.text_input("title:", value="Signaling Role Analysis")
                      #  title=" " #natokaenterreruandplotweightnateviewenakubecomeis、NonedaandgdrawfacewithumakuDisplaysarenot
                        xlabel = st.text_input("Xaxislabel:", value="Outgoing interaction strength" if x_measure == "outdeg" else x_measure)
                        ylabel = st.text_input("Yaxislabel:", value="Incoming interaction strength" if y_measure == "indeg" else y_measure)
                        role_x = st.slider("Fig width:", min_value=1, max_value=20, value=6)
                        role_y = st.slider("Fig height:", min_value=1, max_value=20, value=4)
                    
                    # analysisRunbutton
                    if st.button("Generate signaling role scatter plot"):
                        with st.spinner("Calculating..."):
                            try:
                                from pages.cellchat_vis import netAnalysis_signalingRole_scatter
                                
                                # pathwayParameterSettings
                                pathway_param = None if selected_path == "Aggregate" else [selected_path]
                                
                                # scatter plotthegenerate
                                fig_role = netAnalysis_signalingRole_scatter(
                                    net=result,
                                    signaling=pathway_param,
                                    color_use=st.session_state.get('cell_color_map', None),
                                    slot_name="netP",
                                    group=cell_groups,
                                    weight_MinMax=weight_MinMax,
                                    dot_size=(dot_size_min, dot_size_max),
                                    label_size=label_size,
                                    dot_alpha=dot_alpha,
                                    x_measure=x_measure,
                                    y_measure=y_measure,
                                    xlabel=xlabel,
                                    ylabel=ylabel,
                                    title=title,
                                    do_label=do_label,
                                    width=role_x,
                                    height=role_y,
                                    use_count=use_count,
                                    sources_use=sources_use,
                                    targets_use=targets_use,
                                    sorted_order=sorted_order
                                )
                                
                                # figuretheDisplay
                                st.pyplot(fig_role)
                                
                                # PDFofSaveandDownloadprovideprovide
                                pdf_path = f"{cellchat_temp_dir}/signaling_role_{selected_path}.pdf"
                                fig_role.savefig(pdf_path, bbox_inches='tight')
                                
                                with open(pdf_path, "rb") as pdf_file:
                                    pdf_bytes = pdf_file.read()
                                    st.download_button(
                                        label="signalroleratioscatter plotthePDFwithDownload",
                                        data=pdf_bytes,
                                        file_name=f'cellchat_signaling_role_{selected_path}.pdf',
                                        mime='application/octet-stream'
                                    )
                            except Exception as e:
                                st.error(f"signalroleratioscatter plotError: {str(e)}")
                                st.error(traceback.format_exc())  # tore-subakutheDisplayanddebaguthecontenteasyto
                            

                
                with tabs[6]:
                    st.markdown("#### signalcontributededicatedegreeanalysis")
                    
                    st.info("""
                    koofhi-tomapwithis、specificofCellGroupofsendalsoisreceivesignaltomostmocontributegivedosignalPathwaytheshowshimasu。
                    color-ba-is、CellGroupbetweenofsignalingPathwayofmutualpairalnastrengththetableshimasu（note: valueisrowwaydirectionwithsuke-ringusareteimasu）。
                    uppartofcolor-ba-plotis、hi-tomaptoDisplaysareruallofsignalingPathwaythematchcountshitaCellGroupoftotalsignalingstrengththeshowshimasu。
                    rightsideofgure-ba-plotis、hi-tomaptoDisplaysareruallofCellGroupthematchcountshitasignalingPathwayoftotalsignalingstrengththeshowshimasu。
                    """)
                    from pages.cellchat_vis import netAnalysis_signalingRole_heatmap
                    
                    # pata-nSelect（send/receive）
                    pattern = st.radio(
                        "signalpata-n:",
                        ["outgoing", "incoming"],
                        horizontal=True,
                        help="outgoing: Cellfromofsendsignal, incoming: Cellheofreceivesignal"
                    )
                    
                    # signalPathwayofSelect
                    if 'netP' in result and 'pathways' in result['netP']:
                        pathway_options = result['netP']['pathways']
                        
                        selection_mode = st.radio(
                            "signalPathwaySelect:",
                            ["all", "Select"],
                            horizontal=True
                        )
                        
                        if selection_mode == "Select":
                            selected_pathways = st.multiselect(
                                "analysisdoPathway:",
                                options=sorted(pathway_options),
                                default=sorted(pathway_options)[:min(5, len(pathway_options))],
                                help="analysisdospecificofsignalPathwaytheSelectandplease。"
                            )
                            signaling_param = selected_pathways if selected_pathways else None
                        else:
                            signaling_param = None  # allofPathwaytheanalysis
                            st.info(f"allofPathwaytheanalysisshimasu ({len(pathway_options)}pieceofPathway)")
                    else:
                        st.warning("Resulttopathwayinfoisnot exist。")
                        signaling_param = None
                    
                    # VisualizationSettings
                    with st.expander("VisualizationSettings"):
                        col1, col2 = st.columns(2)
                        
                        with col1:

                            
                            display_numbers = st.checkbox(
                                "Show values in heatmap",
                                value=False,
                                help="hi-tomapofserutocountvaluetheDisplayshimasu"
                            )
                        
                        with col2:
                            cluster_rows = st.checkbox(
                                "rowthekurasutaringu",
                                value=True,
                                help="floorlayeralkurasutaringutheuseandrowthearrangereasonshimasu"
                            )
                            
                            cluster_cols = st.checkbox(
                                "columnthekurasutaringu",
                                value=False,
                                help="floorlayeralkurasutaringutheuseandcolumnthearrangereasonshimasu"
                            )
                        
                        thresh = st.slider(
                            "SignificantnatureThreshold:",
                            min_value=0.01,
                            max_value=0.5,
                            value=0.05,
                            step=0.01,
                            help="SignificantandminasuInteractionofpvalueThreshold"
                        )
                        
                        title = st.text_input(
                            "kasutamutitle:",
                            value="",
                            help="figureofkasutamutitle（emptywhiteofcaseisselfmovegenerate）"
                        )
                        
                        width = st.slider("figureofwidth:", min_value=6, max_value=15, value=10, step=1)
                        height = st.slider("figureofhighsa:", min_value=4, max_value=12, value=8, step=1)
                    
                    # hi-tomapgeneratebutton
                    if st.button("signalcontributededicatedegreehi-tomapthegenerate"):
                        with st.spinner("hi-tomapthegeneratein..."):
                            try:
                                # hi-tomapgenerate
                                # Generate the heatmap
                                fig_heatmap = netAnalysis_signalingRole_heatmap(
                                    net=result,
                                    signaling=signaling_param,
                                    pattern=pattern,
                                    title=title if title else None,
                                    color_heatmap=color_heatmap,
                                    thresh=thresh,
                                    width=width,
                                    height=height,
                                    font_size=10,
                                    cluster_rows=cluster_rows,
                                    cluster_cols=cluster_cols,
                                    display_numbers=display_numbers,
                                    sorted_order=sorted_order,
                                    cmap_name=cmap_name,
                                    sources_use=sources_use,
                                    targets_use=targets_use
                                )
                                
                                # figuretheDisplay
                                st.pyplot(fig_heatmap)
                                
                                # PDFandandSave,Downloadprovideprovide
                                pathway_str = "all_pathways" if signaling_param is None else "_".join(signaling_param[:3]) + (f"_plus{len(signaling_param)-3}" if len(signaling_param) > 3 else "")
                                pdf_path = f"{cellchat_temp_dir}/signaling_role_heatmap_{pattern}_{pathway_str}.pdf"
                                fig_heatmap.savefig(pdf_path, bbox_inches='tight')
                                
                                with open(pdf_path, "rb") as pdf_file:
                                    pdf_bytes = pdf_file.read()
                                    st.download_button(
                                        label="signalcontributededicatedegreehi-tomapthePDFwithDownload",
                                        data=pdf_bytes,
                                        file_name=f'cellchat_signaling_role_heatmap_{pattern}_{pathway_str}.pdf',
                                        mime='application/octet-stream'
                                    )
                            
                            except Exception as e:
                                st.error(f"hi-tomapofgenerateintoErrorisoccurred: {str(e)}")
                                st.error(traceback.format_exc())

                with tabs[7]:
                    st.markdown("#### Network circle plot")
                    
                    # signalSelect（Aggregate andeachpathway）
                    if 'netP' in result and 'pathways' in result['netP']:


                        col1, col2 = st.columns(2)

                        with col1:
                            # Selectlimb: "Aggregate" alsois pieceseppathwayofmulticountSelect
                            option_type = st.radio(
                                "Pathway type:",
                                ["Aggregate", "Specific pathways"],
                                horizontal=True
                            )
                            
                            if option_type == "Aggregate":
                                selected_pathway = "Aggregate"
                                # Interaction measure theSelectdorajiobutton
                                aggregate_measure = st.radio(
                                    "Select measure:",
                                    ["Interaction Weight", "Interaction Number"],
                                    horizontal=True
                                )
                                measure_key = "weight" if aggregate_measure == "Interaction Weight" else "count"
                            else:
                                # multicountofpathwaytheSelectpossibletodo
                                pathways = list(result['netP']['pathways'])
                                selected_pathways = st.multiselect(
                                    "Select pathways (selected pathways are summed):",
                                    sorted(pathways),
                                    default=[pathways[0]] if pathways else []
                                )
                                selected_pathway = selected_pathways  # listandandtransfersu
                                measure_key = "weight"  # specificpathwaywithisnormaltoweight
                            circle_type = st.radio(
                                    "Circle plot type:",
                                    ["Total network", "Individual cell-type"],
                                    horizontal=True
                                )

                        with col2:
                            # DisplaykasutamaizuOption
                            edge_width_max = st.slider("Max edge width:", min_value=1, max_value=20, value=8, step=1)
                            circle_alpha_edge = st.slider("Edge transparency:", min_value=0.1, max_value=1.0, value=0.6, step=0.1)
                            vertex_size_max= st.slider("Max node size:", min_value=1, max_value=15, value=6, step=1)
                            show_vertex = st.checkbox("Node size reflects the number of cells", value=False)

                    else:
                        st.warning("pathwayinfoisuseforwithkimasen")
                        selected_pathway = None
                        measure_key = "weight"


                    
                    # requiredtorespondjite vertex_weight moget（Example：result['adata'].obs fromeachCellofcount）
                    vertex_weight = None
                    if hasattr(result['adata'].obs, 'value_counts') and show_vertex:
                        try:
                            cell_counts = result['adata'].obs[groupby].value_counts()
                          #  st.write(cell_counts)
                            vertex_weight = [cell_counts.get(ct, 1) for ct in result['net']['dimnames'][0]]
                          #  st.write(vertex_weight)
                        except Exception as e:
                            st.warning(f"vertex_weight ofgettolosefailshimashita: {str(e)}")

                    if st.button("Generate circle plot"):
                        if option_type == "Specific pathways" and (not selected_pathways or len(selected_pathways) == 0):
                            st.warning("Please select at least one pathway.")
                        else:
                            try:
                                if option_type == "Aggregate":
                                    # aggregateNetworkofcase
                                    title = f"Cell-Cell {aggregate_measure}"
                                    with st.spinner("aggregatesa-kuruplotthegeneratein..."):
                                        if circle_type == "Total network":
                                            fig_circle = netVisual_circle(
                                                net=result,
                                                title_name=title,
                                                pathway_name=selected_pathway,
                                                measure=measure_key,
                                                sorted_order=sorted_order,
                                                edge_width_max=edge_width_max,
                                                alpha_edge=circle_alpha_edge,
                                                vertex_weight=vertex_weight,
                                                vertex_size_max=vertex_size_max,
                                                color_use=st.session_state.get('cell_color_map', None),  # colormapingutheadd
                                                sources_use=sources_use,
                                                targets_use=targets_use
                                            )
                                        else:
                                            fig_circle = netVisual_circle_individual(
                                                net=result,
                                                measure=measure_key,
                                                title_name=title,
                                                cmap_name=cmap_name,
                                                sorted_order=sorted_order,
                                                edge_width_max=edge_width_max,
                                                alpha_edge=circle_alpha_edge,
                                                vertex_weight=vertex_weight,
                                                vertex_size_max=vertex_size_max,
                                                arrow=True,
                                                color_use=st.session_state.get('cell_color_map', None) 
                                            )
                                else:
                                    # piecesepsignalalsoismulticountsignalofcase
                                    if len(selected_pathways) == 1:
                                        title = f"Cell-Cell Interaction: {selected_pathways[0]}"
                                    else:
                                        if len(selected_pathways) <= 3:
                                            title = f"Combined Pathways: {', '.join(selected_pathways)}"
                                        else:
                                            title = f"Combined {len(selected_pathways)} Pathways"
                                    
                                    with st.spinner(f"sa-kuruplotthegeneratein... ({title})"):
                                        if circle_type == "Total network":
                                            fig_circle = netVisual_circle(
                                                net=result,
                                                title_name=title,
                                                pathway_name=selected_pathway,
                                                measure="weight",  # piecesepofcaseispassnormal weight theuse
                                                cmap_name=cmap_name,
                                                edge_width_max=edge_width_max,
                                                sorted_order=sorted_order,
                                                vertex_weight=vertex_weight,
                                                vertex_size_max=vertex_size_max,
                                                alpha_edge=circle_alpha_edge,
                                                arrow=True,
                                                color_use=st.session_state.get('cell_color_map', None),
                                                sources_use=sources_use,
                                                targets_use=targets_use
                                            )
                                        else:
                                            fig_circle = netVisual_circle_individual(
                                                net=result,
                                                pathway_name=selected_pathway,
                                                measure="weight",  # piecesepofcaseispassnormal weight theuse
                                                title_name=title,
                                                cmap_name=cmap_name,
                                                edge_width_max=edge_width_max,
                                                sorted_order=sorted_order,
                                                vertex_weight=vertex_weight,
                                                vertex_size_max=vertex_size_max,
                                                alpha_edge=circle_alpha_edge,
                                                arrow=True,
                                                color_use=st.session_state.get('cell_color_map', None) 
                                            )
                                st.pyplot(fig_circle)
                                
                                # PDF Save,Download
                                import io
                                pdf_buffer = io.BytesIO()
                                fig_circle.savefig(pdf_buffer, format="pdf", bbox_inches='tight')
                                pdf_data = pdf_buffer.getvalue()
                                
                                # FilenameofSettings
                                if option_type == "Aggregate":
                                    filename = 'cellchat_circle_aggregate.pdf'
                                else:
                                    # multicountpathwayofcaseispathwaycounttheDisplay
                                    if len(selected_pathway) <= 2:
                                        pathway_str = '_'.join(selected_pathway)
                                    else:
                                        pathway_str = f'{len(selected_pathway)}_pathways'
                                    filename = f'cellchat_circle_{pathway_str}_{circle_type.replace(" ", "_")}.pdf'
                                
                                st.download_button(
                                    label="Download circle plot",
                                    data=pdf_data,
                                    file_name=filename,
                                    mime='application/octet-stream'
                                )
                            except Exception as e:
                                st.error(f"sa-kuruplotgenerateError: {str(e)}")
                                import traceback
                                st.error(traceback.format_exc())


                with tabs[8]:
                    st.markdown("#### signalroleratioanalysis")
                    
                    st.info("""
                    signalroleratioanalysisis、communicationNetworktookerueachCell typeof
                    sendperson(Sender)、receiveperson(Receiver)、relationmediateperson(Mediator)、shadoweffectperson(Influencer)andandofroleratiotheshowshimasu。
                    hi-tomapofcoloriseachCell typeofmutualpairalnaweightneeddegreethetableandimasu（eachrowwithmaximumvaluethe1toNormalization）。
                    """)
                    # 2colormutodivkeru
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        if 'netP' in result and 'pathways' in result['netP']:
                            pathways = result['netP']['pathways']
                            selected_path = st.selectbox(
                                "Select pathways for heatmap:",
                                options=["Aggregate"] + sorted(list(pathways)),
                                index=0
                            )
                                                    # debagufortodetailstheOutput
                            print("Pathways:", pathways)
                            print("Pathways type:", type(pathways))
                            print("Pathways length:", len(pathways))
                        show_value_role = st.checkbox("Show values in plot?", value=False)
                        hide_color_bar = not st.checkbox("Show cell color bar?", value=False)
                    with col2:
                        font_size = st.slider("Font size:", min_value=8, max_value=14, value=10)
                    
                        # hi-tomapofsizeSettings
                        width = st.slider("Fig width:", min_value=6, max_value=15, value=10)
                        height = st.slider("Fig height:", min_value=1.0, max_value=8.0, value=3.0, step=0.5)

                        pathways = result['netP']['pathways']


                    if st.button("Generate plot"):
                        with st.spinner("Calculating..."):
                            try:
                                from pages.cellchat_vis import netAnalysis_signalingRole_network
                                
                                # roleratioanalysisofRun
                                fig_role = netAnalysis_signalingRole_network(
                                    result,  # obujiekutoallbodythetransfersu
                                    signaling=selected_path,
                                    font_size=font_size,
                                    width=width,
                                    height=height,
                                    color_heatmap=color_heatmap,
                                    sorted_order=sorted_order,  # Cellorderofpointset
                                    show_value=show_value_role,
                                    cmap_name=cmap_name,
                                    hide_color_bar=hide_color_bar,
                                    color_use=st.session_state.get('cell_color_map', None)
                                )
                                
                                # figuretheDisplay
                                st.pyplot(fig_role)
                                
                                # PDFandandSave
                                pdf_path = f"{cellchat_temp_dir}/signaling_role_{selected_path}.pdf"
                                fig_role.savefig(pdf_path, bbox_inches='tight')
                                
                                # PDFDownloadbuttontheprovideprovide
                                with open(pdf_path, "rb") as pdf_file:
                                    pdf_bytes = pdf_file.read()
                                    st.download_button(
                                        label="Download PDF",
                                        data=pdf_bytes,
                                        file_name=f'cellchat_signaling_role_{selected_path}.pdf',
                                        mime='application/octet-stream'
                                    )
                                    
                            except Exception as e:
                                st.error(f"signalroleratioanalysisError: {str(e)}")
                                st.error(traceback.format_exc())

                # Add the Gene Expression Analysis tab
                with tabs[9]:
                    st.markdown("#### GeneExpressionanalysis")
                                    
                    st.info("""
                    specificofsignalingpathwaytorelatedwaruGeneofExpressiontheCell typebetweenwithVisualizationshimasu。
                    """)

                    if 'cellchatdb' not in locals() or cellchatdb is None:
                        gene_list = result['adata'].var_names.tolist()
                        species_index = check_species_index(gene_list)
                        species = 'human' if species_index == 1 else 'mouse'
                        try:
                            with st.spinner(f'{species} CellChatDBthegetin...'):
                                cellchatdb = get_cellchatdb_from_r(species=species)
                                st.success(f"{species} CellChatDBthecorrectnormaltogetshimashita")
                        except Exception as e:
                            st.error(f"CellChatDBofgetintoErrorisoccurred: {str(e)}")
                            cellchatdb = None


                    # pathwaySelect（Aggregate alsoisspecificpathway）
                    if ('netP' in result) and ('pathways' in result['netP']):
                        pathways = result['netP']['pathways']
                        selected_path_expr = st.selectbox(
                            "Pathways to show gene expression:",
                            options=sorted(list(pathways)),
                            index=0
                        )
                      
                        if st.button("GeneExpressiontheanalysis"):
                            from pages.cellchat_vis import plotGeneExpression
                            with st.spinner("GeneExpressiontheanalysisin..."):
                                # adataisexistatshinotcaseisWarningtheDisplayandprocesstheincut
                                if 'adata' not in result or result['adata'] is None:
                                    st.warning("AnnDataisreadmiintomareteimasen。")
                                else:
                                    try:
                                        fig_expr = plotGeneExpression(
                                            result,
                                            signaling=selected_path_expr,
                                          #  adata=adata,
                                            cellchatdb=cellchatdb,  # CellChatDBthetransfersu
                                            group_by=groupby,
                                            cmap_name=cmap_name
                                        )
                                        
                                        st.pyplot(fig_expr)
                                        
                                        # Save to PDF
                                        pdf_path = f"{cellchat_temp_dir}/gene_expression_{selected_path_expr}.pdf"
                                        fig_expr.savefig(pdf_path, bbox_inches='tight')
                                        
                                        with open(pdf_path, "rb") as pdf_file:
                                            pdf_bytes = pdf_file.read()
                                            st.download_button(
                                                label="GeneExpressionanalysisPDFtheDownload",
                                                data=pdf_bytes,
                                                file_name=f'cellchat_gene_expression_{selected_path_expr}.pdf',
                                                mime='application/octet-stream'
                                            )
                                    except Exception as e:
                                        st.error(f"GeneExpressionanalysisError: {str(e)}")
                    else:
                        st.warning("Resulttopathwayinfoisnot exist。")
                
                # Visualizationsekutionofmostafter、AnalysistheagainRunbuttonofdirectbeforetoorbelowtheadd:
                st.markdown("---")

                # AnalysisResultofSavebutton
                if st.session_state.get('cellchat_res') is not None:
                    if st.button('Prepare CellChat analysis result to download'):
                        try:
                            # Filenameofcreate
                            if uploaded_file:
                                base_filename = sanitize_filename(uploaded_file.name, 20)
                            else:
                                base_filename = "cellchat_result"
                            
                            # signaltaiputextcharcolumnofcreate
                            signal_str = "_".join([s[:3] for s in selected_types])
                            
                            # SaveforDatathelevelprepare
                            import pickle
                            output = pickle.dumps(st.session_state.get('cellchat_res'))
                            
                            # DownloadbuttontheDisplay
                            st.download_button(
                                label="Download the result file",
                                data=output,
                                file_name=f"cellchat_{base_filename}_{signal_str}.pkl",
                                mime="application/octet-stream"
                            )
                            
                            st.success("AnalysisResulttheSaveshimashita。Downloaddotoisbuttonthekurikuandplease。")
                        except Exception as e:
                            st.error(f"SaveintoErrorisoccurred: {str(e)}")


                if st.button('Rerun analysis'):
                    st.session_state.cellchat_res = None
                    st.rerun()

    else:
        st.info("h5adFilekapklFiletheUploadandstartandplease。h5adisAnalysisfor、pkisDisplayfor。")
    