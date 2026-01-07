import streamlit as st
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.packages import importr
import numpy as np
import pandas as pd
import pickle
import os
import tempfile
import traceback

# apurioftaitoruandheda-
st.title("CellChat QS Converter")
st.markdown("CellChat QSFilethePythonusepickleFiletochangechangedotsu-ru")

# Rringboundofinitialize (seshiyonstarttimetoonedegreedakeRun)
@st.cache_resource
def initialize_r_environment():
    pandas2ri.activate()
    numpy2ri.activate()
    try:
        base = importr('base')
        qs = importr('qs')
        # CellChatmoinpo-to
        ro.r('library(CellChat)')
        return True
    except Exception as e:
        st.error(f"Rringboundofinitializetofailfailshimashita: {str(e)}")
        st.error("requirednaRpake-ji(CellChat, qs)isinsuto-rusareteexistkaConfirmshitekudasai。")
        return False

# Rringboundtheinitialize
r_initialized = initialize_r_environment()

# QSFiletheprocprocdorelnum
def process_qs_file(uploaded_file):
    # onetimeFiletoSave
    with tempfile.NamedTemporaryFile(suffix='.qs', delete=False) as tmp_file:
        tmp_file.write(uploaded_file.getvalue())
        tmp_path = tmp_file.name
    
    try:
        # specspecialtextchartheesuke-pu
        tmp_path_escaped = tmp_path.replace('\\', '\\\\')
        
        # mostmosimplepurenaapuro-chi: mazuRwithchangenumtheSettingsshi、soretheuseusedo
        ro.r(f'temp_file_path <- "{tmp_path_escaped}"')
        
        r_code = """
        library(CellChat)
        library(qs)

        # QSFilethereadmiintomu
        tryCatch({
            cellchat_obj <- qs::qread(temp_file_path)
            
            # requirednaDatatheextractout
            extracted_data <- list()

            # 1. metaData
            if (!is.null(cellchat_obj@meta)) {
                extracted_data$meta <- as.data.frame(cellchat_obj@meta)
            }

            # 2. Networkinfoinfo
            if (!is.null(cellchat_obj@net$prob)) {
                # certainraterowcol
                extracted_data$net_prob <- cellchat_obj@net$prob
                
                # nextsourcenamethepieceseptoSave
                if (!is.null(dimnames(cellchat_obj@net$prob)[[1]])) {
                    extracted_data$dim_rows <- dimnames(cellchat_obj@net$prob)[[1]]
                }
                if (!is.null(dimnames(cellchat_obj@net$prob)[[2]])) {
                    extracted_data$dim_cols <- dimnames(cellchat_obj@net$prob)[[2]]
                }
                if (!is.null(dimnames(cellchat_obj@net$prob)[[3]])) {
                    extracted_data$dim_paths <- dimnames(cellchat_obj@net$prob)[[3]]
                }
            }

            if (!is.null(cellchat_obj@net$pval)) {
                extracted_data$net_pval <- cellchat_obj@net$pval
            }

            if (!is.null(cellchat_obj@net$weight)) {
                extracted_data$net_weight <- as.matrix(cellchat_obj@net$weight)
            }
            
            if (!is.null(cellchat_obj@net$count)) {
                extracted_data$net_count <- as.matrix(cellchat_obj@net$count)
            }

            # 3. pasuueiinfoinfo
            if (!is.null(cellchat_obj@netP$pathways)) {
                extracted_data$pathways <- cellchat_obj@netP$pathways
            }

            if (!is.null(cellchat_obj@netP$prob)) {
                extracted_data$netP_prob <- cellchat_obj@netP$prob
            }

            # 4. LRinfoinfo
            if (!is.null(cellchat_obj@LR$LRsig)) {
                extracted_data$lr_sig <- as.data.frame(cellchat_obj@LR$LRsig)
            }

            # returnrival
            extracted_data
        }, error = function(e) {
            list(error = paste("QSFileofLoadingError:", e$message))
        })
        """
        
        # Rko-dotheRun
        result_r = ro.r(r_code)
        
        # Errorchieku
        if isinstance(result_r, ro.vectors.ListVector) and 'error' in result_r.names:
            error_msg = result_r.rx2('error')
            st.error(f"QSFileofprocproctofailfailshimashita: {error_msg}")
            return None

        # Resultofchangechange
        cellchat_data = {}
        for key in result_r.names:
            try:
                # ki-oftypetoyotechangechangewaymethodthechangeeru
                if key in ["net_prob", "net_pval", "netP_prob"]:
                    # RdistcolfromNumPydistcolheofchangechange
                    cellchat_data[key] = np.array(result_r.rx2(key))
                else:
                    # soofotherofobujiekutoispandas2riwithchangechange
                    r_value = result_r.rx2(key)
                    try:
                        cellchat_data[key] = pandas2ri.rpy2py(r_value)
                    except Exception as e2:
                        # changechangetofailfailshitaplacematch、basemainalnawaymethodwithDatathegetget
                        if isinstance(r_value, ro.vectors.StrVector):
                            cellchat_data[key] = list(r_value)
                        elif isinstance(r_value, ro.vectors.FloatVector):
                            cellchat_data[key] = list(r_value)
                        elif isinstance(r_value, ro.vectors.IntVector):
                            cellchat_data[key] = list(r_value)
                        else:
                            cellchat_data[key] = str(r_value)
                st.success(f"{key}thecorrectnormaltoextractoutshimashita")
            except Exception as e:
                st.warning(f"{key}ofchangechangetofailfailshimashita: {e}")
                st.warning(traceback.format_exc())

        # nextsourcenameofarrangeproc
        if all(key in cellchat_data for key in ['dim_rows', 'dim_cols', 'dim_paths']):
            cellchat_data['net_prob_dims'] = {
                'rows': cellchat_data['dim_rows'],
                'cols': cellchat_data['dim_cols'],
                'pathways': cellchat_data['dim_paths']
            }
            # piecesepofnextsourceki-thedeleteremove（Option）
            for key in ['dim_rows', 'dim_cols', 'dim_paths']:
                if key in cellchat_data:
                    del cellchat_data[key]

        # Resultthearrangeshape
        result = format_cellchat_for_visualization(cellchat_data)
        
        # onetimeFileofdeleteremove
        os.unlink(tmp_path)
        
        return result
        
    except Exception as e:
        st.error(f"ProcessingtoErrorisoccurgenshimashita: {e}")
        st.error(traceback.format_exc())
        # onetimeFileofdeleteremovethetrymiru
        try:
            os.unlink(tmp_path)
        except:
            pass
        return None

def format_cellchat_for_visualization(cellchat_data):
    """extractoutshitaDatathecellchat_vis.pyofbijiyuaraize-shiyonusetoarrangeshapedo"""
    # cellchat.pyofOutputshapeformattomatchwaseteDatathestructmakeize
    result = {}
    
    # 1. Networkinfoinfo
    # 3nextsourcedistcolofprobandpvalisexistkaConfirm
    has_3d_prob = 'net_prob' in cellchat_data and isinstance(cellchat_data['net_prob'], np.ndarray) and cellchat_data['net_prob'].ndim == 3
    has_3d_pval = 'net_pval' in cellchat_data and isinstance(cellchat_data['net_pval'], np.ndarray) and cellchat_data['net_pval'].ndim == 3
    
    result['net'] = {}
    
    # probandpvalis3nextsourceofplacematchofmiSettings
    if has_3d_prob:
        result['net']["prob"] = cellchat_data['net_prob']
    else:
        # Dataisnakerebaemptyof3nextsourcedistcolthemakebecome
        if 'net_prob_dims' in cellchat_data and all(key in cellchat_data['net_prob_dims'] for key in ['rows', 'cols']):
            rows = len(cellchat_data['net_prob_dims']['rows'])
            cols = len(cellchat_data['net_prob_dims']['cols'])
            paths = len(cellchat_data['net_prob_dims'].get('pathways', []))
            # fewnakuandmo1tsuofpasuueithekeepproof
            paths = max(1, paths)
            result['net']["prob"] = np.zeros((rows, cols, paths))
        else:
            # nextsourceinfoinfomonotplacematchisminimumsaizuofdistcol
            result['net']["prob"] = np.zeros((1, 1, 1))
    
    if has_3d_pval:
        result['net']["pval"] = cellchat_data['net_pval']
    else:
        # probandsamesaizuofdistcolthemakebecome
        result['net']["pval"] = np.ones_like(result['net']["prob"])
    
    # nextsourcenameofinfoinfotheaddadd
    if 'net_prob_dims' in cellchat_data:
        dims = cellchat_data['net_prob_dims']
        cell_types_rows = dims.get('rows', [])
        cell_types_cols = dims.get('cols', [])
        pathways = dims.get('pathways', [])
        
        # textcharcolrisutotochangechange
        cell_types_rows = [str(ct) for ct in cell_types_rows]
        cell_types_cols = [str(ct) for ct in cell_types_cols]
        pathways = [str(p) for p in pathways]
        
        # minimumlimitofsaizuthekeepproof
        if not cell_types_rows:
            cell_types_rows = [f"Cell_{i}" for i in range(result['net']["prob"].shape[0])]
        if not cell_types_cols:
            cell_types_cols = [f"Cell_{i}" for i in range(result['net']["prob"].shape[1])]
        if not pathways:
            pathways = [f"Pathway_{i}" for i in range(result['net']["prob"].shape[2])]
        
        result['net']['dimnames'] = [
            cell_types_rows,
            cell_types_cols,
            pathways
        ]
    else:
        # nextsourcenameisnotplacematchisdistcolofsaizutomatchwasetegenbecome
        shape = result['net']["prob"].shape
        result['net']['dimnames'] = [
            [f"Cell_{i}" for i in range(shape[0])],
            [f"Cell_{i}" for i in range(shape[1])],
            [f"Pathway_{i}" for i in range(shape[2])]
        ]
    
    # gathercalcshitaweightrowcolthemakebecome
    cell_types = result['net']['dimnames'][0]
    num_cells = len(cell_types)
    
    # allpasuueitowatarumatchcalc（gathercalcNetwork）
    prob_sum = np.sum(result['net']["prob"], axis=2)
    
    # weightrowcolandcountrowcoltheDataFrametochangechange
    result['net']["weight"] = prob_sum
    result['net']["count"] = (prob_sum > 0).astype(int)
    
    # gathercalcNetworkofmidcenternaturetheCalculation（NetworkXtheuseuse）
    import networkx as nx
    
    # initialize
    net_centr = {
        'outdeg': np.zeros(num_cells),
        'indeg': np.zeros(num_cells),
        'outdeg_unweighted': np.zeros(num_cells),
        'indeg_unweighted': np.zeros(num_cells),
        'betweenness': np.zeros(num_cells),
        'page_rank': np.zeros(num_cells),
        'hub': np.zeros(num_cells),
        'authority': np.zeros(num_cells),
        'eigen': np.zeros(num_cells),
        'flowbet': np.zeros(num_cells),
        'info': np.zeros(num_cells)
    }
    
    # Networkisvalid（ejiisexist）placematchofmiCalculation
    if np.any(prob_sum > 0):
        # gurafumakebecome
        G = nx.DiGraph()
        for i in range(num_cells):
            G.add_node(i, name=cell_types[i])
        
        # ejiaddadd
        for i in range(num_cells):
            for j in range(num_cells):
                if prob_sum[i, j] > 0:
                    G.add_edge(i, j, weight=float(prob_sum[i, j]))
        
        # outinnextnum（addweight/nonaddweight）
        for i in range(num_cells):
            # addweightnextnum
            net_centr['outdeg'][i] = sum(G[i][j]['weight'] for j in G.successors(i)) if i in G else 0
            net_centr['indeg'][i] = sum(G[j][i]['weight'] for j in G.predecessors(i)) if i in G else 0
            
            # nonaddweightnextnum
            net_centr['outdeg_unweighted'][i] = G.out_degree(i) if i in G else 0
            net_centr['indeg_unweighted'][i] = G.in_degree(i) if i in G else 0
        
        # weightmiinvnumgurafu（mediumviamidcenternatureuse）
        G_inv = G.copy()
        for u, v, d in G_inv.edges(data=True):
            if d['weight'] > 0:
                d['weight'] = 1.0 / d['weight']
        
        try:
            # mediumviamidcenternature
            betweenness = nx.betweenness_centrality(G_inv, weight='weight')
            for i in range(num_cells):
                net_centr['betweenness'][i] = betweenness.get(i, 0)
                # flowbetisbetweennesswithsubuse
                net_centr['flowbet'][i] = betweenness.get(i, 0)
        except:
            pass
        
        try:
            # PageRank
            page_rank = nx.pagerank(G, weight='weight')
            for i in range(num_cells):
                net_centr['page_rank'][i] = page_rank.get(i, 0)
        except:
            pass
        
        try:
            # HubandAuthority
            hub, authority = nx.hits(G, max_iter=1000)
            for i in range(num_cells):
                net_centr['hub'][i] = hub.get(i, 0)
                net_centr['authority'][i] = authority.get(i, 0)
        except:
            pass
        
        try:
            # fixhavebekutorumidcenternature
            eigen = nx.eigenvector_centrality(G, weight='weight', max_iter=1000)
            for i in range(num_cells):
                net_centr['eigen'][i] = eigen.get(i, 0)
        except:
            pass
        
        # Information Centrality - kasutamuCalculation（nearlike）
        try:
            # infomidcenternatureofsimpleeasynearlike (Normalizationshitadegreeandbetweennessofsetmimatchwase)
            in_degree = np.array([net_centr['indeg'][i] for i in range(num_cells)])
            out_degree = np.array([net_centr['outdeg'][i] for i in range(num_cells)])
            between = np.array([net_centr['betweenness'][i] for i in range(num_cells)])
            
            # Normalization
            norm_in = in_degree / (np.max(in_degree) if np.max(in_degree) > 0 else 1)
            norm_out = out_degree / (np.max(out_degree) if np.max(out_degree) > 0 else 1)
            norm_between = between / (np.max(between) if np.max(between) > 0 else 1)
            
            # infoinfomidcenternaturethematchbecome
            info_score = (norm_in + norm_out + norm_between) / 3
            for i in range(num_cells):
                net_centr['info'][i] = info_score[i]
        except:
            pass
    
    # midcenternaturepointmarkthenettoaddadd
    result['net']['centr'] = net_centr
    
    # 2. pasuueiinfoinfo
    result['netP'] = {}
    
    # pasuueiinfoinfoofaddadd
    if 'pathways' in cellchat_data and cellchat_data['pathways'] is not None:
        result['netP']['pathways'] = cellchat_data['pathways']
    elif 'net_prob_dims' in cellchat_data and 'pathways' in cellchat_data['net_prob_dims']:
        result['netP']['pathways'] = cellchat_data['net_prob_dims']['pathways']
    else:
        # sorewithmoviewtsufromnotplacematchisnet['dimnames']fromgetget
        if 'dimnames' in result['net'] and len(result['net']['dimnames']) > 2:
            result['netP']['pathways'] = result['net']['dimnames'][2]
        else:
            result['netP']['pathways'] = [f"Pathway_{i}" for i in range(result['net']["prob"].shape[2])]
    
    # Networkcertainrateinfoinfoofaddadd
    if 'netP_prob' in cellchat_data and cellchat_data['netP_prob'] is not None:
        result['netP']['prob'] = cellchat_data['netP_prob']
    else:
        # netP_probisnotplacematch、pasuueigoandofgathercalctherowu
        prob = result['net']["prob"]
        num_pathways = len(result['netP']['pathways'])
        
        # eachpasuueioftotalsumtheCalculation
        netP_prob = np.zeros((num_cells, num_cells, num_pathways))
        for k in range(min(prob.shape[2], num_pathways)):
            netP_prob[:, :, k] = prob[:, :, k]
        
        result['netP']['prob'] = netP_prob
    
    # netPuseofmidcenternatureinfoinfotheCalculation
    netP_centr = {}
    
    # eachpasuueitotsuitemidcenternaturetheCalculation
    for k in range(min(result['netP']['prob'].shape[2], len(result['netP']['pathways']))):
        # koofpasuueiofinitialmidcenternaturewordwrite
        netP_centr[k] = {
            'outdeg': np.zeros(num_cells),
            'indeg': np.zeros(num_cells),
            'outdeg_unweighted': np.zeros(num_cells),
            'indeg_unweighted': np.zeros(num_cells),
            'betweenness': np.zeros(num_cells),
            'page_rank': np.zeros(num_cells),
            'hub': np.zeros(num_cells),
            'authority': np.zeros(num_cells),
            'eigen': np.zeros(num_cells),
            'flowbet': np.zeros(num_cells),
            'info': np.zeros(num_cells)
        }
        
        # pasuueirowcol
        pathway_weight = result['netP']['prob'][:, :, k]
        
        # Networkisvalid（ejiisexist）placematchofmiCalculation
        if np.any(pathway_weight > 0):
            # gurafumakebecome
            G_pathway = nx.DiGraph()
            for i in range(num_cells):
                G_pathway.add_node(i, name=cell_types[i])
            
            # ejiaddadd
            for i in range(num_cells):
                for j in range(num_cells):
                    if pathway_weight[i, j] > 0:
                        G_pathway.add_edge(i, j, weight=float(pathway_weight[i, j]))
            
            # outinnextnum（addweight/nonaddweight）
            for i in range(num_cells):
                # addweightnextnum
                netP_centr[k]['outdeg'][i] = sum(G_pathway[i][j]['weight'] for j in G_pathway.successors(i)) if i in G_pathway else 0
                netP_centr[k]['indeg'][i] = sum(G_pathway[j][i]['weight'] for j in G_pathway.predecessors(i)) if i in G_pathway else 0
                
                # nonaddweightnextnum
                netP_centr[k]['outdeg_unweighted'][i] = G_pathway.out_degree(i) if i in G_pathway else 0
                netP_centr[k]['indeg_unweighted'][i] = G_pathway.in_degree(i) if i in G_pathway else 0
            
            # weightmiinvnumgurafu（mediumviamidcenternatureuse）
            G_inv = G_pathway.copy()
            for u, v, d in G_inv.edges(data=True):
                if d['weight'] > 0:
                    d['weight'] = 1.0 / d['weight']
            
            try:
                # mediumviamidcenternature
                betweenness = nx.betweenness_centrality(G_inv, weight='weight')
                for i in range(num_cells):
                    netP_centr[k]['betweenness'][i] = betweenness.get(i, 0)
                    # flowbetisbetweennesswithsubuse
                    netP_centr[k]['flowbet'][i] = betweenness.get(i, 0)
            except:
                pass
            
            try:
                # PageRank
                page_rank = nx.pagerank(G_pathway, weight='weight')
                for i in range(num_cells):
                    netP_centr[k]['page_rank'][i] = page_rank.get(i, 0)
            except:
                pass
            
            try:
                # HubandAuthority
                hub, authority = nx.hits(G_pathway, max_iter=1000)
                for i in range(num_cells):
                    netP_centr[k]['hub'][i] = hub.get(i, 0)
                    netP_centr[k]['authority'][i] = authority.get(i, 0)
            except:
                pass
            
            try:
                # fixhavebekutorumidcenternature
                eigen = nx.eigenvector_centrality(G_pathway, weight='weight', max_iter=1000)
                for i in range(num_cells):
                    netP_centr[k]['eigen'][i] = eigen.get(i, 0)
            except:
                pass
            
            # Information Centrality - kasutamuCalculation
            try:
                # infomidcenternatureofsimpleeasynearlike (Normalizationshitadegreeandbetweennessofsetmimatchwase)
                in_degree = np.array([netP_centr[k]['indeg'][i] for i in range(num_cells)])
                out_degree = np.array([netP_centr[k]['outdeg'][i] for i in range(num_cells)])
                between = np.array([netP_centr[k]['betweenness'][i] for i in range(num_cells)])
                
                # Normalization
                norm_in = in_degree / (np.max(in_degree) if np.max(in_degree) > 0 else 1)
                norm_out = out_degree / (np.max(out_degree) if np.max(out_degree) > 0 else 1)
                norm_between = between / (np.max(between) if np.max(between) > 0 else 1)
                
                # infoinfomidcenternaturethematchbecome
                info_score = (norm_in + norm_out + norm_between) / 3
                for i in range(num_cells):
                    netP_centr[k]['info'][i] = info_score[i]
            except:
                pass
    
    # midcenternatureinfoinfotheaddadd
    result['netP']['centr'] = netP_centr
    
    # 3. ResultDataFramemakebecome
    # Cell typenameandpasuueinamethegetget
    cell_types_rows = result['net']['dimnames'][0]
    cell_types_cols = result['net']['dimnames'][1]
    pathways = result['net']['dimnames'][2]
    
    # prob, pvaldistcolfromInteractionDatatheextractout
    prob = result['net']["prob"]
    pval = result['net']["pval"]
    
    results_data = []
    for i in range(prob.shape[0]):
        for j in range(prob.shape[1]):
            for k in range(min(prob.shape[2], len(pathways))):
                if prob[i, j, k] > 0:
                    # Interactionname
                    interaction_name = pathways[k] if k < len(pathways) else f"interaction_{k}"
                    
                    # Ligand-reseputa-infoinfo（areba）
                    ligand = ""
                    receptor = ""
                    
                    # LRinfoinfogetget
                    if 'lr_sig' in cellchat_data and isinstance(cellchat_data['lr_sig'], pd.DataFrame):
                        lr_df = cellchat_data['lr_sig']
                        if interaction_name in lr_df.index:
                            if 'ligand' in lr_df.columns:
                                ligand = lr_df.loc[interaction_name, 'ligand']
                            if 'receptor' in lr_df.columns:
                                receptor = lr_df.loc[interaction_name, 'receptor']
                    
                    results_data.append({
                        'source': cell_types_rows[i] if i < len(cell_types_rows) else f"Cell_{i}",
                        'target': cell_types_cols[j] if j < len(cell_types_cols) else f"Cell_{j}",
                        'interaction_name': interaction_name,
                        'ligand': str(ligand),
                        'receptor': str(receptor),
                        'prob': float(prob[i, j, k]),
                        'pval': float(pval[i, j, k])
                    })
    
    result['results'] = pd.DataFrame(results_data)
    
    # 4. Networksamari-makebecome
    # Cell typegoandofgathercalc（probofmatchcalc）
    prob = result['net']["prob"]
    strength_matrix = np.sum(prob, axis=2)
    
    # DataFrametochangechange
    strength_df = pd.DataFrame(strength_matrix, 
                              index=cell_types, 
                              columns=cell_types)
    
    # kauntorowcol（nonzeroofInteractionnum）
    count_matrix = np.sum(prob > 0, axis=2)
    count_df = pd.DataFrame(count_matrix,
                           index=cell_types,
                           columns=cell_types)
    
    # sendrecvtrusttotalamount
    outgoing = pd.Series(np.sum(strength_matrix, axis=1), index=cell_types)
    incoming = pd.Series(np.sum(strength_matrix, axis=0), index=cell_types)
    
    # LRcontribgivedegreeofCalculation
    lr_contribution = np.zeros(prob.shape[2])
    for k in range(prob.shape[2]):
        lr_contribution[k] = np.sum(prob[:, :, k])
    
    # NetworkmidcenternatureDatafure-muthemakebecome
    network_centrality_dict = {
        'cell_type': cell_types,
    }
    
    # eachmidcenternaturepointmarkthedeikushiyonaritoaddadd
    for metric in ['outdeg', 'indeg', 'outdeg_unweighted', 'indeg_unweighted', 
                  'betweenness', 'page_rank', 'hub', 'authority', 'eigen', 'flowbet', 'info']:
        network_centrality_dict[metric] = net_centr[metric]
    
    network_centrality = pd.DataFrame(network_centrality_dict)
    
    # koamidcenternaturepointmarkthethanmeanmeanofexistnamebeforetochangefurther（shigunaringuroleratiodivanalyzeuse）
    role_data = {
        'cell_type': cell_types,
        'sender': net_centr['outdeg'],        # sendtrustpersonandshiteofroleratio
        'receiver': net_centr['indeg'],       # recvtrustpersonandshiteofroleratio
        'mediator': net_centr['flowbet'],     # 仲viapersonandshiteofroleratio
        'influencer': net_centr['info']       # shadoweffectpowerandshiteofroleratio
    }
    
    # roleratioDatafure-muthemakebecome
    signaling_role = pd.DataFrame(role_data)
    
    result['network'] = {
        'strength_matrix': strength_df,
        'count_matrix': count_df,
        'outgoing': outgoing,
        'incoming': incoming,
        'lr_contribution': lr_contribution,
        'network_centrality': network_centrality,
        'signaling_role': signaling_role      # shigunaringuroleratiodivanalyzeuse
    }
    
    # 5. metaDatainfoinfo
    if 'meta' in cellchat_data:
        result['meta'] = cellchat_data['meta']
    
    # 6. GroupByinfoinfo（defuorutoSettings）
    result['groupby'] = 'cell.ident'  # defuorutoval
    
    # 7. LRinfoinfo
    if 'lr_sig' in cellchat_data:
        result['lr_sig'] = cellchat_data['lr_sig']
    
    # 8. adatainfoinfo - doubtlikealnaanndatastructmakethemakebecome
    from types import SimpleNamespace
    
    # adataobujiekutothemakebecome
    adata = SimpleNamespace()
    
    # metaDatatheaddadd
    if 'meta' in cellchat_data and isinstance(cellchat_data['meta'], pd.DataFrame):
        adata.obs = cellchat_data['meta']
    else:
        adata.obs = pd.DataFrame(index=cell_types)
    
    # Genenamerisutotheextractout
    gene_names = []
    
    # LRsigfromLigandandreseputa-infoinfotheextractout
    if 'lr_sig' in cellchat_data and isinstance(cellchat_data['lr_sig'], pd.DataFrame):
        lr_df = cellchat_data['lr_sig']
        if 'ligand' in lr_df.columns:
            ligands = lr_df['ligand'].dropna().astype(str).unique()
            gene_names.extend(ligands)
        if 'receptor' in lr_df.columns:
            receptors = lr_df['receptor'].dropna().astype(str).unique()
            gene_names.extend(receptors)
    
    # Genenameisnotplacematchoffuo-rubaku
    if not gene_names:
        gene_names = ["Gene1", "Gene2", "Gene3"]
    
    # weightmultithedeleteremove
    gene_names = sorted(list(set(gene_names)))
    
    # AnnDataofvarbelongnatureandvar_namesthemodel倣
    adata.var_names = pd.Index(gene_names)
    adata.var = pd.DataFrame(index=adata.var_names)
    
    # AnnDataofshapebelongnaturetheaddadd
    adata.shape = (len(adata.obs), len(adata.var_names))
    
    # Xbelongnaturetheaddadd
    adata.X = np.zeros(adata.shape)
    
    result['adata'] = adata
    
    return result

# meinprocproc
if r_initialized:
    # Fileapuro-da-
    uploaded_file = st.file_uploader("CellChat QSFiletheUpload", type=["qs"])
    
    if uploaded_file is not None:
        st.info(f"File「{uploaded_file.name}」theprocprocshiteimasu...")
        
        with st.spinner("QSFiletheProcessing..."):
            result = process_qs_file(uploaded_file)
        
        if result is not None:
            st.success("QSFileofchangechangeisCompleteshimashita！")
            
            # ResultofoverviewneedtheDisplay
            st.subheader("changechangeResultofoverviewneed")
            
            # Cell typeinfoinfoofDisplay
            if 'net' in result and 'dimnames' in result['net'] and len(result['net']['dimnames']) > 0:
                cell_types = result['net']['dimnames'][0]
                cell_count = len(cell_types)
                if cell_count > 0:
                    st.write(f"checkoutsaretaCell typenum: {cell_count}")
                    if cell_count <= 20:
                        st.write(f"Cell type: {', '.join(cell_types)}")
                    else:
                        with st.expander("Cell typeonebrowse"):
                            st.write(", ".join(cell_types))
            
            # pasuueiinfoinfoofDisplay
            if 'netP' in result and 'pathways' in result['netP'] and result['netP']['pathways'] is not None:
                pathways = result['netP']['pathways']
                pathway_count = len(pathways)
                if pathway_count > 0:
                    st.write(f"shigunaruPathwaynum: {pathway_count}")
                    if pathway_count <= 10:
                        st.write(f"shigunaruPathway: {', '.join(pathways)}")
                    else:
                        with st.expander("shigunaruPathwayonebrowse"):
                            st.write(", ".join(pathways))
            
            # InteractioninfoinfoofDisplay
            if 'results' in result and not result['results'].empty:
                st.write(f"extractoutsaretaInteractionnum: {len(result['results'])}")
                with st.expander("InteractionDataofSample（uprank5item）"):
                    st.dataframe(result['results'].head())
            
            # ResultthePickleFiletoSave
            output_filename = os.path.splitext(uploaded_file.name)[0] + "_python.pkl"
            pickle_data = pickle.dumps(result)
            
            # Downloadbotan
            st.download_button(
                label="changechangeResulttheDownload",
                data=pickle_data,
                file_name=output_filename,
                mime="application/octet-stream"
            )
            
            st.info("""
            DownloadshitapickleFileis、cellchat.pywithCellChatofVisualizationtouseusewithkimasu。
            StreamlitapuritoUploaddoka、ko-doinwithpickle.load()theuseuseshitero-doshitekudasai。
            """)
else:
    st.error("Rringboundofinitializetofailfailshitafor、apuritheRunwithkimasen。requirednapake-jiisinsuto-rusareteexistkaConfirmshitekudasai。")