import streamlit as st
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import requests
import io
from urllib.parse import urlencode

def process_gene_list(gene_input):
    """Function to process gene list"""
    # Replace all newlines with spaces
    genes = gene_input.replace('\n', ' ')

    # Replace semicolons with spaces
    genes = genes.replace(';', ' ')

    # Replace consecutive spaces with a single space
    genes = ' '.join(genes.split())

    # Split by spaces to create list
    gene_list = [g.strip() for g in genes.split()]

    # Remove empty elements
    gene_list = [g for g in gene_list if g]

    # Remove duplicates while maintaining order
    seen = set()
    unique_genes = []
    for gene in gene_list:
        if gene not in seen:
            seen.add(gene)
            unique_genes.append(gene)

    return unique_genes

def fetch_string_db_interactions(genes, species=10090, score_threshold=0):
    """Get protein-protein interaction data from StringDB"""
    base_url = "https://string-db.org/api/json"
    endpoint = "network"

    # Set parameters
    params = {
        "identifiers": "\r".join(genes),
        "species": species,
        "required_score": 0,  # Score threshold applied later
        "caller_identity": "your_app_name"
    }

    try:
        # Send API request
        api_url = f"{base_url}/{endpoint}"
        st.write(f"Requesting URL: {api_url}")

        response = requests.post(api_url, data=params)
        response.raise_for_status()

        # Display debug information
        st.write("API Response Status:", response.status_code)
        st.write("API Response Content Type:", response.headers.get('content-type'))

        # Try to parse JSON response
        if response.text.strip():
            try:
                data = response.json()
                if isinstance(data, list):
                    st.write(f"Found {len(data)} interactions")
                    return data
                else:
                    st.warning("Unexpected response format")
                    st.write("Response structure:", type(data))
                    return []
            except ValueError as e:
                st.error(f"Failed to parse JSON response: {str(e)}")
                st.write("Raw response:", response.text[:500])
                return []
        else:
            st.warning("Empty response from STRING-DB")
            return []

    except requests.exceptions.RequestException as e:
        st.error(f"Error fetching data from STRING-DB: {str(e)}")
        if hasattr(response, 'content'):
            st.error(f"Response content: {response.content.decode()[:500]}")
        return []
    except Exception as e:
        st.error(f"Unexpected error: {str(e)}")
        return []

def process_interactions(interactions, G, genes, score_threshold):
    """Process interaction data and add edges"""
    edge_count = 0
    skipped_edges = 0
    all_scores = []

    # Check score distribution
    scores = [float(interaction.get('score', 0)) for interaction in interactions]
    if scores:
        st.write(f"Score distribution - Min: {min(scores)}, Max: {max(scores)}, Mean: {np.mean(scores):.2f}")

    for interaction in interactions:
        source = interaction.get('preferredName_A')
        target = interaction.get('preferredName_B')
        
        try:
            score = float(interaction.get('score', 0))
            all_scores.append(score)
            
            if score >= score_threshold:
                # Add nodes as needed
                if source not in G:
                    G.add_node(source)
                if target not in G:
                    G.add_node(target)
                
                G.add_edge(source, target, weight=score)
                edge_count += 1
            else:
                skipped_edges += 1
        except ValueError:
            skipped_edges += 1
            continue

    st.write(f"Added {edge_count} edges to the network")
    st.write(f"Skipped {skipped_edges} edges due to score threshold")

    # Display histogram of score distribution
    if all_scores:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.hist(all_scores, bins=20, edgecolor='black')
        ax.axvline(x=score_threshold, color='red', linestyle='--', label='Threshold')
        ax.set_xlabel('Interaction Score')
        ax.set_ylabel('Frequency')
        ax.set_title('Distribution of Interaction Scores')
        ax.legend()
        st.pyplot(fig)
        plt.close()

    return edge_count

def visualize_network(G, minimum_score):
    """Function to visualize network"""
    plt.figure(figsize=(15, 10))

    # Calculate node positions with spring layout
    pos = nx.spring_layout(G, k=1.5/np.sqrt(len(G.nodes())), iterations=50)

    # Draw edges
    if G.number_of_edges() > 0:
        edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
        max_weight = max(edge_weights)
        min_weight = min(edge_weights)

        # Normalize edge thickness (range 0.5-5)
        if max_weight != min_weight:
            normalized_weights = [0.5 + 4.5 * (w - min_weight) / (max_weight - min_weight)
                                for w in edge_weights]
        else:
            normalized_weights = [1.0 for _ in edge_weights]

        nx.draw_networkx_edges(G, pos,
                             width=normalized_weights,
                             alpha=0.5,
                             edge_color='gray')

    # Draw nodes
    degrees = dict(G.degree())
    max_degree = max(degrees.values()) if degrees else 1

    # Adjust node size (minimum 1000, maximum 3000)
    node_sizes = [1000 + 2000 * (degrees[node] / max_degree) for node in G.nodes()]

    nx.draw_networkx_nodes(G, pos,
                          node_size=node_sizes,
                          node_color='lightblue',
                          alpha=0.7)

    # Draw labels (adjust font size based on node degree)
    font_sizes = {node: min(8 + 4 * (degrees[node] / max_degree), 14)
                 for node in G.nodes()}
    
    nx.draw_networkx_labels(G, pos, 
                          font_size=8,
                          font_weight='bold')
    
    plt.title("Protein-Protein Interaction Network")
    plt.axis('off')
    
    return plt.gcf()

def find_optimal_threshold(interactions, min_edges=10):
    """Find the optimal score threshold"""
    if not interactions:
        return 400  # Default value

    # Get all available scores
    scores = []
    for interaction in interactions:
        try:
            score = float(interaction.get('score', 0))
            scores.append(score)
        except (ValueError, TypeError):
            continue

    if not scores:
        return 400

    # Sort scores in descending order
    scores.sort(reverse=True)

    # Find the maximum score that guarantees at least min_edges edges
    for score in scores:
        edge_count = sum(1 for s in scores if s >= score)
        if edge_count >= min_edges:
            return score

    return min(scores)

def analyze_network(G):
    """Function to analyze network"""
    analysis_results = {
        "basic_stats": {
            "nodes": G.number_of_nodes(),
            "edges": G.number_of_edges(),
            "average_degree": 2 * G.number_of_edges() / G.number_of_nodes() if G.number_of_nodes() > 0 else 0
        }
    }

    if G.number_of_nodes() > 0:
        # Calculate centrality
        analysis_results["centrality"] = {
            "degree": nx.degree_centrality(G),
            "betweenness": nx.betweenness_centrality(G) if G.number_of_edges() > 0 else {},
            "closeness": nx.closeness_centrality(G) if G.number_of_edges() > 0 else {}
        }

        # Clustering coefficient (only if edges exist)
        if G.number_of_edges() > 0:
            analysis_results["clustering"] = nx.clustering(G)
        else:
            analysis_results["clustering"] = {node: 0 for node in G.nodes()}

    return analysis_results

def main():
    st.title("Protein-Protein Interaction Network Analysis using STRING-db")
    
    # Gene list input
    gene_input = st.text_area(
        "Enter gene symbols (separated by spaces, newlines, or semicolons):",
        height=200,
        help="You can paste gene symbols separated by spaces, newlines, or semicolons"
    )

    # Process gene list
    genes = process_gene_list(gene_input)

    # Display processed gene list
    st.write(f"Number of unique genes: {len(genes)}")
    if st.checkbox("Show processed gene list"):
        st.write(", ".join(genes))
    



    col1, col2 = st.columns(2)
    with col1:
        species = st.selectbox(
            "Species:", 
            options=[(10090, "Mus musculus"), (9606, "Homo sapiens")],
            format_func=lambda x: x[1]
        )
        species_id = species[0]
    
    with col2:
        auto_threshold = st.checkbox("Auto-adjust score threshold", value=True)
        if auto_threshold:
            score_threshold = st.slider(
                "Minimum Score Threshold:", 
                min_value=0.0, max_value=1.0, value=0.4,
                step=0.05,
                help="This will be automatically adjusted if needed"
            )
        else:
            score_threshold = st.slider(
                "Score Threshold:", 
                min_value=0.0, max_value=1.0, value=0.4,
                step=0.05,
                help="Higher values mean more stringent interaction criteria (0-1)"
            )
    
    # Add filtering options
    with st.expander("Advanced Filtering Options"):
        max_edges = st.number_input(
            "Maximum number of edges to show",
            min_value=10,
            max_value=1000,
            value=200,
            help="Limit the number of edges to prevent overcrowded visualization"
        )
        
        show_strongest = st.checkbox(
            "Show only strongest interactions",
            value=True,
            help="If checked, will show only the strongest interactions up to the maximum number"
        )
    
    if st.button("Analyze Network"):
        with st.spinner("Fetching protein interaction data..."):
            # Get all interactions with score threshold 0
            interactions = fetch_string_db_interactions(genes, species_id, 0)
            
            if not interactions:
                st.error("No interactions found for the given parameters.")
                return
            
            # Display sample of interaction data
            with st.expander("Show raw interaction data sample"):
                sample_df = pd.DataFrame(interactions[:5])
                st.dataframe(sample_df)

            # Create network using NetworkX
            G = nx.Graph()

            # Add nodes
            for gene in genes:
                G.add_node(gene)

            # Sort interactions by score
            if show_strongest:
                interactions = sorted(
                    interactions,
                    key=lambda x: float(x.get('score', 0)),
                    reverse=True
                )[:max_edges]

            # Add edges
            edge_count = process_interactions(interactions, G, genes, score_threshold)
            
            if edge_count == 0:
                st.warning("No edges in the network. Trying with lower threshold...")
                if auto_threshold:
                    score_threshold = max(0.0, score_threshold - 0.1)
                    st.info(f"Automatically lowering score threshold to {score_threshold:.2f}")
                    edge_count = process_interactions(interactions, G, genes, score_threshold)
    
            if edge_count > 0:
                try:
                    # Network visualization
                    fig = visualize_network(G, score_threshold)
                    st.pyplot(fig)

                    # Network analysis
                    analysis_results = analyze_network(G)

                    # Display basic statistics
                    st.subheader("Network Statistics")
                    stats = analysis_results["basic_stats"]
                    st.write(f"Number of nodes: {stats['nodes']}")
                    st.write(f"Number of edges: {stats['edges']}")
                    st.write(f"Average degree: {stats['average_degree']:.2f}")

                    # Centrality analysis and display
                    st.subheader("Centrality Analysis")

                    # Top 10 by degree centrality
                    degree_df = pd.DataFrame.from_dict(
                        analysis_results["centrality"]["degree"],
                        orient='index',
                        columns=['Degree Centrality']
                    )
                    degree_df = degree_df.sort_values('Degree Centrality', ascending=False)

                    # Top 10 by betweenness centrality
                    betweenness_df = pd.DataFrame.from_dict(
                        analysis_results["centrality"]["betweenness"],
                        orient='index',
                        columns=['Betweenness Centrality']
                    )
                    betweenness_df = betweenness_df.sort_values('Betweenness Centrality', ascending=False)

                    # Combine centrality results
                    centrality_df = pd.concat([
                        degree_df['Degree Centrality'],
                        betweenness_df['Betweenness Centrality']
                    ], axis=1)
                    
                    st.write("Top 10 genes by centrality measures:")
                    st.dataframe(centrality_df.head(10))
                    
                    # Display clustering coefficient
                    if "clustering" in analysis_results:
                        st.subheader("Clustering Analysis")
                        clustering_df = pd.DataFrame.from_dict(
                            analysis_results["clustering"],
                            orient='index',
                            columns=['Clustering Coefficient']
                        )
                        clustering_df = clustering_df.sort_values('Clustering Coefficient', ascending=False)

                        st.write("Top 10 genes by clustering coefficient:")
                        st.dataframe(clustering_df.head(10))

                    # Export results
                    st.subheader("Export Results")

                    # Prepare edge data
                    edge_df = pd.DataFrame([
                        {
                            'Source': u,
                            'Target': v,
                            'Score': d['weight']
                        }
                        for u, v, d in G.edges(data=True)
                    ])

                    # Prepare node data
                    node_df = centrality_df.reset_index()
                    node_df.columns = ['Gene', 'Degree_Centrality', 'Betweenness_Centrality']

                    # Provide download buttons for data
                    col1, col2 = st.columns(2)
                    with col1:
                        csv_edges = edge_df.to_csv(index=False)
                        st.download_button(
                            label="Download Edge Data",
                            data=csv_edges,
                            file_name="network_edges.csv",
                            mime="text/csv"
                        )
                    
                    with col2:
                        csv_nodes = node_df.to_csv(index=False)
                        st.download_button(
                            label="Download Node Data",
                            data=csv_nodes,
                            file_name="network_nodes.csv",
                            mime="text/csv"
                        )
                
                except Exception as e:
                    st.error(f"Error in network analysis: {str(e)}")
                    import traceback
                    st.error(f"Traceback: {traceback.format_exc()}")
            else:
                st.error("Could not create network with current parameters. Please try adjusting the score threshold.")

if __name__ == "__main__":
    main()