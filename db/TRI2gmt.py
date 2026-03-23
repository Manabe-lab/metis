import pandas as pd
from collections import defaultdict

def tri_to_gmt(input_file, output_file):
    """
    Convert a TRI format TSV file to GMT format.
    Output positive/negative regulatory relationships as separate gene sets.

    Parameters:
    -----------
    input_file : str
        Input TRI file path (TSV format)
    output_file : str
        Output GMT file path
    """
    # Read the TRI file
    tri_df = pd.read_csv(input_file, sep='\t')
    
    # Group targets by source (separately for positive and negative)
    positive_regulons = defaultdict(list)
    negative_regulons = defaultdict(list)
    
    for _, row in tri_df.iterrows():
        source = row['source']
        target = row['target']
        weight = row['weight']
        
        if weight > 0:
            positive_regulons[source].append(target)
        elif weight < 0:
            negative_regulons[source].append(target)
    
    # Write the GMT file
    with open(output_file, 'w') as f:
        # Positive targets
        for source, targets in positive_regulons.items():
            if targets:  # Only output if targets exist
                description = f"{source}_positive_targets"
                target_genes = '\t'.join(targets)
                line = f"{source}_pos\t{description}\t{target_genes}\n"
                f.write(line)
        
        # Negative targets
        for source, targets in negative_regulons.items():
            if targets:  # Only output if targets exist
                description = f"{source}_negative_targets"
                target_genes = '\t'.join(targets)
                line = f"{source}_neg\t{description}\t{target_genes}\n"
                f.write(line)

    # Display statistics
    print(f"Conversion complete. Output file: {output_file}")
    print(f"\nStatistics:")
    print(f"Number of positive regulatory gene sets: {len(positive_regulons)}")
    print(f"Number of negative regulatory gene sets: {len(negative_regulons)}")
    
    if positive_regulons:
        avg_pos = sum(len(targets) for targets in positive_regulons.values()) / len(positive_regulons)
        print(f"Average target genes per positive regulatory set: {avg_pos:.1f}")
    
    if negative_regulons:
        avg_neg = sum(len(targets) for targets in negative_regulons.values()) / len(negative_regulons)
        print(f"Average target genes per negative regulatory set: {avg_neg:.1f}")

# Usage example
if __name__ == "__main__":
    input_file = "TRI.mouse.tsv"
    output_file = "TRI.mouse.gmt"
    tri_to_gmt(input_file, output_file)
    input_file = "TRI.human.tsv"
    output_file = "TRI.human.gmt"
    tri_to_gmt(input_file, output_file)
