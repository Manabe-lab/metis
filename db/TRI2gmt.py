import pandas as pd
from collections import defaultdict

def tri_to_gmt(input_file, output_file):
    """
    TRIフォーマットのTSVファイルをGMTフォーマットに変換する
    Positive/Negative の制御関係を別々の遺伝子セットとして出力する
    
    Parameters:
    -----------
    input_file : str
        入力のTRIファイルパス（TSV形式）
    output_file : str
        出力するGMTファイルパス
    """
    # TRIファイルを読み込む
    tri_df = pd.read_csv(input_file, sep='\t')
    
    # sourceごとにtargetをグループ化（正と負で別々に）
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
    
    # GMTファイルを書き出す
    with open(output_file, 'w') as f:
        # Positive targets
        for source, targets in positive_regulons.items():
            if targets:  # ターゲットが存在する場合のみ出力
                description = f"{source}_positive_targets"
                target_genes = '\t'.join(targets)
                line = f"{source}_pos\t{description}\t{target_genes}\n"
                f.write(line)
        
        # Negative targets
        for source, targets in negative_regulons.items():
            if targets:  # ターゲットが存在する場合のみ出力
                description = f"{source}_negative_targets"
                target_genes = '\t'.join(targets)
                line = f"{source}_neg\t{description}\t{target_genes}\n"
                f.write(line)

    # 統計情報を表示
    print(f"変換が完了しました。出力ファイル: {output_file}")
    print(f"\n統計情報:")
    print(f"Positive 制御の遺伝子セット数: {len(positive_regulons)}")
    print(f"Negative 制御の遺伝子セット数: {len(negative_regulons)}")
    
    if positive_regulons:
        avg_pos = sum(len(targets) for targets in positive_regulons.values()) / len(positive_regulons)
        print(f"Positive 制御の平均ターゲット遺伝子数: {avg_pos:.1f}")
    
    if negative_regulons:
        avg_neg = sum(len(targets) for targets in negative_regulons.values()) / len(negative_regulons)
        print(f"Negative 制御の平均ターゲット遺伝子数: {avg_neg:.1f}")

# 使用例
if __name__ == "__main__":
    input_file = "TRI.mouse.tsv"
    output_file = "TRI.mouse.gmt"
    tri_to_gmt(input_file, output_file)
    input_file = "TRI.human.tsv"
    output_file = "TRI.human.gmt"
    tri_to_gmt(input_file, output_file)
