import pandas as pd
import argparse

def get_tss(path_to_file):
    return pd.read_csv(path_to_file)

def get_methylation(path_to_file):
    df = pd.read_csv(path_to_file)
    df['unity_position'] = df['unity_position'].astype(int)
    return df

def find_hits(df_tss, df_methylation):
    promoter_hits = []

    df_tss['tss_position'] = df_tss['tss_position'].astype(int)

    for _, row in df_tss.iterrows():
        tss = int(row['tss_position'])
        strand = row['strand']
        
        if strand == "+":
            promoter_start = tss - 50
            promoter_end = tss
        elif strand == "-":
            promoter_start = tss
            promoter_end = tss + 50
        else:
            continue
        
        hits = df_methylation[
            (df_methylation['unity_position'] >= promoter_start) &
            (df_methylation['unity_position'] <= promoter_end)
        ]

        for _, hit in hits.iterrows():
            promoter_hits.append({
                'tss_position': tss,
                'strand': strand,
                'unity_position': hit['unity_position']
            })

    return pd.DataFrame(promoter_hits)

def main():
    parser = argparse.ArgumentParser(
        description="""Finds methylation sites in promoter regions"""
    )
    parser.add_argument(
        "-t", "--transcription_start_sites",
        required=True,
        help="Path to a CSV file containing transcription start sites"
    )
    parser.add_argument(
        "-m", "--methylation_sites",
        required=True,
        help="Path to a CSV file containing methylation sites"
    )
    parser.add_argument(
        "-o", "--output",
        default="methylation_in_tss.csv",
        help="Output file name"
    )
    args = parser.parse_args()
    

    df_tss = get_tss(args.transcription_start_sites)
    df_methylation = get_methylation(args.methylation_sites)

    result_df = find_hits(df_tss, df_methylation)
    result_df.to_csv(args.output, index=False)
    print(f"Done! Results written to: {args.output}")

if __name__ == "__main__":
    main()
