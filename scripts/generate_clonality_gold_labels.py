import argparse
from pathlib import Path
import pandas as pd


def _ensure_not_gold_output(output_csv: Path) -> None:
    if output_csv.name == "clonality_gold_labels.csv":
        raise ValueError(
            "Refusing to write pipeline-generated labels to clonality_gold_labels.csv. "
            "That file is reserved for human-reviewed labels."
        )


def generate_labels(ladder_csv: Path, pk_csv: Path, output_csv: Path):
    _ensure_not_gold_output(output_csv)
    dfs = []
    
    if ladder_csv.exists():
        ladder = pd.read_csv(ladder_csv)
        if not ladder.empty:
            l_gold = pd.DataFrame()
            l_gold['artifact_table'] = ladder['artifact_table']
            l_gold['artifact_row_key'] = ladder['artifact_row_key']
            l_gold['label'] = ladder['selected_for_fit'].apply(lambda x: 'pipeline_accept' if str(x).lower() == 'true' else 'pipeline_reject')
            l_gold['label_source'] = 'pipeline_pseudo_label'
            l_gold['reviewer'] = 'pipeline'
            l_gold['reviewed_at_utc'] = pd.Timestamp.utcnow().isoformat()
            l_gold['notes'] = 'Pseudo-label from current pipeline selection; not human gold.'
            dfs.append(l_gold)
            
    if pk_csv.exists():
        pk = pd.read_csv(pk_csv)
        if not pk.empty:
            p_gold = pd.DataFrame()
            p_gold['artifact_table'] = pk['artifact_table']
            p_gold['artifact_row_key'] = pk['artifact_row_key']
            p_gold['label'] = pk['selected'].apply(lambda x: 'pipeline_selected' if str(x).lower() == 'true' else 'pipeline_reject')
            p_gold['label_source'] = 'pipeline_pseudo_label'
            p_gold['reviewer'] = 'pipeline'
            p_gold['reviewed_at_utc'] = pd.Timestamp.utcnow().isoformat()
            p_gold['notes'] = 'Pseudo-label from current pipeline selection; not human gold.'
            dfs.append(p_gold)
            
    if dfs:
        res = pd.concat(dfs, ignore_index=True)
        res.to_csv(output_csv, index=False)
        print(f"Generated {len(res)} gold labels at {output_csv}")
    else:
        print("No candidates found to generate labels from.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ladder-csv', type=Path, required=True)
    parser.add_argument('--pk-csv', type=Path, required=True)
    parser.add_argument('--output-csv', type=Path, required=True, help='Pseudo-label CSV output. Not for human gold labels.')
    args = parser.parse_args()
    generate_labels(args.ladder_csv, args.pk_csv, args.output_csv)
