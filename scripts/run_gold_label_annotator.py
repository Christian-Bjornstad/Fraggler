import argparse
import json
import logging
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import panel as pn
import plotly.graph_objects as go

# Setup Panel context
pn.extension('plotly', sizing_mode="stretch_width")

logging.basicConfig(level=logging.INFO)

# =========================================================
# CONFIG
# =========================================================

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--candidates-dir",
        type=Path,
        required=True,
        help="Directory containing the candidate CSVs and gold labels",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        required=True,
        help="Root directory where the FSA files are stored (e.g. 2025_data)",
    )
    parser.add_argument(
        "--port", type=int, default=5006, help="Port for the panel server"
    )
    return parser

# =========================================================
# FUZZY FILE RESOLUTION
# =========================================================

def resolve_fsa_path(data_dir: Path, identity_key: str) -> Path:
    """Resolve FSA path even if filename has pipeline prefixes or hashes."""
    parts = identity_key.split("::")
    if len(parts) >= 2:
        source_dir = parts[0]
        fname = list(parts)[-1]
    else:
        source_dir = ""
        fname = identity_key
        
    base_dir = data_dir / source_dir
    direct_path = base_dir / fname
    if direct_path.exists():
        return direct_path
        
    # Try stripping 00000_xxxxxxxx_ prefix if it exists
    cleaned_name = re.sub(r'^[0-9]{5}_[a-f0-9]{8}_', '', fname)
    fuzzy_path = base_dir / cleaned_name
    if fuzzy_path.exists():
        return fuzzy_path
        
    # Last ditch: search for anything ending with the cleaned name or similar
    if base_dir.exists():
        for candidate in base_dir.glob("*.fsa"):
            if candidate.name.endswith(cleaned_name) or cleaned_name in candidate.name:
                return candidate
                
    return direct_path

# =========================================================
# LAZY FSA LOADING
# =========================================================

def load_fsa_full_analysis(fsa_path: Path):
    """Load FSA and run full ladder fit (ROX/LIZ) using the official API."""
    import sys
    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
        
    from fraggler import FsaFile
    from core.analysis import analyse_fsa_rox, analyse_fsa_liz
    
    if not fsa_path.exists():
        raise FileNotFoundError(f"Missing FSA file: {fsa_path}")
    
    # Try ROX fit first
    fsa = analyse_fsa_rox(fsa_path, sample_channel="DATA1")
    if fsa is None or not getattr(fsa, "fitted_to_model", False):
        fsa = analyse_fsa_liz(fsa_path, sample_channel="DATA1")
        
    if fsa is None or not getattr(fsa, "fitted_to_model", False):
        raise ValueError(f"Could not fit ladder for {fsa_path.name}")
        
    return fsa

# =========================================================
# DASHBOARD APP
# =========================================================

class GoldLabelAnnotator:
    def __init__(self, candidates_dir: Path, data_dir: Path):
        self.candidates_dir = candidates_dir
        self.data_dir = data_dir
        self.gold_path = candidates_dir / "clonality_gold_labels.csv"
        
        # Load tables
        ladder_path = candidates_dir / "clonality_ladder_candidates.csv"
        pk_path = candidates_dir / "clonality_pk_candidates.csv"
        
        self.ladder_df = pd.read_csv(ladder_path) if ladder_path.exists() else pd.DataFrame()
        self.pk_df = pd.read_csv(pk_path) if pk_path.exists() else pd.DataFrame()
        
        # Normalize Ladder DF for unified curation
        if not self.ladder_df.empty:
            # Map ladder fields to PK fields for the UI
            self.ladder_df = self.ladder_df.rename(columns={
                "candidate_time": "found_bp", 
                "candidate_intensity": "height",
                "ladder": "channel"
            })
            self.ladder_df['marker_name'] = "Ladder Check"
            # Ensure channel is string like "DATA105" or "ROX"
            self.ladder_df['channel'] = self.ladder_df['channel'].fillna("DATA105").astype(str)
            
        # Merge PK and Ladder candidates for unified curation
        self.all_candidates = pd.concat([self.pk_df, self.ladder_df], ignore_index=True)
        # Final safety on channel
        self.all_candidates['channel'] = self.all_candidates['channel'].fillna("UNKNOWN").astype(str)
        
        if not self.gold_path.exists():
            pd.DataFrame(columns=[
                "artifact_table", "artifact_row_key", "label", "label_source", 
                "reviewer", "reviewed_at_utc", "notes"
            ]).to_csv(self.gold_path, index=False)
            
        self.gold_df = pd.read_csv(self.gold_path) if self.gold_path.exists() else pd.DataFrame(columns=["artifact_table", "artifact_row_key", "label", "label_source", "reviewer", "reviewed_at_utc", "notes"])

        # Unified candidates for grouping
        self.all_candidates = pd.concat([self.ladder_df, self.pk_df], ignore_index=True)

        # Create Dir if it doesn't exist
        self.gold_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Prepare Tasks (Group BY FILE identity_key)
        self.tasks = self._build_tasks()
        self.current_idx = 0
        
        # UI Components
        self.title_ui = pn.pane.Markdown("### 🔍 Fraggler Gold Label Annotator (Multi-Trace View)")
        self.progress_ui = pn.pane.Markdown("")
        self.info_ui = pn.pane.Markdown("")
        
        # Container for stacked plots
        self.plot_container = pn.Column(sizing_mode="stretch_width")
        
        # Buttons
        self.btn_accept = pn.widgets.Button(name="Correct (Selected)", button_type="success", sizing_mode="stretch_width")
        self.btn_reject = pn.widgets.Button(name="Reject File", button_type="danger", sizing_mode="stretch_width")
        self.btn_skip = pn.widgets.Button(name="Skip ⏭️", button_type="primary", sizing_mode="stretch_width")
        self.btn_prev = pn.widgets.Button(name="⏪ Previous", button_type="default", sizing_mode="stretch_width")
        
        self.note_input = pn.widgets.TextInput(placeholder="Optional note...")
        
        # Callbacks
        self.btn_accept.on_click(lambda event: self.save_decision('selected_correct'))
        self.btn_reject.on_click(lambda event: self.save_decision('reject'))
        self.btn_skip.on_click(lambda event: self.next_task())
        self.btn_prev.on_click(lambda event: self.prev_task())

        self.layout = pn.Column(
            pn.Row(self.title_ui, self.progress_ui, align="center"),
            self.info_ui,
            pn.Row(self.btn_prev, self.btn_skip, self.note_input),
            self.plot_container,
            pn.Row(self.btn_accept, self.btn_reject),
            sizing_mode="stretch_both"
        )
        
        if not self.tasks:
            self.info_ui.object = "**All candidates reviewed! No pending tasks.**"
            self.btn_accept.disabled = True
            self.btn_reject.disabled = True
            self.btn_skip.disabled = True
        else:
            self.load_task(self.current_idx)
            
    def _build_tasks(self) -> list[dict]:
        tasks = []
        reviewed_keys = set(self.gold_df['artifact_row_key'])
        
        if self.all_candidates.empty:
            return []
            
        # Group by File (identity_key)
        for ident, group in self.all_candidates.groupby('identity_key'):
            # Only curate tasks that clearly point to an FSA file
            if ".fsa" not in str(ident).lower():
                continue
                
            # Filter markers that haven't been reviewed yet
            unreviewed_group = group[~group['artifact_row_key'].isin(reviewed_keys)]
            if unreviewed_group.empty:
                continue
                
            tasks.append({
                "identity_key": ident,
                "markers": group.sort_values(['channel', 'expected_bp']).to_dict('records')
            })
                
        return tasks
        
    def load_task(self, idx: int):
        self.progress_ui.object = f"**{idx + 1} / {len(self.tasks)}** Files Pending"
        if idx >= len(self.tasks) or idx < 0:
            return
            
        task = self.tasks[idx]
        self.plot_container.clear()
        
        try:
            fsa_path = resolve_fsa_path(self.data_dir, task['identity_key'])
            fsa = load_fsa_full_analysis(fsa_path)
            
            # Map of marker name -> rows in CSV
            marker_groups = {}
            for row in task['markers']:
                m = row['marker_name']
                if m not in marker_groups: marker_groups[m] = []
                marker_groups[m].append(row)
                
            # Create a Plot for each marker (Trace)
            for m_name, candidates in marker_groups.items():
                first = candidates[0]
                channel = first['channel']
                is_ladder = "ladder" in str(m_name).lower() or "ladder" in str(first.get('artifact_row_key', '')).lower()
                
                ebp = first.get('expected_bp', np.nan)
                wbp = first.get('window_bp', 20.0)
                if pd.isna(wbp): wbp = 20.0
                
                # HYPER-AGGRESSIVE CHANNEL DISCOVERY
                def find_trace_data(fsa_obj, target_ch):
                    raw_keys = list(fsa_obj.fsa.keys())
                    
                    def normalize(k):
                        if isinstance(k, bytes): k = k.decode('ascii', errors='ignore')
                        return str(k).strip().upper()

                    # 1. Try normalizing the target
                    target_norm = normalize(target_ch)
                    if not target_norm: return np.array([])
                    
                    # Special mappings for size standards
                    SPECIAL_MAPPINGS = {
                        "LIZ": ["DATA5", "DATA105", "LIZ"],
                        "ROX": ["DATA4", "DATA105", "ROX"],
                    }
                    search_targets = [target_norm]
                    if target_norm in SPECIAL_MAPPINGS:
                        search_targets.extend(SPECIAL_MAPPINGS[target_norm])
                    
                    # 2. Try normalized exact match
                    for t in search_targets:
                        for k in raw_keys:
                            if normalize(k) == t:
                                return np.asarray(fsa_obj.fsa[k])
                            
                    # 3. Try number-based match (e.g. "1" or "DATA1")
                    # Extract last integer from target
                    import re
                    match = re.search(r'(\d+)$', target_norm)
                    if match:
                        num = str(int(match.group(1))) # e.g. "105" -> "105"
                        short_num = num[-1] # e.g. "105" -> "5"
                        for k in raw_keys:
                            k_norm = normalize(k)
                            if k_norm == f"DATA{num}" or k_norm == f"DATA{short_num}":
                                return np.asarray(fsa_obj.fsa[k])
                                
                    return np.array([])

                trace_raw = find_trace_data(fsa, channel)
                
                # Dynamic trace extraction
                raw_df = fsa.sample_data_with_basepairs
                bp_arr = raw_df["basepairs"].to_numpy()
                t_arr = raw_df["time"].astype(int).to_numpy()
                
                # Defensive check for trace presence
                if trace_raw.size == 0:
                    logging.warning(f"Trace {channel} not found in {fsa_path.name}")
                    continue # Skip this plot if trace is missing
                
                # Fast baseline for plot
                from core.analysis import estimate_running_baseline
                try:
                    baseline = estimate_running_baseline(trace_raw, bin_size=200, quantile=0.10)
                    y_proc = (trace_raw.astype(float) - baseline).clip(min=0)
                except Exception:
                    y_proc = trace_raw.astype(float)
                
                # Choose X-axis: BP for sample peaks, Time for ladder peaks
                if is_ladder:
                    # For ladders, we show the FULL RAW TRACE for context
                    y_arr = y_proc
                    x_arr = np.arange(len(y_proc))
                    x_label = "Time (Scan Index)"
                else:
                    # Final index safety for sample peaks (show only where we have basepairs)
                    valid_mask = (t_arr >= 0) & (t_arr < len(y_proc))
                    y_arr = y_proc[t_arr[valid_mask]]
                    x_arr = bp_arr[valid_mask]
                    x_label = "Basepairs (bp)"
                
                # Build Figure
                fig = go.Figure()
                fig.add_trace(go.Scatter(
                    x=x_arr, y=y_arr, mode="lines", name=m_name,
                    line=dict(color="#00FFAA" if not is_ladder else "#FFCC00", width=2),
                    fill="tozeroy", fillcolor="rgba(0, 255, 170, 0.05)"
                ))
                
                # Window
                fig.add_vrect(x0=ebp-wbp, x1=ebp+wbp, fillcolor="rgba(255, 255, 255, 0.05)", line_width=0)
                fig.add_vline(x=ebp, line_dash="dash", line_color="#FF0066", opacity=0.5)
                
                # Plot Candidates
                cand_x = []
                cand_y = []
                for cand in candidates:
                    if cand['ok']:
                        is_sel = cand.get('selected', False) or cand.get('selected_for_fit', False)
                        # Ladder candidates store time in 'found_bp' (due to my normalization)
                        cx = cand['found_bp'] if is_ladder else cand.get('found_bp', np.nan)
                        cy = cand['height']
                        
                        color = "#FF3333" if not is_sel else "#FFFF00"
                        fig.add_trace(go.Scatter(
                            x=[cx], y=[cy],
                            mode="markers", marker=dict(size=14 if is_sel else 8, color=color, symbol="star" if is_sel else "circle"),
                            name=f"{cand['artifact_row_key'].split(':')[-1]}"
                        ))
                        if not np.isnan(cx): cand_x.append(cx)

                # Layout Adjustment
                if is_ladder:
                    # For ladders, show whole trace or at least candidate range
                    if cand_x:
                        x_min, x_max = min(cand_x), max(cand_x)
                        margin = (x_max - x_min) * 0.1
                        ax_range = [x_min - margin, x_max + margin]
                    else:
                        ax_range = None
                else:
                    # Zoom to expected BP window
                    if not pd.isna(ebp):
                        ax_range = [ebp-wbp-25, ebp+wbp+25]
                    else:
                        ax_range = None

                fig.update_layout(
                    title=f"<b>{m_name}</b> ({channel}) | {fsa_path.name}",
                    xaxis_title=x_label,
                    xaxis_range=ax_range,
                    template="plotly_dark", height=350 if is_ladder else 300,
                    margin=dict(l=40, r=40, t=40, b=40),
                    showlegend=False
                )
                self.plot_container.append(pn.pane.Plotly(fig))
                
            if len(self.plot_container) == 0:
                available = []
                for k in fsa.fsa.keys():
                    k_str = k.decode('ascii', errors='ignore') if isinstance(k, bytes) else str(k)
                    if "DATA" in k_str or k_str in ["ROX", "LIZ"]:
                        type_str = "bytes" if isinstance(k, bytes) else "str"
                        available.append(f"{k_str} ({type_str})")
                
                # Show what we were looking for
                expected = list(set([str(row['channel']) for row in task['markers']]))
                        
                self.plot_container.append(pn.pane.Markdown(
                    f"### ⚠️ No traces found for this file.\n"
                    f"**Looking for channels:** `{', '.join(expected)}`\n\n"
                    f"**Actual tags available:**\n`{', '.join(available)}`"
                ))
                
            self.info_ui.object = ""
            
        except Exception as e:
            self.info_ui.object = f"**Error loading File:** {e}"
            
    def save_decision(self, label: str):
        task = self.tasks[self.current_idx]
        new_rows = []
        
        timestamp = datetime.now(timezone.utc).isoformat()
        
        # Save decision for EVERY marker in this task
        for row_data in task['markers']:
            # We ONLY label the "SELECTED" candidate as 'selected_correct'
            # Or if rejecting, we reject the whole set for that row key
            if label == 'selected_correct':
                # Only the one the pipeline picked gets the prize
                if row_data.get('selected', False):
                    l = 'selected_correct'
                else:
                    l = 'noise' # Implicitly noise if not selected
            else:
                l = label
                
            new_rows.append({
                "artifact_table": "pk_candidates" if "pk" in str(row_data['artifact_row_key']) else "ladder_candidates",
                "artifact_row_key": row_data["artifact_row_key"],
                "label": l,
                "label_source": "human_expert",
                "reviewer": "christian",
                "reviewed_at_utc": timestamp,
                "notes": self.note_input.value
            })
        
        self.gold_df = pd.concat([self.gold_df, pd.DataFrame(new_rows)], ignore_index=True)
        self.gold_df.to_csv(self.gold_path, index=False)
        
        self.note_input.value = ""
        self.next_task()
        
    def next_task(self):
        if self.current_idx < len(self.tasks) - 1:
            self.current_idx += 1
            self.load_task(self.current_idx)
        else:
            self.progress_ui.object = "**All files complete!**"
            self.plot_container.clear()
            
    def prev_task(self):
        if self.current_idx > 0:
            self.current_idx -= 1
            self.load_task(self.current_idx)

def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    
    def create_app():
        app = GoldLabelAnnotator(args.candidates_dir, args.data_dir)
        return app.layout
    
    pn.serve(create_app, port=0, show=True, title="Fraggler File Annotator")

if __name__ == "__main__":
    main()
