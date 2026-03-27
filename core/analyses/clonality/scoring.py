from __future__ import annotations
import numpy as np
import pandas as pd
from typing import Any

class MLScorer:
    """
    Offline scorer helper for clonality candidate research.

    This module is intentionally not part of the live PK selection path.
    It can be used to compare offline candidate rankings against
    human-reviewed labels while the deterministic runtime path stays primary.
    """
    
    # These were generated from eval_output_human_curated_20260325/pk_weights.json
    INTERCEPT = 0.251561368795845
    COEFFICIENTS = {
        "expected_bp": 0.11430551031876494,
        "window_bp": 0.5743183426461221,
        "search_window_bp": -1.3898395605853142,
        "ok": 0.16903158323156806,
        "found_bp": -0.11172504170071756,
        "delta_bp": -0.22603055707642847,
        "height": 0.00017728339358381073,
        "area": -5.676824438860409e-07,
        "selection_score": -0.33827296851052147,
        "fallback_from_window_bp": 1.0400823269182697,
        "selected": 1.3936213181151473,
        "marker_name_DHJH_D_PK_DATA2_139": 0.11063842564177144,
        "marker_name_DHJH_E_PK_DATA1_109": -0.09113063582030805,
        "marker_name_FR1_PK_DATA1_325": -0.07578435806769084,
        "marker_name_FR2_PK_DATA1_260": -0.043912618555202676,
        "marker_name_FR3_PK_DATA2_145": 0.041420315583865074,
        "marker_name_IGK_PK_DATA1_279": -0.24965563061890086,
        "marker_name_IGK_PK_DATA2_150": -0.3463896441433765,
        "marker_name_KDE_PK_DATA3_287": -0.2270875018326563,
        "marker_name_KDE_PK_DATA3_377": -0.2958401423125043,
        "marker_name_LIZ_Ladder_200": -0.6963255782725891,
        "marker_name_ROX_Ladder_150": -0.04021515692454417,
        "marker_name_ROX_Ladder_280": 1.3684519966040167,
        "marker_name_ROX_Ladder_300": -0.3787542888843909,
        "marker_name_TCRbA_PK_DATA2_265": -0.07472673022712192,
        "marker_name_TCRbB_PK_DATA1_254": 0.21597817602957126,
        "marker_name_TCRbC_PK_DATA2_311": -0.011492133410984946,
        "marker_name_TCRgA_PK1_DATA1_249": 0.048143991609845616,
        "marker_name_TCRgA_PK1_DATA2_212": 0.050616828606292956,
        "marker_name_TCRgA_PK1_LIZ_Ladder_200": -0.06760613995628287,
        "marker_name_TCRgA_PK2_DATA2_163": 0.02333081401649318,
        "marker_name_TCRgA_PK2_LIZ_Ladder_200": -0.035406338215881544,
        "marker_name_TCRgB_PK1_DATA2_115": 0.27418663218214073,
        "marker_name_TCRgB_PK1_LIZ_Ladder_200": 0.21325360628191642,
        "marker_name_TCRgB_PK2_DATA2_178": 0.19387296678105434,
        "marker_name_TCRgB_PK2_LIZ_Ladder_200": 0.26346472208237237,
        "kind_ladder": 0.6268628166232185,
        "kind_sample": -0.45783123127120934,
        "channel_DATA1": -0.1963610728280221,
        "channel_DATA105": -0.3226197273453048,
        "channel_DATA2": 0.2614574711705638,
        "channel_DATA3": -0.5229276425069711,
        "channel_DATA4": 0.9494825304215534,
        "search_mode_fallback": -0.3782037174958675,
        "search_mode_primary": 0.5472353466360266,
        "assay_DHJH_D": -0.26811586580190505,
        "assay_DHJH_E": -0.13134579583150924,
        "assay_FR1": 0.06389136404968589,
        "assay_FR2": 0.09519378819402509,
        "assay_FR3": 0.1810365653346401,
        "assay_IGK": -0.6527429794743941,
        "assay_KDE": -1.1625554982999844,
        "assay_TCRbA": 0.24083326636693758,
        "assay_TCRbB": 0.5193010562440299,
        "assay_TCRbC": 0.3196786107110753,
        "assay_TCRgA": 0.019079158050843214,
        "assay_TCRgB": 0.944777918138888,
        "control_PK": -0.7948254963065366,
        "control_PK1": 0.5185949126551729,
        "control_PK2": 0.44526216877037,
        "sample_kind_control": 0.16903158323156484
    }

    def score_candidate(self, features: dict[str, Any]) -> float:
        """
        Calculates the probability [0, 1] that this candidate is the 'correct' peak
        as determined by the human clinical expert.
        """
        score = self.INTERCEPT
        
        # 1. Continuous Features
        num_cols = ["expected_bp", "window_bp", "search_window_bp", "ok", 
                    "found_bp", "delta_bp", "height", "area", 
                    "selection_score", "fallback_from_window_bp", "selected"]
        
        for col in num_cols:
            val = features.get(col, 0.0)
            if val is None or pd.isna(val):
                val = 0.0
            score += float(val) * self.COEFFICIENTS.get(col, 0.0)
            
        # 2. Categorical Features (One-Hot)
        cat_cols = ["marker_name", "kind", "channel", "search_mode", "assay", "control", "sample_kind"]
        for col in cat_cols:
            val = str(features.get(col, "")).strip()
            if not val:
                continue
            one_hot_key = f"{col}_{val}"
            score += self.COEFFICIENTS.get(one_hot_key, 0.0)
            
        # Sigmoid function
        return 1.0 / (1.0 + np.exp(-score))

    def select_best(self, candidates: list[dict[str, Any]]) -> dict[str, Any] | None:
        """
        Given a list of candidates for a marker, returns the one with the highest ML score.
        """
        if not candidates:
            return None
        
        best_cand = None
        max_prob = -1.0
        
        for cand in candidates:
            # We must make sure the candidate has 'delta_bp' before scoring
            # if it doesn't already (delta_bp = found_bp - expected_bp)
            if 'delta_bp' not in cand and 'found_bp' in cand and 'expected_bp' in cand:
                cand = cand.copy()
                cand['delta_bp'] = abs(float(cand['found_bp']) - float(cand['expected_bp']))
            
            prob = self.score_candidate(cand)
            cand['ml_score'] = prob
            if prob > max_prob:
                max_prob = prob
                best_cand = cand
                
        return best_cand

# Singleton instance for production use
CLONALITY_SCORER = MLScorer()
