"""
Autogluon Model Prediction 
"""

import pandas as pd
from autogluon.tabular import TabularDataset, TabularPredictor
import itertools
import logging

# Configuration
MODE = "medium_quality"
METRICS_LIST = ["ZIP"]
FP_TYPES = ["MACCS", "PubChem", "Substructure"]
MODELS = [
    "WeightedEnsemble_L2", 
    "NeuralNetFastAI_BAG_L1", 
    "NeuralNetTorch_BAG_L1",
    "LightGBM_BAG_L1", 
    "XGBoost_BAG_L1", 
    "LightGBMLarge_BAG_L1", 
    "RandomForestMSE_BAG_L1",
    "LightGBMXT_BAG_L1", 
    "CatBoost_BAG_L1", 
    "ExtraTreesMSE_BAG_L1"
]
TCM_FEAT_PATH = "./ML/tcm/tcm_dat.csv"
OUTPUT_FILE = "./ML/tcm/tcm_pred.csv"

# Logging setup
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_tcm_data(file_path):
    """Load TCM feature data from the specified file path."""
    logger.info(f"Loading TCM feature data from {file_path}")
    return pd.read_csv(file_path)

def model_predict(tcm_feat, metrics, fp_type, model):
    """
    Perform predictions using a pre-trained Autogluon model.

    Args:
        tcm_feat (pd.DataFrame): TCM feature data.
        metrics (str): Metrics used for training the model.
        fp_type (str): Fingerprint type used for training the model.
        model (str): Model to use for prediction.

    Returns:
        pd.DataFrame: Predictions with additional metadata.
    """
    logger.info(f"Predicting with {metrics}, {fp_type}, {model}")
    
    try:
        # Select relevant columns based on fingerprint type
        selected_columns = [col for col in tcm_feat.columns if col.startswith(fp_type)]
        data_Tab = TabularDataset(tcm_feat[selected_columns])

        # Load the pre-trained predictor
        save_path = f"./ML/model/{metrics}_{MODE}/{fp_type}"
        predictor = TabularPredictor.load(save_path)

        # Perform predictions
        preds_df = pd.DataFrame(predictor.predict(data_Tab, model=model))
        preds_df.columns = ['Prob']
        preds_df['metrics'] = metrics
        preds_df['fp_type'] = fp_type
        preds_df['model'] = model
        preds_df = pd.concat([tcm_feat['Id'], preds_df], axis=1)

        return preds_df
    except Exception as e:
        logger.error(f"Error during prediction for {metrics}, {fp_type}, {model}: {e}")
        return None

def main():
    # Load TCM feature data
    tcm_feat = load_tcm_data(TCM_FEAT_PATH)

    # Perform predictions for all combinations of metrics, fingerprint types, and models
    metrics_df_merge = []
    for metrics, fp_type, model in itertools.product(METRICS_LIST, FP_TYPES, MODELS):
        logger.info(f"Processing combination: {metrics}, {fp_type}, {model}")
        metrics_df = model_predict(tcm_feat, metrics, fp_type, model)
        if metrics_df is not None:
            metrics_df_merge.append(metrics_df)

    # Merge all predictions and save to file
    if metrics_df_merge:
        merged_df = pd.concat(metrics_df_merge, axis=0)
        merged_df.to_csv(OUTPUT_FILE, index=False)
        logger.info(f"Predictions saved to {OUTPUT_FILE}")
    else:
        logger.warning("No predictions were generated. Check logs for errors.")

if __name__ == "__main__":
    main()