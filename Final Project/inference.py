import pandas as pd
import joblib
import os
import sys
import argparse

class HDClassifier:
    def __init__(self, model_dir='hd_models'):
        """
        Initialize the classifier by loading the saved artifacts for BOTH models.
        """
        self.model_dir = model_dir
        
        # Define paths to artifacts
        rf_path = os.path.join(model_dir, 'rf_model.pkl')
        svm_path = os.path.join(model_dir, 'svm_model.pkl')
        scaler_path = os.path.join(model_dir, 'scaler.pkl')
        features_path = os.path.join(model_dir, 'feature_names.pkl')
        le_path = os.path.join(model_dir, 'label_encoder.pkl')

        # Load artifacts with error handling
        try:
            self.rf_model = joblib.load(rf_path)
            self.svm_model = joblib.load(svm_path)
            self.scaler = joblib.load(scaler_path)
            self.feature_names = joblib.load(features_path)
            self.label_encoder = joblib.load(le_path)
            print(f"Successfully loaded RF and SVM models from '{model_dir}'")
        except FileNotFoundError as e:
            print(f"Error: Could not load model artifacts. {e}")
            print(f"Ensure the folder '{model_dir}' exists and contains all .pkl files.")
            sys.exit(1)

    def predict(self, input_data):
        """
        Predict HD status using a Weighted Ensemble (SVM dominant).
        Logic: Final_Prob = (0.7 * SVM_Prob) + (0.3 * RF_Prob)
        """
        # 1. Load Data
        if isinstance(input_data, str):
            if not os.path.exists(input_data):
                print(f"Error: Input file '{input_data}' not found.")
                sys.exit(1)
            df = pd.read_csv(input_data)
        else:
            df = input_data.copy()

        # 2. Validate Features
        missing_cols = [col for col in self.feature_names if col not in df.columns]
        if missing_cols:
            print(f"Error: Input data is missing {len(missing_cols)} required genes.")
            sys.exit(1)

        # 3. Align Data
        X = df[self.feature_names]

        print("Running Weighted Ensemble Inference (SVM: 0.7, RF: 0.3)...")

        # --- Random Forest ---
        rf_probs = self.rf_model.predict_proba(X)[:, 1]
        rf_preds = self.rf_model.predict(X)

        # --- SVM ---
        X_scaled = self.scaler.transform(X)
        svm_probs = self.svm_model.predict_proba(X_scaled)[:, 1]
        svm_preds = self.svm_model.predict(X_scaled)

        # --- Weighted Ensemble Logic ---
        # SVM gets 0.7 weight because it has higher accuracy (0.95 vs 0.85)
        ensemble_probs = (0.7 * svm_probs) + (0.3 * rf_probs)
        
        # Final Thresholding (0.5)
        ensemble_pred_idx = (ensemble_probs >= 0.5).astype(int)

        # 5. Decode Results
        final_labels = self.label_encoder.inverse_transform(ensemble_pred_idx)
        rf_labels = self.label_encoder.inverse_transform(rf_preds)
        svm_labels = self.label_encoder.inverse_transform(svm_preds)

        # 6. Results DataFrame
        results = pd.DataFrame({
            'Sample_Index': df.index,
            'Final_Prediction': final_labels,
            'Biomarker_Confidence_Score': ensemble_probs,
            'SVM_Prediction': svm_labels,
            'RF_Prediction': rf_labels
        })
        
        return results

def main():
    parser = argparse.ArgumentParser(description="HD Biomarker Classifier (Weighted Ensemble)")
    parser.add_argument("input_file", help="Path to input CSV.")
    parser.add_argument("--model_dir", default="hd_models", help="Model directory.")
    parser.add_argument("--output", default="weighted_predictions.csv", help="Output CSV path.")

    args = parser.parse_args()

    classifier = HDClassifier(model_dir=args.model_dir)
    results = classifier.predict(args.input_file)

    results.to_csv(args.output, index=False)
    print(f"\nPredictions saved to '{args.output}'")
    print("\n--- Result Preview ---")
    print(results[['Sample_Index', 'Final_Prediction', 'Biomarker_Confidence_Score']].head())

if __name__ == "__main__":
    main()