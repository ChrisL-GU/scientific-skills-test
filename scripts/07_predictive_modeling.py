"""
Predictive Modeling for Disease Classification/Prognosis
Build machine learning models using omics biomarkers
"""

import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (classification_report, confusion_matrix, roc_curve,
                             auc, roc_auc_score, precision_recall_curve)
import warnings
warnings.filterwarnings('ignore')

# Set up paths
OUTPUT_DIR = Path("results/predictions")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def load_omics_expression_data():
    """Load raw omics data for modeling"""

    data_dict = {}

    try:
        # Load RNA-seq counts
        rnaseq = pd.read_csv(Path("data/raw/rnaseq_counts.csv"), index_col=0)
        rnaseq_metadata = pd.read_csv(Path("data/raw/rnaseq_metadata.csv"))
        data_dict['rnaseq'] = (rnaseq.T, rnaseq_metadata)
        print(f"Loaded RNA-seq: {rnaseq.shape}")
    except FileNotFoundError:
        print("RNA-seq data not found")

    try:
        # Load proteomics intensities
        proteomics = pd.read_csv(Path("data/raw/proteomics_intensities.csv"), index_col=0)
        prot_metadata = pd.read_csv(Path("data/raw/proteomics_metadata.csv"))
        data_dict['proteomics'] = (proteomics.T, prot_metadata)
        print(f"Loaded Proteomics: {proteomics.shape}")
    except FileNotFoundError:
        print("Proteomics data not found")

    try:
        # Load metabolomics abundances
        metabolomics = pd.read_csv(Path("data/raw/metabolomics_abundances.csv"), index_col=0)
        met_metadata = pd.read_csv(Path("data/raw/metabolomics_metadata.csv"))
        data_dict['metabolomics'] = (metabolomics.T, met_metadata)
        print(f"Loaded Metabolomics: {metabolomics.shape}")
    except FileNotFoundError:
        print("Metabolomics data not found")

    return data_dict

def build_models(X, y, output_dir):
    """Build and train multiple classification models"""

    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )

    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    print(f"\nTraining set: {X_train.shape}")
    print(f"Test set: {X_test.shape}")

    models = {}
    results = []

    # Model 1: Logistic Regression
    print("\n1. Logistic Regression")
    lr = LogisticRegression(max_iter=1000, random_state=42)
    lr.fit(X_train_scaled, y_train)
    lr_pred = lr.predict(X_test_scaled)
    lr_prob = lr.predict_proba(X_test_scaled)[:, 1]
    lr_score = lr.score(X_test_scaled, y_test)
    lr_auc = roc_auc_score(y_test, lr_prob)
    models['Logistic Regression'] = (lr, lr_prob)
    print(f"Accuracy: {lr_score:.4f}, AUC: {lr_auc:.4f}")
    results.append({
        'model': 'Logistic Regression',
        'accuracy': lr_score,
        'auc': lr_auc
    })

    # Model 2: Random Forest
    print("\n2. Random Forest")
    rf = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
    rf.fit(X_train, y_train)
    rf_pred = rf.predict(X_test)
    rf_prob = rf.predict_proba(X_test)[:, 1]
    rf_score = rf.score(X_test, y_test)
    rf_auc = roc_auc_score(y_test, rf_prob)
    models['Random Forest'] = (rf, rf_prob)
    print(f"Accuracy: {rf_score:.4f}, AUC: {rf_auc:.4f}")
    results.append({
        'model': 'Random Forest',
        'accuracy': rf_score,
        'auc': rf_auc
    })

    # Model 3: Gradient Boosting
    print("\n3. Gradient Boosting")
    gb = GradientBoostingClassifier(n_estimators=100, random_state=42)
    gb.fit(X_train, y_train)
    gb_pred = gb.predict(X_test)
    gb_prob = gb.predict_proba(X_test)[:, 1]
    gb_score = gb.score(X_test, y_test)
    gb_auc = roc_auc_score(y_test, gb_prob)
    models['Gradient Boosting'] = (gb, gb_prob)
    print(f"Accuracy: {gb_score:.4f}, AUC: {gb_auc:.4f}")
    results.append({
        'model': 'Gradient Boosting',
        'accuracy': gb_score,
        'auc': gb_auc
    })

    results_df = pd.DataFrame(results)
    results_df.to_csv(output_dir / "model_performance.csv", index=False)

    return models, X_test_scaled, y_test, results_df

def plot_roc_curves(models, X_test, y_test, output_dir):
    """Plot ROC curves for all models"""

    fig, ax = plt.subplots(figsize=(10, 8))

    for model_name, (model, probs) in models.items():
        if isinstance(model, LogisticRegression):
            y_pred_prob = model.predict_proba(X_test)[:, 1]
        else:
            y_pred_prob = probs

        fpr, tpr, _ = roc_curve(y_test, y_pred_prob)
        roc_auc = auc(fpr, tpr)
        ax.plot(fpr, tpr, linewidth=2, label=f'{model_name} (AUC = {roc_auc:.3f})')

    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, label='Random Classifier')
    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    ax.set_title('ROC Curves - Infection Status Classification', fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'roc_curves.png', dpi=300)
    plt.close()

def plot_confusion_matrices(models, X_test, y_test, output_dir):
    """Plot confusion matrices for all models"""

    n_models = len(models)
    fig, axes = plt.subplots(1, n_models, figsize=(5*n_models, 4))
    if n_models == 1:
        axes = [axes]

    for ax, (model_name, (model, _)) in zip(axes, models.items()):
        if isinstance(model, LogisticRegression):
            y_pred = model.predict(X_test)
        else:
            y_pred = model.predict(X_test)

        cm = confusion_matrix(y_test, y_pred)
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax, cbar=False)
        ax.set_title(f'{model_name}', fontsize=12)
        ax.set_xlabel('Predicted')
        ax.set_ylabel('True')

    plt.tight_layout()
    plt.savefig(output_dir / 'confusion_matrices.png', dpi=300)
    plt.close()

def feature_importance_analysis(models, output_dir):
    """Analyze and plot feature importance from tree-based models"""

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Random Forest feature importance
    rf_model, _ = models.get('Random Forest', (None, None))
    if rf_model is not None:
        importances = rf_model.feature_importances_
        indices = np.argsort(importances)[-15:]
        ax = axes[0]
        ax.barh(range(len(indices)), importances[indices])
        ax.set_yticks(range(len(indices)))
        ax.set_yticklabels([f'Feature {i}' for i in indices])
        ax.set_xlabel('Importance')
        ax.set_title('Random Forest - Top 15 Features')

    # Gradient Boosting feature importance
    gb_model, _ = models.get('Gradient Boosting', (None, None))
    if gb_model is not None:
        importances = gb_model.feature_importances_
        indices = np.argsort(importances)[-15:]
        ax = axes[1]
        ax.barh(range(len(indices)), importances[indices])
        ax.set_yticks(range(len(indices)))
        ax.set_yticklabels([f'Feature {i}' for i in indices])
        ax.set_xlabel('Importance')
        ax.set_title('Gradient Boosting - Top 15 Features')

    plt.tight_layout()
    plt.savefig(output_dir / 'feature_importance.png', dpi=300)
    plt.close()

if __name__ == "__main__":
    print("=" * 60)
    print("Predictive Modeling - Infection Status Classification")
    print("=" * 60)

    # Load data
    print("\nLoading omics data...")
    data_dict = load_omics_expression_data()

    if len(data_dict) == 0:
        print("No data loaded. Generating example data...")
        # Generate synthetic data for demonstration
        np.random.seed(42)
        n_samples = 24
        n_features = 500

        # Features (combined omics)
        X = np.random.randn(n_samples, n_features)

        # Increase feature values for infected samples
        X[:n_samples//2] += np.random.randn(n_samples//2, n_features) * 2

        # Labels
        y = np.array([1]*12 + [0]*12)

        print(f"Generated synthetic data: {X.shape}, Labels: {y.shape}")

    else:
        # Combine omics data
        all_data = []
        for omics_type, (data, metadata) in data_dict.items():
            all_data.append(data)

        X = pd.concat(all_data, axis=1).fillna(0).values
        y = data_dict[list(data_dict.keys())[0]][1]['condition'].map({'Infected': 1, 'Control': 0}).values

    # Build models
    print("\nBuilding predictive models...")
    models, X_test, y_test, results_df = build_models(X, y, OUTPUT_DIR)

    # Plot ROC curves
    print("\nPlotting ROC curves...")
    plot_roc_curves(models, X_test, y_test, OUTPUT_DIR)

    # Plot confusion matrices
    print("Plotting confusion matrices...")
    plot_confusion_matrices(models, X_test, y_test, OUTPUT_DIR)

    # Feature importance
    print("Analyzing feature importance...")
    feature_importance_analysis(models, OUTPUT_DIR)

    print("\n" + "=" * 60)
    print("Predictive modeling complete!")
    print(f"Results saved to {OUTPUT_DIR}")
    print("=" * 60)

    print("\nModel Comparison:")
    print(results_df.to_string(index=False))
