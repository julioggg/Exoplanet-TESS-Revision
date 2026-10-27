"""
============================================================
GENERATE PUBLICATION FIGURES
============================================================

Generates all publication-quality figures used in the
scientific manuscript.

Project:
Multi-objective Optimization for Exoplanet Detection
using TESS Light Curves

Author:
Bryan Cabrera

============================================================
"""
from sklearn.metrics import (
    confusion_matrix,
    ConfusionMatrixDisplay
)
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    roc_curve,
    roc_auc_score
)
import numpy as np
import lightkurve as lk
from pathlib import Path
from astropy.timeseries import BoxLeastSquares
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_config import (
    setup_plot,
    finish_plot,
    get_color
)

# ============================================================
# PATHS
# ============================================================

DATASET = Path("data/dataset/final_dataset.csv")

POSITIVE_BLS = Path(
    "data/features/positives/positive_bls_features.csv"
)

FEATURE_IMPORTANCE = Path(
    "results/feature_importance.csv"
)

PREDICTIONS = Path(
    "results/test_predictions.csv"
)

OUTPUT = Path("results/figures")

OUTPUT.mkdir(
    parents=True,
    exist_ok=True
)

# ============================================================
# LOAD DATA
# ============================================================

print("=" * 60)
print("LOADING DATA")
print("=" * 60)

dataset = pd.read_csv(DATASET)

positive_bls = pd.read_csv(POSITIVE_BLS)

feature_importance = pd.read_csv(FEATURE_IMPORTANCE)

predictions = pd.read_csv(PREDICTIONS)

print()

print(f"Dataset shape              : {dataset.shape}")

print(f"Positive BLS shape         : {positive_bls.shape}")

print(f"Feature importance entries : {len(feature_importance)}")

print(f"Predictions                : {len(predictions)}")

print()

# ============================================================
# FIGURE PLACEHOLDERS
# ============================================================
from matplotlib.patches import FancyBboxPatch

def figure01_pipeline():

    print("[Figure 01] Methodology Pipeline")

    plt.figure(figsize=(16,4))

    ax = plt.gca()
    ax.set_xlim(0,16)
    ax.set_ylim(0,4)
    ax.axis("off")

    labels = [

        "TESS\nLight Curves",

        "Preprocessing\nCleaning\nNormalization",

        "Feature\nExtraction\n(BLS + Statistics)",

        "Multi-objective\nOptimization",

        "Random Forest\nClassification",

        "SHAP\nExplainability"

    ]

    x = [0.5,3.2,6.1,9.0,11.9,14.8]

    colors = [

        "#4C78A8",

        "#72B7B2",

        "#54A24B",

        "#F58518",

        "#E45756",

        "#B279A2"

    ]

    for xi,label,c in zip(x,labels,colors):

        box = FancyBboxPatch(

            (xi,1.3),

            1.6,

            1.2,

            boxstyle="round,pad=0.05",

            linewidth=2,

            facecolor=c,

            edgecolor="black"

        )

        ax.add_patch(box)

        ax.text(

            xi+0.8,

            1.9,

            label,

            ha="center",

            va="center",

            fontsize=11,

            color="white",

            fontweight="bold"

        )

    for i in range(len(x)-1):

        ax.annotate(

            "",

            xy=(x[i+1],1.9),

            xytext=(x[i]+1.6,1.9),

            arrowprops=dict(

                arrowstyle="->",

                lw=2.5

            )

        )

    plt.title(

        "Overall Methodology Pipeline",

        fontsize=18,

        weight="bold"

    )

    plt.savefig(

        OUTPUT/"Figure01_Methodology_Pipeline.png",

        dpi=600,

        bbox_inches="tight"

    )

    plt.savefig(

        OUTPUT/"Figure01_Methodology_Pipeline.pdf",

        bbox_inches="tight"

    )

    plt.close()

    print("Saved: Figure01_Methodology_Pipeline")
# ============================================================
# FIGURE 02
# Dataset Distribution
# ============================================================

def figure02_dataset_distribution():

    print("[Figure 02] Dataset Distribution")

    setup_plot()

    counts = (
        dataset["label"]
        .value_counts()
        .sort_index()
    )

    negatives = counts[0]
    positives = counts[1]

    total = negatives + positives

    percentages = [
        negatives / total * 100,
        positives / total * 100
    ]

    labels = [
        "Negative",
        "Positive"
    ]

    values = [
        negatives,
        positives
    ]

    colors = [
        get_color("orange"),
        get_color("blue")
    ]

    fig, ax = plt.subplots()

    bars = ax.bar(

        labels,

        values,

        color=colors,

        edgecolor="black",

        linewidth=1.2,

        width=0.72

    )

    ax.set_ylabel(
        "Number of Samples",
        fontweight="bold"
    )

    ax.set_title(
        "Dataset Class Distribution",
        fontweight="bold",
        pad=15
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.grid(
        axis="y",
        linestyle="--",
        alpha=0.30
    )

    ax.grid(False, axis="x")

    for i, bar in enumerate(bars):

        height = bar.get_height()

        ax.text(

            bar.get_x() + bar.get_width()/2,

            height + 4,

            f"{height}\n({percentages[i]:.1f}%)",

            ha="center",

            va="bottom",

            fontsize=12,

            fontweight="bold"

        )

    png_file = OUTPUT / "Figure02_Dataset_Distribution.png"

    pdf_file = OUTPUT / "Figure02_Dataset_Distribution.pdf"

    plt.savefig(
        png_file,
        dpi=600
    )

    plt.savefig(
        pdf_file
    )

    plt.close()

    print("   PNG :", png_file)

    print("   PDF :", pdf_file)

    print()

# ============================================================
# FIGURE 03
# Orbital Period Distribution
# ============================================================

def figure03_period_distribution():

    print("[Figure 03] Orbital Period Distribution")

    setup_plot()

    periods = positive_bls["best_period"].dropna()

    mean_period = periods.mean()
    median_period = periods.median()

    fig, ax = plt.subplots()

    ax.hist(
        periods,
        bins=25,
        color=get_color("blue"),
        edgecolor="black",
        alpha=0.80
    )

    ax.axvline(
        mean_period,
        color=get_color("orange"),
        linestyle="--",
        linewidth=2,
        label=f"Mean = {mean_period:.2f} d"
    )

    ax.axvline(
        median_period,
        color=get_color("green"),
        linestyle="-.",
        linewidth=2,
        label=f"Median = {median_period:.2f} d"
    )

    ax.set_xlabel(
        "Orbital Period (days)",
        fontweight="bold"
    )

    ax.set_ylabel(
        "Number of Exoplanets",
        fontweight="bold"
    )

    ax.set_title(
        "Distribution of Detected Orbital Periods",
        fontweight="bold",
        pad=15
    )

    ax.legend()

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.grid(
        axis="y",
        linestyle="--",
        alpha=0.30
    )

    finish_plot()

    png_file = OUTPUT / "Figure03_Orbital_Period_Distribution.png"

    pdf_file = OUTPUT / "Figure03_Orbital_Period_Distribution.pdf"

    plt.savefig(
        png_file,
        dpi=600
    )

    plt.savefig(
        pdf_file
    )

    plt.close()

    print("   PNG :", png_file)
    print("   PDF :", pdf_file)
    print()

# ============================================================
# FIGURE 04
# Transit Depth Distribution
# ============================================================

def figure04_transit_depth_distribution():

    print("[Figure 04] Transit Depth Distribution")

    setup_plot()

    depths = positive_bls["transit_depth"].dropna()

    mean_depth = depths.mean()
    median_depth = depths.median()

    fig, ax = plt.subplots()

    ax.hist(
        depths,
        bins=25,
        color=get_color("purple"),
        edgecolor="black",
        alpha=0.80
    )

    ax.axvline(
        mean_depth,
        color=get_color("orange"),
        linestyle="--",
        linewidth=2,
        label=f"Mean = {mean_depth:.4f}"
    )

    ax.axvline(
        median_depth,
        color=get_color("green"),
        linestyle="-.",
        linewidth=2,
        label=f"Median = {median_depth:.4f}"
    )

    ax.set_xlabel(
        "Transit Depth",
        fontweight="bold"
    )

    ax.set_ylabel(
        "Number of Exoplanets",
        fontweight="bold"
    )

    ax.set_title(
        "Distribution of Transit Depth",
        fontweight="bold",
        pad=15
    )

    ax.legend()

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.grid(
        axis="y",
        linestyle="--",
        alpha=0.30
    )

    finish_plot()

    png_file = OUTPUT / "Figure04_Transit_Depth_Distribution.png"
    pdf_file = OUTPUT / "Figure04_Transit_Depth_Distribution.pdf"

    plt.savefig(
        png_file,
        dpi=600
    )

    plt.savefig(
        pdf_file
    )

    plt.close()

    print("   PNG :", png_file)
    print("   PDF :", pdf_file)
    print()

# ============================================================
# FIGURE 05
# Original vs Preprocessed TESS Light Curve
# ============================================================

def figure05_lightcurve():

    print("[Figure 05] Original vs Preprocessed Light Curve")

    setup_plot()

    fits_file = Path(
        "data/raw/positives/TIC_101145648.fits"
    )

    lc = lk.read(fits_file)

    # -----------------------------
    # Original data
    # -----------------------------

    time = np.array(lc.time.value)

    flux = np.array(lc.flux.value)

    mask = np.isfinite(time) & np.isfinite(flux)

    time_original = time[mask]

    flux_original = flux[mask]

    # -------------------------------------------------
    # Detect temporal gaps
    # -------------------------------------------------

    dt = np.diff(time_original)

    gap_threshold = 0.05      # days

    gap_index = np.where(dt > gap_threshold)[0]

    # -----------------------------
    # Preprocessed data
    # -----------------------------

    median_flux = np.median(flux_original)

    flux_normalized = flux_original / median_flux

    # -----------------------------
    # Figure
    # -----------------------------

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(10,8),
        sharex=True
    )

    # ============================================
    # Original Light Curve
    # ============================================
    # Plot each continuous segment separately

    start = 0

    for idx in gap_index:

        axes[0].plot(

            time_original[start:idx+1],

            flux_original[start:idx+1],

            color=get_color("blue"),

            linewidth=0.5

        )

        start = idx + 1

    axes[0].plot(

        time_original[start:],

        flux_original[start:],

        color=get_color("blue"),

        linewidth=0.5

    )
    

    axes[0].set_ylabel(
        "Flux (e-/s)",
        fontweight="bold"
    )

    axes[0].set_title(
        "Original TESS Light Curve",
        fontweight="bold"
    )

    axes[0].spines["top"].set_visible(False)
    axes[0].spines["right"].set_visible(False)

    axes[0].grid(
        alpha=0.30,
        linestyle="--"
    )

    # ============================================
    # Preprocessed Light Curve
    # ============================================

    start = 0

    for idx in gap_index:

        axes[1].plot(

            time_original[start:idx+1],

            flux_normalized[start:idx+1],

            color=get_color("green"),

            linewidth=0.5

    )

    start = idx + 1

    axes[1].plot(

        time_original[start:],

        flux_normalized[start:],

        color=get_color("green"),

        linewidth=0.5

    )

    axes[1].set_xlabel(
        "Time (BTJD)",
        fontweight="bold"
    )

    axes[1].set_ylabel(
        "Normalized Flux",
        fontweight="bold"
    )
    from matplotlib.ticker import ScalarFormatter

    formatter = ScalarFormatter(useOffset=False)
    formatter.set_scientific(False)

    axes[1].yaxis.set_major_formatter(formatter)

    axes[1].set_title(
        "Preprocessed Light Curve",
        fontweight="bold"
    )

    axes[1].spines["top"].set_visible(False)
    axes[1].spines["right"].set_visible(False)

    axes[1].grid(
        alpha=0.30,
        linestyle="--"
    )

    plt.tight_layout()

    png_file = OUTPUT / "Figure05_Light_Curve.png"

    pdf_file = OUTPUT / "Figure05_Light_Curve.pdf"

    plt.savefig(
        png_file,
        dpi=600
    )

    plt.savefig(
        pdf_file
    )

    plt.close()

    print("   PNG :", png_file)

    print("   PDF :", pdf_file)

    print()

# ==========================================================
# Load an example TESS light curve
# ==========================================================

def load_example_lightcurve():

    """
    Load a representative positive TESS light curve.

    Returns
    -------
    lc : LightCurve
        Preprocessed normalized light curve.
    """

    import lightkurve as lk

    fits_path = (
        "data/raw/positives/"
        "TIC_101145648.fits"
    )

    print()
    print("Loading representative TESS target...")
    print(fits_path)

    lc = lk.read(fits_path)

    lc = (
        lc
        .remove_nans()
        .normalize()
        .flatten(window_length=401)
    )

    return lc

# ==========================================================
# Figure 06
# BLS Periodogram
# ==========================================================

def figure06_bls_periodogram():
    
    print("[Figure 06] BLS Periodogram")

    # ------------------------------------------------------
    # Load the same processed light curve used during
    # feature extraction
    # ------------------------------------------------------
    # ------------------------------------------------------
    # Load representative processed light curve
    # ------------------------------------------------------

    tic = "TIC_101145648"

    print(f"Selected target: {tic}")

    df = pd.read_csv(
        f"data/processed/positives/{tic}_processed.csv"
    )

    time = df["TIME"].values
    flux = df["FLUX_NORM"].values
    print()

    print(df.describe())

    # ------------------------------------------------------
    # Search periods
    # ------------------------------------------------------

    periods = np.linspace(
        0.5,
        25,
        5000
    )

    bls = BoxLeastSquares(
        time,
        flux
    )

    durations = np.array([
        0.05,
        0.10,
        0.15,
        0.20,
        0.25
    ])

    result = bls.power(
        periods,
        durations
    )
    power = result.power
    print()

    print("Maximum power:", np.max(result.power))

    print("Minimum power:", np.min(result.power))

    print("Mean power:", np.mean(result.power))

    best_index = np.argmax(power)

    best_period = result.period[best_index]

    best_power = power[best_index]

    print(f"Best Period : {best_period:.5f} days")
    print(f"BLS Power : {best_power:.6e}")

    # ------------------------------------------------------
    # Plot
    # ------------------------------------------------------

    setup_plot()

    plt.figure(figsize=(10,6))

    plt.plot(
        result.period,
        power,
        color=get_color("blue"),
        linewidth=1.4
    )

    plt.axvline(
        best_period,
        color=get_color("orange"),
        linestyle="--",
        linewidth=2.2,
        label=f"Best Period = {best_period:.3f} d"
    )

    plt.scatter(
        best_period,
        best_power,
        s=80,
        color="red",
        zorder=5
    )

    plt.xlabel("Orbital Period (days)")

    plt.ylabel("BLS Power")

    plt.title(
        "Box Least Squares Periodogram"
    )

    plt.legend()

    finish_plot()

    plt.savefig(
    OUTPUT / "Figure06_BLS_Periodogram.png",
    dpi=600
    )

    plt.savefig(
    OUTPUT / "Figure06_BLS_Periodogram.pdf"
    )

    plt.close()

    print("Saved: Figure06_BLS_Periodogram")

def figure07_folded_transit():

    print("[Figure 07] Folded Transit Curve")

    # ------------------------------------------------------
    # Load representative processed light curve
    # ------------------------------------------------------

    # ------------------------------------------------------
    # Automatically select the best transit candidate
    # ------------------------------------------------------

    features = pd.read_csv(
    "data/features/positives/positive_bls_features.csv"
    )

    # Highest SNR first
    features = features.sort_values(
        "snr_estimate",
        ascending=False
    )

    selected = None

    for _, row in features.iterrows():

        tic = row["tic_id"]

        candidate = Path(
            f"data/processed/positives/{tic}_processed.csv"
        )

        if candidate.exists():

            selected = row
            break

    if selected is None:
        raise FileNotFoundError(
            "No processed candidate was found."
        )

    tic = selected["tic_id"]

    best_period = selected["best_period"]

    print(f"Selected target : {tic}")
    print(f"Best period     : {best_period:.5f} d")
    print(f"SNR estimate    : {selected['snr_estimate']:.2f}")

    df = pd.read_csv(
        f"data/processed/positives/{tic}_processed.csv"
    )

    time = df["TIME"].values
    flux = df["FLUX_NORM"].values

    # ------------------------------------------------------
    # Fold light curve
    # ------------------------------------------------------

    phase = (
        (time / best_period + 0.5) % 1
    ) - 0.5

    order = np.argsort(phase)

    phase = phase[order]
    flux = flux[order]

    # ------------------------------------------------------
    # Bin folded light curve
    # ------------------------------------------------------

    nbins = 120

    bins = np.linspace(
        -0.5,
        0.5,
        nbins + 1
    )

    centers = 0.5 * (bins[:-1] + bins[1:])
    
    digitized = np.digitize(
    phase,
    bins
    )

    binned_flux = []

    for i in range(1, len(bins)):

        mask = (digitized == i)

        if np.sum(mask) > 0:

            binned_flux.append(
                np.mean(flux[mask])
            )

        else:

            binned_flux.append(
                np.nan
            )
    binned_flux = np.array(binned_flux)
    # ------------------------------------------------------
    # Plot
    # ------------------------------------------------------

    setup_plot()

    plt.figure(figsize=(10,6))

    plt.scatter(
        phase,
        flux,
        s=5,
        alpha=0.18,
        color=get_color("blue"),
        label="Observations"
    )

    plt.plot(
        centers,
        binned_flux,
        color=get_color("orange"),
        linewidth=3,
        label="Binned curve"
    )

    plt.axvline(
        0,
        linestyle="--",
        color="gray",
        linewidth=2
    )

    plt.xlabel("Orbital Phase")

    plt.ylabel("Normalized Flux")

    plt.title(
        "Phase-folded Transit Light Curve"
    )

    plt.legend()

    finish_plot()

    plt.savefig(
        OUTPUT / "Figure07_Folded_Transit.png",
        dpi=600
    )

    plt.savefig(
        OUTPUT / "Figure07_Folded_Transit.pdf"
    )

    plt.close()

    print("Saved: Figure07_Folded_Transit")

def figure08_roc_curve():

    print("[Figure 08] ROC Curve")

    # ------------------------------------------------------
    # Load dataset
    # ------------------------------------------------------

    df = pd.read_csv(
        "data/dataset/final_dataset.csv"
    )

    X = df.drop(
        columns=[
            "tic_id",
            "label"
        ]
    )

    y = df["label"]

    # ------------------------------------------------------
    # Train/Test split
    # ------------------------------------------------------

    X_train, X_test, y_train, y_test = train_test_split(

        X,
        y,

        test_size=0.30,

        stratify=y,

        random_state=42
    )

    # ------------------------------------------------------
    # Train Random Forest
    # ------------------------------------------------------

    model = RandomForestClassifier(

        n_estimators=500,

        class_weight="balanced",

        random_state=42,

        n_jobs=-1
    )

    model.fit(
        X_train,
        y_train
    )

    # ------------------------------------------------------
    # Prediction probabilities
    # ------------------------------------------------------

    y_prob = model.predict_proba(X_test)[:,1]

    # ------------------------------------------------------
    # ROC
    # ------------------------------------------------------

    fpr, tpr, _ = roc_curve(
        y_test,
        y_prob
    )

    auc = roc_auc_score(
        y_test,
        y_prob
    )

    print(f"AUC = {auc:.4f}")

    # ------------------------------------------------------
    # Plot
    # ------------------------------------------------------

    setup_plot()

    plt.figure(figsize=(8,7))

    plt.plot(

        fpr,

        tpr,

        color=get_color("blue"),

        linewidth=3,

        label=f"Random Forest (AUC = {auc:.3f})"
    )

    plt.plot(

        [0,1],

        [0,1],

        linestyle="--",

        color=get_color("gray"),

        linewidth=3,

        label="Random classifier"
    )

    plt.xlim(0,1)

    plt.ylim(0,1.02)

    plt.xlabel(
        "False Positive Rate"
    )

    plt.ylabel(
        "True Positive Rate"
    )

    plt.title(
        "Receiver Operating Characteristic (ROC)"
    )

    plt.legend(
        loc="lower right"
    )

    finish_plot()

    plt.savefig(
        OUTPUT / "Figure08_ROC_Curve.png",
        dpi=600
    )

    plt.savefig(
        OUTPUT / "Figure08_ROC_Curve.pdf"
    )

    plt.close()

    print(
        "Saved: Figure08_ROC_Curve"
    )

    print()

def figure09_confusion_matrix():

    print("[Figure 09] Confusion Matrix")

    # --------------------------------------------------
    # Load dataset
    # --------------------------------------------------

    df = pd.read_csv(
        "data/dataset/final_dataset.csv"
    )

    X = df.drop(
        columns=[
            "tic_id",
            "label"
        ]
    )

    y = df["label"]

    # --------------------------------------------------
    # Train/Test split
    # --------------------------------------------------

    X_train, X_test, y_train, y_test = train_test_split(

        X,
        y,

        test_size=0.30,

        stratify=y,

        random_state=42

    )

    # --------------------------------------------------
    # Train Random Forest
    # --------------------------------------------------

    model = RandomForestClassifier(

        n_estimators=500,

        class_weight="balanced",

        random_state=42,

        n_jobs=-1

    )

    model.fit(
        X_train,
        y_train
    )

    y_pred = model.predict(X_test)

    # --------------------------------------------------
    # Confusion Matrix
    # --------------------------------------------------

    cm = confusion_matrix(
        y_test,
        y_pred
    )

    setup_plot()

    fig, ax = plt.subplots(
        figsize=(7,6)
    )

    disp = ConfusionMatrixDisplay(

        confusion_matrix=cm,

        display_labels=[
            "Negative",
            "Positive"
        ]

    )

    disp.plot(

        cmap="Blues",

        colorbar=False,

        ax=ax

    )

    plt.title(
        "Confusion Matrix"
    )

    finish_plot()

    plt.savefig(

        OUTPUT / "Figure09_Confusion_Matrix.png",

        dpi=600

    )

    plt.savefig(

        OUTPUT / "Figure09_Confusion_Matrix.pdf"

    )

    plt.close()

    print("Saved: Figure09_Confusion_Matrix")

def figure10_feature_importance():

    print("[Figure 10] Feature Importance")

    # --------------------------------------------------
    # Load dataset
    # --------------------------------------------------

    df = pd.read_csv(
        "data/dataset/final_dataset.csv"
    )

    X = df.drop(
        columns=[
            "tic_id",
            "label"
        ]
    )

    y = df["label"]

    # --------------------------------------------------
    # Train/Test split
    # --------------------------------------------------

    X_train, X_test, y_train, y_test = train_test_split(

        X,
        y,

        test_size=0.30,

        stratify=y,

        random_state=42

    )

    # --------------------------------------------------
    # Train Random Forest
    # --------------------------------------------------

    model = RandomForestClassifier(

        n_estimators=500,

        class_weight="balanced",

        random_state=42,

        n_jobs=-1

    )

    model.fit(
        X_train,
        y_train
    )

    # --------------------------------------------------
    # Feature importance
    # --------------------------------------------------

    importance = pd.DataFrame({

        "Feature": X.columns,

        "Importance": model.feature_importances_

    })

    importance = importance.sort_values(

        "Importance",

        ascending=False

    )

    # Mostrar las 15 variables más importantes

    importance = importance.head(15)

    feature_names = {
        "n_points": "Number of Samples",
        "snr_estimate": "Transit SNR",
        "median_flux": "Median Flux",
        "best_period": "Orbital Period",
        "mad_flux": "MAD Flux",
        "transit_depth": "Transit Depth",
        "std_flux": "Flux Std. Dev.",
        "amplitude": "Flux Amplitude",
        "bls_power": "BLS Power",
        "mean_flux": "Mean Flux",
        "transit_duration": "Transit Duration"
    }

    importance["Feature"] = importance["Feature"].map(feature_names)
    setup_plot()

    plt.figure(figsize=(9,6))

    plt.barh(

        importance["Feature"],

        importance["Importance"],

        color=get_color("blue")

    )

    plt.gca().invert_yaxis()

    plt.xlabel("Importance")

    plt.ylabel("Feature")

    plt.title("Top Feature Importances")

    finish_plot()

    plt.savefig(

        OUTPUT / "Figure10_Feature_Importance.png",

        dpi=600

    )

    plt.savefig(

        OUTPUT / "Figure10_Feature_Importance.pdf"

    )

    plt.close()

    print("Saved: Figure10_Feature_Importance")

def figure11_precision_recall():

    print("[Figure 11] Precision-Recall Curve")

    from sklearn.model_selection import train_test_split
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import (
        precision_recall_curve,
        average_precision_score
    )

    df = pd.read_csv(
        "data/dataset/final_dataset.csv"
    )

    X = df.drop(
        columns=[
            "tic_id",
            "label"
        ]
    )

    y = df["label"]

    X_train, X_test, y_train, y_test = train_test_split(

        X,

        y,

        test_size=0.30,

        random_state=42,

        stratify=y

    )

    model = RandomForestClassifier(

        n_estimators=500,

        class_weight="balanced",

        random_state=42,

        n_jobs=-1

    )

    model.fit(

        X_train,

        y_train

    )

    y_prob = model.predict_proba(X_test)[:,1]

    precision, recall, _ = precision_recall_curve(

        y_test,

        y_prob

    )

    ap = average_precision_score(

        y_test,

        y_prob

    )

    setup_plot()

    plt.figure(figsize=(7,6))

    plt.plot(

        recall,

        precision,

        color=get_color("blue"),

        linewidth=2.5,

        label=f"Random Forest (AP = {ap:.3f})"

    )

    plt.xlabel("Recall")

    plt.ylabel("Precision")

    plt.title("Precision–Recall Curve for Exoplanet Detection")

    plt.legend()

    finish_plot()

    plt.savefig(

        OUTPUT/"Figure11_Precision_Recall.png",

        dpi=600

    )

    plt.savefig(

        OUTPUT/"Figure11_Precision_Recall.pdf"

    )

    plt.close()

    print(f"Average Precision = {ap:.4f}")

    print("Saved: Figure11_Precision_Recall")

def figure12_shap_summary():

    print("[Figure 12] SHAP Summary Plot")

    import shap

    from sklearn.model_selection import train_test_split

    from sklearn.ensemble import RandomForestClassifier

    # ------------------------------------------------------
    # Load dataset
    # ------------------------------------------------------

    df = pd.read_csv(
        "data/dataset/final_dataset.csv"
    )

    X = df.drop(
        columns=[
            "tic_id",
            "label"
        ]
    )

    y = df["label"]

    # ------------------------------------------------------
    # Train/Test split
    # ------------------------------------------------------

    X_train, X_test, y_train, y_test = train_test_split(

        X,

        y,

        test_size=0.30,

        random_state=42,

        stratify=y

    )

    # ------------------------------------------------------
    # Train Random Forest
    # ------------------------------------------------------

    model = RandomForestClassifier(

        n_estimators=500,

        class_weight="balanced",

        random_state=42,

        n_jobs=-1

    )

    model.fit(

        X_train,

        y_train

    )

    # ------------------------------------------------------
    # SHAP values
    # ------------------------------------------------------

    explainer = shap.TreeExplainer(model)

    shap_values = explainer.shap_values(X_test)

    # ------------------------------------------------------
    # Plot
    # ------------------------------------------------------

    plt.figure(figsize=(10,7))

    if isinstance(shap_values, list):

        values = shap_values[1]

    else:

        values = shap_values[:, :, 1]
    
    feature_names = {
        "n_points": "Number of Samples",
        "snr_estimate": "Transit SNR",
        "median_flux": "Median Flux",
        "best_period": "Orbital Period",
        "transit_depth": "Transit Depth",
        "std_flux": "Flux Std. Dev.",
        "mad_flux": "MAD Flux",
        "amplitude": "Flux Amplitude",
        "bls_power": "BLS Power",
        "mean_flux": "Mean Flux",
        "transit_duration": "Transit Duration"
    }

    X_plot = X_test.rename(columns=feature_names)

    shap.summary_plot(
        values,
        X_test,
        show=False,
        max_display=15
    )

    plt.title("SHAP Summary of Feature Contributions")

    plt.tight_layout()

    plt.savefig(

        OUTPUT / "Figure12_SHAP_Summary.png",

        dpi=600,

        bbox_inches="tight"

    )

    plt.savefig(

        OUTPUT / "Figure12_SHAP_Summary.pdf",

        bbox_inches="tight"

    )

    plt.close()

    print("Saved: Figure12_SHAP_Summary")
# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 60)
    print("PUBLICATION FIGURES")
    print("=" * 60)
    print()
    print("Output folder:")
    print(OUTPUT)
    print()
    figure01_pipeline()
    figure02_dataset_distribution()
    figure03_period_distribution()
    figure04_transit_depth_distribution()
    figure05_lightcurve()
    figure06_bls_periodogram()
    figure07_folded_transit()
    figure08_roc_curve()
    figure09_confusion_matrix()
    figure10_feature_importance()
    figure11_precision_recall()
    figure12_shap_summary()
    print()
    print("Finished.")
if __name__ == "__main__":
    main()