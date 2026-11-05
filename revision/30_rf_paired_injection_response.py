#!/usr/bin/env python3
"""Experiment 30: source-disjoint Random Forest response to existing TESS
physical transit injections. GitHub-only, no HPC/FITS/MAST/holdout photometry.

Uses the already published 447 control-injection feature PAIRS from
experiment 06 and the 42 TICs frozen as its EVALUATION cohort in experiment
08 (126 pairs). The 42 evaluation TICs and the 12 examined pilot TICs are
EXCLUDED from RF fitting and RF ordinary source-label testing.

The original catalogue feature table (369 sources) supplies RF training.
Training and ordinary 30% test use ONLY the other 315 TICs (220 labelled
hosts + 95 catalogue-not-confirmed stars). The 42 paired-injection TICs
remain separate for transfer evaluation.

The eleven legacy RF features are reconstructed exactly by name:
  mean_flux, median_flux, std_flux, mad_flux, amplitude, n_points,
  best_period, transit_duration, transit_depth, bls_power, snr_estimate.
The injection BLS power was computed with a constant dy=noise_sigma_baseline;
convert to UNWEIGHTED units via bls_power*dy**2. Legacy snr_estimate
is abs(transit_depth)/std_flux, NOT Astropy depth_snr.

The full, no-n_points and five-BLS-feature RFs each use n_estimators=500,
class_weight='balanced', random_state=42, and a fixed .5 threshold.
This is an EXPLORATORY TRANSFER EVALUATION: the injected samples come
from originally unconfirmed-catalogue sources (not vetted transit-free
stars), the injection preprocessing and BLS grid differ somewhat from
the original catalogue processing, and binary RF probability is a
CATALOGUE-HOST score, NOT a calibrated transit probability.

Main endpoint: paired injected-minus-control RF output for each source.
Source-TIC clustered CIs, not 126 independent stars.
No previously frozen 18-host final benchmark is accessed or rerun.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (roc_auc_score, average_precision_score,
                             balanced_accuracy_score, f1_score,
                             confusion_matrix)
from sklearn.model_selection import train_test_split

ROOT=Path(__file__).resolve().parents[1]
CAT=ROOT/"data/dataset/final_dataset.csv"
PAIRS=ROOT/"results/injection_full/paired_features.csv"
TRUTH=ROOT/"results/injection_full/injection_truth_and_recovery.csv"
COHORT=ROOT/"results/injection_validation/source_star_split.csv"
PILOT=ROOT/"results/injection_pilot/injection_truth_and_recovery.csv"
OUT=ROOT/"results/rf_paired_injection_response"
SEED=42
BOOT_SEED=20261105
BOOT=3000
FEATURES=["mean_flux","median_flux","std_flux","mad_flux","amplitude",
          "n_points","best_period","transit_duration","transit_depth",
          "bls_power","snr_estimate"]
MODELS={
  "FULL_11":FEATURES,
  "NO_N_POINTS":[x for x in FEATURES if x!="n_points"],
  "BLS_ONLY_5":["best_period","transit_duration","transit_depth",
                 "bls_power","snr_estimate"],
}
RADIUS_BINS=[.034999,.07,.11,.150001]
RADIUS_NAMES=["0.035–0.070","0.070–0.110","0.110–0.150"]

def sha256(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()

def ensure_finite(df,columns,name):
    if not np.isfinite(df[columns].to_numpy(float)).all():
        raise ValueError("Nonfinite "+name+" feature inputs")

def load_inputs():
    original=pd.read_csv(CAT,dtype={"tic_id":str})
    pairs=pd.read_csv(PAIRS,dtype={"tic_id":str,"pair_id":str})
    truth=pd.read_csv(TRUTH,dtype={"tic_id":str,"pair_id":str})
    cohort=pd.read_csv(COHORT,dtype={"tic_id":str})
    pilot=pd.read_csv(PILOT,dtype={"tic_id":str})
    return original,pairs,truth,cohort,pilot

def validate(original,pairs,truth,cohort,pilot):
    if len(original)!=369 or original.tic_id.nunique()!=369:
        raise ValueError("Expected 369 unique original labelled source TICs")
    if not (set(original.label.astype(int))=={0,1} and
            int((original.label==1).sum())==220 and
            int((original.label==0).sum())==149):
        raise ValueError("Original catalogue source labels changed")
    if len(pairs)!=894 or pairs.tic_id.nunique()!=149:
        raise ValueError("Missing original 447 paired injection trials")
    if len(truth)!=447 or truth.pair_id.duplicated().any():
        raise ValueError("Injection truth/provenance changed")
    if cohort.tic_id.nunique()!=137 or len(cohort)!=137:
        raise ValueError("Missing 95/42 fixed injection validation cohort")
    if set(cohort.split)!={"calibration","evaluation"}:
        raise ValueError("Unexpected injection split labels")
    test=set(cohort.loc[cohort.split=="evaluation","tic_id"])
    cal=set(cohort.loc[cohort.split=="calibration","tic_id"])
    pilots=set(pilot.tic_id)
    if len(test)!=42 or len(cal)!=95 or len(pilots)!=12:
        raise ValueError("Expected 42 eval/95 cal/12 pilot source TICs")
    if test&cal or test&pilots or cal&pilots:
        raise ValueError("Leakage between injection source cohorts")
    neg=set(original.loc[original.label==0,"tic_id"])
    if test|cal|pilots != neg:
        raise ValueError("Injection cohorts do not equal original 149 negative-catalogue TICs")
    if (pairs.duplicated(["pair_id","label"]).any() or
        not pairs.groupby("pair_id").agg(n=("label","size"),s=("label","sum")).eval(
            "(n == 2) & (s == 1)"
        ).all()):
        raise ValueError("Every injection must have exactly one control and one injection")
    if not truth.pair_id.isin(pairs.pair_id).all():
        raise ValueError("Injection truth has unpaired entries")
    if not pairs.groupby("pair_id").tic_id.nunique().eq(1).all():
        raise ValueError("One pair spans multiple TICs")
    for col in ("n_points","baseline_days","noise_sigma_baseline"):
        pivot=pairs.pivot(index="pair_id",columns="label",values=col)
        if not np.allclose(pivot[0].astype(float),pivot[1].astype(float),
                           rtol=0,atol=1e-9):
            raise ValueError("Control/injected mismatch in "+col)
    for col in FEATURES:
        if col not in original.columns:
            raise ValueError("Original RF feature missing: "+col)
    for col in [*FEATURES[:10],"bls_depth_snr","noise_sigma_baseline"]:
        if col not in pairs.columns:
            raise ValueError("Paired injected feature missing: "+col)
    if "radius_ratio" not in truth.columns:
        raise ValueError("No physical injection radius provenance")
    return test,cal,pilots

def harmonize_weighting(pairs):
    """Restore original Astropy unweighted BLS 'power' from dy-weighted BLS.

    A CONSTANT sigma supplied for every cadence multiplies BLS log-likelihood
    power by 1/sigma**2. The reconstructed power is in the original
    dimensionless FLUX**2 units, but original BLS grids/preprocessing differ.
    """
    out=pairs.copy()
    sigma=out.noise_sigma_baseline.astype(float)
    if not (np.isfinite(sigma).all() and (sigma>0).all()):
        raise ValueError("Invalid BLS dy in injection feature provenance")
    out["bls_power_weighted"]=out.bls_power.astype(float)
    out["bls_power"]=(out.bls_power.astype(float)*sigma**2)
    out["snr_estimate"]=np.abs(out.transit_depth.astype(float)) / (
        out.std_flux.astype(float)
    )
    if (out.std_flux.astype(float)<=0).any():
        raise ValueError("Bad standard deviation for legacy SNR")
    ensure_finite(out,FEATURES,"injection")
    if (out.bls_power<0).any():
        raise ValueError("Invalid unweighted BLS power")
    return out

def prepare(original,pairs,truth,test,cal,pilots):
    # Fit on 315 catalogue sources, with NO injection-evaluation/pilot TIC.
    pool=original[~original.tic_id.isin(test|pilots)].copy()
    if (len(pool)!=315 or pool.tic_id.nunique()!=315 or
        not pool.tic_id.isin(test|pilots).eq(False).all()):
        raise ValueError("Unexpected excluded source sizes")
    ensure_finite(pool,FEATURES,"original")
    # Exactly the same source-specific split for all RF model variants.
    idx=np.arange(len(pool))
    train_i,normal_i=train_test_split(
        idx,test_size=.30,random_state=SEED,stratify=pool.label
    )
    train=pool.iloc[train_i].copy()
    normal=pool.iloc[normal_i].copy()
    if set(train.tic_id)&set(normal.tic_id) or (
        set(normal.tic_id)&test or set(train.tic_id)&test
    ):
        raise ValueError("Source leaked between normal and injection eval")
    if not set(normal.label)=={0,1}:
        raise ValueError("Normal heldout has only one catalogue class")
    paired=pairs[pairs.tic_id.isin(test)].copy()
    truth_test=truth[truth.tic_id.isin(test)].copy()
    if (len(paired)!=252 or len(truth_test)!=126
        or paired.tic_id.nunique()!=42):
        raise ValueError("Expected 126 source-held-out physical injection pairs")
    if set(paired.pair_id)!=set(truth_test.pair_id):
        raise ValueError("Physical injection pairs and truth do not match")
    for tic,g in paired.groupby("tic_id"):
        if len(g)!=6 or (g.label==0).sum()!=3 or (g.label==1).sum()!=3:
            raise ValueError("Incomplete 3-injection source pair: "+tic)
    prepared=harmonize_weighting(paired)
    return train,normal,prepared,truth_test

def bootstrap_tic(frame,metric,boot=BOOT,seed=BOOT_SEED):
    # One per-TIC average for point estimates and intervals.
    averages=frame.groupby("tic_id")[metric].mean().sort_index().to_numpy(float)
    if len(averages)!=42 or not np.isfinite(averages).all():
        raise ValueError("Bootstrap requires 42 TIC clusters")
    rng=np.random.default_rng(seed)
    sample=rng.integers(0,len(averages),size=(boot,len(averages)))
    draws=averages[sample].mean(axis=1)
    return [float(averages.mean()),
            float(np.quantile(draws,.025)),
            float(np.quantile(draws,.975))]

def classify_and_compare(train,normal,features,truth,model_name,boot):
    sub=features.copy()
    raw=pd.Series(index=sub.index,dtype=float)
    mdl=RandomForestClassifier(
        n_estimators=500,class_weight="balanced",
        random_state=SEED,n_jobs=2
    )
    mdl.fit(train[MODELS[model_name]],train.label.astype(int))
    p=mdl.predict_proba(normal[MODELS[model_name]])[:,1]
    y=normal.label.astype(int).to_numpy()
    pred=(p>=.5).astype(int)
    tn,fp,fn,tp=confusion_matrix(y,pred,labels=[0,1]).ravel()
    normal_score={
        "model":model_name,"n_features":len(MODELS[model_name]),
        "train_n_source_TICs":len(train),
        "normal_test_n_source_TICs":len(normal),
        "normal_catalogue_AUC":float(roc_auc_score(y,p)),
        "normal_catalogue_average_precision":float(average_precision_score(y,p)),
        "normal_catalogue_F1":float(f1_score(y,pred,zero_division=0)),
        "normal_catalogue_balanced_accuracy":float(balanced_accuracy_score(y,pred)),
        "normal_TN":int(tn),"normal_FP":int(fp),
        "normal_FN":int(fn),"normal_TP":int(tp)
    }
    sub["rf_probability"]=mdl.predict_proba(
        sub[MODELS[model_name]]
    )[:,1]
    c=sub[sub.label==0].set_index("pair_id")
    z=sub[sub.label==1].set_index("pair_id")
    if not c.index.equals(z.index):
        z=z.loc[c.index]
    if not (c.tic_id.astype(str).equals(z.tic_id.astype(str))):
        raise ValueError("Paired TIC order mismatch")
    df=pd.DataFrame({
        "model":model_name,"pair_id":c.index.astype(str),
        "tic_id":c.tic_id.to_numpy(str),
        "control_probability":c.rf_probability.to_numpy(float),
        "injected_probability":z.rf_probability.to_numpy(float),
        "n_points_original":c.n_points.to_numpy(float),
        "n_points_injected":z.n_points.to_numpy(float),
        "bls_power_unweighted_control":c.bls_power.to_numpy(float),
        "bls_power_unweighted_injected":z.bls_power.to_numpy(float),
    })
    df=df.merge(truth[["pair_id","tic_id","radius_ratio","true_period_days",
                       "period_recovered_direct"]],
                on=["pair_id","tic_id"],validate="one_to_one")
    if len(df)!=126:
        raise ValueError("Expected precisely 126 paired predictions")
    df["delta_probability"]=(
        df.injected_probability-df.control_probability
    )
    df["control_RF_positive"]=df.control_probability>=.5
    df["injected_RF_positive"]=df.injected_probability>=.5
    df["flip_negative_to_positive"]=(
        ~df.control_RF_positive & df.injected_RF_positive
    )
    df["flip_positive_to_negative"]=(
        df.control_RF_positive & ~df.injected_RF_positive
    )
    df["positive_delta"]=df.delta_probability>0.
    df["radius_bin"]=pd.cut(
        df.radius_ratio,RADIUS_BINS,labels=RADIUS_NAMES,
        include_lowest=True
    ).astype(str)
    if df.radius_bin.eq("nan").any():
        raise ValueError("Injection radius outside original physical design")
    if not np.allclose(df.n_points_original,df.n_points_injected,
                       atol=0,rtol=0):
        raise ValueError("n_points source shortcut changed across pair")
    del_ci=bootstrap_tic(df,"delta_probability",boot,BOOT_SEED)
    plus_ci=bootstrap_tic(df,"positive_delta",boot,BOOT_SEED+1)
    flip_ci=bootstrap_tic(df,"flip_negative_to_positive",boot,BOOT_SEED+2)
    control_ci=bootstrap_tic(df,"control_RF_positive",boot,BOOT_SEED+3)
    inject_ci=bootstrap_tic(df,"injected_RF_positive",boot,BOOT_SEED+4)
    summary={
        "model":model_name,
        "n_injection_pairs":len(df),
        "n_distinct_evaluation_stars":df.tic_id.nunique(),
        "mean_delta_RF_probability":del_ci[0],
        "mean_delta_cluster_CI95_low":del_ci[1],
        "mean_delta_cluster_CI95_high":del_ci[2],
        "median_pair_delta":float(df.delta_probability.median()),
        "positive_delta_fraction":plus_ci[0],
        "positive_delta_cluster_CI95_low":plus_ci[1],
        "positive_delta_cluster_CI95_high":plus_ci[2],
        "control_RF_positive_fraction":control_ci[0],
        "injected_RF_positive_fraction":inject_ci[0],
        "negative_to_positive_flip_fraction":flip_ci[0],
        "negative_to_positive_flip_CI95_low":flip_ci[1],
        "negative_to_positive_flip_CI95_high":flip_ci[2],
        "positive_to_negative_flips":int(df.flip_positive_to_negative.sum()),
        "direct_BLS_recovery_fraction_in_same_pairs":float(
            df.period_recovered_direct.mean()
        ),
        "note":"RF scores are catalogue-host outputs, not transit probability"
    }
    return normal_score,df,summary

def figure_plot(paired,out):
    main=paired[paired.model=="FULL_11"].copy()
    fig,ax=plt.subplots(figsize=(7.5,5),layout="constrained")
    bins=[main[main.radius_bin==k].delta_probability for k in RADIUS_NAMES]
    ax.boxplot(bins,tick_labels=RADIUS_NAMES,showfliers=True)
    ax.axhline(0,ls="--",color="black",lw=1)
    ax.set(xlabel="Injected radius ratio (planet / assumed solar-like star)",
           ylabel="RF injected minus control host-label probability",
           title="Paired RF transfer response to physical transit injections")
    fig.savefig(out/"figure_rf_paired_probability_by_radius.png",dpi=220)
    fig.savefig(out/"figure_rf_paired_probability_by_radius.pdf")
    plt.close(fig)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--outdir",type=Path,default=OUT)
    ap.add_argument("--bootstrap",type=int,default=BOOT)
    args=ap.parse_args()
    if not 300<=args.bootstrap<=30000:
        raise ValueError("Bootstrap replicates must be 300..30000")
    original,pairs,truth,cohort,pilot=load_inputs()
    test,cal,pilots=validate(original,pairs,truth,cohort,pilot)
    train,normal,paired,truth_test=prepare(
        original,pairs,truth,test,cal,pilots
    )
    args.outdir.mkdir(parents=True,exist_ok=True)
    scores=[]; predictions=[]; effects=[]
    for model_name in MODELS:
        norm,pred,effect=classify_and_compare(
            train,normal,paired,truth_test,model_name,args.bootstrap
        )
        scores.append(norm)
        predictions.append(pred)
        effects.append(effect)
        print(
            f"[PASS] {model_name} normal AUC={norm['normal_catalogue_AUC']:.3f} "
            f"paired mean shift={effect['mean_delta_RF_probability']:.4f} "
            f"on {effect['n_distinct_evaluation_stars']} withheld TICs",
            flush=True
        )
    normal_frame=pd.DataFrame(scores)
    pair_frame=pd.concat(predictions,ignore_index=True)
    effects_frame=pd.DataFrame(effects)
    normal_frame.to_csv(args.outdir/"normal_catalogue_test_metrics.csv",index=False)
    pair_frame.to_csv(args.outdir/"paired_126_injection_RF_predictions.csv",index=False)
    effects_frame.to_csv(args.outdir/"source_clustered_RF_effect_summary.csv",index=False)
    radius=pair_frame.groupby(["model","radius_bin"],observed=True).agg(
        n_pairs=("pair_id","size"),
        n_stars=("tic_id","nunique"),
        mean_delta_probability=("delta_probability","mean"),
        fraction_positive_delta=("positive_delta","mean"),
        control_RF_positive=("control_RF_positive","mean"),
        injected_RF_positive=("injected_RF_positive","mean"),
        flip_negative_to_positive=("flip_negative_to_positive","sum"),
        direct_BLS_period_recovered=("period_recovered_direct","sum"),
    ).reset_index()
    radius.to_csv(args.outdir/"paired_RF_effect_by_radius.csv",index=False)
    figure_plot(pair_frame,args.outdir)
    manifest={
      "experiment":30,
      "method":"Exploratory source-disjoint RF transfer evaluation on published paired physical injections",
      "source":"Pre-existing experiment 06/08 447 paired trials; 42 TICs/126 paired trials held for evaluation",
      "n_original_RF_sources":369,
      "n_removed_from_RF_fit_and_RF_normal_test":54,
      "n_RF_source_pool_after_exclusions":315,
      "n_RF_training_sources":len(train),
      "n_RF_normal_testing_sources":len(normal),
      "n_injected_testing_source_TICs":42,
      "n_paired_injection_evaluation_trials":126,
      "models":MODELS,
      "model_params":{"n_estimators":500,"class_weight":"balanced",
                      "random_state":42,"threshold":.5},
      "first_random_source_partition": "remove 42 prior injection-eval and 12 pilot stars; stratified 70/30 of remaining 315",
      "same_RF_per_model_evaluates_both_normal_and_injections":True,
      "original_RF_feature_powers": "Legacy original RF model had Astropy unweighted BLS power",
      "reconstructed_BLS_power": "original paired dy-weighted power * noise_sigma_baseline**2",
      "legacy_RF_SNR_formula":"abs(transit_depth)/std_flux",
      "bootstrap_unit":"stellar TIC (42); all three injections retained together",
      "bootstrap_replicates":args.bootstrap,
      "locked_18_final_holdout_not_accessed":True,
      "source_hashes":{str(f.relative_to(ROOT)):sha256(f)
                       for f in [CAT,PAIRS,TRUTH,COHORT,PILOT]},
      "limitations":[
        "The paired injection sources come from the catalogue-negative group, not verified transit-free stars.",
        "The RF was originally built to discriminate confirmed-host catalogue labels, not injection presence.",
        "Injection BLS preprocessing and trial-period grid differ from those used to build original RF features; power is rescaled for unit comparability, not complete pipeline equivalence.",
        "The original RF sampling/cadence shortcut can persist even after physical injections.",
        "RF probabilities are not calibrated astrophysical transit-detection probabilities.",
        "Normal-test AUC is for catalogue labels among other stars; it is not an injection-detection AUC.",
        "Only 42 source clusters, not 126 independent stars, define the injection transfer population.",
        "This is a separate experiment from the 108 TGLC injections in manuscript section 4.7 and from the 18-host frozen final test.",
        "All results postdate the 18-host frozen benchmark and cannot be used to change its algorithm."
      ]
    }
    (args.outdir/"experiment_protocol.json").write_text(
      json.dumps(manifest,indent=2,ensure_ascii=False)+"\n",encoding="utf8"
    )
    mainrow=effects_frame.set_index("model").loc["FULL_11"]
    typical=normal_frame.set_index("model").loc["FULL_11"]
    report=[
      "# Experiment 30: Random Forest response to physical transit injections",
      "",
      "EXPLORATORY complementary injection-transfer analysis. This is not a rerun "
      "of the 108 TGLC-host injections or the locked 18-host final evaluation.",
      "",
      f"Original source catalogue: 369 unique TICs; 42 injection test TICs "
      f"and 12 pilot TICs excluded before RF model fitting.",
      f"Remaining RF source pool: 315; training {len(train)}; normal test {len(normal)}.",
      "Independent paired transfer test: 42 stars, 126 synthetic control/injected pairs.",
      "",
      "## Ordinary catalogued-label test (same trained models)",
      normal_frame.round(5).to_string(index=False),
      "",
      "## Paired physical-injection effect (42 source TIC bootstrap clusters)",
      effects_frame.round(5).to_string(index=False),
      "",
      "## By injected radius-ratio bin",
      radius.round(5).to_string(index=False),
      "",
      "## Interpretation",
      f"The full 11-feature RF obtained a catalogue-label AUC of "
      f"{typical.normal_catalogue_AUC:.3f} on its separate normal heldout test, "
      f"but on 126 paired injected-minus-original samples its mean score shift was "
      f"{mainrow.mean_delta_RF_probability:+.4f} (source bootstrap 95% interval "
      f"{mainrow.mean_delta_cluster_CI95_low:+.4f} to "
      f"{mainrow.mean_delta_cluster_CI95_high:+.4f}). "
      "A difference in these endpoints is not a direct difference of RF AUC.",
      "",
      *["- "+z for z in manifest["limitations"]]
    ]
    (args.outdir/"RF_injection_transfer_report.md").write_text(
       "\n".join(report)+"\n",encoding="utf8"
    )
    print("\n"+"\n".join(report),flush=True)

if __name__=="__main__":
    main()
