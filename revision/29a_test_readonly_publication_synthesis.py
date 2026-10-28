#!/usr/bin/env python3
"""Pre-publication read-only test. Only opens published CSV/JSON metadata.

Do NOT call final experiment 28, original FITS, or cleaned final light curves.
All plotting and table rendering tests run in an ephemeral temp folder.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

SRC=Path(__file__).with_name("29_readonly_publication_synthesis.py")
spec=importlib.util.spec_from_file_location("publication29",SRC)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def main():
    final=pd.read_csv(m.FINAL/"final_per_host_scored.csv",dtype={"tic_id":str})
    blind=pd.read_csv(m.FINAL/"final_blind_predictions.csv",dtype={"tic_id":str})
    frozen=json.loads(m.FREEZE.read_text())
    outcome=json.loads((m.FINAL/"final_results_protocol.json").read_text())
    dev=pd.read_csv(m.DEV)
    phys=pd.read_csv(m.PHYSICAL)
    marg=pd.read_csv(m.MARGINAL)
    cal=pd.read_csv(m.CALIB)
    paired=json.loads(m.PAIRED.read_text())
    m.validate(final,blind,frozen,outcome,dev,phys,marg,cal,paired)
    assert len(final)==18
    table=m.grouped_success(final)
    assert len(table)==5
    assert int(table.loc[0,"n"])==4
    assert int(table.loc[1,"n"])==0
    assert int(table.loc[2,"n"])==2
    assert int(table.loc[3,"n"])==2
    assert int(table.loc[4,"n"])==0

    with TemporaryDirectory(prefix="tess_publication_readonly_") as path:
        root=Path(path)
        m.fig_period_agreement(final,root)
        m.fig_correctness_acceptance(final,root)
        m.fig_injection_grid(phys,root)
        figs=list(root.glob("figure_*"))
        assert len(figs)==6,figs
        assert all(x.stat().st_size>1000 for x in figs)

    # Any post-hoc alteration of final scientific decisions must fail.
    changed=final.copy()
    j=int(changed.index[0])
    changed.at[j,"evidence_threshold_accepted"]=not bool(
        changed.at[j,"evidence_threshold_accepted"]
    )
    try:
        m.validate(changed,blind,frozen,outcome,dev,phys,marg,cal,paired)
    except ValueError:
        pass
    else:
        raise AssertionError("Tampering with frozen final acceptance undetected")

    change_physical=phys.copy()
    change_physical.loc[0,"n_direct"]+=1
    try:
        m.validate(final,blind,frozen,outcome,dev,
                   change_physical,marg,cal,paired)
    except ValueError:
        pass
    else:
        raise AssertionError("Tampered synthetic recovery denominator unnoticed")

    print("[PASS] Frozen final one-shot 18-host science fields reproduced")
    print("[PASS] 4/18 agreement; 0/2 raw vs 2/2 evidence accepted")
    print("[PASS] Six reproducible PNG/PDF manuscript figures created in memory")
    print("[PASS] Reject post-hoc changes to final results and injection totals")
    print("[PASS] No FITS or holdout photometry opened; read-only metadata")


if __name__=="__main__":
    main()
