"""Stable full/last-window and top-candidate hotspot table adapters."""

from pathlib import Path
import pandas as pd

SCRIPT_ID="PY-130"

def summarize_full_last(interactions_csv, metadata_csv, total_poses, last_start, chunksize=250_000):
    cols=["population_state","molecule_id","residue_name","original_residue_number","source_step"]
    parts=[]
    for chunk in pd.read_csv(interactions_csv,usecols=cols,chunksize=chunksize):
        chunk=chunk.dropna(subset=cols).copy(); chunk["source_step"]=chunk.source_step.astype(int)
        chunk["residue"]=chunk.residue_name.astype(str)+chunk.original_residue_number.astype(int).astype(str)
        parts.append(chunk[["population_state","molecule_id","source_step","residue"]].drop_duplicates())
    unique=pd.concat(parts,ignore_index=True).drop_duplicates()
    full=unique.groupby(["residue","population_state"]).molecule_id.nunique().rename("poses").reset_index()
    full["prevalence"]=full.apply(lambda r:r.poses/total_poses[r.population_state],axis=1)
    recent=unique[unique.apply(lambda r:r.source_step>=last_start[r.population_state],axis=1)]
    meta=pd.read_csv(metadata_csv,usecols=["population_state","molecule_id","source_step","pose_read_status"])
    meta=meta[meta.pose_read_status.eq("valid") & meta.apply(lambda r:int(r.source_step)>=last_start[r.population_state],axis=1)]
    denominators={state:int(meta.loc[meta.population_state.eq(state),"molecule_id"].nunique()) for state in ("PPS","PR")}
    if any(v<=0 for v in denominators.values()): raise ValueError(f"Empty last-window denominator: {denominators}")
    last=recent.groupby(["residue","population_state"]).molecule_id.nunique().rename("poses").reset_index()
    last["prevalence"]=last.apply(lambda r:r.poses/denominators[r.population_state],axis=1)
    return full,last,denominators

def clean_residue(series): return series.astype(str).str.replace(r"\.[A-Za-z0-9]+$","",regex=True)

def load_top_candidates(summary_csv,transitions_csv):
    summary=pd.read_csv(summary_csv); transitions=pd.read_csv(transitions_csv)
    summary["own_state"]=summary.own_state.astype(str); transitions["residue"]=clean_residue(transitions.residue)
    transitions["molecule_id"]=transitions.molecule_id.astype(str)
    denominators=summary.groupby("own_state").molecule_id.nunique().astype(int).to_dict(); frames=[]
    for condition,flag in (("own","own_present"),("opposite","opposite_present")):
        present=transitions[transitions[flag].astype(str).str.lower().eq("true")]
        counts=present[["own_state","molecule_id","residue"]].drop_duplicates().groupby(["residue","own_state"]).molecule_id.nunique().rename("poses").reset_index()
        counts["prevalence"]=counts.apply(lambda r:r.poses/denominators[r.own_state],axis=1)
        frames.append(counts.rename(columns={"own_state":"population_state"}).assign(condition=condition))
    own_typed=transitions[transitions.own_present.astype(str).str.lower().eq("true")].drop_duplicates(["own_state","molecule_id","residue","interaction_type"])
    typed=own_typed.groupby(["residue","own_state","interaction_type"]).molecule_id.nunique().rename("poses").reset_index().rename(columns={"own_state":"population_state"})
    typed["prevalence"]=typed.apply(lambda r:r.poses/denominators[r.population_state],axis=1)
    return frames[0],frames[1],typed,denominators

