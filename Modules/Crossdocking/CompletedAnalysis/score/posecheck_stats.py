"""Summarize PoseCheck populations and deterministic docking top-1% cohorts."""
from __future__ import annotations
import argparse, math
from pathlib import Path
import pandas as pd
from Modules.Analysis_Figen_block.common.output_naming import allocate_versioned_group
from Modules.Analysis_Figen_block.common.paths import validate_results_root
from Modules.Analysis_Figen_block.common.provenance import update_manifest
from Modules.Analysis_Figen_block.common.run_context import create_run_context, generate_run_id

SCRIPT_ID="PY-162"; STATES=("PPS","PR")
REQUIRED={"population_state","pose_index","molecule_id","docking_score","heavy_atom_count","rotatable_bond_count","clash_count","clashes_per_heavy_atom","strain_energy","strain_per_rotatable_bond"}
QUALITY=("clash_count","clashes_per_heavy_atom","strain_energy","strain_per_rotatable_bond")
SPECS=[("pose_count",None,"count_rows","count"),("docking_score_mean","docking_score","mean","Glide score"),("docking_score_median","docking_score","median","Glide score"),("docking_score_q05","docking_score","q05","Glide score"),("docking_score_q95","docking_score","q95","Glide score"),("heavy_atom_count_mean","heavy_atom_count","mean","atoms"),("heavy_atom_count_median","heavy_atom_count","median","atoms"),("rotatable_bond_count_mean","rotatable_bond_count","mean","bonds"),("rotatable_bond_count_median","rotatable_bond_count","median","bonds")]
for base,column,unit in (("clash_count","clash_count","count"),("clashes_per_heavy_atom","clashes_per_heavy_atom","ratio"),("strain_energy","strain_energy","PoseCheck strain energy"),("strain_per_rotatable_bond","strain_per_rotatable_bond","ratio")):
    SPECS.extend((f"{base}_{s}",column,s,unit) for s in ("mean","median","q05","q25","q75","q95","max"))
SPECS.append(("strain_per_rotatable_bond_missing","strain_per_rotatable_bond","count_missing","count"))

def _stat(series,stat):
    values=pd.to_numeric(series,errors="coerce")
    if stat=="count_missing": return int(values.isna().sum())
    if stat.startswith("q"): return values.quantile(int(stat[1:])/100)
    return getattr(values,stat)()

def build_posecheck_stats(quality,*,top_fraction=.01):
    missing=sorted(REQUIRED-set(quality.columns))
    if missing: raise ValueError(f"Missing required columns: {missing}")
    if set(quality.population_state.dropna().unique())!=set(STATES): raise ValueError("Expected exactly PPS and PR population states")
    if not 0<top_fraction<=1: raise ValueError("top_fraction must be in (0, 1]")
    tops=[]; selections=[]
    for state in STATES:
        frame=quality.loc[quality.population_state==state].copy().sort_values(["docking_score","pose_index","molecule_id"],kind="mergesort")
        n=math.ceil(len(frame)*top_fraction); selected=frame.head(n).copy(); selected.insert(0,"top1pct_rank",range(1,n+1)); tops.append(selected)
        selections.append({"population":state,"full_population_n":len(frame),"selected_n":n,"selected_fraction":n/len(frame),"best_docking_score":selected.docking_score.min(),"top1pct_boundary_docking_score":selected.docking_score.max(),"next_unselected_docking_score":frame.iloc[n].docking_score if n<len(frame) else pd.NA,"selection_rule":"ceil(N*fraction), ascending docking_score; deterministic tie-break by pose_index then molecule_id"})
    top=pd.concat(tops,ignore_index=True); cohorts={"full_population":quality,"docking_top1pct":top}; long=[]
    for scope,cohort in cohorts.items():
        groups={s:cohort.loc[cohort.population_state==s] for s in STATES}
        for metric,column,stat,unit in SPECS:
            row={"scope":scope,"metric":metric,"unit":unit}
            for state,frame in groups.items(): row[state]=len(frame) if stat=="count_rows" else _stat(frame[column],stat)
            long.append(row)
    long=pd.DataFrame(long); summary=pd.DataFrame([{"metric":metric,**{f"{state}_{suffix}":long.loc[(long.scope==scope)&(long.metric==metric),state].iloc[0] for state in STATES for suffix,scope in (("full","full_population"),("T1","docking_top1pct"))}} for metric,*_ in SPECS])
    relationships=[]
    for scope,cohort in cohorts.items():
        for method in ("spearman","pearson"):
            for metric in QUALITY:
                row={"scope":scope,"correlation_method":method,"pose_quality_metric":metric}
                for state in STATES:
                    pair=cohort.loc[cohort.population_state==state,["docking_score",metric]].apply(pd.to_numeric,errors="coerce").dropna(); row[state]=pair.docking_score.corr(pair[metric],method=method); row[f"{state}_n"]=len(pair)
                relationships.append(row)
    return summary,pd.DataFrame(relationships),pd.DataFrame(selections),top

def run_posecheck_stats(*,source_csv,run_id,results_root=None,resume=False,top_fraction=.01,upstream_run_ids=()):
    context=create_run_context(("analysis","posecheck_stats"),run_id,results_root=results_root,resume=resume); config={"top_fraction":top_fraction}
    try:
        tables=build_posecheck_stats(pd.read_csv(source_csv,low_memory=False),top_fraction=top_fraction); paths=allocate_versioned_group(context.tables_dir,("pose_quality_summary_full_vs_top1pct.csv","docking_pose_quality_relationship_full_vs_top1pct.csv","top1pct_selection_summary.csv","pose_quality_top1pct_selected_poses.csv"))
        for table,path in zip(tables,paths): table.to_csv(path,index=False)
        update_manifest(context,script_id=SCRIPT_ID,status="completed",config=config,inputs=(Path(source_csv),),outputs=paths,upstream_run_ids=upstream_run_ids)
    except Exception:
        update_manifest(context,script_id=SCRIPT_ID,status="failed",config=config,upstream_run_ids=upstream_run_ids); raise
    return context,tables,paths

def main(argv=None):
    p=argparse.ArgumentParser(description="Summarize full and docking top-1% PoseCheck cohorts"); p.add_argument("--source",required=True,type=Path); p.add_argument("--top-fraction",type=float,default=.01); p.add_argument("--run-id"); p.add_argument("--results-root",type=Path); p.add_argument("--resume",action="store_true"); p.add_argument("--upstream-run-id",action="append",default=[]); a=p.parse_args(argv)
    c,_,paths=run_posecheck_stats(source_csv=a.source,top_fraction=a.top_fraction,run_id=a.run_id or generate_run_id("PR_PPS","posecheck_stats"),results_root=validate_results_root(a.results_root) if a.results_root else None,resume=a.resume,upstream_run_ids=tuple(a.upstream_run_id)); print(f"PoseCheck statistics run: {c.run_dir}"); [print(f"Table: {x}") for x in paths]; return 0
if __name__=="__main__": raise SystemExit(main())
