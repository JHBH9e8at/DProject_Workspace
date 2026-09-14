"""Run one analysis workflow or an ordered JSON plan from one block entry point."""
from __future__ import annotations
import argparse, json
from pathlib import Path
from .workflows.workflow_runner import WORKFLOWS

def run_plan(path: str|Path) -> list[dict]:
    source=Path(path).resolve(strict=True); payload=json.loads(source.read_text(encoding="utf-8"))
    commands=payload.get("commands") if isinstance(payload,dict) else None
    if not isinstance(commands,list) or not commands: raise ValueError("Plan must contain a non-empty commands list")
    completed=[]
    for index,command in enumerate(commands,1):
        if not isinstance(command,dict): raise ValueError(f"Command {index} must be an object")
        name=command.get("workflow"); arguments=command.get("args",[])
        if name not in WORKFLOWS: raise ValueError(f"Command {index} has unknown workflow: {name!r}")
        if not isinstance(arguments,list) or not all(isinstance(value,str) for value in arguments): raise ValueError(f"Command {index} args must be a list of strings")
        print(f"[{index}/{len(commands)}] {name}")
        code=WORKFLOWS[name](arguments)
        if code not in (None,0): raise RuntimeError(f"Workflow {name} returned {code}")
        completed.append({"index":index,"workflow":name,"status":"completed"})
    return completed

def main(argv=None):
    parser=argparse.ArgumentParser(description="Run the Analysis/Figure block")
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--workflow",choices=sorted(WORKFLOWS)); mode.add_argument("--plan",type=Path)
    parser.add_argument("workflow_args",nargs=argparse.REMAINDER)
    args=parser.parse_args(argv)
    if args.plan:
        if args.workflow_args: parser.error("workflow arguments cannot be used with --plan")
        completed=run_plan(args.plan); print(json.dumps({"status":"completed","commands":completed},indent=2)); return 0
    workflow_args=args.workflow_args[1:] if args.workflow_args[:1]==["--"] else args.workflow_args
    return WORKFLOWS[args.workflow](workflow_args)

if __name__=="__main__": raise SystemExit(main())
