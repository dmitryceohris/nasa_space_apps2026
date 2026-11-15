"""Fetch Gaia DR3 reference-star catalogues for the five primary targets."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'src')]
from processing.references import fetch_references
from processing.targets import load_targets

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'cache/references')
    args=p.parse_args(); results=[]
    for t in load_targets().values():
        try:
            path=fetch_references(t,args.output)
            result={'target':t.id,'status':'downloaded','manifest':str(path.with_suffix('.json'))}
        except Exception as error:
            result={'target':t.id,'status':'unavailable','error':str(error)}
        results.append(result);print(result,flush=True)
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'status.json').write_text(json.dumps(results,indent=2)+'\n')
