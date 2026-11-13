"""Replay 106 type II witnesses without modifying the shipped reference reports."""
from pathlib import Path
import argparse,importlib.util,json,shutil

ROOT=Path(__file__).resolve().parents[1]

def main():
    if not __debug__:raise SystemExit('Do not use Python -O')
    ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args()
    a.output_dir.mkdir(parents=True,exist_ok=False)
    paper=ROOT/'census'
    for name in ('preparation.json','results.json','input.scan'):shutil.copy2(paper/'data/type_ii'/name,a.output_dir/name)
    spec=importlib.util.spec_from_file_location('portable_type_ii',paper/'tools/type_ii_census.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    module.PAPER=paper;module.OUT=a.output_dir
    module.verify(native=False)
    result=json.loads((a.output_dir/'replay.json').read_text())
    assert result['verified_count']==106 and result['target_complete'] and result['missing']==[]
    print('PASS: 106 type II witnesses, independent region reader and syndrome checks.')

if __name__=='__main__':main()
