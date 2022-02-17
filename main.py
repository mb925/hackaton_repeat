import os
import gzip
import json
import pandas as pd
import argparse
import config as cfg
import numpy as np

if __name__ == '__main__':

    # execute for each uniprot,pdb



    parser = argparse.ArgumentParser(description="map pdb on uniprot through seqres")
    parser.add_argument('--pdb', '-p', help="identifier of the source pdb")
    parser.add_argument('--chain', '-c', help="identifier of the source chain")
    parser.add_argument('--start', '-s', help="index of annotation start")
    parser.add_argument('--end', '-e', help="index of annotation end")
    parser.add_argument('--uniprot', '-u', help="identifier of the target uniprot")
    parser.add_argument('--mappings', '-m', help="path to seqres directory")
    parser.add_argument('--output', '-o', help="path to output file")

    args = parser.parse_args()

    pdb, chain = args.pdb, args.chain
    uniprot = args.uniprot
    seqres_dir = args.mappings
    seqres_path = os.path.join(seqres_dir, ''.join(pdb[1:3]), '%s_seqres.mjson.gz' % pdb)

    if not os.path.isfile(seqres_path):
        raise FileNotFoundError('could not find mappings file %s' % seqres_path)

    data = []
    with gzip.open(seqres_path, 'rt') as f:
        for line in f:
            line = json.loads(line)
            data.append(line)

    ds = pd.DataFrame(data)

    # ds = (['pdb_id', 'pdb_chain', 'uniprot_id']).apply(lambda x: findResidues(x))
    ds = ds[ds.pdb_id == pdb]
    ds = ds[ds.pdb_chain == chain]
    ds = ds[ds.uniprot_id == uniprot]

    ds = ds.sort_values("seqres_index")

    start = ds[ds.pdb_residue_id == args.start].index[0]
    end = ds[ds.pdb_residue_id == args.end].index[-1]

    ds = ds.iloc[start:(end + 1), :]
    ds.to_csv(args.output, sep="\t")





