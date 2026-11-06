import os
import glob
import subprocess
import time
import sys
import re

from pyrosetta import *
from pyrosetta.rosetta.core.scoring.hbonds import HBondSet, fill_hbond_set
import pyrosetta.rosetta.protocols.minimization_packing
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover

def get_sequence_from_chainA(pdb_path):
    result = subprocess.run(['/wynton/home/kortemme/rdalal/tools/get_chainA.sh', pdb_path], capture_output=True, text=True, check=True)
    residues = result.stdout.strip().split('\n')
    return residues

def design_pdbs_fixed_designs_by_motif(rfdiffusionfiles):
    motif = "TDSNSIT"
    print(f"Fixing motif: {motif}")
    for pdb_path in rfdiffusionfiles:
        pose = pose_from_file(pdb_path)
        full_seq = pose.sequence()

        matches = [m.start() for m in re.finditer(f"(?={motif})", full_seq)]

        for m in matches:
            indices = list(range(m + 1, m + len(motif) + 1))
            print(f"Motif '{motif}' found at positions: {indices}")
            add_fixed_resis_label(indices, pdb_path)

def add_fixed_resis_label(fixresis, pdb):
    remarks = []
    for position in fixresis:
        remarks.append(f"REMARK PDBinfo-LABEL:{position: >5} FIXED")
    remarks_str = '\n'.join(remarks)
    with open(pdb, 'a') as f:
        f.write("\n" + remarks_str + "\n")

def submit_design(roundnum, pdbdir):
    outdir = f"{pdbdir}/design_round{roundnum}"
    os.makedirs(outdir, exist_ok=True)
    
    if roundnum == 1:
        cmd = (f"python /wynton/home/kortemme/rdalal/dl_binder_design/mpnn_fr/dl_interface_design.py -pdbdir {pdbdir}/ "
               f"-outpdbdir {outdir} "
               f"-seqs_per_struct {'20' if roundnum == 1 else '1'} "
               f"-temperature 0.1 "
               f"-relax_cycles 0 "
               f"-debug "
               f"-checkpoint_name {pdbdir}/check.point "
               f"-checkpoint_path /wynton/home/kortemme/rdalal/tools/ProteinMPNN/soluble_model_weights/v_48_020.pt")
    else:
        cmd = (f"python /wynton/home/kortemme/rdalal/dl_binder_design/mpnn_fr/dl_interface_design.py -pdbdir {pdbdir}/design_round{roundnum - 1}/ "
               f"-outpdbdir {outdir} "
               f"-seqs_per_struct {'1' if roundnum == 1 else '1'} "
               f"-temperature 0.1 "
               f"-relax_cycles 0 "
               f"-debug "
               f"-checkpoint_name {pdbdir}/check.point "
               f"-checkpoint_path /wynton/home/kortemme/rdalal/tools/ProteinMPNN/soluble_model_weights/v_48_020.pt")
