"""
Phase 3.7 — Parse PLIP XML reports and summarize interactions.
"""
import os
import xml.etree.ElementTree as ET
import pandas as pd
from collections import defaultdict

LIGANDS = ['EAI045', 'EAI001', 'decoy_0086', 'decoy_0454']
BASE = "data/md/plip/results"

# Interaction types PLIP detects
INTERACTION_TAGS = {
    'hydrophobic_interactions': 'Hydrophobic',
    'hydrogen_bonds': 'HBond',
    'water_bridges': 'WaterBridge',
    'salt_bridges': 'SaltBridge',
    'pi_stacks': 'PiStack',
    'pi_cation_interactions': 'PiCation',
    'halogen_bonds': 'Halogen',
    'metal_complexes': 'Metal',
}

all_interactions = []

for LIG in LIGANDS:
    xml_path = f"{BASE}/{LIG}/{LIG}_dry_report.xml"
    if not os.path.exists(xml_path):
        print(f"❌ {LIG}: no XML at {xml_path}")
        continue
    
    print(f"\n{'='*60}")
    print(f"  {LIG}")
    print(f"{'='*60}")
    
    tree = ET.parse(xml_path)
    root = tree.getroot()
    
    # Find the binding site
    for site in root.findall('.//bindingsite'):
        # Get ligand info
        lig_id = site.find('.//ligand')
        lig_name = lig_id.get('name') if lig_id is not None else 'UNK'
        
        # For each interaction type
        interactions_found = {}
        for tag, label in INTERACTION_TAGS.items():
            for interaction in site.findall(f'.//{tag}/*'):
                # Extract residue info
                resnr  = interaction.findtext('resnr', '?')
                restype = interaction.findtext('restype', '?')
                reschain = interaction.findtext('reschain', '?')
                
                key = f"{restype}{resnr}"
                all_interactions.append({
                    'Ligand': LIG,
                    'Residue': key,
                    'Chain': reschain,
                    'Interaction': label,
                })
                interactions_found[label] = interactions_found.get(label, 0) + 1
        
        # Print summary
        print(f"  Ligand: {lig_name}")
        print(f"  Interactions detected:")
        for label, count in sorted(interactions_found.items()):
            print(f"    {label:15s}: {count}")

# Save full table
df = pd.DataFrame(all_interactions)
os.makedirs("data/md/plip/results", exist_ok=True)
df.to_csv("data/md/plip/results/plip_interactions.csv", index=False)

# Build pivot: which residues interact with which ligands
if not df.empty:
    matrix = df.groupby(['Residue', 'Interaction', 'Ligand']).size().unstack(fill_value=0)
    matrix.to_csv("data/md/plip/results/plip_matrix.csv")
    
    print(f"\n{'='*60}")
    print(f"  PLIP ANALYSIS COMPLETE")
    print(f"{'='*60}")
    print(f"  Total interactions logged: {len(df)}")
    print(f"  Unique residues involved:  {df['Residue'].nunique()}")
    print(f"  Full matrix: data/md/plip/results/plip_matrix.csv")
else:
    print(f"\n⚠️  No interactions found in any system")
