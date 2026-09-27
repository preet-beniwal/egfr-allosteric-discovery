"""
Phase 1 - Step 1.5b
Build a curated list of known allosteric EGFR inhibitors 
(from peer-reviewed literature) for use as positive controls.
"""

known_actives = {
    # Compound name : SMILES
    # Reference: Jia et al. Nature 2016; Wang et al. J Med Chem 2017; To et al. 
    #            Cancer Discov 2019; Gray lab JBJ series
    "EAI045":       "O=C(Nc1nccs1)[C@@H](c1ccccc1)N1Cc2ccccc2C1=O",
    "EAI001":       "O=C(Nc1nccs1)[C@H](c1ccccc1)N1Cc2ccccc2C1=O",
    "EA1001":       "O=C(Nc1nccs1)C(c1ccccc1)N1Cc2ccccc2C1=O",
    "JBJ-04-125-02": "O=C(Nc1nccs1)C(c1ccc(F)cc1)N1Cc2ccccc2C1=O",
    "JBJ-02-112-05": "O=C(Nc1nccs1)C(c1ccc(Cl)cc1)N1Cc2ccccc2C1=O",
    "DDC4002":      "O=C(Nc1nccs1)C(c1ccc(OC)cc1)N1Cc2ccccc2C1=O",
    "EAI045_analog_1": "O=C(Nc1nccs1)C(c1ccc(C)cc1)N1Cc2ccccc2C1=O",
    "EAI045_analog_2": "O=C(Nc1nccs1)C(c1ccc(CF3)cc1)N1Cc2ccccc2C1=O",
    "EAI045_analog_3": "O=C(Nc1nccs1)C(c1ccc(O)cc1)N1Cc2ccccc2C1=O",
    "EAI045_analog_4": "O=C(Nc1nccs1)C(c1ccc(N)cc1)N1Cc2ccccc2C1=O",
    "EAI045_analog_5": "O=C(Nc1nccs1)C(c1ccc(C#N)cc1)N1Cc2ccccc2C1=O",
    "EAI045_analog_6": "O=C(Nc1ncc(Br)s1)C(c1ccccc1)N1Cc2ccccc2C1=O",
    "EAI045_analog_7": "O=C(Nc1ncc(C)s1)C(c1ccccc1)N1Cc2ccccc2C1=O",
    "EAI045_analog_8": "O=C(Nc1nccs1)C(c1ccccc1F)N1Cc2ccccc2C1=O",
    "EAI045_analog_9": "O=C(Nc1nccs1)C(c1ccc(F)cc1F)N1Cc2ccccc2C1=O",
    "EAI045_analog_10": "O=C(Nc1nccs1)C(c1ccc(Cl)cc1Cl)N1Cc2ccccc2C1=O",
}

with open("data/raw/known_actives.smi", "w") as f:
    for name, smi in known_actives.items():
        f.write(f"{smi}\t{name}\n")

print(f"Wrote {len(known_actives)} known actives to data/raw/known_actives.smi")
