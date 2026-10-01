#!/bin/bash
MDP="/home/preet/projects/egfr-allosteric/data/md/mdp"
BASE="/home/preet/projects/egfr-allosteric/data/md/systems"
LOG="/home/preet/projects/egfr-allosteric/logs/auto_chain.log"

mkdir -p /home/preet/projects/egfr-allosteric/logs
exec >> "$LOG" 2>&1

echo ""
echo "############ AUTO-CHAIN STARTED: $(date) ############"

while pgrep -f "mdrun.*nvt" > /dev/null; do sleep 60; done
echo "[$(date)] NVT complete. Launching NPT..."

for LIG in EAI045 EAI001 decoy_0086 decoy_0454; do
    cd "$BASE/${LIG}"
    gmx grompp -f "$MDP/npt.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt \
        -p topol.top -o npt.tpr -maxwarn 2 > grompp_npt.log 2>&1
    if [ ! -f npt.tpr ]; then echo "❌ $LIG npt grompp"; continue; fi
    nohup gmx mdrun -v -deffnm npt -ntmpi 1 -ntomp 1 > mdrun_npt.log 2>&1 &
    echo "[$(date)] ✅ $LIG NPT launched"
done

while pgrep -f "mdrun.*npt" > /dev/null; do sleep 60; done
echo "[$(date)] NPT complete. Launching production MD..."

for LIG in EAI045 EAI001 decoy_0086 decoy_0454; do
    cd "$BASE/${LIG}"
    gmx grompp -f "$MDP/md.mdp" -c npt.gro -t npt.cpt \
        -p topol.top -o md.tpr -maxwarn 2 > grompp_md.log 2>&1
    if [ ! -f md.tpr ]; then echo "❌ $LIG md grompp"; continue; fi
    nohup gmx mdrun -v -deffnm md -ntmpi 1 -ntomp 1 > mdrun_md.log 2>&1 &
    echo "[$(date)] ✅ $LIG production MD launched"
done

echo "[$(date)] Auto-chain complete."
