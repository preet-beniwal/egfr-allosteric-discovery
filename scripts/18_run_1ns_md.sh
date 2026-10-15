#!/bin/bash
# Sequential 1 ns production MD — one system at a time, 6 cores each
# Total estimated time: ~1.5 days

MDP="/home/preet/projects/egfr-allosteric/data/md/mdp/md_1ns.mdp"
BASE="/home/preet/projects/egfr-allosteric/data/md/systems"
LOG="/home/preet/projects/egfr-allosteric/logs/1ns_production.log"

mkdir -p /home/preet/projects/egfr-allosteric/logs
exec >> "$LOG" 2>&1

echo ""
echo "##########################################################"
echo "  1 ns PRODUCTION MD STARTED: $(date)"
echo "##########################################################"

for LIG in EAI045 EAI001 decoy_0086 decoy_0454; do
    echo ""
    echo "[$(date)] ===== Starting: $LIG ====="
    cd "$BASE/${LIG}" || continue
    
    # Skip if already done
    if [ -f md_1ns.gro ]; then
        echo "[$(date)] ⏭️  $LIG already complete — skipping"
        continue
    fi
    
    # Build TPR
    gmx grompp -f "$MDP" \
        -c npt.gro \
        -t npt.cpt \
        -p topol.top \
        -o md_1ns.tpr \
        -maxwarn 2 > grompp_1ns.log 2>&1
    
    if [ ! -f md_1ns.tpr ]; then
        echo "[$(date)] ❌ $LIG grompp failed"
        tail -10 grompp_1ns.log
        continue
    fi
    echo "[$(date)] ✅ TPR built"
    
    # Run MD with all 6 cores
    gmx mdrun -v -deffnm md_1ns -ntmpi 1 -ntomp 6 > mdrun_1ns.log 2>&1
    
    if [ -f md_1ns.gro ]; then
        echo "[$(date)] ✅ $LIG production MD complete"
    else
        echo "[$(date)] ❌ $LIG production MD failed"
    fi
done

echo ""
echo "[$(date)] ##########################################################"
echo "[$(date)]  ALL 4 SYSTEMS COMPLETE"
echo "[$(date)] ##########################################################"
