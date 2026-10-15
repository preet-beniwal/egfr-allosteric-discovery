#!/bin/bash
# Resume 1 ns production MD:
#  - Skip completed systems
#  - Resume interrupted systems from checkpoint
#  - Start fresh for un-run systems
#
# Uses -cpi to continue from .cpt if present.

MDP="/home/preet/projects/egfr-allosteric/data/md/mdp/md_1ns.mdp"
BASE="/home/preet/projects/egfr-allosteric/data/md/systems"
LOG="/home/preet/projects/egfr-allosteric/logs/1ns_recovery.log"

mkdir -p /home/preet/projects/egfr-allosteric/logs
exec >> "$LOG" 2>&1

echo ""
echo "##########################################################"
echo "  1 ns PRODUCTION MD — RECOVERY: $(date)"
echo "##########################################################"

for LIG in EAI045 EAI001 decoy_0086 decoy_0454; do
    echo ""
    echo "[$(date)] ===== Processing: $LIG ====="
    cd "$BASE/${LIG}" || continue
    
    # Skip if already complete
    if [ -f md_1ns.gro ]; then
        echo "[$(date)] ✅ $LIG already complete — skipping"
        continue
    fi
    
    # Determine if we need to resume or start fresh
    if [ -f md_1ns.cpt ]; then
        echo "[$(date)] 🔄 $LIG has checkpoint — resuming"
        RESUME_FLAG="-cpi md_1ns.cpt"
        APPEND_FLAG="-append"
    else
        echo "[$(date)] 🆕 $LIG fresh start"
        RESUME_FLAG=""
        APPEND_FLAG=""
    fi
    
    # Build TPR if missing
    if [ ! -f md_1ns.tpr ]; then
        echo "[$(date)] Building TPR..."
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
    fi
    
    # Run MD with all 6 cores, resume if applicable
    echo "[$(date)] Running mdrun (resume=$([ -n "$RESUME_FLAG" ] && echo yes || echo no))..."
    gmx mdrun -v -deffnm md_1ns \
        $RESUME_FLAG $APPEND_FLAG \
        -ntmpi 1 -ntomp 6 > mdrun_1ns.log 2>&1
    
    if [ -f md_1ns.gro ]; then
        echo "[$(date)] ✅ $LIG production MD complete"
    else
        echo "[$(date)] ❌ $LIG production MD failed"
    fi
done

echo ""
echo "[$(date)] ##########################################################"
echo "[$(date)]  ALL SYSTEMS COMPLETE"
echo "[$(date)] ##########################################################"
