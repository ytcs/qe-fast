#!/usr/bin/env bash
# after the pre-registered queue: regression of the corrected binary on B1/B2 smoke, then the labelled diagnostics
G=/root/gen; Q=/root/qefast/qe-fast; P=$Q/qefast/examples/eph-tc/bench/phq.sh; B2=/root/qefast/bin-gen2
source /root/qefast/env.sh; cd $G/p1
while ! grep -q ALL_DONE $G/all.log; do sleep 5; done
echo "== bin-gen2 regression $(date -u +%H:%M:%S)"
( source $G/env/max.env; bash $P scfB1 r1 $B2 1 3 ) 2>&1 | tail -1
for i in 1 2 3; do cmp r1/MgB2.dyn$i t1g/MgB2.dyn$i && echo "B1 q$i max: bin-gen2 identical to bin-gen"; done
( source $G/env/max.env; export QEF_STH_CBSR_RESERVE=1000; bash $P scfB1 r2 $B2 1 3 ) 2>&1 | tail -1
for i in 1 2 3; do cmp r2/MgB2.dyn$i t1cg/MgB2.dyn$i && echo "B1 q$i forced decline: identical to CG"; done
grep -h 'qef cbsr auto' r1/ph.out r2/ph.out | cut -c1-260
( source $G/env/max.env; bash $P scfB2 r3 $B2 1 2 ) 2>&1 | tail -1
grep -h 'qef cbsr auto' r3/ph.out | cut -c1-260
rm -rf r?/_ph0
echo "== ablation G6 $(date -u +%H:%M:%S)"; bash $G/diag.sh G6
echo "== shared gauge $(date -u +%H:%M:%S)"; bash $G/sg.sh G6 safe fast max; bash $G/sg.sh G3 safe fast max; bash $G/sg.sh G5 safe fast max
echo POST_DONE
