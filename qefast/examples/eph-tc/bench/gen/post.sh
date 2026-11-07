#!/usr/bin/env bash
# Copyright (C) 2026 The qe-fast authors
# SPDX-License-Identifier: GPL-2.0-or-later
# after the queue: regression on B1/B2 smoke, then the diagnostics
G=/root/gen; Q=/root/qefast/qe-fast; P=$Q/qefast/examples/eph-tc/bench/phq.sh; B2=$Q/build-gpu/bin
source /root/qefast/env.sh; cd $G/p1
while ! grep -q ALL_DONE $G/all.log; do sleep 5; done
echo "== regression $(date -u +%H:%M:%S)"
( export QEF_PRESET=max; bash $P scfB1 r1 $B2 1 3 ) 2>&1 | tail -1
for i in 1 2 3; do cmp r1/MgB2.dyn$i t1g/MgB2.dyn$i && echo "B1 q$i max: identical to the earlier run"; done
( export QEF_PRESET=max; export QEF_STH_CBSR_RESERVE=1000; bash $P scfB1 r2 $B2 1 3 ) 2>&1 | tail -1
for i in 1 2 3; do cmp r2/MgB2.dyn$i t1cg/MgB2.dyn$i && echo "B1 q$i forced decline: identical to CG"; done
grep -h 'qef cbsr auto' r1/ph.out r2/ph.out | cut -c1-260
( export QEF_PRESET=max; bash $P scfB2 r3 $B2 1 2 ) 2>&1 | tail -1
grep -h 'qef cbsr auto' r3/ph.out | cut -c1-260
rm -rf r?/_ph0
echo "== ablation G6 $(date -u +%H:%M:%S)"; bash $G/diag.sh G6
echo "== shared gauge $(date -u +%H:%M:%S)"; bash $G/sg.sh G6 safe fast max; bash $G/sg.sh G3 safe fast max; bash $G/sg.sh G5 safe fast max
echo POST_DONE
