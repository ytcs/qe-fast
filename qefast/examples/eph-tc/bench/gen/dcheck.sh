#!/usr/bin/env bash
# Copyright (C) 2026 The qe-fast authors
# SPDX-License-Identifier: GPL-2.0-or-later
# Fallback checks D1-D3: ph.x only. D1, D2: off vs max, q1-q2. D3: max to q2 under a time limit, then safe if the check chose CG.
G=/root/gen; Q=/root/qefast/qe-fast; B=$Q/build-gpu/bin; P=$Q/qefast/examples/eph-tc/bench/phq.sh; PS=$Q/qefast/pseudo/sg15
source /root/qefast/env.sh
pwr() { ( export QEF_PRESET=off; cd $1; OMP_NUM_THREADS=8 mpirun --bind-to none -np 1 -x OMP_NUM_THREADS -x QEF_PRESET $B/pw.x -in scf.in > scf.out 2>&1 ); }
for D in D1 D2 D3; do
  S=$G/run/$D/scf; rm -rf $G/run/$D; mkdir -p $S; cp $G/in/$D/scf.in $G/in/$D/ph.in $S/; sed -i "s|@PSEUDO@|$PS|" $S/scf.in
  t0=$(date +%s); pwr $S; echo "$D scf rc=$? $(( $(date +%s)-t0 )) s; $(grep -E 'Sym|total magnetization' $S/scf.out | tail -2 | tr -s ' ' | tr '\n' ';' | cut -c1-120)"
  if [ $D = D3 ]; then
    ( export QEF_PRESET=max; timeout -k 20 2400 bash $P $S $G/run/D3/max $B ); rc=$?
    echo "D3 max rc=$rc"; grep -E 'qef cbsr auto|sth_engine: (needs|irrep)|Error in|CG kept' $G/run/D3/max/ph.out | cut -c1-400 | head -6
    pkill -9 -f '[p]h\.x'
    if grep -q 'qef cbsr auto: CG' $G/run/D3/max/ph.out; then
      ( export QEF_PRESET=safe; timeout -k 20 2400 bash $P $S $G/run/D3/safe $B ); echo "D3 safe rc=$?"
      for f in $G/run/D3/max/x.dyn[0-9]*; do cmp $f $G/run/D3/safe/$(basename $f) && echo "D3 $(basename $f) identical to safe"; done
    fi
  else
    ( export QEF_PRESET=off; bash $P $S $G/run/$D/off $B ); echo "$D off rc=$?"
    ( export QEF_PRESET=max; bash $P $S $G/run/$D/max $B ); echo "$D max rc=$?"
    for f in $G/run/$D/off/x.dyn[0-9]*; do cmp $f $G/run/$D/max/$(basename $f) && echo "$D $(basename $f) identical"; done
    python3 $Q/qefast/examples/eph-tc/bench/cmpgate.py $G/run/$D/off $G/run/$D/max x | tail -4
    echo "$D qef lines in max:"; grep -E 'qef|sth_engine|Error' $G/run/$D/max/ph.out | sed 's/[0-9.]\+/#/g' | sort | uniq -c | sort -rn | head -12 | cut -c1-130
  fi
  rm -rf $G/run/$D/*/_ph0 $G/run/$D/*/x.save
done
echo DCHECK_DONE
