#!/usr/bin/env bash
# Copyright (C) 2026 The qe-fast authors
# SPDX-License-Identifier: GPL-2.0-or-later
# chain.sh <material> <preset> <time limit s>: one held-out chain. Writes <run>/<material>/<preset>/ with
# timings.tsv, status (exit code of run-chain.sh; 124 = timeout), peakmem (MiB, nvidia-smi every 2 s).
M=$1; S=$2; LIM=$3
[ "$M" = G8 ] && LIM=5400
G=/root/gen; W=$G/run/$M/$S; Q=/root/qefast/qe-fast
source /root/qefast/env.sh; export QEF_PRESET=$S
[ "$M" = G1 ] && export QEF_RUN_NOEPW=1
# G8's EPW is OOM-killed on a 50 GB pod; if safe dies the same way, the later chains skip it
if [ "$M" = G8 ] && { [ "$S" = fast ] || [ "$S" = max ]; } && [ "$(cat $G/run/G8/safe/status)" != 0 ]; then export QEF_RUN_NOEPW=1; fi
rm -rf $W; mkdir -p $W
( while :; do nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits; sleep 2; done ) > $W/mem.log 2>/dev/null &
MP=$!
t0=$(date +%s)
flock /root/qefast/TIMING.lock timeout -k 20 $LIM bash $Q/qefast/examples/eph-tc/bench/run-chain.sh $G/in/$M $Q/build-gpu/bin $W 8 8 > $W/run.log 2>&1
rc=$?
t1=$(date +%s)
kill $MP 2>/dev/null
if [ $rc -ne 0 ]; then
  pkill -9 -f '[p]h\.x' ; pkill -9 -f '[p]w\.x'; pkill -9 -f '[e]pw\.x'; echo quit | nvidia-cuda-mps-control 2>/dev/null; sleep 3
fi
echo $rc > $W/status; sort -n $W/mem.log | tail -1 > $W/peakmem
du -sm $W | cut -f1 > $W/disk
# keep the disk bounded
rm -rf $W/_ph0 $W/x.epmatwp $W/*.epb* $W/*.ephmat* 2>/dev/null   # save/, x.save, x.mmn, x.amn stay: EPW-only diagnostics on a shared gauge
echo "$(date -u +%H:%M:%S) $M $S rc=$rc wall=$((t1-t0)) peak=$(cat $W/peakmem) MiB"
