#!/usr/bin/env bash
# chain.sh <material> <set> <time limit s>: one held-out chain (groups/compute/GEN-LOG.md). Writes <run>/<material>/<set>/ with
# timings.tsv, status (exit code of run-fork.sh; 124 = timeout), peakmem (MiB, nvidia-smi every 2 s).
M=$1; S=$2; LIM=$3
[ "$M" = G8 ] && LIM=5400   # GEN-LOG amendment 05:04 UTC
G=/root/gen; W=$G/run/$M/$S; Q=/root/rtsc/qe-rtsc
source /root/rtsc/env.sh; source $G/env/$S.env
[ "$M" = G1 ] && export RTSC_RUN_NOEPW=1
rm -rf $W; mkdir -p $W
( while :; do nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits; sleep 2; done ) > $W/mem.log 2>/dev/null &
MP=$!
t0=$(date +%s)
flock /root/rtsc/TIMING.lock timeout -k 20 $LIM bash $Q/rtsc/bench/run-fork.sh $G/in/$M /root/rtsc/bin-gen $W 8 8 > $W/run.log 2>&1
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
