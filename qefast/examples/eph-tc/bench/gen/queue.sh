#!/usr/bin/env bash
# the pre-registered order; all-off first, then limit = min(2700, 2 x all-off + 300)
G=/root/gen
for M in "$@"; do
  if [ "$M" != G1 ] && [ ! -f $G/in/$M/epw.in ]; then echo "$(date -u +%H:%M:%S) SKIP $M: no epw.in (preparation failed; see prep.log)"; continue; fi
  A=2700
  bash $G/chain.sh $M off 2700
  if [ "$(cat $G/run/$M/off/status)" = 0 ]; then
    T=$(awk -F'\t' '$1=="total"{print int($2)}' $G/run/$M/off/timings.tsv); A=$((2*T+300)); [ $A -gt 2700 ] && A=2700
  fi
  for S in safe fast max; do bash $G/chain.sh $M $S $A; done
done
echo QUEUE_DONE "$@"
