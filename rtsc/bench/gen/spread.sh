#!/usr/bin/env bash
# the pre-registered spread: 3 repeats of G2's all-off ph.x, on the preparation scf
G=/root/gen; Q=/root/rtsc/qe-rtsc
source /root/rtsc/env.sh; source $G/env/alloff.env
sed -i "s|@PSEUDO@|$Q/rtsc/pseudo/sg15|" $G/prep/G2/scf/ph.in 2>/dev/null
cp $G/in/G2/ph.in $G/prep/G2/scf/
for i in 1 2 3; do bash $Q/rtsc/bench/phq.sh $G/prep/G2/scf $G/run/spread$i /root/rtsc/bin-gen > /dev/null 2>&1; rm -rf $G/run/spread$i/_ph0 $G/run/spread$i/x.save; done
cat $G/run/spread?/ph.time | tr '\n' ' ' | awk '{m=($1+$2+$3)/3; hi=$1;lo=$1; for(i=2;i<=3;i++){if($i>hi)hi=$i;if($i<lo)lo=$i}; printf "SPREAD G2 all-off ph.x: %s %s %s s; (max-min)/mean %.2f %%\n",$1,$2,$3,100*(hi-lo)/m}'
