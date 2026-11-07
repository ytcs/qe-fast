#!/usr/bin/env bash
# Copyright (C) 2026 The qe-fast authors
# SPDX-License-Identifier: GPL-2.0-or-later
G=/root/gen; cd $G
while ! grep -q PREP_DONE prep.log; do sleep 5; done
bash queue.sh G1 G2
bash spread.sh
bash queue.sh G3 G4 G6 G5 G7 G8
bash dcheck.sh
echo ALL_DONE
