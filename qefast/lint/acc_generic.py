# Copyright (C) 2026 The qe-fast authors
# SPDX-License-Identifier: GPL-2.0-or-later
"""Flag generic calls that nvfortran may resolve to the device routine without the code asking for it.

Inside an OpenACC data region, nvfortran resolves a generic call whose argument is named in a data clause to the
specific with DEVICE dummies (invfft_y_gpu, mp_sum_*_gpu, ...), even outside host_data. Code written for the host
copy then silently works on the device copy. Example: lr_calc_dens transformed psic on the device and left the host
copy in G-space, dropping the ultrasoft term on GPU builds. Such calls should sit inside an explicit
'!$acc host_data use_device(...)', with the host/device updates around it spelled out.

  acc_generic.py [files...]     (default: every tracked .f90/.F90 file); exit 1 on findings
"""
import re
import subprocess
import sys

GENERICS = r'invfft|fwfft|mp_sum|mp_bcast|mp_max|mp_min'
CLAUSES = r'present_or_copyin|present_or_copyout|present_or_copy|copyin|copyout|copy|create|present|no_create'
data = re.compile(r'^\s*!\$acc\s+data\b(.*)', re.I)
end_data = re.compile(r'^\s*!\$acc\s+end\s+data', re.I)
host_data = re.compile(r'^\s*!\$acc\s+host_data', re.I)
end_host_data = re.compile(r'^\s*!\$acc\s+end\s+host_data', re.I)
unit = re.compile(r'^\s*(?:END\s+)?(?:SUBROUTINE|FUNCTION|CONTAINS)\b', re.I)
call = re.compile(r'^\s*(?:IF\s*\(.*\)\s*)?CALL\s+(' + GENERICS + r')\s*\((?:\s*[\'"]\w+[\'"]\s*,)?\s*([A-Za-z_]\w*)',
                  re.I)


def names(clause_text):
    """Variable names in the data clauses of one directive."""
    out = set()
    for m in re.finditer(r'\b(?:' + CLAUSES + r')\s*\(', clause_text, re.I):
        depth, j, buf = 1, m.end(), ''
        while j < len(clause_text) and depth:
            ch = clause_text[j]
            depth += (ch == '(') - (ch == ')')
            if depth == 1 or (depth == 0 and ch == ')'):
                buf += ch if depth else ''
            elif depth > 1:
                buf += ' '
            j += 1
        out |= {v.lower() for v in re.findall(r'(?:^|,)\s*([A-Za-z_]\w*)', buf.replace('(', ' '))}
    return out


def scan(path):
    found = []
    try:
        lines = open(path, errors='replace').read().split('\n')
    except OSError:
        return found
    stack, in_hd, i = [], 0, 0
    while i < len(lines):
        line = lines[i]
        if unit.match(line):
            stack, in_hd = [], 0
        m = data.match(line)
        if m:
            text = m.group(1)
            while text.rstrip().endswith('&') and i + 1 < len(lines):
                i += 1
                text = text.rstrip()[:-1] + lines[i].replace('!$acc', '').lstrip().lstrip('&')
            stack.append(names(text))
        elif end_data.match(line) and stack:
            stack.pop()
        elif host_data.match(line):
            in_hd += 1
        elif end_host_data.match(line):
            in_hd = max(0, in_hd - 1)
        else:
            c = call.match(line)
            if c and not in_hd and any(c.group(2).lower() in s for s in stack):
                found.append(f'{path}:{i + 1}: {line.strip()}')
        i += 1
    return found


files = sys.argv[1:] or subprocess.run(['git', 'ls-files', '*.f90', '*.F90'], capture_output=True,
                                        text=True, check=True).stdout.split()
hits = [h for f in files for h in scan(f)]
print('\n'.join(hits) if hits else f'acc_generic: {len(files)} files, no generic call on a data-clause array '
                                   'outside host_data')
sys.exit(1 if hits else 0)
