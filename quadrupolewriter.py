#  This file is part of vlx-onee-ints.
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2024-2025 Xin Li
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from quadrupole import Quadrupole
from intsutils import apply_hrr_a, apply_hrr_b, simplify_coef
from intswriter import write_integrals


def write_quadrupole(ab, flag=''):

    comp_dict_1 = {
        's': [],
        'p': ['a0'],
        'd': ['a0', 'a1'],
        'f': ['a0', 'a1', 'a2'],
    }

    comp_dict_2 = {
        's': [],
        'p': ['b0'],
        'd': ['b0', 'b1'],
        'f': ['b0', 'b1', 'b2'],
    }

    quadrupole_ab = Quadrupole(comp_dict_1[ab[0]], comp_dict_2[ab[1]], 'm', 'n')

    flip_ab = False

    # apply HRR on a

    coefs, eris = quadrupole_ab.apply_hrr_a()

    # apply HRR on b

    coefs, eris = apply_hrr_b(coefs, eris)

    # simplify coefficients

    final_list = []
    for ind, (c, e) in enumerate(zip(coefs, eris)):
        final_list.append((simplify_coef(c), ind, e))

    coefs = []
    eris = []
    for (c, ind, e) in sorted(final_list):
        coefs.append(c)
        eris.append(e)

    # reorganize terms

    final_ints_lines = write_integrals(coefs, eris, indent=20)

    # finalize the result

    final_ints_str = ''

    for line in final_ints_lines:

        if '_M_00_m ( S S )' in line:
            line = line.replace('_M_00_m ( S S )', 'PC[m] * (')
        elif '_M_00_n ( S S )' in line:
            line = line.replace('_M_00_n ( S S )', 'PC[n] * (')
        elif '_Q_00_m_n ( S S )' in line:
            line = line.replace('_Q_00_m_n ( S S )', '(PC[n] * PC[m] + 0.5 / (a_i + a_j) * delta[n][m]) * (')
        elif '( S S )' in line:
            line = line.replace('( S S )', '(')

        if line.endswith('__nonewline__'):
            final_ints_str += line.split('__nonewline__')[0]
        else:
            final_ints_str += line + '\n'

    final_ints_str = final_ints_str.replace(r'(1.0 * (-1.0))', '(-1.0)')
    final_ints_str = final_ints_str.replace(r'(1.0 * (-2.0))', '(-2.0)')
    final_ints_str = final_ints_str.replace(r'(1.0 * (-3.0))', '(-3.0)')

    return final_ints_str
