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


def apply_gradient(coefs, eris, flag='', grad_symbol='n'):

    assert flag in ['a', 'b']

    new_coefs, new_eris = [], []

    for c, e in zip(coefs, eris):

        if flag == 'a':
            c_s, e_s = e.apply_gradient_a(grad_symbol=grad_symbol)
        elif flag == 'b':
            c_s, e_s = e.apply_gradient_b(grad_symbol=grad_symbol)

        for c2, e2 in zip(c_s, e_s):
            if c == '1':
                new_coefs.append(f'{c2}')
            elif c2 == '1':
                new_coefs.append(f'{c}')
            else:
                new_coefs.append(f'{c} * {c2}')
            new_eris.append(e2)

    return new_coefs, new_eris


def apply_hrr_a(coefs, eris):

    return apply_hrr(coefs, eris, flag='a')


def apply_hrr_b(coefs, eris):

    return apply_hrr(coefs, eris, flag='b')


def apply_hrr_a_once(coefs, eris):

    return apply_hrr(coefs, eris, flag='a_once')


def apply_hrr_b_once(coefs, eris):

    return apply_hrr(coefs, eris, flag='b_once')


def apply_hrr(coefs, eris, flag=''):

    new_coefs, new_eris = [], []

    for c, e in zip(coefs, eris):

        if flag == 'a':
            c_s, e_s = e.apply_hrr_a()
        elif flag == 'b':
            c_s, e_s = e.apply_hrr_b()
        elif flag == 'a_once':
            c_s, e_s = e.apply_hrr_a_once()
        elif flag == 'b_once':
            c_s, e_s = e.apply_hrr_b_once()

        for c2, e2 in zip(c_s, e_s):
            if c == '1':
                new_coefs.append(f'{c2}')
            elif c2 == '1':
                new_coefs.append(f'{c}')
            else:
                new_coefs.append(f'{c} * {c2}')
            new_eris.append(e2)

    return new_coefs, new_eris


def simplify_coef(c):

    new_coef_terms = []

    numerator = []
    denominator = []
    sign = 1
    denom_2 = 1

    for t in c.split('*'):
        term = t.strip()
        if term == '1':
            continue
        elif term == '(-1)':
            sign *= -1
        elif '/' in term:
            x, y = term.split('/')
            if x != '1':
                numerator.append(x)
            if y == '2':
                denom_2 *= 2
            else:
                denominator.append(y)
        #elif term == '(-2.0)' or term.startswith('one'):
        elif term == '(-2.0)' or term == '-2':
            new_coef_terms.append('(-2.0)')
        elif term.isdigit():
            new_coef_terms.append(term)
        elif term.startswith('delta_'):
            content = term.split('_')
            assert content[0] == 'delta'
            new_term = '_'.join(content[:1] + sorted(content[1:]))
            new_coef_terms.append(new_term)
        elif (term.startswith('A') or term.startswith('C')
              or term.startswith('P') or term.startswith('Q')
              or term.startswith('G')):
            new_coef_terms.append(term)
        elif term in ['ksi', 'zeta', 'a_i', 'a_j', 'S1', 'zeta_c']:
            new_coef_terms.append(term)
        else:
            print()
            print('c:   ', c)
            print('term:', term)
            print()
            assert False

    while True:
        found_x_y = False
        ix, jy = -1, -1
        for i, x in enumerate(numerator):
            for j, y in enumerate(denominator):
                if x == y:
                    found_x_y = True
                    ix = i
                    jy = j
                    break
        if not found_x_y:
            break
        else:
            numerator.pop(ix)
            denominator.pop(jy)

    new_coef_terms = sorted(new_coef_terms)
    denominator = sorted(denominator)
    numerator = sorted(numerator)

    if denominator:
        if len(denominator) > 1:
            denominator_str = '( ' + ' * '.join(denominator) + ' )'
        else:
            denominator_str = ' * '.join(denominator)
        if numerator:
            if len(numerator) > 1:
                numerator_str = '( ' + ' * '.join(numerator) + ' )'
            else:
                numerator_str = ' * '.join(numerator)
            new_coef_terms = [numerator_str + ' / ' + denominator_str
                              ] + new_coef_terms
        else:
            new_coef_terms = ['1 / ' + denominator_str] + new_coef_terms

    if denom_2 > 1:
        new_coef_terms = [f'{sign*1/denom_2}'] + new_coef_terms
    elif sign < 0:
        new_coef_terms = [f'{sign}'] + new_coef_terms

    return ' * '.join(new_coef_terms).replace(' * 1 / ', ' / ')
