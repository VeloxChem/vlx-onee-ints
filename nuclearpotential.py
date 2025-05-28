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

from intsutils import apply_hrr_a_once, apply_hrr_b_once
from overlap import Overlap


class NuclearPotential:

    def __init__(self, cart_a, cart_b, bf_order=0):

        assert isinstance(bf_order, int)
        assert bf_order >= 0

        self.cart_a = list(cart_a)
        self.cart_b = list(cart_b)
        self.bf_order = bf_order

    def __repr__(self):

        label_a = self.get_angmom_str(self.cart_a)
        label_b = self.get_angmom_str(self.cart_b)

        return f'_A_00 ( {label_a} {label_b} )^{self.bf_order}'

    @staticmethod
    def get_angmom_str(cart_comp):

        labels = 'SPDFGHI'
        if len(cart_comp) == 0:
            return labels[len(cart_comp)]
        else:
            return labels[len(cart_comp)] + '_' + '_'.join(cart_comp)

    @property
    def La(self):

        return len(self.cart_a)

    @property
    def Lb(self):

        return len(self.cart_b)

    def apply_gradient_a(self, grad_symbol='m'):

        # Note: This only contains the orbital part of gradient and does not
        # include the operator part of gradient. The latter can be obtained
        # from ElectricField integral.

        # Eq.(5), Obara-Saika JCP 1986

        coef_s = ['2 * a_i']
        eri_s = [
            NuclearPotential(self.cart_a + [grad_symbol], self.cart_b, self.bf_order)
        ]

        for ind_a in range(len(self.cart_a)):

            coef_s.append(f'(-1) * delta_{self.cart_a[ind_a]}_{grad_symbol}')

            new_cart_a = list(self.cart_a)
            new_cart_a.pop(ind_a)
            eri_s.append(
                NuclearPotential(new_cart_a, self.cart_b, self.bf_order))

        return coef_s, eri_s

    def apply_gradient_b(self, grad_symbol='n'):

        # Note: This only contains the orbital part of gradient and does not
        # include the operator part of gradient. The latter can be obtained
        # from ElectricField integral.

        # Eq.(5), Obara-Saika JCP 1986

        coef_s = ['2 * a_j']
        eri_s = [
            NuclearPotential(self.cart_a, self.cart_b + [grad_symbol], self.bf_order)
        ]

        for ind_b in range(len(self.cart_b)):

            coef_s.append(f'(-1) * delta_{self.cart_b[ind_b]}_{grad_symbol}')

            new_cart_b = list(self.cart_b)
            new_cart_b.pop(ind_b)
            eri_s.append(
                NuclearPotential(self.cart_a, new_cart_b, self.bf_order))

        return coef_s, eri_s

    def apply_hrr_a_once(self):

        if self.La == 0:

            coef_s = ['1']
            eri_s = [NuclearPotential(self.cart_a, self.cart_b, self.bf_order)]

        else:

            # Eq.(A19), Obara-Saika JCP 1986

            coef_s = [f'PA_{self.cart_a[-1]}']
            eri_s = [
                NuclearPotential(self.cart_a[:-1], self.cart_b, self.bf_order)
            ]

            coef_s.append(f'(-1) * PC_{self.cart_a[-1]}')
            eri_s.append(
                NuclearPotential(self.cart_a[:-1], self.cart_b,
                                 self.bf_order + 1))

            for ind_a in range(len(self.cart_a) - 1):

                coef_s.append(
                    f'1/2 * 1/S1 * delta_{self.cart_a[ind_a]}_{self.cart_a[-1]}'
                )

                new_cart_a = list(self.cart_a[:-1])
                new_cart_a.pop(ind_a)
                eri_s.append(
                    NuclearPotential(new_cart_a, self.cart_b, self.bf_order))

                coef_s.append(
                    f'(-1) * 1/2 * 1/S1 * delta_{self.cart_a[ind_a]}_{self.cart_a[-1]}'
                )
                eri_s.append(
                    NuclearPotential(new_cart_a, self.cart_b,
                                     self.bf_order + 1))

            for ind_b in range(len(self.cart_b)):

                coef_s.append(
                    f'1/2 * 1/S1 * delta_{self.cart_b[ind_b]}_{self.cart_a[-1]}'
                )

                new_cart_b = list(self.cart_b)
                new_cart_b.pop(ind_b)
                eri_s.append(
                    NuclearPotential(self.cart_a[:-1], new_cart_b,
                                     self.bf_order))

                coef_s.append(
                    f'(-1) * 1/2 * 1/S1 * delta_{self.cart_b[ind_b]}_{self.cart_a[-1]}'
                )
                eri_s.append(
                    NuclearPotential(self.cart_a[:-1], new_cart_b,
                                     self.bf_order + 1))

        return coef_s, eri_s

    def apply_hrr_b_once(self):

        if self.Lb == 0:

            coef_s = ['1']
            eri_s = [NuclearPotential(self.cart_a, self.cart_b, self.bf_order)]

        else:

            # Eq.(A19), Obara-Saika JCP 1986

            coef_s = [f'PB_{self.cart_b[-1]}']
            eri_s = [
                NuclearPotential(self.cart_a, self.cart_b[:-1], self.bf_order)
            ]

            coef_s.append(f'(-1) * PC_{self.cart_b[-1]}')
            eri_s.append(
                NuclearPotential(self.cart_a, self.cart_b[:-1],
                                 self.bf_order + 1))

            for ind_b in range(len(self.cart_b) - 1):

                coef_s.append(
                    f'1/2 * 1/S1 * delta_{self.cart_b[ind_b]}_{self.cart_b[-1]}'
                )

                new_cart_b = list(self.cart_b[:-1])
                new_cart_b.pop(ind_b)
                eri_s.append(
                    NuclearPotential(self.cart_a, new_cart_b, self.bf_order))

                coef_s.append(
                    f'(-1) * 1/2 * 1/S1 * delta_{self.cart_b[ind_b]}_{self.cart_b[-1]}'
                )
                eri_s.append(
                    NuclearPotential(self.cart_a, new_cart_b,
                                     self.bf_order + 1))

            for ind_a in range(len(self.cart_a)):

                coef_s.append(
                    f'1/2 * 1/S1 * delta_{self.cart_a[ind_a]}_{self.cart_b[-1]}'
                )

                new_cart_a = list(self.cart_a)
                new_cart_a.pop(ind_a)
                eri_s.append(
                    NuclearPotential(new_cart_a, self.cart_b[:-1],
                                     self.bf_order))

                coef_s.append(
                    f'(-1) * 1/2 * 1/S1 * delta_{self.cart_a[ind_a]}_{self.cart_b[-1]}'
                )
                eri_s.append(
                    NuclearPotential(new_cart_a, self.cart_b[:-1],
                                     self.bf_order + 1))

        return coef_s, eri_s

    def apply_hrr_a(self):

        coefs, eris = self.apply_hrr_a_once()

        while True:

            done = True
            for e in eris:
                if e.La > 0:
                    done = False
                    break

            if done:
                break

            coefs, eris = apply_hrr_a_once(coefs, eris)

        return coefs, eris

    def apply_hrr_b(self):

        coefs, eris = self.apply_hrr_b_once()

        while True:

            done = True
            for e in eris:
                if e.Lb > 0:
                    done = False
                    break

            if done:
                break

            coefs, eris = apply_hrr_b_once(coefs, eris)

        return coefs, eris
