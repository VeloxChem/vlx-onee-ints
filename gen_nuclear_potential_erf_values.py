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

from nuclearpotentialwriter import write_nuclear_potential

multi = {'s': '', 'p': ' * 3', 'd': ' * 6', 'f': ' * 10'}
denom = {'s': '', 'p': ' / 3', 'd': ' / 6', 'f': ' / 10'}
rem = {'s': '', 'p': ' % 3', 'd': ' % 6', 'f': ' % 10'}

angmoms = 'spdf'

ab_list = []
for i in range(len(angmoms)):
    for j in range(i, len(angmoms)):
        ab_list.append(angmoms[i] + angmoms[j])

for ab in ab_list:
    a, b = ab

    text = f"""
    // {a.upper()}-{b.upper()} block

    #pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
    for (int ij = 0; ij < {a}{b}_prim_pair_count; ij++)
    {{
        const auto thread_id = omp_get_thread_num();

        const auto i = std::get<0>(pair_inds_{a}{b}[ij]);
        const auto j = std::get<1>(pair_inds_{a}{b}[ij]);

        const auto a_i = {a}_prim_info[i{denom[a]} + {a}_prim_count * 0];
        const auto c_i = {a}_prim_info[i{denom[a]} + {a}_prim_count * 1];
        const auto x_i = {a}_prim_info[i{denom[a]} + {a}_prim_count * 2];
        const auto y_i = {a}_prim_info[i{denom[a]} + {a}_prim_count * 3];
        const auto z_i = {a}_prim_info[i{denom[a]} + {a}_prim_count * 4];

        const auto a_j = {b}_prim_info[j{denom[b]} + {b}_prim_count * 0];
        const auto c_j = {b}_prim_info[j{denom[b]} + {b}_prim_count * 1];
        const auto x_j = {b}_prim_info[j{denom[b]} + {b}_prim_count * 2];
        const auto y_j = {b}_prim_info[j{denom[b]} + {b}_prim_count * 3];
        const auto z_j = {b}_prim_info[j{denom[b]} + {b}_prim_count * 4];

"""

    for x, m, n in zip([a, b], 'ab', 'ij'):
        if x == 'p':
            text += f'        const auto {m}0 = {n} % 3;\n'
        elif x == 'd':
            text += f'        const auto {m}0 = d_cart_inds[{n} % 6][0];\n'
            text += f'        const auto {m}1 = d_cart_inds[{n} % 6][1];\n'
        elif x == 'f':
            text += f'        const auto {m}0 = f_cart_inds[{n} % 10][0];\n'
            text += f'        const auto {m}1 = f_cart_inds[{n} % 10][1];\n'
            text += f'        const auto {m}2 = f_cart_inds[{n} % 10][2];\n'
        text += '\n'

    if a == 's':
        text += f'        const auto i_cgto = {a}_prim_aoinds[i{denom[a]}];\n'
    else:
        text += f'        const auto i_cgto = {a}_prim_aoinds[(i{denom[a]}) + {a}_prim_count * (i{rem[a]})];\n'

    if b == 's':
        text += f'        const auto j_cgto = {b}_prim_aoinds[j{denom[b]}];\n'
    else:
        text += f'        const auto j_cgto = {b}_prim_aoinds[(j{denom[b]}) + {b}_prim_count * (j{rem[b]})];\n'

    text += """
        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto inv_aij = 1.0 / (a_i + a_j);

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI * inv_aij, 1.5) * std::exp(-a_i * a_j * inv_aij * r2_ij);

"""

    if a == 'p' or a == 'd' or a == 'f':
        text += f'        const auto PA_0 = (a_j * inv_aij) * rij[a0];\n'
    if a == 'd' or a == 'f':
        text += f'        const auto PA_1 = (a_j * inv_aij) * rij[a1];\n'
    if a == 'f':
        text += f'        const auto PA_2 = (a_j * inv_aij) * rij[a2];\n'
    text += '\n'

    if b == 'p' or b == 'd' or b == 'f':
        text += f'        const auto PB_0 = (-a_i * inv_aij) * rij[b0];\n'
    if b == 'd' or b == 'f':
        text += f'        const auto PB_1 = (-a_i * inv_aij) * rij[b1];\n'
    if b == 'f':
        text += f'        const auto PB_2 = (-a_i * inv_aij) * rij[b2];\n'
    text += '\n'

    max_bf_order = write_nuclear_potential(f'{a}{b}', 'bf_order')

    text += f'        // product of density and cart-sph transformation coefficients\n'
    text += f'\n'
    text += f'        double dens_coef_prod = 0.0;\n'
    text += f'\n'

    if a != 's':
        text += f'        for (const auto& i_cgto_sph_ind_coef : cart_sph_{a}[i_cgto])\n'
    text += '        {\n'
    if a != 's':
        text += '            auto i_cgto_sph = i_cgto_sph_ind_coef.first;\n'
        text += '            auto i_coef_sph = i_cgto_sph_ind_coef.second;\n'
    else:
        text += '            auto i_cgto_sph = i_cgto;\n'
        text += '            double i_coef_sph = 1.0;\n'
    text += '\n'

    if b != 's':
        text += f'            for (const auto& j_cgto_sph_ind_coef : cart_sph_{b}[j_cgto])\n'
    text += '            {\n'
    if b != 's':
        text += '                auto j_cgto_sph = j_cgto_sph_ind_coef.first;\n'
        text += '                auto j_coef_sph = j_cgto_sph_ind_coef.second;\n'
    else:
        text += '                auto j_cgto_sph = j_cgto;\n'
        text += '                double j_coef_sph = 1.0;\n'

    if a == b:
        D_sym_factor = '((i == j) ? Dij : (Dij + Dji))'
    else:
        D_sym_factor = '(Dij + Dji)'

    text += f"""
                auto coef_sph = i_coef_sph * j_coef_sph;

                auto Dij = D[i_cgto_sph * naos + j_cgto_sph];
                auto Dji = D[j_cgto_sph * naos + i_cgto_sph];

                double D_sym = {D_sym_factor};

                dens_coef_prod += coef_sph * D_sym;
            }}
        }}

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto rho = a_i + a_j;

        const auto V_ij_00 = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j) * S_ij_00;

        for (int c = 0; c < npoints; c++)
        {{
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];

            const double PC[3] = {{(a_i * x_i + a_j * x_j) * inv_aij - x_c,
                                  (a_i * y_i + a_j * y_j) * inv_aij - y_c,
                                  (a_i * z_i + a_j * z_j) * inv_aij - z_c}};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double d2 = 1.0;

            if ((omega != nullptr) && (omega[c] != 0.0)) d2 = omega[c] * omega[c] / (rho + omega[c] * omega[c]);

            double F{max_bf_order}_t[{max_bf_order+1}];

            onee::computeBoysFunction(F{max_bf_order}_t, rho * d2 * r2_PC, {max_bf_order}, boys_func_table.data(), boys_func_ft.data());

            if ((omega != nullptr) && (omega[c] != 0.0))
            {{
                const double sqrt_d2 = std::sqrt(d2);

                F{max_bf_order}_t[0] *= sqrt_d2;
"""

    for omega_bf_order in range(1, max_bf_order + 1):
        d2_list = ['d2' for tmp in range(omega_bf_order)]
        d2_rhs = ' * '.join(d2_list)
        text += f"""
                F{max_bf_order}_t[{omega_bf_order}] *= {d2_rhs} * sqrt_d2;
"""

    text += f"""
            }}

            // Note: minus sign from electron charge

            double npot_val = (-1.0) * V_ij_00 * (
"""

    text += write_nuclear_potential(f'{a}{b}')

    text += f"""
            );

            npot_values_omp[thread_id][c] += npot_val * dens_coef_prod;
        }}
    }}
"""

    print(text)
