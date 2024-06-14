from nuclearpotentialwriter import write_nuclear_potential
from electricfieldwriter import write_electric_field

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

    #pragma omp parallel for schedule(static, PAD_SIZE)
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

    text += f"""
        const auto i_atom = cart_ao_to_atom_ids[i_cgto];
        const auto j_atom = cart_ao_to_atom_ids[j_cgto];

        const double rij[3] = {{x_j - x_i, y_j - y_i, z_j - z_i}};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

"""

    if a == 'p' or a == 'd' or a == 'f':
        text += f'        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];\n'
    if a == 'd' or a == 'f':
        text += f'        const auto PA_1 = (a_j / (a_i + a_j)) * rij[a1];\n'
    if a == 'f':
        text += f'        const auto PA_2 = (a_j / (a_i + a_j)) * rij[a2];\n'
    text += '\n'

    if b == 'p' or b == 'd' or b == 'f':
        text += f'        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];\n'
    if b == 'd' or b == 'f':
        text += f'        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];\n'
    if b == 'f':
        text += f'        const auto PB_2 = (-a_i / (a_i + a_j)) * rij[b2];\n'
    text += '\n'

    max_bf_order = write_nuclear_potential(f'{a}{b}', 'bf_order')
    max_bf_order += 1

    text += f"""
        // J. Chem. Phys. 84, 3963-3974 (1986)

        for (int c = 0; c < npoints; c++)
        {{
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto q_c = points_info[c + npoints * 3];

            const double PC[3] = {{(a_i * x_i + a_j * x_j) / (a_i + a_j) - x_c,
                                  (a_i * y_i + a_j * y_j) / (a_i + a_j) - y_c,
                                  (a_i * z_i + a_j * z_j) / (a_i + a_j) - z_c}};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double F{max_bf_order}_t[{max_bf_order+1}];

            onee::computeBoysFunction(F{max_bf_order}_t, (a_i + a_j) * r2_PC, {max_bf_order}, boys_func_table.data(), boys_func_ft.data());

            for (int m = 0; m < 3; m++)
            {{
                const auto PA_m = (a_j / (a_i + a_j)) * rij[m];

                // Note: minus sign from electron charge

                double V_ij_for_i = (-1.0) * q_c * (
"""

    text += write_nuclear_potential(f'{a}{b}', 'gradient_a')

    text += f"""
                );

                // Note: minus sign from force (electric field) -> gradient conversion

                double V_ij_for_c = (-1.0) * q_c * (
"""

    text += write_electric_field(f'{a}{b}', 'm')

    text += f"""
                );

                double grad_i = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j) * S_ij_00 * V_ij_for_i;
                double grad_c = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j) * S_ij_00 * V_ij_for_c;
                double grad_j = -grad_c - grad_i;

"""

    if a != 's':
        text += f'                for (const auto& i_cgto_sph_ind_coef : cart_sph_{a}[i_cgto])\n'
    text += '                {\n'
    if a != 's':
        text += '                    auto i_cgto_sph = i_cgto_sph_ind_coef.first;\n'
        text += '                    auto i_coef_sph = i_cgto_sph_ind_coef.second;\n'
    else:
        text += '                    auto i_cgto_sph = i_cgto;\n'
        text += '                    double i_coef_sph = 1.0;\n'
    text += '\n'

    if b != 's':
        text += f'                    for (const auto& j_cgto_sph_ind_coef : cart_sph_{b}[j_cgto])\n'
    text += '                    {\n'
    if b != 's':
        text += '                        auto j_cgto_sph = j_cgto_sph_ind_coef.first;\n'
        text += '                        auto j_coef_sph = j_cgto_sph_ind_coef.second;\n'
    else:
        text += '                        auto j_cgto_sph = j_cgto;\n'
        text += '                        double j_coef_sph = 1.0;\n'

    if a == b:
        D_sym_factor = '((i == j) ? Dij : (Dij + Dji))'
    else:
        D_sym_factor = '(Dij + Dji)'

    text += f"""
                        auto coef_sph = i_coef_sph * j_coef_sph;

                        auto Dij = D[i_cgto_sph * naos + j_cgto_sph];
                        auto Dji = D[j_cgto_sph * naos + i_cgto_sph];

                        double D_sym = {D_sym_factor};

                        V_grad_omp[thread_id].row(i_atom)[m] += grad_i * coef_sph * D_sym;
                        V_grad_omp[thread_id].row(c)[m]      += grad_c * coef_sph * D_sym;
                        V_grad_omp[thread_id].row(j_atom)[m] += grad_j * coef_sph * D_sym;
                    }}
                }}
            }}
        }}
    }}
"""

    print(text)
