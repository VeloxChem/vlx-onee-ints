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

    text += """
        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

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

    max_bf_order = write_electric_field(f'{a}{b}', 'm', 'bf_order')

    text += f"""
        // J. Chem. Phys. 84, 3963-3974 (1986)

        double V_ij = 0.0;

        for (int c = 0; c < ndipoles; c++)
        {{
            const auto x_c   = dipoles_info[c + ndipoles * 0];
            const auto y_c   = dipoles_info[c + ndipoles * 1];
            const auto z_c   = dipoles_info[c + ndipoles * 2];
            const auto x_dip = dipoles_info[c + ndipoles * 3];
            const auto y_dip = dipoles_info[c + ndipoles * 4];
            const auto z_dip = dipoles_info[c + ndipoles * 5];

            const double dip[3] = {{x_dip, y_dip, z_dip}};

            const double PC[3] = {{(a_i * x_i + a_j * x_j) / (a_i + a_j) - x_c,
                                  (a_i * y_i + a_j * y_j) / (a_i + a_j) - y_c,
                                  (a_i * z_i + a_j * z_j) / (a_i + a_j) - z_c}};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double F{max_bf_order}_t[{max_bf_order+1}];

            onee::computeBoysFunction(F{max_bf_order}_t, (a_i + a_j) * r2_PC, {max_bf_order}, boys_func_table.data(), boys_func_ft.data());

            for (int m = 0; m < 3; m++)
            {{
                // Note: minus sign from electric field - electric dipole interaction

                V_ij += (-1.0) * dip[m] * (
"""

    text += write_electric_field(f'{a}{b}', 'm')

    text += f"""
                );
            }}
        }}

        efield[ij] = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j) * S_ij_00 * V_ij;
    }}

    for (int ij = 0; ij < {a}{b}_prim_pair_count; ij++)
    {{
        const auto i = std::get<0>(pair_inds_{a}{b}[ij]);
        const auto j = std::get<1>(pair_inds_{a}{b}[ij]);

"""

    if a == 's':
        text += f'        const auto i_cgto = {a}_prim_aoinds[i{denom[a]}];\n'
    else:
        text += f'        const auto i_cgto = {a}_prim_aoinds[(i{denom[a]}) + {a}_prim_count * (i{rem[a]})];\n'

    if b == 's':
        text += f'        const auto j_cgto = {b}_prim_aoinds[j{denom[b]}];\n'
    else:
        text += f'        const auto j_cgto = {b}_prim_aoinds[(j{denom[b]}) + {b}_prim_count * (j{rem[b]})];\n'

    text += f"""
        // Cartesian to spherical
"""

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

    if_str = 'if (i != j) ' if a == b else ''

    text += f"""
                auto coef_sph = i_coef_sph * j_coef_sph;

                V.row(i_cgto_sph)[j_cgto_sph] += efield[ij] * coef_sph;

                {if_str}V.row(j_cgto_sph)[i_cgto_sph] += efield[ij] * coef_sph;
            }}
        }}
    }}
"""

    print(text)
