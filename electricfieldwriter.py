from electricfield import ElectricField
from intsutils import apply_hrr_a, apply_hrr_b, simplify_coef
from intsutils import apply_hrr_a_once, apply_hrr_b_once
from intswriter import write_integrals


def write_electric_field(ab, mu, flag=''):

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

    efield_ab = ElectricField(comp_dict_1[ab[0]], comp_dict_2[ab[1]], 0, mu)

    flip_ab = False

    # apply HRR on a

    coefs, eris = efield_ab.apply_hrr_a()

    # apply HRR on b

    coefs, eris = apply_hrr_b(coefs, eris)

    # apply HRR on a and b once more, to make sure that all ElectricField (_A_mu)
    # integrals are converted to NuclearPotential (_A_00) integrals

    coefs, eris = apply_hrr_a_once(coefs, eris)
    coefs, eris = apply_hrr_b_once(coefs, eris)

    # simplify coefficients

    final_list = []
    for ind, (c, e) in enumerate(zip(coefs, eris)):
        final_list.append((simplify_coef(c), ind, e))

    coefs, eris = [], []
    for (c, ind, e) in sorted(final_list):
        coefs.append(c)
        eris.append(e)

    # find out maximum order of Boys function

    max_bf_order = write_integrals(coefs, eris, flag='bf_order')

    if flag == 'bf_order':
        return max_bf_order

    # reorganize terms

    final_ints_lines = write_integrals(coefs, eris, indent=20)

    # finalize the result

    final_ints_str = ''

    for line in final_ints_lines:

        if '( S S )^' in line:
            assert '_A_00 ( S S )^' in line
            bf_order = int(line.split('( S S )^')[1].split()[0])
            old_keyword = f'_A_00 ( S S )^{bf_order}'
            new_keyword = f'F{max_bf_order}_t[' + str(bf_order) + '] * ('
            line = line.replace(old_keyword, new_keyword)

        if line.endswith('__nonewline__'):
            final_ints_str += line.split('__nonewline__')[0]
        else:
            final_ints_str += line + '\n'

    final_ints_str = final_ints_str.replace(r'(1.0 * (-1.0))', '(-1.0)')
    final_ints_str = final_ints_str.replace(r'(1.0 * (-2.0))', '(-2.0)')
    final_ints_str = final_ints_str.replace(r'(1.0 * (-3.0))', '(-3.0)')

    return final_ints_str
