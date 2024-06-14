from kineticenergy import KineticEnergy
from intsutils import apply_hrr_a, apply_hrr_b, simplify_coef
from intswriter import write_integrals


def write_kinetic_energy(ab, flag=''):

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

    kin_ab = KineticEnergy(comp_dict_1[ab[0]], comp_dict_2[ab[1]])

    flip_ab = False

    # if requested, apply gradient on a
    # apply HRR on a

    if flag == 'gradient_a':
        coefs, eris = kin_ab.apply_gradient_a()
        coefs, eris = apply_hrr_a(coefs, eris)
    else:
        coefs, eris = kin_ab.apply_hrr_a()

    # apply HRR on b

    coefs, eris = apply_hrr_b(coefs, eris)

    # simplify coefficients

    final_list = []
    for ind, (c, e) in enumerate(zip(coefs, eris)):
        final_list.append((simplify_coef(c), ind, e))

    coefs, eris = [], []
    for (c, ind, e) in sorted(final_list):
        coefs.append(c)
        eris.append(e)

    # reorganize terms

    final_ints_lines = write_integrals(coefs, eris, indent=16)

    # finalize the result

    final_ints_str = ''

    for line in final_ints_lines:

        if '_K_00 ( S S )' in line:
            line = line.replace('_K_00 ( S S )', 'K_ij_00 * (')
        elif '( S S )' in line:
            # Note: lift ksi outside of the bracket
            line = line.replace('( S S )', '__ksi__ * (')

        if line.endswith('__nonewline__'):
            final_ints_str += line.split('__nonewline__')[0]
        else:
            # Note: lift ksi outside of the bracket
            line = line.replace(' * ksi', '')
            line = line.replace('ksi * ', '')
            line = line.replace('__ksi__', 'ksi')
            final_ints_str += line + '\n'

    final_ints_str = final_ints_str.replace(r'(1.0 * (-1.0))', '(-1.0)')
    final_ints_str = final_ints_str.replace(r'(1.0 * (-2.0))', '(-2.0)')
    final_ints_str = final_ints_str.replace(r'(1.0 * (-3.0))', '(-3.0)')

    return final_ints_str
