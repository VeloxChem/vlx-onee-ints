import subprocess


def get_output(command):

    try:
        output = subprocess.check_output(command)
    except subprocess.CalledProcessError:
        return ''
    return output.decode('utf-8')


pycodes_str = """
gen_angular_momentum.py
gen_electric_dipole_gradient.py
gen_electric_dipole.py
gen_electric_field.py
gen_electric_field_values.py
gen_kinetic_energy_gradient.py
gen_kinetic_energy.py
gen_linear_momentum.py
gen_nuclear_potential_gradient.py
gen_nuclear_potential.py
gen_overlap_gradient.py
gen_overlap.py
"""

references_str = """
reference_output/ref_angular_momentum.txt
reference_output/ref_electric_dipole_gradient.txt
reference_output/ref_electric_dipole.txt
reference_output/ref_electric_field.txt
reference_output/ref_electric_field_values.txt
reference_output/ref_kinetic_energy_gradient.txt
reference_output/ref_kinetic_energy.txt
reference_output/ref_linear_momentum.txt
reference_output/ref_nuclear_potential_gradient.txt
reference_output/ref_nuclear_potential.txt
reference_output/ref_overlap_gradient.txt
reference_output/ref_overlap.txt
"""

pycodes = pycodes_str.split()
references = references_str.split()

for code, ref in zip(pycodes, references):
    output = get_output(['python3', code])
    with open(ref, 'r') as fh:
        ref_output = fh.readlines()
    ref_output = ''.join(ref_output)
    assert output == ref_output
