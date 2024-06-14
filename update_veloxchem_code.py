import subprocess


def get_output(command):

    try:
        output = subprocess.check_output(command)
    except subprocess.CalledProcessError:
        return ''
    return output.decode('utf-8')


pycodes_str = """
gen_angular_momentum.py
gen_linear_momentum.py
gen_electric_dipole.py
gen_electric_dipole_gradient.py
gen_electric_field.py
gen_electric_field_values.py
gen_nuclear_potential.py
gen_nuclear_potential_gradient.py
gen_kinetic_energy.py
gen_kinetic_energy_gradient.py
gen_overlap.py
gen_overlap_gradient.py
"""

references_str = """
reference_output/ref_angular_momentum.txt
reference_output/ref_linear_momentum.txt
reference_output/ref_electric_dipole.txt
reference_output/ref_electric_dipole_gradient.txt
reference_output/ref_electric_field.txt
reference_output/ref_electric_field_values.txt
reference_output/ref_nuclear_potential.txt
reference_output/ref_nuclear_potential_gradient.txt
reference_output/ref_kinetic_energy.txt
reference_output/ref_kinetic_energy_gradient.txt
reference_output/ref_overlap.txt
reference_output/ref_overlap_gradient.txt
"""

vlxcodes_str = """
/home/xinli/gitlab/VeloxChem.simd_master/src/onee_ints/AngularMomentumIntegrals.cpp
/home/xinli/gitlab/VeloxChem.simd_master/src/onee_ints/LinearMomentumIntegrals.cpp
/home/xinli/gitlab/VeloxChem.simd_master/src/onee_ints/ElectricDipoleIntegrals.cpp
/home/xinli/gitlab/VeloxChem.simd_master/src/onee_ints/ElectricDipoleGradient.cpp
/home/xinli/gitlab/VeloxChem.simd_master/src/onee_ints/ElectricFieldIntegrals.cpp
/home/xinli/gitlab/VeloxChem.simd_master/src/onee_ints/ElectricFieldValues.cpp
/home/xinli/gitlab/VeloxChem.simd_master/src/onee_ints/NuclearPotentialIntegrals.cpp
/home/xinli/gitlab/VeloxChem.simd_master/src/onee_ints/NuclearPotentialGradient.cpp
/home/xinli/gitlab/VeloxChem.simd_master/src/onee_ints/KineticEnergyIntegrals.cpp
/home/xinli/gitlab/VeloxChem.simd_master/src/onee_ints/KineticEnergyGradient.cpp
/home/xinli/gitlab/VeloxChem.simd_master/src/onee_ints/OverlapIntegrals.cpp
/home/xinli/gitlab/VeloxChem.simd_master/src/onee_ints/OverlapGradient.cpp
"""

pycodes = pycodes_str.split()
references = references_str.split()
vlxcodes = vlxcodes_str.split()

for code, ref, vlxcode in zip(pycodes, references, vlxcodes):
    output = get_output(['python3', code])
    with open(ref, 'w') as fh:
        fh.write(output)
    with open(vlxcode, 'r') as fh:
        vlxlines = fh.readlines()
    with open(vlxcode, 'w') as fh:
        auto_gen_flag = False
        code_written = False
        for line in vlxlines:
            if '// auto-generated code begins here' in line:
                fh.write(line)
                auto_gen_flag = True
            if '// auto-generated code ends here' in line:
                auto_gen_flag = False
            if not auto_gen_flag:
                fh.write(line)
            else:
                if not code_written:
                    fh.write(output)
                    code_written = True
