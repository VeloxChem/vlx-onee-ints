from intsutils import apply_hrr_a_once, apply_hrr_b_once
from electricdipole import ElectricDipole


class AngularMomentum:

    def __init__(self, cart_a, cart_b, mu):

        self.cart_a = list(cart_a)
        self.cart_b = list(cart_b)
        self.mu = mu

    def __repr__(self):

        label_a = self.get_angmom_str(self.cart_a)
        label_b = self.get_angmom_str(self.cart_b)

        return f'_AM_00 ( {label_a} {label_b} )'

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

    def apply_hrr_a(self):

        # rewrite AngularMomentum integrals in ElectricDipole integrals
        # using Eq.(7) and (5), Obara-Saika JCP 1986
        # also making use of antisymmetry of angular momentum

        m1 = f'{self.mu}1'
        m2 = f'{self.mu}2'

        coefs, eris = [], []

        coefs.append('2 * a_i')
        eris.append(ElectricDipole(self.cart_a + [m2], self.cart_b, m1))

        for ind_a in range(len(self.cart_a)):
            coefs.append(f'(-1) * delta_{m2}_{self.cart_a[ind_a]}')
            new_cart_a = list(self.cart_a)
            new_cart_a.pop(ind_a)
            eris.append(ElectricDipole(new_cart_a, self.cart_b, m1))

        coefs.append('(-1) * 2 * a_i')
        eris.append(ElectricDipole(self.cart_a + [m1], self.cart_b, m2))

        for ind_a in range(len(self.cart_a)):
            coefs.append(f'delta_{m1}_{self.cart_a[ind_a]}')
            new_cart_a = list(self.cart_a)
            new_cart_a.pop(ind_a)
            eris.append(ElectricDipole(new_cart_a, self.cart_b, m2))

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

        # rewrite AngularMomentum integrals in ElectricDipole integrals
        # using Eq.(7) and (5), Obara-Saika JCP 1986

        m1 = f'{self.mu}1'
        m2 = f'{self.mu}2'

        coefs, eris = [], []

        coefs.append('(-1) * 2 * a_j')
        eris.append(ElectricDipole(self.cart_a, self.cart_b + [m2], m1))

        for ind_b in range(len(self.cart_b)):
            coefs.append(f'delta_{m2}_{self.cart_b[ind_b]}')
            new_cart_b = list(self.cart_b)
            new_cart_b.pop(ind_b)
            eris.append(ElectricDipole(self.cart_a, new_cart_b, m1))

        coefs.append('2 * a_j')
        eris.append(ElectricDipole(self.cart_a, self.cart_b + [m1], m2))

        for ind_b in range(len(self.cart_b)):
            coefs.append(f'(-1) * delta_{m1}_{self.cart_b[ind_b]}')
            new_cart_b = list(self.cart_b)
            new_cart_b.pop(ind_b)
            eris.append(ElectricDipole(self.cart_a, new_cart_b, m2))

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
