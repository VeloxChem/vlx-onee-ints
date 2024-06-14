from intsutils import apply_hrr_a_once, apply_hrr_b_once
from overlap import Overlap


class LinearMomentum:

    def __init__(self, cart_a, cart_b, mu):

        self.cart_a = list(cart_a)
        self.cart_b = list(cart_b)
        self.mu = mu

    def __repr__(self):

        label_a = self.get_angmom_str(self.cart_a)
        label_b = self.get_angmom_str(self.cart_b)

        return f'_LM_00_{self.mu} ( {label_a} {label_b} )'

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

        # Eq.(7) and (5), Obara-Saika JCP 1986
        # also making use of antisymmetry of linear momentum

        coefs = ['2 * a_i']
        eris = [Overlap(self.cart_a + [self.mu], self.cart_b)]

        for ind_a in range(len(self.cart_a)):
            coefs.append(f'(-1) * delta_{self.mu}_{self.cart_a[ind_a]}')
            new_cart_a = list(self.cart_a)
            new_cart_a.pop(ind_a)
            eris.append(Overlap(new_cart_a, self.cart_b))

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

        # Eq.(7) and (5), Obara-Saika JCP 1986

        coefs = ['(-1) * 2 * a_j']
        eris = [Overlap(self.cart_a, self.cart_b + [self.mu])]

        for ind_b in range(len(self.cart_b)):
            coefs.append(f'delta_{self.mu}_{self.cart_b[ind_b]}')
            new_cart_b = list(self.cart_b)
            new_cart_b.pop(ind_b)
            eris.append(Overlap(self.cart_a, new_cart_b))

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
