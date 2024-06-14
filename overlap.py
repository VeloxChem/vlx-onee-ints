from intsutils import apply_hrr_a_once, apply_hrr_b_once


class Overlap:

    def __init__(self, cart_a, cart_b):

        self.cart_a = list(cart_a)
        self.cart_b = list(cart_b)

    def __repr__(self):

        label_a = self.get_angmom_str(self.cart_a)
        label_b = self.get_angmom_str(self.cart_b)

        return f'( {label_a} {label_b} )'

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

    def apply_gradient_a(self):

        # Eq.(5), Obara-Saika JCP 1986

        coef_s = ['2 * a_i']
        eri_s = [Overlap(self.cart_a + ['m'], self.cart_b)]

        for ind_a in range(len(self.cart_a)):

            coef_s.append(f'(-1) * delta_{self.cart_a[ind_a]}_m')

            new_cart_a = list(self.cart_a)
            new_cart_a.pop(ind_a)
            eri_s.append(Overlap(new_cart_a, self.cart_b))

        return coef_s, eri_s

    def apply_hrr_a_once(self):

        if self.La == 0:

            coef_s = ['1']
            eri_s = [Overlap(self.cart_a, self.cart_b)]

        else:

            # Eq.(A2), Obara-Saika JCP 1986

            coef_s = [f'PA_{self.cart_a[-1]}']
            eri_s = [Overlap(self.cart_a[:-1], self.cart_b)]

            for ind_a in range(len(self.cart_a) - 1):

                coef_s.append(
                    f'1/2 * 1/S1 * delta_{self.cart_a[ind_a]}_{self.cart_a[-1]}'
                )

                new_cart_a = list(self.cart_a[:-1])
                new_cart_a.pop(ind_a)
                eri_s.append(Overlap(new_cart_a, self.cart_b))

            for ind_b in range(len(self.cart_b)):

                coef_s.append(
                    f'1/2 * 1/S1 * delta_{self.cart_b[ind_b]}_{self.cart_a[-1]}'
                )

                new_cart_b = list(self.cart_b)
                new_cart_b.pop(ind_b)
                eri_s.append(Overlap(self.cart_a[:-1], new_cart_b))

        return coef_s, eri_s

    def apply_hrr_b_once(self):

        if self.Lb == 0:

            coef_s = ['1']
            eri_s = [Overlap(self.cart_a, self.cart_b)]

        else:

            # Eq.(A2), Obara-Saika JCP 1986

            coef_s = [f'PB_{self.cart_b[-1]}']
            eri_s = [Overlap(self.cart_a, self.cart_b[:-1])]

            for ind_b in range(len(self.cart_b) - 1):

                coef_s.append(
                    f'1/2 * 1/S1 * delta_{self.cart_b[ind_b]}_{self.cart_b[-1]}'
                )

                new_cart_b = list(self.cart_b[:-1])
                new_cart_b.pop(ind_b)
                eri_s.append(Overlap(self.cart_a, new_cart_b))

            for ind_a in range(len(self.cart_a)):

                coef_s.append(
                    f'1/2 * 1/S1 * delta_{self.cart_a[ind_a]}_{self.cart_b[-1]}'
                )

                new_cart_a = list(self.cart_a)
                new_cart_a.pop(ind_a)
                eri_s.append(Overlap(new_cart_a, self.cart_b[:-1]))

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
