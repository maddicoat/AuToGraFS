import numpy as np

class ClusterBase:
    def get_layer_distance(self, miller, layers=1):
        """Returns the distance between planes defined by the given miller
        index.
        """
        n = self.miller_to_direction(miller)
        d1 = d2 = 0.0

        d = np.abs(np.sum(n * self.lattice_basis, axis=1))
        mask = np.greater(d, 1e-10)
        if mask.sum() > 0:
            d1 = np.min(d[mask])

        if len(self.atomic_basis) > 1:
            atomic_basis = np.dot(self.atomic_basis, self.lattice_basis)
            d = np.sum(n * atomic_basis, axis=1)
            s = np.sign(d)
            d = np.abs(d)
            mask = np.greater(d, 1e-10)
            if mask.sum() > 0:
                d2 = np.min(d[mask])
                s2 = s[mask][np.argmin(d[mask])]

        if d2 > 1e-10:
            if s2 < 0 and d1 - d2 > 1e-10:
                d2 = d1 - d2
            elif s2 < 0 and d2 - d1 > 1e-10:
                d2 = 2 * d1 - d2
            elif s2 > 0 and d2 - d1 > 1e-10:
                d2 = d2 - d1

            if np.abs(d1 - d2) < 1e-10:
                ld = np.array([d1])
            elif np.abs(d1 - 2 * d2) < 1e-10:
                ld = np.array([d2])
            else:
                assert d1 > d2, "Something is wrong with the layer distance."
                ld = np.array([d2, d1 - d2])
        else:
            ld = np.array([d1])

        if len(ld) > 1:
            if layers < 0:
                ld = np.array([-ld[1], -ld[0]])
                layers *= -1

            map = np.arange(layers - (layers % 1), dtype=int) % len(ld)
            r = ld[map].sum() + (layers % 1) * ld[np.abs(map[-1] - 1)]
        else:
            r = ld[0] * layers

        return r

    def miller_to_direction(self, miller, norm=True):
        """Returns the direction corresponding to a given Miller index."""
        d = np.dot(miller, self.resiproc_basis)
        if norm:
            d = d / np.linalg.norm(d)
        return d

