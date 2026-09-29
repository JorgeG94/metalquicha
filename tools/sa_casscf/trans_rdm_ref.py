import numpy as np
from pyscf.fci import direct_spin1, cistring

NORB = 4
NCHOL = 3

h1e = np.zeros((NORB, NORB))
for q in range(NORB):
    for p in range(NORB):
        h1e[p, q] = -1.0/(p + q + 2)

b = np.zeros((NORB, NORB, NCHOL))
for l in range(NCHOL):
    for q in range(NORB):
        for p in range(NORB):
            b[p, q, l] = 1.0/(p + q + l + 3)

eri = np.einsum('pql,rsl->pqrs', b, b)

nelec = (2, 2)
na = cistring.num_strings(NORB, 2)
nb = cistring.num_strings(NORB, 2)

bra = np.zeros((na, nb))
ket = np.zeros((na, nb))
for ia in range(na):
    for ib in range(nb):
        bra[ia, ib] = np.sin(0.3*(ia + 1) + 0.7*(ib + 1))
        ket[ia, ib] = np.sin(0.5*(ia + 1) - 0.2*(ib + 1) + 0.9)

dm1, dm2 = direct_spin1.trans_rdm12(bra, ket, NORB, nelec)


def our_dm1(p, q):
    return dm1[q - 1, p - 1]


def our_dm2(p, q, r, s):
    return dm2[p - 1, q - 1, r - 1, s - 1]


for (p, q) in [(1, 1), (2, 2), (1, 2), (2, 1), (3, 4), (4, 1)]:
    print("dm1", p, q, repr(our_dm1(p, q)))

for (p, q, r, s) in [(1, 1, 1, 1), (2, 3, 4, 1), (1, 2, 2, 1),
                      (3, 1, 2, 4), (4, 4, 1, 1), (2, 1, 3, 4)]:
    print("dm2", p, q, r, s, repr(our_dm2(p, q, r, s)))
