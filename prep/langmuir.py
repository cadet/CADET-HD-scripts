#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt

ka = 1.144
kd = 2.0e-3
qm = 4.88

cin = 7.14e-3

max_c = 1e-1

def q(c): 
    value =  ka * c * qm / (kd + ka*c)
    return value


cs = np.linspace(0,max_c, 1000)
qs = [ q(c) for c in cs ]



with plt.style.context(['science']):
    fig, ax = plt.subplots()
    ax.plot(cs, qs)
    # ax.vlines(cin, 0, qm, ls='dashed', color='red', label='inlet conc.')
    # ax.hlines(qm, 0, max_c, ls='dashed', color='black', label='qmax')
    ax.set(xlabel='c', ylabel='q')
    # plt.show()
    plt.savefig('langmuir.pdf')



