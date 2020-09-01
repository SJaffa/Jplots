"""
Created on Thu Jul 27 17:09:25 2017

@author: sjaffa
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import momentsmod as mm

testlist = [#'diag-fil','diag-fil2',
    'disk', 'disk-cc',#'peak',
    #'noisy', 'dot',
    'fil-thin', 'fil-thick2',#'fil-thin2',
    'ellipse', 'ellipse-cc',#'disk-cc2',
    #'half-ring-thick','half-ring-thin',
    'ring-thick', 'ring-thin',
    'ring-lumpy', #'asym_cc','fil-thin_dot',
    ]

fig = plt.figure(1)
gridspec.GridSpec(5, 5)
fig.set_size_inches(w=10, h=10)

bigax = plt.subplot2grid((5, 5), (0, 0), colspan=4, rowspan=4)
c = 0
loclist = [(0,4), (1,4), (2,4), (3,4),
           (4,0), (4,1), (4,2), (4,3), (4,4)]

for n in testlist:

    g = mm.testdata2d(100, shape=n)
    smallax = plt.subplot2grid((5,5), loclist[c])
    im = smallax.imshow(np.sqrt(g.T), interpolation='None',
                        origin='lower', cmap='Greys', vmin=0)
    plt.setp(smallax.get_xticklabels(), visible=False)
    plt.setp(smallax.get_yticklabels(), visible=False)
    smallax.text(5, 85, "%i"%c, color='k', fontsize=14)

    
    J1,J2,com,t1 = mm.moments_2d(g)


    bigax.plot(J1, J2, marker='o', label=n, color='k')
    bigax.text(J1, J2+0.05, "%i"%c, color='k', fontsize=14)
    print()
    c += 1

bigax=mm.plot_moments_axes(bigax, text=True)

plt.savefig('proof.pdf')
