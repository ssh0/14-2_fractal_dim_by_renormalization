#! /usr/bin/env python
# -*- coding:utf-8 -*-
#
# written by Shotaro Fujimoto, July 2014.

from Tkinter import *
import numpy as np
import matplotlib.pyplot as plt


class Percolation:

    def __init__(self, L=16, p=0.5927):
        self.L = L  # lattice size
        self.p = p
        self.sub = None
        self.lattice = np.zeros([self.L, self.L], dtype=bool)

    def percolate(self):
        if self.sub is None or not self.sub.winfo_exists():
            lattice = self.lattice
            rn = np.random.random([self.L, self.L])
            lattice[rn < p] = True
            lattice[rn >= p] = False
            self.lattice = lattice

    def labeling(self):
        label = np.zeros([self.L + 2, self.L + 2], dtype=int)
        n = 1
        r = range(1, self.L + 1)
        for i in r:
            for j in r:
                if self.lattice[i - 1][j - 1]:
                    nn = []
                    if label[i - 1][j] > 0:
                        nn.append(label[i - 1][j])
                    if label[i][j - 1] > 0:
                        nn.append(label[i][j - 1])
                    if len(nn) > 0:
                        label[i][j] = min(nn)
                    else:
                        label[i][j] = n
                        n += 1
        tag = range(1, n + 1)

        for i in reversed(r):
            for j in reversed(r):
                if label[i][j] > 0:
                    nn = []
                    if label[i + 1][j] > 0:
                        nn.append(label[i + 1][j])
                    if label[i][j + 1] > 0:
                        nn.append(label[i][j + 1])
                    nn.append(label[i][j])
                    min_tag = min(nn)
                    nn = set([x for x in nn if x != min_tag])
                    for t in nn:
                        tag[t - 1] = min_tag
                        label[label == t] = tag[t - 1]

        self.lattice = label[1:-1, 1:-1]
        left = set(self.lattice[0])
        right = set(self.lattice[self.L - 1])
        top = set([self.lattice[t][0] for t in range(self.L)])
        bottom = set([self.lattice[t][self.L - 1] for t in range(self.L)])
        self.ptag = (left.intersection(right)
                     | top.intersection(bottom)) - set([0])
        if len(self.ptag) == 0:
            self.percolate()
            self.labeling()
        return self.lattice, self.ptag

    def renormalization(self, b=2):
        if self.L % b != 0:
            raise ValueError("lattice cannot be divided by scale factor b")

        lattice = np.zeros([self.L, self.L])
        lattice[self.lattice == list(self.ptag)[0]] = 1
        rlattice = np.zeros([self.L / b, self.L / b])
        for i in range(self.L / 2):
            ic = 2 * i
            for j in range(self.L / 2):
                jc = 2 * j
                if lattice[ic, jc] * lattice[ic, jc + 1] == 1 or \
                        lattice[ic + 1, jc] * lattice[ic + 1, jc + 1] == 1:
                    rlattice[i, j] = 1

        self.rlattice = rlattice * list(self.ptag)[0]

        M = np.sum(lattice)
        rM = np.sum(rlattice)
        M_2 = M * M
        rM_2 = rM * rM
        return M_2, rM_2

    def draw_canvas(self, rect, L):
        default_size = 640  # default size of canvas
        r = int(default_size / (2 * L))
        fig_size = 2 * r * L
        margin = 10
        sub = Tk()

        sub.title('figure  ' + '(p=%s)' % str(self.p))
        self.canvas = Canvas(sub, width=fig_size + 2 * margin,
                             height=fig_size + 2 * margin)
        self.canvas.create_rectangle(margin, margin,
                                     fig_size + margin, fig_size + margin,
                                     outline='black', fill='white')
        self.canvas.pack()
        c = self.canvas.create_rectangle
        colors = ['blue', 'green', 'red', 'purple']
        colordict = dict(zip(list(self.ptag),
                             colors * (int(len(self.ptag) / len(colors)) + 1)
                             )
                         )

        nonzero_rect = np.nonzero(rect)
        for m, n in zip(nonzero_rect[0], nonzero_rect[1]):
            if rect[m][n] in self.ptag:
                c(2 * m * r + margin, 2 * n * r + margin,
                  2 * (m + 1) * r + margin, 2 * (n + 1) * r + margin,
                  outline='', fill=colordict[rect[m][n]])
            else:
                c(2 * m * r + margin, 2 * n * r + margin,
                  2 * (m + 1) * r + margin, 2 * (n + 1) * r + margin,
                  outline='', fill='black')

        sub.mainloop()

if __name__ == '__main__':
    L = 16
    p = 0.5927
    per = Percolation(L, p)
    b = 2
    trial = 10

    M_2, rM_2 = [], []
    for t in range(trial):
        per.percolate()
        per.labeling()
        m_2, rm_2 = per.renormalization(b)
        M_2.append(m_2)
        rM_2.append(rm_2)

#    per.draw_canvas(per.lattice, L)
#    per.draw_canvas(per.rlattice, L/b)
    ave_M_2 = np.average(M_2)
    ave_rM_2 = np.average(rM_2)
    D = np.log(ave_M_2 / ave_rM_2) / (2 * np.log(b))

    print 'D = %f' % D
