#!/usr/bin/python

import matplotlib.pyplot as plt
import sys

fin = open(sys.argv[1], 'r')

lines = [line.split() for line in fin]
i = []
real = []
imag = []

for u,v,w in lines:
    i.append(float(u))
    real.append(float(v))
    imag.append(float(w))

plt.plot(i, real, '.-', label="real")
plt.plot(i, imag, '.-', label="imag")
plt.legend()
plt.show()
