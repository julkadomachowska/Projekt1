# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 11:26:17 2022

@author: Julka
"""

import numpy as np
from geo_v1 import *
geo = Transformacje(model = 'wgs84')

plik = 'wsp_inp.txt'
#odczyt z pliku
tablica = np.genfromtxt(plik, delimiter = ',', skip_header = 4)

w, r = np.shape(tablica)

wynik = np.zeros((w,10))

i = 0

#tutaj trzeba zrobić pętlę i przeliczenie do wszystkich układóW 
for wiersz in tablica:
     fi, lam, h = geo.xyz2plh(wiersz[0], wiersz[1], wiersz[2])
     x, y, z = geo.flh2XYZ(fi, lam, h)
     xgk, ygk = geo.fl2xy(fi, lam)
     x2000, y2000 = geo.u2000(xgk, ygk)
     x1992, y1992 = geo.u1992(xgk, ygk)
     
     wynik[i, 0] = fi
     wynik[i, 1] = lam
     wynik[i, 2] = h
     wynik[i, 3] = x
     wynik[i, 4] = y
     wynik[i, 5] = z
     wynik[i, 6] = x2000
     wynik[i, 7] = y2000
     wynik[i, 8] = x1992
     wynik[i, 9] = y1992
     
     
     i+=1
print(wynik)
    
np.savetxt('wsp_out.txt', wynik, delimiter = ',', fmt = ['%10.7f', '%10.7f','%10.3f','%10.3f','%10.3f','%10.3f','%10.3f','%10.3f','%10.3f','%10.3f'], header = 'zamiana wspolrzednych geodezyjnych \\ Julka Domachowska')