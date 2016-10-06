# -*- coding: utf-8 -*-
import numpy as np
import pylab as pl
import matplotlib as mp

# parametry
n=10 # liczba modulow(sprezynek)
k=80 # stala sprezystosci
m=0.01 # masa jednego modulu(sprezynki)
g=9.81 # przyspieszenie ziemskie
d=0.1 # dlugosc swobodna jednego modulu(sprezynki)
dt = 0.0002 # krok czasowy
t=8 # czas symulacji
e=0 # wspolczynnik elastycznosci; Jezeli e=0 to mamy do czynienia ze zderzeniem calkowicie nieelastycznym; jezeli e=1 to
#mamy do czynienia ze zderzeniem w pelni elastycznym


"""
symulacja - funkcja symulujaca ruch sprężyny 
Input : parametry n, k, m, g, d
Output : tablica zawierajaca dane o ruchu spadajacej sprezyny w czasie, zaleznosc przyspieszenia srodka masy od czasu, 
zaleznosc predkosci srodka masy od czasu, zaleznosc polozenia srodka masy od czasu, zaleznosc dlugosci sprezyny od czasu.
tablica jest postaci [0-polozenie, 1-predkosc, 2-przyspieszenie][numer sprezynki 0, ..., n-1][czas] """
def symulacja(n, k, m, g, d, e):
    ## I Stan poczatkowy(sprezyna jest zawieszona nieruchomo)
    
    xl=[] # odleglosci pomiedzy sprezynkami
    sl=[] # polozenia
    
    #obliczam odleglosci miedzy modulami(sprezynkami)
    for i in range(1,n+1):  #przechodzi po 1, 2, ..., n
        xl.append( m*g*( n-i+1 )/k + d)
        
    #obliczam polozenia modulow(sprezynek)
    odleglosc=0 # zmienna pomocnicza
    for i in range(1,n+1): #przechodzi po 1, 2, ..., n
        odleglosc+=xl[i-1]
        sl.append(odleglosc)
    
    ## II Ruch ukladu(sprezyna spada swobodnie)
    
    tablica = np.zeros((3,n,int(t/dt))) # tablica[zmienna 0-sl, 1-vl , 2-al][numer sprezynki 0, ..., n-1][czas]
    
    # wypelniam tabelice wartosciami polozen moduow dla chwili 0 
    tablica[0][0][0]=d # sprezynka polozona najwyzej ma dlugosc d (nie ma sily do gory)
    for i in range(1,n):
        tablica[0][i][0]=tablica[0][i-1][0]+xl[i]
        
    # wypelniam tablice wartosciami przyspieszen modulow w chwili 0
    # przyspieszenie dla l = 1
    tablica[2][0][0]=k*(tablica[0][1][0]-tablica[0][0][0]-d)/m + g  
    # przyspieszenie dla l = 1, 2, ..., n-1
    for i in range(1,n-1):
        tablica[2][i][0]=k*(tablica[0][i+1][0]-2*tablica[0][i][0]+tablica[0][i-1][0])/m + g    
    # przyspieszenie dla l = n
    tablica[2][n-1][0]=-k*(tablica[0][n-1][0]-tablica[0][n-2][0]-d)/m + g
    
    u=False #zmienna pomocnicza do sprawdzenia warunku czy dlugosc sprezynki nie jest ujemna

    #wypelniam tablice wartosciami predkosci, polozen i przyspieszen w czasie od 0 do t z krokiem dt
    liczba_zderzen=0 #liczba zderzen pomiedzy modulami
    for t_ruch in range(1,int(t/dt)):
        for i in range(n):
            # predkosci dla l = 1, ..., n
            tablica[1][i][t_ruch] = tablica[2][i][t_ruch-1] * dt + tablica[1][i][t_ruch-1]
            
        # polozenia dla l = 1
        tablica[0][0][t_ruch] = (tablica[1][0][t_ruch-1] + tablica[2][0][t_ruch-1] * dt / 2) * dt + tablica[0][0][t_ruch-1]
        
        for i in range(1,n):    
            # polozenia dla l = 2, ..., n
            tablica[0][i][t_ruch] = (tablica[1][i][t_ruch-1] + tablica[2][i][t_ruch-1] * dt / 2) * dt + tablica[0][i][t_ruch-1]
            
            # sprawdzam warunek czy dlugosc sprezynki nie jest ujemna
            # obliczam predkosci modulow po zdezeniu za pomoca funkcji predkosc_zderzenie
            if tablica[0][i][t_ruch]<=tablica[0][i-1][t_ruch]:
                tablica[1][i][t_ruch]=(tablica[1][i][t_ruch]+tablica[1][i-1][t_ruch])/2
                tablica[1][i-1][t_ruch]=tablica[1][i][t_ruch]
                tablica[0][i][t_ruch]=tablica[0][i-1][t_ruch] # Wprowadzam zmiane w polozeniu na niekorzysc ukladu w dol(dla sprezynki o nizszym polozeniu)  
                # fala idzie szybciej do dolu
                liczba_zderzen=liczba_zderzen+1
                u=True
            
        # obliczam przyspieszenie dla l=1
        tablica[2][0][t_ruch] = k * (tablica[0][1][t_ruch-1]-tablica[0][0][t_ruch-1]-d)/m + g
        
        # obliczam przyspieszenie dla l = 2, ..., n-1
        for i in range(1,n-1):
            tablica[2][i][t_ruch] = k * (tablica[0][i+1][t_ruch-1] - 2*tablica[0][i][t_ruch-1]+tablica[0][i-1][t_ruch-1])/m + g
            
        # obliczam przyspieszenie dla l = n
        tablica[2][n-1][t_ruch] = - k * (tablica[0][n-1][t_ruch-1]-tablica[0][n-2][t_ruch-1]-d)/m + g
    
    # wyswietla komunikat o ujemnej dlugosci sprezynki
    if u:
        print "Uwaga ujemna dlugosc sprezynki"
        print "dla k i d: " + str(k) + " " + str(d)
        print "liczba zderzen wynosi: " + str(liczba_zderzen) 
        
    return tablica
    
"""
predkosc_zderzenie - funkcja obiczajaca predkosci sprezynek po zderzeniu(kiedy zajdzie przypadek xl<=0)
Input : predkosci modulow przed zderzeniem
Output : predkosci modulow po zdrzeniu
"""
def predkosc_zderzenie(v1,v2,e):
    v2k=(v1+v2)/2
    v1k=v2k
    
    #v2k=0.5*v1*(1+e) + 0.5*v2*(1-e)
    #v1k=0.5*v1*(1-e) + 0.5*v2*(1+e)
    return (v1k,v2k)

"""
przysp_sr_masy - funkcja obliczajaca przyspieszenie srodka masy
Input : tablica[0-polozenie, 1-predkosc, 2-przyspieszenie][numer sprezynki 0, ..., n-1][czas]
Output : przyspieszenie srodka masy w kolejnych krokach czasowych od 0 do t
        """
def  przysp_sr_masy(tablica):
    #obliczam przyspieszenie srodka masy sprezyny
    przysp_sr_masy=[] # tablica zawierajaca wartosci przyspieszen w kolejnych krokach czasowych
    for i in range(int(t/dt)):
        suma=0
        for j in range(n):
            suma+=tablica[2][j][i]
        przysp_sr_masy.append(suma/n)
    return przysp_sr_masy
    
    
"""
predk_sr_masy - funkcja predkosc srodka masy
Input : tablica[0-polozenie, 1-predkosc, 2-przyspieszenie][numer sprezynki 0, ..., n-1][czas]
Output : predkosc w kolejnych krokach czasowych od 0 do t            
"""
def predk_sr_masy(tablica):
    #obliczam predkosc srodka masy
    predk_sr_masy=[] # tablica zawierajaca wartosci predkosci w kolejnych krokach czasowych
    for i in range(int(t/dt)):
        suma=0
        for j in range(n):
            suma+=tablica[1][j][i]
        predk_sr_masy.append(suma/n)
    return predk_sr_masy
        
        
"""
pol_sr_masy - funkcja zwracajaca polozenie srodka masy
Input : tablica[0-polozenie, 1-predkosc, 2-przyspieszenie][numer sprezynki 0, ..., n-1][czas]
Output : polozenie w kolejnych krokach czasowych od 0 do t 
"""
def pol_sr_masy(tablica):
    #obliczam polozenie srodka masy
    polozenie_sr_masy=[] # tablica zawierajaca wartosci polozen w kolejnych krokach czasowych
    czas=[] # tablica pomocicza krokow czasowych
    for i in range(int(t/dt)):
        czas.append(i)
        suma=0
        for j in range(n):
            suma+=tablica[0][j][i]
        polozenie_sr_masy.append(suma/n)
    return polozenie_sr_masy
        
        
        
"""
pol_najniz funkcja obliczajaca polozenie najnizszego modulu
Input : tablica[0-polozenie, 1-predkosc, 2-przyspieszenie][numer sprezynki 0, ..., n-1][czas]
Output : polozenie najnizszego modulu w kolejnych krokach czasowych od 0 do t 
"""
def pol_najniz(tablica):
    #obliczam polozenie modulu(sprezynki) najnizszego
    polozenie_ostatniego=[] # tablica zawierajaca wartosci polozen najnizszego modulu w kolejnych krokach czasowych
    for i in range(int(t/dt)):
        polozenie_ostatniego.append(tablica[0][n-1][i])
    return polozenie_ostatniego
        
        
    
"""dl_sprezyny - funkcja wyznaczajaca polozenie najnizszego modulu i dlugosc sprezyny
Input : tablica[0-polozenie, 1-predkosc, 2-przyspieszenie][numer sprezynki 0, ..., n-1][czas]
Output : dlugosc sprezyny w kolejnych krokach czasowych od 0 do t
"""
def dl_sprezyny(tablica):
    #obliczam dlugosc sprezyny
    dlugosc_sprezyny=[] # tablica wartowsci dlugosci sprezyny w kolejnych krokach czasowych
    for i in range(int(t/dt)):
            dlugosc_sprezyny.append(tablica[0][n-1][i]-tablica[0][0][i])
    return dlugosc_sprezyny
    
    
mp.rc('font', family='DejaVu Sans')

"""
tablica=symulacja(n,k,m,g,d,e)


#Wykres zaleznosci predkosci srodka masy od czasu
fig = pl.figure()
pl.plot(predk_sr_masy(tablica))
pl.xlabel('t [ms]')
pl.ylabel('predkosc srodka masy [m/s]')
fig.savefig('/home/pawel/Documents/licencjat/symulacja_kod/predk_sr.jpg')
pl.show()
pl.clf()


#Wykres zaleznosci polozenia srodka masy od czasu 
fig = pl.figure()
pl.plot(pol_sr_masy(tablica))
pl.xlabel('t [ms]')
pl.ylabel('polozenie srodka masy [m]')
pl.xlim([0,1000])
fig.savefig('/home/pawel/Documents/licencjat/symulacja_kod/polozenie_sr.jpg')
pl.show()
pl.clf()

"""

#Wykres zaleznosci przyspieszenia srodka masy od czasu
tablica=symulacja(n,k,m,g,d,e)
fig = pl.figure()
pl.plot(przysp_sr_masy(tablica))
pl.xlabel('t [ms]')
pl.ylabel('przyspieszenie srodka masy [m/s2]')
fig.savefig('/home/pawel/Documents/licencjat/symulacja_kod/przyspiesz_sr.jpg')
pl.show()
pl.clf()


"""

# Wykres zaleznosci polozenia dolnego modulu i dlugosci sprezyny od czasu
fig = pl.figure()
pl.plot(dl_sprezyny(tablica))
pl.xlabel('t [ms]')
pl.ylabel('polozenie dolnej sprezynki(zielony) i dlugosc sprezyny(nieb)')
pl.plot(pol_najniz(tablica))
fig.savefig('/home/pawel/Documents/licencjat/symulacja_kod/polozenie_ostat.jpg')
pl.clf()


# Wykresy zaleznosci dlugosci sprezyny od czasu dla k = 50, 70, 90
fig = pl.figure()
k=50
tablica=symulacja(n,k,m,g,d)

pl.plot(dl_sprezyny(tablica))
pl.xlabel('t [ms]')
pl.ylabel('dlugosc sprezyny[m], niebieski k=50, zielony k=70, czerwony k=90')

k=70
tablica=symulacja(n,k,m,g,d)


pl.plot(dl_sprezyny(tablica))

k=90
tablica=symulacja(n,k,m,g,d)

pl.plot(dl_sprezyny(tablica))
fig.savefig('/home/pawel/Documents/licencjat/symulacja_kod/dlugosci.jpg')
pl.clf()
pl.show()

#Wykresy predkosci ostatniego ciezarka dla 5 różnych k
fig = pl.figure()
k=65
pl.xlabel(u"t [ms]")
pl.ylabel(u"v(t) dla dolnego ciężarka [m/s]")
#pl.xlim([0,1000])
tablica=symulacja(n,k,m,g,d,e)

pl.plot(tablica[1][n-1][:], label="k=65")

k=80
tablica=symulacja(n,k,m,g,d,e)

pl.plot(tablica[1][n-1][:], label="k=80")

k=95
tablica=symulacja(n,k,m,g,d,e)

pl.plot(tablica[1][n-1][:], label="k=95")

k=110
tablica=symulacja(n,k,m,g,d,e)

pl.plot(tablica[1][n-1][:], label="k=110")

k=135
tablica=symulacja(n,k,m,g,d,e)

pl.plot(tablica[1][n-1][:], label="k=135")

pl.legend(bbox_to_anchor=(1.1, 0.97), loc=1, borderaxespad=0.)

fig.savefig('/home/pawel/Documents/licencjat/symulacja_kod/predkosci_ostatniego_k.jpg')
pl.show()
pl.clf()
"""
"""
#Znajdowanie czasu charakterystycznego sprezyny
tablica=symulacja(n,k,m,g,d,e)
t_ruch=0
while(tablica[1][n-1][t_ruch]<0.1):
    t_ruch+=1
print "#Znajdowanie czasu charakterystycznego sprezyny"
print "dla k i d: " + str(k) + " " + str(d)
print "czas po ktorym dolny spodek ruszyl"
print t_ruch*dt
"""

"""
#Wykresy predkosci ostatniego ciezarka dla 5 roznych d
print "#Wykresy predkosci ostatniego ciezarka dla 5 roznych d"
fig = pl.figure()
k=150
d=0.001

tablica=symulacja(n,k,m,g,d,e)
pl.xlabel(u"t [ms]")
pl.ylabel(u"v(t) dla dolnego ciężarka [m/s]")

pl.plot(tablica[1][n-1][:], label="d=0.23[m]")

d=0.005
tablica=symulacja(n,k,m,g,d,e)
pl.plot(tablica[1][n-1][:], label="d=0.25[m]")

d=0.025
tablica=symulacja(n,k,m,g,d,e)
pl.plot(tablica[1][n-1][:], label="d=0.27[m]")

d=0.125
tablica=symulacja(n,k,m,g,d,e)
pl.plot(tablica[1][n-1][:], label="d=0.31[m]")

d=0.625
tablica=symulacja(n,k,m,g,d,e)
pl.plot(tablica[1][n-1][:], label="d=0.33[m]")


pl.legend(bbox_to_anchor=(1.1, 0.97), loc=1, borderaxespad=0.)

pl.plot()

fig.savefig('/home/pawel/Documents/licencjat/symulacja_kod/predkosci_ostatniego_d.jpg')
pl.show()
pl.clf()
"""


"""
#Wykres zaleznosci dlugosci pierwszego modulu od czasu
fig = pl.figure()
k=80
d=0.12

tablica=symulacja(n,k,m,g,d,e)
pl.plot(tablica[1][n-1][:])
pl.show()
"""

#print predkosc_zderzenie(2,1,e)
