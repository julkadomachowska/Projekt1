from math import sin, cos, sqrt, atan, atan2, degrees, radians, pi, tan, acos
from numpy import deg2rad, matrix, sqrt, transpose

class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2


#Wspolrzedne geocentryczne XYZ --> phi, lam, h
   
    def xyz2plh(self, X, Y, Z, output = 'dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        lat
            [stopnie dziesiętne] - szerokość geodezyjna
        lon
            [stopnie dziesiętne] - długośc geodezyjna.
        h : TYPE
            [metry] - wysokość elipsoidalna
        output [STR] - optional, defoulf 
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """
        r   = sqrt(X**2 + Y**2)           # promień
        lat_prev = atan(Z / (r * (1 - self.ecc2)))    # pierwsze przybliilizenie
        lat = 0
        while abs(lat_prev - lat) > 0.000001/206265:    
            lat_prev = lat
            N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
            h = r / cos(lat_prev) - N
            lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lon = atan(Y/X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2);
        h = r / cos(lat) - N       
        if output == "dec_degree":
            return degrees(lat), degrees(lon), h 
        elif output == "dms":
            lat = self.deg2dms(degrees(lat))
            lon = self.deg2dms(degrees(lon))
            return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - output format not defined")
            

#Definicja N
    
    def Np(self, f):
        N = self.a / sqrt(1-self.ecc2*(sin(f)**2))
        return(N)
     
# Phi, lam, h --> XYZ
        
    def flh2XYZ(self, f,l,h):
        N = self.Np(f)
        X = (N+h)*cos(f)*cos(l)
        Y = (N+h)*cos(f)*sin(l)
        Z = (N*(1-self.ecc2)+h)*sin(f)
    
        return(X,Y,Z)

#definicja sigmy

    def sigma(self, f):

        A0 = 1-(self.ecc2/4)-(3/64)*(self.ecc2**2)-(5/256)*(self.ecc2**3);
        A2 = (3/8)*(self.ecc2 + (self.ecc2**2)/4 + (15/128)*(self.ecc2**3));
        A4 = (15/256)*(self.ecc2**2 + 3/4*(self.ecc2**3));
        A6 = (35/3072)*self.ecc2**3;
        si = self.a*(A0*f - A2*sin(2*f) + A4*sin(4*f) - A6*sin(6*f));
    
        return(si)

#odcechowanie wspolrzednych

    def fl2xy(self,f,l,L0=19):

        b2 = (self.a**2)*(1-self.ecc2)
        ep2 = (self.a**2-b2)/b2
        t = tan(f)
        n2 = ep2*(cos(f)**2)
        N = self.Np(f)
        si = self.sigma(f)
        dL = l - L0
    
        xgk = si + (dL**2/2)*N*sin(f)*cos(f)*(1 + (dL**2/12)*cos(f)**2*(5 - t**2 + 9*n2 + 4*n2**2) + (dL**4/360)*cos(f)**4*(61 - 58*t**2 + t**4 + 14*n2 - 58*n2*t**2))
        ygk = dL*N*cos(f)*(1 + (dL**2/6)*cos(f)**2*(1 - t**2 + n2) + (dL**4/120)*cos(f)**4*(5 - 18*t**2 + t**4 + 14*n2 - 58*n2*t**2))
    
        return(xgk,ygk)

#zamiana na wspolrzedne 2000

    def u2000(self,xgk,ygk,L0=21):

        m2000 = 0.999923;
    
        x = xgk * m2000;
        y = ygk * m2000 + L0*(180/pi/3)* 1000000 + 500000;
    
       # print("2000: x=",format(x, '10.3f'), "y= ",format(y, '10.3f'))
    
    
        return(x,y)

#zamiana na wspolrzedne 1992

    def u1992(self,xgk,ygk,L0=19):

        m92 = 0.9993;
    
        x = xgk * m92 - 5300000;
        y = ygk * m92 + 500000;

        #print("1992: x=",format(x, '10.3f'), "y= ",format(y, '10.3f'))
    
    
        return(x,y)
#zamiana na NEU
    
    def s_A_z2neu (self,s, A, z):
        A = deg2rad(A)
        z = deg2rad(z)
        n = s*sin(z)*cos(A)
        e = s*sin(z)*sin(A)
        u = s*cos(z)

        return (n, e, u)
    
# Wyznaczanie kata azymutu i odleglosci 2D
    def azytmut2D(self, fa, la, fb, lb, a, e2):
    
        b = self.a*sqrt(1-self.ecc2);
        f = 1 - self.b/self.a;
    
        dl = lb - la;
        Ln = dl;
    
        U_a = atan((1-f)*tan(fa));
        U_b = atan((1-f)*tan(fb));
    
        while 1:
            L = Ln;
        
            ssigma = sqrt((cos(U_b) * sin(L))**2 + (cos(U_a) * sin(U_b) - sin(U_a) * cos(U_b) * cos(L))**2)
            csigma = (sin(U_a) * sin(U_b)) + (cos(U_a) * cos(U_b) * cos(L))
            sigma = atan(ssigma/csigma)
    
            salpha = (cos(U_a) * cos(U_b) * sin(L))/ssigma
            c2alpha = 1 - (salpha)**2
        
            c2sigmam = csigma - ((2*sin(U_a)*sin(U_b))/ c2alpha)
        
            C = ((f/16)*c2alpha)*(4 + f*(4 - 3*c2alpha))
        
            Ln = dl + (1 - C) * f * salpha * (sigma + C*ssigma*(c2sigmam + C * csigma*(-1 + 2*(c2sigmam)**2)))
        
            if abs(Ln - L) < 0.000001/206265:
                break


    
        u2 = ((a**2 - b**2)/b**2)*c2alpha
    
        A = 1 + u2/16384 *(4096 + u2*(-768 + u2*(320 - 175*u2)))
        B = u2/1024 * (256 + u2*(-128 + u2* (74 - 47*u2)))
    
        dsigma = B * ssigma *(c2sigmam + (1/4)*B*(csigma*(-1 + 2*(c2sigmam)**2)- 1/6*B*c2sigmam*(-3 + 4*(ssigma)**2*(-3 + 4*(c2sigmam)**2))))
    
        sab = b * A*(sigma - dsigma)
    
        Aab = atan((cos(U_b) * sin(L))/(cos(U_a)*sin(U_b) - sin(U_a) * cos(U_b) * cos(L)))
        Aba = atan((cos(U_a) * sin(L))/(-sin(U_a)*cos(U_b) + cos(U_a)*sin(U_b)*cos(L))) + pi
    
        if Aba > 2*pi:
            Aba = Aba - 2*pi


        if Aab < 0:
            Aab = Aab + 2*pi
        
        return(sab, Aab, Aba)

# Odleglosc 3D
    def odleglosc_3D(self,Xa,Xb,Ya,Yb,Za,Zb):
        
        delta_x = Xb - Xa
        delta_y = Yb - Ya
        delta_z = Zb - Za
        odleglosc = sqrt(delta_x**2+delta_y**2+delta_z**2)
        
        return(odleglosc)
    
#kat elewacji
    def az_el(self, fia, lama, ha, fib, lamb, hb):  

            Na = self.a/(sqrt(1-self.ecc2*(sin(fia)**2)))
            Nb = self.a/(sqrt(1-self.ecc2*(sin(fib)**2)))
        
            wRr = matrix([[(Na+ha)*cos(fia)*cos(lama)],
                            [(Na+ha)*cos(fia)*sin(lama)],
                            [(Na*(1-self.ecc2)+ha)*sin(fia)]])
            print(wRr)
            wRs = matrix([[(Nb+hb)*cos(fib)*cos(lamb)],
                            [(Nb+hb)*cos(fib)*sin(lamb)],
                            [(Nb*(1-self.ecc2)+hb)*sin(fib)]])
        
            R = wRs - wRr
            dlR = sqrt(R[0, 0]**2 + R[1, 0]**2 + R[2, 0]**2)
        
            wR = matrix([[R[0, 0]/dlR],
                            [R[1, 0]/dlR],
                            [R[2, 0]/dlR]])
        
            u = matrix([[cos(fia)*cos(lama)],
                           [cos(fia)*sin(lama)],
                           [sin(fia)]])
        
            n = matrix([[-sin(fia)*cos(lama)],
                           [-sin(fia)*sin(lama)],
                           [cos(fia)]])
        
            e = matrix([[-sin(lama)],
                           [cos(lama)],
                           [0]])
        
            alfa = atan((transpose(wR)*e)/(transpose(wR)*n))*180/pi+180
            print(alfa)
            azymut = Transformacje.st2sms(alfa)
            
        
            z = acos(transpose(u)*wR)*180/pi
            e = 90 - z
            E = Transformacje.st2sms(e)
            
            print('Azymut = ', azymut)
            print('Kąt elewacji = ', E)
        
            return azymut, E





        
if __name__ == "__main__":
    # utworzenie obiektu
    geo = Transformacje(model = "wgs84")
    # dane XYZ geocentryczne
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    phi, lam, h = geo.xyz2plh(X, Y, Z)
    print(phi, lam, h)
    phi, lam, h = geo.xyz2plh(X, Y, Z)
    print(phi, lam, h)
    
