import numpy as np
from math import inf, sqrt

#Objetos

#Esfera - Centro e raio 

class Sphere:
    def __init__(self, center, radius, ka, kd, ks, exp, kr, kt, refraction_index):
        self.center = np.array(center)
        self.radius = radius
        #Coeficientes
        self.ka = ka
        self.kd = kd
        self.ks = ks
        self.exp = exp
        self.kr = kr
        self.kt = kt
        self.refraction_index = refraction_index

    def coloring(self, RGB_color):
        self.color = np.array(RGB_color)

    def getColor(self):
        return self.color
    
    #Calcula a normal para ajudar a pegar a interseção
    def getNormal(self, P):
        normal_vector = np.array(P - self.center) 
        normal_vector = ( normal_vector/ np.linalg.norm(normal_vector))
        return normal_vector


    def raysphere_intersect(center, radius, ray_origin, ray_direction):
        I = center - ray_origin
        tmin = np.dot(I, ray_direction)
        d_power_2 = np.dot(I,I) - (tmin**2)
        if d_power_2 <= (radius**2):
            delta = sqrt(radius**2 - d_power_2)
            t1 = tmin - delta
            t2 = tmin + delta
            if t1 < 0:
                if t2 < 0:
                    t = inf
                else:
                    t = t2
            else:
                t = t1
        else:
            t = inf
        return t     

    def __str__(self):
        return "Sphere"



#Plano - Vetor normal e ponto contido no plano

class Plane:
    def __init__(self, normal_vector, P_point, ka, kd, ks, exp, kr, kt, refraction_index):
        self.normal_vector = np.array(normal_vector)
        self.P_point = np.array(P_point)
        #Coeficientes
        self.ka = ka
        self.kd = kd
        self.ks = ks
        self.exp = exp
        self.kr = kr
        self.kt = kt
        self.refraction_index = refraction_index
    
    def coloring(self, RGB_color):
        self.color = np.array(RGB_color)

    def getColor(self):
        return self.color
    
    def getNormal(self, P):
        return self.normal_vector

    def rayplane_intersect(normal_vector, P_point, ray_origin, ray_direction, epsilon = 1e-6, infinito = inf ):
        v = np.dot(ray_direction, normal_vector)

	
        if abs(v) > epsilon:

            w = P_point - ray_origin 
            h = np.dot(w, normal_vector) 
            t = h / v
            
            if t < 0:
                t = infinito
        else:
            t = infinito

        return t

    def __str__(self):
        return "Plane"
    
#Triangulo - 3 pontos

class Triangle:
    def __init__(self, A_point, B_point, C_point, ka, kd, ks, exp, kr, kt, refraction_index):
        self.A_point = np.array(A_point)
        self.B_point = np.array(B_point)
        self.C_point = np.array(C_point)
        #Coeficientes
        self.ka = ka
        self.kd = kd
        self.ks = ks
        self.exp = exp
        self.kr = kr
        self.kt = kt
        self.refraction_index = refraction_index
    
    def coloring(self, RGB_color):
        self.color = np.array(RGB_color)
    
    def getColor(self):
        return self.color
    
    def getNormal(self,P):
        AB = self.B_point - self.A_point
        AC = self.C_point - self.A_point
        normal_vector = np.cross(AB,AC) 
        normal_vector = ( normal_vector/ np.linalg.norm(normal_vector))
        return normal_vector


    def raytriangle_intersect(A_point, B_point, C_point, ray_origin, ray_direction):
        #Pré processamento        
        u = np.array(B_point - A_point)
        v = np.array(C_point - A_point)
        n = np.cross(u,v)
        n_ = n/ np.linalg.norm(n)

        t = Plane.rayplane_intersect(n_, A_point, ray_origin, ray_direction, epsilon = 1e-6, infinito = inf)

        projuv = (np.dot(u,v)/np.dot(v,v)) * v
        projvu = (np.dot(v,u)/np.dot(u,u)) * u
        
        hb = u - projuv
        hc = v - projvu
        hb = hb/(np.dot(hb,hb))
        hc = hc/(np.dot(hc,hc))
        
        #Cálculo da interseção
       
        if t  < inf:
            P_point = ray_origin + t * ray_direction

            AP = P_point - A_point

            beta = np.dot(AP, hb)
            gama =  np.dot(AP, hc)
            alpha = 1 - (beta + gama)

            if alpha < 0 or beta < 0 or gama < 0:
                t = inf
        return t
    
    def __str__(self):
        return "Triangle"