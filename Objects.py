import numpy as np
from math import inf, sqrt
#Objetos

#Esfera - Centro e raio
class Sphere:
    def __init__(self, center, radius, ka, kd, ks, exp):
        self.center = np.array(center)
        self.radius = radius
        self.ka = ka
        self.kd = kd
        self.ks = ks
        self.exp = exp

    def coloring(self, RGB_color):
        self.color = np.array(RGB_color)

    def getColor(self):
        return self.color
    
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
    
    
    def __init__(self, normal_vector, P_point, ka, kd, ks, exp):
        self.normal_vector = np.array(normal_vector)
        self.P_point = np.array(P_point)
        self.ka = ka
        self.kd = kd
        self.ks = ks
        self.exp = exp
    
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
    def __init__(self, A_point, B_point, C_point, ka, kd, ks, exp):
        self.A_point = np.array(A_point)
        self.B_point = np.array(B_point)
        self.C_point = np.array(C_point)
        self.ka = ka
        self.kd = kd
        self.ks = ks
        self.exp = exp
    
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

       
        if t  < inf:
            P_point = ray_origin + t * ray_direction

            AP = P_point - A_point

            beta = np.dot(AP, hb)
            gama =  np.dot(AP, hc)
            alpha = 1 - (beta + gama)

            if alpha < 0 or beta < 0 or gama < 0:
                #print('intersect')
                t = inf
        return t
        
        #print('alpha: ' + str(alpha))
        #print('beta: ' + str(beta))
        #print('gama: ' + str(gama))
 

        """
        AB = B_point - A_point               # Oriented segment A to B
        AC = C_point - A_point               # Oriented segment A to C
        ray_direction_ = ray_direction/np.linalg.norm(ray_direction)
        n = np.cross(AB, AC)     # Normal vector
        n_ = n/np.linalg.norm(n) # Normalized normal
        # Using the point A to find d
        d = - np.dot(n_, A_point)
        # Finding parameter t
        t = - (np.dot(n_, ray_origin) + d)/np.dot(n_, ray_direction_)
        # Finding P
        P = ray_origin + t * ray_direction_
        # Get the resulting vector for each vertex
        # following the construction order
        Pa = np.dot(np.cross(B_point - A_point, P - A_point), n_)
        Pb = np.dot(np.cross(C_point - B_point, P - B_point), n_)
        Pc = np.dot(np.cross(A_point - C_point, P - C_point), n_)
        if(t < 0):
            return
        elif(Pa < 0 and Pb < 0 and Pc < 0):
            return 
        else:
            return t """
    
    def __str__(self):
        return "Triangle"