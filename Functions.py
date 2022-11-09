from math import inf,sqrt
import numpy as np
from Objects import Plane, Sphere, Triangle
from main import *

#Função de normalização
def normalize(vector):
    return vector / np.linalg.norm(vector)

#Trace - retorna o objeto mais próximo e o parâmetro de distancia mínima
def trace(objects, ray_origin, ray_direction):
    distances = []
    for object in objects:
        if object.__str__() == "Plane": 
            distances.append(Plane.rayplane_intersect( object.normal_vector, object.P_point , ray_origin, ray_direction))
        if object.__str__() == "Sphere":
            distances.append(Sphere.raysphere_intersect( object.center, object.radius, ray_origin, ray_direction))
        if object.__str__() == "Triangle":
            distances.append(Triangle.raytriangle_intersect( object.A_point, object.B_point, object.C_point, ray_origin, ray_direction))
    nearest_object = None
    min_distance = np.inf
    for index, distance in enumerate(distances):
        if distance and distance < min_distance:
            min_distance = distance
            nearest_object = objects[index]
    return min_distance, nearest_object



def reflect(l_vector,n_vector):
    return 2 * np.dot(n_vector,l_vector) * n_vector - l_vector


#Iluminação
def phong(objects, nearest_object, intersection_point, Ca, lights, w0_vector, normal_vector):

    cd = nearest_object.getColor()
    ka = nearest_object.ka
    kd = nearest_object.kd
    ks = nearest_object.ks
    q = nearest_object.exp
    
    #Componente ambiental
    Cp = np.multiply(ka*cd , Ca)
    
    for c,L in lights:
        l_vector = normalize(np.array(L - intersection_point))
        r_vector = reflect(l_vector, normal_vector)
        shift_point = intersection_point + 10E-5 * l_vector
        tmin,nearest_object =  trace(objects, shift_point, l_vector)
        
        if not nearest_object or (np.dot(l_vector, (L-shift_point) ) < tmin):
            if np.dot(normal_vector, l_vector) > 0:
                #Componente difusa
                Cp += np.multiply((kd * cd), (np.dot(normal_vector, l_vector) * c))
            if np.dot(w0_vector, r_vector) > 0:
                #Componente especular
                Cp += ks * (np.dot(w0_vector, r_vector) ** q) * c
    return Cp

#Cálculo da refração
def refract(nearest_object, w0_vector, normal_vector):
    costheta = np.dot(normal_vector, w0_vector)
    refraction_index = nearest_object.refraction_index
    if costheta < 0:
        normal_vector = -normal_vector
        refraction_index = 1/refraction_index
        costheta = -costheta
    delta = 1 - (1/(refraction_index**2))*(1-costheta**2)
    if delta < 0:
        return
    
    return - ((1/refraction_index) * w0_vector) - (sqrt(delta) - (1/refraction_index)*costheta)*normal_vector
    

#Função de cast
def cast(objects, ray_origin, ray_direction, background_color, Ca, lights, ttl):
    color = background_color
    #Parametro de ditancia minima e objeto mais próximo
    tmin,nearest_object = trace(objects, ray_origin, ray_direction)
    #Se tiver objeto mais próximo
    if nearest_object:
        kt = nearest_object.kt
        kr = nearest_object.kr
        #Ponto de interseção
        intersection_point = ray_origin + (tmin * ray_direction)
        
        w0_vector = -ray_direction
        normal_vector = nearest_object.getNormal(intersection_point)
        color = phong(objects, nearest_object, intersection_point, Ca, lights, w0_vector, normal_vector)
        
        #Usa o max_depth(ttl)
        if ttl > 0:
            r = reflect(w0_vector, normal_vector)
            shift_intersection_point_r = intersection_point + 10e-5  * r
            try:
                if kt > 0:
                    t = refract(nearest_object, w0_vector, normal_vector)
                    shift_intersection_point_t = intersection_point + 10e-5 * t
                    color += kt * cast(objects, shift_intersection_point_t, t, background_color, Ca, lights, ttl -1)
                if kr > 0:
                    color += kr * cast(objects, shift_intersection_point_r, r, background_color, Ca, lights, ttl -1)
            except:
                color += cast(objects, shift_intersection_point_r, r, background_color, Ca, lights, ttl-1)  
    return color





def image(objects, E_point, L_point, up_vector, BC_RGB, height, width, d, s, Ca, lights, max_depth):
    #Vetores w/ u/ up_vector
    w_vector = normalize(np.array(E_point - L_point))
    u_vector = normalize(np.cross(up_vector, w_vector))
    v_vector = np.cross(w_vector, u_vector)
    v_res = height
    h_res = width
    dist = d
    pixel_side = s
    
    #Calculando centro da tela
    C_point = np.array(E_point - (dist*w_vector))
    
    #Inicializando matriz Q
    Q_matrix = np.zeros((v_res, h_res, 3))
    image = np.zeros((v_res, h_res, 3))
    c = np.zeros((v_res, h_res, 3))

    #Cálculo do Q00 -- ponto esquerda/cima
    Q_matrix[0,0] = np.array((C_point + 1/2 * pixel_side * (v_res - 1) * v_vector) - 
    (1/2 * pixel_side * (h_res - 1) * u_vector))


    for i in range(v_res):
        for j in range(h_res):
            #Cálculo do Qij
            Q_matrix[i,j] = Q_matrix[0,0] + (pixel_side * j * u_vector) - (pixel_side * i * v_vector)
            ray_direction = normalize(Q_matrix[i, j] - E_point)
            aux = cast(objects, E_point, ray_direction, BC_RGB, Ca, lights, max_depth)

            #Tratando o Overflow
            aux = aux/max(*aux,1)
            image[i , j] = aux 
            
    return image