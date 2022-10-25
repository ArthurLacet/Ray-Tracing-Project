from math import inf
import numpy as np
from Objects import Plane, Sphere, Triangle
from main import *


#Função de normalização
def normalize(vector):
    return vector / np.linalg.norm(vector)




def trace(objects, ray_origin, ray_direction):
    distances = []
    #print('opa')
    for object in objects:
        if object.__str__() == "Plane": 
            #print('O: ' + str(ray_origin))
            #print('D: ' + str(ray_direction))
            distances.append(Plane.rayplane_intersect( object.normal_vector, object.P_point , ray_origin, ray_direction))
        if object.__str__() == "Sphere":
            #print('O: ' + str(ray_origin))
            #print('D: ' + str(ray_direction))
            distances.append(Sphere.raysphere_intersect( object.center, object.radius, ray_origin, ray_direction))
        if object.__str__() == "Triangle":
            #print('A: ' + str(object.A_point) + ' B: '+ str(object.B_point) + ' C: ' +  str(object.C_point))
            #print('O: ' + str(ray_origin))
            #print('D: ' + str(2*ray_direction)) 
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

def phong(objects, nearest_object, intersection_point, Ca, lights, w0_vector, normal_vector):
    cd = nearest_object.getColor()
    #print(cd)
    ka = nearest_object.ka
    kd = nearest_object.kd
    ks = nearest_object.ks
    q = nearest_object.exp
    
    Cp = np.multiply(ka*cd , Ca)
    
    for c,L in lights:
        l_vector = normalize(np.array(L - intersection_point))
        r_vector = reflect(l_vector, normal_vector)
        shift_point = intersection_point + 10e-5 * l_vector
        tmin,nearest_object =  trace(objects, intersection_point, l_vector)
        
        if not nearest_object or (np.dot(l_vector, (L-shift_point) ) < tmin):
            if np.dot(normal_vector, l_vector) > 0:
                Cp += np.multiply((kd * cd), (np.dot(normal_vector, l_vector) * c))
            
            if np.dot(w0_vector, r_vector) > 0:
                Cp += ks * (np.dot(w0_vector, r_vector) ** q) * c
    return Cp



def cast(objects, ray_origin, ray_direction, background_color, Ca, lights):
    color = background_color
    #print(f'objects:{objects}')
    tmin,nearest_object = trace(objects, ray_origin, ray_direction)
    #print(f'nearest:{nearest_object}')
    if nearest_object:
        intersection_point = ray_origin + (tmin * ray_direction)
        w0_vector = -ray_direction
        normal_vector = nearest_object.getNormal(intersection_point)
        color = phong(objects, nearest_object, intersection_point, Ca, lights, w0_vector, normal_vector)
    #print(color)    
    return color





def image(objects, E_point, L_point, up_vector, BC_RGB, height, width, d, s, Ca, lights):
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

    #Cálculo do Q00
    Q_matrix[0,0] = np.array((C_point + 1/2 * pixel_side * (v_res - 1) * v_vector) - 
    (1/2 * pixel_side * (h_res - 1) * u_vector))


    for i in range(v_res):
        for j in range(h_res):
            #Cálculo do Qij
            Q_matrix[i,j] = Q_matrix[0,0] + (pixel_side * j * u_vector) - (pixel_side * i * v_vector)
            ray_direction = normalize(Q_matrix[i, j] - E_point)
            aux = cast(objects, E_point, ray_direction, BC_RGB, Ca, lights)
            aux = aux/max(*aux, 1)
            #print(aux)
            image[i , j] = aux
            
    return image