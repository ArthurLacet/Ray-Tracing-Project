from math import inf
import numpy as np
from Objects import Plane, Sphere, Triangle
from main import *


#Função de normalização
def normalize(vector):
    return vector / np.linalg.norm(vector)

def trace(objects, ray_origin, ray_direction):
    distances = []
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



def cast(objects, ray_origin, ray_direction, background_color):
    color = background_color
    tmin,nearest_object = trace(objects, ray_origin, ray_direction)
    if nearest_object:
        color = nearest_object.getColor()
    return color





def image(objects, E_point, L_point, up_vector, BC_RGB, height, width, d, s):
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
            image[i, j] = cast(objects, E_point, ray_direction, BC_RGB)

    return image / 255