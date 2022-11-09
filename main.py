from math import inf
import numpy as np
import matplotlib.pyplot as plt
from Objects import Plane, Sphere, Triangle
from Functions import *
import random

def readline(Type): return map(Type, input().split())

if __name__ == "__main__":
    ########### -- INPUT -- ###########
    height, width = readline(int)
    s, d = readline(float)
    E_x, E_y, E_z = readline(float)
    L_x, L_y, L_z = readline(float)
    UP_x, UP_y, UP_z = readline(float)
    BC_R, BC_G, BC_B = readline(int)
    qtd_obj = int(input())
    ###################################
    E_point = np.array([E_x, E_y, E_z])
    L_point = np.array([L_x, L_y, L_z])
    up_vector = np.array([UP_x, UP_y, UP_z])
    BC_RGB = np.array([BC_R, BC_G, BC_B])

    objects = []


    for i in range(qtd_obj):
        object_definition = str(input())
        
        if '/' in object_definition:
            RGB_plane, plane_definition = object_definition.split(' / ')

            RGB_plane = RGB_plane.split(' ')
            plane_definition = plane_definition.split(' ')
            
            R_plane = int(RGB_plane[0])
            G_plane = int(RGB_plane[1])
            B_Plane = int(RGB_plane[2])
            

            P_point_x = float(plane_definition[0])
            P_point_y = float(plane_definition[1])
            P_point_z = float(plane_definition[2])
            N_vector_x = float(plane_definition[3])
            N_vector_y = float(plane_definition[4])
            N_vector_z = float(plane_definition[5])

            plane = Plane([N_vector_x, N_vector_y, N_vector_z], [P_point_x, P_point_y, P_point_z])
            plane.coloring([R_plane, G_plane, B_Plane])
            print([R_plane, G_plane, B_Plane])
            objects.append(plane)

        if '*' in object_definition:
            RGB_sphere, sphere_definition = object_definition.split(' * ')
            
            RGB_sphere = RGB_sphere.split(' ')
            sphere_definition = sphere_definition.split(' ')

            R_sphere = int(RGB_sphere[0])
            G_sphere = int(RGB_sphere[1])
            B_sphere = int(RGB_sphere[2])

            C_point_x = float(sphere_definition[0])
            C_point_y = float(sphere_definition[1])
            C_point_z = float(sphere_definition[2])
            radius = float(sphere_definition[3])

            sphere = Sphere([C_point_x, C_point_y , C_point_z], radius)
            sphere.coloring([R_sphere, G_sphere, B_sphere])
            objects.append(sphere)

        if '>' in object_definition:
            RGB_triangle, triangle_definition = object_definition.split(' > ')

            RGB_triangle = RGB_triangle.split(' ')
            triangle_definition = triangle_definition.split(' ')
            
            R_triangle = int(RGB_triangle[0])
            G_triangle = int(RGB_triangle[1])
            B_triangle = int(RGB_triangle[2])
            
            
            A_point_x = float(triangle_definition[0])
            A_point_y = float(triangle_definition[1])
            A_point_z = float(triangle_definition[2])

            B_point_x = float(triangle_definition[3])
            B_point_y = float(triangle_definition[4])
            B_point_z = float(triangle_definition[5])

            C_point_x = float(triangle_definition[6])
            C_point_y = float(triangle_definition[7])
            C_point_z = float(triangle_definition[8])
            #print('A: ' + str([A_point_x, A_point_y, A_point_z]) + ' B: '+ str([B_point_x, B_point_y, B_point_z]) + ' C: ' +  str([C_point_x, C_point_y, C_point_z]))


            triangle = Triangle([A_point_x, A_point_y, A_point_z], [B_point_x, B_point_y, B_point_z], [C_point_x, C_point_y, C_point_z])
            #print([R_triangle, G_triangle, B_triangle])
            triangle.coloring([R_triangle, G_triangle, B_triangle])

            objects.append(triangle)

    
    imagem = image(objects, E_point, L_point, up_vector, BC_RGB, height, width, d, s)
    plt.imsave("./image.png", imagem)