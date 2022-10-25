from math import inf
import numpy as np
import matplotlib.pyplot as plt
from Objects import Plane, Sphere, Triangle
from Functions import *

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
            color_coeficients_plane, plane_definition = object_definition.split(' / ')

            color_coeficients_plane = color_coeficients_plane.split(' ')
            plane_definition = plane_definition.split(' ')

            #print(color_coeficients_plane)
            
            Cd_R_plane = int(color_coeficients_plane[0])
            Cd_G_plane = int(color_coeficients_plane[1])
            Cd_B_plane = int(color_coeficients_plane[2])

            ka = float(color_coeficients_plane[3])
            kd = float(color_coeficients_plane[4])
            ks = float(color_coeficients_plane[5])
            exp = float(color_coeficients_plane[6])
            

            P_point_x = float(plane_definition[0])
            P_point_y = float(plane_definition[1])
            P_point_z = float(plane_definition[2])
            N_vector_x = float(plane_definition[3])
            N_vector_y = float(plane_definition[4])
            N_vector_z = float(plane_definition[5])

            plane = Plane([N_vector_x, N_vector_y, N_vector_z], [P_point_x, P_point_y, P_point_z], ka, kd, ks, exp)
            plane.coloring([Cd_R_plane/255, Cd_G_plane/255, Cd_B_plane/255])
            print(f'Cd_RGB: {Cd_R_plane, Cd_G_plane, Cd_B_plane}')
            objects.append(plane)

        if '*' in object_definition:
            color_coeficients_sphere, sphere_definition = object_definition.split(' * ')

            #print(color_coeficients_sphere)
            #print(sphere_definition)
            
            color_coeficients_sphere = color_coeficients_sphere.split(' ')
            sphere_definition = sphere_definition.split(' ')

            Cd_R_sphere = int(color_coeficients_sphere[0])
            Cd_G_sphere = int(color_coeficients_sphere[1])
            Cd_B_sphere = int(color_coeficients_sphere[2])

            #print(f'Cd_RGB: {Cd_R_sphere, Cd_G_sphere, Cd_B_sphere}')
            
            ka = float(color_coeficients_sphere[3])
            kd = float(color_coeficients_sphere[4])
            ks = float(color_coeficients_sphere[5])
            exp = float(color_coeficients_sphere[6])

            #print(f'ka, kd, ks, exp: {ka, kd, ks, exp}')

            C_point_x = float(sphere_definition[0])
            C_point_y = float(sphere_definition[1])
            C_point_z = float(sphere_definition[2])
            radius = float(sphere_definition[3])

            #print(f'Center: {[C_point_x, C_point_y , C_point_z]}')
            #print(f'Radius: {radius}')

            sphere = Sphere([C_point_x, C_point_y , C_point_z], radius, ka, kd, ks, exp)
            sphere.coloring([Cd_R_sphere/255, Cd_G_sphere/255, Cd_B_sphere/255])
            objects.append(sphere)

        if '>' in object_definition:
            color_coeficients_triangle, triangle_definition = object_definition.split(' > ')

            color_coeficients_triangle = color_coeficients_triangle.split(' ')
            triangle_definition = triangle_definition.split(' ')
            
            Cd_R_triangle = int(color_coeficients_triangle[0])
            Cd_G_triangle = int(color_coeficients_triangle[1])
            Cd_B_triangle = int(color_coeficients_triangle[2])

            ka = float(color_coeficients_triangle[3])
            kd = float(color_coeficients_triangle[4])
            ks = float(color_coeficients_triangle[5])
            exp = float(color_coeficients_triangle[6])
            
            
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


            triangle = Triangle([A_point_x, A_point_y, A_point_z], [B_point_x, B_point_y, B_point_z], [C_point_x, C_point_y, C_point_z], ka, kd, ks, exp)
            #print([R_triangle, G_triangle, B_triangle])
            triangle.coloring([Cd_R_triangle/255, Cd_G_triangle/255, Cd_B_triangle/255])

            objects.append(triangle)
    
    Ca_R, Ca_G, Ca_B = readline(float)
    Ca = np.array([Ca_R/255, Ca_G/255, Ca_B/255])
    qtd_lights = int(input())

    lights = []

    for i in range(qtd_lights):
         c_R, c_G, c_B, L_x, L_y, L_z = readline(float)
         #print([I_R, I_G, I_B])
         lights.append((np.array((c_R/255, c_G/255, c_B/255)), np.array((L_x, L_y, L_z))))

    image = image(objects, E_point, L_point, up_vector, BC_RGB, height, width, d, s, Ca, lights)
    plt.imsave("./image.png", image)