from math import inf
import numpy as np
import matplotlib.pyplot as plt
from Objects import Plane, Sphere, Triangle
from Functions import *

def readline(Type): return map(Type, input().split())

if __name__ == "__main__":
    #Recebendo input
    height, width = readline(int)
    s, d = readline(float)
    E_x, E_y, E_z = readline(float)
    L_x, L_y, L_z = readline(float)
    UP_x, UP_y, UP_z = readline(float)
    BC_R, BC_G, BC_B = readline(int)
    max_depth = int(input())
    qtd_obj = int(input())

    #Organizando dados estratégicos
    E_point = np.array([E_x, E_y, E_z])
    L_point = np.array([L_x, L_y, L_z])
    up_vector = np.array([UP_x, UP_y, UP_z])
    BC_RGB = np.array([BC_R, BC_G, BC_B])
   
    #Inicializa array de objetos vazio
    objects = []

    #Para cada objeto, coloca-se os dados no devido lugar e ao fim adiciona ao array de objetos 
    for i in range(qtd_obj):
        
        #Lê-se os dados
        object_definition = str(input())
        
        #PLANO
        if '/' in object_definition:
            color_coeficients_plane, plane_definition = object_definition.split(' / ')

            color_coeficients_plane = color_coeficients_plane.split(' ')
            plane_definition = plane_definition.split(' ')
            
            #RGB plane
            Cd_R_plane = int(color_coeficients_plane[0])
            Cd_G_plane = int(color_coeficients_plane[1])
            Cd_B_plane = int(color_coeficients_plane[2])
            
            #Coeficientes
            ka = float(color_coeficients_plane[3])
            kd = float(color_coeficients_plane[4])
            ks = float(color_coeficients_plane[5])
            exp = float(color_coeficients_plane[6])
            kr = float(color_coeficients_plane[7])
            kt = float(color_coeficients_plane[8])
            refraction_index = float(color_coeficients_plane[9])
            
            #Ponto P
            P_point_x = float(plane_definition[0])
            P_point_y = float(plane_definition[1])
            P_point_z = float(plane_definition[2])
            
            #Normal Vector
            N_vector_x = float(plane_definition[3])
            N_vector_y = float(plane_definition[4])
            N_vector_z = float(plane_definition[5])
            
            #Cria-se objeto plano
            plane = Plane([N_vector_x, N_vector_y, N_vector_z], [P_point_x, P_point_y, P_point_z], ka, kd, ks, exp, kr, kt, refraction_index)
            
            #Coloca-se a cor dele
            plane.coloring([Cd_R_plane/255, Cd_G_plane/255, Cd_B_plane/255])
            
            #Adiciona ao array de objetos 
            objects.append(plane)

        #ESFERA
        if '*' in object_definition:
            color_coeficients_sphere, sphere_definition = object_definition.split(' * ')
            color_coeficients_sphere = color_coeficients_sphere.split(' ')
            sphere_definition = sphere_definition.split(' ')

            #RGB esfera
            Cd_R_sphere = int(color_coeficients_sphere[0])
            Cd_G_sphere = int(color_coeficients_sphere[1])
            Cd_B_sphere = int(color_coeficients_sphere[2])
            
            #Coeficientes
            ka = float(color_coeficients_sphere[3])
            kd = float(color_coeficients_sphere[4])
            ks = float(color_coeficients_sphere[5])
            exp = float(color_coeficients_sphere[6])
            kr = float(color_coeficients_sphere[7])
            kt = float(color_coeficients_sphere[8])
            refraction_index = float(color_coeficients_sphere[9])

            #Centro
            C_point_x = float(sphere_definition[0])
            C_point_y = float(sphere_definition[1])
            C_point_z = float(sphere_definition[2])
            
            #Raio
            radius = float(sphere_definition[3])
            
            #Cria objeto esfera
            sphere = Sphere([C_point_x, C_point_y , C_point_z], radius, ka, kd, ks, exp, kr, kt, refraction_index)
            
            #Coloca-se a cor dela
            sphere.coloring([Cd_R_sphere/255, Cd_G_sphere/255, Cd_B_sphere/255])
            
            #Adiciona ao array de objetos
            objects.append(sphere)

        #TRIANGULO
        if '>' in object_definition:
            color_coeficients_triangle, triangle_definition = object_definition.split(' > ')

            color_coeficients_triangle = color_coeficients_triangle.split(' ')
            triangle_definition = triangle_definition.split(' ')
            
            #RGB triangle
            Cd_R_triangle = int(color_coeficients_triangle[0])
            Cd_G_triangle = int(color_coeficients_triangle[1])
            Cd_B_triangle = int(color_coeficients_triangle[2])

            #Coeficientes
            ka = float(color_coeficients_triangle[3])
            kd = float(color_coeficients_triangle[4])
            ks = float(color_coeficients_triangle[5])
            exp = float(color_coeficients_triangle[6])
            kr = float(color_coeficients_triangle[7])
            kt = float(color_coeficients_triangle[8])
            refraction_index = float(color_coeficients_triangle[9])
            
            #Ponto A
            A_point_x = float(triangle_definition[0])
            A_point_y = float(triangle_definition[1])
            A_point_z = float(triangle_definition[2])

            #Ponto B
            B_point_x = float(triangle_definition[3])
            B_point_y = float(triangle_definition[4])
            B_point_z = float(triangle_definition[5])

            #Ponto C
            C_point_x = float(triangle_definition[6])
            C_point_y = float(triangle_definition[7])
            C_point_z = float(triangle_definition[8])

            #Cria-se o triângulo
            triangle = Triangle([A_point_x, A_point_y, A_point_z], [B_point_x, B_point_y, B_point_z], [C_point_x, C_point_y, C_point_z], ka, kd, ks, exp, kr, kt, refraction_index)
            triangle.coloring([Cd_R_triangle/255, Cd_G_triangle/255, Cd_B_triangle/255])
            #Adiciona ao array de objetos
            objects.append(triangle)

    #RGB ambient
    Ca_R, Ca_G, Ca_B = readline(float)
    Ca = np.array([Ca_R/255, Ca_G/255, Ca_B/255])
    
    #Quantidade de fontes de luz
    qtd_lights = int(input())

    lights = []

    for i in range(qtd_lights):
         c_R, c_G, c_B, L_x, L_y, L_z = readline(float)
         lights.append((np.array((c_R/255, c_G/255, c_B/255)), np.array((L_x, L_y, L_z))))

    #Obtenção da imagem
    image = image(objects, E_point, L_point, up_vector, BC_RGB, height, width, d, s, Ca, lights, max_depth)

    plt.imsave("./image.png", image)

