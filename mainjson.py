import json
import numpy as np
import matplotlib.pyplot as plt
from Objects import Plane, Sphere
from Functions import image, cast

path='eclipse.json'
with open(path) as f:
    data = json.load(f)
    v_res = data['v_res']
    h_res = data['h_res']
    s = data['square_side']
    d = data['dist']
    E = np.array(data['eye'])
    L = np.array(data['look_at'])
    up = np.array(data['up'])
    background_color = np.array(data['background_color'])/255
    objects = []
    for obj in data['objects']:
        if 'plane' in obj:
            plane = Plane(np.array(obj['plane']['sample']), np.array(obj['plane']['normal']))
            plane.coloring(np.array(obj['color'])/255)
            plane.set_illumination(obj['ka'], obj['kd'], obj['ks'], obj['exp'])
            objects.append(plane)
        if 'sphere' in obj:
            sphere = Sphere(np.array(obj['sphere']['center']), 
            np.array(obj['sphere']['radius']),
            obj['ka'], 
            obj['kd'],
            obj['ks'],
            obj['exp'])
            sphere.coloring(np.array(obj['color'])/255)
            objects.append(sphere)
    ambient_light = np.array(data['ambient_light'])/255
    lights = []
    for light in data['lights']:
        lights.append((np.array(light['intensity'])/255, np.array(light['position'])))
    img = image(objects, E, L, up, background_color, v_res, h_res, d, s, ambient_light, lights)
    plt.imshow(img)
    plt.show()