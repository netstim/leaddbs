import bpy
import csv
import sys
import json

print("Starting Blender Script")

# Function to read vertices from CSV file
def read_vertices(file):
    vertices = []
    with open(file, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            vertices.append((float(row[0]), float(row[1]), float(row[2])))
    return vertices

# Function to read faces from CSV file
def read_faces(file):
    faces = []
    with open(file, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            faces.append([int(index) for index in row])
    return faces

# Function to read lines from CSV file
def read_lines(file):
    lines = []
    with open(file, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            lines.append((float(row[0]), float(row[1]), float(row[2])))
    return lines

# Function to read surface data from CSV files
def read_surface(fileX, fileY, fileZ):
    X, Y, Z = [], [], []
    with open(fileX, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            X.append([float(value) for value in row])
    with open(fileY, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            Y.append([float(value) for value in row])
    with open(fileZ, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            Z.append([float(value) for value in row])
    return X, Y, Z

# Function to read streamtube data from CSV files
def read_streamtube(fileX, fileY, fileZ, fileU, fileV, fileW):
    X, Y, Z, U, V, W = [], [], [], [], [], []
    with open(fileX, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            X.append([float(value) for value in row])
    with open(fileY, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            Y.append([float(value) for value in row])
    with open(fileZ, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            Z.append([float(value) for value in row])
    with open(fileU, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            U.append([float(value) for value in row])
    with open(fileV, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            V.append([float(value) for value in row])
    with open(fileW, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            W.append([float(value) for value in row])
    return X, Y, Z, U, V, W

# Create mesh for patch
def create_patch(vertices, faces):
    mesh = bpy.data.meshes.new("PatchMesh")
    obj = bpy.data.objects.new("PatchObject", mesh)
    bpy.context.collection.objects.link(obj)
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)
    mesh.from_pydata(vertices, [], faces)
    mesh.update()

# Create lines for plot3
def create_plot3(lines):
    for i in range(len(lines) - 1):
        bpy.ops.mesh.primitive_cylinder_add(radius=0.02, depth=1, location=lines[i])
        cylinder = bpy.context.object
        cylinder.rotation_euler = (lines[i+1][0] - lines[i][0], lines[i+1][1] - lines[i][1], lines[i+1][2] - lines[i][2])

# Create surface
def create_surface(X, Y, Z):
    vertices = [(X[i][j], Y[i][j], Z[i][j]) for i in range(len(X)) for j in range(len(X[0]))]
    faces = [(i + j*len(X), (i+1) + j*len(X), (i+1) + (j+1)*len(X), i + (j+1)*len(X)) for i in range(len(X)-1) for j in range(len(X[0])-1)]
    create_patch(vertices, faces)

# Create streamtubes
def create_streamtube(X, Y, Z, U, V, W):
    for i in range(len(X)):
        for j in range(len(X[0])):
            bpy.ops.mesh.primitive_uv_sphere_add(radius=0.1, location=(X[i][j], Y[i][j], Z[i][j]))
            sphere = bpy.context.object
            sphere.data.materials.append(create_material("Blue", (0, 0, 1, 1)))

# Create a material with color
def create_material(name, color):
    mat = bpy.data.materials.new(name)
    mat.diffuse_color = color
    return mat

# Add light source
def add_light(location, type='SUN', energy=3):
    bpy.ops.object.light_add(type=type, location=location)
    light = bpy.context.object
    light.data.energy = energy

# Read command-line arguments
plot_type = sys.argv[-1]

try:
    print(f"Running plot type: {plot_type}")

    # Delete default cube
    if 'Cube' in bpy.data.objects:
        bpy.data.objects['Cube'].select_set(True)
        bpy.ops.object.delete()

    if plot_type == 'patch':
        vertices = read_vertices('vertices.csv')
        faces = read_faces('faces.csv')
        create_patch(vertices, faces)
    elif plot_type == 'plot3':
        lines = read_lines('lines.csv')
        create_plot3(lines)
    elif plot_type == 'surf':
        X, Y, Z = read_surface('X.csv', 'Y.csv', 'Z.csv')
        create_surface(X, Y, Z)
    elif plot_type == 'streamtube':
        X, Y, Z, U, V, W = read_streamtube('X.csv', 'Y.csv', 'Z.csv', 'U.csv', 'V.csv', 'W.csv')
        create_streamtube(X, Y, Z, U, V, W)
    elif plot_type == 'light':
        with open('lights.json', 'r') as f:
            lights = json.load(f)
        for light in lights:
            add_light(location=light['position'], type=light.get('type', 'SUN'), energy=light.get('energy', 3))
    elif plot_type == 'camlight':
        with open('camlights.json', 'r') as f:
            camlights = json.load(f)
        for camlight in camlights:
            add_light(location=camlight['position'], type=camlight.get('type', 'SUN'), energy=camlight.get('energy', 3))

    # Set camera view
    bpy.ops.object.camera_add(location=(0, 0, 10))
    camera = bpy.context.object
    camera.rotation_euler = (0, 0, 0)
    bpy.context.scene.camera = camera

    print("Finished executing plot type: ", plot_type)

except Exception as e:
    print(f"Error during execution: {e}")