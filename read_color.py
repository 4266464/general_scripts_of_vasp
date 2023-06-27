color_list=[]
with open("color","r") as f:
    for i in f.readlines():
        color_list+=[float(j) for j in i.split(" ") if j]

range_x,range_y,range_z=(72,72,72)
color_dict={}
for x in range(range_x):
    for y in range(range_y):
        for z in range(range_z):
            color_dict[(x,y,z)]=color_list[x*72*72+y*72+z]

dx=vx/a*range_x
dy=vy/a*range_x
dz=vz/a*range_x
new_color_dict={}
for x in range(range_x):
    for y in range(range_y):
        for z in range(range_z):
            new_color_dict[(x,y,z)]=new_color_dict[(x+dx,y+dy,z+dz)]




color_list.
print(max(color_list),min(color_list))
