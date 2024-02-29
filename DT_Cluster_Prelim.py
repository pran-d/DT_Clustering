import numpy as np 
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import math

def Area(vertex1, vertex2, vertex3):
    x1, y1 = vertex1
    x2, y2 = vertex2
    x3, y3 = vertex3
    area = abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2)
    return area

def find(point, clusterSet):
  clusterNum = -1
  counter = 0
  for cluster in clusterSet:
    if point in cluster:
      clusterNum = counter
      break
    counter=counter+1
  return clusterNum

def notIn(point, clusterSet):
  return (find(point, clusterSet)==-1)

def Shortest_edge(x1, y1, x2, y2, x3, y3):
  # Calculate the lengths of the sides
    a = math.sqrt((x2 - x3)**2 + (y2 - y3)**2)
    b = math.sqrt((x1 - x3)**2 + (y1 - y3)**2)
    c = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)
    if(a<=b and a<=c):
      ord = [0, 1, 2]
    elif(b<=c and b<=a):
      ord = [1, 2, 0]
    else:
      ord = [2, 1, 0]
    return (max(a,b,c), min(a,b,c), ord)

def sortingKey(simplex):
  return (Area(points[simplex[0]], points[simplex[1]], points[simplex[2]]))

def distance(pointA, pointB):
  return math.sqrt(math.pow(pointA[0]-pointB[0],2)+math.pow(pointA[1]-pointB[1],2))

points=np.array([[-8.0, -7.44], [4.08, -4.56], [-5.7, -7.55], [-0.7, 4.66], [-2.67, 0.28], [4.08, 2.11], [2.57, 4.71], [9.71, 7.66], [-0.36, 0.59], [-2.82, 9.74], [-8.37, 4.65], [-1.58, -1.64], [-9.63, 3.13], [-5.08, 5.53], [-7.19, -2.16], [-6.99, 1.45], [-1.09, -1.65], [6.51, 8.23], [0.98, -9.26], [1.15, 4.46], [-2.91, -7.77], [-1.04, 3.24], [6.13, -5.65], [6.87, 9.14], [-1.56, -3.97], [6.84, 8.35], [9.1, -1.15], [9.57, 3.86], [9.46, -1.78], [1.64, -4.83], [-5.52, 1.38], [1.29, 8.79], [-6.31, -5.5], [-9.18, 5.79], [5.42, -4.99], [4.69, -9.09], [3.36, -3.46], [-3.24, -5.74], [-1.53, -3.78], [-2.84, -2.13], [-2.92, -2.67], [-5.55, -9.94], [-7.74, 1.56], [-6.57, 5.43], [-8.23, -1.82], [5.92, -5.66], [0.66, 0.62], [8.51, 8.77], [6.43, -1.14], [-6.43, 8.43], [-9.47, 7.46], [-3.32, 1.51], [-7.74, -8.6], [-1.76, 3.49], [5.08, -0.95], [0.21, 0.52], [0.32, 9.04], [8.37, -7.84], [1.94, 1.7], [-6.25, -0.03], [-3.16, -5.27], [7.84, 0.75], [-3.3, -1.72], [-7.14, 6.15], [-3.23, -5.86], [-9.05, -5.86], [-8.45, 6.57], [-2.78, -5.8], [2.07, -5.92], [-9.64, 4.77], [4.64, 5.25], [-0.03, -9.48], [-1.16, 3.66], [7.43, -2.7], [2.23, -0.38], [-1.19, 3.72], [-0.34, -8.43], [-1.52, -0.35], [6.0, 1.57], [1.82, 5.28], [-2.14, 0.86], [0.4, -1.07], [-1.28, -9.78], [-8.49, -9.45], [-6.74, -7.37], [-3.55, -5.35], [0.53, -1.71], [0.26, -6.41], [-5.31, 6.91], [-4.48, 3.41], [-0.07, 6.11], [-9.02, 7.58], [0.28, 6.08], [-5.0, 3.41], [2.25, -5.38], [5.27, 0.11], [-2.33, -8.08], [6.89, 1.2], [-4.78, -8.73], [-1.26, 5.27], [4.83, -6.44], [-0.46, 1.4], [4.34, 0.35], [5.56, 8.09], [-8.73, 6.96], [-0.45, -3.13], [-9.4, 3.5], [9.48, -5.65], [-1.13, -2.5], [2.76, -9.52], [1.48, 3.37], [1.5, -2.37], [2.53, 2.1], [4.57, -7.67], [-7.03, -8.12], [4.99, 5.85], [-6.67, -2.05], [4.36, 8.52], [3.23, -5.28], [-7.91, 3.08], [5.8, -8.72], [-4.05, -9.9], [-8.1, -3.37], [-2.35, -2.08], [-8.51, 6.49], [-5.66, 3.22], [-8.9, 5.38], [1.48, 5.76], [5.48, 6.85], [9.83, 0.18], [-2.85, 4.96], [6.91, -0.91], [-1.52, 1.98], [-5.96, 3.41], [-8.64, -1.79], [-5.44, -8.32], [5.49, 8.93], [3.76, -9.83], [-7.48, 2.81], [-5.28, 6.7], [-7.21, -6.48], [4.97, 1.23], [-1.59, -9.16], [4.0, 3.19], [0.66, 1.16], [-5.99, -2.54], [4.79, -8.54], [6.98, -6.99], [-0.63, -4.52], [-3.03, 4.75], [1.9, 0.79], [-7.68, -9.01], [9.91, 5.73], [9.43, -9.49], [7.2, -8.86], [0.43, -5.35], [2.21, -3.46], [0.67, -6.32], [-4.2, 5.41], [5.62, 7.99], [-6.54, 2.14], [-6.43, 7.29], [3.18, 1.62], [6.32, 6.65], [5.87, -5.64], [-3.04, 3.2], [5.45, -3.47], [-7.99, 0.84], [1.06, -4.24], [7.44, 8.94], [8.93, 1.51], [1.85, -6.46], [-6.34, 7.17], [9.89, 3.63], [9.3, -2.15], [-7.11, 6.6], [9.68, 3.71], [6.6, 5.84], [-5.64, -6.25], [-7.15, 9.5], [9.65, 7.0], [-6.83, 6.73], [5.12, 9.7], [-0.45, 9.69], [-0.77, 3.45], [-2.48, -1.14], [-1.76, -6.23], [-3.04, -5.03], [5.94, -8.32], [-9.95, 0.62], [-0.99, 0.21], [6.94, -7.77], [2.19, -1.44], [-7.45, -3.01], [9.83, -8.22], [-2.76, 8.49], [-7.9, -3.81], [9.33, 1.8], [-8.35, -4.49], [-5.62, 6.69]])

tri=Delaunay(points)
clusters=[]
MECs=[]
areas=[]
total_area=0
cluster_areas = {}
mec_areas = []

# sort the triangles by increasing area
tri.simplices=sorted(tri.simplices, key=sortingKey)

for triangle in tri.simplices:
  currArea = Area(points[triangle[0]], points[triangle[1]], points[triangle[2]])
  areas.append(currArea)

# calculating average area of largest 5% of triangles 
for i in range(1, math.ceil(0.05*len(areas))):
  total_area = total_area+areas[-i]

area_lim = total_area/(0.05*len(areas))
sl_ratio = 0.12
breaking_ratio = 0.05
merging_ratio = 2

for idx,triangle in enumerate(tri.simplices):
  area=areas[idx]
  longest_edge, shortest_edge, order=Shortest_edge(points[triangle[0], 0], points[triangle[0],1], points[triangle[1],0], points[triangle[1],1], points[triangle[2],0], points[triangle[2],1])
  
  # for a small triangle
  if(area<area_lim):
    clusterNum0 = find(triangle[0], clusters)
    clusterNum1 = find(triangle[1], clusters)
    clusterNum2 = find(triangle[2], clusters)

    # if all three vertices are not added to any cluster yet
    if clusterNum0==-1 and clusterNum1==-1 and clusterNum2==-1:
      # if the shortest edge is much smaller than longest edge, add the point opp to shortest edge to one cluster, and endpoints of shortest edge to one cluster
      if(shortest_edge/longest_edge<sl_ratio):
        clusters.append([triangle[order[0]]])
        clusters.append([triangle[order[1]], triangle[order[2]]])
      # or else if triangle is small and no edge is much longer, add all three points to one cluster 
      else:
        clusters.append([triangle[0], triangle[1], triangle[2]])
    
    # if all three vertices are added to cluster already
    elif not(clusterNum0==-1 or clusterNum1==-1 or clusterNum2==-1):
      continue

    # any other case
    else:
      # if point 1 already added, add whichever point not added to cluster
      if(not clusterNum0==-1):
        if(clusterNum1==-1):
          # if thin triangle, add faraway point to diff cluster
          if(shortest_edge/longest_edge<sl_ratio and order[0]==1):
            clusters.append([triangle[1]])
          else:
            clusters[clusterNum0].append(triangle[1])
        if(clusterNum2==-1):
          if(shortest_edge/longest_edge<sl_ratio and order[0]==2):
            clusters.append([triangle[2]])
          else:
            clusters[clusterNum0].append(triangle[2])
      # if point 1 not yet added
      else:
        # if point 2 not yet added, then point 3 must have already been added
        if(clusterNum1==-1):
          # if triangle is thin, add the faraway vertex to one cluster and the remaining vertices to one cluster
          if(shortest_edge/longest_edge<sl_ratio):
            if(order[0]==0):
              clusters.append([triangle[0]])
              clusters[clusterNum2].append(triangle[1])
            elif(order[0]==1):
              clusters[clusterNum2].append(triangle[0])
              clusters.append([triangle[1]])
            else:
              clusters.append([triangle[0], triangle[1]])
          else:
            clusters[clusterNum2].append(triangle[0])
            clusters[clusterNum2].append(triangle[1])
        # if point 3 not yet added, then point 2 must have been added
        elif(clusterNum2==-1):
          if(shortest_edge/longest_edge<sl_ratio):
            if(order[0]==0):
              clusters.append([triangle[0]])
              clusters[clusterNum1].append(triangle[2])
            elif(order[0]==2):
              clusters[clusterNum2].append(triangle[0])
              clusters.append([triangle[2]])
            else:
              clusters.append([triangle[0], triangle[2]])
          else:
            clusters[clusterNum1].append(triangle[0])
            clusters[clusterNum1].append(triangle[2])
        # if point 2 and 3 both have been added
        else:
          if(shortest_edge/longest_edge<sl_ratio):
            if(order[0]==1):
              clusters[clusterNum2].append(triangle[0])
            elif(order[0]==2):
              clusters[clusterNum1].append(triangle[0])
            else:
              clusters.append([triangle[0]])
          else:
            if(distance(points[triangle[0]],points[triangle[1]])<distance(points[triangle[0]],points[triangle[2]])):
              clusters[clusterNum1].append(triangle[0])
            else:
              clusters[clusterNum2].append(triangle[0])
  
  # for a big triangle, automatically add all three to different clusters
  else:
    for i in range(3):
      if(notIn(triangle[i], clusters)):
        clusters.append([triangle[i]])

colors=['#4C2B9F', '#78290C', '#CA223D', '#C0A364', '#2C9645', '#C785F1', '#6807B7', '#78F358', '#A15DF7', '#511544', '#BFC1AF', '#2CFA45', '#8D2E0A', '#AA4658', '#5EE75C', '#711D25', '#149CA2', '#A90F54', '#0EA111', '#C14018', '#7BABF5', '#9287D4', '#969077', '#E9E628', '#632C88', '#9161A6', '#549B58', '#C38CA6', '#01FED3', '#B59929', '#28833F', '#064012', '#7ACDC6', '#84DD3D', '#8BDF23', '#CA5112', '#7E9FC3', '#540F61', '#E7D76A', '#880F50', '#D1CEA4', '#416609', '#1754A3', '#B3C5E6', '#B45B3A', '#54935C', '#981A43', '#C9A86B', '#610553', '#739A64', '#CA6354', '#DCE2EF', '#9F848C', '#6E3DFE', '#665164', '#7A30C5', '#851034', '#7078CB', '#8FAB86', '#BA5500', '#1235A0', '#34B204', '#7F937C', '#522CC9', '#DEE80D', '#D6D5FD', '#2B3918', '#581997', '#FE5878', '#0A9AAC', '#DDB1CB', '#2CE5DD', '#7E4F34', '#3B6824', '#B43761', '#3EEE83', '#11A600', '#616DFB', '#A7D121', '#D451FD', '#846595', '#A7438D', '#BAAF55', '#CDBAB9', '#249203', '#013D1B', '#8AF1F9', '#BEC070', '#E89E39', '#EFCE05', '#B54C5D', '#2A448F', '#2E705E', '#E312CD', '#B99450', '#98496D', '#D1BD8C', '#65E2B8', '#E67FFD', '#3D978A']

# plotting MECs
import smallestenclosingcircle
from collections import defaultdict
cluster_areas = defaultdict(lambda:0,cluster_areas)

for idx,triangle in enumerate(tri.simplices):
  clusterNum0 = find(triangle[0], clusters)
  clusterNum1 = find(triangle[1], clusters)
  clusterNum2 = find(triangle[2], clusters)
  if(clusterNum0 == clusterNum1 and clusterNum1 == clusterNum2 and not clusterNum0==-1 ):
    cluster_areas[clusterNum0] = cluster_areas[clusterNum0] + areas[idx]

for idx,cluster in enumerate(clusters):
  pointset=[]
  for point in cluster:
    pointset.append(points[point])
  mec = smallestenclosingcircle.make_circle(pointset)
  MECs.append(mec)
  radius = mec[2]  # radius of the circle
  mec_areas.append(math.pi*radius**2)

 # breaking the clusters with large MEC/area ratio
for idx,cluster in enumerate(clusters):
  if(cluster_areas[idx]!=0 and cluster_areas[idx]/mec_areas[idx]<breaking_ratio):
    endpts = []
    for point in cluster:
      if abs(distance(points[point], MECs[idx][:2])-MECs[idx][2])<0.000001:
        endpts.append(point) 
    if(len(endpts)==2):
      new_cluster=[]
      for pt in cluster:
        if(pt==endpts[0] or distance(points[pt], points[endpts[0]])<distance(points[pt], points[endpts[1]])):
          new_cluster.append(pt)
      for pt in new_cluster:
        cluster.remove(pt)
      clusters.append(new_cluster)
# plt.show()

new_mecs = [] 
# breaking large MEC
for idx,cluster in enumerate(clusters):
  pointset=[]
  for point in cluster:
    pointset.append(points[point])
  mec = smallestenclosingcircle.make_circle(pointset)
  new_mecs.append(mec)
  radius = mec[2]  # radius of the circle
  mec_areas.append(math.pi*radius**2)

# if any 2 new_mecs are overlapping a lot, merge them and add
# or else add them as it is
new_clusters=[]
for idx1, cluster1 in enumerate(clusters):
  for idx2, cluster2 in enumerate(clusters):
    if(idx1==idx2):
      continue
    circle1 = new_mecs[idx1]
    circle2 = new_mecs[idx2]  
    if(distance(circle1[:2], circle2[:2])+min(circle1[2],circle2[2])<merging_ratio*max(circle1[2],circle2[2])):
      found1 = find(cluster1[0], new_clusters) 
      found2 = find(cluster2[0], new_clusters)
      if(found1==-1 and found2==-1):
        new_clusters.append(cluster1+cluster2) 
      elif(found1==-1 and found2!=-1):
        new_clusters[found2] = new_clusters[found2]+cluster1
      elif(found1!=-1 and found2==-1):
        new_clusters[found1] = new_clusters[found1]+cluster2

for idx, cluster in enumerate(clusters):
  if notIn(cluster[0], new_clusters):
    new_clusters.append(cluster)

# after merging overlapping MEC
f = open("clusters.txt", "w")
plt.figure("After merging overlapping MECs")
plt.gca().set_aspect('equal')
for idx,cluster in enumerate(new_clusters):
  pointset=[]
  for point in cluster:
    f.write("["+str(points[point,0])+","+str(points[point,1])+"] ")
    pointset.append(points[point])
    plt.scatter(points[point,0], points[point,1], color=colors[idx%len(colors)])
  mec = smallestenclosingcircle.make_circle(pointset)
  f.write("["+str(mec[0])+","+str(mec[1])+","+str(mec[2])+"]\n")
  radius = mec[2]  # radius of the circle
  mec_areas.append(math.pi*radius**2)
  num_points = 100  # number of points on the circle
  angles = np.linspace(0, 2*np.pi, num_points)  # angles for each point
  x = mec[0] + radius * np.cos(angles)  # x-coordinates
  y = mec[1] + radius * np.sin(angles)  # y-coordinates
  plt.plot(x,y)
  
f.close()  
plt.show()
