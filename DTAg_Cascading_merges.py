from collections import defaultdict
from collections import Counter
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import math
import random
import smallestenclosingcircle


def largest_jump_value(sorted_data):
    # if len(sorted_data) < 2:
    #     return None  # There are not enough data points to calculate a jump
    # largest_jump = (sorted_data[1] - sorted_data[0])
    # jump_value = sorted_data[1]

    # for i in range(round(len(sorted_data)*0.1), round(len(sorted_data)*0.8)):
    #     diff = (sorted_data[i + 1] - sorted_data[i])
    #     if diff > largest_jump:
    #         largest_jump = diff
    #         jump_value = sorted_data[i]
    # if (equispaced == 0):
    #     fraction = 0.08
    # else:
    #     fraction = 0.05
    fraction = spacingLevel*0.05
    sorted_data = list(set(sorted_data))
    sorted_data.sort()
    return sorted_data[round(len(sorted_data)*fraction)]


def Area(vertex1, vertex2, vertex3):
    x1, y1 = vertex1
    x2, y2 = vertex2
    x3, y3 = vertex3
    area = abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2)
    return area


def MaxSide(vertex1, vertex2, vertex3):
    return (max(distance(vertex1, vertex2), distance(vertex1, vertex3), distance(vertex2, vertex3)))


def find(point, clusterSet):
    clusterNum = -1
    counter = 0
    for cluster in clusterSet:
        if point in cluster:
            clusterNum = counter
            break
        counter = counter+1
    return clusterNum


def notIn(point, clusterSet):
    return (find(point, clusterSet) == -1)


def Shortest_edge(x1, y1, x2, y2, x3, y3):
    # Calculate the lengths of the sides
    a = math.sqrt((x2 - x3)**2 + (y2 - y3)**2)
    b = math.sqrt((x1 - x3)**2 + (y1 - y3)**2)
    c = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)
    if (a <= b and a <= c):
        ord = [0, 1, 2]
    elif (b <= c and b <= a):
        ord = [1, 2, 0]
    else:
        ord = [2, 1, 0]
    return (max(a, b, c), min(a, b, c), ord)


def sortingKey(simplex):
    return (-1*MaxSide(points[simplex[0]], points[simplex[1]], points[simplex[2]]))


def distance(pointA, pointB):
    return math.sqrt(math.pow(pointA[0]-pointB[0], 2)+math.pow(pointA[1]-pointB[1], 2))


# def mergingStep(clusters, MECs, merging_ratio):
#     new_clusters = []
#     for idx1, cluster1 in enumerate(clusters):
#         circle1 = MECs[idx1]
#         best_overlap_ratio = -1
#         best_overlap_idx = -1
#         best_c2c = -1
#         for idx2, cluster2 in enumerate(clusters):
#             if (idx1 == idx2 or MECs[idx2][2] == 0):
#                 continue
#             circle2 = MECs[idx2]
#             overlap_ratio = (distance(circle1[:2], circle2[:2])+min(
#                 circle1[2], circle2[2]))/max(circle1[2], circle2[2]) if circle1[2] != 0 else -1
#             if (overlap_ratio != -1 and overlap_ratio < merging_ratio):
#                 if (overlap_ratio < best_overlap_ratio or best_overlap_ratio == -1):
#                     best_overlap_ratio = overlap_ratio
#                     best_overlap_idx = idx2
#             elif (overlap_ratio == -1):
#                 if (distance(circle1[:2], circle2[:2]) < best_c2c or best_c2c == -1):
#                     best_c2c = distance(circle1[:2], circle2[:2])
#                     best_overlap_idx = idx2
#         if (best_overlap_idx != -1):
#             found1 = find(cluster1[0], new_clusters)
#             found2 = find(clusters[best_overlap_idx][0], new_clusters)
#             if (found1 == -1 and found2 == -1):
#                 new_clusters.append(cluster1+clusters[best_overlap_idx])
#             elif (found1 == -1 and found2 != -1):
#                 new_clusters[found2] = new_clusters[found2]+cluster1
#             elif (found1 != -1 and found2 == -1):
#                 new_clusters[found1] = new_clusters[found1] + \
#                     clusters[best_overlap_idx]
#     for cluster in clusters:
#         if notIn(cluster[0], new_clusters):
#             new_clusters.append(cluster)
#     return new_clusters

def mergingStep(clusters, MECs, merging_ratio):
    new_clusters = []
    for idx1, cluster1 in enumerate(clusters):
        circle1 = MECs[idx1]
        found1 = find(cluster1[0], new_clusters)
        tbm = []
        for idx2, cluster2 in enumerate(clusters):
            if (idx1 == idx2):
                continue
            found2 = find(cluster2[0], new_clusters)
            circle2 = MECs[idx2]
            if (found2 == -1 and distance(circle1[:2], circle2[:2]) < max(circle1[2], circle2[2])):
                tbm.append(idx2)
        premerged = []
        for idx in tbm:
            premerged = premerged+clusters[idx]
        if (found1 == -1):
            new_clusters.append(premerged+cluster1)
        elif (found1 != -1):
            new_clusters[found1] = new_clusters[found1] + premerged
    for cluster in clusters:
        if notIn(cluster[0], new_clusters):
            new_clusters.append(cluster)
    print(len(new_clusters))
    return new_clusters


def mergingStep2(clusters, MECs, merging_ratio):
    new_clusters = []
    for idx1, cluster1 in enumerate(clusters):
        circle1 = MECs[idx1]
        found1 = find(cluster1[0], new_clusters)
        if (circle1[2] == 0):
            continue
        tbm = []
        for idx2, cluster2 in enumerate(clusters):
            if (idx1 == idx2):
                continue
            found2 = find(cluster2[0], new_clusters)
            circle2 = MECs[idx2]
            if (found2 == -1 and distance(circle1[:2], circle2[:2]) < max(circle1[2], circle2[2])):
                tbm.append(idx2)
            else:
                if (found2 == -1 and circle2[2] == 0 and distance(circle1[:2], circle2[:2]) < merging_ratio*2*circle1[2]):
                    tbm.append(idx2)
                elif (found2 == -1 and circle2[2] > 0 and (distance(circle1[:2], circle2[:2])/(circle1[2]+circle2[2])) < merging_ratio):
                    tbm.append(idx2)
        premerged = []
        for idx in tbm:
            premerged = premerged+clusters[idx]
        if (found1 == -1):
            new_clusters.append(premerged+cluster1)
        if (found1 != -1):
            new_clusters[found1] = new_clusters[found1] + premerged
    for cluster in clusters:
        if notIn(cluster[0], new_clusters):
            new_clusters.append(cluster)
    print(len(new_clusters))
    return new_clusters


def computeMEC(clusters):
    MECs = []
    mec_areas = []
    for idx, cluster in enumerate(clusters):
        pointset = []
        for point in cluster:
            pointset.append(points[point])
        mec = smallestenclosingcircle.make_circle(pointset)
        MECs.append(mec)
        radius = mec[2]  # radius of the circle
        mec_areas.append(math.pi*radius**2)
    return (MECs, mec_areas)


def displayMEC(clusters, MECs, title):
    plt.figure(title)
    plt.gca().set_aspect('equal')
    for idx, cluster in enumerate(clusters):
        for point in cluster:
            plt.scatter(points[point, 0], points[point, 1],
                        color=colors[idx % len(colors)])
        radius = MECs[idx][2]  # radius of the circle
        num_points = 100  # number of points on the circle
        angles = np.linspace(0, 2*np.pi, num_points)  # angles for each point
        x = MECs[idx][0] + radius * np.cos(angles)  # x-coordinates
        y = MECs[idx][1] + radius * np.sin(angles)  # y-coordinates
        plt.plot(x, y)


def duplicateCheck(data):
    # Flatten the list of lists into a single list
    flat_data = [num for sublist in data for num in sublist]

    # Count the occurrences of each number
    number_counts = Counter(flat_data)

    # Find the numbers that are repeated
    repeated_numbers = [num for num,
                        count in number_counts.items() if count > 1]

    print("Repeated numbers:", repeated_numbers)
    ###


def DT_Clustering(points, sl_ratio, breaking_ratio, merging_ratio, spacingLevel):
    tri = Delaunay(points)
    clusters = []
    areas = []
    maxSideLengths = []
    cluster_areas = {}

    # sort the triangles by increasing side length
    tri.simplices = sorted(tri.simplices, key=sortingKey)

    for triangle in tri.simplices:
        areas.append(Area(points[triangle[0]],
                          points[triangle[1]], points[triangle[2]]))
        maxSideLengths.append(
            MaxSide(points[triangle[0]], points[triangle[1]], points[triangle[2]]))

    # calculating average length of longest 30% of triangles
    # for i in range(0, math.ceil(0.3*len(maxSideLengths))):
    #     total_length = total_length+maxSideLengths[i]

    # len_lim = total_length/(0.3*len(maxSideLengths))
    sl_ratio = 0.2

    len_lim = largest_jump_value(maxSideLengths)

    plt.figure("long sides")
    for idx, triangle in enumerate(tri.simplices[:]):
        if (maxSideLengths[idx] > len_lim):
            plt.plot([points[triangle[0], 0], points[triangle[1], 0], points[triangle[2], 0]], [
                points[triangle[0], 1], points[triangle[1], 1], points[triangle[2], 1]], color="black")
            longest_edge, shortest_edge, order = Shortest_edge(
                points[triangle[0], 0], points[triangle[0], 1], points[triangle[1], 0], points[triangle[1], 1], points[triangle[2], 0], points[triangle[2], 1])
            clusterNum0 = find(triangle[order[0]], clusters)
            clusterNum1 = find(triangle[order[1]], clusters)
            clusterNum2 = find(triangle[order[2]], clusters)
            if (clusterNum0 == -1):
                clusters.append([triangle[order[0]]])
            if (shortest_edge/longest_edge < sl_ratio):
                if (clusterNum1 != -1 and clusterNum2 != -1):
                    continue
                if (clusterNum1 != -1 and clusterNum2 == -1):
                    clusters[clusterNum1].append(triangle[order[2]])
                elif (clusterNum2 != -1 and clusterNum1 == -1):
                    clusters[clusterNum2].append(triangle[order[1]])
                else:
                    clusters.append([triangle[order[1]], triangle[order[2]]])
            else:
                if (clusterNum1 == -1):
                    clusters.append([triangle[order[1]]])
                if (clusterNum2 == -1):
                    clusters.append([triangle[order[2]]])
        else:
            plt.plot([points[triangle[0], 0], points[triangle[1], 0], points[triangle[2], 0]], [
                points[triangle[0], 1], points[triangle[1], 1], points[triangle[2], 1]], color="red")
            clusterNum0 = find(triangle[0], clusters)
            clusterNum1 = find(triangle[1], clusters)
            clusterNum2 = find(triangle[2], clusters)

            # if all three vertices are not added to any cluster yet
            if clusterNum0 == -1 and clusterNum1 == -1 and clusterNum2 == -1:
                clusters.append([triangle[0], triangle[1], triangle[2]])
            # if all three vertices are added to cluster already
            elif (clusterNum0 != -1 and clusterNum1 != -1 and clusterNum2 != -1):
                continue
            # any other case
            else:
                # if point 1 already added, add whichever point not added to cluster
                if (clusterNum0 != -1):
                    if (clusterNum1 == -1):
                        clusters[clusterNum0].append(triangle[1])
                    if (clusterNum2 == -1):
                        clusters[clusterNum0].append(triangle[2])
                # if point 1 not yet added
                else:
                    # if point 2 not yet added, then point 3 must have already been added
                    if (clusterNum1 == -1 and clusterNum2 != -1):
                        clusters[clusterNum2].append(triangle[0])
                        clusters[clusterNum2].append(triangle[1])
                    # if point 3 not yet added, then point 2 must have been added
                    elif (clusterNum2 == -1 and clusterNum1 != -1):
                        clusters[clusterNum1].append(triangle[0])
                        clusters[clusterNum1].append(triangle[2])
                    # if point 2 and 3 both have been added
                    else:
                        if clusterNum1 == clusterNum2:
                            clusters[clusterNum2].append(triangle[0])
                        else:
                            clusters.append(clusters[clusterNum1] +
                                            clusters[clusterNum2]+[triangle[0]])
                            clusters.pop(clusterNum1)
                            clusterNum2 = find(triangle[2], clusters)
                            clusters.pop(clusterNum2)

    # colors=['#4C2B9F', '#78290C', '#CA223D', '#C0A364', '#2C9645', '#C785F1', '#6807B7', '#78F358', '#A15DF7', '#511544', '#BFC1AF', '#2CFA45', '#8D2E0A', '#AA4658', '#5EE75C', '#711D25', '#149CA2', '#A90F54', '#0EA111', '#C14018', '#7BABF5', '#9287D4', '#969077', '#E9E628', '#632C88', '#9161A6', '#549B58', '#C38CA6', '#01FED3', '#B59929', '#28833F', '#064012', '#7ACDC6', '#84DD3D', '#8BDF23', '#CA5112', '#7E9FC3', '#540F61', '#E7D76A', '#880F50', '#D1CEA4', '#416609', '#1754A3', '#B3C5E6', '#B45B3A', '#54935C', '#981A43', '#C9A86B', '#610553', '#739A64', '#CA6354', '#DCE2EF', '#9F848C', '#6E3DFE', '#665164', '#7A30C5', '#851034', '#7078CB', '#8FAB86', '#BA5500', '#1235A0', '#34B204', '#7F937C', '#522CC9', '#DEE80D', '#D6D5FD', '#2B3918', '#581997', '#FE5878', '#0A9AAC', '#DDB1CB', '#2CE5DD', '#7E4F34', '#3B6824', '#B43761', '#3EEE83', '#11A600', '#616DFB', '#A7D121', '#D451FD', '#846595', '#A7438D', '#BAAF55', '#CDBAB9', '#249203', '#013D1B', '#8AF1F9', '#BEC070', '#E89E39', '#EFCE05', '#B54C5D', '#2A448F', '#2E705E', '#E312CD', '#B99450', '#98496D', '#D1BD8C', '#65E2B8', '#E67FFD', '#3D978A']

    cluster_areas = defaultdict(lambda: 0, cluster_areas)

    for idx, triangle in enumerate(tri.simplices):
        clusterNum0 = find(triangle[0], clusters)
        clusterNum1 = find(triangle[1], clusters)
        clusterNum2 = find(triangle[2], clusters)
        if (clusterNum0 == clusterNum1 and clusterNum1 == clusterNum2 and not clusterNum0 == -1):
            cluster_areas[clusterNum0] = cluster_areas[clusterNum0] + areas[idx]

    MECs, mec_areas = computeMEC(clusters)
    displayMEC(clusters, MECs, "Initial clustering")
    print(clusters)
    duplicateCheck(clusters)

    # # breaking the clusters with large MEC/area ratio
    # for idx, cluster in enumerate(clusters):
    #     if (cluster_areas[idx] != 0 and cluster_areas[idx]/mec_areas[idx] < breaking_ratio):
    #         endpts = []
    #         for point in cluster:
    #             if abs(distance(points[point], MECs[idx][:2])-MECs[idx][2]) < 0.000001:
    #                 endpts.append(point)
    #         if (len(endpts) == 2):
    #             new_cluster = []
    #             for pt in cluster:
    #                 if (pt == endpts[0] or distance(points[pt], points[endpts[0]]) < distance(points[pt], points[endpts[1]])):
    #                     new_cluster.append(pt)
    #             for pt in new_cluster:
    #                 cluster.remove(pt)
    #             clusters.append(new_cluster)
    ###

    # MECs, mec_areas = computeMEC(clusters)
    # plt.figure("Post breakdown")
    # plt.gca().set_aspect('equal')
    # displayMEC(clusters, MECs)

    # if any 2 new_mecs are overlapping a lot, merge them and add
    # or else add them as it is
    merge_counter = 0
    while (True):
        merge_counter += 1
        new_clusters = mergingStep2(clusters, MECs, merging_ratio)
        if (len(new_clusters) == len(clusters)):
            displayMEC(clusters, MECs, f'Post Merging{merge_counter}')
            break
        clusters = new_clusters
        MECs, mec_areas = computeMEC(clusters)
    plt.show()
    result = []
    for idx, point in enumerate(points):
        result.append([point[0], point[1], find(idx, clusters)])

    return result


# # -------------------------------------#
# # Random seed for reproducibility
# np.random.seed(90)

# # Number of clusters and points in each cluster
# num_clusters = random.randint(6, 10)
# num_clusters = 14
# points_per_cluster = np.random.randint(10, 16, num_clusters)

# # Generate the dataset
# dataset = []
# for i in range(num_clusters):
#     center_x = np.random.uniform(-20, 20)
#     center_y = np.random.uniform(-20, 20)
#     cluster_points = []
#     for _ in range(points_per_cluster[i]):
#         point_x = round(np.random.normal(center_x, 1.0), 2)
#         point_y = round(np.random.normal(center_y, 1.0), 2)
#         cluster_points.append([point_x, point_y])
#     dataset.extend(cluster_points)

# # Convert the dataset to a numpy array
# points = np.array(dataset)
# # -------------------------------------#


# ---------------------------------#
csv_file_path = 'Testcase16.csv'

# Initialize an empty list to store the points
points = []

# Open the CSV file and read its contents
with open(csv_file_path, 'r') as csv_file:
    for line in csv_file:
        # Split each line using spaces as the delimiter
        values = line.strip().split("\t")
        # Assuming each line contains two values for x and y coordinates, change as needed
        x, y = map(float, values)  # Convert values to float
        points.append([x, y])

# Convert the list of points to a NumPy array
points = np.array(points)
# --------------------------------#

colors = ['#4C2B9F', '#78290C', '#CA223D', '#C0A364', '#2C9645', '#C785F1', '#6807B7', '#78F358', '#A15DF7', '#511544', '#BFC1AF', '#2CFA45', '#8D2E0A', '#AA4658', '#5EE75C', '#711D25', '#149CA2', '#A90F54', '#0EA111', '#C14018', '#7BABF5', '#9287D4', '#969077', '#E9E628', '#632C88', '#9161A6', '#549B58', '#C38CA6', '#01FED3', '#B59929', '#28833F', '#064012', '#7ACDC6', '#84DD3D', '#8BDF23', '#CA5112', '#7E9FC3', '#540F61', '#E7D76A', '#880F50', '#D1CEA4', '#416609', '#1754A3', '#B3C5E6', '#B45B3A', '#54935C', '#981A43', '#C9A86B', '#610553', '#739A64',
          '#CA6354', '#DCE2EF', '#9F848C', '#6E3DFE', '#665164', '#7A30C5', '#851034', '#7078CB', '#8FAB86', '#BA5500', '#1235A0', '#34B204', '#7F937C', '#522CC9', '#DEE80D', '#D6D5FD', '#2B3918', '#581997', '#FE5878', '#0A9AAC', '#DDB1CB', '#2CE5DD', '#7E4F34', '#3B6824', '#B43761', '#3EEE83', '#11A600', '#616DFB', '#A7D121', '#D451FD', '#846595', '#A7438D', '#BAAF55', '#CDBAB9', '#249203', '#013D1B', '#8AF1F9', '#BEC070', '#E89E39', '#EFCE05', '#B54C5D', '#2A448F', '#2E705E', '#E312CD', '#B99450', '#98496D', '#D1BD8C', '#65E2B8', '#E67FFD', '#3D978A']

######## TUNABLE PARAMETERS #########
sl_ratio = 0.3
breaking_ratio = 0.15
merging_ratio = 0.75
spacingLevel = 19.61 # higher spacingLevel => more points clustered together
plotMerges = 0
#####################################

result_data = DT_Clustering(points[:], sl_ratio, breaking_ratio,
                            merging_ratio, spacingLevel)

output_csv_file = 'outcome.csv'
with open(output_csv_file, 'w') as output_file:
    # Write the header
    # output_file.write("x,y,cluster_id\n")
    for row in result_data:
        output_file.write(f"{row[0]}, {row[1]}, {row[2]}\n")
