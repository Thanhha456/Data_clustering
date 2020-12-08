"""
implement two methods for computing closest pairs and two methods for clustering data.
"""
from timeit import default_timer as timer
import random
import alg_cluster
import math
import matplotlib.pyplot as plt
import urllib.request as urllib2
import alg_clusters_matplotlib
# Code to load data tables
# URLs for cancer risk data tables of various sizes
# Numbers indicate number of counties in data table

DIRECTORY = "http://commondatastorage.googleapis.com/codeskulptor-assets/"
DATA_3108_URL = DIRECTORY + "data_clustering/unifiedCancerData_3108.csv"
DATA_896_URL = DIRECTORY + "data_clustering/unifiedCancerData_896.csv"
DATA_290_URL = DIRECTORY + "data_clustering/unifiedCancerData_290.csv"
DATA_111_URL = DIRECTORY + "data_clustering/unifiedCancerData_111.csv"


def load_data_table(data_url):
    """
    Import a table of county-based cancer risk data
    from a csv format file
    """
    data_file = urllib2.urlopen(data_url)
    data = data_file.read()
    data = data.decode()
    data_lines = data.split('\n')
    print("Loaded", len(data_lines), "data points")
    data_tokens = [line.split(',') for line in data_lines]
    return [[tokens[0], float(tokens[1]), float(tokens[2]), int(tokens[3]), float(tokens[4])]
            for tokens in data_tokens]

def slow_closest_pair(cluster_list):
    """
    Compute the distance between the closest pair of clusters in a list (slow)
    Input: cluster_list is the list of clusters
    Output: tuple of the form (dist, idx1, idx2) where the centers of the clusters
    cluster_list[idx1] and cluster_list[idx2] have minimum distance dist.
    This function should implement the brute-force closest pair method described in SlowClosestPair
    A set P of (≥2) points whose i-th point,pi, is a pair(xi,yi).
    :param cluster_list:
    :return:A tuple(d,i,j)where d is the smallest pairwise distance of points in P,
    and i,j are the indices of two points whose distance is d
    """
    num = len(cluster_list)
    dist, idx_i, idx_j = float("inf"), -1, -1

    if num <= 1:
        return dist, idx_i, idx_j
    elif num == 2:
        dist = cluster_list[0].distance(cluster_list[1])
        return dist, 0, 1
    else:
        dist_list = []
        for dummy_x in range(num):
            for dummy_y in range(dummy_x + 1, num):
                dist_xy = cluster_list[dummy_x].distance(cluster_list[dummy_y])
                dist_list.append((dist_xy, dummy_x, dummy_y))
                # dist, idx_i, idx_j = min((dist, idx_i, idx_j), (dist_xy, dummy_x, dummy_y))

        dist_min, x_min, y_min = min(dist_list)
        dist, idx_i, idx_j = min((dist, idx_i, idx_j), (dist_min, x_min, y_min))
        return dist, idx_i, idx_j


def fast_closest_pair(cluster_list):
    """
    Compute the distance between the closest pair of clusters in a list (fast)

    Input: cluster_list is list of clusters SORTED such that horizontal positions of their
    centers are in ascending order

    Output: tuple of the form (dist, idx1, idx2) where the centers of the clusters
    cluster_list[idx1] and cluster_list[idx2] have minimum distance dist.

     A set P of (≥2) points whose i-th point, pi, is a pair(xi,yi), sorted in
     non decreasing order of their horizontal (x) coordinates.
    :param cluster_list:
    :return:A tuple(d,i,j)where d is the smallest pairwise distance of the points in P,
     and i,j are the indices of two points whose distance is d.
    """

    cluster_list.sort(key=lambda cluster: cluster.horiz_center())
    num = len(cluster_list)
    # base case
    if num <= 3:
        return slow_closest_pair(cluster_list)

    else:
        mid = int(math.ceil((num / 2.)))
        left_set = cluster_list[:mid]  # left_set and right_set are sorted too
        right_set = cluster_list[mid:]

        dist_l, idx_xl, idx_yl = fast_closest_pair(left_set)
        dist_r, idx_xr, idx_yr = fast_closest_pair(right_set)
        dist, idx_i, idx_j = min((dist_l, idx_xl, idx_yl), (dist_r, idx_xr + mid, idx_yr + mid))

        # centre line of strip
        mid_strip = (cluster_list[mid - 1].horiz_center() + cluster_list[mid].horiz_center()) / 2.
        dist, idx_i, idx_j = min((dist, idx_i, idx_j), closest_pair_strip(cluster_list, mid_strip, dist))
    return dist, idx_i, idx_j


def closest_pair_strip(cluster_list, horiz_center, half_width):
    """
     Helper function to compute the closest pair of clusters in a vertical strip

    Input: cluster_list is a list of clusters produced by fast_closest_pair
    horiz_center is the horizontal position of the strip's vertical center line
    half_width is the half the width of the strip (i.e; the maximum horizontal distance
    that a cluster can lie from the center line)

    Output: tuple of the form (dist, idx1, idx2) where the centers of the clusters
    cluster_list[idx1] and cluster_list[idx2] lie in the strip and have minimum distance dist.

    A set P of points whose i-th point,pi, is a pair(xi,yi);mid and w, both of which are real numbers.
    :param cluster_list: 
    :param horiz_center: 
    :param half_width: 
    :return: A tuple(d,i,j) where d is the smallest pairwise distance of points in P
    whose horizontal (x) coordinates are within w from mid.
    """
    # Let S be a list of the set{i:|xi−mid|< w}
    set_strip = []
    for cluster in cluster_list:
        if abs(cluster.horiz_center() - horiz_center) < half_width:
            set_strip.append(cluster)
    set_strip.sort(key=lambda cluster: cluster.vert_center())
    num = len(set_strip)
    dist, idx_i, idx_j = float("inf"), -1, -1
    if num <= 1:
        return float("inf"), -1, -1
    elif num < 3:
        for idx_u in range(num - 1):
            dist_u = set_strip[idx_u].distance(set_strip[idx_u + 1])
            idx_i = cluster_list.index(set_strip[idx_u])
            idx_j = cluster_list.index(set_strip[idx_u + 1])
            dist, idx_i, idx_j = min((dist, idx_i, idx_j), (dist_u, min(idx_i, idx_j), max(idx_i, idx_j)))
        return dist, idx_i, idx_j

    else:
        for idx_u in range(num - 1):
            for idx_v in range(idx_u + 1, min(idx_u + 3, num - 1) + 1):
                dist_uv = set_strip[idx_u].distance(set_strip[idx_v])
                idx_i_uv = cluster_list.index(set_strip[idx_u])
                idx_j_uv = cluster_list.index(set_strip[idx_v])
                dist, idx_i, idx_j = min((dist, idx_i, idx_j),
                                         (dist_uv, min(idx_i_uv, idx_j_uv), max(idx_i_uv, idx_j_uv)))
        return dist, idx_i, idx_j


def hierarchical_clustering(cluster_list, num_clusters):
    """
    Compute a hierarchical clustering of a set of clusters
    Note: the function may mutate cluster_list

    Input: A set P of points whose i-th point,pi, is a pair(xi,yi);k, the desired number of clusters.
    Output: A set C of k clusters that provides a clustering of the points in P.
    """
    cluster_list_copy = list(cluster_list)

    if len(cluster_list) <= num_clusters:
        return cluster_list
    while len(cluster_list) > num_clusters:
        cluster_list_copy.sort(key=lambda cluster: cluster.horiz_center())
        dummy, cluster_i, cluster_j = fast_closest_pair(cluster_list)
        cluster_list[cluster_i].merge_clusters(cluster_list[cluster_j])
        cluster_list.remove(cluster_list[cluster_j])

    return cluster_list

######################################################################################################################

def kmeans_clustering(cluster_list, num_clusters, num_iterations):
    """
    Compute the k-means clustering of a set of clusters
    Note: the function may not mutate cluster_list

    Input: List of clusters, integers number of clusters and number of iterations
    Output: List of clusters whose length is num_clusters
    Input: A set P of points whose i-th point,pi, is a pair(xi,yi);k, the desired number of clusters;q, a number of iterations.
    Output: A set C of k clusters that provides a clustering of the points in P
    """
    # position initial clusters at the location of clusters with largest populations
    cluster_list_copy = sorted(cluster_list,
                             reverse = True,
                             key=lambda cluster: cluster.total_population())
    cluster_list_copy = cluster_list_copy[: num_clusters]
    cluster_cent = [(cluster.horiz_center(), cluster.vert_center()) for cluster in cluster_list_copy]
    result = []
    #clustering to k initial centers adjusting the centers after each iteration
    for dummy_q in range(num_iterations):
        #Initialize k empty sets C1,...,Ck
        k_clusters = []
        for dummy_k in range(num_clusters):
            k_clusters.append(alg_cluster.Cluster(set(), 0, 0, 0, 0))
        for idx_j in range(len(cluster_list)):
            # defining the closest k center and add the cluster to it
            dist_list = []
            for idx_k in range(num_clusters):
                center_x, center_y = cluster_cent[idx_k]
                dist = cluster_list[idx_j].distance(
                        alg_cluster.Cluster(set(), center_x, center_y, 0, 0))
                dist_list.append((dist, idx_k))
            dummy_k, idx = min(dist_list)
            k_clusters[idx].merge_clusters(cluster_list[idx_j])
        result = k_clusters
        #update the new center of k clusters
        cluster_cent = [(k_clusters[idx_f].horiz_center(), k_clusters[idx_f].vert_center()) for idx_f in range(num_clusters)]
    return result

##################################################################################################################################

def  gen_random_clusters(num_clusters):
     """
     Create a list of clusters where each cluster in this list corresponds to one randomly generated point in the square with corners (±1,±1).

     """
     cluster_list = []
     for dummy in range(num_clusters):
         x_coordinate = random.choice([-1, 1]) * random.random()
         y_coordinate = random.choice([-1, 1]) * random.random()
         cluster_list.append(alg_cluster.Cluster(set(), x_coordinate, y_coordinate, 0, 0))
     return cluster_list

def plot_running_time(num_clusters):
    """
    Compute the running times of the functions slow_closest_pair and fast_closest_pair for lists of clusters of size 2 to 200.
     The horizontal axis for your plot should be the the number of initial clusters 
     while the vertical axis should be the running time of the function in seconds.
    :param num_clusters: 
    :return: 
    """
    slow_running = []
    fast_running = []
    for dummy_i in range(2, num_clusters):
         cluster_list = gen_random_clusters(dummy_i)
         start = timer()
         fast_closest_pair(cluster_list)
         end = timer()
         fast_running.append((end - start))
         
         start = timer()
         slow_closest_pair(cluster_list)
         end = timer()
         slow_running.append((end - start))
 #
    plt.plot(range(2, num_clusters), fast_running)
    plt.plot(range(2, num_clusters), slow_running)
    plt.xlabel("num clusters")
    plt.ylabel("running time in seconds")
    plt.title("Running time slow closest pair vs fast closest pair.")
    plt.legend(["fast closest pair", "slow closest pair"])
    plt.show()

#######################################################################################################
#Given a cluster C, its error is the sum of the squares of the distances from each county in the cluster
# to the cluster's center, weighted by each county's population. If pi is the position of the county and wi
# is its population, the cluster's error is:
#    error(C)=∑pi∈Cwi(dpic)2
#where c is the center of the cluster C.
#Given a list of clusters L, the distortion of the clustering is the sum of the errors associated with its clusters.

def compute_distortion(data_table, cluster_list):
    """
    take a list of clusters and use cluster_error to compute its distortion of two clustering:
    hierarchical clustering and  k-means clustering of the 111 counties k = 9, m = 5
    :param cluster_list:
    :return:
    """
    distortion = 0
    for cluster in cluster_list:
        distortion += cluster.cluster_error(data_table)
    return distortion

def distortion_of_hierarchical_clustering(data_table):
    """
    Compute the distortion of the list of clusters produced by hierarchical clustering on the 111, 290, and 896 county data sets,
    respectively, where the number of output clusters ranges from 6 to 20 (inclusive)
    You should remember that you can use the hierarchical cluster of size 20 to compute the hierarchical clustering of size 19 and so on
    :param data_table:
    :return:
    """
    max_range = 20
    min_range = 6
    singleton_list = []
    for line in data_table:
        singleton_list.append(alg_cluster.Cluster(set([line[0]]), line[1], line[2], line[3], line[4]))
    if len(singleton_list) < max_range:
        return None
    cluster_list = hierarchical_clustering(singleton_list, max_range)
    distortion = compute_distortion(data_table, cluster_list)
    distortion_list = [distortion]
    num = max_range
    while num >= min_range + 1:
        cluster_list = hierarchical_clustering(singleton_list, num)
        cluster_list.sort(key=lambda cluster: cluster.horiz_center())
        dummy, idx_x, idx_y = fast_closest_pair(cluster_list)
        cluster_list[idx_x].merge_clusters(cluster_list[idx_y])
        distortion += cluster_list[idx_x].cluster_error(data_table)
        distortion_list.append(distortion)
        num -= 1
    return distortion_list

def distortion_of_kmeans_clustering(data_table):
    """
    Compute the distortion of the list of clusters produced by kmeans clustering on the 111, 290, and 896 county data sets,
    respectively, where the number of output clusters ranges from 6 to 20 (inclusive)
    :param data_table:
    :return:
    """
    num_iritations = 5
    singleton_list = []
    for line in data_table:
        singleton_list.append(alg_cluster.Cluster(set([line[0]]), line[1], line[2], line[3], line[4]))
    distortion_list = []
    for num in range(20, 5, -1):
        cluster_list = kmeans_clustering(singleton_list,num, num_iritations)
        distortion = compute_distortion(data_table, cluster_list)
        distortion_list.append(distortion)
    return distortion_list

#####################################################################
# Code to load cancer data, compute a clustering and
# visualize the results


# def run_example():
#     """
#     Load a data table, compute a list of clusters and
#     plot a list of clusters
#
#     Set DESKTOP = True/False to use either matplotlib or simplegui
#     """
#     data_table = load_data_table(DATA_3108_URL)
#     singleton_list = []
#     for line in data_table:
#         singleton_list.append(alg_cluster.Cluster(set([line[0]]), line[1], line[2], line[3], line[4]))
    num_clusters = 16
    # cluster_list = sequential_clustering(singleton_list, num_clusters)
    # print("Displaying", len(cluster_list), "sequential clusters")
    #
    # cluster_list = alg_project3_solution.hierarchical_clustering(singleton_list, num_clusters)
    # print("Displaying", len(cluster_list), "hierarchical clusters")
    #
    # cluster_list = alg_project3_solution.kmeans_clustering(singleton_list, num_clusters, 5)
    # print("Displaying", len(cluster_list), "k-means clusters")

    # draw the clusters using matplotlib or simplegui
    #
    # if DESKTOP:
    #     # alg_clusters_matplotlib.plot_clusters(data_table, cluster_list, False)
    #      alg_clusters_matplotlib.plot_clusters(data_table, cluster_list, True)  #add cluster centers

    # else:
    #     alg_clusters_simplegui.PlotClusters(data_table, cluster_list)  # use toggle in GUI to add cluster centers

# run_example()
# data_table = load_data_table(DATA_111_URL)
# print(distortion_of_hierarchical_clustering(data_table))
# print(distortion_of_kmeans_clustering(data_table))
# data_table1 = load_data_table(DATA_111_URL)
# data_table2 = load_data_table(DATA_290_URL)
# data_table3 = load_data_table(DATA_896_URL)
# distortions_y1 = distortion_of_hierarchical_clustering(data_table1)
# distortions_y2 = distortion_of_hierarchical_clustering(data_table2)
# distortions_y3 = distortion_of_hierarchical_clustering(data_table3)
# distortions_x = range(20, 5, -1)
# distortions_z1 = distortion_of_kmeans_clustering(data_table1)
# distortions_z2 = distortion_of_kmeans_clustering(data_table2)
# distortions_z3 = distortion_of_kmeans_clustering(data_table3)
# #
# plt.plot(distortions_x, distortions_y1)
# plt.plot(distortions_x, distortions_z1)
# plt.title("Distortions of the hierarchical_clustering vs kmeans_clustering method on DATA_111_URL")

# # #
# plt.plot(distortions_x, distortions_y2)
# plt.plot(distortions_x, distortions_z2)
# plt.title("Distortions of the hierarchical_clustering vs kmeans_clustering method on DATA_290_URL")
# #
# plt.plot(distortions_x, distortions_y3)
# plt.plot(distortions_x, distortions_z3)
# plt.title("Distortions of the hierarchical_clustering vs kmeans_clustering method on DATA_896_URL")
# #
# plt.xlabel("number of clusters")
# plt.ylabel("distortions")
# plt.legend(["hierarchical clustering","kmeans clustering"])
# plt.show()

# plot_running_time(200)


# print(gen_random_clusters(10))
# def run_example():
#     """
#     Load a data table, compute a list of clusters
#     """
#     cluster_list = load_data_table(DATA_111_URL)
# print(fast_closest_pair([alg_cluster.Cluster(set([]), 0, 0, 1, 0), alg_cluster.Cluster(set([]), 1, 0, 1, 0)]))
# print(closest_pair_strip([alg_cluster.Cluster(set([]), 0, 0, 1, 0), alg_cluster.Cluster(set([]), 1, 0, 1, 0), alg_cluster.Cluster(set([]), 2, 0, 1, 0), alg_cluster.Cluster(set([]), 3, 0, 1, 0)], 1.5, 1.0))
# cluster_list = [alg_cluster.Cluster(set([0]), 3, 4, 1, 1), alg_cluster.Cluster(set([0]), 5,1, 1, 1),alg_cluster.Cluster(set([0]), 7, 9, 1, 1),alg_cluster.Cluster(set([0]), 5, 8, 1, 1) ]
# # print(hierarchical_clustering(cluster_list, 2))
# print(fast_closest_pair([alg_cluster.Cluster(set([]), 0.11, 0.75, 1, 0), alg_cluster.Cluster(set([]), 0.62, 0.86, 1, 0), alg_cluster.Cluster(set([]), 0.65, 0.68, 1, 0), alg_cluster.Cluster(set([]), 0.68, 0.48, 1, 0), alg_cluster.Cluster(set([]), 0.7, 0.9, 1, 0), alg_cluster.Cluster(set([]), 0.79, 0.18, 1, 0)]))
# print(slow_closest_pair([alg_cluster.Cluster(set([]), 0.11, 0.75, 1, 0), alg_cluster.Cluster(set([]), 0.62, 0.86, 1, 0), alg_cluster.Cluster(set([]), 0.65, 0.68, 1, 0), alg_cluster.Cluster(set([]), 0.68, 0.48, 1, 0), alg_cluster.Cluster(set([]), 0.7, 0.9, 1, 0), alg_cluster.Cluster(set([]), 0.79, 0.18, 1, 0)]))
# print(kmeans_clustering([alg_cluster.Cluster(set([1]), 0.11, 0.75, 1, 1), alg_cluster.Cluster(set([2]), 0.62, 0.86, 1, 3), alg_cluster.Cluster(set([3]), 0.65, 0.68, 1, 5), alg_cluster.Cluster(set([4]), 0.68, 0.48, 1, 4), alg_cluster.Cluster(set([5]), 0.7, 0.9, 1, 2), alg_cluster.Cluster(set([6]), 0.79, 0.18, 1, 7)], 3, 2))
# print(hierarchical_clustering([alg_cluster.Cluster(set([1]), 0.11, 0.75, 1, 0), alg_cluster.Cluster(set([2]), 0.62, 0.86, 1, 0), alg_cluster.Cluster(set([3]), 0.65, 0.68, 1, 0), alg_cluster.Cluster(set([4]), 0.68, 0.48, 1, 0), alg_cluster.Cluster(set([5]), 0.7, 0.9, 1, 0), alg_cluster.Cluster(set([6]), 0.79, 0.18, 1, 0)],3))
