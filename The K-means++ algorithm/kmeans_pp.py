import math
import numpy as np
import pandas as pd
import sys
import mykmeanssp as km


def initialize_centroids(data_points_np, k):
    np.random.seed(0)
    index_arr = np.empty((k))
    #Choose one center uniformly at random among the data points.
    N,d=data_points_np.shape
    rand = np.random.choice(N) #4654
    index_arr[0] = rand
    first_center = data_points_np[rand]
    centroides = np.empty((k,d-1))
    centroides[0] = first_center[1:]
    #the eclidean_distance array
    eclidean_distance = np.empty(N)#4465465
    for i in range(1, k):
        for j in range (0,N):
            data = data_points_np[j][1:]
            eclidean_distance[j] =compute_min_d(centroides,i,data)
        eclidean_distance /= np.sum(eclidean_distance)
        #Choose one new data point at random as a new center
        index = np.arange(N)  # Array containing numbers from 0 to 10
        chosen_index = np.random.choice(index, p=eclidean_distance)
        index_arr[i] = chosen_index
        centroides[i] = data_points_np[chosen_index][1:]
    return (centroides,index_arr)


def compute_min_d(centroides,last,element):
    min = float('inf')
    for i in range (0,last):
        distance = np.linalg.norm(centroides[i] - element)
        if (distance<min):
            min = distance
    return min




def p(matrix):
    for i in range(0, len(matrix.astype(int))):
        print(matrix.astype(int)[i], end="")
        if i + 1 != len(matrix.astype(int)):
            print(",", end="")
    print()


def printMat(arr):
    i = 0
    while i < len(arr):
        j = 0
        while j < len(arr[0]):
            print(np.round(arr[i][j], 4), end="")
            if j + 1 != len(arr[0]):
                print(",", end="")
            j += 1
        print()
        i += 1


def kmeansplusplus(k, max_iterations, epsilon, file1, file2):
    # combine both input files by inner join using the first column in each file as a key
    df1 = pd.read_csv(file1, header=None)
    df2 = pd.read_csv(file2, header=None)

    # Rename the first column to 'key'
    df1.rename(columns={df1.columns[0]: "key"}, inplace=True)
    df2.rename(columns={df2.columns[0]: "key"}, inplace=True)

    data_points = pd.merge(df1, df2, on="key", how='inner')
    # sort the data points by the ’key’ in ascending order.
    data_points.sort_values("key", inplace=True)
    input_matrix = data_points.values
    input_array = input_matrix.tolist()
    data = []
    for point in input_array:
        data.append(point[1:])
    # n := number of lines of an input file = number of vectors = len(input_data).
    n = len(data)
    # d := number of columns of an input file.
    d = len(data[0])
    # Check if the data is correct
    if (k >= n) or (k <= 1) or (int(k) != k):
        print("Invalid number of clusters!")
        exit(1)
    if (max_iterations <= 1) or (max_iterations >= 1000) or (int(max_iterations) != max_iterations):
        print("Invalid maximum iteration!")
        exit(1)
    # centroids µ1, µ2, ..., µK ∈ R^d where 1 < K < N.
    centroids, index_arr = initialize_centroids(input_matrix, k)
    matrix = km.fit(k, max_iterations, epsilon, n, d, data, centroids.tolist())
    p(index_arr)
    printMat(np.array(matrix))


input_argv = sys.argv
input_argc = len(sys.argv)
try:
    k = int(input_argv[1])


except ValueError:
    print("Invalid number of clusters!")
#input :K = 3, iter = 100, eps=0.01 , file 1 file 2
try:
    if (input_argc == 6):
        iteration = int(input_argv[2])
        eps = float(input_argv[3])
        kmeansplusplus(k, iteration, eps, (input_argv[4]), (input_argv[5]))


    elif (input_argc == 5):
        eps = float(input_argv[2])
        kmeansplusplus(int(input_argv[1]), 300, eps, (input_argv[3]), (input_argv[4]))

    else:
        print("An Error Has Occured")

except ValueError:
    print("An Error Has Occured")
