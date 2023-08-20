import sys

def kmeans(k ,num_of_iter ,file):
    try:
        num_from_file=[]
        clusters = []
        centroides = []
        with open(file) as file:
            for line in file:
                arr = line.rstrip().split(",")
                for i in range(0, len(arr)):
                    arr[i] = float(arr[i])
                num_from_file.append(arr)
        N = len(num_from_file)
        error = False
        if k<=1 or k>=N or k!=int(k):
            print("Invalid number of clusters!")
            error = True
        if num_of_iter<=1 or num_of_iter>=1000 or num_of_iter!=int(num_of_iter) :
            print("Invalid maximum iteration!")
            error = True
        if (error):
                return

        for i in range(0,k):
            centroides.append(num_from_file[i])
            
        initialize_clusters(clusters,k)
        curr_num_of_iter = 0
        while(curr_num_of_iter < num_of_iter):
            curr_num_of_iter += 1
            initialize_clusters(clusters,k)
            refill_cluster(num_from_file,centroides,clusters)
            stop = update_centroides(centroides,clusters)
            #print(" cluster = ",clusters ," centroides = ", centroides)
            if (stop):
                break
            
        #print centroides
        for i in range (0,len(centroides)):
            for j in range (0,len(centroides[i])):
                if (j < len(centroides[i])-1):
                    print("%.4f" %centroides[i][j],end=",")
                else:
                     print("%.4f" %centroides[i][j])
                     
        
    except ValueError :
        print("An Error Has Occured")

        
                
    
def initialize_clusters(clusters,k):
    clusters.clear()
    for i in range (0,k):
        clusters.append([])


def refill_cluster(num_from_file,centroides,clusters):
    positive_infinity = float('inf')
    for i in range (0 , len(num_from_file)):#iterate through the elemnts in num from file
        min_distance = positive_infinity
        index_of_min_distance = -1
        for j in range (0 , len(centroides)):#iterate through centoides
            curr_min = compute_distance(num_from_file[i],centroides[j])
            if (curr_min < min_distance):
                min_distance = curr_min
                index_of_min_distance = j
        clusters[index_of_min_distance].append(num_from_file[i])
    
    
def compute_distance(vector1 ,vector2):
    tmp = 0
    for i in range (0 , len(vector1)):
        tmp += pow((vector1[i] - vector2[i]),2)
    tmp = pow(tmp,0.5)
    return tmp

def update_centroides(centroides,clusters):
    stop = True
    for j in range (0 , len(centroides)):#iterate through centoides
        new_centroid =[]
        for t in range (0,len(centroides[j])): #iterate through coordinates
            tmp = 0
            for i in range (0,len(clusters[j])):#iterate through elemnts in cluster
                tmp += clusters[j][i][t]
            tmp = tmp/(len(clusters[j]))
            new_centroid.append(tmp)
        if ( compute_distance(new_centroid,centroides[j]) >= 0.001 ):
            stop = False
        centroides[j] = new_centroid
    return stop

input_argv=sys.argv
input_argc=len(sys.argv)
try:
    k=int(input_argv[1])
    
except ValueError :
        print("Invalid number of clusters!")
        
try:
    if (input_argc==4):
        iteration = int(input_argv[2])
        kmeans(k,iteration,(input_argv[3]))

    elif (input_argc==3):
        kmeans(int(input_argv[1]),200,(input_argv[2]))
        
    else:
        print("An Error Has Occured")
        
except ValueError :
        print("Invalid maximum iteration!")
