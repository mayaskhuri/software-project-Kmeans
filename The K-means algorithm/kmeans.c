#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double compute_distance(double* vector1, double* vector2, int d);
void refill_cluster(double** centroids, int** clusters, double** num_from_file, int K, int N, int d);
int update_centroides(int** clusters, double** centroids, double** num_from_file, int K, int d);
void append_to_cluster_j(int** clusters, int j, int index);
int** initialize_clusters(int K);
void kmeans(int K, int iter, char* input);

int main(int argc, char** argv)
{
    int k = strtol(argv[1], NULL, 10);
    if (argc >= 4)
    {
        int iter = strtol(argv[2], NULL, 10);
        kmeans(k, iter, argv[3]);
    }
    else
    {
        kmeans(k, 200, argv[2]);
    }
    return 0;

}

void kmeans(int K, int iter, char* input)
{
    double** centroids =NULL;
    double** num_from_file = NULL;
    int i=0;
    int j=0;
    int d = 0;
    int s = 0;
    int e = 0;
    int num_of_curr_iter = 0;
    FILE* file_p;
    char* line = NULL;
    ssize_t current;
    int stop=0;
    size_t length = 0;
    int num_lines = 0;
    int dimensions = 1;
    int N = 0;
    int** clusters=NULL;
    int invalid_k = 0, invalid_iter = 0;

    file_p = fopen(input, "r");
    if (file_p == NULL)
    {
        perror("An Error Has Occurred");
        exit(1);
    }

    while ((current = getline(&line, &length, file_p)) != -1)
    {
        num_lines++;
        if (num_lines == 1)
        {
            int j = 0;
            while (line[j] != '\0')
            {
                if (line[j] == ',')
                {
                    dimensions++;
                }
                j++;
            }
        }
    }

    
    N = num_lines;
    fseek(file_p, 0, SEEK_SET);

    if (K >= N || floor(K) != K || K < 1)
    {
        fprintf(stderr, "Invalid number of clusters!\n");
        invalid_k = 1;
    }
    if (iter >= 1000 || iter < 1 || floor(iter) != iter)
    {
        fprintf(stderr, "Invalid maximum iteration!\n");
        invalid_iter = 1;
    }
    if (invalid_k || invalid_iter)
    {
        exit(1);
    }
    fseek(file_p, 0, SEEK_SET);
    num_from_file= (double**)malloc(num_lines * sizeof(double*));
    if(num_from_file == NULL){
        perror("An Error Has Occurred");
        exit(1);
    }
    for ( i = 0; i < num_lines; i++)
    {
        double* arr = (double*)malloc(dimensions * sizeof(double));
        if(arr == NULL){
            perror("An Error Has Occurred");
            exit(1);
        }
        current = getline(&line, &length, file_p);
        s=0;
        e=0;
        j=0;
        while (e <= current)
        {
            if (line[e] == ',' || line[e] == '\0')
            {
                line[e] = '\0';
                arr[j++] = atof(&line[s]);
                s = e + 1;
            }
            
            e++;
        }
        num_from_file[i] = arr;
    }


    centroids= (double**)malloc(K * sizeof(double*));
    if(centroids == NULL){
        perror("An Error Has Occurred");
        exit(1);
    }
    for ( i = 0; i < K; i++)
    {
        centroids[i] = num_from_file[i];
    }

    d = dimensions;
    num_of_curr_iter = 0;
    
    clusters= initialize_clusters(K);
    while (num_of_curr_iter < iter)
    {  
        stop=0; 
        refill_cluster(centroids, clusters, num_from_file, K, N, dimensions);
        stop = update_centroides(clusters, centroids, num_from_file, K, dimensions);
        if (stop)
        {
            break;
        }
        num_of_curr_iter++;

        for ( i = 0; i < K; i++)
        {
            free(clusters[i]);
        }
        clusters = initialize_clusters(K);
        if(clusters == NULL){
            perror("An Error Has Occurred");
            exit(1);
    }
    }

    fclose(file_p);
    if (line)
    {
        free(line);
    }

    for ( i = 0; i < K; i++)
    {
        for ( j = 0; j < d; j++)
        {
            if(j<d-1)
                printf("%0.4f,", centroids[i][j]);
            else
                printf("%0.4f", centroids[i][j]);
        }
        printf("\n");
    }

    for ( i = 0; i < num_lines; i++)
    {
        free(num_from_file[i]);
    }
    free(num_from_file);
    free(centroids);
}

void append_to_cluster_j(int** clusters, int j, int index)
{
    int* tmp=NULL;
    clusters[j][0]++;
    tmp = (int*)realloc(clusters[j], (clusters[j][0] + 1) * sizeof(int));
    if (tmp != NULL)
    {
        clusters[j] = tmp;
        clusters[j][clusters[j][0]] = index;
    }
    else{
        perror("An Error Has Occurred");
        exit(1);
    }


}

int** initialize_clusters(int K)
{
    int** clusters;
    int i=0;
    int** tmp = (int**)calloc(K, sizeof(int*));
    int* tmp_1=NULL;
    if (tmp != NULL)
    {
        clusters = tmp;
        for ( i = 0; i < K; i++)
        {
             tmp_1 = (int*)calloc(1, sizeof(int));
            if (tmp_1 != NULL)
            {
                clusters[i] = tmp_1;
                clusters[i][0] = 0;
            }
            else{
                perror("An Error Has Occurred");
                exit(1);
    }
        }
    }
    else{
        perror("An Error Has Occurred");
        exit(1);
    }
    tmp_1[0]=0;
    tmp[0][0]=0;
    return clusters;
}

int update_centroides(int** clusters, double** centroids, double** num_from_file, int K, int d)
{
    int stop = 1;
    int m=0;
    int j=0,i=0;
    double* tmp=NULL;
    double* old_centroid = NULL;
    old_centroid=centroids[0];
    for ( i = 0; i < K; i++)
    {
        old_centroid = centroids[i];
        tmp = (double*)calloc(d, sizeof(double));
        if (tmp != NULL)
        {
            int num_of_elemnt_in_cluster = clusters[i][0];
            for ( j = 1; j <= num_of_elemnt_in_cluster; j++)
            {
                double* vector = num_from_file[clusters[i][j]];
                for ( m = 0; m < d; m++)
                {
                    tmp[m] += vector[m];
                }
            }
            for ( j = 0; j < d; j++)
            {
                tmp[j] = tmp[j] * (1.0 / num_of_elemnt_in_cluster);
            }
            if (compute_distance(tmp, old_centroid, d) >= 0.001)
            {
                stop = 0;
            }
            centroids[i] = tmp;
        }
        else{
            perror("An Error Has Occurred");
            exit(1);
         }
    }
    return stop;
}

void refill_cluster(double** centroids, int** clusters, double** num_from_file, int K, int N, int d)
{
    int i=0,j=0;
    for ( i = 0; i < N; i++)
    {
        double* vector = num_from_file[i];
        int cluster_index = 0;
        double min_distance = INFINITY;
        for ( j = 0; j < K; j++)
        {
            double distance = compute_distance(vector, centroids[j], d);
            if (distance < min_distance)
            {
                min_distance = distance;
                cluster_index = j;
            }
        }
        append_to_cluster_j(clusters, cluster_index, i);
    }
}

double compute_distance(double* vector1, double* vector2, int d)
{
    int i=0;
    double distance = 0.0;
    for ( i = 0; i < d; i++)
    {
        distance += pow(vector1[i] - vector2[i], 2);
    }
    return sqrt(distance);
}
