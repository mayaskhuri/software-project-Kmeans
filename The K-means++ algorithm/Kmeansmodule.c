#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "math.h"
#include <float.h>
#include <stdbool.h>
#include <math.h>
#include <Python.h>
#define PY_SSIZE_T_CLEAN

int find_closest_cluster_recursive(double *data_point, double **cluster_centers, int num_clusters, int num_dimensions);
double** create_empty_matrix(int r,int c);
double** initialize_clusters(double** mat,int k,int d);
double compute_distance(double* vector1, double* vector2, int d);
double*** create_empty(int k,int* points_num,int d);
double** update_centroids(double*** groups, int* num_points, int num_clusters, int num_dimensions);
double** kmeans(double** mat,double** curr_clusters,int n,int k,int d,int iter,double eps);
int validate_iter(int iter);
int validate_k(int k,int n);
void set_clusters(PyObject *ce, double **clusters, int k, int d) ;
void set_matrix(PyObject *inputt, double **mat, int n, int d) ;
double** allocate_matrix(int r, int c) ;
void update_points_counts(double** mat, double** curr_clusters, int* points_num, int* ind_vec, int n, int k, int d) ;
int validate_allocation(int* array);
void update_points_group(double** data_points, double** cluster_centers, int* cluster_indices, double*** point_groups, int num_points, int num_clusters, int dimensions) ;
int find_closest_cluster(double *data_point, double **cluster_centers, int num_clusters, int num_dimensions);
int iter;
int d;
int n;
int k;

static PyObject* fit(PyObject *self, PyObject *args){
    PyObject *inputt;
    PyObject *ce;
    double **clusters;
    double eps;
    PyObject *ret;
    double **mat=NULL;
    int i=0;
    int j=0;
    if (!PyArg_ParseTuple(args, "iidiiOO", &k, &iter, &eps, &n, &d, &inputt, &ce)) {
        return NULL;
    }
    if (!PyList_Check(inputt) || !PyList_Check(ce)){
        return NULL;
    }
    mat=create_empty_matrix(n,d);
    clusters = create_empty_matrix(k, d);
    ret = PyList_New(k);
    if (ret == NULL) {
        return NULL;
    }
    if (mat == NULL) {
        printf("An Error Occured");
        return NULL;
    }
    set_matrix(inputt,mat,n,d);
    if (clusters == NULL) {
        printf("An Error Occured");
        return NULL;
    }
    set_clusters(ce,clusters,k,d);
    clusters=kmeans(mat,clusters,n,k,d,iter,eps);
    while (i < k) {
        PyObject* l = PyList_New(d);
        if (l == NULL) {
            return NULL;
        }
        j = 0;
        while (j < d) {
            PyObject* item = PyFloat_FromDouble(clusters[i][j]);
            if (item != NULL) {
                PyList_SetItem(l, j, item);
                j++;
            }
            else{
                return NULL;
            }
        }
        if (PyList_SetItem(ret, i, l) == 0) {
            i++;
        }
        else{
            return NULL;
        }
    }
    for (i=0;i<n;i++){
        free(mat[i]);
    }
    free(mat);
    for (i=0;i<k;i++){
        free(clusters[i]);
    }
    free(clusters);
    return Py_BuildValue("O", ret);
}
static PyMethodDef myMethods[] = {
        { "fit",
        (PyCFunction)fit, METH_VARARGS, PyDoc_STR("a") },
        { NULL, NULL, 0, NULL }
};
static struct PyModuleDef mykmeanssp = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",
        "",
        -1,
        myMethods
};
PyMODINIT_FUNC PyInit_mykmeanssp(void) {
    PyObject *m;

    m = PyModule_Create(&mykmeanssp);
    if (!m) {
        return NULL;
    }
    return m;
}



void set_matrix(PyObject *inputt, double **mat, int n, int d) {
    int i = 0;
    while (i < n) {
        PyObject *l = PyList_GetItem(inputt, i);
        int j = 0;
        while (j < d) {
            PyObject *item = PyList_GetItem(l, j);
            double obj = PyFloat_AsDouble(item);
            mat[i][j] = obj;
            j++;
        }
        i++;
    }
}


void set_clusters(PyObject *ce, double **clusters, int k, int d) {
    int i, j;
    for (i = 0; i < k; i++) {
        j = 0;
        PyObject *r = PyList_GetItem(ce, i);
        while (j < d) {
            double obj = PyFloat_AsDouble(PyList_GetItem(r, j));
            clusters[i][j] = obj;
            j++;
        }
    }
}


double** create_empty_matrix(int r, int c) {
    double** matrix = (double**)malloc(r * sizeof(double*));
    if (matrix != NULL) {
        matrix=allocate_matrix(r,c);
    }
    else{
        return NULL;
    }
    return matrix;
}

double** allocate_matrix(int r, int c) {
    int i = 0, j;
    double** m = (double**)malloc(r * sizeof(double*));
    if (m != NULL) {
        while (i < r) {
        m[i] = (double*)malloc(c * sizeof(double));
        if (m[i] == NULL) {
            // Free previously allocated memory
            for (j = i - 1; j >= 0; j--) {
                free(m[j]);
            }
            free(m);
            return NULL;
        }
        
        i++;
    }
    }
    else{
        return NULL;
    }
    return m;
}


double** initialize_clusters(double** mat, int k, int d) {
    double** clusters = create_empty_matrix(k, d);
    if (clusters != NULL) {
        int i, j;
        i = 0;
        while (i < k) {
            for (j = 0; j < d; j++) {
                clusters[i][j] = mat[i][j];
            }
            i++;
        }
    }
    else{
        return NULL;
    }

    return clusters;
}

int validate_k(int k,int n){
    if (k>=n||k<=1 || k!=floor(k)){
        printf("Invalid number of clusters!");
        printf("\n");
        return 0;
    }
    return 1;
}

int validate_iter(int iter){
    if (iter<=1 || iter>=1000 || (iter!=floor(iter))){
        printf("Invalid maximum iteration!");

        printf("\n");
        return 0;
    }
    return 1;
}

int validate_allocation(int* array){
    if(array==NULL){
        printf("An Error Has Occurred");
        printf("\n");
        return 0;
    }
    return 1;
}

void update_points_counts(double **mat, double **curr_clusters, int *points_num, int *in, int n, int k, int d) {
  for (int x_i = 0; x_i < n; x_i++) {
    int min_index = find_closest_cluster(mat[x_i], curr_clusters, k, d);
    points_num[min_index]++;
    in[min_index]++;
  }
}

int find_closest_cluster_recursive(double *data_point, double **cluster_centers, int num_clusters, int num_dimensions) {
  if (num_clusters == 1) {
    return 0; 
  } else {
    double distance = compute_distance(data_point, cluster_centers[num_clusters - 1], num_dimensions);
    int closest_cluster = find_closest_cluster_recursive(data_point, cluster_centers, num_clusters - 1, num_dimensions);
    double closest_distance = compute_distance(data_point, cluster_centers[closest_cluster], num_dimensions);

    if (distance < closest_distance) {
      return num_clusters - 1; 
    } else {
      return closest_cluster;
    }
  }
}


void update_points_group(double **data_points, double **cluster_centers, int *cluster_indices, double ***point_groups, int num_points, int num_clusters, int num_dimensions) {
  int point_index, cluster_index;
  for (point_index = 0; point_index < num_points; point_index++) {
    cluster_index = find_closest_cluster(data_points[point_index], cluster_centers, num_clusters, num_dimensions);
    point_groups[cluster_index][cluster_indices[cluster_index] - 1] = data_points[point_index];
    cluster_indices[cluster_index]--;
  }
}

int find_closest_cluster(double *data_point, double **cluster_centers, int num_clusters, int num_dimensions) {
  if (num_clusters == 1) {
    return 0;  
  } else {
    double distance = compute_distance(data_point, cluster_centers[num_clusters - 1], num_dimensions);
    int closest_cluster = find_closest_cluster(data_point, cluster_centers, num_clusters - 1, num_dimensions);
    double closest_distance = compute_distance(data_point, cluster_centers[closest_cluster], num_dimensions);

    if (distance < closest_distance) {
      return num_clusters - 1; 
    } else {
      return closest_cluster; 
    }
  }
}


double** update_centroids(double*** groups, int* num_points, int num_clusters, int num_dimensions) {
    int i = 0, q, l;
    double average, sum;
    double** centroids;
    centroids = (double**)calloc(num_clusters, sizeof(double*));
    if (centroids != NULL) {
        while (i < num_clusters) {
        centroids[i] = (double*)calloc(num_dimensions, sizeof(double));
        if (centroids[i] == NULL) {
            printf("An Error Has Occurred\n");
            return NULL;
        }
        
        q = 0;
        while (q < num_dimensions) {
            sum = 0.0;
            l = 0;
            while (l < num_points[i]) {
                sum += groups[i][l][q];
                l++;
            }
            
            average = sum / num_points[i];
            centroids[i][q] = average;
            
            q++;
        }
        
        i++;
    }
    }
    else{
        printf("An Error Has Occurred\n");
        return NULL;

    }
    return centroids;
}


double*** create_empty(int k, int* counts, int d)
{
    int i = 0, j, n;
    double*** res = (double***)calloc(k, sizeof(double**));
    if (res == NULL) {
        return NULL;
    }
    
    while (i < k) {
        n = counts[i];
        res[i] = (double**)calloc(n, sizeof(double*));
        if (res[i] != NULL) {
            j = 0;
            while (j < n) {
                res[i][j] = (double*)calloc(d, sizeof(double));
                if (res[i][j] == NULL) {
                    return NULL;
                }
                j++;
            }
        }
        else {
            return NULL;
        }
        i++;
    }
    
    return res;
}

double compute_distance(double* vector1, double* vector2, int d)
{
    int i=0;
    double distance = 0.0;
    for( i = 0; i < d; i++)
    {
        distance += pow(vector1[i] - vector2[i], 2);
    }
    return sqrt(distance);
}


double** kmeans(double** data, double** initial_centroids, int numPoints, int numClusters, int numDimensions, int maxIterations, double epsilon)
{
    double** prevCentroids;
    int iterationCount;
    double*** pointGroups = NULL;
    int* pointCounts;
    int* indexArray;
    bool isFirstIteration;
    int i, j;
    
    if (validate_k(numClusters, numPoints) == 0) {
        return NULL;
    }
    
    if (validate_iter(maxIterations) == 0) {
        return NULL;
    }
    
    prevCentroids = create_empty_matrix(numClusters, numDimensions);
    
    if (prevCentroids == NULL) {
        printf("An Error Has Occurred\n");
        return NULL;
    }
    
    isFirstIteration = true;
    iterationCount = 0;
    pointCounts = (int*)calloc(numClusters, sizeof(int));
    if (validate_allocation(pointCounts) == 0) {
        return NULL;
    }
    
    indexArray = (int*)calloc(numClusters, sizeof(int));
    if (validate_allocation(indexArray) == 0) {
        return NULL;
    }
    
    for (iterationCount = 0; iterationCount < maxIterations && (!isFirstIteration || iterationCount == 0); iterationCount++) {
        isFirstIteration = false;
        prevCentroids = initial_centroids;
        update_points_counts(data, initial_centroids, pointCounts, indexArray, numPoints, numClusters, numDimensions);
        pointGroups = create_empty(numClusters, pointCounts, numDimensions);
        if (pointGroups == NULL) {
            printf("An Error Has Occurred\n");
            return NULL;
        }
        
        int pointIndex = 0;
        while (pointIndex < numPoints) {
            double minDistance = -1;
            int minIndex = 0;
            
            for (int clusterIndex = 0; clusterIndex < numClusters; clusterIndex++) {
                double currDistance = compute_distance(data[pointIndex], initial_centroids[clusterIndex], numDimensions);
                if ((currDistance < minDistance) || (minDistance == -1)) {
                    minDistance = currDistance;
                    minIndex = clusterIndex;
                }
            }
            
            int x = indexArray[minIndex] - 1;
            pointGroups[minIndex][x] = data[pointIndex];
            indexArray[minIndex] = indexArray[minIndex] - 1;
            pointIndex++;
        }
        
        initial_centroids = update_centroids(pointGroups, pointCounts, numClusters, numDimensions);
        free(pointCounts);
        pointCounts = (int*)calloc(numClusters, sizeof(int));
        if (pointCounts == NULL) {
            printf("An Error Has Occurred\n");
            return NULL;
        }
    }
    
    for (i = 0; i < numClusters; i++) {
        for (j = 0; j < pointCounts[i]; j++) {
            free(pointGroups[i][j]);
        }
        free(pointGroups[i]);
    }
    free(pointGroups);
    free(pointCounts);
    
    for (i = 0; i < numClusters; i++) {
        free(prevCentroids[i]);
    }
    
    free(prevCentroids);
    
    return initial_centroids;
}

