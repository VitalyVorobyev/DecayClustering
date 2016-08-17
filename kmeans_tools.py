# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 15:57:54 2016

@author: Vitaly Vorobyev
"""
import numpy  as np
import sys
import time
import os
import datetime
from scipy.sparse import csr_matrix
from sklearn.metrics import pairwise_distances
import matplotlib.pyplot as plt
import decay_tree_tools as dtt
import evtpdlparser as pdl

pdl_weights = {
 pdl.pname_id_dict['Upsilon(5S)'] : 50, # 9000553
 pdl.pname_id_dict['Upsilon(4S)'] : 50, # 300553
 pdl.pname_id_dict['Upsilon(3S)'] : 10, # 200553
 pdl.pname_id_dict['Upsilon(2S)'] :  7, # 100553
 pdl.pname_id_dict['Upsilon']     :  7, # 553
 pdl.pname_id_dict['B0']          : 10, # 511
 pdl.pname_id_dict['B+']          : 10, # 521
 pdl.pname_id_dict['B_s0']        :  7, # 531
 pdl.pname_id_dict['D+']          :  7, # 411
 pdl.pname_id_dict['D0']          :  7, # 421
 pdl.pname_id_dict['D_s+']        :  7, # 431
 pdl.pname_id_dict['J/psi']       :  7, # 443
 pdl.pname_id_dict['psi(2S)']     :  5, # 100443
 pdl.pname_id_dict['D*+']         :  4, # 413
 pdl.pname_id_dict['D*0']         :  4, # 423
 pdl.pname_id_dict['D_s*+']       :  4, # 433
 pdl.pname_id_dict['psi(3770)']   :  4, # 30443
 pdl.pname_id_dict['K_S0']        :  4, # 310
 pdl.pname_id_dict['K_L0']        :  4, # 130
 pdl.pname_id_dict['K+']          :  4, # 321
 pdl.pname_id_dict['pi+']         :  2, # 211
 pdl.pname_id_dict['pi0']         :  2, # 111
 pdl.pname_id_dict['p+']          :  2, # 2212
 pdl.pname_id_dict['n0']          :  2, # 2112
 pdl.pname_id_dict['e+']          :  2, # 11
 pdl.pname_id_dict['mu+']         :  2, # 13
 pdl.pname_id_dict['tau+']        :  2, # 15
}

def particle_weight(idhep):
    """ Retrieves particle weight for a given particle ID
    Args:
        idhep : particle EvtGen ID
    Returns:
        Integer weight.
        If EvtGen ID is not in the evt.pdl table returns 1
    """
    if idhep in pdl_weights:
        return pdl_weights[idhep]
    else:
        return 1

def particle_index(idhep):
    """ Retrieves index of particle in evt.pdl table for a given EvtGen ID
    Args:
        idhep : particle EvtGen ID
    Returns:
        Integer index.
        If EvtGen ID is not in the evt.pdl table returns 0
    """
    if idhep in pdl.pid_tuple:
        return pdl.pid_tuple.index(idhep)
    else:
        return 0

def decay_code(idhep,daF):
    """ Retrieves 1d NumPy arrays describing the decay chain
    Args:
        idhep: tuple of EvtGen ID
        daF: tuple of indeces corresponding to the first daughters
    Returns:
        1d NumPy arrays of integers of zise 539
    """
    code = np.zeros(pdl.numpdl)
    for i in range(0,len(idhep)):
        if daF[i] != -1:
            code[particle_index(idhep[i])] += particle_weight(abs(idhep[i]))
    return code

def decay_codes(idhep,daF,daL,pids):
    """
    """
    decay_trees = dtt.subdecays(idhep,daF,daL,pids)
    codes = []
    for tree in decay_trees:
        codes.append(decay_code(idhep,tree[0]))
    return tuple(codes)

def get_csr_matrix(idhepl,daFl,daLl,nevt=0):
    """
    
    """
    codes = []
    for i in xrange(nevt):
        if 911 in idhepl[i]:
            N = idhepl[i].index(911)
        else:
            N = len(idhepl[i])

        idhep = idhepl[i][:N]
        daF   = daFl[i][:N]
        daL   = daLl[i][:N]
        for code in decay_codes(idhep,daF,daL,dtt.b_mesons):
            codes.append(code)
        if (i % 1000) == 0:
            print("%d events..." % i)
            sys.stdout.flush()
        if nevt != 0:
            if i == nevt:
                break
    return csr_matrix(np.array(codes))

def tree_to_csr_matrix(tree,nevt=0):
    """
    
    """
    codes = []
    i = 0
    for evt in tree:
        i += 1
        idhep = tuple(tree.gen_idhep)
        if 911 in idhep:
            N = idhep.index(911)
        else:
            N = len(idhep)

        idhep = idhep[:N]
        daF   = tuple(tree.gen_daF)[:N]
        daL   = tuple(tree.gen_daL)[:N]
        for code in decay_codes(idhep,daF,daL,dtt.b_mesons):
            codes.append(code)
        if (i % 1000) == 0:
            print("%d events..." % i)
            sys.stdout.flush()
        if nevt != 0:
            if i == nevt:
                break
    return csr_matrix(np.array(codes))
    
def save_sparse_csr(filename,array):
    """ Saves sparse matrix to .npz file
    """
    np.savez(filename,data = array.data ,indices=array.indices,
             indptr =array.indptr, shape=array.shape )

def load_sparse_csr(filename):
    """ Loads sparse matrix from .npz file
    Args:
        filename: file name
    Returns:
        csr_matrix
    """
    loader = np.load(filename)
    return csr_matrix((  loader['data'], loader['indices'], loader['indptr']),
                         shape = loader['shape'])

def b_decay_codes(tree,evt):
    """ Retrieves 
    """
    idhep, daF,daL = dtt.gen_table(tree,evt)
    return decay_codes(idhep,daF,daL,dtt.b_mesons)

def get_random_centroids(data, k, seed=None):
    """ Randomly choose k data points as initial centroids
    """
    if seed is not None: # useful for obtaining consistent results
        np.random.seed(seed)
    n = data.shape[0] # number of data points    
    # Pick K indices from range [0, N).
    rand_indices = np.random.randint(0, n, k)
    # Keep centroids as dense format, as many entries will be nonzero due to averaging.
    # As long as at least one document in a cluster contains a word,
    # it will carry a nonzero weight in the TF-IDF vector of the centroid.
    centroids = data[rand_indices,:].toarray()
    return centroids

def get_kpp_centroids(data, k, seed=None):
    """ Use k-means++ to initialize a good set of centroids
    """
    if seed is not None: # useful for obtaining consistent results
        np.random.seed(seed)
    centroids = np.zeros((k, data.shape[1]))
    
    # Randomly choose the first centroid.
    # Since we have no prior knowledge, choose uniformly at random
    idx = np.random.randint(data.shape[0])
    centroids[0] = data[idx,:].toarray()
    # Compute distances from the first centroid chosen to all the other data points
    distances = pairwise_distances(data, centroids[0:1], metric='euclidean').flatten()

    for i in xrange(1, k):
        # Choose the next centroid randomly, so that the probability for each data point to be chosen
        # is directly proportional to its squared distance from the nearest centroid.
        # Roughtly speaking, a new centroid should be as far as from ohter centroids as possible.
        idx = np.random.choice(data.shape[0], 1, p=distances/sum(distances))
        centroids[i] = data[idx,:].toarray()
        # Now compute distances from the centroids to all data points
        distances = np.min(pairwise_distances(data, centroids[0:i+1], metric='euclidean'),axis=1)
    
    return centroids

def assign_clusters(data, centroids):  
    """
    
    """
    # Compute distances between each data point and the set of centroids:
    # Fill in the blank (RHS only)
    distances_from_centroids = pairwise_distances(data, centroids, metric='euclidean')
    
    # Compute cluster assignments for each data point:
    # Fill in the blank (RHS only)
    cluster_assignment = cluster_assignment = np.argmin(distances_from_centroids,axis=1)
    
    return cluster_assignment
    
def revise_centroids(data, k, cluster_assignment):
    """
    
    """
    new_centroids = []
    for i in xrange(k):
        # Select all data points that belong to cluster i. Fill in the blank (RHS only)
        member_data_points = data[cluster_assignment == i,:]
        # Compute the mean of the data points. Fill in the blank (RHS only)
        centroid = member_data_points.mean(axis=0)
        
        # Convert numpy.matrix type to numpy.ndarray type
        centroid = centroid.A1
        new_centroids.append(centroid)
    new_centroids = np.array(new_centroids)
    
    return new_centroids

def cluster_heterogeneity(data,centroid):
    distances = pairwise_distances(data, [centroid], metric='euclidean')
    squared_distances = distances**2
    return(np.sum(squared_distances))
    
def ith_cluster_heterogeneity(data,centroids,cluster_assignment,i):
    # Select all data points that belong to cluster i. Fill in the blank (RHS only)
    d = data[cluster_assignment==i, :]
    if d.shape[0] == 0:
        return -1
    else:
        return cluster_heterogeneity(d,centroids[i])
    
def compute_heterogeneity(data, k, centroids, cluster_assignment, for_each_cluster=None):
    """ Computes sum of squared distances
    Args:
        data:
        k:
        centroids:
        cluster_assignment:
        for_each_cluster:
    Returns:
        Sum of squared distances to the nearest centroid summed up over
        all points. If for_each_cluster is specified, heterogeneity list 
        for each cluster is stored in for_each_cluster
    """
    if for_each_cluster is not None:
        for_each_cluster = []
    heterogeneity = 0.0
    for i in xrange(k):
        Hi = ith_cluster_heterogeneity(data,centroids, cluster_assignment,i)
        heterogeneity += Hi
        if for_each_cluster is not None:
                for_each_cluster.append(Hi)
    return heterogeneity

def kmeans(data, k, initial_centroids, maxiter, record_heterogeneity=None, verbose=False):
    """ This function runs k-means on given data and initial set of centroids.
       maxiter: maximum number of iterations to run.
       record_heterogeneity: (optional) a list, to store the history of heterogeneity as function of iterations
                             if None, do not store the history.
       verbose: if True, print how many data points changed their cluster labels in each iteration
    """
    centroids = initial_centroids[:]
    prev_cluster_assignment = None
    
    for itr in xrange(maxiter):        
        if verbose:
            print(itr)
        
        # 1. Make cluster assignments using nearest centroids
        cluster_assignment = assign_clusters(data,centroids)
            
        # 2. Compute a new centroid for each of the k clusters, averaging all data points assigned to that cluster.
        centroids = revise_centroids(data, k, cluster_assignment)
            
        # Check for convergence: if none of the assignments changed, stop
        if prev_cluster_assignment is not None and \
          (prev_cluster_assignment==cluster_assignment).all():
            break
        
        # Print number of new assignments 
        if prev_cluster_assignment is not None:
            num_changed = np.sum(prev_cluster_assignment!=cluster_assignment)
            if verbose:
                print('    {0:5d} elements changed their cluster assignment.'.format(num_changed))
                sys.stdout.flush()
        
        # Record heterogeneity convergence metric
        if record_heterogeneity is not None:
            score = compute_heterogeneity(data, k, centroids, cluster_assignment)
            record_heterogeneity.append(score)
        
        prev_cluster_assignment = cluster_assignment[:]
        
    return centroids, cluster_assignment
    
def plot_heterogeneity(heterogeneity, k):
    """
    
    """
    plt.figure(figsize=(7,4))
    plt.plot(heterogeneity, linewidth=4)
    plt.xlabel('# Iterations')
    plt.ylabel('Heterogeneity')
    plt.title('Heterogeneity of clustering over time, K={0:d}'.format(k))
    plt.rcParams.update({'font.size': 16})
    plt.tight_layout()

def get_datetime_dirname():
    return os.path.join(os.getcwd(),'plots/', datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    
def save_gvs(gvs,nclusters_to_show = 0):
    if len(gvs) == 0:
        return

    path = get_datetime_dirname();
    print(path)
    
    for c in range(len(gvs)):
        save_gv_list(path,gvs[c],c, c<nclusters_to_show)
            
def save_gv_list(path,gvt,cluster,show_flag = False):
    if len(gvt) == 0:
        return

    for nn in range(len(gvt)):
        gv = gvt[nn]
        gv.render(path + '/b_decay_tree_Cluster' + str(cluster) + '_nn' + str(nn), view=( (nn == 0) & (show_flag) ))

def visualize_clusters(tree, matrix, centroids, cluster_assignment, k, nclusters_to_show = 8, ntrees_to_save=0, collect_gvt = None):
    """
    
    """
    bcount = np.bincount(cluster_assignment)
    selected_cluster_idxs = tuple(bcount.argsort()[::-1])
    path = get_datetime_dirname();
    print('==========================================================')
    # Visualize each cluster c
    for i in xrange(k):
        cidx = selected_cluster_idxs[i]
        # Cluster heading
        print('Cluster {0:d}    '.format(i)),
        # Print top 7 particle types with largest weights in the cluster
        idx = centroids[cidx].argsort()[::-1]
        for j in xrange(7): # Print each word along with the TF-IDF weight
            print('{0:s}:{1:.3f}'.format(pdl.pname_tuple[idx[j]], centroids[cidx,idx[j]])),
        print('')

        gv_trees = []
#        if ntrees_to_save>0:
            # Compute distances from the centroid to all data points in the cluster,
            # and compute nearest neighbors of the centroids within the cluster.
        distances = pairwise_distances(matrix, centroids[cidx].reshape(1, -1), metric='euclidean').flatten()
        distances[cluster_assignment!=cidx] = float('inf') # remove non-members from consideration
        nearest_neighbors = distances.argsort()
#        print("%d minimal distances to centroid:" % ntrees_to_save)
        for nn in xrange(ntrees_to_save):
            idx = nearest_neighbors[nn]
#            print(distances[idx])
            evt= int(idx / 2)
            gv_trees.append(dtt.b_meson_gv_decay_trees(tree,evt)[idx % 2])
        if(collect_gvt is not None):
            collect_gvt[i] = tuple(gv_trees)
        save_gv_list(path,gv_trees,i, i<nclusters_to_show)
        print("Cluster size: %d" % bcount[cidx])
        Hi = ith_cluster_heterogeneity(matrix,centroids,cluster_assignment,cidx)
        Hi = np.sqrt(Hi)/bcount[cidx]
        print("Cluster heterogeneity: %lf" % (Hi))
        print('==========================================================')
#    save_gvs(collect_gvt,nclusters_to_show)

#def boxplot_heterogenious(heterogeneity,heterogeneity_smart):
#    """
#    
#    """
#    plt.figure(figsize=(8,5))
#    plt.boxplot([heterogeneity.values(), heterogeneity_smart.values()], vert=False)
#    plt.yticks([1, 2], ['k-means', 'k-means++'])
#    plt.rcParams.update({'font.size': 16})
#    plt.tight_layout()

#def kmeans_multiple_runs(data, k, maxiter, num_runs, seed_list=None, verbose=False):
#    """
#    Args:
#        data: sparse matrix N*M, where N is number of events and M = 539 is evt.pdl table size
#        k: number of clusters
#        maxiter: max number of iterations
#        seed_list: self-explanatory
#        verbose: self-explanatory
#    Returns:
#        Centroids and cluster assignments that minimize heterogeneity
#    """
#    heterogeneity = {}
#
#    min_heterogeneity_achieved = float('inf')
#    final_centroids = None
#    final_cluster_assignment = None
#    
#    for i in xrange(num_runs):
#        
#        # Use UTC time if no seeds are provided 
#        if seed_list is not None: 
#            seed = seed_list[i]
#            np.random.seed(seed)
#        else: 
#            seed = int(time.time())
#            np.random.seed(seed)
#        
#        # Use k-means++ initialization
#        initial_centroids = get_kpp_centroids(data, k, seed)
#        
#        # Run k-means
#        centroids, cluster_assignment = kmeans(data, k, initial_centroids, maxiter=400,
#                                           record_heterogeneity=None, verbose=False)
#        
#        # To save time, compute heterogeneity only once in the end
#        heterogeneity[seed] = compute_heterogeneity(data, k, centroids, cluster_assignment)
#        
#        if verbose:
#            print('seed={0:06d}, heterogeneity={1:.5f}'.format(seed, heterogeneity[seed]))
#            sys.stdout.flush()
#        
#        # if current measurement of heterogeneity is lower than previously seen,
#        # update the minimum record of heterogeneity.
#        if heterogeneity[seed] < min_heterogeneity_achieved:
#            min_heterogeneity_achieved = heterogeneity[seed]
#            final_centroids = centroids
#            final_cluster_assignment = cluster_assignment
#    
#    # Return the centroids and cluster assignments that minimize heterogeneity.
#    return final_centroids, final_cluster_assignment

#def plot_k_vs_heterogeneity(k_values, heterogeneity_values):
#    plt.figure(figsize=(7,4))
#    plt.plot(k_values, heterogeneity_values, linewidth=4)
#    plt.xlabel('K')
#    plt.ylabel('Heterogeneity')
#    plt.title('K vs. Heterogeneity')
#    plt.rcParams.update({'font.size': 16})
#    plt.tight_layout()

#def bipartition(cluster, maxiter=400, num_runs=4, seed=None):
#    '''cluster: should be a dictionary containing the following keys
#                * dataframe: original dataframe
#                * matrix:    same data, in matrix format
#                * centroid:  centroid for this particular cluster'''
#    
#    data_matrix = cluster['matrix']
#    dataframe   = cluster['dataframe']
#    
#    # Run k-means on the data matrix with k=2. We use scikit-learn here to simplify workflow.
#    kmeans_model = KMeans(n_clusters=2, max_iter=maxiter, n_init=num_runs, random_state=seed, n_jobs=-1)    
#    kmeans_model.fit(data_matrix)
#    centroids, cluster_assignment = kmeans_model.cluster_centers_, kmeans_model.labels_
#    
#    # Divide the data matrix into two parts using the cluster assignments.
#    data_matrix_left_child, data_matrix_right_child = data_matrix[cluster_assignment==0], \
#                                                      data_matrix[cluster_assignment==1]
#    
#    # Divide the dataframe into two parts, again using the cluster assignments.
#    cluster_assignment_sa = graphlab.SArray(cluster_assignment) # minor format conversion
#    dataframe_left_child, dataframe_right_child     = dataframe[cluster_assignment_sa==0], \
#                                                      dataframe[cluster_assignment_sa==1]
#        
#    
#    # Package relevant variables for the child clusters
#    cluster_left_child  = {'matrix': data_matrix_left_child,
#                           'dataframe': dataframe_left_child,
#                           'centroid': centroids[0]}
#    cluster_right_child = {'matrix': data_matrix_right_child,
#                           'dataframe': dataframe_right_child,
#                           'centroid': centroids[1]}
#    
#    return (cluster_left_child, cluster_right_child)
