import numpy as np 
import eam_calculator as eamc
import matplotlib.pyplot as plt
def calculate_distance(rij):
  """
  Args:
    rij (np.array): vector distance between particle i and j, shape(1, 3)
  Return:
    float: abs distance between i and j
  """
  return np.sum(rij ** 2) ** 0.5

def my_pair_correlation(dists, natom, nbins, dr, lbox):
  """ Calculate the pair correlation function g(r).

  the first two steps are provided:
  histogram is a tuple (counts, bin_edges)
  and r gives the centers of the bins
  please see documentation on numpy.histogram() for more details
  
  your task is to transform histogram[0] into g(r)

  Args:
    dists (np.array): 1d array of pair distances
    natom (int): number of atoms
    nbins (int): number of bins to histogram
    dr (float): size of bins
    lbox (float): side length of cubic box
  Return:
    array of shape (nbins,): the pair correlation g(r)
  """
  histogram = np.histogram(dists, bins=nbins, range=(0, nbins*dr))
  r = (histogram[1] + dr/2)[:-1] # centers of the bins
  hi = []
  n_pairs = len(dists)
  for i in r:
      hi.append((4 / 3 * np.pi) * ((i + dr / 2) ** 3 - (i - dr / 2) ** 3) * (  n_pairs / lbox ** 3))
  return histogram[0] / np.array(hi), r
  pass

# get distance table
def get_distance_table(position, particle_number, box_size):
  distance_table = np.zeros((particle_number, particle_number))
  for i in range(particle_number):
   for j in range(particle_number):
     distance_table[i][j] = calculate_distance(eamc.minimum_image(position[i] - position[j], box_size))
  return distance_table

def exclude_replicated_pair(distance_table):
  pair_dis = []
  for i in range(distance_table.shape[0]):
    for j in range(i + 1, distance_table.shape[1]):
      pair_dis.append(distance_table[i][j])
  return pair_dis

def get_pair_correlation(T, position_arr, N, L):
    n_bins = 20
    gr_total = np.array([])
    bin_centers = np.array([])
    count = 0
    for i in position_arr:
        count+=1
        dis_table = get_distance_table(i, N, L)
        pair_distance = exclude_replicated_pair(dis_table)
        gr, bin_centers = my_pair_correlation(pair_distance, N, n_bins, L/(2 * n_bins), L)
        if len(gr_total) == 0:
            gr_total = gr.copy()
        else:
            gr_total += gr.copy()
    
    gr_avg = gr_total / len(position_arr)

    return gr_avg