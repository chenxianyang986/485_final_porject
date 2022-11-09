import numpy as np
def my_calc_rhok(kvecs, pos):
  """ Calculate the fourier transform of particle density.

  Args:
    kvecs (np.array): array of k-vectors, shape (nk, ndim)
    pos (np.array): particle positions, shape (natom, ndim)
  Return:
    array of shape (nk,): fourier transformed density rho_k
  """
  rho_k = []
  for i in range(kvecs.shape[0]):
    sum = 0
    for rj in pos:
      sum += complex(np.cos(np.dot(np.array(rj), np.array(kvecs[i]))), - np.sin(np.dot(np.array(rj), np.array(kvecs[i]))))
    rho_k.append(sum)
  return rho_k
  

def my_calc_sk(kvecs, pos):
  """ Calculate the structure factor S(k).

  Args:
    kvecs (np.array): array of k-vectors, shape (nk, ndim)
    pos (np.array): particle positions, shape (natom, ndim)
  Return:
    array of shape (nk,): structure factor s(k)
  """
  rho_k = my_calc_rhok(kvecs, pos)
  sk = []
  for i in rho_k:
    sk.append(1 / pos.shape[0] * (i.real ** 2 + i.imag ** 2))
  return sk
  pass

def my_legal_kvecs(maxn, lbox):
  """ Calculate k vectors commensurate with a cubic box.

  Consider only k vectors in the all-positive octant of reciprocal space.

  Args:
    maxn : maximum value for nx, ny, nz; maxn+1 is number of k-points along each axis
    lbox : side length of cubic cell

  Return:
    array of shape (nk, ndim): collection of k vectors
  """
  kvecs = []
  kdim = []
  for i in range(maxn + 1):
    kdim.append(2 * np.pi / lbox * i)
  for i in range(maxn + 1):
      for j in range(maxn + 1):
          for k in range(maxn + 1):
              kvecs.append([kdim[i], kdim[j], kdim[k]])
  return np.array(kvecs)

def structure_factor_plot_helper(kvecs, sk_list):
  ''' Select unique k_vectors and plot corresonding sk_list 
  for plotting 
  Args:
    kvecs: collections of kvecs (nk, ndim)
    sk_list: array of shape nk 
  '''
  kmags = [np.linalg.norm(kvec) for kvec in kvecs]
  sk_arr = np.array(sk_list)

  unique_kmags = np.unique(kmags)
  unique_sk = np.zeros(len(unique_kmags))

  for iukmag in range(len(unique_kmags)):
    kmag = unique_kmags[iukmag]
    idx2avg = np.where(kmags == kmag)
    unique_sk[iukmag] = np.mean(sk_arr[idx2avg])
  return unique_kmags, unique_sk