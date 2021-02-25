import numpy as np
import scipy.sparse as sp

# A function to zoom in on the polytope
def gmscale(A, iprint, scltol):

   #--------------------------------------------------------------------------------------
   #
   # This function is a translation of a matlab cobra script you may find here:
   # https://github.com/opencobra/cobratoolbox/blob/master/src/analysis/subspaces/gmscale.m 
   #
   # USAGE:
   # cscale, rscale = gmscale(A, iprint, scltol)

   # INPUTS:
   #     A:          is a (m*n) sparse matrix 
   #     scltol:      should be in the range (0.0, 1.0).
   #              Typically `scltol` = 0.9.  A bigger value like 0.99 asks
   #              gmscale` to work a little harder (more passes).   

   # OUTPUTS:
   #    cscale, rscale:    column vectors of column and row scales such that
   #                       `R` (inverse) `A` `C` (inverse) should have entries near 1.0,
   #                        where `R= diag(rscale)`, `C = diag(cscale)`.
   #
   #--------------------------------------------------------------------------------------
   
   m = A.shape[0] ; n = A.shape[1]
   A = np.abs(A)
   A_t = A.T     ## We will work with the transpose matrix as numpy works on row based
   maxpass = 10
   aratio = 1e+50
   damp = 1e-4
   small = 1e-8
   rscale = np.ones((m,1))
   cscale = np.ones((n,1))

   #----------------------------------------------------
   # Main loop
   #----------------------------------------------------

   for npass in range(maxpass):
      
      rscale[rscale == 0] = 1
      sparse = sp.csr_matrix(1./rscale)
      r_diagonal_elements = np.asarray(sparse.todense())      
      Rinv = diags(r_diagonal_elements[0].T, shape = (m,m)).toarray()
      
      SA = np.dot(Rinv, A)
      SA_T = SA.T
      
      J, I = SA_T.nonzero()
      V = SA.T[SA.T > 0]
      invSA = sp.csr_matrix( (1./V, (J, I)), shape = (n, m) ).T
      
      cmax = np.max(SA, axis = 0)
      cmin = np.max(invSA, axis = 0).data
      cmin = 1./ (cmin + 2.2204e-16)
      
      sratio = np.max(cmax/cmin)
      
      if npass > 0:
         c_product = np.multiply(np.max((np.array((cmin.data, damp*cmax))), axis =0), cmax)
         cscale = np.sqrt(c_product)

      check = aratio*scltol     
      
      if npass >= 2 and sratio >= check:
         break

      if npass == maxpass:
         break
      
      aratio = sratio
      
      # Set new row scales for the next pass.
      
      cscale[cscale == 0] = 1
      sparse = sp.csr_matrix(1./cscale)
      c_diagonal_elements = np.asarray(sparse.todense())
      
      Cinv = diags(c_diagonal_elements[0].T, shape = (n,n)).toarray()      
      SA = np.dot(A, Cinv)
      SA_T = SA.T
      
      J, I = SA_T.nonzero()
      V = SA.T[SA.T > 0]
      invSA = sp.csr_matrix( (1./V, (J, I)), shape = (n, m) ).T
      
      rmax = np.max(SA, axis = 1)
      rmin = np.max(invSA, axis = 1).data
      tmp = rmin + 2.2204e-16
      rmin = 1.0/tmp
      
      r_product = np.multiply(np.max((np.array((rmin.data, damp*rmax))), axis =0), rmax)
      rscale = np.sqrt(r_product)      
      
   #----------------------------------------------------
   # End of main loop
   #----------------------------------------------------

   # Reset column scales so the biggest element in each scaled column will be 1.
   # Again, allow for empty rows and columns.
   rscale[rscale == 0] = 1
   sparse = sp.csr_matrix(1./rscale)
   r_diagonal_elements = np.asarray(sparse.todense())      
   Rinv = diags(r_diagonal_elements[0].T, shape = (m,m)).toarray()

   SA = np.dot(Rinv, A)
   SA_T = SA.T  
   J, I = SA_T.nonzero()
   V = SA.T[SA.T > 0]
   
   cscale = np.max(SA, axis = 0)
   cscale[cscale == 0] = 1
   
   if iprint > 0 :
   
      rmin = np.amin(rscale) ; imin = np.where(rscale == np.amin(rscale))
      cmin = np.amin(cscale) ; jmin = np.where(cscale == np.amin(cscale))
      
      rmax = np.amax(rscale) ; imax = np.where(rscale == np.amax(rscale))
      cmax = np.amax(cscale) ; jmax = np.where(cscale == np.amax(cscale))
      
   return cscale, rscale
     

     