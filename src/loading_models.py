import json
import scipy.io
import numpy as np



## Read a Bigg file and get the necessary A and b
# The .json format case
def read_json_file(input_file):

   with open(input_file, 'r') as f:

      distros_dict = json.load(f)

      reactions_list = distros_dict['reactions']

      metabolites = []
      reactions = []

      for reaction in reactions_list:

         metabolites_dic = reaction['metabolites']
         reaction_name = reaction['id']
         reactions.append(reaction_name)

         for metabolite in metabolites_dic.keys():
            if metabolite not in metabolites:
               metabolites.append(metabolite)

      list_of_reaction_lists = []
      vector_of_ubs = []
      vector_of_lbs = []

      for reaction in reactions_list:

         ub = float(reaction['upper_bound']) ; vector_of_ubs.append(ub)
         lb = float(reaction['lower_bound']) ; vector_of_lbs.append(lb)

         metabolites_dic = reaction['metabolites']
         reaction_vector = []

         for term in range(len(metabolites)):

            metabolite = metabolites[term]

            if metabolite in metabolites_dic.keys():

               reaction_vector.append(metabolites_dic[metabolite])
            else:
               reaction_vector.append(0)

         list_of_reaction_lists.append(reaction_vector)

   f.close()

   # Build function's output; first the A matrix
   n = len(list_of_reaction_lists)
   A = np.zeros((2*n, n), dtype=np.float)
   A[0:n] = np.eye(n)
   A[n:] -=  np.eye(n,n, dtype=np.float)

   # Now, the b vector
   vector_of_lbs = [-x for x in vector_of_lbs]
   b = np.asarray(vector_of_ubs + vector_of_lbs)

   # The Aeq matrix
   Aeq = np.asarray(list_of_reaction_lists)
   Aeq = np.transpose(Aeq)

   # And the beq vector
   m = len(metabolites)
   beq = np.zeros(m)
   
   # Make everything C contigeous
   A = np.asarray(A, dtype = 'float')
   A = np.ascontiguousarray(A, dtype='float')
   b = np.asarray(b, dtype = 'float')
   b = np.ascontiguousarray(b, dtype='float')
   Aeq = np.asarray(Aeq, dtype = 'float')
   Aeq = np.ascontiguousarray(Aeq, dtype='float')
   beq = np.asarray(beq, dtype = 'float')
   beq = np.ascontiguousarray(beq, dtype='float')

   return A, b, Aeq, beq, metabolites, reactions

