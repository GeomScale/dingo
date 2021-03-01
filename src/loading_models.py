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


# The .mat format case
def read_mat_file(input_file):

   data_from_mat = scipy.io.loadmat(input_file)

   species_name = ''
   for key in data_from_mat.keys():
      if key[0] != "_":
         species_name = key

   species = data_from_mat[species_name]
   list_of_lists = species.tolist()

   counter = 0

   metabolites = []

   for element in list_of_lists[0][0]:

      if counter == 0:

         m =len(element)

         for i in element:

            metabolite = i[0][0]

            if metabolite not in metabolites:
               metabolites.append(metabolite)

      if counter == 7:
         reactions_list = element.tolist()
         reactions = [reaction[0][0] for reaction in reactions_list]
      if counter == 11:
         lb_tmp = element
         lb_tmp = lb_tmp.tolist()
      if counter == 12:
         ub_tmp = element
      if counter == 10:
         Aeq = element

      counter += 1

   # Build function's output; first the A matrix
   n = len(ub_tmp)
   A = np.zeros((2*n, n), dtype=np.float)
   A[0:n] = np.eye(n)
   A[n:] -=  np.eye(n,n, dtype=np.float)

   # Now, the b vector
   ub = [i[0] for i in ub_tmp]
   lb = [-x[0] for x in lb_tmp]
   b = np.asarray(ub + lb)

   # The Aeq matrix
   Aeq = np.asarray(Aeq)

   # And the beq vector
   beq = np.zeros(m)

   return A, b, Aeq, beq, metabolites, reactions
