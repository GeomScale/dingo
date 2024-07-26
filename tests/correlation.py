
from dingo.illustrations import plot_corr_matrix
from dingo.utils import correlated_reactions
from dingo import MetabolicNetwork, PolytopeSampler
import numpy as np


#model = MetabolicNetwork.from_json('ext_data/e_coli_core.json')
#reactions = model.reactions

#sampler = PolytopeSampler(model)
#steady_states = sampler.generate_steady_states()

#plot_corr_matrix(steady_states, reactions, color="RdYlBu")
#correlated_reactions(steady_states, reactions, pearson_cutoff = 0.80, n = 7)


# Creation of 4 possible copulas classifying 2 reactions as
# Positive, Negative or No Correlated


# Positive Correlation
data = [
    [2, 2, 0, 1, 1, 1, 1, 1, 1],
    [2, 2, 2, 0, 1, 1, 1, 1, 1],
    [0, 2, 2, 2, 0, 1, 1, 1, 1],
    [1, 0, 2, 2, 2, 0, 1, 1, 1],
    [1, 1, 0, 2, 2, 2, 0, 1, 1],
    [1, 1, 1, 0, 2, 2, 2, 0, 1],
    [1, 1, 1, 1, 0, 2, 2, 2, 0],
    [1, 1, 1, 1, 1, 0, 2, 2, 2],
    [1, 1, 1, 1, 1, 1, 0, 2, 2]
     
]


# Negative Correlation
data = [
    [1, 1, 1, 1, 1, 1, 0, 2, 2],
    [1, 1, 1, 1, 1, 0, 2, 2, 2],
    [1, 1, 1, 1, 0, 2, 2, 2, 0],
    [1, 1, 1, 0, 2, 2, 2, 0, 1],
    [1, 1, 0, 2, 2, 2, 0, 1, 1],
    [1, 0, 2, 2, 2, 0, 1, 1, 1],
    [0, 2, 2, 2, 0, 1, 1, 1, 1],
    [2, 2, 2, 0, 1, 1, 1, 1, 1],
    [2, 2, 0, 1, 1, 1, 1, 1, 1]
     
]


# No Correlation
data = [
    [2, 2, 1, 1, 1, 1, 1, 1, 1],
    [2, 2, 2, 1, 1, 1, 1, 1, 1],
    [1, 2, 2, 2, 1, 1, 1, 1, 1],
    [1, 1, 2, 2, 2, 1, 1, 1, 1],
    [1, 1, 1, 2, 2, 2, 1, 1, 1],
    [1, 1, 1, 1, 2, 2, 2, 1, 1],
    [1, 1, 1, 1, 1, 2, 2, 2, 1],
    [1, 1, 1, 1, 1, 1, 2, 2, 2],
    [1, 1, 1, 1, 1, 1, 1, 2, 2]
     
]


# No Correlation
data = [
    [2, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 2, 0, 0, 0, 1, 1, 1, 1],
    [0, 0, 2, 0, 0, 0, 1, 1, 1],
    [0, 0, 0, 2, 0, 0, 0, 1, 1],
    [1, 0, 0, 0, 2, 0, 0, 0, 1],
    [1, 1, 0, 0, 0, 2, 0, 0, 0],
    [1, 1, 1, 0, 0, 0, 2, 0, 0],
    [1, 1, 1, 1, 0, 0, 0, 2, 0],
    [1, 1, 1, 1, 1, 0, 0, 0, 2]
     
]


# Convert nested list to numps 2D array
copula1 = np.array(data)
# Flip array upside down to get the other diagonal starting in the top left corner
copula2 = np.flipud(copula1)

print("Copula 1:\n" , copula1 , end = "\n\n")
print("Copula 2:\n" , copula2 , end = "\n\n")

# Get the dimensions of the copula
rows, cols = copula1.shape

# Variable to store the sum of values across the 1st diagonal
cop_a_diag_sum = 0
# Variable to store the count of values across the 1st diagonal
cop_a_diag_count = 0
# Variable to store the sum of values across the edges (non 1st diagonal values)
cop_a_edge_sum = 0
# Variable to store the count of values across the edges (non 1st diagonal counts)
cop_a_edge_count = 0

# Variable to store the sum of values across the 2st diagonal
cop_b_diag_sum = 0
# Variable to store the count of values across the 2st diagonal
cop_b_diag_count = 0
# Variable to store the sum of values across the edges (non 2nd diagonal values)
cop_b_edge_sum = 0
# Variable to store the count of values across the edges (non 2nd diagonal counts)
cop_b_edge_count = 0

# Iterate every copula's row
for row in range(rows):
    # Iterate every copula's column
    for col in range(cols):
        # Find combination of rows and column near the diagonal
        if ((row-col >= -0.2*rows) & (row-col <= 0.2*rows)):
            # Sum the values of combinations that met the criteria
            cop_a_diag_sum = cop_a_diag_sum + copula1[row][col]
            # Count the occurences of combinations that met the criteria
            cop_a_diag_count += 1
            
            # Do the same for the flipped copula (2nd diagonal)
            cop_b_diag_sum = cop_b_diag_sum + copula2[row][col]
            cop_b_diag_count += 1
            
        else:
            # Do the same for values not belonging in the 1st diagonal
            cop_a_edge_sum = cop_a_edge_sum + copula1[row][col]
            cop_a_edge_count += 1
            
            # Do the same for values not belonging in the 2nd diagonal
            cop_b_edge_sum = cop_b_edge_sum + copula2[row][col]
            cop_b_edge_count += 1
            
            
# Calculate the fraction of mass in the diagonal to the mass in the edges  (1st diagonal)
# Add a value of 1e-9 to avoid inf values
diagonal_a_avg = ((cop_a_diag_sum / cop_a_diag_count + 1e-9) / 
                  (cop_a_edge_sum / cop_a_edge_count + 1e-9))

# Calculate the fraction of mass in the diagonal to the mass in the edges  (2nd diagonal)
# Add a value of 1e-9 to avoid inf values
diagonal_b_avg = ((cop_b_diag_sum / cop_b_diag_count + 1e-9) / 
                  (cop_b_edge_sum / cop_b_edge_count + 1e-9))

# If the fraction between the 2 diagonals is similar ==> No Correlation
if( ((diagonal_a_avg / diagonal_b_avg) < 2) and ((diagonal_b_avg / diagonal_a_avg) < 2) ):
    print("No Correlation")
else:
    # If 1st diagonal has a higher fraction
    if diagonal_a_avg > diagonal_b_avg:
        print("Positive Correlation")
    # If 2nd diagonal has a higher fraction        
    else:
        print("Negative Correlation")