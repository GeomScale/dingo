import pickle 
import numpy as np
import sys
import hopsy

def binary_search_ess(left_index: int, right_index: int, unified_chain):

    previous_left_index  = -1
    previous_right_index = -1

    while True:

        print("left_index: ", left_index)
        print("right_index: ", right_index)
        print("\n\n~~~~~~\n\n")

        if  right_index - left_index == 1:
            index = right_index
        else:
            index = left_index + int((right_index - left_index) / 2)


        index = index - index % 5


        chains_on_index = np.array(np.split(unified_chain[:index], 5))

        chains_on_left_index = np.array(np.split(unified_chain[:left_index], 5))

        ess_value_left = hopsy.ess(chains_on_left_index)
        ess_value = hopsy.ess(chains_on_index)
        
        if ess_value_left.min() > 1000:
            return index
            break
        elif ess_value.min() > 1000:
            right_index = index
        else:
            left_index = index

        if left_index == previous_left_index and right_index == previous_right_index:
            oops = np.array(np.split(unified_chain[:right_index], 5))
            print(hopsy.ess(oops).min())
            return oops
        else:
            previous_left_index = left_index
            previous_right_index = right_index 
        
with open(sys.argv[1], "rb") as f:
    samples = pickle.load(f)    

unified_chain = np.reshape(samples,  (samples.shape[0] * samples.shape[1], samples.shape[2]))    

num_of_steps = int(sys.argv[2])

left  = (samples.shape[0]-1) * num_of_steps
right = samples.shape[0] * num_of_steps 

binary_search_ess(left, right, unified_chain)

