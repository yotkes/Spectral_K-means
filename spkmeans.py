import spkmeansmodule as spm
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()

# Receiving inputs from user
parser.add_argument("k", type=int, help="Number of centroids - integer >= 1")
parser.add_argument("goal", type=str, help="Indicates specific procedure")
parser.add_argument("file_name", type=str, help="Name of file")
args = parser.parse_args()
k = args.k
goal = args.goal
filename = args.file_name

df = pd.read_csv(filename, header=None) # Converting file into a data frame
N = len(df) # File's number of rows (# of observations)
d = len(df.columns) # File's number of columns (vector's dimension)

# Converting goal's value into an integer
if(goal == "spk"):
    goal = 0
elif (goal == "wam"):
    goal = 1
elif (goal == "ddg"):
    goal = 2
elif (goal == "lnorm"):
    goal = 3
elif (goal == "jacobi"):
    goal = 4
else:
    print("Invalid Input!")
    exit(0)

# Checks for invalid input
if ((k>=N) or (d<=0) or (k<0) or (N<=0)):
    print("Invalid Input!")
    exit(0)

# Method that calculates the difference between an observation (a vector) and a centroid
def calculate_difference(vector, centroid):
    sum = 0
    for j in range(k):
        sum += ((vector[j] - centroid[j])**2)
    return sum

# Method that checks the closest centroid to a specific observation (vector)
def check_min_distance(vector, centroids):
    dist = calculate_difference(vector,centroids[0])
    for cluster in range(1, len(centroids)):
        new_dist = calculate_difference(vector, centroids[cluster])
        if (new_dist<dist):
            dist = new_dist
    return dist

# First Phase of K-means ++ Algorithm - computing k initial centroids
def k_means_pp(obs):
    np.random.seed(0)
    rand = np.random.choice(N, 1)
    indexes = [rand[0]]
    centroids = [obs[rand[0]]]
    for j in range(1,k):
        dists = [check_min_distance(obs[i], centroids) for i in range(N)]
        sums = sum(dists)
        probs = [dists[i]/sums for i in range(N)]
        rand = np.random.choice(N,1,p=probs)
        centroids.append(obs[rand[0]])
        indexes.append(rand[0])

    rep = ""
    for i in indexes:
        rep += str(i)
        rep += ","
    print(rep[:-1])
    return indexes

obs = df.to_numpy() # Converting the data frame into a NumPy array
to_c_obs = obs.tolist() # Returns the array as a deep nested list of Python scalars.
to_send = [k, N, d, goal, to_c_obs] # List of inputs sent to C
list_mat_c = spm.fit(to_send) # Receiving matrix as output - defined by goal's value
print_mat = np.array(list_mat_c) # Converting list into a NumPy array
list_mat_dimc = print_mat.shape[1]  # Deriving matrix's column dimension


# For SPK procedure - receives k-clusters matrix 
if(goal == 0):   # goal = spk
    if(k == 0):
        k = list_mat_dimc  # Updating k according to T's column dimension
    first = k_means_pp(print_mat)  # Receiving k initial centroids
    to_send = [k, N, first, list_mat_c] # List of inputs sent to C
    list_final_c = spm.kmeans_pp(to_send) # Receiving list of k clusters (after operating K-means ++)
    list_final_c = np.array(list_final_c) # Converting list into a NumPy array
    print_mat = list_final_c # Assigning printing variable

# Printing loop for every necessary output (for each goal)
i = 0
for element in print_mat.flat:
    i=i+1
    if(element < 0):
        if(-element < 0.00005):  # Matrix's element is negative and its absolute value is smaller than the threshold
            print("0.0000", end = '')
        else:
            print(f'{element:.4f}', end = '') # Formatting list's elements values to 4 decimal points
    else:
        print(f'{element:.4f}', end = '') # Formatting list's elements values to 4 decimal points
    if i == list_mat_dimc:  # End of row
        i=0
        print("\n", end='')
        continue
    print(",", end='')