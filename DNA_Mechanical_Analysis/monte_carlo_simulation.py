import pandas as pd
import numpy as np
import random
from itertools import islice
import matplotlib.pyplot as plt

phase_factor= (145*np.pi)/10
lam=4.46
Number_iterations = 1000

stiff_table=pd.read_csv("stiffness_intrinsic_values.csv")

###Generate the random dna sequence
def random_dna_sequence(length):
    return ''.join(random.choice('ACTG') for _ in range(length))

# print (random_dna_sequence(147))

def mutated_sequence_gen(s):
    location = random.randrange(0,147)
    mutation_letter=random.choice('ACTG')
    s_list = list(s)
    s_list[location]= mutation_letter
    return ''.join(s_list), location, mutation_letter

def theta(steps):
    theta_tilt = lam*np.sin(((2*np.pi*steps)/10)-phase_factor)
    theta_roll = lam*np.cos(((2*np.pi*steps)/10)-phase_factor)
    theta_twist = 35.575
    return theta_tilt, theta_roll, theta_twist


def K_intrinsic(dinucl):
    parameter_values = stiff_table[stiff_table["dinucleotide"] == dinucl]

    return parameter_values["theta_tilt_dash"].values[0], parameter_values["K_tilt"].values[0], \
           parameter_values["theta_roll_dash"].values[0], parameter_values["K_roll"].values[0]


def T_dinucl_step(step, dinucleotide):

    beta = 1 # We take just constant because the E value in T is already KbT units.
                #beta value other than means we are taking different temperature than E value

    ### get the step values of tilt and roll
    theta_t, theta_r, theta_tw = theta(step)

    #### get the intrinsic values for a basepair step
    theta_t_dash, K_tilt, theta_r_dash, K_roll = K_intrinsic(dinucleotide)
#     print (theta_t - theta_t_dash)

    E_tilt = 0.5*K_tilt*((theta_t - theta_t_dash)**2)

    E_roll = 0.5*K_roll*((theta_r - theta_r_dash)**2)

    E = E_tilt + E_roll
    T = np.exp(-beta*E)

    return T

def window(seq, n=2):
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result


def next_sequence(e_delta, m_sequence, o_sequence):
    if e_delta <= 0:
        return e_delta, m_sequence

    elif np.exp(-e_delta) >= random.random():
        return e_delta, m_sequence
    else:
        return -99999, o_sequence


def gen_delta_E(m_sequence, all_steps, p):
    mutated_sequence_list = ["".join(x) for x in window(m_sequence, 2)]

    if p == 0:
        E_mutated = T_df.loc[T_df['step'] == p, mutated_sequence_list[p]].values[0]
        E_sequence = T_df.loc[T_df['step'] == p, all_steps[p]].values[0]

    elif p == 146:
        E_mutated = T_df.loc[T_df['step'] == p - 1, mutated_sequence_list[p - 1]].values[0]
        E_sequence = T_df.loc[T_df['step'] == p - 1, all_steps[p - 1]].values[0]

    else:

        E_mutated = T_df.loc[T_df['step'] == p - 1, mutated_sequence_list[p - 1]].values[0] + \
                    T_df.loc[T_df['step'] == p, mutated_sequence_list[p]].values[0]

        E_sequence = T_df.loc[T_df['step'] == p - 1, mutated_sequence_list[p - 1]].values[0] + \
                     T_df.loc[T_df['step'] == p, mutated_sequence_list[p]].values[0]

    d_E = E_mutated - E_sequence
    return d_E

########################################## DEFINE THE ENERGY MATRIX

all_dinucleotides = ["AA", "AT", "AC", "AG",
                     "TA", "TT", "TC", "TG",
                     "CA", "CT", "CC", "CG",
                     "GA", "GT", "GC", "GG"]
data = []
temp_dict = dict()
# step =1

for step in range(0, 146):
    for di in all_dinucleotides:
        temp_dict[di] = T_dinucl_step(step, di)
        temp_dict["step"] = step
    data.append(temp_dict)
    temp_dict = dict()

T_df = pd.DataFrame(data)

################################################

def stable_sequence():
    Energy_bucket = []
    seq_number = []
    ind = 0
    E = 0

    for i in range(0, Number_iterations):

        if i == 0:

            sequence = random_dna_sequence(147)
            all_steps = ["".join(x) for x in window(sequence, 2)]

            count = 0
            E = 0
            for s in all_steps:
                E = E + T_dinucl_step(count, s)
                count = count + 1

            mutated_sequence, position, letter = mutated_sequence_gen(sequence)

            delta_E = gen_delta_E(mutated_sequence, all_steps, position)

            DELTA_E, sequence = next_sequence(delta_E, mutated_sequence, sequence)


        else:

            all_steps = ["".join(x) for x in window(sequence, 2)]

            mutated_sequence, position, letter = mutated_sequence_gen(sequence)

            delta_E = gen_delta_E(mutated_sequence, all_steps, position)

            DELTA_E, sequence = next_sequence(delta_E, mutated_sequence, sequence)

        if DELTA_E != -99999:
            E = E + DELTA_E
            Energy_bucket.append(E)
            seq_number.append(ind)
            ind = ind + 1
    return sequence, Energy_bucket, seq_number


if __name__ == "__main__":


    stable_seq, Energy_bucket, seq_number = stable_sequence()

    import seaborn as sns
    fig = plt.figure(figsize=(20, 10))

    # sns.lineplot(data=Energy_bucket).
    sns.scatterplot(x=seq_number, y=Energy_bucket)
    plt.show()