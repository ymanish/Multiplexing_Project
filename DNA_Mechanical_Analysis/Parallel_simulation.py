import time
import monte_carlo_simulation
import concurrent.futures

f = open('sequences_energy.txt', 'w')

def simulation():

    seq, E_bucket, s_number = monte_carlo_simulation.stable_sequence()
    return seq+"_"+str(E_bucket[-1])+"\n"


if __name__ == "__main__":

    start = time.perf_counter()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(10)
        pool = [executor.submit(simulation) for i in iter_seq]
        for i in concurrent.futures.as_completed(pool):
            print(f'Return Value: {i.result()}')
            f.write(i.result())

    f.close()
    end = time.perf_counter()
    print(f'Finished in {round(end - start, 2)} second(s)')