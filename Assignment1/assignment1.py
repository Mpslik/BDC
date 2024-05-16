import multiprocessing as mp

if __name__ == "__main__":
    cpus = mp.cpu_count()
    with mp.Pool(cpus) as pool:
        results = pool.map()
    print(results)