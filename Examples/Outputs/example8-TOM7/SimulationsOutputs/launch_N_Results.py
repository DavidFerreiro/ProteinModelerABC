import Functions
import multiprocessing

if __name__ == "__main__":
  a_pool = multiprocessing.Pool(500)
  result = a_pool.map(Functions.SummaryResultN, range(20000, 30000), chunksize=1)
  a_pool.close()
  a_pool.join()
