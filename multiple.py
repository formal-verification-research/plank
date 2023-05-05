import os
from main import main
from combine_data import combineAll
from parameter_vault import NUMBER_OF_RUNS
import threading


# The file to run multiple simulations in a row and save them to different folders
# Change the number of run in parameter_vault --- NUMBER_OF_RUNS

def multiple():
    file = r"C:\Users\cassa\Desktop\Results\test"
    os.chdir(file)
    threads = []
    for i in range(NUMBER_OF_RUNS):
        # try:
        #     os.mkdir(f"{i}")
        # except FileExistsError:
        #     pass
        # # os.chdir(f"{i}")
        # threads.append(threading.Thread(target=main, args=(f"{i}",)))
        os.mkdir(f"{i}")
        os.chdir(f"{i}")
        main()

        # print("made a new file")
        os.chdir(file)
    # for thread in threads:
    #     thread.start()
    #
    # for thread in threads:
    #     thread.join()
    combineAll()
if __name__=="__main__":
    multiple()




