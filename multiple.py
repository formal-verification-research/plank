import os
from main import main
from combine_data import combineAll
from parameter_vault import NUMBER_OF_RUNS

# The file to run multiple simulations in a row and save them to different folders
# Change the number of run in parameter_vault --- NUMBER_OF_RUNS
def multiple():
    print("running multiple")
    print("why are you not running")
    file = r"C:\Users\cassa\Desktop\Results\test"
    os.chdir(file)
    for i in range(NUMBER_OF_RUNS):
        os.mkdir(f"{i}")
        os.chdir(f"{i}")
        main()
        print("made a new file")
        os.chdir(file)

    combineAll()

multiple()

