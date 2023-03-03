import os
from main import main
from combine_data import combineData

# The file to run multiple simulations in a row and save them to different folders
def multiple():
    file = r"C:\Users\cassa\Desktop\Results\test"
    os.chdir(file)
    for i in range(3):
        os.mkdir(f"{i}")
        os.chdir(f"{i}")
        main()
        os.chdir(file)
    combineData()

multiple()

