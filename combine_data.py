# this file combines all the collected data from each run and makes graphs

import pandas as pd
import glob
import sys


def combineData():
    file = r"C:\Users\cassa\Desktop\Results\test"
    files = (glob.glob(file+"/**/Density"))
    tempfile = open(file, "w")
    number = 0
    for file in files:

        file = pd.read_csv(file, delim_whitespace=True)

        f = pd.DataFrame(file)
        with pd.ExcelWriter(r"C:\Users\cassa\Desktop\Results\test\test.xlsx", mode="a", engine="openpyxl") as writer:
            f.to_excel(writer, sheet_name=str(number))
        number += 1

def combine2():
    file = r"C:\Users\cassa\Desktop\Results\test"
    files = (glob.glob(file + "/**/Density"))
    # tempfile = open(file, "w")
    allData = {}
    fileIndex = 0
    for fl in files:
        openFile = open(fl, "r")
        for line in openFile.readlines():
            print(f"Read line {line}")
            line = line.split("\t")
            if len(line) != 2:
                print(f"Line {line} is misformatted! Ignoring line.", file=sys.stderr)
                continue
            time = line[0]
            density = line[1]
            if not time in allData:
                allData[time] = ["" for i in range(fileIndex + 1)]
                allData[time][fileIndex] = density
            else:
                currentRow = allData[time]
                while len(currentRow) <= fileIndex:
                    currentRow.append("")
                currentRow[fileIndex] = density
        fileIndex += 1
        openFile.close()
    headers = ["Time (Hours)"]
    for f in files:
        headers.append(str(f))
    data = [headers]
    for time, densityRow in allData.items():
        currentRow = [time]
        for density in densityRow:
            currentRow.append(density)
        data.append(currentRow)
    frame = pd.DataFrame(data)
    with pd.ExcelWriter(r"C:\Users\cassa\Desktop\Results\test\test_other.xlsx", mode="w", engine="openpyxl") as writer:
        frame.to_excel(writer, sheet_name="Aggregate Data")






# combineData()
combine2()
