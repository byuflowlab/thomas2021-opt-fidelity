import numpy as np
import pickle 
import pandas as pd
import re

from scipy.stats.stats import energy_distance

def pickle_to_csv(case="low_ti",layout="",n="4"):
    if case == "low-ti":
        case2 = "low_ti"
    if case == "high-ti":
        case2 = "high_ti"
        case3 = "high_TI"
    if layout == "-opt":
        layout2 = "_opt"
    else:
        layout2 = ""

    pathname = "/Users/jaredthomas/OneDrive - Brigham Young University/Documents/Jared/School/PhD/Data/thomas2021-opt-fidelity/les-data/"+case3+"_"+n+"/"
    filename = "saved_sowfa_"+case2+"_byu_38turb"+layout2+n+".p"
    data = pickle.load(open(pathname+filename, "rb"),encoding='ASCII')
    
    # set regex to search for in case name to get wind direction 
    pattern = re.compile(r'(?<=wd)\d+')

    # initialize arrays
    turbinepower = np.zeros((12,38))
    winddirections = np.zeros(12) 
    
    for i in np.arange(0,12):
        turbinepower[i, :] = np.asarray(data["turbine_power"][i])*1E3
        winddirections[i] = float(pattern.findall(data["case_name"][i])[0])

    np.savetxt("turbine-power-"+case+"-les"+layout+n+".txt", np.c_[turbinepower.T], header=np.array2string(winddirections))



if __name__ == "__main__":
    pickle_to_csv(case="high-ti", layout="-opt", n="4")
    # pickle_to_csv(case="low-ti", layout="")
    # pickle_to_csv(case="low-ti", layout="-opt")
