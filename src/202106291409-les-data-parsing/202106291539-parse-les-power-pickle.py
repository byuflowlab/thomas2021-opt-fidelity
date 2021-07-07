import numpy as np
import pickle 
import pandas as pd
import re

def pickle_to_csv(case="low_ti",layout=""):
    filename="/Users/jaredthomas/OneDrive - BYU/Documents/Jared/School/PhD/Data/thomas2021-opt-fidelity/saved_sowfa_"+case+"_byu_38turb"+layout+".p"
    data = pickle.load(open(filename, "rb"),encoding='ASCII')
    
    # set regex to search for in case name to get wind direction 
    pattern = re.compile(r'(?<=wd)\d+')

    # initialize arrays
    turbinepower = np.zeros((12,38))
    winddirections = np.zeros(12) 
    
    for i in np.arange(0,12):
        turbinepower[i, :] = np.asarray(data["turbine_power"][i])*1E3
        winddirections[i] = float(pattern.findall(data["case_name"][i])[0])

    np.savetxt("turbine-power-low-ti"+layout+".txt", np.c_[turbinepower.T], header=np.array2string(winddirections))



if __name__ == "__main__":
    pickle_to_csv(layout="_opt")
    # pickle_to_csv(layout="")
