import numpy as np
import pickle 
import pandas as pd
import re

def pickle_to_csv():
    filename="/Users/jaredthomas/OneDrive - BYU/Documents/Jared/School/PhD/Data/thomas2021-opt-fidelity/saved_sowfa_low_ti_byu_38turb.p"
    data = pickle.load(open(filename, "rb"),encoding='ASCII')
    
    # set regex to search for in case name to get wind direction 
    pattern = re.compile(r'(?<=wd)\d+')

    # initialize arrays
    turbinepower = np.zeros((12,38))
    winddirections = np.zeros(12) 
    
    for i in np.arange(0,12):
        turbinepower[i, :] = np.asarray(data["turbine_power"][i])
        winddirections[i] = float(pattern.findall(data["case_name"][i])[0])

    np.savetxt("turbine-power-low-ti.txt", np.c_[turbinepower.T], header=np.array2string(winddirections))



if __name__ == "__main__":
    pickle_to_csv()
