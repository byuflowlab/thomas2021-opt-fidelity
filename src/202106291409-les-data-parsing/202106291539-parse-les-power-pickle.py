import numpy as np
import pickle 
import pandas as pd

def pickle_to_csv():
    filename="/Users/jaredthomas/OneDrive - BYU/Documents/Jared/School/PhD/Data/thomas2021-opt-fidelity/saved_sowfa_low_ti_byu_38turb.p"
    data = pickle.load(open(filename, "rb"),encoding='ASCII')
    print(type(data))
    data.head()
    data.to_csv("turbine-power-low-ti.txt")
    print(data["turbine_power"])



if __name__ == "__main__":
    pickle_to_csv()
