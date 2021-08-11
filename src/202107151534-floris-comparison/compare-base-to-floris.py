# Copyright 2021 NREL

# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License. You may obtain a copy of
# the License at http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations under
# the License.

# See https://floris.readthedocs.io for documentation


import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import floris.tools as wfct
import floris.tools.cut_plane as cp
import floris.tools.wind_rose as rose
import floris.tools.power_rose as pr
import floris.tools.visualization as vis
from floris.tools.optimization.scipy.yaw_wind_rose import YawOptimizationWindRose


# Instantiate the FLORIS object
file_dir = os.path.dirname(os.path.abspath(__file__))
fi = wfct.floris_interface.FlorisInterface(
    os.path.join(file_dir, "example_input.json")
)

# Define wind farm coordinates and layout
wf_coordinate = [39.8283, -98.5795]

# Below minimum wind speed, assumes power is zero.
minimum_ws = 3.0

# Set wind farm to N_row x N_row grid with constant spacing
# (2 x 2 grid, 5 D spacing)
D = 126.4
layoutdata = np.loadtxt("../inputfiles/farms/layout_38turb_round.txt")*D
layout_x = layoutdata[:, 0]
layout_y = layoutdata[:, 1]

N_turb = len(layout_x)

fi.reinitialize_flow_field(
    layout_array=(layout_x, layout_y), wind_direction=[270.0], wind_speed=[8.0]
)
fi.calculate_wake()

# ================================================================================
print("Plotting the FLORIS flowfield...")
# ================================================================================

# Initialize the horizontal cut
hor_plane = fi.get_hor_plane(height=fi.floris.farm.turbines[0].hub_height)

# Plot and show
fig, ax = plt.subplots()
wfct.visualization.visualize_cut_plane(hor_plane, ax=ax)
ax.set_title("Baseline flow for U = 8 m/s, Wind Direction = 270$^\\circ$")

# ================================================================================
print("Importing wind rose data...")
# ================================================================================

# Create wind rose object and import wind rose dataframe using WIND Toolkit
# HSDS API. Alternatively, load existing file with wind rose information.
calculate_new_wind_rose = False

wind_rose = rose.WindRose()

if calculate_new_wind_rose:

    wd_list = np.arange(0, 360, 5)
    ws_list = np.arange(0, 26, 1)

    df = wind_rose.import_from_wind_toolkit_hsds(
        wf_coordinate[0],
        wf_coordinate[1],
        ht=100,
        wd=wd_list,
        ws=ws_list,
        limit_month=None,
        st_date=None,
        en_date=None,
    )

else:
    # df = wind_rose.load(
    #     os.path.join(file_dir, "windtoolkit_geo_center_us.p")
    # )
    winddata = np.loadtxt("../inputfiles/wind/windrose_nantucket_12dir.txt")
    df = pd.DataFrame({
            'wd': winddata[:,0],
            "ws": winddata[:,1]*0 + 8.055,
            "freq_val": winddata[:,2]
            })
                                                    
print(df.head())
print(df.wd)
wind_rose.num_wd = len(df.wd)
wind_rose.num_ws = 1,
wind_rose.wd_step = df.wd[1]-df.wd[0],
wind_rose.ws_step = 0,
wind_rose.wd = df.wd
wind_rose.ws = df.ws
wind_rose.df = df

print('Set specific model parameters on the current wake model:\n')
alphastar = 2.32
betastar = 0.154
k1 = 0.3837
k2 = 0.003678

alpha = alphastar/4.0
beta = betastar/2.0

ti = 0.061

params = {
    'Wake Velocity Parameters': {'alpha': alpha},
    'Wake Velocity Parameters': {'beta': beta},
    'Wake Velocity Parameters': {'ka': k1},
    'Wake Velocity Parameters': {'kb': k2},
    'Wake Velocity Parameters': {"calculate_VW_velocities": False},
    'Wake Velocity Parameters': {"use_yaw_added_recovery": False},
    'Wake Deflection Parameters': {'alpha': alpha},
    'Wake Deflection Parameters': {'beta': beta},
    'Wake Deflection Parameters': {'ka': k1},
    'Wake Deflection Parameters': {'kb': k2},
    'Wake Turbulence Parameters': {'ti_initial': ti}
}

print('\n')

fi.set_model_parameters(params)

print('Observe that the requested paremeters changes have been made:\n')
model_params = fi.get_model_parameters()
print(model_params)

fi.show_model_parameters()
print('\n')

# plot wind rose
wind_rose.plot_wind_rose()

# =============================================================================
print("Finding power with and without wakes in FLORIS...")
# =============================================================================
print(df.ws)
print(df.freq_val)
AEP = fi.get_farm_AEP(df.wd, df.ws, df.freq_val, limit_ws=False)

# Instantiate the Optimization object
# Note that the optimization is not performed in this example.
# yaw_opt = YawOptimizationWindRose(fi, df.wd, df.ws, minimum_ws=minimum_ws)

# # Determine baseline power with and without wakes
# df_base = yaw_opt.calc_baseline_power()

# # Initialize power rose
# case_name = "Round Wind Farm"
# power_rose = pr.PowerRose()
# power_rose.make_power_rose_from_user_data(
#     case_name, df, df_base["power_no_wake"], df_base["power_baseline"]
# )

# # Display AEP analysis
# fig, axarr = plt.subplots(2, 1, sharex=True, figsize=(6.4, 6.5))
# power_rose.plot_by_direction(axarr)
# power_rose.report()
print("AEP: ", AEP*1E-9)

# plt.show()