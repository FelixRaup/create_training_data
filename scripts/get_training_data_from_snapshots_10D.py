### --------- loading packages / modules ------------ ###
# import h5py
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import joblib
import math

from flchem import flchem_evolve, flchem_evolve_mp, flchem_init
from timeit import default_timer as timer
from multiprocessing import Pool, freeze_support, set_start_method

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler


# load functions from external script for cleaner script
from pychem_toolkit import get_config_from_snapshot, pychem_data_generator, balance_by_histogram
from pychem_toolkit import get_random_config, transform, inverse_transform


### ------------- settings --------------------------- ###


N = int(1e3)
stand_scaler = True
include_existing_sample = False
training_folder = "../data/test" # "data_hist_balanced_10D_longtime"
input_size = 10
save_physical = True
plot = True

balance = False
target_count = 10000
min_threshold = 2000

consider_deltas = False

snapshots_path = "../snapshots/"

if not os.path.exists(training_folder):
    os.makedirs(training_folder)


### ---------- constants for pychem ------------- ###

kyr = 3.1556926e10

# totaltime, timestep = 3.1556926e10, 3.1556926e10
abundc, abundo = 1.4e-4, 3.2e-4
dust_to_gas_ratio, dust_temp = 1e0, np.nan
G0, cosmic_ray_ion_rate = 1.7, 3e-17
NH_ext, Z_atom = 1e20, 1.
dl, divv = 1e15, 0.
abcp = np.nan
G0RP = 0.0

ch_muC = 12.011
ch_mf_scale = 1.0 + abundc * ch_muC

eps = 1e-40

### ------- input parameter setting -------- ###

from dataset_config_king import input_range, input_scaling, output_scaling


### ------------- CREATE NEW TRAINING POINTS -------------- ###

X_physical = get_config_from_snapshot(snapshots_path, N)


X_physical, y_physical = pychem_data_generator(
    X_physical = X_physical, 
    filter_data = False,
    abundance_thresholds = [1e-10, 1e-15, 1e-21],
    low_threshold = np.nan)
if consider_deltas:
    y_physical[:, 0] = y_physical[:, 0] - X_physical[:, 6]  # abh2
    y_physical[:, 1] = y_physical[:, 1] - X_physical[:, 7]  # abhp
    y_physical[:, 2] = y_physical[:, 2] - X_physical[:, 8]  # abco
    y_physical[:, 4] = y_physical[:, 4] - X_physical[:, 5]  # eint


delta_eint = y_physical[:, 4]

print(np.percentile(
    np.abs(y_physical),
    [1, 10, 25, 50, 75, 90, 99]
))




# load existing sample and combine
if include_existing_sample:
    
    data = np.load(training_folder + "/raw_sample.npz")
    X_old = data["X_physical"]
    y_old = data["y_physical"]

    print(X_physical.shape)

    X_physical = np.concatenate( (X_physical, X_old), axis=0 )
    y_physical = np.concatenate( (y_physical, y_old), axis=0 )

    print(X_physical.shape)


print("input shape = " , X_physical.shape)


### -------- SCALING OF SAMPLE (e.g. log10) -------- ###

X_scaled = transform(X_physical, input_scaling)
y_scaled = transform(y_physical, output_scaling)


# for i, name in enumerate(output_scaling.keys()):
#     print(
#         f"{i:2d} {name:5s} | "
#         f"physical: {y_physical[0, i]: .6e} | "
#         f"scaled: {y_scaled[0, i]: .6e}"
#     )

# exit()


### -------------- BALANCE SAMPLE ------------------- ### 

if balance:
    X_scaled_balanced, y_scaled_balanced = balance_by_histogram(X_scaled, y_scaled, bins=50, target_count=target_count, min_threshold=min_threshold, dimensions=[0, 1, 2, 4])
else:
    X_scaled_balanced, y_scaled_balanced = X_scaled, y_scaled


print(X_scaled.shape)
print(X_scaled_balanced.shape)


X_balanced = inverse_transform(
    X_scaled_balanced,
    input_scaling
)

y_balanced = inverse_transform(
    y_scaled_balanced,
    output_scaling
)


### ---------------- PLOT SAMPLE ----------------- ###


# print(y_physical == y_balanced)

# exit()

if plot:

    datasets = [X_physical, X_scaled, X_scaled_balanced, X_balanced]
    titles = ["X physical", "X scaled", "X scaled balanced", "X physical balanced"]
    quantity = ["fshield_H2", "fshield_CO", "AV_mean", "chi_mean", "dens", "eint", "abh2","abhp","abco","time"]

    fig = plt.figure(figsize=(14, 12))
    outer = gridspec.GridSpec(2, 2, wspace=0.25, hspace=0.25)

    for idx, data in enumerate(datasets):
        inner = gridspec.GridSpecFromSubplotSpec(
            data.shape[1], 1,
            subplot_spec=outer[idx],
            hspace=0.1
        )
        
        for i in range(data.shape[1]):
            ax = plt.Subplot(fig, inner[i])
            ax.hist(data[:, i], bins=50, edgecolor='black')
            ax.set_ylabel(quantity[i])
            ax.grid(True, alpha=0.5)
            ax.set_ylim(bottom=0)
            ax.set_xticks([])
            ax.set_xticklabels([])

            
            fig.add_subplot(ax)

        fig.axes[-data.shape[1]].set_title(titles[idx])

    plt.tight_layout()
    plt.show()

    fig.savefig(training_folder + "/input_sampling")



    datasets = [y_physical, y_scaled, y_scaled_balanced, y_balanced]
    titles = ["y physical", "y scaled", "y scaled balanced", "y physical balanced"]
    quantity = ["h2", "hp", "co", "T_dust", "eint"]

    fig = plt.figure(figsize=(14, 12))
    outer = gridspec.GridSpec(2, 2, wspace=0.25, hspace=0.25)

    for idx, data in enumerate(datasets):
        inner = gridspec.GridSpecFromSubplotSpec(
            data.shape[1], 1,
            subplot_spec=outer[idx],
            hspace=0.1
        )
        
        for i in range(data.shape[1]):
            ax = plt.Subplot(fig, inner[i])
            ax.hist(data[:, i], bins=50, edgecolor='black')
            ax.set_ylabel(quantity[i])
            ax.grid(True, alpha=0.5)
            
            # ax.set_xticks([])
            # ax.set_xticklabels([])
            # ax.set_yscale('log')
            ax.set_ylim(bottom=0)
 
            fig.add_subplot(ax)

        fig.axes[-data.shape[1]].set_title(titles[idx])

    plt.tight_layout()
    plt.show()

    fig.savefig(training_folder + "/output_sampling")



### -------------- SAVE DATA ------------------------ ### 



np.savez_compressed(
    training_folder + "/raw_sample.npz",
    X_physical=X_balanced,
    y_physical=y_balanced
)
print(f"balanced raw sample saved in {training_folder}")


X_train_phys, X_test_phys, y_train_phys, y_test_phys = train_test_split(
    X_balanced, y_balanced, test_size=0.2, random_state=42
)

X_train_phys_mins = X_train_phys.min(axis=0)
y_train_phys_mins = y_train_phys.min(axis=0)

X_train_phys_max = X_train_phys.max(axis=0)
y_train_phys_max = y_train_phys.max(axis=0)


if save_physical:

    X_train_phys, X_test_phys, y_train_phys, y_test_phys = train_test_split(
        X_balanced, y_balanced, test_size=0.2, random_state=42
    )

    N = X_train_phys.shape[0]
    new_col =  np.arange(1, N+1).reshape(N, 1)
    training_data = np.hstack((new_col, X_train_phys, y_train_phys))

    N = X_test_phys.shape[0]
    new_col =  np.arange(1, N+1).reshape(N, 1)
    test_data = np.hstack((new_col, X_test_phys, y_test_phys))

    if not os.path.exists(training_folder):
        os.makedirs(training_folder)

    np.savez_compressed(
        training_folder + "/training_data_physical.npz",
        data=training_data,
    )

    np.savez_compressed(
        training_folder + "/test_data_physical.npz",
        data=test_data,
    )

    exit()