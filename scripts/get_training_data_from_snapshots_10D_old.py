### --------- loading packages / modules ------------ ###
import h5py
import numpy as np
import os
import matplotlib.pyplot as plt
import joblib
import matplotlib.gridspec as gridspec

from flchem import flchem_evolve, flchem_evolve_mp, flchem_init
from timeit import default_timer as timer
from multiprocessing import Pool, freeze_support

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from pychem_toolkit import get_config_from_snapshot, pychem_data_generator, balance_by_histogram
from pychem_toolkit import get_random_config, transform, inverse_transform

### ------------- settings --------------------------- ###

snapshots_path = "../snapshots/"


N = int(1e2)
stand_scaler = True
include_existing_sample = False
training_folder = "data_snapshot_10D"
input_size = 10
save_physical = True
plot = True

if not os.path.exists(training_folder):
    os.makedirs(training_folder)


eps = 1e-40

param_log_scales = {
    "fshield_H2": True,
    "fshield_CO": True,
    "AV_mean": True,
    "chi_mean": False,
    "dens": True,
    "eint": True,
    "abh2": True,
    "abhp": True,
    "abco": True,
    "time": False,
    }



### ------------- CREATE NEW TRAINING POINTS -------------- ###

X_physical = get_config_from_snapshot(snapshots_path, N)


X_physical, y_physical = pychem_data_generator(
    X_physical = X_physical, 
    filter_data = False,
    abundance_thresholds = [1e-10, 1e-15, 1e-21],
    low_threshold = np.nan)


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

X_log_idx = [i for i, v in enumerate(param_log_scales.values()) if v]
print(X_log_idx)
y_log_idx = [1, 4] # linear dust temp

X_physical_log = X_physical.copy()
X_physical_log[:,X_log_idx] = np.log10(X_physical[:,X_log_idx] + eps) 

y_physical_log = y_physical.copy()
y_physical_log[:,y_log_idx] = np.log10(y_physical[:,y_log_idx] + eps)

y_physical_log[:, 0] = y_physical[:, 0]**0.2
y_physical_log[:, 2] = y_physical[:, 2]**0.2



### -------------- BALANCE SAMPLE ------------------- ### 

X_physical_balance_log, y_physical_balance_log = balance_by_histogram(X_physical_log, y_physical_log, bins=50, target_count=10000, min_threshold=2000, dimensions=[0, 1, 2, 4])

print(X_physical_log.shape)
print(X_physical_balance_log.shape)


X_physical_balance = X_physical_balance_log.copy()
X_physical_balance[:,X_log_idx] = 10**(X_physical_balance_log[:,X_log_idx]) - eps

# reverse scaling 
y_physical_balance = y_physical_balance_log.copy()
y_physical_balance[:,y_log_idx] = 10**(y_physical_balance_log[:,y_log_idx]) - eps

y_physical_balance[:, 0] = y_physical_balance_log[:,0]**5
y_physical_balance[:, 2] = y_physical_balance_log[:,2]**5



### ---------------- PLOT SAMPLE ----------------- ###

if plot:

    datasets = [y_physical, y_physical_log, y_physical_balance_log, y_physical_balance]
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
            
            ax.set_xticks([])
            ax.set_xticklabels([])
            # ax.set_yscale('log')
            ax.set_ylim(bottom=0)
 
            fig.add_subplot(ax)

        fig.axes[-data.shape[1]].set_title(titles[idx])

    plt.tight_layout()
    plt.show()

    fig.savefig(training_folder + "/output_sampling")



np.savez_compressed(
    training_folder + "/raw_sample.npz",
    X_physical=X_physical_balance,
    y_physical=y_physical_balance
)
print(f"balanced raw sample saved in {training_folder}")


X_train_phys, X_test_phys, y_train_phys, y_test_phys = train_test_split(
    X_physical_balance, y_physical_balance, test_size=0.2, random_state=42
)

X_train_phys_mins = X_train_phys.min(axis=0)
y_train_phys_mins = y_train_phys.min(axis=0)

X_train_phys_max = X_train_phys.max(axis=0)
y_train_phys_max = y_train_phys.max(axis=0)



if save_physical:

    X_train_phys, X_test_phys, y_train_phys, y_test_phys = train_test_split(
        X_physical_balance, y_physical_balance, test_size=0.2, random_state=42
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






