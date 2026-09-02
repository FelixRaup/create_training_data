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

N = int(1e5)
stand_scaler = True
include_existing_sample = False
training_folder = "../data/data_hist_balanced_10D_test"
input_size = 10
save_physical = True
plot = True

consider_deltas = False


### --- data generator settings ---------- ###

from dataset_config_4 import input_range, input_scaling, output_scaling

filter_data = False
abundance_thresholds = [1e-12, 1e-15, 1e-21]
# abundance_thresholds = [1e-40, 1e-40, 1e-40]
# min_output_values = [np.nan]
min_output_values = [1e-12, 1e-12, 1e-20, 1e-20, 1e-20]
cut_1e40s = True


### ------ data balancing settings ------- ###

balance = True
target_count = 100000
min_threshold = 20000
filter_dimensions = [0, 1, 2, 4]


if not os.path.exists(training_folder):
    os.makedirs(training_folder)



### ------------- CREATE NEW TRAINING POINTS -------------- ###

X_physical = get_random_config(
    param_ranges=input_range,
    scalings=input_scaling,
    N=N,
)

X_physical, y_physical = pychem_data_generator(
    X_physical = X_physical, 
    filter_data = filter_data,
    abundance_thresholds = abundance_thresholds,
    min_output_values = min_output_values,
    cut_1e40s = cut_1e40s,
    )


if consider_deltas:
    y_physical[:, 0] = y_physical[:, 0] - X_physical[:, 6]  # abh2
    y_physical[:, 1] = y_physical[:, 1] - X_physical[:, 7]  # abhp
    y_physical[:, 2] = y_physical[:, 2] - X_physical[:, 8]  # abco
    y_physical[:, 4] = y_physical[:, 4] - X_physical[:, 5]  # eint



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
    X_scaled_balanced, y_scaled_balanced = balance_by_histogram(X_scaled, y_scaled, bins=50, target_count=target_count, min_threshold=min_threshold, dimensions=filter_dimensions)
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



### ----------- SAVE LATENT TRAINING DATA AND SCALING PARAMETERS ------- ###

X_train, X_test, y_train, y_test = train_test_split(
    X_physical_balance_log, y_physical_balance_log, test_size=0.2, random_state=42
)

if (stand_scaler == True):

    scaler_X = StandardScaler()
    scaler_y = StandardScaler()

    X_train_latent = scaler_X.fit_transform(X_train)
    X_test_latent  = scaler_X.transform(X_test)

    y_train_latent = scaler_y.fit_transform(y_train)
    y_test_latent  = scaler_y.transform(y_test)


# X_train_log = X_train_phys.copy()
# X_test_log = X_test_phys.copy()
# y_train_log = y_train_phys.copy()
# y_test_log = y_test_phys.copy()

# X_train_phys_mins = X_train_phys.min(axis=0)
# y_train_phys_mins = y_train_phys.min(axis=0)

# X_train_phys_max = X_train_phys.max(axis=0)
# y_train_phys_max = y_train_phys.max(axis=0)

# print("X_train mins = ", X_train_phys_mins)
# print("y_train mins = ", y_train_phys_mins)

# print("X_train max = ", X_train_phys_max)
# print("y_train max = ", y_train_phys_max)


# if (stand_scaler == True):

#     scaler_X = StandardScaler()
#     scaler_y = StandardScaler()


#     X_train_latent = scaler_X.fit_transform(X_train_log)
#     X_test_latent  = scaler_X.transform(X_test_log)

#     y_train_latent = scaler_y.fit_transform(y_train_log)
#     y_test_latent  = scaler_y.transform(y_test_log)





### ------------------ save data and scalers ---------------------- ###

N = X_train_latent.shape[0]
new_col =  np.arange(1, N+1).reshape(N, 1)
training_data = np.hstack((new_col, X_train_latent, y_train_latent))

# print(training_data.shape)

N = X_test_latent.shape[0]
new_col =  np.arange(1, N+1).reshape(N, 1)
test_data = np.hstack((new_col, X_test_latent, y_test_latent))

# print(test_data.shape)
# print(np.max(training_data[:1,:]))
# print(np.max(test_data[:1,:]))


if not os.path.exists(training_folder):
    os.makedirs(training_folder)

np.savez_compressed(
    training_folder + "/training_data.npz",
    data=training_data,
)

np.savez_compressed(
    training_folder + "/test_data.npz",
    data=test_data,
)



print(f"training-data saved in '{training_folder}'")


### save scalers and log info in txt file ###


# ======================================================
# Extract scaler parameters
# ======================================================

X_mean = scaler_X.mean_
X_std  = scaler_X.scale_

y_mean = scaler_y.mean_
y_std  = scaler_y.scale_

# ======================================================
# Build binary log masks
# ======================================================

X_log_mask = np.zeros(X_mean.shape[0], dtype=int)
y_log_mask = np.zeros(y_mean.shape[0], dtype=int)

X_log_mask[X_log_idx] = 1
y_log_mask[y_log_idx] = 1




# ======================================================
# Save everything
# ======================================================

with open(training_folder + "/scalers.txt", "w") as f:

    # --------------------------------------------------
    # X scaler
    # --------------------------------------------------

    f.write("X_LOG\n")
    f.write(" ".join(map(str, X_log_mask)) + "\n")

    f.write("X_MIN\n")
    f.write(" ".join(map(str, X_train_phys_mins)) + "\n")

    f.write("X_MAX\n")
    f.write(" ".join(map(str, X_train_phys_max)) + "\n")

    f.write("X_MEAN\n")
    f.write(" ".join(map(str, X_mean)) + "\n")

    f.write("X_STD\n")
    f.write(" ".join(map(str, X_std)) + "\n")

    # --------------------------------------------------
    # Y scaler
    # --------------------------------------------------

    f.write("Y_LOG\n")
    f.write(" ".join(map(str, y_log_mask)) + "\n")

    f.write("Y_MIN\n")
    f.write(" ".join(map(str, y_train_phys_mins)) + "\n")

    f.write("Y_MAX\n")
    f.write(" ".join(map(str, y_train_phys_max)) + "\n")

    f.write("Y_MEAN\n")
    f.write(" ".join(map(str, y_mean)) + "\n")

    f.write("Y_STD\n")
    f.write(" ".join(map(str, y_std)) + "\n")


print("Export complete -> scalers.txt")



###########################################
########### CHECK SKALING #################
###########################################


def load_scalers_txt(filename):

    data = {}

    with open(filename, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    i = 0
    while i < len(lines):
        key = lines[i]
        values = np.array(list(map(float, lines[i + 1].split())))
        data[key] = values
        i += 2

    return data



def latent_to_physical(X_latent, scaler_file):

    params = load_scalers_txt(scaler_file)

    X_mean = params["X_MEAN"]
    X_std = params["X_STD"]
    X_log_mask = params["X_LOG"].astype(bool)

    # Undo StandardScaler
    X_log = X_latent * X_std + X_mean

    # Undo log10 transform
    X_phys = X_log.copy()

    X_phys[:, X_log_mask] = (
        10.0 ** X_log[:, X_log_mask]
    ) - eps

    return X_phys

# ======================================================
# Verify inverse transform
# ======================================================

def check_inverse_scaling(
        X_train_phys,
        X_train_latent,
        scaler_file,
        atol=1e-12,
        rtol=1e-8):

    X_recovered = latent_to_physical(
        X_train_latent,
        scaler_file
    )

    abs_err = np.max(np.abs(X_recovered - X_train_phys))

    rel_err = np.max(
        np.abs(X_recovered - X_train_phys)
        / np.maximum(np.abs(X_train_phys), 1e-30)
    )

    print("\nInverse scaling test")
    print("--------------------")
    print("max abs error =", abs_err)
    print("max rel error =", rel_err)

    passed = np.allclose(
        X_recovered,
        X_train_phys,
        atol=atol,
        rtol=rtol
    )

    print("PASS =", passed)

    return passed, X_recovered

passed, X_recovered = check_inverse_scaling(
    X_train_phys,
    X_train_latent,
    training_folder + "/scalers.txt"
)

