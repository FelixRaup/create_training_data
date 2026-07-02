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




### ------------- settings --------------------------- ###


N = int(1e5)
stand_scaler = True
include_existing_sample = True
training_folder = "data_hist_balanced_10D_large"
input_size = 10
save_physical = True
plot = True


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

param_ranges = {
    "fshield_H2": (1e-10, 1),
    "fshield_CO": (1e-10, 1),
    "AV_mean": (1e-5, 1e3), # 1e-3
    "chi_mean": (1e-5, 1),
    "dens": (1e-26, 1e-16),
    "eint": (1e-15, 1e-5),
    "abh2": (1e-10, 0.5),
    "abhp": (1e-15, 1),
    "abco": (1e-21, abundc),
    "time": (0.1 * kyr, 20 * kyr)
    }


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

param_alphas = {
    "fshield_H2": 0.1,
    "fshield_CO": 0.1,
    "AV_mean": 1,
    "chi_mean": 0,
    "dens": 1,
    "eint": 1,
    "abh2": 0.1,
    "abhp": 1,
    "abco": 0.1,
    "time": 1,
    }


### ------------------------------------------------ ### 



if not os.path.exists(training_folder):
    os.makedirs(training_folder)

def get_random_arr(range = (1,10), N=int(1e3), log=True):

    if (log==True):
        return 10**np.random.uniform(np.log10(range[0]), np.log10(range[1]), N)
    
    else:
        return np.random.uniform(range[0], range[1], N)
    

def get_random_arr_log_uniform(range=(1, 10), N=int(1e3), alpha=1.0):
    """
    alpha=0.0 → purely linear
    alpha=1.0 → purely logarithmic
    """
    u = np.random.uniform(0, 1, N)  # single draw

    # Linear and log mappings of the same u
    linear = range[0] + u * (range[1] - range[0])
    log    = range[0] * (range[1] / range[0]) ** u

    # alpha in log space
    return np.exp((1 - alpha) * np.log(linear) + alpha * np.log(log))

def get_variable_base_samples(range_val=(1, 10), N=1000, b=10.0):
    """
    Forward Function: Generates samples using a specific log-base 'b'.
    b > 1: Clustered towards the lower bound (like alpha > 0)
    """
    u = np.random.uniform(0, 1, N)
    r0, r1 = range_val
    
    # Scale u using the variable base b
    return r0 + (r1 - r0) * (b**u - 1) / (b - 1)

def invert_variable_base(x, range_val=(1, 10), b=10.0):
    """
    Inverse Function: Instantly maps physical values (x) back to 
    the uniform latent space [0, 1] using the variable base 'b'.
    """
    r0, r1 = range_val
    
    # Pure algebraic inversion (vectorized for speed)
    numerator = np.log(1 + ((x - r0) / (r1 - r0)) * (b - 1))
    denominator = np.log(b)
    
    return numerator / denominator
    
# def balance_by_histogram(X, y, bins=50, target_count=500):

#     indices = np.arange(len(y))

#     for dim in [0,1,2,3,4]: # [1,0,2,4]: 

#         print(f"Filtering with respect to dimension {dim}")

#         y_ref = y[indices, dim]

#         counts, bin_edges = np.histogram(y_ref, bins=bins)

#         nonzero_counts = counts[counts > 0]
#         if len(nonzero_counts) == 0:
#             continue

#         min_count = np.min(nonzero_counts)
#         #target_count = int(500) # int(min_count * multiple)

#         new_indices = []

#         for i in range(bins):
#             bin_mask = (y_ref >= bin_edges[i]) & (y_ref < bin_edges[i+1])
#             bin_indices = indices[np.where(bin_mask)[0]]

#             if len(bin_indices) == 0:
#                 continue

#             if len(bin_indices) > target_count:
#                 chosen = np.random.choice(bin_indices, target_count, replace=False)
#             else:
#                 chosen = bin_indices

#             new_indices.append(chosen)

#         if len(new_indices) == 0:
#             continue

#         indices = np.concatenate(new_indices)

#     return X[indices], y[indices]



# def balance_by_histogram(X, y, bins=50, target_count=500, min_threshold=50):
#     indices = np.arange(len(y))
#     dimensions = [0, 1, 2, 4]

#     for dim in dimensions:
#         print(f"Filtering with respect to dimension {dim}")

#         # WICHTIG: Wir betrachten die Werte der AKTUELL noch übrig gebliebenen Punkte
#         y_ref = y[indices, dim]
#         counts, bin_edges = np.histogram(y_ref, bins=bins)

#         nonzero_counts = counts[counts > 0]
#         if len(nonzero_counts) == 0:
#             continue

#         new_indices = []

#         for i in range(bins):
#             # Maske für den aktuellen Bin innerhalb der noch vorhandenen Indizes
#             bin_mask = (y_ref >= bin_edges[i]) & (y_ref < bin_edges[i+1])
#             bin_indices = indices[np.where(bin_mask)[0]]

#             current_bin_size = len(bin_indices)
#             if current_bin_size == 0:
#                 continue

#             # --- DYNAMISCHER SCHUTZ ---
#             # Wenn der Bin JETZT SCHON kleiner oder gleich dem Threshold ist, 
#             # ODER wenn das Filtern auf target_count ihn unter den Threshold drücken würde:
#             # Behalte ALLE Punkte dieses Bins unberührt!
#             if current_bin_size <= min_threshold:
#                 chosen = bin_indices
#             elif current_bin_size > target_count:
#                 # Wir filtern nur, wenn nach dem Filtern immer noch mindestens `min_threshold` übrig bleiben.
#                 # Falls target_count < min_threshold definiert wurde, sichern wir hier ab, 
#                 # dass wir maximal auf `min_threshold` runtergehen, niemals tiefer.
#                 effective_target = max(target_count, min_threshold)
                
#                 chosen = np.random.choice(bin_indices, effective_target, replace=False)
#             else:
#                 # Bin liegt zwischen min_threshold und target_count -> komplett behalten
#                 chosen = bin_indices

#             new_indices.append(chosen)

#         if len(new_indices) == 0:
#             continue

#         indices = np.concatenate(new_indices)

#     return X[indices], y[indices]

# import numpy as np




def balance_by_histogram(X, y, bins=50, target_count=500, min_threshold=50, max_iterations=10, dimensions=[0, 1, 2, 3, 4]):
    all_dims = [0, 1, 2, 3, 4] 
    indices = np.arange(len(y))
    
    # Pre-Computing: Kanten für alle Dimensionen berechnen
    bin_edges_all = []
    for dim in all_dims:
        _, edges = np.histogram(y[:, dim], bins=bins)
        bin_edges_all.append(edges)

    print(f"Start-Datenpunkte: {len(indices)}")
    print(f"Aktive Dimensionen (Filtern + Schutz): {dimensions}")

    for iteration in range(max_iterations):
        points_to_remove = set()
        
        # 1. Aktuellen Zustand berechnen
        current_counts = []
        point_bin_assignments = np.zeros((len(y), 5), dtype=int)
        
        for dim in all_dims:
            edges = bin_edges_all[dim]
            bin_idx = np.clip(np.digitize(y[:, dim], edges) - 1, 0, bins - 1)
            point_bin_assignments[:, dim] = bin_idx
            
            counts, _ = np.histogram(y[indices, dim], bins=edges)
            current_counts.append(counts)

        # 2. Finde zu volle Bins NUR in den aktiven Dimensionen
        candidates_to_remove = []
        for dim in dimensions: # <--- Ignorierte Dimensionen triggern hier keinen Filter
            counts = current_counts[dim]
            overpopulated_bins = np.where(counts > target_count)[0]
            
            for b in overpopulated_bins:
                in_bin_mask = (point_bin_assignments[:, dim] == b)
                available_in_bin = np.intersect1d(indices, np.where(in_bin_mask)[0])
                
                excess = len(available_in_bin) - target_count
                if excess > 0:
                    chosen_candidates = np.random.choice(available_in_bin, excess, replace=False)
                    candidates_to_remove.extend(chosen_candidates)
        
        candidates_to_remove = list(set(candidates_to_remove))
        np.random.shuffle(candidates_to_remove)

        # 3. DER GUARD-CHECK (Veto NUR noch für aktive Dimensionen!)
        for p in candidates_to_remove:
            p_bins = point_bin_assignments[p]
            
            veto = False
            for dim in dimensions: # <--- GEÄNDERT: Nur noch aktive Dimensionen dürfen ein Veto einlegen!
                assigned_bin = p_bins[dim]
                if current_counts[dim][assigned_bin] <= min_threshold:
                    veto = True
                    break
            
            if not veto:
                points_to_remove.add(p)
                # Wir tracken den Verlust trotzdem in allen Dimensionen, damit die Counts stimmen
                for dim in all_dims:
                    assigned_bin = p_bins[dim]
                    current_counts[dim][assigned_bin] -= 1

        if len(points_to_remove) == 0:
            print(f"Konvergenz erreicht in Iteration {iteration}. Keine weiteren sicheren Löschungen möglich.")
            break
            
        indices = np.array([i for i in indices if i not in points_to_remove])
        print(f"Iteration {iteration+1}: {len(points_to_remove)} Punkte gelöscht. Verbleibend: {len(indices)}")

    return X[indices], y[indices]


def get_random_config(param_ranges, log_scales):

    X_phys = np.zeros((N, input_size), dtype=np.float64)

    X_phys[:, 0] = get_random_arr(param_ranges["fshield_H2"], N, log_scales["fshield_H2"])
    X_phys[:, 1] = get_random_arr(param_ranges["fshield_CO"], N, log_scales["fshield_CO"])
    X_phys[:, 2] = get_random_arr(param_ranges["AV_mean"], N, log_scales["AV_mean"])
    X_phys[:, 3] = get_random_arr(param_ranges["chi_mean"], N, log_scales["chi_mean"])
    X_phys[:, 4] = get_random_arr(param_ranges["dens"], N, log_scales["dens"])
    X_phys[:, 5] = get_random_arr(param_ranges["eint"], N, log_scales["eint"])
    X_phys[:, 6] = get_random_arr(param_ranges["abh2"], N, log_scales["abh2"])
    X_phys[:, 7] = get_random_arr(param_ranges["abhp"], N, log_scales["abhp"])
    X_phys[:, 8] = get_random_arr(param_ranges["abco"], N, log_scales["abco"])
    X_phys[:, 9] = get_random_arr(param_ranges["time"], N, log_scales["time"])

    hydrogen_mask = (2 * X_phys[:, 6] + X_phys[:, 7] < 1 ) # filter for total hydrogen > 1
    print(N - hydrogen_mask.sum())

    return X_phys[hydrogen_mask]


def get_random_config_log_uniform(param_ranges, param_alphas):

    X_phys = np.zeros((N, input_size), dtype=np.float64)

    X_phys[:, 0] = get_random_arr_log_uniform(param_ranges["fshield_H2"], N, param_alphas["fshield_H2"])
    X_phys[:, 1] = get_random_arr_log_uniform(param_ranges["fshield_CO"], N, param_alphas["fshield_CO"])
    X_phys[:, 2] = get_random_arr_log_uniform(param_ranges["AV_mean"], N, param_alphas["AV_mean"])
    X_phys[:, 3] = get_random_arr_log_uniform(param_ranges["chi_mean"], N, param_alphas["chi_mean"])
    X_phys[:, 4] = get_random_arr_log_uniform(param_ranges["dens"], N, param_alphas["dens"])
    X_phys[:, 5] = get_random_arr_log_uniform(param_ranges["eint"], N, param_alphas["eint"])
    X_phys[:, 6] = get_random_arr_log_uniform(param_ranges["abh2"], N, param_alphas["abh2"])
    X_phys[:, 7] = get_random_arr_log_uniform(param_ranges["abhp"], N, param_alphas["abhp"])
    X_phys[:, 8] = get_random_arr_log_uniform(param_ranges["abco"], N, param_alphas["abco"])
    X_phys[:, 9] = get_random_arr_log_uniform(param_ranges["time"], N, param_alphas["time"])

    hydrogen_mask = (2 * X_phys[:, 6] + X_phys[:, 7] < 1 ) # filter for total hydrogen > 1
    print(N - hydrogen_mask.sum())

    return X_phys[hydrogen_mask]


def PychemDataGenerator(X_phys): # X_phys = (10 x N) array


    if __name__ == '__main__':

        try:
            set_start_method('fork')
        except RuntimeError:
            pass

        freeze_support()

        flchem_init(abundc, abundo, dust_to_gas_ratio,
                    cosmic_ray_ion_rate, NH_ext, Z_atom)

        fshield_H2_arr = X_phys[:,0]
        fshield_CO_arr = X_phys[:,1]
        AV_mean_arr = X_phys[:,2]
        chi_mean_arr = X_phys[:,3]
        dens_arr = X_phys[:,4]
        eint_arr = X_phys[:,5]
        abh2_arr = X_phys[:,6]
        abhp_arr = X_phys[:,7]
        abco_arr = X_phys[:,8]
        time_arr = X_phys[:,9]

        args = (time_arr, time_arr, G0, G0RP,
                fshield_H2_arr, fshield_CO_arr, AV_mean_arr, chi_mean_arr,
                dl, divv, dens_arr, eint_arr,
                abh2_arr, abhp_arr, abco_arr, abcp
                )

        print(args)

        results = flchem_evolve_mp(*args)

        print("RESULTS: ")
        print(results)

        # Filtering
        solver_error = results[-3] != 999

        num_zero = (results[1] == 1e-40) | (results[3] == 1e-40)
        low_abund = (results[1] < 1e-20) | (results[2] < 1e-20) | (results[3] < 1e-25)


        valid = (~solver_error) & (~low_abund) # & (~num_zero) # & (~low_abund)

        num_valid = np.sum(valid)
        print(f"Invalid points: {len(dens_arr) - num_valid}(solver error: {np.sum(solver_error)} num zeros: {np.sum(num_zero)}, abundances < 1e-20 : {np.sum(low_abund)})")

        
        X_out = np.zeros((num_valid, 10), dtype=np.float64)
        X_out = X_phys[valid]

        y_out = np.zeros((num_valid, 5), dtype=np.float64)
        # y[:,0] = np.maximum(results[1][valid], 1e-10) # abh2
        # y[:, 1] = np.maximum(results[2][valid], 1e-20) # abhp
        # y[:, 2] = np.maximum(results[3][valid], 1e-20) # abco
        # y[:, 3] = results[6][valid] # dust temp
        # y[:, 4] = results[7][valid] # internal energy

        print("T min: ", np.min(results[6][valid]))
        print("T max: ", np.max(results[6][valid]))

        y_out[:, 0] = results[1][valid] # abh2
        y_out[:, 1] = results[2][valid] # abhp
        y_out[:, 2] = results[3][valid] # abco
        y_out[:, 3] = np.maximum(np.minimum(results[6][valid], 40), 20) # dust temp
        y_out[:, 4] = results[7][valid] # internal energy

    # Combine all snapshots

    print(f"\nTotal training points: {X_out.shape[0]}")

    return X_out, y_out




### ------------- CREATE NEW TRAINING POINTS -------------- ###

X_phys = get_random_config(param_ranges, param_log_scales)
# X_phys = get_random_config(param_ranges, param_alphas)

X_physical, y_physical = PychemDataGenerator(X_phys)


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

# np.savez_compressed(
#     training_folder + "/raw_sample.npz",
#     X_physical=X_physical,
#     y_physical=y_physical
# )



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

### -------- SAVE POINTS WITH LOW HP ------------ ###

print("X_scaled shape: ", X_physical_log.shape)
X_physical_log_low_hp, y_physical_log_low_hp = X_physical_log[X_physical_log[:,1]<-10],  y_physical_log[y_physical_log[:,1]<-10]
print("X_scaled_low_hp shape: ", X_physical_log_low_hp.shape)

### -------------- BALANCE SAMPLE ------------------- ### 

# X_physical_balance_log, y_physical_balance_log = balance_by_histogram(X_physical_log, y_physical_log, bins=50, target_count=1000)
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


# datasets = [X_physical, X_physical_log, X_physical_balance_log, X_physical_balance]
# titles = ["X physical", "X physical log", "X physical balance log", "X balanced"]
# quantity = ["fshield_H2","fshield_CO","AV_mean", "chi_mean", "dens", "eint", "abh2","abhp","abco","time"]


# fig = plt.figure(figsize=(14, 12))
# outer = gridspec.GridSpec(2, 2, wspace=0.25, hspace=0.25)

# for idx, data in enumerate(datasets):
#     inner = gridspec.GridSpecFromSubplotSpec(
#         data.shape[1], 1,
#         subplot_spec=outer[idx],
#         hspace=0.1
#     )
    
#     for i in range(data.shape[1]):
#         ax = plt.Subplot(fig, inner[i])
#         ax.hist(data[:, i], bins=50, edgecolor='black')
#         ax.set_ylabel(quantity[i])
#         ax.grid(True, alpha=0.5)
#         ax.set_ylim(bottom=0)
#         ax.set_xticks([])
#         ax.set_xticklabels([])

        
#         fig.add_subplot(ax)

#     fig.axes[-data.shape[1]].set_title(titles[idx])

# plt.tight_layout()
# plt.show()

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

