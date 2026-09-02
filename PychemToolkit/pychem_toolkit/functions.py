import os
import h5py
import numpy as np

from flchem import flchem_evolve, flchem_evolve_mp, flchem_init
from timeit import default_timer as timer
from multiprocessing import Pool, freeze_support


### ---------- constants for pychem ------------- ###

kyr = 3.1556926e10

# totaltime, timestep = 3.1556926e10, 3.1556926e10
abundc, abundo = 1.4e-4, 3.2e-4
dust_to_gas_ratio, dust_temp = 1e0, 10.
G0, cosmic_ray_ion_rate = 1.7, 3e-17
NH_ext, Z_atom = 1e20, 1.
dl, divv = 1e15, 0.
abcp = np.nan
G0RP = 0.0

ch_muC = 12.011
ch_mf_scale = 1.0 + abundc * ch_muC

eps = 1e-38


def get_config_from_snapshot(snapshot_folder, N):
    snapshot_files = [f for f in os.listdir(snapshot_folder) if os.path.isfile(os.path.join(snapshot_folder, f))]
    print(snapshot_files)

    X_phys = []

    for snapshot_file in snapshot_files:
        print(f"\n Processing {snapshot_file}...")

        with h5py.File(snapshot_folder + snapshot_file, "r") as f:

            fshield_H2 = f["cdh2"][...]
            fshield_CO = f["cdco"][...]
            AV_mean = f["cdto"][...]
            chi_mean = f["chid"][...]
            dens = f["dens"][...]
            eint = f["eint"][...]
            abh2 = f["ih2 "][...]
            abhp = f["ihp "][...]
            abco = f["ico "][...]

            rng = np.random.random(size=dens.shape)
            mask = rng < N / dens.size

            # Apply mask
            fshield_H2_arr = fshield_H2[mask]
            fshield_CO_arr = fshield_CO[mask]
            AV_mean_arr = AV_mean[mask]
            chi_mean_arr = chi_mean[mask]
            dens_arr = dens[mask]
            eint_arr = eint[mask] * dens_arr
            abh2_arr = abh2[mask] * ch_mf_scale * 0.5
            abhp_arr = abhp[mask] * ch_mf_scale
            abco_arr = abco[mask] * ch_mf_scale / ch_muC

            time_arr = np.random.uniform(0.1 * kyr, 20.0 * kyr, size=dens_arr.shape)

            X_phys.append(np.column_stack([
                fshield_H2_arr, 
                fshield_CO_arr, 
                AV_mean_arr, 
                chi_mean_arr, 
                dens_arr, 
                eint_arr, 
                abh2_arr, 
                abhp_arr, 
                abco_arr, 
                time_arr
            ]))

        print(f"Sampled cells: {len(dens_arr)}")

    return np.vstack(X_phys)


def pychem_data_generator(X_physical, filter_data, abundance_thresholds, min_output_values, cut_1e40s): # X_phys = (10 x N) array

    # if __name__ == '__main__':

    freeze_support()

    flchem_init(abundc, abundo, dust_to_gas_ratio,
                cosmic_ray_ion_rate, NH_ext, Z_atom)

    fshield_H2_arr = X_physical[:,0]
    fshield_CO_arr = X_physical[:,1]
    AV_mean_arr = X_physical[:,2]
    chi_mean_arr = X_physical[:,3]
    dens_arr = X_physical[:,4]
    eint_arr = X_physical[:,5]
    abh2_arr = X_physical[:,6]
    abhp_arr = X_physical[:,7]
    abco_arr = X_physical[:,8]
    time_arr = X_physical[:,9]

    args = (time_arr, time_arr, G0, G0RP,
            fshield_H2_arr, fshield_CO_arr, AV_mean_arr, chi_mean_arr,
            dl, divv, dens_arr, eint_arr,
            abh2_arr, abhp_arr, abco_arr, abcp
            )


    results = flchem_evolve_mp(*args)

    # Filtering
    solver_error = results[-3] != 999

    num_zero = (results[1] == 1e-40) | (results[3] == 1e-40)
    low_abund = (results[1] < abundance_thresholds[0]) | (results[2] < abundance_thresholds[1]) | (results[3] < abundance_thresholds[2]) # ensure single precision


    valid = (~solver_error)

    if (filter_data == True):
        valid = valid & (~low_abund)
        print("abundances below threshold cutted out of sample")

    # if (filter_data == True):
    #     valid = (~solver_error) & (~low_abund) 
    # else: 
    #     valid = (~solver_error)

    if (cut_1e40s == True):
        valid = valid & (~num_zero)
        print("1e-40 abundances cutted out of sample")

    num_valid = np.sum(valid)

    print(f"Invalid points: {len(dens_arr) - num_valid}(solver error: {np.sum(solver_error)} num zeros: {np.sum(num_zero)}, abundances < threshold : {np.sum(low_abund)})")

    
    X_out = np.zeros((num_valid, 10), dtype=np.float64)
    X_out = X_physical[valid]

    y_out = np.zeros((num_valid, 5), dtype=np.float64)

    print("T min: ", np.min(results[6][valid]))
    print("T max: ", np.max(results[6][valid]))

    y_out[:, 0] = results[1][valid] # abh2
    y_out[:, 1] = results[2][valid] # abhp
    y_out[:, 2] = results[3][valid] # abco
    y_out[:, 3] = results[6][valid] # dust temp
    y_out[:, 4] = results[7][valid] # internal energy


    if not np.isnan(min_output_values[0]):
        y_out[:, 0] = np.maximum(y_out[:, 0], min_output_values[0])
        y_out[:, 1] = np.maximum(y_out[:, 1], min_output_values[1])
        y_out[:, 2] = np.maximum(y_out[:, 2], min_output_values[2])
        y_out[:, 3] = np.maximum(y_out[:, 3], min_output_values[3])
        y_out[:, 4] = np.maximum(y_out[:, 4], min_output_values[4])

    # Combine all snapshots

    print(f"\nTotal training points: {X_out.shape[0]}")

    print("TEST")

    return X_out, y_out


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



# def transform(X, scaling):
#     """
#     Apply feature transformations.

#     scaling:
#         "linear"
#         "log"
#         "asinh"
#         "power:<exponent>"  e.g. "power:0.2"
#     """
#     X = X.copy()

#     for i, mode in enumerate(scaling.values()):

#         if mode == "linear":
#             continue

#         elif mode == "log":
#             X[:, i] = np.log10(X[:, i] + eps)
        
#         elif mode == "asinh": 
#             X[:, i] = np.arcsinh(X[:, i])
#             print(f"Transformed column {i} using asinh scaling.")

#         elif mode.startswith("power:"):
#             exponent = float(mode.split(":")[1])
#             X[:, i] = X[:, i] ** exponent

#         else:
#             raise ValueError(f"Unknown scaling mode: {mode}")

#     return X

def transform(X, scaling):
    """
    Apply feature transformations.

    scaling:
        "linear"
        "log"
        "asinh:<scale>"
        "power:<exponent>"
    """
    X = X.copy()

    for i, mode in enumerate(scaling.values()):

        if mode == "linear":
            continue

        elif mode == "log":
            X[:, i] = np.log10(X[:, i] + eps)

        elif mode.startswith("asinh:"):
            scale = float(mode.split(":")[1])
            X[:, i] = np.arcsinh(X[:, i] / scale)

        elif mode.startswith("power:"):
            exponent = float(mode.split(":")[1])
            X[:, i] = X[:, i] ** exponent

        else:
            raise ValueError(f"Unknown scaling mode: {mode}")

    return X


# def inverse_transform(X, scaling):
#     """
#     Reverse feature transformations.
#     """
#     X = X.copy()

#     for i, mode in enumerate(scaling.values()):

#         if mode == "linear":
#             continue

#         elif mode == "log":
#             X[:, i] = 10 ** X[:, i] - eps

#         elif mode == "asinh": 
#             X[:, i] = np.sinh(X[:, i])

#         elif mode.startswith("power:"):
#             exponent = float(mode.split(":")[1])
#             X[:, i] = X[:, i] ** (1 / exponent)

#         else:
#             raise ValueError(f"Unknown scaling mode: {mode}")

#     return X


def inverse_transform(X, scaling):
    """
    Reverse feature transformations.

    scaling:
        "linear"
        "log"
        "asinh:<scale>"
        "power:<exponent>"
    """
    X = X.copy()

    for i, mode in enumerate(scaling.values()):

        if mode == "linear":
            continue

        elif mode == "log":
            X[:, i] = 10 ** X[:, i] - eps

        elif mode.startswith("asinh:"):
            scale = float(mode.split(":")[1])
            X[:, i] = scale * np.sinh(X[:, i])

        elif mode.startswith("power:"):
            exponent = float(mode.split(":")[1])
            X[:, i] = X[:, i] ** (1 / exponent)

        else:
            raise ValueError(f"Unknown scaling mode: {mode}")

    return X



def get_random_arr(value_range=(1, 10), N=1000, scaling="log"):
    """
    Draw N random samples according to the specified scaling.

    Parameters
    ----------
    value_range : tuple(float, float)
        (min, max) range.
    N : int
        Number of samples.
    scaling : {"linear", "log", "power:<exponent>"}
        Sampling strategy.
    """
    low, high = value_range

    if scaling == "linear":
        return np.random.uniform(low, high, N)

    elif scaling == "log":
        return 10 ** np.random.uniform(np.log10(low), np.log10(high), N)

    elif scaling.startswith("power:"):
        exponent = float(scaling.split(":")[1])

        low_t = low ** exponent
        high_t = high ** exponent

        return np.random.uniform(low_t, high_t, N) ** (1 / exponent)

    else:
        raise ValueError(
            f"Unknown scaling '{scaling}'. "
            "Use 'linear', 'log', or 'power:<exponent>'."
        )


def get_random_config(param_ranges, scalings, N):
    """
    Generate random physical configurations.

    Returns
    -------
    ndarray
        Shape (M, len(param_ranges)), where M <= N after filtering.
    """
    parameter_names = list(param_ranges.keys())

    X_phys = np.empty((N, len(parameter_names)), dtype=np.float64)

    for i, name in enumerate(parameter_names):
        X_phys[:, i] = get_random_arr(
            value_range=param_ranges[name],
            N=N,
            scaling=scalings[name],
        )

    # Keep only physically valid hydrogen abundances
    h2_idx = parameter_names.index("abh2")
    hp_idx = parameter_names.index("abhp")

    hydrogen_mask = (2 * X_phys[:, h2_idx] + X_phys[:, hp_idx]) < 1

    print(f"Rejected {N - hydrogen_mask.sum()} samples")

    return X_phys[hydrogen_mask]