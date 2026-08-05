import numpy as np
import matplotlib.pyplot as plt

training_folder = "data_hist_balanced_10D_nolims_2"

data = np.load(training_folder + "/training_data_physical.npz")
training_data = data["data"]

data = np.load(training_folder + "/test_data_physical.npz")
test_data = data["data"]


total_sample = np.concatenate( (training_data, test_data), axis=0 )


X_physical = total_sample[:,1:11]
y_physical = total_sample[:,11:]


print("input shape = " , X_physical.shape)
print("output = " , y_physical.shape)


# ----------- check min max from sampel --------------- #


quantity = ["fshield_H2","fshield_CO","AV_mean", "chi_mean", "dens", "eint", "abh2", "abhp", "abco", "time"]
data = X_physical 

print(f"{'Quantity':<10} {'Min':>15} {'Max':>15}")
print("-" * 42)

for i, name in enumerate(quantity):
    col_min = np.min(data[:, i])
    col_max = np.max(data[:, i])
    print(f"{name:<10} {col_min:15.6e} {col_max:15.6e}")


quantity = ["h2", "hp", "co", "T_dust", "eint"]
data = y_physical 

print(f"{'Quantity':<10} {'Min':>15} {'Max':>15}")
print("-" * 42)

for i, name in enumerate(quantity):
    col_min = np.min(data[:, i])
    col_max = np.max(data[:, i])
    print(f"{name:<10} {col_min:15.6e} {col_max:15.6e}")


# mask_h2 = np.isclose(y_physical[:,0], 1e-40, rtol=0, atol=1e-45)

# print(y_physical[:,0])

# print(np.sum(mask_h2))


# plt.hist(np.log10(X_physical[:,6]), bins= 50)
# # plt.hist(np.log10(X_physical[mask_h2,a]), bins= 50)
# plt.show()

# plt.hist(np.log10(y_physical[:,0]), bins= 50)
# plt.show()


# ------------ plot sample distribution ---------------- #

title = "X physical"
data = X_physical
quantity = ["fshield_H2","fshield_CO","AV_mean", "chi_mean", "dens", "eint", "abh2","abhp","abco","time"]


fig, axes = plt.subplots(
    nrows=data.shape[1],
    ncols=1,
    figsize=(8, 10),
    sharex=False
)

# Ensure axes is iterable if there is only one column
if data.shape[1] == 1:
    axes = [axes]

for i, ax in enumerate(axes):
    ax.hist(data[:, i], bins=50, edgecolor="black")
    ax.set_ylabel(quantity[i])
    ax.grid(True, alpha=0.5)
    ax.set_ylim(bottom=0)

    if i < data.shape[1] - 1:
        ax.set_xticklabels([])

axes[0].set_title(title)

plt.tight_layout()
plt.savefig(training_folder + "/output_sampling")
plt.show()


title = "y physical"
data = y_physical
quantity = ["h2", "hp", "co", "T_dust", "eint"]


fig, axes = plt.subplots(
    nrows=data.shape[1],
    ncols=1,
    figsize=(8, 10),
    sharex=False
)

# Ensure axes is iterable if there is only one column
if data.shape[1] == 1:
    axes = [axes]

for i, ax in enumerate(axes):
    ax.hist(data[:, i], bins=50, edgecolor="black")
    ax.set_ylabel(quantity[i])
    ax.grid(True, alpha=0.5)
    ax.set_ylim(bottom=0)

    if i < data.shape[1] - 1:
        ax.set_xticklabels([])

axes[0].set_title(title)

plt.tight_layout()
plt.savefig(training_folder + "/output_sampling")
plt.show()