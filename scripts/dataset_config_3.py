kyr = 3.1556926e10
abundc = 1.4e-4

input_range = {
    "fshield_H2": (1e-10, 1),
    "fshield_CO": (1e-10, 1),
    "AV_mean": (1e-2, 1e3),
    "chi_mean": (1e-5, 1),
    "dens": (1e-26, 1e-16),
    "eint": (1e-15, 1e-5),
    "abh2": (1e-10, 0.5),
    "abhp": (1e-15, 1),
    "abco": (1e-21, abundc),
    "time": (0.1 * kyr, 20 * kyr),
}

input_scaling = {
    "fshield_H2": "log",
    "fshield_CO": "log",
    "AV_mean": "log",
    "chi_mean": "linear",
    "dens": "log",
    "eint": "log",
    "abh2": "log",
    "abhp": "log",
    "abco": "log",
    "time": "linear",
}

output_scaling = {
    "abh2": "log",
    "abhp": "log",
    "abco": "log",
    "tdus": "linear",
    "eint": "log",
}