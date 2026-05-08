# Standard Atmosphere
# Thu May  7
# Jacob Birch


# %% Initializing

import math
from dataclasses import dataclass

# CONSTANTS

T0 = 288.15  # Sea-level temperature, K
P0 = 101325.0  # Sea-level pressure, Pa
g0 = 9.80665  # Standard gravity, m/s²
R = 287.058  # Specific gas constant for dry air, J/(kg·K)

# Each row: (base altitude m, lapse rate K/m)
# Lapse rate = 0 → isothermal layer (exponential decay)
LAYERS = [
    (-0.0065),  # Troposphere
    (0.0000),  # Tropopause (isothermal)
    (0.0010),  # Stratosphere 1
    (0.0028),  # Stratosphere 2
    (0.0000),  # Stratopause (isothermal)
    (-0.0028),  # Mesosphere 1
    (-0.0020),  # Mesosphere 2
    (None),  # Upper boundary sentinel
]

BOUNDS = [0, 11000, 20000, 32000, 47000, 51000, 71000]


@dataclass
class AtmProps:
    T: float
    P: float
    rho: float

    @property  # K to C
    def T_C(self):
        return self.T - 273.15

    @property  # rho to slug
    def rho_slug(self):
        return self.rho / 515.379

    def __repr__(self):
        return (
            f"AtmProps(T={self.T:.2f} K, "
            f"P={self.P:.2f} Pa, "
            f"rho={self.rho:.4f} kg/m³)"
        )


# %% Standard ATM function


def atm(alt: float, alt0=0) -> AtmProps:

    # error checking -----
    if not isinstance(alt, (float, int)):
        raise TypeError(f"{alt} is not a number.")
    elif alt < 0 or alt > 86000:
        raise ValueError(f"Altitude {alt} m is outside valid range")
    # ----------

    if 0 <= alt <= BOUNDS[1]:
        T = T0 + LAYERS[0] * (alt - alt0)
        P = P0 * (T / T0) ** (g0 / (LAYERS[0] * R))
        rho = P / (R * T)

    elif BOUNDS[1] < alt <= BOUNDS[2]:
        T = T0 + LAYERS[0] * (11000 - alt0)
        P = P0 * (T / T0) ** (g0 / (LAYERS[0] * R))
        P = P * math.exp((g0 / (R * T)) * (alt - 11000))
        rho = P / (T * R)

    elif BOUNDS[2] < alt <= BOUNDS[3]:
        T = T0 + LAYERS[0] * (11000 - alt0)
        T1 = T
        P = P0 * (T / T0) ** (g0 / (LAYERS[0] * R))
        P = P * math.exp((g0 / (R * T)) * (20000 - 11000))
        T = T + LAYERS[2] * (alt - 20000)
        P = P * (T / T1) ** (g0 / (LAYERS[2] * R))
        rho = P / (R * T)

    elif BOUNDS[3] < alt <= BOUNDS[4]:
        T = T0 + LAYERS[0] * (11000 - alt0)
        T1 = T
        P = P0 * (T / T0) ** (g0 / (LAYERS[0] * R))
        P = P * math.exp((g0 / (R * T)) * (20000 - 11000))
        T = T + LAYERS[2] * (32000 - 20000)
        T2 = T
        P = P * (T / T1) ** (g0 / (LAYERS[2] * R))
        T = T + LAYERS[3] * (alt - 32000)
        P = P * (T / T2) ** (g0 / (LAYERS[3] * R))
        rho = P / (R * T)

    elif BOUNDS[4] < alt <= BOUNDS[5]:
        T = T0 + LAYERS[0] * (11000 - alt0)
        T1 = T
        P = P0 * (T / T0) ** (g0 / (LAYERS[0] * R))
        P = P * math.exp((g0 / (R * T)) * (20000 - 11000))
        T = T + LAYERS[2] * (32000 - 20000)
        P = P * (T / T1) ** (g0 / (LAYERS[2] * R))
        T2 = T
        T = T + LAYERS[3] * (47000 - 32000)
        P = P * (T / T2) ** (g0 / (LAYERS[3] * R))
        P = P * math.exp((g0 / (R * T)) * (alt - 47000))
        rho = P / (R * T)

    elif BOUNDS[5] < alt <= BOUNDS[6]:
        T = T0 + LAYERS[0] * (11000 - alt0)
        T1 = T
        P = P0 * (T / T0) ** (g0 / (LAYERS[0] * R))
        P = P * math.exp((g0 / (R * T)) * (20000 - 11000))
        T = T + LAYERS[2] * (32000 - 20000)
        P = P * (T / T1) ** (g0 / (LAYERS[2] * R))
        T2 = T
        T = T + LAYERS[3] * (47000 - 32000)
        T3 = T
        P = P * (T / T2) ** (g0 / (LAYERS[3] * R))
        P = P * math.exp((g0 / (R * T)) * (51000 - 47000))
        T = T + LAYERS[5] * (alt - 51000)
        P = P * (T / T3) ** (g0 / (LAYERS[5] * R))
        rho = P / (R * T)

    elif alt > BOUNDS[6]:
        T = T0 + LAYERS[0] * (11000 - alt0)
        T1 = T
        P = P0 * (T / T0) ** (g0 / (LAYERS[0] * R))
        P = P * math.exp((g0 / (R * T)) * (20000 - 11000))
        T = T + LAYERS[2] * (32000 - 20000)
        P = P * (T / T1) ** (g0 / (LAYERS[2] * R))
        T2 = T
        T = T + LAYERS[3] * (47000 - 32000)
        T3 = T
        P = P * (T / T2) ** (g0 / (LAYERS[3] * R))
        P = P * math.exp((g0 / (R * T)) * (51000 - 47000))
        T = T + LAYERS[5] * (71000 - 51000)
        T4 = T
        P = P * (T / T3) ** (g0 / (LAYERS[5] * R))
        T = T + LAYERS[6] * (alt - 71000)
        P = P * (T / T4) ** (g0 / (LAYERS[6] * R))
        rho = P / (R * T)

    else:
        raise ValueError(f"Altitude {alt} m is out of range.")

    return AtmProps(T, P, rho)


# %% Main

if __name__ == "__main__":

    try:
        atm_test = float(input("Enter altitude in m: "))
    except ValueError:
        print("Altitude must be a number.")
    else:
        try:
            print(atm(atm_test))
        except ValueError as e:
            print(e)
