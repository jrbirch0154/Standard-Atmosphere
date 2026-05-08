# Standard Atmosphere Calculator

A Python implementation of the ISA (International Standard Atmosphere) model, covering altitudes from 0 to 86,000 m (sea level to the mesosphere).

## Usage

**Interactive:**
```
python atm.py
Enter altitude in m: 10000
AtmProps(T=223.15 K, P=26436.88 Pa, rho=0.4127 kg/m³)
```

**As a module:**
```python
from atm import atm

props = atm(10000)
print(props.T)      # Temperature, K
print(props.P)      # Pressure, Pa
print(props.rho)    # Density, kg/m³

# Bonus properties
print(props.T_C)        # Temperature, °C
print(props.rho_slug)   # Density, slug/ft³
```

## Atmosphere Layers

| Layer | Altitude (m) | Lapse Rate (K/m) |
|---|---|---|
| Troposphere | 0 – 11,000 | −0.0065 |
| Tropopause | 11,000 – 20,000 | 0 (isothermal) |
| Stratosphere 1 | 20,000 – 32,000 | +0.0010 |
| Stratosphere 2 | 32,000 – 47,000 | +0.0028 |
| Stratopause | 47,000 – 51,000 | 0 (isothermal) |
| Mesosphere 1 | 51,000 – 71,000 | −0.0028 |
| Mesosphere 2 | 71,000 – 86,000 | −0.0020 |

## Constants

| Symbol | Value | Description |
|---|---|---|
| T₀ | 288.15 K | Sea-level temperature |
| P₀ | 101,325 Pa | Sea-level pressure |
| g₀ | 9.80665 m/s² | Standard gravity |
| R | 287.058 J/(kg·K) | Specific gas constant, dry air |

## Requirements

- Python 3.7+
- Scipy
