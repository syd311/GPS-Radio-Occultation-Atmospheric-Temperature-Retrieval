import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy import integrate
import os
# creating data 
def generate_standard_atmosphere_data():
    # Grid: Surface to 60km
    z = np.linspace(0, 60, 500) 
    Re = 6371.0 # Earth Radius
    g = 9.80665
    R = 287.05
    
    # A. Defining "True" Temperature (Standard Atmosphere)
    # Lapse rate 6.5 K/km up to 11km, then constant
    T_true = 288.15 - 6.5 * z
    T_true[z > 11] = 216.65 # Stratosphere (Isothermal)
    
    # B. Deriving "True" Pressure (Hydrostatic)
    # dP = -rho * g * dz = -(P/RT) * g * dz  => dP/P = -(g/RT) dz
    P_true = np.zeros_like(T_true)
    P_true[0] = 1013.25 # Surface hPa
    # Integrating upwards
    for i in range(len(z)-1):
        dz = (z[i+1] - z[i]) * 1000 # meters
        T_mean = 0.5 * (T_true[i] + T_true[i+1])
        # Hydrostatic step
        P_true[i+1] = P_true[i] * np.exp(- (g * dz) / (R * T_mean))
        
    # C. Deriving "True" Refractivity
    # N = 77.6 * P / T
    N_true = 77.6 * P_true / T_true
    
    # D. Forward Abel Transform: N -> Bending Angle
    # alpha(a) = -2a * integral( dN/dr / sqrt(r^2 - a^2) )
    print("   > Calculating synthetic bending angles...")
    bending_angle = np.zeros_like(z)
    
    # r = Re + z
    r = Re + z
    # dN/dr (gradient of refractivity)
    dN_dr = np.gradient(N_true, z) # N-units per km
    
    for i in range(len(z)-1):
        a = r[i]
        # Integrating from current layer upwards
        current_r = r[i:]
        current_dN = dN_dr[i:]
        
        # Kernel: 1 / sqrt(r^2 - a^2)
        # Nudging 'a' slightly to avoid div/0
        kernel = 1.0 / np.sqrt(current_r**2 - (a - 1e-5)**2)
        
        # Integrand: -2 * a * dN/dr * kernel * 1e-6 
        integrand = -2 * a * (current_dN * 1e-6) * kernel
        
        # Integrate
        bending_angle[i] = integrate.trapezoid(integrand, current_r)

    # E. Adding Noise
    np.random.seed(99) # Fixed seed for reproducibility
    noise = np.random.normal(0, 1.0e-6, size=len(z)) # 1 microradian noise
    bending_angle += noise
    
    # Saving the data to NetCDF
    ds = xr.Dataset(
        data_vars={
            'bending_angle': (['height'], bending_angle),
            'impact_param': (['height'], r),
            'altitude': (['height'], z),
            'temp_true': (['height'], T_true) # Saving truth for plotting
        },
        coords={'height': z}
    )
    os.makedirs('data', exist_ok=True)
    ds.to_netcdf('data/simulation_data.nc')
    return 'data/simulation_data.nc'

def abel_inverse(impact, alpha):
    n_levels = len(impact)
    ref = np.zeros(n_levels)
    # Abel Inversion
    for i in range(n_levels - 1):
        a = impact[i]
        # Integration upwards
        x = impact[i:]
        val_alpha = alpha[i:]
        # Nudging to handle singularity
        x[0] = a + 1e-5
        # Kernel for inversion
        kernel = 1.0 / np.sqrt(x**2 - a**2)
        integrand = val_alpha * kernel
        # Integrate
        integral = integrate.trapezoid(integrand, x)
        ref[i] = np.exp(integral / np.pi)
        
    return (ref - 1) * 1e6
# Retrieving Temperature from Refractivity
def retrieve_T(altitude, N):
    g = 9.80665
    R = 287.05
    k = 77.6
    
    # 1. Density
    N_safe = np.maximum(N, 0.1)
    rho = (N_safe * 100) / (k * R)
    
    # 2. Pressure (Top-Down)
    P = np.zeros_like(rho)
    P[-1] = 0.1 # Top boundary guess
    
    # Integrating backwards
    for i in range(len(rho)-2, -1, -1):
        dz = (altitude[i+1] - altitude[i]) * 1000.0
        rho_avg = 0.5 * (rho[i] + rho[i+1])
        dP = rho_avg * g * dz
        P[i] = P[i+1] + dP / 100.0 # to hPa
        
    # 3. Temp
    return k * P / N_safe

def run_full_simulation():
    # 1. Generating Data
    filepath = generate_standard_atmosphere_data()
    ds = xr.open_dataset(filepath)
    
    # 2. Extracting Variables
    z = ds['altitude'].values
    imp = ds['impact_param'].values
    b_ang = ds['bending_angle'].values
    t_true = ds['temp_true'].values
    
    # 3. Running Retrieval
    my_N = abel_inverse(imp, b_ang)

    my_T = retrieve_T(z, my_N)
    
    # 4. Plotting Results
    print("Step 4: Plotting Results...")
    plt.figure(figsize=(9, 11))
    
    # Truth (Red)
    plt.plot(t_true, z, 'r--', linewidth=2, label='True Atmosphere (Input)')
    
    # Retrieval (Blue)
    # Plotting from 2km to 40km to avoid surface noise and top boundary error
    mask = (z > 2) & (z < 40)
    plt.plot(my_T[mask], z[mask], 'b-', linewidth=2.5, label='Retrieved Profile (Output)')
    
    plt.title('End-to-End Simulation: GPS-RO Temperature Retrieval', fontsize=14)
    plt.xlabel('Temperature (K)', fontsize=12)
    plt.ylabel('Altitude (km)', fontsize=12)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(fontsize=12)
    plt.xlim(200, 300)
    plt.ylim(0, 45)
    
    plt.savefig('result.png')
if __name__ == "__main__":
    run_full_simulation()