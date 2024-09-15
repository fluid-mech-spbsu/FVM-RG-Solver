import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from tqdm import tqdm

import warnings
warnings.filterwarnings("ignore")

# Constants
R = 8.3144598
kB = 1.38064852e-23
N_A = 6.02214129e23
hc = 6.62559e-34 * 2.99792458e8

# Gas data 
gases = {
    "argon": {
        "molecule": False,
        "gamma": 1.667,
        "molarMass": 0.039948, # kg/mol
        "mass": 6.633521356992e-26, # kg
        "speed_of_sound_room": 322.6 , # m/s
    },
    "methane": {
        "molecule": True,
        "gamma": 1.3084,
        "molarMass": 0.016043, # kg/mol
        "mass": 2.663732314e-26,    # kg
        "speed_of_sound_room": 450.06, # m/s
        "number_of_atoms": 4,
        "rotational_levels": 3,
        "spectral_data": np.array([302550, 158270, 315680, 136740]), # m^-1
        "degeneracy": np.array([1, 2, 3, 3]),
        "D_diss": 3668582.3189, # m^-1, converted from 438.86 kJ/mol https://www.weizmann.ac.il/oc/martin/tools/hartree.html
        "max vibrational levels": [9, 17, 9, 20] # harmonic oscillator 
    }
}

gas = "methane"

molecule = gases[gas]["molecule"]
gamma = gases[gas]["gamma"]
molarMass = gases[gas]["molarMass"]
mass = gases[gas]["mass"]
speed_of_sound_room = gases[gas]["speed_of_sound_room"]

ground_st_energy = 0
possible_inds = []
ground_state_only = False

# Calculate vibrational energy levels
if molecule:
    om_e = gases[gas]["spectral_data"]
    ds = gases[gas]["degeneracy"]
    D_diss = gases[gas]["D_diss"]
    max_vibr_lvls = gases[gas]["max vibrational levels"]
    rot_lvls = gases[gas]["rotational_levels"]

    es = hc * om_e
    gr_st_energy = sum(hc * (om_e * ds) / 2)  # Ground state energy

    if ground_state_only:
        possible_inds.append([0] * gases[gas]["number_of_atoms"])
    else:
        indices = [0] * gases[gas]["number_of_atoms"]

        for i1 in range(max_vibr_lvls[0]):
            indices[0] = i1
            for i2 in range(max_vibr_lvls[1]):
                indices[1] = i2
                for i3 in range(max_vibr_lvls[2]):
                    indices[2] = i3
                    for i4 in range(max_vibr_lvls[3]):
                        indices[3] = i4

                        e = (
                            om_e[0] * (i1 + ds[0] / 2.) +
                            om_e[1] * (i2 + ds[1] / 2.) +
                            om_e[2] * (i3 + ds[2] / 2.) +
                            om_e[3] * (i4 + ds[3] / 2.)
                        )

                        # If vibrational energy is within the dissociation energy, add the combination
                        if e <= D_diss:
                            possible_inds.append(indices.copy())


# Initial data
Ma = 3.8
T_left = 300 # K
pressure = 100 # Pa

# Calculated parameters
velocity_left = Ma * speed_of_sound_room
density_left = pressure * molarMass / (R * T_left)

################################################################################################

print("Initial conditions:")
print("v_0 = ", velocity_left)
print("rho_0 = ", density_left)
print("T_0 = ", T_left)
print("p_0 = ", pressure)
print("______________________________________")

################################################################################################

def solver_approx(velocity_left, density_left, T_left):
    # Rankine-Hugoniot boundary conditions
    density_right = ((gamma + 1) * pow(Ma,2)) / (2 + (gamma - 1) * pow(Ma,2)) * density_left
    velocity_right = velocity_left * density_left / density_right 

    pressure_left = R * T_left * density_left / molarMass 
    pressure_right = (pow(Ma,2) * 2 *gamma - (gamma - 1)) / (gamma + 1) * pressure_left
    T_right = pressure_right / (density_right * R / molarMass)

    return velocity_right, density_right, T_right

print("Getting answer via approximate solver:")
ans2 = solver_approx(velocity_left, density_left, T_left)
print("v_n = ", ans2[0])
print("rho_n = ", ans2[1])
print("T_n = ", ans2[2])
print("p_n = ", ans2[1] * R * ans2[2] / molarMass)
print("______________________________________")


################################################################################################


def solver(velocity_left, density_left, T_left):
    
    v_0, rho_0, T_0 = velocity_left, density_left, T_left
    
    def func(x):

        if molecule:
            Zvibr_0, Zvibr_1 = 0, 0
            for inds in possible_inds:
                s = ((inds[1]+1)*(inds[2]+1)*(inds[2]+2)*(inds[3]+1)*(inds[3]+2)/4)
                e_0 = sum(es*inds)
                Zvibr_0 += s*np.exp(-e_0/(kB*T_0))
                Zvibr_1 += s*np.exp(-e_0/(kB*x[2]))
            
            Uvibr_0, Uvibr_1 = 0, 0
            for inds in possible_inds:
                s = ((inds[1]+1)*(inds[2]+1)*(inds[2]+2)*(inds[3]+1)*(inds[3]+2)/4)
                e_0 = sum(es*inds)
                Uvibr_0 += s * e_0 * np.exp(-e_0/(kB*T_0))
                Uvibr_1 += s * e_0 * np.exp(-e_0/(kB*x[2]))
            
            Uvibr_0 = Uvibr_0/Zvibr_0 + gr_st_energy
            # Uvibr_0 = 0 # only ground state case
            Uvibr_1 = Uvibr_1/Zvibr_1 + gr_st_energy
            # Uvibr_1 = 0 # only ground state case

            return [
            x[0] * x[1] - v_0 * rho_0, # x[0] - velocity, x[1] - density, x[2] - temperature
            x[1] * np.power(x[0],2) + R * x[2] * x[1] / molarMass 
            - rho_0 * np.power(v_0,2) - R * T_0 * rho_0 / molarMass,
            x[1] * x[0] * (3/2 * R * x[2]  / molarMass +  rot_lvls / 2 * kB * x[2] / mass 
                           + Uvibr_1 / mass + np.power(x[0],2) / 2 + R * x[2] * x[1] / (molarMass * x[1]))
            - rho_0 * v_0 * ( 3/2 * R * T_0 / molarMass + rot_lvls / 2 * kB * T_0 / mass 
                             + Uvibr_0 / mass + np.power(v_0,2) / 2 + R * T_0 * rho_0 / (molarMass * rho_0))
            ]
        else: # monatomic gas
            return [
            x[0] * x[1] - v_0 * rho_0, # x[0] - velocity, x[1] - density, x[2] - temperature
            x[1] * np.power(x[0],2) + R*x[2]*x[1]/molarMass - rho_0 * np.power(v_0,2) - R*T_0*rho_0/molarMass,
            5 * (rho_0 / x[1] - 1 / 5) * R * x[1] * x[2] / molarMass - 4 * R * T_0 * rho_0 / molarMass
            # x[1] * x[0] * (3/2*R*x[2]/molarMass + np.power(x[0],2)/2 + R*x[2]*x[1]/(molarMass*x[1]))
            # - rho_0 * v_0 * (3/2*R*T_0/molarMass + np.power(v_0,2)/2 + R*T_0*rho_0/(molarMass*rho_0))
            ]

        
    
    vs, rhos, Ts = [], [], []
    
    for x in tqdm(range(1, 200)):
        cur_ans = fsolve(func, [x, x/100, x])
        if any(np.isclose(func(cur_ans), [0.0,0.0,0.0])):
            vs.append(cur_ans[0])
            rhos.append(cur_ans[1])
            Ts.append(cur_ans[2])
    ans = [np.median(vs), np.median(rhos), np.median(Ts)]
    
    print("Must be near zero values =", func(ans))

    return ans

print("Getting answer via numeric scipy.fsolver:")
ans = solver(velocity_left, density_left, T_left)
print("v_n = ", ans[0])
print("rho_n = ", ans[1])
print("T_n = ", ans[2])
print("p_n = ", ans[1]*R*ans[2]/molarMass)
print("______________________________________")