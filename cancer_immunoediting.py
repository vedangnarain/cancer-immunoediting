"""
ENGR-E 599: Computational Bioengineering

Cancer Immunoediting: Simulating metabolic competition in the tumor microenvironment

@author: Vedang Narain

4 December 2018
"""

#==============================================================================
# NOTES
#==============================================================================

"""
— Acronyms: EV = Estimated Value
            DV = Default Value
            CTI = Cytotoxic Immune

— The model below is based on the model proposed by Kareva and Berezovskaya (2015).

  I. Kareva and F. Berezovskaya, "Cancer immunoediting: A process driven by
  metabolic competition as a predator–prey–shared resource type model", Journal
  of Theoretical Biology, vol. 380, pp. 463-472, 2015.
"""

#==============================================================================
# LIBRARIES
#==============================================================================

import matplotlib.pyplot as plt
import numpy as np
import tellurium as te

#==============================================================================
# FUNCTIONS
#==============================================================================

# defines function to plot year marker and legend (does not account for leap days)


def marker(years):  # accepts number of years (integer)
    plt.axvline(x=years*365, color='black', linestyle=':', alpha=0.5)
    plt.legend()

#==============================================================================
# MODEL
#==============================================================================

microenvironment = te.loada('''

// defines population dynamics as ODEs
Ta' = (ra*Ta*(G/(xT+bT*G))) - (ea*(Ta*I)/(I+s*(Ta+Tg))) - (ua*Ta)
Tg' = (rg*Tg*(G/(xT+bT*G))) - (eg*(Tg*I)/(I+s*(Ta+Tg))) - (ug*Tg)
G' = (G0-uG*G) - ((da*Ta+dg*Tg)*(G/(xT+bT*G))) - (dI*I*(G/(xI+bI*G)))
I' = (i0*(Ta+Tg)) - (uIC*I) + (rI*(I*kI/(I+kI))*(G/(xI+bI*G))*((I*(Ta+Tg))/(I+s*(Ta+Tg))))
C' = -uC*C
uIC' = (C*uIC) + ((uI-uIC)/p)

// defines chemotherapeutic intervention as events for different species
CI: at (time>4000): C = dose
CTa: at (time>4000): Ta = Ta * (1-(e/100))
CTg: at (time>4000): Tg = Tg * (1-(e/100))

// initializes species (all values are estimations)
T = 1                       // unit: vol            ,   EV: T(t) ≥ 0
Ta = 0.1 * T                // unit: vol            ,
Tg = T - Ta                 // unit: vol            ,   Note: T = Ta + Tg
G = 10                      // unit: moles          ,   EV: G(t) ≥ 0
I = 1                       // unit: vol            ,   EV: I(t) ≥ 0
C = 0                       // unit: vol
ra = 0.1                    // unit: 1/day          ,   DV: 0.431
rg = ra * k1 ; k1 = 1.9     // unit: 1/day          ,   EV: k1 ∈ [0.5, 5.0]
ua = 0.01                   // unit: 1/day
ug = ua * k2 ; k2 = 3       // unit: 1/day          ,   EV: k2 ∈ [1.0, 10.0]
ea = 0.1                    // unit: 1/day
eg = ea - k3 ; k3 = 0       // unit: 1/day          ,   EV: k3 ∈ [-0.1, 0.1]
s = 0.7                     // unit: n/a            ,   DV: 0.618
bT = 0.9                    // unit: n/a
bI = 0.9                    // unit: n/a
G0 = 1                      // unit: mol/day        ,   EV: G0 ∈ [0.01, 2.06]
uG = 0.01                   // unit: 1/day
da = 0.1                    // unit: mol/vol/day    ,   EV: da ∈ [0.048, 0.5]
dg = 2 * da                 // unit: mol/vol/day    ,   EV: dg ∈ [0.1268, 0.5262]
dI = dg                     // unit: mol/vol/day    ,   EV: dI ≈ dg
i0 = 0.0003                 // unit: vol/day        ,   DV: 0.0002, EV: i0 ∈ [0.001, 0.100]
uI = 0.01                   // unit: 1/day          ,   DV: 0.02
rI = 0.01                   // unit: 1/vol/day      ,   DV: 0.303
xT = 1                      // unit: moles
xI = 1                      // unit: moles
kI = 100                    // unit: vol
uC = 1                      // unit: 1/day
uIC = uI                    // unit: 1/day
p = 28                      // unit: days
e = 90                      // unit: %
dose = 1.0                  // unit: n/a
''')

#==============================================================================
# SIMULATION
#==============================================================================

# runs simulation
simulation_i0 = microenvironment.simulate(0, 15000, 100000)

# plots results of simulation
plt.figure(1)
plt.subplot(3, 1, 1)
plt.plot(simulation_i0[:, 0], simulation_i0[:, 1]+simulation_i0[:, 4], label='T')
plt.plot(simulation_i0[:, 0], simulation_i0[:, 4], label='T$_g$', linestyle='--')
plt.plot(simulation_i0[:, 0], simulation_i0[:, 1], label='T$_a$', linestyle='--')
plt.title('i$_0$ = %1.5f mm$^3$/day | dose = %1.1f units'
          % (microenvironment.i0, microenvironment.dose))
plt.ylabel('volume of tumor cells (mm$^3$)')
marker(30)
plt.subplot(3, 1, 2)
plt.plot(simulation_i0[:, 0], simulation_i0[:, 2], label='G')
plt.ylabel('amount of glucose (mol)')
marker(30)
plt.subplot(3, 1, 3)
plt.plot(simulation_i0[:, 0], simulation_i0[:, 3], label='I')
plt.xlabel('time (days)')
plt.ylabel('volume of CTI cells (mm$^3$)')
marker(30)
plt.subplots_adjust(top=2)
plt.show()

#==============================================================================
# PARAMETER SCAN
#==============================================================================

# acquires time array
parascan = simulation_i0[:, 0]

# initializes label list
label_list = []

# plots parameter scan
for dose in np.arange(0, 2.1, 0.2):
    microenvironment.reset()
    microenvironment.dose = dose
    label_list.append('dose = '+str(round(dose, 3)))  # adjusts decimal places
    simulation_T = microenvironment.simulate(0, 15000, 100000)
    Tsum = simulation_T[:, 1] + simulation_T[:, 4]
    parascan = np.vstack([parascan, Tsum])
parascan = np.transpose(parascan)
te.plotArray(parascan, labels=label_list, xlabel='time (days)',
             ylabel='volume of tumor cells (mm$^3$)', title='i$_0$ = %1.5f mm$^3$/day'
             % (microenvironment.i0))
