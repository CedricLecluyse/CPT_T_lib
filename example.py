"""

@author: Cedric Lecluyse
"""

from tparameter_lib import Tparameter_lib as tpr
import numpy as np 
import pandas as pd


# =============================================================================
# Import all the components L / C / R / Cm needed to analyze the capacitive wireless power transfer system
# L= compensation inductance
# C= compensation capacitance
# R= load resistance
# Cm= mutual coupling capacitance (Total capacitance of the coupler)
# import L / C / R /Cm as a dataframe, for C and L: First column= primary value, Second column= loss resistance
# =============================================================================
f= np.arange(400e3, 10.05e6, 5e3)
L= pd.DataFrame({'L': [1e-6]*len(f), 'RL':[0]*len(f)}, index=f)
C= pd.DataFrame({'C': [1000e-12]*len(f), 'RC':[0]*len(f)}, index=f)
Load= 10
R= pd.DataFrame({'L': [0]*len(f), 'R':[Load]*len(f)}, index=f)
Cm= pd.DataFrame({'C': [100e-12]*len(f), 'RC':[0]*len(f)}, index=f)
V_in=48  # peak voltage

# =============================================================================
# Import library file
# =============================================================================
T= tpr(f)

# =============================================================================
# Simulate CPT system
# a : 0= primary L-compensation/ 1= Secondary L_compensation / 2= double sided L-compensation / 
#    3= double sided LC-compensation / 4= double sided LCL-compensation
# 
# Output values:
# 'S_in', 'P_in', 'V_in', 'I_in', 'Z_in', 'V_out', 'I_out', 'S_out',
#        'P_out', 'eff', 'T_out', 'T_in', 'T_imp', 'V1', 'I1', 'V2', 'I2', 'DV',
#        'mat_loss'
# 
#
# =============================================================================

df_sys= T.CPT_System(V_in,L1=L, L2=L, C1=C, C2=C, Cm=Cm, Load= R, a=3)
df_sys.P_out.plot()
