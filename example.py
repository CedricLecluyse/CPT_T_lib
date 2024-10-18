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

# Optimized system:
df_Popt= pd.DataFrame(index=f)
Pmax= pd.DataFrame(index=f)
df_EffOpt= pd.DataFrame(index=f)
eff_max= pd.DataFrame(index=f)

df_Popt,df_EffOpt= T.CPT_optimize(Cmin=100, Cmax=1100, DC=100, Lmin=1, Lmax=21, DL=1, Vin=V_in, Cm=Cm, Load=Load, a=3, L=L, C=C)
df_Popt.columns.names= ['Capacitance', 'Inductance']
df_EffOpt.columns.names= ['Capacitance', 'Inductance']

Ltemp=[]
Ctemp=[]
df_P= df_Popt.max().unstack(level='Inductance')
df_Eff= df_EffOpt.max().unstack(level='Inductance')
for i in df_P.columns:
   Ltemp.append(int(i[2:]))
for j in df_P.index:
    Ctemp.append(int(j[2:]))

df_P.index= Ctemp
df_P.columns= Ltemp
df_Eff.index= Ctemp
df_Eff.columns= Ltemp

df_P= df_P.sort_index()
df_P= df_P.sort_index(axis=1)
df_Eff= df_Eff.sort_index()
df_Eff= df_Eff.sort_index(axis=1)

'''
Give equal weight to Pout and Efficiency. This is done by dividing Pout by the maximum value
in the data frame. The same goes for efficiency. The optimum point between Pout and Efficiency
can now be found by summing the two data frames and dividing by 2. This gives maximum power at maximum efficiency.
'''
df_Popt=df_Popt[df_Popt>100]
df_Popt.fillna(0, inplace=True)
wheight_Pout= df_Popt/(df_Popt.max().max())
wheight_eff= abs(df_EffOpt)/(abs(df_EffOpt).max().max())
#%%
opti= wheight_Pout+wheight_eff
opti_2= opti.loc[:2000000]

opti_max= opti.max()
opti_max.fillna(0, inplace=True)
opti_max= opti_max.unstack(level='Inductance')
print('Optimal combination: '+str(opti_max.stack().idxmax()) + 'No frequency restriction')
opti_max2= opti_2.max()
opti_max2.fillna(0, inplace=True)
opti_max2= opti_max2.unstack(level='Inductance')
print('Optimal combination: '+str(opti_max2.stack().idxmax()) + 'Max 2MHz')

