# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 16:48:15 2023

@author: Cedric Lecluyse
Library file based on unpulished IEEE ECCE paper
"""
import numpy as np
import pandas as pd

class Tparameter_lib:
    
    def __init__(self, f):
        self.f= f
        
    def coup_cap(self, c):
        C= pd.DataFrame(index=c.index)
        RC= pd.DataFrame(index=c.index)        
        C['C']= c['C']
        RC['RC']=c['RC']
        
        
        '''
        Coupler transformation matrix T_m
        --       --
        | 1   Z_m |
        | 0    1  |
        --       --
        '''
        
        Z_m= [0]*len(self.f) #-1j/(2*np.pi*f*C) 
        T_m= [0]*len(self.f)

        for i in range(len(self.f)):

            Z_m[i]= (RC[RC.columns[-1]].iloc[RC.index.get_loc(self.f[i],method='nearest')])-1j/(2*np.pi*self.f[i]*C[C.columns[-1]].iloc[C.index.get_loc(self.f[i],method='nearest')])
            T_m[i]= [[1, Z_m[i]],[0,1]]

        return T_m
    
    

    def comp_L(self,l):
        L= pd.DataFrame(index=l.index)
        RL= pd.DataFrame(index=l.index)
        L['L']=l['L']
        RL['RL']=l['RL']
        
        '''
        L-compensation transformation matrix T_L
         --       --
         | 1   Z_L |
         | 0    1  |
         --       --
        '''
        
        Z_L= [0]*len(self.f)
        T_L= [0]*len(self.f)
        for i in range(len(self.f)):
            Z_L[i]= RL[RL.columns[0]].iloc[RL.index.get_loc(self.f[i],method='nearest')]+1j*(2*np.pi*self.f[i]*L[L.columns[0]].iloc[L.index.get_loc(self.f[i],method='nearest')])
            T_L[i]= np.array([[1, Z_L[i]],[0,1]])
        
        return T_L
    

    def comp_C(self, c):
        C= pd.DataFrame(index=c.index)
        RC= pd.DataFrame(index=c.index)
        C['C']= c['C']
        RC['RC']=c['RC']
        '''
        C-compensation transformation matrix T_L
         --          --
         | 1        0 |
         | 1/Z_C    1 |
         --          --
        '''
        
        Z_RC= [0]*len(self.f)
        T_RC= [0]*len(self.f)
        
        for i in range(len(self.f)):
            Z_RC[i]= RC[RC.columns[0]].iloc[RC.index.get_loc(self.f[i],method='nearest')] - 1j/(2*np.pi*self.f[i]*C[C.columns[0]].iloc[C.index.get_loc(self.f[i],method='nearest')])
            T_RC[i]= np.array([[1, 0] , [1/Z_RC[i],1]])
            
        return T_RC  
    
    def input_Z(self,T_L1=None, T_L2=None, T_L1_2=None, T_L2_2=None, T_C1=None, T_C2=None, T_C1_2=None, T_C2_2=None, T_m=None, Load=None, a=None):
        '''
            Parameters
            ----------
            T_L1 : T-matrix Primary L
            T_L2 : T-matrix Secondary L
            T_L1_2:T-matrix Primary L #2 (LCL or LCLC)
            T_L2_2: T-matrix Secondary L #2 (LCL or LCLC)
            T_C1 : T-matrix Primary C
            T_C2 : T-matrix Secondary C
            T_C1_2 : T-matrix Primary C #2 (LCLC)
            T_C2_2 : T-matrix Secondary C #2 (LCLC)
            T_m : T-matrix coupling capacitor
            Load: Impedance of the load
            a : 0= primary L-compensation/ 1= Secondary L_compensation / 2= double sided L-compensation / 3= double sided LC-compensation / 4= double sided LCL-compensation
    
            Returns
            -------
            - T-matrix of the black box CPT impedance
            - Input impedance of the complete system, including load
    
        '''
        T_imp=[0]*len(self.f)
        Z_in=[0]*len(self.f)
        
        for i in range(len(self.f)):
            R_L= Load[Load.columns[1]].iloc[Load.index.get_loc(self.f[i],method='nearest')]+1j*(2*np.pi*self.f[i]*Load[Load.columns[0]].iloc[Load.index.get_loc(self.f[i],method='nearest')])
            if a==0:
                T_imp[i]= T_L1[i]@T_m[i] #single sided L-compensation Primary
            elif a==1: 
                T_imp[i]= T_m[i]@T_L2[i] #single sided L-compensation Secondary
            elif a==2:
                T_imp[i]= T_L1[i]@T_m[i]@T_L2[i] #double sided L compensation
            elif a==3:
                T_imp[i]= T_L1[i]@T_C1[i]@T_m[i]@T_C2[i]@T_L2[i]     #double sided LC compensation
            elif a==4:
                T_imp[i]= T_L1[i]@T_C1[i]@T_L1_2[i]@T_m[i]@T_L2[i]@T_C2[i]@T_L2_2[i] #double sides LCL-compensation   
            else:
                  print(" Give a number: 0 - 4") 
                    
           
            Z_in[i]=  (T_imp[i][0][0]*R_L+T_imp[i][0][1])/(T_imp[i][1][0]*R_L+T_imp[i][1][1])
            # Z_in[i]=  (T_imp[i][0][0]*Load-T_imp[i][0][1])/(T_imp[i][1][0]*Load-T_imp[i][1][1])
            
        return T_imp, Z_in
    
    
    def output_matrix(self,T_imp, V, Z_in, Pin):
        '''
        
        Feel free to add cases e.g. LC-LCL / Change else statement!

        Parameters
        ----------
        T_imp : T-matrix of the black box CPT impedance (see input_Z)
        V_in : RMS value of input voltage
        Z_in : Input impedance of the complete system, including load (see input_Z)

        Returns
        -------
        T_out : [V_out I_out]
        S_in : Schijnbaar vermogen
        T_in : [V_in I_in]
        '''
        V_in= [V]*len(self.f)
        P_in=[Pin]*len(self.f)
        I_in=[0]*len(self.f)
        S_in=[0]*len(self.f)
        T_in=[0]*len(self.f)
        T_out=[0]*len(self.f)
        
        
        
        
        
        
        for i in range(len(self.f)):
            
            I_in[i]=V_in[i]/Z_in[i]
            # if abs(I_in[i]) > abs(P_in[i]/V_in[i]):
            #     I_in[i]= P_in[i]/V_in[i]
                
            S_in[i]= V_in[i]*I_in[i]
            T_in[i]=[[V_in[i]], [I_in[i]]] 
            
            try:
                T_out[i]=np.linalg.inv(T_imp[i])@T_in[i]
                
            except:
                T_out[i]=[[0,0],[0,0]]
                print("error")
                return T_out, S_in, T_in
        
        return T_out ,  S_in , T_in
        
    
    def V1_Voltage(self,T_in=None, T_L1=None, T_L2=None, T_L1_2=None, T_L2_2=None, T_C1=None, T_C2=None, T_C1_2=None, T_C2_2=None, T_m=None, T_Load=None, a=None):
        '''
        
        Feel free to add cases e.g. LC-LCL / Change else statement! 
        
        Parameters
        ----------
        V_out : Output voltage CPT system (see output matrix)
        I_out : Output current CPT system (see output matrix)
        T_L1 : T-matrix Primary L
        T_L2 : T-matrix Secondary L
        T_L1_2:T-matrix Primary L #2 (LCL or LCLC)
        T_L2_2: T-matrix Secondary L #2 (LCL or LCLC)
        T_C1 : T-matrix Primary C
        T_C2 : T-matrix Secondary C
        T_C1_2 : T-matrix Primary C #2 (LCLC)
        T_C2_2 : T-matrix Secondary C #2 (LCLC)
        T_m : T-matrix coupling capacitor
        Load: Impedance of the load
        a : 0= primary L-compensation/ 1= Secondary L_compensation / 2= double sided L-compensation / 3= double sided LC-compensation / 4= double sided LCL-compensation
    
        Returns
        -------
        - Returns voltage V1 at primary plates
        - Returns current I1 at primary plates        
        '''
        # df_Zcpt= pd.DataFrame(index=self.f)
        # df_Tcpt= pd.DataFrame(index=self.f)
        # df_out= pd.DataFrame(index=self.f)
        
        # df_Tcpt['Tcpt'], df_Zcpt['Zcpt']= self.CPT_System(L1= L1, L2=L2, L1_2= L1_2, L2_2= L2_2, C1= C1,C2=C2, Cm=Cm, Load=Load, a=a)
        
        # T_out, S_in, T_in= self.output_matrix(df_Tcpt[df_Tcpt.columns[0]], V_in, df_Zcpt[ df_Zcpt.columns[0]])
        
        # df_out['V_out']=[item[0][0] for item in T_out]
        # df_out['I_out']=[item[1][0] for item in T_out]
      
        
        V_1= [0]*len(self.f)
        I_1= [0]*len(self.f)
        T_imp=[0]*len(self.f)
        T_1=[0]*len(self.f)


        # T_m= self.coup_cap(Cm)
        
        if a== 0:
            # T_L1= self.comp_L(L1)
            for i in range(len(self.f)):
                T_imp[i]= T_L1[i]
                
        elif a== 1:
            for i in range(len(self.f)):
                T_imp[i]= np.identity(2)
                
        elif a==2:
            # T_L1= self.comp_L(L1)
            for i in range(len(self.f)):
                T_imp[i]= T_L1[i]
                
        elif a==3:
            # T_L1= self.comp_L(L1)
            # T_C1= self.comp_C(C1)
            for i in range(len(self.f)):
                T_imp[i]=T_L1[i]@T_C1[i]
            
        elif a==4:
            # T_L1= self.comp_L(L1)
            # T_L1_2=self.comp_L(L1_2)
            # T_C1= self.comp_C(C1)
            for i in range(len(self.f)):
               T_imp[i]=T_L1[i]@T_C1[i]@T_L1_2[i]
                    
                
            
        else:
            print(" Give a number: 0 - 4") 

        
        for i in range(len(self.f)):    
            T_1[i]= np.linalg.inv(T_imp[i])@T_in[i]
            
        return T_1
    
    def V2_Voltage(self,T_in=None, T_L1=None, T_L2=None, T_L1_2=None, T_L2_2=None, T_C1=None, T_C2=None, T_C1_2=None, T_C2_2=None, T_m=None, T_Load=None, a=None):
        '''
        
        Feel free to add cases e.g. LC-LCL / Change else statement! 
        
        Parameters
        ----------
        V_out : Output voltage CPT system (see output matrix)
        I_out : Output current CPT system (see output matrix)
        T_L1 : T-matrix Primary L
        T_L2 : T-matrix Secondary L
        T_L1_2:T-matrix Primary L #2 (LCL or LCLC)
        T_L2_2: T-matrix Secondary L #2 (LCL or LCLC)
        T_C1 : T-matrix Primary C
        T_C2 : T-matrix Secondary C
        T_C1_2 : T-matrix Primary C #2 (LCLC)
        T_C2_2 : T-matrix Secondary C #2 (LCLC)
        T_m : T-matrix coupling capacitor
        Load: Impedance of the load
        a : 0= primary L-compensation/ 1= Secondary L_compensation / 2= double sided L-compensation / 3= double sided LC-compensation / 4= double sided LCL-compensation
    
        Returns
        -------
        - Returns voltage V2 at Secondary plates
        - Returns current I2 at Secondary plates        
        '''
    
        V_2= [0]*len(self.f)
        I_2= [0]*len(self.f)
        T_imp=[0]*len(self.f)
        T_2=[0]*len(self.f)
        
        # T_m= self.coup_cap(Cm)
        
        if a== 0:
            # T_L1= self.comp_L(L1)
            for i in range(len(self.f)):
                T_imp[i]= T_L1[i]@T_m[i]
                
        elif a== 1:
            for i in range(len(self.f)):
                T_imp[i]= T_m[i]
                
        elif a==2:
            # T_L1= self.comp_L(L1)
            for i in range(len(self.f)):
                T_imp[i]= T_L1[i]@T_m[i]
                
        elif a==3:
            # T_L1= self.comp_L(L1)
            # T_C1= self.comp_C(C1)
            for i in range(len(self.f)):
                T_imp[i]=T_L1[i]@T_C1[i]@T_m[i]
        
        elif a==4:
            # T_L1= self.comp_L(L1)
            # T_L1_2=self.comp_L(L1_2)
            # T_C1= self.comp_C(C1)
            for i in range(len(self.f)):
               T_imp[i]=T_L1[i]@T_C1[i]@T_L1_2[i]@T_m[i]
                
                
        else:
            print(" Give a number: 0 - 4") 

        
        for i in range(len(self.f)):    
            T_2[i]= np.linalg.inv(T_imp[i])@T_in[i]
            
        return T_2
    
    def CPT_System(self,V_in=None, L1=None, L2=None, L1_2=None, L2_2=None, C1=None, C2=None, C1_2=None, C2_2=None, Cm=None, Load=None, P_in=None, a=None, opt=False):
        '''
        a : 0= primary L-compensation/ 1= Secondary L_compensation / 2= double sided L-compensation / 3= double sided LC-compensation / 4= double sided LCL-compensation
    
        '''        
        df_output=pd.DataFrame(index=self.f)
        
        
        T_Cm= self.coup_cap(Cm)
        if a== 0:
            T_L1= self.comp_L(L1)
            T_L2=None
            T_C1=None
            T_C2=None
            T_L1_2= None
            T_L2_2= None
            T_imp, Z_in= self.input_Z(T_L1=T_L1,T_m= T_Cm, Load=Load, a=a)
        elif a== 1:
            T_L1=None
            T_L2= self.comp_L(L2)
            T_C1=None
            T_C2=None
            T_L1_2= None
            T_L2_2= None
            T_imp, Z_in= self.input_Z(T_L2=T_L2, T_m= T_Cm, Load=Load, a=a)
        elif a==2:
            T_L1= self.comp_L(L1)
            T_L2= self.comp_L(L2)
            T_C1=None
            T_C2=None
            T_L1_2= None
            T_L2_2= None
            T_imp, Z_in= self.input_Z(T_L1=T_L1, T_L2=T_L2, T_m= T_Cm, Load=Load, a=a)
        elif a==3:
            T_L1= self.comp_L(L1)
            T_L2= self.comp_L(L2)
            T_C1= self.comp_C(C1)
            T_C2= self.comp_C(C2)
            T_L1_2= None
            T_L2_2= None
            T_imp, Z_in= self.input_Z(T_L1= T_L1, T_L2=T_L2, T_C1= T_C1,T_C2=T_C2,T_m= T_Cm, Load=Load, a=a)
        elif a==4:
            T_L1= self.comp_L(L1)
            T_L1_2=self.comp_L(L1_2)
            T_L2= self.comp_L(L2)
            T_L2_2=self.comp_L(L2_2)
            T_C1= self.comp_C(C1)
            T_C2= self.comp_C(C2)
            T_imp, Z_in= self.input_Z(T_L1= T_L1, T_L2=T_L2, T_L1_2= T_L1_2, T_L2_2= T_L2_2, T_C1= T_C1,T_C2=T_C2, T_m= T_m, Load=Load, a=a)
        else:
                print(" Give a number: 0 - 4") 
        
        
        
        
        T_out, S_in, T_in= self.output_matrix(T_imp, V_in, Z_in, P_in)
        df_output['S_in']= S_in
        df_output['P_in']= np.real(df_output['S_in'])/2
        df_output['V_in']= pd.DataFrame([item[0][0] for item in T_in], index=self.f)
        df_output['I_in']= pd.DataFrame([item[1][0] for item in T_in], index=self.f)
        df_output['Z_in']= Z_in
        df_output['V_out']=[item[0][0] for item in T_out]
        df_output['I_out']=[item[1][0] for item in T_out]
        df_output['S_out']= df_output['V_out']* np.conjugate(df_output['I_out'])/2
        df_output['P_out']= np.real(df_output['S_out'])
        df_output['eff']= df_output['P_out']/df_output['P_in']
        
        df_output['T_out']= T_out
        df_output['T_in']= T_in
        df_output['T_imp']= T_imp
        
        if opt==False:
            T1= self.V1_Voltage(T_in, T_L1= T_L1, T_L2=T_L2, T_C1=T_C1, T_C2=T_C2 , T_m=T_Cm,  a=a)
            df_output['V1']=[item[0][0] for item in T1]
            df_output['I1']=[item[1][0] for item in T1]
            T2= self.V2_Voltage(T_in, T_L1= T_L1, T_L2=T_L2, T_C1=T_C1, T_C2=T_C2 , T_m=T_Cm,  a=a)
            df_output['V2']=[item[0][0] for item in T2]
            df_output['I2']=[item[1][0] for item in T2]
            df_output['DV']= df_output['V1']-df_output['V2']
            df_output['mat_loss']= df_output['DV']*np.conjugate(df_output['I2'])
        
        
        
        return df_output
        

            

            
            
            
            
            
            
