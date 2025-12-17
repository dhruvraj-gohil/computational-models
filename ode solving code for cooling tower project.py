import numpy as np
import pandas as pd

#constants
rho_a=1.184      # kg/m³, air density
c_pa=1005.0      # J/kg·K, specific heat of air
c_pw=4180.0      # J/kg·K, specific heat of water
h_fg=2.26e6      # J/kg, latent heat of vaporization
p_atm=101325.0   # Pa, atmospheric pressure

#Saturation vapor pressure (Tetens formula)
def p_sat_tetens(T_C):
    return 610.78*np.exp((17.269*T_C)/(T_C+237.3))

#Humidity ratio 
def Y_from_RH_T(RH, T_C):
    e_s=p_sat_tetens(T_C)
    e=RH*e_s
    return 0.622*e / (p_atm-e)

#Saturation humidity ratio at Tw
def Y_sat(Tw_C):
    e_s=p_sat_tetens(Tw_C)
    return 0.622*e_s / (p_atm-e_s)

#Cooling tower
def simulate_tower(n_baffles ,theta_deg):
    #Input parameters
    m_dot_w=1.0          # kg/s water
    m_dot_a=2.0          # kg/s dry air
    T_w_in=40.0          # °C
    T_a_in=25.0          # °C
    RH_a_in=0.50         # relative humidity
    H=20               # tower height
    A_cross=20      # m² cross-section area
    h0=9             # W/m²K
    N=10                # segments

    #Baffle parameters
    alpha1=0.7
    alpha2=0.15
    alpha3=0.01

    dz=H / N
    #Initial conditions
    T_w=np.zeros(N+1)
    T_a=np.zeros(N+1)
    Y=np.zeros(N+1)
    T_w[0]=T_w_in
    T_a[0]=T_a_in
    Y[0]=Y_from_RH_T(RH_a_in, T_a_in)

    #convective coefficient for baffles
    h_mult=(1.0+alpha1*np.sin(np.radians(theta_deg)))*(1.0+alpha2*n_baffles-alpha3*n_baffles*n_baffles)
    h_local=h0*h_mult
    k_c=0.0074 #  m/s 
    q_s_arr, q_l_arr, evap_arr=[], [], []
    #Finite difference method
    for i in range(N):
        Y_sat_i=Y_sat(T_w[i])
        m_evap=k_c*A_cross*(Y_sat_i-Y[i])
        q_s=h_local*A_cross*(T_w[i]-T_a[i])
        q_l=m_evap*h_fg

        T_w[i+1]=T_w[i]-(h_local*A_cross*(T_w[i]-T_a[i])/(m_dot_w*c_pw))*dz
        Y[i+1]=Y[i]+(k_c*A_cross*(Y_sat(T_w[i])-Y[i]) / m_dot_a)*dz
        T_a[i+1]=T_a[i]+((k_c*A_cross*h_fg*(Y[i+1]-Y[i]))+(h_local*A_cross*(T_w[i]-T_a[i]))*dz)/(m_dot_a*c_pa)
     
        q_s_arr.append(q_s)
        q_l_arr.append(q_l)
        evap_arr.append(m_evap)
    #Summarize results
    Q_total=sum(q_s_arr)+sum(q_l_arr)  # W
    evap_total=sum(evap_arr)             # kg/s
    T_w_out=T_w[-1]
    T_a_out=T_a[-1]
    #efficiency
    eff=abs(100.0*(T_w_in-T_w_out)/(T_w_in-T_a_in))
    return {
        "theta": theta_deg,
        "n" :n_baffles,
        "T_w_out_C": round(T_w_out, 3),
        "T_a_out_C": round(T_a_out, 3),
        "Q_kW": round(Q_total / 1000.0, 3),
        "Evap_kg_s": round(evap_total, 5),
        "Efficiency_%": round(eff, 4)
    }
results=[]
results.append(simulate_tower(0,0))
for n in range(1, 7):  # n = 0 to 6
    for theta in [0, 15, 30, 45, 60, 75,90,105,120,135,150]:
        results.append(simulate_tower(n, theta))
df = pd.DataFrame(results)
print(df.to_string(index=False))
