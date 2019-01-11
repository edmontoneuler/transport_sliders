import numpy as np

#Generic power-law transport integral 
def base_transport(m,E):
    t = np.zeros_like(E)
    for k in range(len(E)):
        if E[k] >0:
            t[k] = (E[k])**m
        else:
            t[k]=0
    return t

def transport(Egrid, vbands = [np.array([1, 1e8, 0])], cbands = [np.array([1, 2e8, 0.5])] ):
    """Returns the transport distribution on Egrid defined by """
    num_vbands = len(vbands)
    num_cbands = len(cbands)
    t = np.zeros_like(Egrid)
    N = len(t)
    
    for k in range(num_cbands):
        powerlaw = cbands[k][0]
        scale = cbands[k][1]
        gap = cbands[k][2]
        temp_t = np.zeros_like(Egrid)
        for i in range(N):
            if Egrid[i] - gap > 0:
                temp_t[i] = (Egrid[i]-gap)**powerlaw
        t = t + scale*temp_t
        
    for k in range(num_vbands):
        powerlaw = vbands[k][0]
        scale = vbands[k][1]
        gap = vbands[k][2]
        temp_t = np.zeros_like(Egrid)
        for i in range(N):
            if Egrid[i] - gap >0:
                temp_t[i] = (Egrid[i] - gap)**powerlaw
        t = t+ scale*np.flip(temp_t, 0) #Assumes symmetric energy grid 
    
    return t 

def fermi(E, mu, T):
    #Fermi Window (Zeroth Order)
    kB_eV = 8.617332385e-5
    ones = np.ones_like(E)
    return (1/(kB_eV*T))*np.exp((E-mu*ones)/(kB_eV*T))/(ones + np.exp((E-mu*ones)/(kB_eV*T)) )**2

def fermi1(E,mu,T):
    #First Order Fermi Window
    kB_eV = 8.617332385e-5
    ones = np.ones_like(E)
    return fermi(E,mu, T)*(E-mu*ones)/(kB_eV*T) 

def fermi2(E,mu,T):
    #Second Order Fermi Window
    kB_eV = 8.617332385e-5
    ones = np.ones_like(E)
    return fermi(E,mu,T)*((E-mu*ones)/(kB_eV*T))**2

def moments(transport, E, T_vec):
    """Calculates the thermoelectric moments
    of a provided transport distribution as a function
    of chemical potential and energy"""
    N = len(E)
    Nt = len(T_vec)
    moment0 = np.zeros([N, Nt])
    moment1 = np.zeros([N, Nt])
    moment2 = np.zeros([N, Nt])

    #Thermoelectric moment calculation
    for k in range(len(E)):
        percent = 100*k/N
        print("%.2f" % percent,'%', end="\r")
        for j in range(len(T_vec)):
            f0 = fermi(E, E[k], T_vec[j])
            f1 = fermi1(E, E[k], T_vec[j])
            f2 = fermi2(E, E[k], T_vec[j])
        
            moment0[k, j] = np.trapz(f0*transport, E)
            moment1[k, j] = np.trapz(f1*transport, E)
            moment2[k, j] = np.trapz(f2*transport, E)
    
    return moment0, moment1, moment2

def shift_fermi_int(integral, E, gap):
    dE = E[2] - E[1]
    index_shift = int(gap/dE)
    N = len(E)
    cond_int = np.zeros_like(integral)
    for k in range(1,N):
        i = N-k
        if i > index_shift:
            cond_int[i, :] = integral[i-index_shift, :]
        else:
            cond_int[i, :] = 0
    return cond_int

def shift_transport(t, E, gap):
    dE = E[2] =E[1]
    index_shift = int(gap/dE)
    N = len(E)
    shifted_t = np.zeros_like(t)
    for k in range(1, N):
        i = N-k
        if E[i] > gap:
            shifted_t[i] = t[i-index_shift]
        else:
            shifted_t[i] = 0
    return shifted_t

def get_index(value, vector):
    """Returns the index corresponding to the component 
    of the provided vector whose value is closest to the
    provided value (assumes uniform monotonic growth with index)"""

    array = np.array(vector)
    idx = (np.abs(array-value)).argmin()
    return idx


def clean_image(image, E, T_vec, plotlims=[-2.0, 2.0, 300, 800]):
    """Plots a function of mu and T within the plot limits provided"""
    for k in range(len(E)):
        if E[k] >plotlims[0]:
            low_E_index = k 
            break
    for k in range(len(E)):
        if E[k] >= plotlims[1]:
            high_E_index = k
            break
    for k in range(len(T_vec)):
        if T_vec[k] >= plotlims[2]:
            low_T_index = k
            break
    for k in range(len(T_vec)):
        if T_vec[k] > plotlims[3]:
            high_T_index = k 
            break
    trimmed_image = image[low_E_index:high_E_index,low_T_index:high_T_index]
    cleaned_image = np.flip(trimmed_image.T, 0) #Provides desired orientation
    
    return cleaned_image


