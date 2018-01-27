from numpy import sqrt
def E2V(E):
    """
      Takes a neutron energy in meV and converts it to velocity in m/s 
    """
# for energy in mev returns velocity in m/s
    return sqrt(E/5.227e-6)

def V2E(V):
    """
      Takes a neutron velocity in m/s and converts it to energy in meV
    """
# for  v in  m/s returns energy in meV
    return 5.227e-6*V*V
    
def E2K(E):
    """
        Takes a neutron of E mev and converts it to  k ^A-1 
    """
    return sqrt(E/2.0723)
def V2lambda(V):
    """
        Takes a neutron velocity V in m/s and converts it to lambda in angstroms
    """
    return(3956/V)
