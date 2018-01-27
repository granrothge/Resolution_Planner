# A class that defines a slit pack object
import numpy as np
class Slit_pack(object):
    #class that holds slit package parameters and functions
    def __init__(self,d,r,name,R=0.05):
        self.d=d #distance between slits in m
        self.R=R #radius of slit package in m
        self.r=r #radius of curvature of slits
        self.name=name # name of slit package
    def __repr__(self):
        msg="%s\n distance between slits: %g \n radius of slit package: %g \n slit radius of curvature: %g \n" % (self.name,self.d,self.R,self.r)
        return msg
    def Beta(self,nu,v0,v):
        return (2.0*np.pi*nu/self.d*(self.R*self.R))*abs(1/v-1/v0)    
    def Fermi_T(self,nu,v):
          #transmission for neutrons of velovity v through a fermi chopper rotating at nu with  where optimal transmission is at v0  
          v0in=self.v0(nu)
          beta=self.Beta(nu,v0in,v)
          out=np.zeros(len(beta))
          blog=beta<=0.25
          out[blog]=1-8.0/3.0*beta[blog]*beta[blog]
          blog=(beta>0.25)&(beta<1.0)
          out[blog]=16.0/3.0*np.sqrt(beta[blog])-8*beta[blog]+8.0/3.0*beta[blog]*beta[blog]
          return out  
    def v0(self,nu):
        #optimal velocity for a chopper rotating at nu
        return self.r*4.0*np.pi*nu
    def Fermi_dt(self,nu,alpha=0.0):
       """
       nu is the frequency of the Fermi chopper in Hz
       alpha is the angular divergence of the beam perpendicular to the slit package
       """
       return 1/(4.0*np.pi*nu)*(self.d/self.R+2.0*alpha)

#definition of default slit packages
SEQ_100=Slit_pack(0.00203,0.58,'SEQ 100')
SEQ_700=Slit_pack(0.00356,1.53,'SEQ 700')
ARCS_100=Slit_pack(0.00152,0.58,'ARCS 100')
ARCS_300=Slit_pack(0.00305,1.00,'ARCS 300')
ARCS_700_2=Slit_pack(0.00152,1.53,'ARCS 700 2')
ARCS_700_3=Slit_pack(0.00356,1.53,'ARCS 700 3')    
