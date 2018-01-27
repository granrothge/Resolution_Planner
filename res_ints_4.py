#python routines for estimating energy resolution and intensity
#G. E. Granroth
#Updated 2-19-2013 to include tube efficiency
import sys
#sys.path.append('/SNS/users/19g/SEQUOIA/commissioning/python')
from unit_convert import E2V,E2K
from numpy import pi, log,exp,sqrt,tanh,linspace, radians,zeros
from slit_pack import Slit_pack
from scipy.interpolate import interp1d
from pylab import figure, plot, subplot, show, xlabel, ylabel, title
from UB import Bmat_gen,gen_rec_latt,Bmat

class Chopper_spec(object):
     """
     class to define a chopper spectrometer
     w is the width fo the detector in m
     h is the height of the detector in m
     L is a three vector all in m
     	L[0]= moderator to chopper distance
        L[1]= chopper to sample distance
        L[3]= sample to detector distance
     slit_pack is and instance of the slit_pack class
     sw is the sample width in m
     sh is the sample height in m
     He_press is the He pressure in ATM
     mod_file is a moderator file provided by the neutronics group
     """
     
     def __init__(self,instr_name,L,slit_pack,hphilims,vphilims,w=0.0254,h=0.01,sw=0.05,sh=0.05,He_press=10,He_T=300.0,mod_file='source_sct521_bu_17_1.dat'):   
         self.instr_name=instr_name
         self.L=L
         self.w=w
         self.h=h
         self.sw=sw
         self.sh=sh
         self.hphilims=hphilims
         self.vphilims=vphilims
         self.slit_pack=slit_pack
         self.He_press=He_press
         self.He_T=He_T
         self.pix_vol=pi*self.w*self.w/4.0*self.h
         self.num_density=self.He_press*7.336e26/self.He_T
         self.mod_file=mod_file
         self.I_func=read_mod_file(mod_file)
     def domega_in(self,Ei,Ef,nu):
         """
         """
         dtm=H2Omod_dt(Ei)
         dtc=self.slit_pack.Fermi_dt(nu)
         dtd=det_dt(self.w,Ef)
         return domega(Ei,Ef,self.L,dtm,dtc,dtd)
     def dE(self,Ei,nu):
         """
         """
         vi=E2V(Ei)
         dtc=self.slit_pack.Fermi_dt(nu)
         return 2*5.227e-6*vi**3.0*dtc/self.L[0]
     def flux(self,filename,Ei,Ef,nu,He_flag=0):
         "calculate the flux on sample given an Ei and and chopper frequency"
         v=E2V(Ei)
	     #print v
	     #v0=self.slit_pack.v0(nu)
	     #print v0
         T=self.slit_pack.Fermi_T(nu,v)
         #I_func=read_mod_file(filename)
         I_return=self.I_func(Ei/1000.0)*self.dE(Ei,nu)*self.sw*self.sh/((self.L[0]+self.L[1])**2.0)*T
         if He_flag:
	         I_return=I_return*(1-exp(-1.1734E-21/(E2V(Ef))*self.num_density*self.w))
         return I_return


def domega(Ei,Ef,L,dtm,dtc,dtd):
     """
     provides the energy resolution of a direct chopper spectrometer given
     Ei: incident energy in meV
     Ef: fineal energy in meV
     A three element tuple
     L[0]= moderator to chopper distance
     L[1]= chopper to sample distance
     L[3]= sample to detector distance
     dtm = moderator pulse width (s)
     dtc= chopper pulse width (s)
     dtd= detector time uncertainty (s)
     """
     mn=1.674e-5/1.602 #mass of the neutron in meV
     vi=E2V(Ei)
     vf=E2V(Ef)
     return mn*sqrt(((vi**3.0)/L[0]+(vf**3.0)*L[1]/L[0]/L[2])**2.0*(dtm**2.0)+((vi**3.0)/L[0]+(vf**3.0)*(L[1]+L[0])/L[0]/L[2])**2.0*(dtc**2.0)+((vf**3.0)/L[2])**2.0*(dtd)**2.0)
   
def H2Omod_dt(E):
    """
    returns the time width of the neutron distribution as a function of energy
    E(meV)	
    """
    E=E
    x=log(E)
    p=[-0.4494,-0.046,4.3672,0.8530,3.7389,0.1271]
    y=exp(m1tanhm2(x,p))
    return y*1e-6
    
def det_dt(w,Ef):
    """
    w is the detector width in (m)
    Ef is the final energy in meV
    """
    return w/E2V(Ef)
def m1tanhm2(x,p):
    """
    """
    m=[p[0],p[1]]
    #m[0]=p[0]
    #m[1]=p[1]
    x0=p[2]
    w=p[3]
    y0=p[4]
    A=p[5]
    return (1+tanh((x-x0)/w))/2.0*m[0]*x+(1-tanh((x-x0)/w))/2.0*m[1]*x+y0+A*tanh((x-x0)/w)
def read_mod_file(filename):
    """
    a function to read a moderator file from the neutronics group
    """
    fid=open(filename)
    dattmp=fid.readlines()
    fid.close()
    idx=0
    E=[]
    flux=[]
    while '#' in dattmp[idx]:
       idx=idx+1
    while  not ('#' in dattmp[idx]):
       #print dattmp[idx]
       tmp1=dattmp[idx].split()
       if len(tmp1)>0:
          E.append(eval(tmp1[0]))
          flux.append(eval(tmp1[2]))
       idx=idx+1 
    flux_func=interp1d(E,flux,kind='linear')               
    return flux_func
def plot_flux(nu,Ei,Ef,Spec,He_flag=1):
     """
      plot_flux(nu,Ei,Ef,Spec)
      give a range of chopper frequencies (nu) (Hz)
      an incident energy Ei (meV) and a final energy Ef (meV)
      and an instance of a spectrometer class (examples are given in the bottome of this file)
      plot a number proportional to flux and the resolution as a function of chopper frequency  
     """
     dw=Spec.domega_in(Ei,Ef,nu)
     I=zeros(len(nu))
     #for idx in range(len(nu)):
      # I.append(Spec.flux('source_sct521_bu_17_1.dat',Ei,Ef,nu[idx],He_flag=He_flag))
     I=Spec.flux('source_sct521_bu_17_1.dat',Ei,Ef,nu,He_flag=He_flag)
     figure()
     subplot(2,1,1)
     plot(nu,I,'bo')
     ylabel('I (arb. units)')
     title('$E_i$ = %g (meV),$E_f$ = %g (meV), $\hbar\omega$ = %g (meV)\n Slit_Pack:%s '%(Ei,Ef,Ei-Ef,Spec.slit_pack.name))
     subplot(2,1,2)
     plot(nu,dw,'bo')
     ylabel('$d(\hbar\omega)$ (meV)')
     xlabel('$\\nu$ (Hz)')
     show()
def plot_res_omega(nu,Ei,omega,Spec):
    """
        plot_res_omega(nu,Ei,omega,Spec)
        Plot the energy resolution as a function of omega in meV
	nu is the Fermi chopper speed
	Ei is the incident energy
	omega is a numpy array of energy transfers
	Spec is one of the spectrometers defined at the end of the file.	
	
    """     
    Ef=Ei-omega
    dw=Spec.domega_in(Ei,Ef,nu)
    figure()
    plot(omega,dw,'bo')
    ylabel('$d(\hbar\omega)$ (meV)')
    xlabel('$\hbar\omega$ (meV)')
    title('$E_i$ = %d meV $\\nu$ = %d Hz \n SlitPack:%s'% (Ei,nu,Spec.slit_pack.name)) 
    show()
    return [omega,dw]
    
def plot_qrange(Ei, wmin,spec,UB=[[1,0,0],[0,1,0],[0,0,1]]):
   """
    plot_qrange(Ei, wmin,spec,UB)
    given an Ei and a minimum energy transfer (wmin) for a given spectrometer (spec) and with a 
    crystal parameters and orientation defined by UB, plot the Q ranges accesible by the instrument
    predefined values for several chopper spectrometers are given at the end of this file.
   """
   ki=E2K(Ei)
   omega=linspace(wmin,Ei*0.9,100)
   Ef=Ei-omega
   kf=E2K(Ef)
   hphilims=radians(spec.hphilims)
   vphilims=radians(spec.vphilims)
   Qxmax=-kf*sin(hphilims[1])
   Qxmin=-kf*sin(hphilims[0])
   Qxmin2=-kf*sin(hphilims[2])
   Qxmax2=-kf*sin(hphilims[3])
   Qymax=-kf*sin(vphilims[1])
   Qymin=-kf*sin(vphilims[0])
   Qymin2=-kf*sin(vphilims[2])
   Qymax2=-kf*sin(vphilims[3])
   Qzmax=ki-kf*cos(hphilims[1])
   Qzmin=ki-kf*cos(hphilims[0])
   Qzmax2=ki-kf*cos(hphilims[3])
   Qzmin2=ki-kf*cos(hphilims[2])
   Qmins=array([Qxmin,Qymin,Qzmin])
   Qmins2=array([Qxmin2,Qymin2,Qzmin2])
   Qmaxs=array([Qxmax,Qymax,Qzmax])
   Qmaxs2=array([Qxmax2,Qymax2,Qzmax2])
   hklmins=dot(UB,Qmins)
   hklmaxs=dot(UB,Qmaxs)
   hklmins2=dot(UB,Qmins2)
   hklmaxs2=dot(UB,Qmaxs2)
   figure()
   hold('on')
   xlbs=['$Q_x$','$Q_y$','$Q_z$']
   for idx in range(3):
     subplot(2,2,idx+1)
     plot_qlims(hklmins,hklmaxs,hklmins2,hklmaxs2,omega,idx)
     xlabel(xlbs[idx])
   subplot(2,2,4)
   abs_tt=abs(array(hphilims))
   tthetamin=min(abs_tt)
   tthetamax=max(abs_tt)
   Qmin=sqrt(ki*ki+kf*kf-2.*ki*kf*cos(tthetamin))
   Qmax=sqrt(ki*ki+kf*kf-2.*ki*kf*cos(tthetamax))
   plot(Qmin,omega,'r')
   plot(Qmax,omega,'b')
   xlabel('|Q|')
   ylabel('$\omega$')
   show()
        
def plot_qlims(mins,maxs,mins2,maxs2,omega,idx):
   """
   """
   plot(mins[idx,:],omega,'b')
   plot(maxs[idx,:],omega,'b')
   plot(mins2[idx,:],omega,'r')
   plot(maxs2[idx,:],omega,'r')
   ylabel('$\omega$')
   
    
       
       

#define default slit packages
SEQ_100=Slit_pack(0.00203,0.58,'SEQ-100-2.03-AST')
SEQ_700=Slit_pack(0.00356,1.53,'SEQ-700-3.56-AST')
ARCS_100=Slit_pack(0.00152,0.58,'ARCS 100')
ARCS_300=Slit_pack(0.00305,1.00,'ARCS 300')
ARCS_700_2=Slit_pack(0.00152,1.53,'ARCS 700 2')
ARCS_700_3=Slit_pack(0.00356,1.53,'ARCS 700 3')
ARCS_700_sf=Slit_pack(0.0005,1.53,'ARCS-700-0.5-AST')
SEQ_1000=Slit_pack(0.0015,1.83,'SEQ-1000-1.5-AST')
SEQUOIA=Chopper_spec('SEQUOIA',[18.0,2.0,5.5],SEQ_100,[2.1,60.0,-5.3,-30.0],[6.7,18.0,-7.5,-18.0])
SEQUOIA_sloppy=Chopper_spec('SEQUOIA',[18.0,2.0,5.5],SEQ_700,[2.1,60.0,-5.3,-30.0],[6.7,18.0,-7.5,-18.0])
SEQUOIA_700_superfine=Chopper_spec('SEQUOIA',[18.0,2.0,5.5],ARCS_700_sf,[2.1,60.0,-5.3,-30.0],[6.7,18.0,-7.5,-18.0])
SEQUOIA_1000=Chopper_spec('SEQUOIA',[18.0,2.0,5.5],SEQ_1000,[2.1,60.0,-5.3,-30.0],[6.7,18.0,-7.5,-18.0])
ARCS=Chopper_spec('ARCS',[11.6,2.0,3.0],ARCS_100,[2.1,90.0,-5.3,-30.0],[6.7,30.0,-7.5,-30.0])
ARCS_700_fine=Chopper_spec('ARCS',[11.6,2.0,3.0],ARCS_700_2,[2.1,90.0,-5.3,-30.0],[6.7,30.0,-7.5,-30.0])
ARCS_700_sloppy=Chopper_spec('ARCS',[11.6,2.0,3.0],ARCS_700_3,[2.1,90.0,-5.3,-30.0],[6.7,30.0,-7.5,-30.0])
ARCS_700_superfine=Chopper_spec('ARCS',[11.6,2.0,3.0],ARCS_700_sf,[2.1,90.0,-5.3,-30.0],[6.7,30.0,-7.5,-30.0])
