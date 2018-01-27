def Bmat_gen(latt_con,latt_ang):
     """function to generate B Matrix given the lattice constants
     and the lattice angles
     inputs:
        lattice constants in angstroms
        lattice angles in degrees
     """
     from numpy import pi
     latt_ang=latt_ang*pi/180.0
     reclist=gen_rec_latt(latt_con,latt_ang)
     return Bmat(latt_con,reclist[0],latt_ang,reclist[1])

def gen_rec_latt(latt_con,latt_ang,times2pi=1):
     """ function to lattice parameters of the reciprocal lattice given the real 
     space lattice parameters
     inputs:
     lattice constants in angstroms
     lattice angles in radians
     """
     from numpy import prod,sum,cos,sin,arccos,zeros,sqrt,pi     
     rec_latt=zeros(3)
     rec_ang=zeros(3)
     #print latt_ang
     Vol=prod(latt_con)*sqrt((1.0+2.0*prod(cos(latt_ang))-sum((cos(latt_ang))**2)))
     print (Vol)
     rec_latt[0]=latt_con[1]*latt_con[2]*sin(latt_ang[0])/Vol
     rec_latt[1]=latt_con[0]*latt_con[2]*sin(latt_ang[1])/Vol
     rec_latt[2]=latt_con[0]*latt_con[1]*sin(latt_ang[2])/Vol
     rec_ang[0]=arccos((cos(latt_ang[1])*cos(latt_ang[2])-cos(latt_ang[0]))/abs(sin(latt_ang[1])*sin(latt_ang[2])))
     #print rec_ang[0]
     rec_ang[1]=arccos((cos(latt_ang[0])*cos(latt_ang[2])-cos(latt_ang[1]))/abs(sin(latt_ang[0])*sin(latt_ang[2])))
     rec_ang[2]=arccos((cos(latt_ang[0])*cos(latt_ang[1])-cos(latt_ang[2]))/abs(sin(latt_ang[0])*sin(latt_ang[1])))
     if times2pi:
       rec_latt=rec_latt*2.0*pi
     return [rec_latt, rec_ang]
def Bmat(a,b,alpha,beta):
     """ given the real space and reciprocal space lattice parameters,  
     determine the B matrix
     """
     from numpy import cos, sin, pi,zeros
     B=zeros([3,3])
     B[0][0]=b[0]
     B[0][1]=b[1]*cos(beta[2])
     B[0][2]=b[2]*cos(beta[1])
     B[1][1]=b[1]*sin(beta[2])
     B[1][2]=-b[2]*sin(beta[1])*cos(alpha[0])
     B[2][2]=2.0*pi/a[2]
     return B
def unit_vec_omega(omega):
    """
    """
    from numpy import sin,cos
    return [cos(radians(omega)),0,-sin(radians(omega))]
    
def two_vec2U(hkl1,hkl2,B,rot1,rot2,twotheta1,twotheta2,ccwrot=0,debug=0):
     """
     """
     from numpy.linalg import norm
     from numpy import cross
     omega1=rot2omega(rot1,twotheta1,ccwrot)
     omega2=rot2omega(rot2,twotheta2,ccwrot)
     uw1=unit_vec_omega(omega1)
     uw2=unit_vec_omega(omega2)
     twmat=transpose(twovec2_3unit(uw1,uw2))
     qc1=dot(B,hkl1)
     qc2=dot(B,hkl2)
     tcmat=transpose(twovec2_3unit(qc1,qc2))
     # because tcmat is an orthonormal matrix, the transpose is the inverse
     U=dot(twmat,transpose(tcmat))
     if debug:
       print ('omega1= %g omega2=%g \n' %(omega1,omega2))
       print ("uw1\n")
       print (uw1)
       print ("\nuw2\n")
       print (uw2)
       print ("\n")
       print ("twmat\n")
       print (twmat)
       print ("\n")
       print ("tcmat\n")
       print (tcmat)
       print ("\n")
     return U
def twovec2_3unit(prim,sec):
    """
    take a primary and a secondary vector and generate 3 orthogonal unit vectors
    out1 is along the primary 
    out2 is in the plane and orthogonal to out1
    out3 is orthogonal to the plane
    """
    from numpy.linalg import norm
    from numpy import cross
    out1=prim/norm(prim)
    out3=cross(prim,sec)/norm(cross(prim,sec))
    out2=cross(out3,out1)
    return[out1,out2,out3]

def rot2omega(rot,twotheta,ccwrot):
    """
    sign dependent determination of omega from rot and twotheta
    """
    if ccwrot:
       return(rot-twotheta/2.0)
    else:
       return(-rot-twotheta/2.0)
       
def omega2rot(omega,twotheta,ccwrot):
    """
    sign dependent determination of omega from rot and twotheta
    """
    if ccwrot:
       return(omega+twotheta/2.0)
    else:
       return(-omega-twotheta/2.0)       
def mslice_psi_2_rot(psi,U,B,hkl,twotheta):

    """
    """
    from numpy.linalg import norm
    UB=dot(U,B)
    Qw=dot(UB,hkl)
    uw=Qw/norm(Qw)
    omega=degrees(arccos(-uw[0]))+psi
    return omega2rot(omega,twotheta,0)
     
    
          
     
