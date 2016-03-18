import numpy as np
from matplotlib import pyplot as plt
plt.style.use('ggplot')
D=np.arange(0,5,0.1)
Dmm=D
Dcm=D*0.1
N0=0.08*1.e6*1.e-1 #m**{-3} mm^{-1}
R=1.
theLambda=41*R**(-0.21)
curve1=N0*np.exp(-theLambda*Dcm)
R=5.
theLambda=41*R**(-0.21)
curve2=N0*np.exp(-theLambda*Dcm)
R=25.
theLambda=41*R**(-0.21)
curve3=N0*np.exp(-theLambda*Dcm)
fig, ax = plt.subplots(1,1,figsize=(10,8))
ax.semilogy(D,curve1,label='1 mm/hr')
ax.semilogy(D,curve2,'r-',label='5 mm/hr')
ax.semilogy(D,curve3,'g-',label='25 mm/hr')
ax.set_xlabel('Drop diameter (mm)')
ax.set_ylabel('$n(D) m^{-3} mm^{-1}$')
ax.set_title('Marshall Palmer distribution for three rain rates')
ax.legend(loc='best')

int_D = (np.sum((Dcm[:-1]*N0*np.exp(Dcm[:-1]*-theLambda)*np.diff(Dcm))))/(np.sum(N0*
np.exp(Dcm[:-1]*-theLambda)*np.diff(Dcm)))

print("mean diameter of M-P distribution is {}".format(int_D))
print("lambda^-1 = {}".format(1/theLambda))

grav = 9.8 # m/s^2
rho_l = 1000 # kg/m^3
rho_air = 1.225 # kg/m^3
Dnew = np.arange(0,5,0.1)
Dcm_new = Dnew*0.1
fallspeed = np.sqrt((rho_l/rho_air)*grav*Dcm_new[:-1])
R = 15 # mm/hr
theLambda=41*R**(-0.21)

rainrate = np.sum((N0*np.exp(Dcm_new[:-1]*-theLambda))*((np.pi*(Dcm_new[:-1])**3)/6)*
fallspeed*np.diff(Dcm_new))

print(rainrate)
