#Preambule
#   TYPE %matplotlib qt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from numpy import ma

#%matplotlib.use("Agg")

#Domain and parameters
#fig=plt.figure(figsize=(10,6), dpi=100)
fig = plt.figure(num = 0, figsize = (12, 8))#, dpi = 100)
fig.suptitle("Silvernet", fontsize=12)
ax01 = plt.subplot2grid((3, 3), (0, 0))
ax02 = plt.subplot2grid((3, 3), (0, 1))
ax03 = plt.subplot2grid((3, 3), (0, 2))
ax04 = plt.subplot2grid((3, 3), (1, 2))
ax05 = plt.subplot2grid((3, 3), (2, 2))
ax = plt.subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)
#ax04 = ax03.twinx()
xvec=[0,1,2]
yvec=[0,1,2]

plt.xlim([0,2])
plt.ylim([0,2])
plt.xlabel('Zonal distance',fontsize=13)
plt.ylabel('Meridional distance',fontsize=13)
ax.set_title('Model Horizontal Field')
ax01.set_title('Radiation')
ax01.set_ylim(-750,1000)
ax01.set_xlim(0,100)
ax02.set_title('Temperature')
ax02.set_ylim(240,330)
ax02.set_xlim(0,100)
ax03.set_title('Net surface energy balance')
ax03.set_ylim(-300,500)
ax03.set_xlim(0,100)
ax04.set_title('Unknown')
ax05.set_title('Unknown')

#Constants
albedo=[0.1,0.1,0.1,0.3,0.3,0.3,0.3,0.3,0.3]
S=1366.
sigma=5.670*10**(-8)
emis_a=0.7
emis_s=1.0
tau_a=0.8
alb_a=1-tau_a
apt=1.25
Lv=2500000.
Cp=1000.
Cpw=4181.
Cps=800.
gamma=Cp/Lv
Rd=287
Rv=461.5
eps=Rd/Rv
densw=999.7
denss=2000.
densa=1.2
lat=45. # in Summer!
keff=1. #0.3 dry, 2.2 wet
u2=3.
omega=2.*np.pi/(24.)
OMEGA=2.*np.pi/(24.*3600.)
K=0.5*10**6. #0.24 fresh, 0.74 saturated
Kw=0.14*10**6.
Depth=0.5 #of subsurface layer
f=7.2921*10.**(-5.)*np.sin(lat*2.*np.pi/360.)

#Per grid
LENGTH=1000000.
Thickness=0.25
Thicknessa1=1000.
Thicknessa2=1000.
VOLUME=LENGTH**2*Thickness
VOLUMEa1=LENGTH**2*Thicknessa1
VOLUMEa2=LENGTH**2*Thicknessa2
AREA=LENGTH**2

#Initial conditions
U=np.array([1.,1.,1.,1.,1.,1.,1.,1.,1.])
V=np.array([1.,1.,1.,1.,1.,1.,1.,1.,1.])
W=[0,0,0,0,0,0,0,0,0]
P=[1000,1000,1000,1000,1000,1000,1000,1000,1000]
Ts=[280,280,280,283,280,280,280,290,280]
Ta=[296,295,290,290,291,290,293,290,294]
Ts2=[278,278,278,278,278,278,278,278,278]#5 m onder de grond
Ta=np.array(Ts)-10.
Ta2=np.array(Ta)-10.
rho=densa
rh=25.

#Assumptions:
#--- 

#Fluxes
def Qs(a,time): #SWR
    return max(0,np.sin(time*2.*np.pi/24.))*(1.-a)*S*(tau_a)
    
def Ql(Tatm,Tsurf): #LWR: Stefon-Boltzmann
    return emis_s*0.5*emis_a*sigma*Tatm**4.-emis_s*sigma*Tsurf**4.
    
def Qg(Tsurf,Tsurf2): #GHF
    return -keff*(Tsurf-Tsurf2)/(0.2)
        
def Qe(a,Tatm,Tsurf,time,Tsurf2): #LHF: Priest-Taylor
    qs=611.2*np.exp(17.67*((Tatm-273.)/(Tatm-273.+243.5)))
    #es=611.*np.exp(17.2694*((Tatm-273.16)/(Tatm-35.86)))
    scc=eps*Lv*qs/(Rd*Tatm**2)
    A=6.12
    m=7.59
    Tn=240.73
    es=A*10**((m*(Tatm-273))/(Tatm-273+Tn))/100.
    ea=0.5*es
    #scc=2*10**(-11)*np.exp(0.0597*Tatm)
    return (0.408*scc*(-(Qs(a,time)+Ql(Tatm,Tsurf))+Qg(Tsurf,Tsurf2))+gamma*900.*u2*(es-ea)/(Tsurf+273.))/(scc+gamma*(1.+0.34*u2))
    #return apt*scc*(-(Qs(a,time)+Ql(Tatm,Tsurf))+Qg(Tsurf,Tsurf2))/(scc+gamma)
    
def Qh(Tatm,Tsurf): #SHF
    return rho*Cp*(Tatm-Tsurf)/rh
    
def EBAL(a,Tatm,Tsurf,time,Tsurf2):
    return Qs(a,time)+Ql(Tatm,Tsurf)+Qe(a,Tatm,Tsurf,time,Tsurf2)+Qg(Tsurf,Tsurf2)+Qh(Tatm,Tsurf)

def ABAL(a,Tatm,Tsurf,time,Tatm2,Tsurf2):
    return emis_a*(emis_s*sigma*Tsurf**4.-sigma*Tatm**4.-Qh(Tatm,Tsurf)-Qe(a,Tatm,Tsurf,time,Tsurf2)+0.5*emis_a*sigma*Tatm2**4.)
   
def ABAL2(a,Tatm,Tsurf,time,Tatm2,Tsurf2):
    #return emis_a*(emis_s*sigma*Tsurf**4.-sigma*Tatm**4.-Qh(Tatm,Tsurf)-Qe(a,Tatm,Tsurf,time)+0.5*emis_a*sigma*Tatm2**4.)

    return emis_a*(0.5*emis_a*sigma*Tatm**4.-sigma*Tatm2**4.-Qh(Tatm,Tsurf)*(1-emis_a)+emis_s*sigma*Tsurf**4.*(1-emis_a)-Qe(a,Tatm,Tsurf,time,Tsurf2)*(1-emis_a))
   
# Data Placeholders
yp1=np.zeros(0)
yp2=np.zeros(0)
yp3=np.zeros(0)
yp4=np.zeros(0)
yp5=np.zeros(0)
yp6=np.zeros(0)
yp7=np.zeros(0)
yp8=np.zeros(0)
yp9=np.zeros(0)
t=np.zeros(0)
x = 0.0
Uold=np.array(U)
Vold=np.array(V)
DPDX=np.array([0,0,0,0,0,0,0,0,0])
DPDY=np.array([0,0,0,0,0,0,0,0,0])
v_g=3.
u_g=3.

#Variables
TA=Ta
TA2=Ta2
TS=Ts
Evec=np.zeros(9)
Avec=np.zeros(9)
xvec,yvec=np.meshgrid(xvec,yvec)

#Radiation plot
LIJN1, = ax01.plot(t,yp1,'b-', label="QS")
LIJN2, = ax01.plot(t,yp2,'r-', label="QL")
LIJN3, = ax01.plot(t,yp3,'k-', label="QE")
LIJN4, = ax01.plot(t,yp4,'y-', label="QG")
LIJN5, = ax01.plot(t,yp5,'g-', label="QH")

ax01.legend([LIJN1,LIJN2,LIJN3,LIJN4,LIJN5], [LIJN1.get_label(),LIJN2.get_label(),LIJN3.get_label(),LIJN4.get_label(),LIJN5.get_label()],bbox_to_anchor=(0., 1.02, 1., .102), loc=2,
           ncol=5, mode="expand", borderaxespad=3.,fontsize=9,prop={'size':7})

#Temperature plot
LIJN6, = ax02.plot(t,yp6,'k-', label="Ta")
LIJN7, = ax02.plot(t,yp7,'b-', label="Ts")
LIJN9, = ax02.plot(t,yp9,'g-', label="Ta2")
ax02.legend([LIJN6,LIJN7,LIJN9], [LIJN6.get_label(),LIJN7.get_label(),LIJN9.get_label()],bbox_to_anchor=(0., 1.02, 1., .102), loc=2,
           ncol=5, mode="expand", borderaxespad=3.,fontsize=9,prop={'size':7})

#Net surface radiation
LIJN8, = ax03.plot(t,yp8,'k-', label="E_bal")
ax03.legend([LIJN8], [LIJN8.get_label()],bbox_to_anchor=(0., 1.02, 1., .102), loc=2,
           ncol=5, mode="expand", borderaxespad=3.,fontsize=9,prop={'size':7})


#Timeloop
def plotter(vec):
    return np.array(np.transpose([vec[0:3],vec[3:6],vec[6:9]]))
    
def init():
    imobj.set_data(Ta)
    time_text.set_text('time = 0.0')

    return imobj , time_text
    
def animate(self):
    global data
    global LIJN1
    global LIJN2
    global LIJN3
    global LIJN4
    global LIJN5
    global LIJN6
    global LIJN7
    global LIJN8
    global LIJN9
    global t
    global x
    global yp1
    global yp2
    global yp3
    global yp4
    global yp5
    global yp6
    global yp7
    global yp8
    global yp9
    global qk
    Evec=[]
    Avec=[]
    Avec2=[]
    QS=[]
    QL=[]
    QE=[]
    QG=[]
    QH=[]
    
    for j in range(0,9):
        Evec.append(EBAL(albedo[j],Ta[j],Ts[j],x,Ts2[j]))
        Avec.append(ABAL(albedo[j],Ta[j],Ts[j],x,Ta2[j],Ts2[j]))
        Avec2.append(ABAL2(albedo[j],Ta[j],Ts[j],x,Ta2[j],Ts2[j]))
        QS.append(Qs(albedo[j],x))
        QL.append(Ql(Ta[j],Ts[j]))
        QE.append(Qe(albedo[j],Ta[j],Ts[j],x,Ts2[j]))
        QG.append(Qg(Ts[j],Ts2[j]))
        QH.append(Qh(Ta[j],Ts[j]))
    for j in range(0,3):
        Ta[j]=Ta[j]+Avec[j]*AREA*3600./(Cp*VOLUMEa1*densa)
        Ta2[j]=Ta2[j]+Avec2[j]*AREA*3600./(Cp*VOLUMEa2*densa)
        Ts[j]=Ts[j]+Evec[j]*AREA*3600./(Cpw*VOLUME*densw)
        Ts2[j]=300+5*np.exp(np.sqrt(OMEGA/2./Kw)*Depth)*np.cos(omega*x+np.sqrt(OMEGA/2./Kw)*Depth)
    for j in range(3,9):
        Ta[j]=Ta[j]+Avec[j]*AREA*3600./(Cp*VOLUMEa1*densa)
        Ta2[j]=Ta2[j]+Avec2[j]*AREA*3600./(Cp*VOLUMEa2*densa)
        Ts[j]=Ts[j]+Evec[j]*AREA*3600./(Cps*VOLUME*denss)
        Ts2[j]=300+5*np.exp(np.sqrt(OMEGA/2./K)*Depth)*np.cos(omega*x+np.sqrt(OMEGA/2./K)*Depth)
    

    
    for j in range(0,6):
        DPDX[j]=5*np.cos(omega*x)+v_g*densa*f
    for j in range(6,9):
        DPDX[j]=0.5*5*np.cos(omega*x)+v_g*densa*f #0=v_g=B/(rho*f)
    for j in range(0,9):
        DPDY[j]=u_g*densa*f
    
    for j in range(0,9):
        U[j]=Uold[j]+1.*((-1./densa)*DPDX[j]+f*Vold[j])
    for j in range(0,9):
        V[j]=Vold[j]-1.*(f*Uold[j]+(-1./densa)*DPDY[j])
        
    for j in range(0,9):#####################################Hier onder kijken!
        U[j]=0.0001*(OMEGA*np.sin(omega*(x+0.25*2*np.pi/omega))-f*np.sin(3600*f*x+np.pi/2.))/(densa*(f**2.-OMEGA**2.))
    for j in range(0,9):
        V[j]=0.0001*f*(np.cos(omega*x+np.pi/2.)-np.cos(3600*f*x+np.pi/2.))/(densa*(f**2.-OMEGA**2.))
    
    
    yp1=np.append(yp1,np.mean(QS))
    yp2=np.append(yp2,np.mean(QL))
    yp3=np.append(yp3,np.mean(QE))
    yp4=np.append(yp4,np.mean(QG))
    yp5=np.append(yp5,np.mean(QH))
    yp6=np.append(yp6,np.mean(Ta))
    yp7=np.append(yp7,np.mean(Ts))
    yp8=np.append(yp8,np.mean(Evec))
    yp9=np.append(yp9,np.mean(Ta2))
    t=np.append(t,x)

    x += 1.0    
    
    LIJN1.set_data(t,yp1)
    LIJN2.set_data(t,yp2)
    LIJN3.set_data(t,yp3)
    LIJN4.set_data(t,yp4)
    LIJN5.set_data(t,yp5)
    LIJN6.set_data(t,yp6)
    LIJN7.set_data(t,yp7)
    LIJN8.set_data(t,yp8)
    LIJN9.set_data(t,yp9)
    
    time_text.set_text('time = %.1f' % self )
    SB_text.set_text('SB = %.1f' % np.mean(Evec) )
    AB_text.set_text('AB = %.1f' % np.mean(Avec) )
    AB2_text.set_text('A2B = %.1f' % np.mean(Avec2) )
    new=plotter(Ts)
    imobj.set_data(new)
       
    imobj.set_zorder(0)
    Quivers.set_UVC(plotter(U),plotter(V))
    
    ax1=plt.axhline(xmin=0,xmax=1,y=0,linewidth=3,color='k')
    ax2=plt.axvline(ymin=0,ymax=1,x=0,linewidth=3,color='k')
    ax3=plt.axhline(xmin=0,xmax=1,y=2,linewidth=3,color='k')
    ax4=plt.axvline(ymin=0,ymax=1,x=2,linewidth=3,color='k')
    ax5=plt.axvline(ymin=0,ymax=1,x=0.5,linewidth=3,color='k',ls='--')
    #hoi=plt.contour([0,1,2],[0,1,2],plotter(Ta),15,animated=True)
    
    
    if x >= 0.:
        LIJN1.axes.set_xlim(x-100.+1.0,x+1.0)
        LIJN2.axes.set_xlim(x-100.+1.0,x+1.0)
        LIJN3.axes.set_xlim(x-100.+1.0,x+1.0)
        LIJN4.axes.set_xlim(x-100.+1.0,x+1.0)
        LIJN5.axes.set_xlim(x-100.+1.0,x+1.0)
        LIJN6.axes.set_xlim(x-100.+1.0,x+1.0)
        LIJN7.axes.set_xlim(x-100.+1.0,x+1.0)
        LIJN8.axes.set_xlim(x-100.+1.0,x+1.0)
        LIJN9.axes.set_xlim(x-100.+1.0,x+1.0)
    
    qk = plt.quiverkey(Quivers, 0.1, 0.7, 2, r'$2 \frac{m}{s}$', labelpos='W',
                       fontproperties={'weight': 'bold'})
    
    return imobj , Quivers,time_text,ax1,ax2,ax3,ax4,ax5,SB_text,AB_text,AB2_text,LIJN1,LIJN2,LIJN3,LIJN4,LIJN5,LIJN6,LIJN7,LIJN8,LIJN9,qk


def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

data = list(Ts)
#imobj = plt.contourf(xvec,yvec,plotter(data),animated=True,aspect=0.5,extent=[0, 2.0, 0.0, 2.0], alpha=1.0, zorder=1)
imobj=plt.imshow(plotter(data), cmap='jet', animated=True,aspect=0.5,extent=[0, 2.0, 0.0, 2.0], alpha=1.0, zorder=1)
plt.colorbar(orientation='horizontal')
plt.clim(240,310)
X,Y=np.meshgrid([0,1,2],[2,1,0])
Quivers=plt.quiver(plotter(U),plotter(V),units='inches')
qk = plt.quiverkey(Quivers, 0.1, 0.7, 2, r'$2 \frac{m}{s}$', labelpos='W',
                   fontproperties={'weight': 'bold'})

l, r, b, top = plt.axis()
dx, dy = r - l, top - b
plt.axis([l - 0.05*dx, r + 0.05*dx, b - 0.05*dy, top + 0.05*dy])

line_simu, = plt.plot([], [],"r--", lw=2, markersize=4 , label = "Some curve" ,  zorder= 1 )
time_text = plt.text(0.1, 0.1, '', zorder=10)
SB_text = plt.text(0.1, 1.9, '', zorder=10)
AB_text = plt.text(0.1, 1.8, '', zorder=10)
AB2_text = plt.text(0.1, 1.7, '', zorder=10)

anim = animation.FuncAnimation(fig, animate,  frames=range(1000), interval=200,blit=True,repeat=False)