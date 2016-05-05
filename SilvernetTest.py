#Preambule
#   TYPE %matplotlib qt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#Domain and parameters
#fig=plt.figure(figsize=(10,6), dpi=100)
fig = plt.figure(num = 0, figsize = (12, 8))#, dpi = 100)
#plt.suptitle("Silvernet", fontsize=12)
#ax=plt.subplots()
#ax = plt.subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)
#ax04 = ax03.twinx()
xvec=[0,1,2]
yvec=[0,1,2]

plt.xlim([0,2])
plt.ylim([0,2])
plt.xlabel('Zonal distance',fontsize=13)
plt.ylabel('Meridional distance',fontsize=13)
plt.title('Model Horizontal Field')

#Constants
albedo=[0.1,0.1,0.1,0.3,0.3,0.3,0.3,0.3,0.3]
S=1366.-200
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
lat=45 # in Summer!

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
U=[0,0,0,0,0,0,0,0,0]
V=[0,0,0,0,0,0,0,0,0]
W=[0,0,0,0,0,0,0,0,0]
Ts=[280,280,280,283,280,280,280,290,280]
Ta=[296,295,290,290,291,290,293,290,294]
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
    
def Qg(Tatm,Tsurf): #GHF
    return -10.
        
def Qe(a,Tatm,Tsurf,time): #LHF: Priest-Taylor
    es=611.2*np.exp(17.67*((Tatm-273.)/(Tatm-273.+243.5)))
    #es=611.*np.exp(17.2694*((Tatm-273.16)/(Tatm-35.86)))
    scc=eps*Lv*es/(Rd*Tatm**2)
    #scc=2*10**(-11)*np.exp(0.0597*Tatm)
    return apt*scc*(-(Qs(a,time)+Ql(Tatm,Tsurf))+Qg(Tatm,Tsurf))/(scc+gamma)
    
def Qh(Tatm,Tsurf): #SHF
    return rho*Cp*(Tatm-Tsurf)/rh
    
def EBAL(a,Tatm,Tsurf,time):
    #return Qs(a,time)+Ql(Tatm,Tsurf)+Qe(a,Tatm,Tsurf,time)+Qg(Tatm,Tsurf)+Qh(Tatm,Tsurf)
    return 10.

def ABAL(a,Tatm,Tsurf,time,Tatm2):
    #return emis_a*(emis_s*sigma*Tsurf**4.-sigma*Tatm**4.-Qh(Tatm,Tsurf)-Qe(a,Tatm,Tsurf,time)+0.5*emis_a*sigma*Tatm2**4.)
    return 10.
   
def ABAL2(a,Tatm,Tsurf,time,Tatm2):
    #return emis_a*(emis_s*sigma*Tsurf**4.-sigma*Tatm**4.-Qh(Tatm,Tsurf)-Qe(a,Tatm,Tsurf,time)+0.5*emis_a*sigma*Tatm2**4.)
    return 10.
    #return emis_a*(0.5*emis_a*sigma*Tatm**4.-sigma*Tatm2**4.-Qh(Tatm,Tsurf)*(1-emis_a)+emis_s*sigma*Tsurf**4.*(1-emis_a)-Qe(a,Tatm,Tsurf,time)*(1-emis_a))
   
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

#Variables
TA=Ta
TA2=Ta2
TS=Ts
Evec=np.zeros(9)
Avec=np.zeros(9)
xvec,yvec=np.meshgrid(xvec,yvec)

#Radiation plot

#Timeloop
def plotter(vec):
    return np.array(np.transpose([vec[0:3],vec[3:6],vec[6:9]]))
    
def init():
    imobj.set_data(Ta)
    time_text.set_text('time = 0.0')

    return imobj , time_text
    
def animate(self):
    global data
    global t
    global x
    
    t=np.append(t,x)

    x += 1.0    
    
    
    time_text.set_text('time = %.1f' % self )
    new=plotter(Ts)
    imobj.set_data(new)
       
    imobj.set_zorder(0)
    
    #ax1=plt.axhline(xmin=0,xmax=1,y=0,linewidth=3,color='k')
    #ax2=plt.axvline(ymin=0,ymax=1,x=0,linewidth=3,color='k')
    #ax3=plt.axhline(xmin=0,xmax=1,y=2,linewidth=3,color='k')
    #ax4=plt.axvline(ymin=0,ymax=1,x=2,linewidth=3,color='k')
    #ax5=plt.axvline(ymin=0,ymax=1,x=0.5,linewidth=3,color='k',ls='--')
    #hoi=plt.contour([0,1,2],[0,1,2],plotter(Ta),15,animated=True)
        
    return imobj , time_text#,ax1,ax2,ax3,ax4,ax5#,SB_text,AB_text,AB2_text

#def forceAspect(ax,aspect=1):
#    im = ax.get_images()
#    extent =  im[0].get_extent()
#    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

data = list(Ts)
#imobj = plt.contourf(xvec,yvec,plotter(data),animated=True,aspect=0.5,extent=[0, 2.0, 0.0, 2.0], alpha=1.0, zorder=1)
imobj=plt.imshow(plotter(data), cmap='jet', animated=True,aspect=0.5,extent=[0, 2.0, 0.0, 2.0], alpha=1.0, zorder=1)

plt.colorbar(orientation='horizontal')
plt.clim(240,310)

line_simu, = plt.plot([], [],"r--", lw=2, markersize=4 , label = "Some curve" ,  zorder= 1 )
time_text = plt.text(0.1, 0.1, '', zorder=10)
#SB_text = plt.text(0.1, 1.9, '', zorder=10)
#AB_text = plt.text(0.1, 1.8, '', zorder=10)
#AB2_text = plt.text(0.1, 1.7, '', zorder=10)

anim = animation.FuncAnimation(fig, animate,  frames=range(10), interval=50,blit=True,repeat=False)