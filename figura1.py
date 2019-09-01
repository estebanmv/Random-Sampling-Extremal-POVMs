import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

xZ = np.linspace(0,np.pi,200)
yZ = 1-np.absolute(np.cos(xZ))

#ext1 = pd.read_csv('Qubit-VT1ExtCC.dat', header=None, sep=',')
#ext1.columns=['x','y']
ext10 = pd.read_csv('Qubit-VT10ExtCC.dat', header=None, sep=',')
ext10.columns=['x','y']
int1000 = pd.read_csv('Qubit-VT1000IntCC.dat', header=None, sep=',')
int1000.columns=['x','y']
#naimvt10 = pd.read_csv('QubitNaimark-testsVT10.dat', header=None, sep=',')
#naimvt10.columns=['x','y']

ext1Av = pd.read_csv('Qubit-VT1ExtCCAv.dat', header=None, sep=',')
ext1Av.columns = ['x','y']
nymr10Av = pd.read_csv('QubitNaimark-VT10Av.dat', header=None, sep=',')
nymr10Av.columns = ['x','y']
int1000Av = pd.read_csv('Qubit-VT1000IntCCAv.dat', header=None, sep=',')
int1000Av.columns = ['x','y']

xnymr10 = []
ynymr10 = []

for i in range(1,16):
    ynymr10.append(float(np.mean(nymr10Av.y[(i-1)*10:i*10])))
    xnymr10.append(float(np.mean(nymr10Av.x[(i-1)*10:i*10])))

xext1 = []
yext1 = []

for i in range(1,16):
    yext1.append(float(np.mean(ext1Av.y[(i-1)*10:i*10])))
    xext1.append(float(np.mean(ext1Av.x[(i-1)*10:i*10])))

xint1000 = []
yint1000 = []

for i in range(1,16):
    yint1000.append(float(np.mean(int1000Av.y[(i-1)*10:i*10])))
    xint1000.append(float(np.mean(int1000Av.x[(i-1)*10:i*10])))

plt.figure(figsize= (2.5,2))
plt.plot(xZ,yZ, "b-",linewidth=0.7)
#plt.plot(ext1.x,ext1.y,'co',markersize=3.5)
plt.plot(xext1,yext1,'co',markersize=3.5)
plt.plot(ext10.x,ext10.y,'rD',markersize=3)
#plt.plot(naimvt10.x,naimvt10.y,'gs',markersize=2)
plt.plot(xnymr10,ynymr10,'gs',markersize=2)
#plt.plot(int1000.x,int1000.y,'m^',markersize=2)
plt.plot(xint1000,yint1000,'m^',markersize=2)
plt.tick_params(labelsize=7,width=0.5)
plt.savefig("Qubit1-10ext.pdf")

vtresm = pd.read_csv('Qubit-VTErrtimes.dat',header=None)
tvtresm = pd.read_csv('tiemposQubitRSMErrtimes.dat',header=None)
vtresm.columns = ['x','y']
vtnymr = pd.read_csv('QubitNaimark-VTErrtimes.dat',header=None)
tvtnymr = pd.read_csv('tiemposQubitRSMNeymartimes.dat',header=None)
vtnymr.columns = ['x','y']

plt.figure(figsize= (2.5,2))
x = [1,5,10,15,20,25,30,35,40,45,50] 
plt.loglog(tvtresm,1-vtresm.y,"b^",markersize=2)
plt.tick_params(labelsize=7,width=0.5)
plt.loglog(tvtnymr,1-vtnymr.y,"gs",markersize=2)
plt.savefig('Tiempos-ErrorNymrCCAlg.pdf')
