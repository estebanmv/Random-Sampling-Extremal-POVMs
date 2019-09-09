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
#plt.plot(xnymr10,ynymr10,'gs',markersize=2)
#plt.plot(int1000.x,int1000.y,'m^',markersize=2)
plt.plot(xint1000,yint1000,'m^',markersize=2)
plt.tick_params(labelsize=7,width=0.5)
plt.savefig("Qubit1-10ext.pdf")

vtresm30 = pd.read_csv('Qubit-VTErrtimesAvrg30.dat',header=None)
tvtresm30 = pd.read_csv('tiemposQubitRSMErrtimesAvrg30.dat',header=None)
vtresm30.columns = ['x','y']
vtresm20 = pd.read_csv('Qubit-VTErrtimesAvrg20.dat',header=None)
tvtresm20 = pd.read_csv('tiemposQubitRSMErrtimesAvrg20.dat',header=None)
vtresm20.columns = ['x','y']
vtresm = pd.read_csv('Qubit-VTErrtimesAvrg.dat',header=None)
tvtresm = pd.read_csv('tiemposQubitRSMErrtimesAvrg.dat',header=None)
vtresm.columns = ['x','y']
vtnymr = pd.read_csv('QubitNaimark-VTErrtimesAvrg.dat',header=None)
tvtnymr = pd.read_csv('tiemposQubitRSMNeymartimesAvrg.dat',header=None)
vtnymr.columns = ['x','y']
vtnymr20 = pd.read_csv('QubitNaimark-VTErrtimesAvrg20.dat',header=None)
tvtnymr20 = pd.read_csv('tiemposQubitRSMNeymartimesAvrg20.dat',header=None)
vtnymr20.columns = ['x','y']
vtcc = pd.read_csv('Qubit-VTErrtimesIntCCAvrg.dat',header=None)
tvtcc = pd.read_csv('tiemposQubitCCAvrg.dat',header=None)
vtcc.columns = ['x','y']

xvtresm30 = []
yvtresm30 = []

for i in range(1,10):
    yvtresm30.append((3/6)*(1-float(np.mean(vtresm30.y[(i-1)*30:i*30]))))
    xvtresm30.append((3/6)*(float(np.mean(tvtresm30[(i-1)*30:i*30]))))

xvtresm20 = []
yvtresm20 = []

for i in range(1,10):
    yvtresm20.append((2/6)*(1-float(np.mean(vtresm20.y[(i-1)*20:i*20]))))
    xvtresm20.append((2/6)*(float(np.mean(tvtresm20[(i-1)*20:i*20]))))

xvtresm10 = []
yvtresm10 = []

for i in range(1,10):
    yvtresm10.append((1/6)*(1-float(np.mean(vtresm.y[(i-1)*10:i*10]))))
    xvtresm10.append((1/6)*(float(np.mean(tvtresm[(i-1)*10:i*10]))))

xvtresm = [sum(x) for x in zip(xvtresm20, xvtresm10, xvtresm30)]
yvtresm = [sum(x) for x in zip(yvtresm20, yvtresm10, yvtresm30)]

xvtnymr10 = []
yvtnymr10 = []

for i in range(1,18):
    yvtnymr10.append((1/3)*(1-float(np.mean(vtnymr.y[(i-1)*10:i*10]))))
    xvtnymr10.append((1/3)*(float(np.mean(tvtnymr[(i-1)*10:i*10]))))

xvtnymr20 = []
yvtnymr20 = []

for i in range(1,18):
    yvtnymr20.append((2/3)*(1-float(np.mean(vtnymr20.y[(i-1)*20:i*20]))))
    xvtnymr20.append((2/3)*(float(np.mean(tvtnymr20[(i-1)*20:i*20]))))

xvtnymr = [sum(x) for x in zip(xvtnymr20, xvtnymr10)]
yvtnymr = [sum(x) for x in zip(yvtnymr20, yvtnymr10)]

xvtcc = []
yvtcc = []

for i in range(1,18):
    yvtcc.append(1-float(np.mean(vtcc.y[(i-1)*10:i*10])))
    xvtcc.append(float(np.mean(tvtcc[(i-1)*10:i*10])))

plt.figure(figsize= (2.5,2))
x = [1,5,10,15,20,25,30,35,40,45,50] 
#plt.loglog(tvtresm,1-vtresm.y,"b^",markersize=2)
plt.tick_params(labelsize=7,width=0.5)
#plt.loglog(xvtnymr,yvtnymr,"gs",markersize=2)
plt.loglog(xvtresm,yvtresm,"rd",markersize=2)
plt.loglog(xvtcc,yvtcc,"m^",markersize=2)
#plt.plot(xvtresm,yvtresm,"b^",markersize=2)
#plt.savefig('Tiempos-ErrorNymrCCAlg.pdf')
plt.savefig('Tiempos-ErrorNymrCCAlgAvrg.pdf')
