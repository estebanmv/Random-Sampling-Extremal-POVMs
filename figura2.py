import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

DAnsatz = pd.read_csv('AnsatzDisplacedThermal.dat', header=None, sep='\t')
DAncpmt = pd.read_csv('AnsatzDisplacedThermalcmpmt.dat',header=None, sep='\t')

Dprueba = pd.read_csv('CC-DTherGaussian/DisplacedThermalSqrt-redo.dat', header=None, sep=',')
Dprueba.columns=['x','y']
Dprueba.x = Dprueba.x.astype(float)
Dprueba.y = Dprueba.y.astype(float)

#DAnsatz.columns=['x','y']
df = DAncpmt.append(DAnsatz,ignore_index=True)
df.columns=['x','y']

plt.figure(figsize= (2.5,2))
#plt.figure(figsize=(8,8))
plt.plot(np.sqrt(Dprueba.x[:9]),Dprueba.y[:9],'ro', markersize=2)
plt.plot(np.sqrt(df.x[:29]),df.y[:29],'b-',linewidth=0.7)
#plt.tick_params(labelsize=7,width=0.5)
plt.tick_params(labelsize=5,width=0.5)
plt.savefig('DisplacedThermal.pdf')
#plt.show()

CohAnsatz = pd.read_csv('coherent_Numerical_ansatz.dat', header=None, sep='\t')
CohAnsatz.columns = ['x', 'y']
CAlg = pd.read_csv('CC-Coherent/ExtremalCoherent-redo.dat', header=None, sep=',')
CAlg.columns = ['x', 'y']
CAlg9 = pd.read_csv('CC-Coherent/ExtremalCoherent-redoDH9.dat', header=None, sep=',')
CAlg9.columns = ['x', 'y']

plt.figure(figsize= (2.5,2))
plt.plot(CohAnsatz.x,CohAnsatz.y, "b-",linewidth=0.7)
plt.plot(np.sqrt(CAlg.x[:8]), CAlg.y[:8], "ro", markersize=2)
plt.plot(np.sqrt(CAlg9.x[:8]), CAlg9.y[:8], "go", markersize=2)
plt.tick_params(labelsize=5,width=0.5)
#plt.xticks([0.,0.4,0.8,1.2])
#plt.yticks([1.5,2.5,3.5,4.5,5.5])
plt.savefig('Coherent.pdf')
#.iloc[:86]
#.iloc[:86]
