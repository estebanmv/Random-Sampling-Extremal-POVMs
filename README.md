# RANDOM SAMPLING OF EXTREMAL POVMS.
This is an implementation of a random sample plus a decomposition into extremal POVMs with an 
algorithm for maximizing the Van Trees Information.  The system at hand is
a mix of a coherent state plus a thermal state, where a phase is induced. We can also do the calculation for the qubit case.
The algorithm can be found in  "Decomposition of any quantum measurement into extremals", G. Sentís et. al.  J. Phys. A: Mat. 
Th. Vol 46  Num 37. This implementation was done by Esteban Martínez Vargas 2017-18: Esteban.Martinez@uab.cat

Change the path to the Mathematica command, for Linux and Mac systems the path normally is the following:

#!/usr/bin/env wolframscript 

This program requires the Quantum.m package by Carlos Pineda (carlospgmat03), which is included. I leave
a link to his GitHub: 
https://github.com/carlospgmat03/libs

   Usage examples:
   
   ./RandomExtremalPOVMs.wl --Samplings 100 --HilbertDim 6 --Outcomedim 8 --Temperature 1*^-4
   
   ./RandomExtremalPOVMs.wl -o Qubit --EtaAngle 3.1416

The program will print the solution to the screen.

The options are:

  -n  The average number of photons of the field.

  -mix  The weight of the pure state (-mix 1 produces a pure state, -mix 0 the maximally mixed state)

  -T  The temperature.

  -s  The number of random samplings.

  -o  Option of the problem to calculate: CohPlusTherGamma, CohPlusTherGaussian or Qubit.
  This includes the Cooherent and Thermal mix using a Gamma distribution, using a Gaussian distribution
  and the Qubit case.

  -h  The angle eta in the case of choosing the Qubit case.

  -hD The dimension of the Hilbert space to use. As n grows a larger Hilbert space is needed.

  -od The number of outcomes of the POVM.

There are several defaults for the options which can be seen in the main running file.
