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

  -T  The temperature.

  -s  The number of random samplings.

  -o  Option of the problem to calculate: CohPlusTher or Qubit.

  -h  The angle eta in the case of choosing the Qubit case.

  -hD The dimension of the Hilbert space to use. As n grows a larger Hilbert space is needed.

  -od The number of outcomes of the POVM.

There are several defaults for the options.

To reproduce the results in the article:
For the Qubit case, put the flag -o Qubit, and vary flag --EtaAngle from 0 to Pi.
For the Phase Estimation case, put the flag -o CohPlusTher that is default. Also, put
the temperature -T 10^-3, the mixing constant --MixConstant 0.5, the samplings -s 80,
the --HilbertDim 7 and the -- Outcomedim 10. For the pure state the MixingConstant is 1.
The squared norm of alpha is given by --MeanPhotonNumb, wich can be varied to reproduce 
the graph. The graphs are given with the norm of alpha.
