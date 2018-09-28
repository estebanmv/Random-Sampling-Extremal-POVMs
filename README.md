# RANDOM SAMPLING OF EXTREMAL POVMS.
This is an implementation of a random sample plus a decomposition into extremal POVMs with an 
algorithm for maximizing the Van Trees Information.  The system at hand is
a mix of a coherent state plus a thermal state, where a phase is induced. We can also do the calculation for the qubit case.
The algorithm can be found in  "Decomposition of any quantum measurement into extremals", G. Sentís et. al.  J. Phys. A: Mat. 
Th. Vol 46  Num 37.
This implementation was done by Esteban Martínez Vargas 2017-18: Esteban.Martinez@uab.cat

Change the path to the Mathematica command, for Linux and Mac systems the path normally is the following:
#!/usr/bin/env wolframscript 

This program requires the Quantum.m package by Carlos Pineda (carlospgmat03), which is included. I leave
a link to his GitHub: 
https://github.com/carlospgmat03/libs

   Usage examples:
   
   ./RandomExtremalPOVMs.wl --Samplings 100 --HilbertDim 6 --Outcomedim 8 --Temperature 1*^-4
   
   ./RandomExtremalPOVMs.wl -o Qubit --EtaAngle 3.1416

The program will print the solution to the screen.
