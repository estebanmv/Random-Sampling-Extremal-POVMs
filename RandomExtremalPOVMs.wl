#!/home/estebanmv/Software/Wolfram/Mathematica/12.0/Executables/wolframscript -script
(*#!/usr/bin/env wolframscript*)
(* RANDOM SAMPLING OF EXTREMAL POVMS.

This is an implementation of a random sample plus a decomposition into extremal POVMs with an 
algorithm for maximizing the Van Trees Information.  The system at hand is
a mix of a coherent state plus a thermal state, where a phase is induced. We can also do the calculation for the qubit case.
The algorithm can be found in  "Decomposition of any quantum measurement into extremals", G. Sentís et. al.  J. Phys. A: Mat. 
Th. Vol 46  Num 37.

This implementation was done by Esteban Martínez Vargas 2017-18: Esteban.Martinez@uab.cat

*)

(* Change the path to the Mathematica command, for UNIX-related systems the path normally is the following:
   #!/usr/bin/env wolframscript *)

Needs["extremallibrary`"]

(*It requires the Quantum.m package by Carlos Pineda (carlospgmat03): 
https://github.com/carlospgmat03/libs
*)
Needs["Quantum`"]

(*
   Usage examples:
   ./RandomExtremalPOVMs.wl --Samplings 100 --HilbertDim 6 --Outcomedim 8 --Temperature 1*^-4

   ./RandomExtremalPOVMs.wl -o Qubit --EtaAngle 3.1416

*)

(* NOTE This program gives 2 warnings that I have not been able to eliminate. One is for Integral 
convergence and the other has to do with finding a Linear Program solution. Neither of them affect directly the result. 
*)
Off[NIntegrate::ncvb]
Off[NIntegrate::slwcon]
Off[Infinity::indet]
(*Defaults {{{*)
{option = "CohPlusTherGamma", Samplings = 150, MeanPhotonNumb = 0.5, Temperature = 1*^-3, MixConstant = 0.5, EtaAngle = \[Pi]/2, 
HilbertDim = 7, Outcomedim = 10, WriteDirectory = "."}
(*}}}*)
(*Flags {{{*)
counterCommandLine = 2;
While[counterCommandLine <= Length[$ScriptCommandLine],
 Switch[$ScriptCommandLine[[counterCommandLine]],
  "-n" | "--MeanPhotonNumb" , 
    MeanPhotonNumb = ToExpression[$ScriptCommandLine[[counterCommandLine+1]]]; counterCommandLine += 2,
  "-T" | "--Temperature" , 
    Temperature = ToExpression[$ScriptCommandLine[[counterCommandLine+1]]]; counterCommandLine += 2,
  "-s" | "--Samplings", 
    Samplings = ToExpression[$ScriptCommandLine[[counterCommandLine+1]]]; counterCommandLine += 2,
  "--WriteDirectory", 
    WriteDirectory = $ScriptCommandLine[[counterCommandLine+1]]; counterCommandLine += 2,
  "-o" | "--option", 
    option = $ScriptCommandLine[[counterCommandLine+1]]; counterCommandLine += 2,
  "-h" | "--EtaAngle", 
    EtaAngle = ToExpression[$ScriptCommandLine[[counterCommandLine+1]]]; counterCommandLine += 2,
  "-hD" | "--HilbertDim", 
    HilbertDim = ToExpression[$ScriptCommandLine[[counterCommandLine+1]]]; counterCommandLine += 2,
  "-od" | "--Outcomedim", 
    Outcomedim = ToExpression[$ScriptCommandLine[[counterCommandLine+1]]]; counterCommandLine += 2,
  "-mix" | "--MixConstant", 
    MixConstant = ToExpression[$ScriptCommandLine[[counterCommandLine+1]]]; counterCommandLine += 2,
  _, 
    Print["Error in parameter line. Option \"", 
      $ScriptCommandLine[[counterCommandLine]] , "\" not found."];  Exit[1]
 ]
];
(*}}}*)
Switch[option,
  "test",
    Print["Testing the option mode"];
    ComplexCoherent = Sqrt[MeanPhotonNumb];
    FieldFrequency = (kBoltzmann/hBar)*Log[Abs[ComplexCoherent]^(-2) + 1]*Temperature; 
    Print[LaddUp[5,1]];
    (*Print[Displacement[5,1+I]];
    Print[DisplacedThermal[5,\[Pi],1,FieldFrequency,Temperature];*)
    ,
  "CohPlusTherGamma",
(* Initialize lists {{{*)
vanmuchos = {};
maximalist = {};
ListaCocientes = {};
Optimal = {};
ComplexCoherent = Sqrt[MeanPhotonNumb];
FieldFrequency = (kBoltzmann/hBar)*Log[Abs[ComplexCoherent]^(-2) + 1]*Temperature; 
(*}}}*)
(*{{{*)   
   Timing[Do[
      (*Initialize several variables.*)
      VT = {};
      prob = 0;
      Sol = {};
      VanTrees = 0;
      H = CUEMember[Outcomedim*(HilbertDim)];
      unPOVM = Table[POVM[g,HilbertDim,Outcomedim], {g, 1, Outcomedim}];
      
      (*The algorithm runs until the probability of obtaining the current solution is almost 1.*)
     
     While[prob < 1.,
    
        (*Decomposition of the original randomly produced POVM into the matrix A.*)
        
       {A,b} = AConstruction[unPOVM,HilbertDim,Outcomedim];
       
        (*Linear Programming to find a solution. We find the optimum with a random vector, therefore, 
        it is a random element of the polytope. We need to add a Break because sometimes a solution 
        cannot be found.*)
        
       Quiet[ Check[ Sol = LinearProg[A,b]; ,Break[]]];
        
        (*With the result of the Linear Program we construct the extremal POVM.*)
       
       Extremal = BuildExtremal[Sol,unPOVM];
        
        (*We calculate the constant probability p.*)

       {prob,avector} = CalculateP[Sol,unPOVM];

        (*We construct the auxiliar POVM and recall it unPOVM.*)
        
       otroPOVM = AuxiliarSol[prob, avector, Sol, unPOVM]; 

        unPOVM = otroPOVM;
        
        (*We calculate the Van Trees Information with the POVM found in this iteration.*)

   AppendTo[VT, Chop[NIntegrate[
            FisherCoherentPlusTher[Extremal,ThetaPhase, ComplexCoherent,HilbertDim,FieldFrequency,Temperature,MixConstant]*
            Gammapdf[4, 1.5, ThetaPhase], {ThetaPhase, 0, Infinity}], 10^-9] + Sobra];
   
   (*If the Van Trees Information obtained in this iteration is larger, we pick that one.*)

        ranks = Table[MatrixRank[Extremal[[m]]], {m, 1, Length[Extremal]}]; 
        If[Max[VT] > VanTrees, Optimal = Extremal; 
         VanTrees = Max[VT];];];
         
         AppendTo[maximalist, VanTrees];
      
      (*Sampling number view.*)
      
      PrintTemporary[\[Kappa]];
  
  , {\[Kappa], 
       Samplings}];][[1]] (*>>> "./tiemposCoherentplusthermalExtremales.dat";*) (*For checking the computation times.*)
  (*}}}*)
  
  (*At the end print the maximum value obtained.*)
  Print["Max{Van Trees} = ",Max[maximalist]];
  ,
  "CohPlusTherGaussian",
(* Initialize lists {{{*)
vanmuchos = {};
maximalist = {};
ListaCocientes = {};
Optimal = {};
ComplexCoherent = Sqrt[MeanPhotonNumb];
FieldFrequency = (kBoltzmann/hBar)*Log[Abs[ComplexCoherent]^(-2) + 1]*Temperature; 
(*}}}*)
(*{{{*)   
   Timing[Do[
      (*Initialize several variables.*)
      VT = {};
      prob = 0;
      Sol = {};
      VanTrees = 0;
      H = CUEMember[Outcomedim*(HilbertDim)];
      unPOVM = Table[POVM[g,HilbertDim,Outcomedim], {g, 1, Outcomedim}];
      
      (*The algorithm runs until the probability of obtaining the current solution is almost 1.*)
     
     While[prob < 1.,
    
        (*Decomposition of the original randomly produced POVM into the matrix A.*)
        
       {A,b} = AConstruction[unPOVM,HilbertDim,Outcomedim];
       
        (*Linear Programming to find a solution. We find the optimum with a random vector, therefore, 
        it is a random element of the polytope. We need to add a Break because sometimes a solution 
        cannot be found.*)
        
       Quiet[ Check[ Sol = LinearProg[A,b]; ,Break[]]];
        
        (*With the result of the Linear Program we construct the extremal POVM.*)
       
       Extremal = BuildExtremal[Sol,unPOVM];
        
        (*We calculate the constant probability p.*)

       {prob,avector} = CalculateP[Sol,unPOVM];

        (*We construct the auxiliar POVM and recall it unPOVM.*)
        
       otroPOVM = AuxiliarSol[prob, avector, Sol, unPOVM]; 

        unPOVM = otroPOVM;
        
        (*We calculate the Van Trees Information with the POVM found in this iteration.*)

   AppendTo[VT, 
         Chop[NIntegrate[
            FisherCoherentPlusTher[Extremal,Theta, ComplexCoherent,HilbertDim,FieldFrequency,Temperature,MixConstant]*
            Gaussianpdf[Theta, \[Pi]/4, \[Pi]], {Theta, 0, 2*\[Pi]}], 10^-9] + SobraGaussian];
        
        (*If the Van Trees Information obtained in this iteration is larger, we pick that one.*)

        ranks = Table[MatrixRank[Extremal[[m]]], {m, 1, Length[Extremal]}]; 
        If[Max[VT] > VanTrees, Optimal = Extremal; 
         VanTrees = Max[VT];];];
         
         AppendTo[maximalist, VanTrees];
      
      (*Sampling number view.*)
      
      PrintTemporary[\[Kappa]];
  
  , {\[Kappa], 
       Samplings}];][[1]] (*>>> "./tiemposCoherentplusthermalExtremales.dat";*) (*For checking the computation times.*)
   
  (*}}}*)
  
  (*At the end print the maximum value obtained.*)
  Print["Max{Van Trees} = ",Max[maximalist]];
  ,
  "DispTherGaussian",
(* Initialize lists {{{*)
vanmuchos = {};
maximalist = {};
ListaCocientes = {};
Optimal = {};
ComplexCoherent = Sqrt[MeanPhotonNumb];
FieldFrequency = (kBoltzmann/hBar)*Log[Abs[ComplexCoherent]^(-2) + 1]*Temperature; 
(*}}}*)
(*{{{*)   
   Timing[Do[
      (*Initialize several variables.*)
      VT = {};
      prob = 0;
      Sol = {};
      VanTrees = 0;
      H = CUEMember[Outcomedim*(HilbertDim)];
      unPOVM = Table[POVM[g,HilbertDim,Outcomedim], {g, 1, Outcomedim}];
      
      (*The algorithm runs until the probability of obtaining the current solution is almost 1.*)
     
     While[prob < 1.,
    
        (*Decomposition of the original randomly produced POVM into the matrix A.*)
        
       {A,b} = AConstruction[unPOVM,HilbertDim,Outcomedim];
       
        (*Linear Programming to find a solution. We find the optimum with a random vector, therefore, 
        it is a random element of the polytope. We need to add a Break because sometimes a solution 
        cannot be found.*)
        
       Quiet[ Check[ Sol = LinearProg[A,b]; ,Break[]]];
        
        (*With the result of the Linear Program we construct the extremal POVM.*)
       
       Extremal = BuildExtremal[Sol,unPOVM];
        
        (*We calculate the constant probability p.*)

       {prob,avector} = CalculateP[Sol,unPOVM];

        (*We construct the auxiliar POVM and recall it unPOVM.*)
        
       otroPOVM = AuxiliarSol[prob, avector, Sol, unPOVM]; 

        unPOVM = otroPOVM;
        
        (*We calculate the Van Trees Information with the POVM found in this iteration.*)

   AppendTo[VT, Chop[NIntegrate[
            FisherDispTher[Extremal,Theta, ComplexCoherent*(1+I)/Sqrt[2],HilbertDim,FieldFrequency,Temperature]*
            Gaussianpdf[Theta, \[Pi]/4, \[Pi]], {Theta, 0, 2*\[Pi]}], 10^-9] + SobraGaussian];
        
        (*If the Van Trees Information obtained in this iteration is larger, we pick that one.*)

        ranks = Table[MatrixRank[Extremal[[m]]], {m, 1, Length[Extremal]}]; 
        If[Max[VT] > VanTrees, Optimal = Extremal; 
         VanTrees = Max[VT];];];
         
         AppendTo[maximalist, VanTrees];
      
      (*Sampling number view.*)
      
      PrintTemporary[\[Kappa],"  ",Max[maximalist]];
  
  , {\[Kappa], 
       Samplings}];][[1]] (*>>> "./tiemposCoherentplusthermalExtremales.dat";*) (*For checking the computation times.*)
  (*}}}*)
  
    Print["For Complex^2 = ",MeanPhotonNumb,"Max{Van Trees} = ",Max[maximalist]];

    {MeanPhotonNumb,Max[maximalist]} >>> "./DisplacedThermalSqrt-tests.dat";
  ,
  "Qubit",
(* Initialize lists {{{*)
vanmuchos = {};
maximalist = {};
ListaCocientes = {};
Optimal = {};
(*Here, the Dimensions have to be fixed. *)
HilbertDim = 2;
Outcomedim = 4;
(*Quiet[Infinity::indet];*)
(*}}}*) 
(*{{{*)   
Timing[Do[VT = {};
      prob = 0;
      Sol = {};
      VanTrees = 0;
      H = CUEMember[Outcomedim*(HilbertDim)];
      unPOVM = Table[POVM[g,HilbertDim,Outcomedim], {g, 1, Outcomedim}];
      (*The algorithm runs until the probability of obtaining the current solution is almost 1.*)
        While[prob < 1.,
        
        (*Decomposition of the original randomly produced POVM into the matrix A.*)
        
       {A,b} = AConstruction[unPOVM,HilbertDim,Outcomedim];
       
        (*Linear Programming to find a solution. We find the optimum with a random vector, therefore, 
        it is a random element of the polytope. We need to add a Break because sometimes a solution 
        cannot be found.*)
        
       Quiet[ Check[ Sol = LinearProg[A,b]; ,Break[]]];
        
        (*With the result of the Linear Program we construct the extremal POVM.*)
       
       Extremal = BuildExtremal[Sol,unPOVM];
        
        (*We calculate the constant probability p.*)

       {prob,avector} = CalculateP[Sol,unPOVM];

        (*We construct the auxiliar POVM and recall it unPOVM.*)
        
       otroPOVM = AuxiliarSol[prob, avector, Sol, unPOVM]; 

        unPOVM = otroPOVM;
        
        (*We calculate the Van Trees Information with the POVM found in this iteration.*)

       AppendTo[VT,Chop[NIntegrate[(FisherQubit[Extremal,ThetaAngle,EtaAngle])/(2\[Pi]),{ThetaAngle,0,2\[Pi]}]]];
        
        ranks = Table[
          MatrixRank[Extremal[[m]]], {m, 1, Length[Extremal]}]; 
        If[Max[VT] > VanTrees, Optimal = Extremal; 
         VanTrees = Max[VT];]; 
         ]; 
         
         AppendTo[maximalist, VanTrees];
      
      , {\[Kappa], 
       Samplings}];][[1]] (*>>> "./tiemposCoherentplusthermalExtremales.dat";*)
  (*}}}*)
  
  Print["Max{Van Trees} = ",Max[maximalist]];
    ,
  _,
    Print["Option \"", option, "\" not found"]
  ]
