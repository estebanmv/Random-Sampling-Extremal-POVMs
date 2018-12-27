#!/usr/bin/env wolframscript 
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
    ,
  "CohPlusTherGamma",
(* Initialize lists {{{*)
vanmuchos = {};
maximalist = {};
ListaCocientes = {};
Optimal = {};
ComplexCoherent = Sqrt[MeanPhotonNumb];
FieldFrequency = (kBoltzmann/hBar)*Log[ComplexCoherent^(-2) + 1]*Temperature; 
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
        A = {};
        vectA = {};
        dimtot = Outcomedim;
        lenghtunpovm = Length[unPOVM];
        ranks = {};
        Module[{n, r}, 
         For[n = 1, n <= lenghtunpovm, n++, 
           rank = MatrixRank[unPOVM[[n]]];
           If[rank > 1, dimtot = dimtot + rank - 1;
            For[i = 1, i <= rank, i++, vectA = {};
             
             For[r = 1, r <= (HilbertDim^2 - 1), r++, 
              AppendTo[vectA, 
                Chop[0.5*
                  Tr[projector[unPOVM,n, i].(GellMann[HilbertDim][[r]])]]];];
             AppendTo[vectA, 1];
             AppendTo[A, vectA];];, If[rank == 1, vectA = {};
              
              For[r = 1, r <= (HilbertDim^2 - 1), r++, 
               AppendTo[vectA, 
                 Chop[0.5*Tr[POVMUnit[unPOVM,n].(GellMann[HilbertDim][[r]])]]];];
              AppendTo[vectA, 1];
              AppendTo[A, vectA];];]];];
        A = Transpose[A];
        b = Table[0, {n, 1, HilbertDim^2 - 1}]~Join~{HilbertDim};
        RandV = RandomReal[{0, 1}, Dimensions[A][[2]]];
        bprime = Table[{b[[i]], 0}, {i, 1, Length[b]}];
        
        (*Linear Programming to find a solution. 
        We find the optimum with a random vector, therefore, 
        it is a random element of the polytope.*)
        (* Print[Abs[1-prob]];
         Print[MatrixRank[A]]; *)
        
        Quiet[Check[Sol = 
          LinearProgramming[RandV, A, bprime, Method -> "Simplex"],
         Break[]]];
        
        (*With the result of the Linear Program we construct the extremal POVM.*)
        ExtPOVM = {};
        ranksum = 0;
        rank1case = 0;
        prevcases = 0;
        Module[{r, u}, 
         For[r = 1, r <= lenghtunpovm, r++, 
           rank = MatrixRank[unPOVM[[r]]];
           
           If[rank > 1, 
            
            For[u = 1, u <= rank, u++, 
             AppendTo[ExtPOVM, 
              projector[unPOVM,r, u]*Sol[[prevcases + u]]]];
            ranksum = ranksum + rank;, 
            If[rank == 1, 
              AppendTo[ExtPOVM, POVMUnit[unPOVM,r]*Sol[[prevcases + 1]]];
              rank1case = rank1case + 1];];
           prevcases = ranksum + rank1case;];];
        Extremal = {};
        For[u = 1, u <= Length[ExtPOVM], u++, 
         If[MatrixRank[Chop[ExtPOVM[[u]]]] != 0, 
          AppendTo[Extremal, ExtPOVM[[u]]]]];
        avector = {};
        auxa = {};
        
        (*We construct the auxiliar POVM and recall it unPOVM.*)
     
           Module[{\[Nu]}, 
         For[\[Nu] = 1, \[Nu] <= lenghtunpovm, \[Nu]++, 
           rank = MatrixRank[unPOVM[[\[Nu]]]];
           
           If[rank > 1, 
            For[\[Eta] = 1, \[Eta] <= rank, \[Eta]++, 
             AppendTo[avector, Chop[eivalues[unPOVM,\[Nu], \[Eta]]]];], 
            If[rank == 1, 
              AppendTo[avector, Chop[Tr[unPOVM[[\[Nu]]]]]]];]];];
        ListaCocientes = 
         Table[If[Sol[[\[Sigma]]] == 0, 0, 
           avector[[\[Sigma]]]/Sol[[\[Sigma]]]], {\[Sigma], 1, 
           Length[avector]}];
        prob = Min[Select[ListaCocientes, # > 0. &]];
        If[Chop[Abs[1. - prob], 10^-5] == 0, Break[], 
         SolAux = (1./(1. - prob))*(avector - prob*Sol);];
        
        otroPOVM = {};
        prevcases = 0;
        ranksum = 0;
        rank1case = 0;
        Module[{k, u}, 
         For[k = 1, k <= lenghtunpovm, k++, 
           rank = MatrixRank[unPOVM[[k]]];
           
           If[rank > 1, 
            For[u = 1, u <= rank, u++, 
             AppendTo[otroPOVM, 
              projector[unPOVM,k, u]*SolAux[[prevcases + u]]]];
            ranksum = ranksum + rank;, 
            If[rank == 1, 
              AppendTo[otroPOVM, 
               POVMUnit[unPOVM,k]*SolAux[[prevcases + 1]]];
              rank1case = rank1case + 1];];
           prevcases = rank1case + ranksum;];];
        unPOVM = otroPOVM;
        
        (*We calculate the Van Trees Information with the POVM found in this iteration.*)
        
   AppendTo[VT, 
         Chop[NIntegrate[
            FisherCoherentPlusTher[Extremal,ThetaPhase, ComplexCoherent,HilbertDim,FieldFrequency,Temperature,MixConstant]*
            Gammapdf[4, 1.5, ThetaPhase], {ThetaPhase, 0, Infinity}], 10^-9] + Sobra];
        
        ranks = Table[
          MatrixRank[Extremal[[m]]], {m, 1, Length[Extremal]}]; 
        If[Max[VT] > VanTrees, Optimal = Extremal; 
         VanTrees = Max[VT];];] AppendTo[maximalist, VanTrees];
      (*Sampling number view.*)
      
      (*Export["./kappa_samplingnum.dat", \[Kappa]];*)
      PrintTemporary[\[Kappa]];
      , {\[Kappa], 
       Samplings}];][[1]] (*>>> 
  "./tiemposCoherentplusthermalExtremales.dat";*)
  (*}}}*)
(*PutAppend[Max[maximalist], WriteDirectory<>"/VanCohplusthervalues.dat"];*) (*This file contains the Maximum value of the Van Trees Information found.*)
  Print["Max{Van Trees} = ",Max[maximalist]];
  ,
  "CohPlusTherGaussian",
(* Initialize lists {{{*)
vanmuchos = {};
maximalist = {};
ListaCocientes = {};
Optimal = {};
ComplexCoherent = Sqrt[MeanPhotonNumb];
FieldFrequency = (kBoltzmann/hBar)*Log[ComplexCoherent^(-2) + 1]*Temperature; 
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
        A = {};
        vectA = {};
        dimtot = Outcomedim;
        lenghtunpovm = Length[unPOVM];
        ranks = {};
        Module[{n, r}, 
         For[n = 1, n <= lenghtunpovm, n++, 
           rank = MatrixRank[unPOVM[[n]]];
           If[rank > 1, dimtot = dimtot + rank - 1;
            For[i = 1, i <= rank, i++, vectA = {};
             
             For[r = 1, r <= (HilbertDim^2 - 1), r++, 
              AppendTo[vectA, 
                Chop[0.5*
                  Tr[projector[unPOVM,n, i].(GellMann[HilbertDim][[r]])]]];];
             AppendTo[vectA, 1];
             AppendTo[A, vectA];];, If[rank == 1, vectA = {};
              
              For[r = 1, r <= (HilbertDim^2 - 1), r++, 
               AppendTo[vectA, 
                 Chop[0.5*Tr[POVMUnit[unPOVM,n].(GellMann[HilbertDim][[r]])]]];];
              AppendTo[vectA, 1];
              AppendTo[A, vectA];];]];];
        A = Transpose[A];
        b = Table[0, {n, 1, HilbertDim^2 - 1}]~Join~{HilbertDim};
        RandV = RandomReal[{0, 1}, Dimensions[A][[2]]];
        bprime = Table[{b[[i]], 0}, {i, 1, Length[b]}];
        
        (*Linear Programming to find a solution. 
        We find the optimum with a random vector, therefore, 
        it is a random element of the polytope.*)
        (* Print[Abs[1-prob]];
         Print[MatrixRank[A]]; *)
        
        Quiet[Check[Sol = 
          LinearProgramming[RandV, A, bprime, Method -> "Simplex"],
         Break[]]];
        
        (*With the result of the Linear Program we construct the extremal POVM.*)
        ExtPOVM = {};
        ranksum = 0;
        rank1case = 0;
        prevcases = 0;
        Module[{r, u}, 
         For[r = 1, r <= lenghtunpovm, r++, 
           rank = MatrixRank[unPOVM[[r]]];
           
           If[rank > 1, 
            
            For[u = 1, u <= rank, u++, 
             AppendTo[ExtPOVM, 
              projector[unPOVM,r, u]*Sol[[prevcases + u]]]];
            ranksum = ranksum + rank;, 
            If[rank == 1, 
              AppendTo[ExtPOVM, POVMUnit[unPOVM,r]*Sol[[prevcases + 1]]];
              rank1case = rank1case + 1];];
           prevcases = ranksum + rank1case;];];
        Extremal = {};
        For[u = 1, u <= Length[ExtPOVM], u++, 
         If[MatrixRank[Chop[ExtPOVM[[u]]]] != 0, 
          AppendTo[Extremal, ExtPOVM[[u]]]]];
        avector = {};
        auxa = {};
        
        (*We construct the auxiliar POVM and recall it unPOVM.*)
     
           Module[{\[Nu]}, 
         For[\[Nu] = 1, \[Nu] <= lenghtunpovm, \[Nu]++, 
           rank = MatrixRank[unPOVM[[\[Nu]]]];
           
           If[rank > 1, 
            For[\[Eta] = 1, \[Eta] <= rank, \[Eta]++, 
             AppendTo[avector, Chop[eivalues[unPOVM,\[Nu], \[Eta]]]];], 
            If[rank == 1, 
              AppendTo[avector, Chop[Tr[unPOVM[[\[Nu]]]]]]];]];];
        ListaCocientes = 
         Table[If[Sol[[\[Sigma]]] == 0, 0, 
           avector[[\[Sigma]]]/Sol[[\[Sigma]]]], {\[Sigma], 1, 
           Length[avector]}];
        prob = Min[Select[ListaCocientes, # > 0. &]];
        If[Chop[Abs[1. - prob], 10^-5] == 0, Break[], 
         SolAux = (1./(1. - prob))*(avector - prob*Sol);];
        
        otroPOVM = {};
        prevcases = 0;
        ranksum = 0;
        rank1case = 0;
        Module[{k, u}, 
         For[k = 1, k <= lenghtunpovm, k++, 
           rank = MatrixRank[unPOVM[[k]]];
           
           If[rank > 1, 
            For[u = 1, u <= rank, u++, 
             AppendTo[otroPOVM, 
              projector[unPOVM,k, u]*SolAux[[prevcases + u]]]];
            ranksum = ranksum + rank;, 
            If[rank == 1, 
              AppendTo[otroPOVM, 
               POVMUnit[unPOVM,k]*SolAux[[prevcases + 1]]];
              rank1case = rank1case + 1];];
           prevcases = rank1case + ranksum;];];
        unPOVM = otroPOVM;
        
        (*We calculate the Van Trees Information with the POVM found in this iteration.*)
        
   AppendTo[VT, 
         Chop[NIntegrate[
            FisherCoherentPlusTher[Extremal,Theta, ComplexCoherent,HilbertDim,FieldFrequency,Temperature,MixConstant]*
            Gaussianpdf[Theta, \[Pi]/4, \[Pi]], {Theta, 0, 2*\[Pi]}], 10^-9] + SobraGaussian];
        
        ranks = Table[
          MatrixRank[Extremal[[m]]], {m, 1, Length[Extremal]}]; 
        If[Max[VT] > VanTrees, Optimal = Extremal; 
         VanTrees = Max[VT];];] AppendTo[maximalist, VanTrees];
      (*Sampling number view.*)
      
      (*Export["./kappa_samplingnum.dat", \[Kappa]];*)
      PrintTemporary[\[Kappa]];
      , {\[Kappa], 
       Samplings}];][[1]] (*>>> 
  "./tiemposCoherentplusthermalExtremales.dat";*)
  (*}}}*)
(*PutAppend[Max[maximalist], WriteDirectory<>"/VanCohplusthervalues.dat"];*) (*This file contains the Maximum value of the Van Trees Information found.*)
  Print["Max{Van Trees} = ",Max[maximalist]];
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
        A = {};
        vectA = {};
        dimtot = Outcomedim;
        lenghtunpovm = Length[unPOVM];
        ranks = {};
        Module[{n, r}, 
         For[n = 1, n <= lenghtunpovm, n++, 
           rank = MatrixRank[unPOVM[[n]]];
           If[rank > 1, dimtot = dimtot + rank - 1;
            For[i = 1, i <= rank, i++, vectA = {};
             
             For[r = 1, r <= (HilbertDim^2 - 1), r++, 
              AppendTo[vectA, 
                Chop[0.5*
                  Tr[projector[unPOVM,n, i].(GellMann[HilbertDim][[r]])]]];];
             AppendTo[vectA, 1];
             AppendTo[A, vectA];];, If[rank == 1, vectA = {};
              
              For[r = 1, r <= (HilbertDim^2 - 1), r++, 
               AppendTo[vectA, 
                 Chop[0.5*Tr[POVMUnit[unPOVM,n].(GellMann[HilbertDim][[r]])]]];];
              AppendTo[vectA, 1];
              AppendTo[A, vectA];];]];];
        A = Transpose[A];
        b = Table[0, {n, 1, HilbertDim^2 - 1}]~Join~{HilbertDim};
        RandV = RandomReal[{0, 1}, Dimensions[A][[2]]];
        bprime = Table[{b[[i]], 0}, {i, 1, Length[b]}];
        
        (*Linear Programming to find a solution. 
        We find the optimum with a random vector, therefore, 
        it is a random element of the polytope.*)
        
        Quiet[Check[Sol = 
          LinearProgramming[RandV, A, bprime, Method -> "Simplex"],
         Break[]]];
        
        (*With the result of the Linear Program we construct the extremal POVM.*)
        ExtPOVM = {};
        ranksum = 0;
        rank1case = 0;
        prevcases = 0;
        Module[{r, u}, 
         For[r = 1, r <= lenghtunpovm, r++, 
           rank = MatrixRank[unPOVM[[r]]];
           
           If[rank > 1, 
            
            For[u = 1, u <= rank, u++, 
             AppendTo[ExtPOVM, 
              projector[unPOVM,r, u]*Sol[[prevcases + u]]]];
            ranksum = ranksum + rank;, 
            If[rank == 1, 
              AppendTo[ExtPOVM, POVMUnit[unPOVM,r]*Sol[[prevcases + 1]]];
              rank1case = rank1case + 1];];
           prevcases = ranksum + rank1case;];];
        Extremal = {};
        For[u = 1, u <= Length[ExtPOVM], u++, 
         If[MatrixRank[Chop[ExtPOVM[[u]]]] != 0, 
          AppendTo[Extremal, ExtPOVM[[u]]]]];
        avector = {};
        auxa = {};
        
        (*We construct the auxiliar POVM and recall it unPOVM.*)
     
           Module[{\[Nu]}, 
         For[\[Nu] = 1, \[Nu] <= lenghtunpovm, \[Nu]++, 
           rank = MatrixRank[unPOVM[[\[Nu]]]];
           
           If[rank > 1, 
            For[\[Eta] = 1, \[Eta] <= rank, \[Eta]++, 
             AppendTo[avector, Chop[eivalues[unPOVM,\[Nu], \[Eta]]]];], 
            If[rank == 1, 
              AppendTo[avector, Chop[Tr[unPOVM[[\[Nu]]]]]]];]];];
        ListaCocientes = 
         Table[If[Sol[[\[Sigma]]] == 0, 0, 
           avector[[\[Sigma]]]/Sol[[\[Sigma]]]], {\[Sigma], 1, 
           Length[avector]}];
        prob = Min[Select[ListaCocientes, # > 0. &]];
        If[Chop[Abs[1. - prob], 10^-9] == 0,Break[], 
         SolAux = (1./(1. - prob))*(avector - prob*Sol);];
        
        otroPOVM = {};
        prevcases = 0;
        ranksum = 0;
        rank1case = 0;
        Module[{k, u}, 
         For[k = 1, k <= lenghtunpovm, k++, 
           rank = MatrixRank[unPOVM[[k]]];
           
           If[rank > 1, 
            For[u = 1, u <= rank, u++, 
             AppendTo[otroPOVM, 
              projector[unPOVM,k, u]*SolAux[[prevcases + u]]]];
            ranksum = ranksum + rank;, 
            If[rank == 1, 
              AppendTo[otroPOVM, 
               POVMUnit[unPOVM,k]*SolAux[[prevcases + 1]]];
              rank1case = rank1case + 1];];
           prevcases = rank1case + ranksum;];];
        unPOVM = otroPOVM;
        
        (*We calculate the Van Trees Information with the POVM found \
in this iteration.*)
        
       AppendTo[VT,Chop[NIntegrate[(FisherQubit[Extremal,ThetaAngle,EtaAngle])/(2\[Pi]),{ThetaAngle,0,2\[Pi]}]]];
        
        ranks = Table[
          MatrixRank[Extremal[[m]]], {m, 1, Length[Extremal]}]; 
        If[Max[VT] > VanTrees, Optimal = Extremal; 
         VanTrees = Max[VT];];] AppendTo[maximalist, VanTrees];
      (*Sampling number view.*)
      
      (*Export["./kappa_samplingnum.dat", \[Kappa]];*)
      , {\[Kappa], 
       Samplings}];][[1]] (*>>> 
  "./tiemposCoherentplusthermalExtremales.dat";*)
  (*}}}*)
(*PutAppend[Max[maximalist], WriteDirectory<>"/VanQubitvalues.dat"];*) (*This file contains the Maximum value of the Van Trees Information found.*)
  Print["Max{Van Trees} = ",Max[maximalist]];
    ,
  _,
    Print["Option \"", option, "\" not found"]
  ]
