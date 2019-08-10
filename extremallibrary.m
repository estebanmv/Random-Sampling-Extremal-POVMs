BeginPackage["extremallibrary`"]
Get["./Quantum.m"]

Coherent::usage = "Defines a finite approximation to a coherent state. Requires a complex number as first input and an
integer for the dimension to use."
PhaseNumberOperatorMatrix::usage = "Induces a phase depending in the number of the Fock space. Requires the phase angle and
the dimension of the state it is acting on."
DPhaseNumberOperatorMatrix::usage = "Derivative of the previous function."
PsiTheta::usage = "Gives the number operator exponential applied to the coherent state. Requires the same as Coherent plus the phase."
PsiThetaNormalized::usage = "Normalized version of the previous state."
PsiThetaDM::usage = "Density matrix of the normalized state."
expectationN::usage = "Expectation value of the number of photons of the field in terms of the frequency and temperature."
Thermal::usage = "Thermal state that depends on the dimension that one is approximating to, the field frequency and the temperature."
ThermalMod::usage = "Modifiel Thermal state by the phase operation."
kBoltzmann::usage = "Boltzmann constant."
hBar::usage = "Value of the Planck constant divided by 2Pi."
QubitState::usage = "The state of a Qubit given by two angles in the Bloch sphere: Eta and Theta."
QBDM::usage = "Density matrix correspondent to the Qubit."
DilatedQubit::usage = "Generates a Qubit with an ancilla."
CoherentPlusThermalState::usage = "The convex sum of a Coherent plus a thermal state modified by the phase operator."
GellMann::usage = "Yields a table with all the Gell-Mann Matrices of a given dimension. You must input the dimension. This
function was found in the blog: https://blog.suchideas.com/2016/04/sun-gell-mann-matrices-in-mathematica/ ."
CanonicalProjector::usage = "Canonical basis for the matrix space of dimension dim."
POVM::usage = "Generates a random POVM of dimension HilbertDim with Outcomedim the number of outcomes."
Gammapdf::usage = "We use the Gamma distribution for our calculation."
Sobra::usage = "The Fisher Information correspondent to the a priori distribution."
eivect::usage = "Eigenvector j of unPOVM element i."
eivalues::usage = "Eigenvalue j of unPOVM element i."
projector::usage = "Projector to eigenvector i,j."
POVMUnit::usage = "The normalized version of POVM."
PderivCoherentPlusTher::usage = "Derivative of the probability distribution given the extremal POVM for Coherent plus Thermal case. Requires the Phase, the number of the POVM element and a complex number for the Coherent state."
FisherCoherentPlusTher::usage = "Fisher Information given the extremal POVM for the Coherent plus Thermal case. Requires the Phase and the initial Complex number of the coherent state."
PderivQubit::usage = "Derivative of the probability distribution given the extremal POVM for Qubit case. 
Requires two angles in the Bloch sphere."
FisherQubit::usage = "Fisher Information given the extremal POVM for the Qubit case. Requires two angles."
Gaussianpdf::usage = "The Gaussian distribution probability density function."
NormCnsnt::usage = "Normalization constant for Gaussianpdf."
SobraGaussian::usage = "Fisher information for a priori Gaussian distribution."

LaddUp::usage = "Anhilation operator in Fock basis."
LaddDwn::usage = "Creation operator in Fock basis."
Displacement::usage = "Displacement operator for dimension n and complex Alpha."
DisplacedThermal::usage = "Displaced Thermal state with phase variation."
PderivDispTher::usage = "Derivative of the probability distribution given the extremal POVM for Displaced Thermal case. Requires the Phase, the number of the POVM element and a complex number for the Coherent state."
FisherDispTher::usage = "Fisher Information given the extremal POVM for the Displaced Thermal case. Requires the Phase and the initial Complex number of the coherent state."
AConstruction::usage = "Construct the A matrix and the b vector for the Linear program."
LinearProg::usage = "Using the LinearProgramming routine calculate a random feasible solution of the linear program defined by A and b from AConstruction."
BuildExtremal::usage = "Builds extremal solution from the results of the Linear Program."
CalculateP::usage = "Calculates the probability p that separates the extremal POVM from the auxiliar POVM."
AuxiliarSol::usage = "Calculates the auxiliar (residual) solution by extracting the extremal POVM from the original one."

PderivDilatedQubit::usage = "Pderiv for Naimark Qubit."
FisherDilatedQubit::usage = "Fisher for Naimark Qubit."
UnitaryProjectors::usage = "Projectors from unitary matrix."

Begin["Private`"]      
(*Done {{{*)

kBoltzmann = 1.38*^-23;
hBar = (6.626*^-34)/(2 \[Pi]);

Coherent[ComplexCoherent_, Dimension_] := Exp[-(1/2)*((Norm[ComplexCoherent])^2)] Sum[((ComplexCoherent^(n - 1))/Sqrt[((n - 1)!)])*UnitVector[Dimension + 1, n], {n, 1, Dimension + 1}];
PhaseNumberOperatorMatrix[Theta_, Dimension_] := Table[Chop[KroneckerDelta[i, j]*Exp[I*Theta*(i - 1)]], {i, Dimension + 1}, {j, Dimension + 1}];
PsiTheta[ComplexCoherent_, Theta_, Dimension_] := PhaseNumberOperatorMatrix[Theta, Dimension].Coherent[ComplexCoherent, Dimension];
PsiThetaNormalized[ComplexCoherent_, Theta_, Dimension_] := PsiTheta[ComplexCoherent, Theta, Dimension]/(PsiTheta[Conjugate[ComplexCoherent], -Theta, Dimension].PsiTheta[ComplexCoherent, Theta, Dimension]);
PsiThetaDM[ComplexCoherent_, Theta_, Dimension_] := TensorProduct[PsiTheta[ComplexCoherent, Theta, Dimension],PsiTheta[Conjugate[ComplexCoherent], -Theta, Dimension]]/Tr[TensorProduct[PsiTheta[ComplexCoherent, Theta, Dimension],PsiTheta[Conjugate[ComplexCoherent], -Theta, Dimension]]];
expectationN[FieldFrequency_,  Temperature_] := 1/(Exp[(hBar*FieldFrequency)/(kBoltzmann*Temperature)] - 1);
Thermal[Dimension_, FieldFrequency_, Temperature_] := Sum[TensorProduct[UnitVector[Dimension + 1, n],UnitVector[Dimension + 1, n]]*expectationN[FieldFrequency, Temperature]^n/(1 + expectationN[FieldFrequency, Temperature])^(n + 1), {n, Dimension}];
ThermalMod[Dimension_, FieldFrequency_, Temperature_, Theta_] := (PhaseNumberOperatorMatrix[-Theta, Dimension].Thermal[Dimension, FieldFrequency, Temperature]).PhaseNumberOperatorMatrix[Theta, Dimension];
QubitState[Eta_, Theta_] := Exp[-I*Theta/2]*Cos[Eta/2]*UnitVector[2, 1] + Exp[I*Theta/2]*Sin[Eta/2]*UnitVector[2, 2];
QBDM[Eta_, Theta_] := TensorProduct[QubitState[Eta, Theta],QubitState[Eta, -Theta]];
CoherentPlusThermalState[Dimension_, Theta_, ComplexCoherent_, FieldFrequency_, Temperature_, MixConstant_] := MixConstant*PsiThetaDM[ComplexCoherent, Theta, Dimension] + (1 - MixConstant)*ThermalMod[Dimension, FieldFrequency, Temperature, Theta];
GellMann[Dimension_] := Flatten[Table[(*Symmetric case*)
    SparseArray[{{j, k} -> 1, {k, j} -> 1}, {Dimension, Dimension}], {k, 2, Dimension}, {j,  1, k - 1}], 1]~Join~Flatten[
   Table[(*Antisymmetric case*)
    SparseArray[{{j, k} -> -I, {k, j} -> +I}, {Dimension, Dimension}], {k, 2, Dimension}, {j, 1, k - 1}], 1]~Join~
    Table[(*Diagonal case*)
      Sqrt[2/l/(l + 1)] SparseArray[Table[{j, j} -> 1, {j, 1, l}]~Join~{{l + 1, l + 1} -> -l}, {Dimension, Dimension}], {l, 1, Dimension - 1}];
CanonicalProjector[iterator_,Outcomedim_] := TensorProduct[UnitVector[Outcomedim, iterator],UnitVector[Outcomedim, iterator]];
POVM[Iterator_,HilbertDim_,Outcomedim_,HRandomMatrix_] := If[IntegerQ[Iterator] == True,  
    X = {};
      Y = {};
      For[s = 0, s < HilbertDim, s++,
        Y = {};
        For[j = 0, j < HilbertDim, j++,
          AppendTo[Y, 
            Chop[Sum[(Transpose[
            Conjugate[HRandomMatrix]].(KroneckerProduct[IdentityMatrix[HilbertDim], 
                          CanonicalProjector[Iterator,Outcomedim]]).HRandomMatrix)[[i + s*Outcomedim, i + j*Outcomedim]], {i, 1, Outcomedim}]]]];
        AppendTo[X, Y];]; Return[X/Outcomedim]];
Gammapdf[Alpha_, Beta_, x_] := Exp[-x/Beta]*x^(Alpha - 1)*Beta^(-Alpha)/Gamma[Alpha];
Sobra = NIntegrate[Power[D[Gammapdf[4, 1.5, x], x], 2]/Gammapdf[4, 1.5, x], {x, 0, Infinity}];
eivect[unPOVM_,i_, j_] := Eigenvectors[unPOVM[[i]]][[j]];
eivalues[unPOVM_,i_, j_] := Eigenvalues[unPOVM[[i]]][[j]];
projector[unPOVM_,i_, j_] := TensorProduct[eivect[unPOVM,i, j],Conjugate[eivect[unPOVM,i, j]]];
POVMUnit[unPOVM_,Iterator_] := unPOVM[[Iterator]]/Tr[unPOVM[[Iterator]]];
PderivCoherentPlusTher[Extremal_,ThetaPhase_, Iterator_, ComplexCoherent_,HilbertDim_,FieldFrequency_,Temperature_,MixConstant_] := 
Module[{DerivativeVar},D[Tr[CoherentPlusThermalState[HilbertDim -1, DerivativeVar, ComplexCoherent, FieldFrequency, Temperature, MixConstant].Extremal[[Iterator]]], DerivativeVar] /.DerivativeVar -> ThetaPhase];
FisherCoherentPlusTher[Extremal_,ThetaPhase_, ComplexCoherent_,HilbertDim_,FieldFrequency_,Temperature_,MixConstant_] :=
Sum[((PderivCoherentPlusTher[Extremal,ThetaPhase, Iterator, ComplexCoherent,HilbertDim,FieldFrequency,Temperature,MixConstant])^2)/Tr[CoherentPlusThermalState[HilbertDim -1, ThetaPhase, ComplexCoherent, FieldFrequency, Temperature, MixConstant].Extremal[[Iterator]]], 
{Iterator, 1, Length[Extremal]}];
PderivQubit[Extremal_,ThetaAngle_, Iterator_, EtaAngle_] := Module[{DerivativeVar}, D[Tr[QBDM[EtaAngle, DerivativeVar].Extremal[[Iterator]]], DerivativeVar] /. {DerivativeVar -> ThetaAngle}];
FisherQubit[Extremal_,ThetaAngle_, EtaAngle_] := Sum[((PderivQubit[Extremal,ThetaAngle, Iterator, EtaAngle])^2)/Tr[QBDM[EtaAngle, ThetaAngle].Extremal[[Iterator]]], {Iterator, 1, Length[Extremal]}];
NormCnsnt = 1/NIntegrate[Exp[-(Theta - \[Pi])^2/(2*(\[Pi]/4)^2)], {Theta, 0, 2*\[Pi]}];
Gaussianpdf[Theta_, Sigma_, Gamma_] := NormCnsnt*Exp[-(Theta - Gamma)^2/(2*Sigma^2)];
SobraGaussian = NIntegrate[Power[D[Gaussianpdf[Theta, \[Pi]/4, \[Pi]], Theta], 2]/Gaussianpdf[Theta,\[Pi]/4, \[Pi]], {Theta, 0, 2*\[Pi]}];
DPhaseNumberOperatorMatrix[\[Theta]_,m_]:=Table[Chop[KroneckerDelta[i,j]*(i-1)*I*Exp[I*\[Theta]*(i-1)]],{i,m+1},{j,m+1}];
LaddUp[n_, \[Alpha]_] := Array[\[Alpha]*Sqrt[#1 - 1]*KroneckerDelta[#1 + 1, #2] &, {n, n}];
LaddDwn[n_, \[Alpha]_] := Array[\[Alpha]*Sqrt[#1 - 2]*KroneckerDelta[#1 - 1, #2] &, {n, n}];
Displacement[n_, \[Alpha]_] := N[MatrixExp[LaddDwn[n, \[Alpha]] - LaddUp[n, \[Alpha]\[Conjugate]]]];
DisplacedThermal[m_, \[Theta]_, \[Beta]_, \[Nu]_, T_] := PhaseNumberOperatorMatrix[-\[Theta], m].Displacement[m + 1, \[Beta]].Thermal[m, \[Nu], T].Displacement[m + 1, -\[Beta]].PhaseNumberOperatorMatrix[\[Theta], m];

PderivDispTher[Extremal_,ThetaPhase_, Iterator_, ComplexCoherent_,HilbertDim_,FieldFrequency_,Temperature_] := Module[{DerivativeVar}, D[Tr[DisplacedThermal[HilbertDim, DerivativeVar, ComplexCoherent, FieldFrequency, Temperature].Extremal[[Iterator]]], DerivativeVar] /. {DerivativeVar -> ThetaPhase}];

FisherDispTher[Extremal_,ThetaPhase_, ComplexCoherent_,HilbertDim_,FieldFrequency_,Temperature_] :=Sum[((PderivDispTher[Extremal,ThetaPhase, Iterator, ComplexCoherent,HilbertDim-1,FieldFrequency,Temperature])^2)/Tr[DisplacedThermal[HilbertDim -1, ThetaPhase, ComplexCoherent, FieldFrequency, Temperature].Extremal[[Iterator]]], {Iterator, 1, Length[Extremal]}];

AConstruction[unPOVM_,HilbertDim_,Outcomedim_]:= Do[Module[{n, r, i}, 
        A = {};
        vectA = {};
        dimtot = Outcomedim;
        lenghtunpovm = Length[unPOVM];
        ranks = {};
         For[n = 1, n <= lenghtunpovm, n++, 
           rank = MatrixRank[unPOVM[[n]]];
           If[rank > 1, dimtot = dimtot + rank - 1;
            For[i = 1, i <= rank, i++, vectA = {};
             
             For[r = 1, r <= (HilbertDim^2 - 1), r++, 
              AppendTo[vectA, 
                Chop[
                  Tr[projector[unPOVM,n, i].(GellMann[HilbertDim][[r]])]]];];
             AppendTo[vectA, 1];
             AppendTo[A, vectA];];, If[rank == 1, vectA = {};
              
              For[r = 1, r <= (HilbertDim^2 - 1), r++, 
               AppendTo[vectA, 
                 Chop[Tr[POVMUnit[unPOVM,n].(GellMann[HilbertDim][[r]])]]];];
              AppendTo[vectA, 1];
              AppendTo[A, vectA];];]];
        A = Transpose[A];
        b = Table[0, {n, 1, HilbertDim^2 - 1}]~Join~{HilbertDim};
        ];
        Return[{A,b}],1];

LinearProg[A_,b_]:=Do[
        bprime = Table[{b[[i]], 0}, {i, 1, Length[b]}];
        RandV = RandomReal[{0, 1}, Dimensions[A][[2]]];

        Sol = LinearProgramming[RandV, A, bprime, Method -> "Simplex"];
         Return[Sol],1];

BuildExtremal[Sol_,unPOVM_]:= Do[
        Module[{r, u}, 
       ExtPOVM = {};
        ranksum = 0;
        rank1case = 0;
        prevcases = 0;
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
           prevcases = ranksum + rank1case;];
        Extremal = {};
        For[u = 1, u <= Length[ExtPOVM], u++, 
         If[MatrixRank[Chop[ExtPOVM[[u]]]] != 0, 
          AppendTo[Extremal, ExtPOVM[[u]]]]];
          ];
          Return[Extremal],1];

CalculateP[Sol_,unPOVM_]:= 
        Do[avector = {};
        auxa = {};
           Module[{\[Nu]}, 
         For[\[Nu] = 1, \[Nu] <= lenghtunpovm, \[Nu]++, 
           rank = MatrixRank[unPOVM[[\[Nu]]]];
           
           If[rank > 1, 
            For[\[Eta] = 1, \[Eta] <= rank, \[Eta]++, 
             AppendTo[avector, Chop[eivalues[unPOVM,\[Nu], \[Eta]]]];], 
            If[rank == 1, 
              AppendTo[avector, Chop[Tr[unPOVM[[\[Nu]]]]]]];]];];
        ListaCocientes = {};
              ListaCocientes = 
         Table[If[Sol[[\[Sigma]]] == 0, 0, 
           avector[[\[Sigma]]]/Sol[[\[Sigma]]]], {\[Sigma], 1, 
           Length[avector]}];
        prob = Min[Select[ListaCocientes, # > 0. &]];
        Return[{prob, avector}],1];

AuxiliarSol[prob_, avector_, Sol_, unPOVM_] := 
        Do[If[Chop[Abs[1. - prob], 10^-5] == 0, Break[], 
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
           Return[otroPOVM],1]; 


DilatedQubit[EtaAngle_, ThetaAngle_,RandVDM_] := KroneckerProduct[QBDM[EtaAngle, ThetaAngle], RandVDM];        

PderivDilatedQubit[ThetaAngle_, Iterator_, EtaAngle_, RandVDM_,unPOVM_] := Module[{DerivativeVar}, D[Tr[DilatedQubit[EtaAngle, DerivativeVar,RandVDM].unPOVM[[Iterator]]], DerivativeVar] /. {DerivativeVar -> ThetaAngle}];
FisherDilatedQubit[ThetaAngle_, EtaAngle_,RandVDM_,unPOVM_] := Sum[((PderivDilatedQubit[ThetaAngle, Iterator, EtaAngle,RandVDM,unPOVM])^2)/Tr[DilatedQubit[EtaAngle, ThetaAngle,RandVDM].unPOVM[[Iterator]]], {Iterator, 1, 4}];

UnitaryProjectors[Iterator_,G_] := TensorProduct[G[[All, Iterator]],Conjugate[G][[All, Iterator]]];
         
         (*}}}*)

End[]
EndPackage[]
