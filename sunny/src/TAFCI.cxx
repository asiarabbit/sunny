/**
  SUNNY project, Anyang Normal University, IMP-CAS
  \file TAFCI.cxx
  \class TAFCI
  \brief Global class for full configuration interaction calculations, including
  but not limited to read-in of user-input single-particle states, generation of
  many-body basis, establishment and diagonalization of the Hamiltonian matrix.
  \author SUN Yazhou, aisa.rabbit@163.com
  \date Created: 2020/02/12
  \date Last modified: 2020/02/12 by SUN Yazhou
  \copyright 2020 SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include "TAFCI.h"
#include "TAHamiltonian.h"
#include "TAMathFCI.h"

TAFCI *TAFCI::kInstance = nullptr;

TAFCI::TAFCI() : fHamiltonian(nullptr){}

TAFCI::~TAFCI(){}

TAFCI *TAFCI::Instance(){
  if(!kInstance) kInstance = new TAFCI();
  return kInstance;
} // end of member function Instance

void TAFCI::Go(){
  if(!fHamiltonian) fHamiltonian = TAHamiltonian::Instance();

  fHamiltonian->InitializeCoefficient(); // DEBUG
  TAMatrix2D &H = fHamiltonian->Matrix();
//  for(int i = H.nrow(); i--;) H[i][i] -= 16.;
  H.Print();

  const int n = H.ncol();
  TAMatrix2D v(n);
  v = {1,3,1};
  //for(int i = n; i--;) v[i][0] = 1.;

  // solve the dominant eigenvalue and eigenvector using power method //
//  TAMathFCI::EigenPower(H, v);

  // solve all of the eigenvalues using QR method //
  // firstly, let's implement QR factorization for H
  //TAMatrix2D Q(n, n), R(n, n);
  //TAMathFCI::QR(H, Q, R);

  // QR method to calculate all the eigenvalues of a matrix
  TAMatrix2D vv(n);
  TAMatrix2D A(n, n);
  A = {
    0, 2, 0,
    2, 1, 2,
    0, 2, 1
  };
//  TAMathFCI::EigenQR(H, vv);
//  A.Print(); // DEBUG
//  TAMathFCI::EigenQR(A, vv);
  // Jacobi method to calculate all the eigenvalues and eigenvectors of a matrix
  TAMatrix2D P(n, n);
  TAMathFCI::EigenJacobi(A, P, vv);
} // end of member function Go
