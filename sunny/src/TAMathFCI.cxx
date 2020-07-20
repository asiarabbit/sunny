/**
  SUNNY Project, Anyang Normal University, IMP-CAS
  \file TAMathFCI.C
  \class TAMathFCI
  \brief A collection of math functions.
  \author SUN Yazhou, asia.rabbit@163.com
  \date Created: 2020/02/01
  \date Last modified: 2020/02/29
  \copyright 2020, SUN Yazhou
	\copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#include <cstring>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "TAMathFCI.h"
#include "TAException.h"

using std::max_element;
using std::cout;
using std::endl;

double TAMathFCI::sign(double c){
  if(c >= 0) return 1.;
  else return -1.;
}
int TAMathFCI::Binomial(int n, int m){
  if(m > n)
    TAException::Error("TAMathFCI",
      "Binomial: n: %d is larger than m: %d!", n, m);
  return Factorial(n) / (Factorial(m) * Factorial(n-m));
}

/// \return n!
int TAMathFCI::Factorial(int n){
  if(n < 0) TAException::Error("TAMathFCI", "Factorial: n: %d is minus", n);

  return n <= 1 ? 1 : n * Factorial(n-1);
}

/// return the dominant eigenvalue using power method
double TAMathFCI::EigenPower(const TAMatrix2D &ma, TAMatrix2D &v){
  if(!ma.IsSquare())
    TAException::Error("TAMathFCI",
      "EigenPower: Input matrix is not square.");
  if(!v.IsVector()){
    TAException::Error("TAMathFCI", "EigenPower:\
Input vector is not a vector.");
  }
  //in case that the eigenvalue is zero
  TAMatrix2D vzero(v); // then the eigenvector is v itself

  double mv, mv_old = 1E200; // the maximum element in u

  while(fabs(mv - mv_old) > 1E-4){
    static int round = 0; // DEBUG
    v = ma*v;
    mv_old = mv;
    mv = **max_element(v.cv(0).begin(), v.cv(0).end(),
      [](double *a, double *b){ return fabs(*a) < fabs(*b); } );
    if(0. != mv) v /= mv;

    cout << "Round: " << round++ << endl; // DEBUG
    cout << "After renormalization, max_v: " << mv << endl; // DEBUG
    v.Print(); // DEBUG
    getchar(); // DEBUG
  } // end of while

  cout << "The dominant eigenvalue is \033[33;1m" << mv << "\033[0m\n";
  cout << "The corresponding eigenvector is " << endl;
  if(mv != 0) v.Print();
  else vzero.Print();

  return mv;
} // end of member function EigenPower

void TAMathFCI::GramSchmidt(const TAMatrix2D &A, TAMatrix2D &Q){
  const int nc = A.ncol(), nr = A.nrow();
  if(nc > nr){
    TAException::Error("TAMathFCI",
      "GramSchmidt: Vectors to orthogonalize are larger than A's rank.");
  }
  if(Q.ncol() != nc || Q.nrow() != nr) Q.Resize(nr, nc);

  for(int i = 0; i < nc; i++){
    Q.cv(i) = A.cv(i);
    for(int j = 0; j < i; j++){
      Q.cv(i) -= Q.cv(j)*A.cv(i) * Q.cv(j); // subtract qj's component
    } // end for over qj
    Q.cv(i).normalize();
  } // end loop over qi
} // end of member function GramSchmidt

// QR factorization: A=QR, R=Q^(-1)*A=Q^(T)*A
void TAMathFCI::QR(const TAMatrix2D &A, TAMatrix2D &Q, TAMatrix2D &R){
  GramSchmidt(A, Q); // implement Gram-Schmidt orthogonalization
  R = Q.Transpose()*A;
} // end of member function QR

/// solve all the eigenvalues of a matrix using QR method
/// \param v: to store all the eigenvalues
void TAMathFCI::EigenQR(const TAMatrix2D &A, TAMatrix2D &v){
  const int n = A.nrow(), nc = A.ncol();
  if(nc != n){
    TAException::Error("TAMathFCI",
      "EigenQR: Input matrix is not square, nc: %d, nr: %d", nc, n);
  }

  if(!v.IsVector() || v.ncol() < n) v.Resize(n, 1);
  TAMatrix2D Q(n, n), R(n, n), Ak(A);
  TAMatrix2D Qv(n, n); Qv = 1.; // DEBUG
  double epsilon = 1E200; // the sum of elements below the diagonal of Ak
  while(epsilon > 1E-2){
    static int round = 0; // DEBUG
    cout << "Round: " << round++ << " ____________" << endl; // DEBUG
    // the core algorithm //
    QR(Ak, Q, R); // QR factorization: An=Q_n * R_n
    Ak = R*Q; // A_(n+1) = R_n * Q_n
    Qv *= Q;

    // calculate epsilon //
    epsilon = 0.;
    for(int i = 0; i < n; i++)
      for(int j = 0; j < i; j++) epsilon += fabs(Ak[i][j]);

    Q.Print(); // DEBUG
    R.Print(); // DEBUG
    (Q*R).Print(); // DEBUG
    (Q*Q.Transpose()).Print(); // DEBUG
    Ak.Print(); // DEBUG
    Qv.Print(); // DEBUG
    cout << "epsilon: " << epsilon << endl << endl; // DEBUG
    getchar(); // DEBUG
  }
  // output the diagonal elements to v as the obtained eigenvalue set
  for(int i = 0; i < n; i++) v[i][0] = Ak[i][i];

  Ak.Print(); // DEBUG
  v.Print(); // DEBUG
  Qv.Print(); // DEBUG
}

// solve all the eigenvalues (in v) and eigenvectors (in Q) using Jacobi method
void TAMathFCI::EigenJacobi(const TAMatrix2D &A, TAMatrix2D &P, TAMatrix2D &v){
  if(!A.IsSymmetric())
    TAException::Error("TAMathFCI",
      "EigenJacobi: Input matrix is not symmetric");
  const int n = A.ncol();
  if(P.nrow() != n || P.ncol() != n) P.Resize(n, n);
  if(!v.IsVector() || v.ncol() < n) v.Resize(n, 1);

  // A_(k+1) = Rk^T*A_k*Rk; P = R1*R2* ... *R(k-2)*R(k-1)*Rk; P^(-1)AP = D
  TAMatrix2D R(n, n), Ak(A); // R: the rotation matrix
  P = 1.; // initialize to unit matrix
  P.Print(); // DEBUG
  Ak.Print(); // DEBUG
  double mx = 1E20; int p, q; // mx = Ak[p][q] = max(|a_ij|): pivotal element
  while(mx > 1E-3){
    // single out the pivotal element //
    static int round = 0; // DEBUG
    mx = 0;
    double tmp = 0;
    for(int i = 1; i < n; i++) for(int j = 0; j < i; j++)
      if((tmp = fabs(Ak[i][j])) > mx){ p = i; q = j; mx = tmp; }
    cout << "round: " << round++ << endl; // DEBUG
    cout << " mx: " << mx << " A[" << p+1 << "][" << q+1 << "]: "; // DEBUG
    cout << Ak[p][q] << endl; // DEBUG

    // implement the rotation //
    // prepare the rotation matrix -----
    R = 1.; // initialize to unit matrix
    double c = (Ak[p][p] - Ak[q][q]) / (2.*Ak[p][q]); // 1/tan(2*theta)
    double t = sign(c) / (fabs(c) + sqrt(c*c + 1.)); // tan(theta)
    double sqrtt = sqrt(1. + t*t);
    double sinT = t / sqrtt, cosT = 1. / sqrtt;
    R[p][p] = R[q][q] = cosT; R[p][q] = -sinT; R[q][p] = sinT;
    // apply the rotation to Ak
    Ak = R.Transpose()*Ak*R;
    P *= R; // accumulate the rotation action
    cout << "c: " << c << " t: " << t << endl; // DEBUG
    cout << "sinT: " << sinT << " cosT: " << cosT << endl; // DEBUG
    R.Print(); // DEBUG
    P.Print(); // DEBUG
    Ak.Print(); // DEBUG
    getchar(); // DEBUG
  } // end while
  // output the result
  for(int i = 0; i < n; i++) v[i][0] = Ak[i][i];
  Ak.Print(); // DEBUG
  P.Print(); // DEBUG
  v.Print(); // DEBUG
} // end
