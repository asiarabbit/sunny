/**
  SUNNY Project, Anyang Normal University, IMP-CAS
  \file TAMathFCI.h
  \class TAMathFCI
  \brief A collection of math functions.
  \author SUN Yazhou, asia.rabbit@163.com
  \date Created: 2020/02/01
  \date Last modified: 2020/03/01
  \copyright 2020, SUN Yazhou
  \copyright SUNNY project, Anyang Normal University, IMP-CAS
*/

#ifndef _TAMathFCI_h_
#define _TAMathFCI_h_

#include "TAMatrix.h"

class TABit;

class TAMathFCI{
public:
  ~TAMathFCI();

  static constexpr double Pi(){ return 3.14159265358979323846; }
  /// golden cut ratio
	static constexpr double Alpha(){ return 0.61803398874989484820458683436565; }
  /// rad per degree
	static constexpr double DEGREE(){ return 0.01745329251994329547; }
	static constexpr double Sqrt3(){ return 1.73205080756887719318; }

  /// sign of a number
  static double sign(double c);
  /// \return calculate combination number
  /// \param n, m: n should not be smaller than m
  static int Binomial(int n, int m);
  /// \return n!
  static int Factorial(int n);

  //////////// linear algebra operations ///////////////////
  /// solve dominant eigenvalue using power method \retval dominant eigenvalue
  /// \param v the initial vector, and would converge to the eigenvector
  static double EigenPower(const TAMatrix2D &ma, TAMatrix2D &v);

  /// Gram-Schmidt orthogonalization of A. A must be nonsingular
  /// \retval Q  the resulting orthogonal matrix
  /// \param norm  whether to normalize the orthogonalized vector
  static void GramSchmidt(const TAMatrix2D &A, TAMatrix2D &Q);
  /// implement QR factorization: A=QR. A must be nonsingular
  /// Q is orthogonal, and R is uppper triangular
  static void QR(const TAMatrix2D &A, TAMatrix2D &Q, TAMatrix2D &R);
  /// solve all the eigenvalues of a matrix using QR method
  /// \param v: to store all the eigenvalues
  static void EigenQR(const TAMatrix2D &A, TAMatrix2D &v);
  /// solve all the eigenvalues and eigenvectors of A using Jacobi method
  /// A must be symmetric
  /// \param P: column vectors represent eigenvector
  /// \param v: stores eigenvalues corresponding to the eigenvectors in P
  static void EigenJacobi(const TAMatrix2D &A, TAMatrix2D &P, TAMatrix2D &v);
};

#endif
