#ifndef LINALG_FILE_COCG_CXX

#include "CoCg.hxx"

namespace linalg
{

  //! Applique le preconditionneur z = M r
  template<class T>
  void IdentityPreconditioner<T>::Solve(const Vector<T>& r, Vector<T>& z)
  {
    z = r;
  }

  
  //! Resout le systeme lineaire A x = b avec la methode du gradient conjugue
  /*!
    \param[in] A matrice associee au systeme lineaire a resoudre
    \param[inout] x en entree "initial guess", en sortie la solution
    \param[in] b second membre
    \param[in] prec preconditioneur a utiliser
    \param[in] epsilon critere d'arret
    \param[in] nb_iter_max nombre maximal d'iterations
    Cette version du gradient conjugue (appele COCG) converge pour 
    des matrices reelles symetriques ou complexes symetriques.
    Elle ne fonctionne pas pour des matrices non-symetriques ou complexes hermitiennes
   */
  template <class T>
  void ConjugateGradient(const VirtualMatrix<T>& A, Vector<T>& x, const Vector<T>& b,
			 VirtualPreconditioner<T>& prec, double epsilon, int nb_iter_max)
  {  
    T rho(1), rho_1(1);
    T alpha, beta, delta;
    Vector<T> p(b), q(b), r(b), z(b);
    
    // we compute the initial residual r = b - Ax
    r = b;
    A.MltAdd(-1.0, x, r);
    
    int nb_iter = 0;
    // Loop until the stopping criteria are satisfied
    while ((Norm2(r)/Norm2(b) > epsilon) && (nb_iter < nb_iter_max))
      {
	// Preconditioning z = M^{-1} r
	prec.Solve(r, z);
	
	rho = DotProd(r, z);
      
	if (nb_iter == 0)
	  p = z;
	else
	  {
	    // p = beta*p + z  where  beta = rho_i/rho_{i-1}
	    beta = rho / rho_1;
	    p *= beta;
	    Add(T(1), z, p);
	  }
	
	// matrix vector product q = A*p
	A.Mlt(p, q);
	delta = DotProd(p, q);
	alpha = rho / delta;
	
	// x = x + alpha*p  and r = r - alpha*q  where alpha = rho/(bar(p),q)
	Add(alpha, p, x);
	Add(-alpha, q, r);
	
	rho_1 = rho;
	
	nb_iter++;
	if (nb_iter%10 == 0)
	  cout << "Residu at iteration " << nb_iter << " = " << Norm2(r)/Norm2(b) << endl;
      }
  }

}

#define LINALG_FILE_COCG_CXX
#endif
