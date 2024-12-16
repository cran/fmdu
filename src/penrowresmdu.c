//
// Copyright (c) 2020 Frank M.T.A. Busing (e-mail: busing at fsw dot leidenuniv dot nl)
// FreeBSD or 2-Clause BSD or BSD-2 License applies, see Http://www.freebsd.org/copyright/freebsd-license.html
// This is a permissive non-copyleft free software license that is compatible with the GNU GPL. 
//

#include "fmdu.h"

double penrowresmdu( const size_t n, const size_t m, double** delta, const size_t p, const size_t h, double** q, double** b, double** y, int** fy, double** d, const double rlambda, const double llambda, const double glambda, const size_t MAXITER, const double FCRIT, size_t* lastiter, double* lastdif, const bool echo )
// Function penrowresmdu() performs penalized row restricted multidimensional unfolding.
{
  const double EPS = DBL_EPSILON;                                          // 2.2204460492503131e-16
  const double TOL = sqrt( EPS );                                          // 1.4901161193847656e-08
  const double CRIT = sqrt( TOL );                                         // 0.00012207031250000000
  const double TINY = pow( 10.0, ( log10( EPS ) + log10( TOL ) ) / 2.0 );  // 1.8189894035458617e-12
  const double INVTINY = 1.0 / TINY;

  // allocate memory
  double** x = getmatrix( n, p, 0.0 );
  double** imb = getmatrix( n, m, 0.0 );
  double** xtilde = getmatrix( n, p, 0.0 );
  double** ytilde = getmatrix( m, p, 0.0 );
  double** qtrq = getmatrix( h, h, 0.0 );
  double** hhh = getmatrix( h, h, 0.0 );
  double** hhm = getmatrix( h, m, 0.0 );
  double** hhp = getmatrix( h, p, 0.0 );
  double** hmp = getmatrix( m, p, 0.0 );
  double* hh = getvector( h, 0.0 );

  // initialization
  double wr = ( double ) ( m );
  double wc = ( double ) ( n );
  for ( size_t i = 1; i <= h; i++ ) {
    for ( size_t j = 1; j <= h; j++ ) {
      double work = 0.0;
      for ( size_t k = 1; k <= n; k++ ) work += q[k][i] * wr * q[k][j];
      hhh[i][j] = qtrq[i][j] = work;
    }
    hhh[i][i] = qtrq[i][i] += rlambda;
  }

  for ( size_t k = 1; k <= h; k++ ) {
    double work = 0.0;
    for ( size_t i = 1; i <= n; i++ ) work += q[i][k];
    for ( size_t j = 1; j <= m; j++ ) hhm[k][j] = work; 
  }
  int nfy = 0;
  for ( size_t j = 1; j <= m; j++ ) for ( size_t k = 1; k <= p; k++ ) nfy += fy[j][k];

  // update distances and calculate normalized stress
  dgemm( false, false, n, p, h, 1.0, q, b, 0.0, x );
  euclidean2( n, p, x, m, y, d );
  double fridge = 0.0;
  double flasso = 0.0;
  double fgroup = 0.0;
  for ( size_t i = 1; i <= h; i++ ) for ( size_t j = 1; j <= p; j++ ) fridge += b[i][j] * b[i][j];
  for ( size_t i = 1; i <= h; i++ ) for ( size_t j = 1; j <= p; j++ ) flasso += fabs( b[i][j] );
  for ( size_t i = 1; i <= h; i++ ) {
    double work = 0.0;
    for ( size_t j = 1; j <= p; j++ ) work += b[i][j] * b[i][j];
    fgroup += sqrt( work );
  }
  double fold = rlambda * fridge + llambda * flasso + glambda * fgroup;
  for ( size_t i = 1; i <= n; i++ ) {
    for ( size_t j = 1; j <= m; j++ ) {
      const double work = delta[i][j] - d[i][j];
      fold += work * work;
    }
  }
  double fnew = 0.0;

  // echo intermediate results
  if ( echo == true ) echoprogress( 0, fold, fold, fold );

  // start main loop
  size_t iter = 0;
  for ( iter = 1; iter <= MAXITER; iter++ ) {

    // compute B matrix
    for ( size_t i = 1; i <= n; i++ ) {
      for ( size_t j = 1; j <= m; j++ ) imb[i][j] = ( d[i][j] < TINY ? 0.0 : delta[i][j] / d[i][j] );
    }

    // compute preliminary updates: xtilde and xtilde
    for ( size_t i = 1; i <= n; i++ ) {
      double rsb = 0.0;
      for ( size_t k = 1; k <= m; k++ ) rsb += imb[i][k];
      for ( size_t j = 1; j <= p; j++ ) {
        double work = 0.0;
        for ( size_t k = 1; k <= m; k++ ) work += imb[i][k] * y[k][j];
        xtilde[i][j] = rsb * x[i][j] - work;
      }
    }
    for ( size_t i = 1; i <= m; i++ ) {
      double csb = 0.0;
      for ( size_t k = 1; k <= n; k++ ) csb += imb[k][i];
      for ( size_t j = 1; j <= p; j++ ) {
        double work = 0.0;
        for ( size_t k = 1; k <= n; k++ ) work += imb[k][i] * x[k][j];
        ytilde[i][j] = csb * y[i][j] - work;
      }
    }

    // update b
    for ( size_t i = 1; i <= h; i++ ) {
      double work = 0.0;
      for ( size_t j = 1; j <= p; j++ ) work += b[i][j] * b[i][j];
      work = sqrt( work );
      hh[i] = 0.5 * glambda * ( work < TINY ? INVTINY : 1.0 / work );
    }
    dgemm( false, false, h, p, m, 1.0, hhm, y, 0.0, hhp );
    dgemm( true, false, h, p, n, 1.0, q, xtilde, 1.0, hhp );
    for ( size_t k = 1; k <= p; k++ ) {
      dcopy( h * h, &qtrq[1][1], 1, &hhh[1][1], 1 );
															 
      for ( size_t i = 1; i <= h; i++ ) hhh[i][i] += 0.5 * llambda * ( fabs( b[i][k] ) < TINY ? INVTINY : 1.0 / fabs( b[i][k] ) );
      for ( size_t i = 1; i <= h; i++ ) hhh[i][i] += hh[i];
      inverse( h, hhh );
      for ( size_t i = 1; i <= h; i++ ) {
        double work = 0.0;
        for ( size_t j = 1; j <= h; j++ ) work += hhh[i][j] * hhp[j][k];
        b[i][k] = work;
      }
    }

    // update x      
    dgemm( false, false, n, p, h, 1.0, q, b, 0.0, x );

    // update y
    for ( size_t k = 1; k <= p; k++ ) {
      double work = 0.0;
      for ( size_t i = 1; i <= n; i++ ) work += x[i][k];
      for ( size_t j = 1; j <= m; j++ ) hmp[j][k] = work; 
    }
    for ( size_t i = 1; i <= m; i++ ) {
      for ( size_t j = 1; j <= p; j++ ) if ( fy[i][j] == 0 ) y[i][j] = ( ytilde[i][j] + hmp[i][j] ) / wc;
    }

    // update distances and calculate normalized stress
    euclidean2( n, p, x, m, y, d );
    fridge = flasso = fgroup = 0.0;
    for ( size_t i = 1; i <= h; i++ ) for ( size_t j = 1; j <= p; j++ ) fridge += b[i][j] * b[i][j];
    for ( size_t i = 1; i <= h; i++ ) for ( size_t j = 1; j <= p; j++ ) flasso += fabs( b[i][j] );
    for ( size_t i = 1; i <= h; i++ ) {
      double work = 0.0;
      for ( size_t j = 1; j <= p; j++ ) work += b[i][j] * b[i][j];
      fgroup += sqrt( work );
    }
    fnew = rlambda * fridge + llambda * flasso + glambda * fgroup;
    for ( size_t i = 1; i <= n; i++ ) {
      for ( size_t j = 1; j <= m; j++ ) {
        double work = delta[i][j] - d[i][j];
        fnew += work * work;
      }
    }

    // echo intermediate results
    if ( echo == true ) echoprogress( iter, fold, fold, fnew ); 

    // check convergence
    ( *lastdif ) = fold - fnew;
    if ( ( *lastdif ) <= -1.0 * CRIT ) break;
    double fdif = 2.0 * ( *lastdif ) / ( fold + fnew );
    if ( fdif <= FCRIT ) break;
    fold = fnew;
  }
  ( *lastiter ) = iter;

  // de-allocate memory
  freematrix( x );
  freematrix( imb );
  freematrix( xtilde );
  freematrix( ytilde );
  freematrix( qtrq );
  freematrix( hhh );
  freematrix( hhm );
  freematrix( hhp );
  freematrix( hmp );
  freevector( hh );

  return( fnew );
} // penrowresmdu

void Cpenrowresmdu( int* rn, int* rm, double* rdelta, int* rp, int* rh, double* rq, double* rb, double* ry, int* rfy, double* rd, double* rrlambda, double* rllambda, double* rglambda, int* rmaxiter, double* rfdif, double* rfvalue, int* recho )
// Function Cpenrowresmdu() performs penalized row restricted multidimensional unfolding.
{
  // transfer to C
  const size_t n = *rn;
  const size_t m = *rm;
  const size_t h = *rh;
  const size_t p = *rp;
  const size_t MAXITER = *rmaxiter;
  double** delta = getmatrix( n, m, 0.0 );
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) delta[i][j] = rdelta[k];
  double** q = getmatrix( n, h, 0.0 );
  for ( size_t j = 1, k = 0; j <= h; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) q[i][j] = rq[k];
  double** b = getmatrix( h, p, 0.0 );
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= h; i++, k++ ) b[i][j] = rb[k];
  double** y = getmatrix( m, p, 0.0 );
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= m; i++, k++ ) y[i][j] = ry[k];
  int** fy = getimatrix( m, p, 0 );
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= m; i++, k++ ) fy[i][j] = rfy[k];
  double** d = getmatrix( n, m, 0.0 );
  const double rlambda = *rrlambda;
  const double llambda = *rllambda;
  const double glambda = *rglambda;
  const double FCRIT = *rfdif;
  const bool echo = ( *recho ) != 0;

  // run function
  size_t lastiter = 0;
  double lastdif = 0.0;
  double fvalue = penrowresmdu( n, m, delta, p, h, q, b, y, fy, d, rlambda, llambda, glambda, MAXITER, FCRIT, &lastiter, &lastdif, echo );

  // transfer to R
  for ( size_t j = 1, k = 0; j <= h; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) rq[k] = q[i][j];
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= h; i++, k++ ) rb[k] = b[i][j];
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= m; i++, k++ ) ry[k] = y[i][j];
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) rd[k] = d[i][j];
  ( *rmaxiter ) = ( int ) ( lastiter );
  ( *rfdif ) = lastdif;
  ( *rfvalue ) = fvalue;

  // de-allocate memory
  freematrix( delta );
  freematrix( q );
  freematrix( b );
  freematrix( y );
  freeimatrix( fy );
  freematrix( d );

} // Cpenrowresmdu
