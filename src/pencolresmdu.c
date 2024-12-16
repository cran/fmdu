//
// Copyright (c) 2020 Frank M.T.A. Busing (e-mail: busing at fsw dot leidenuniv dot nl)
// FreeBSD or 2-Clause BSD or BSD-2 License applies, see Http://www.freebsd.org/copyright/freebsd-license.html
// This is a permissive non-copyleft free software license that is compatible with the GNU GPL. 
//

#include "fmdu.h"

double pencolresmdu( const size_t n, const size_t m, double** delta, const size_t p, double** x, int** fx, const size_t h, double** q, double** b, double** d, const double rlambda, const double llambda, const double glambda, const size_t MAXITER, const double FCRIT, size_t* lastiter, double* lastdif, const bool echo )
// Function pencolresmdu() performs penalized column restricted multidimensional unfolding.
{
  const double EPS = DBL_EPSILON;                                          // 2.2204460492503131e-16
  const double TOL = sqrt( EPS );                                          // 1.4901161193847656e-08
  const double CRIT = sqrt( TOL );                                         // 0.00012207031250000000
  const double TINY = pow( 10.0, ( log10( EPS ) + log10( TOL ) ) / 2.0 );  // 1.8189894035458617e-12
  const double INVTINY = 1.0 / TINY;

  // allocate memory
  double** y = getmatrix( m, p, 0.0 );
  double** imb = getmatrix( n, m, 0.0 );
  double** xtilde = getmatrix( n, p, 0.0 );
  double** ytilde = getmatrix( m, p, 0.0 );
  double** qtcq = getmatrix( h, h, 0.0 );
  double** hhh = getmatrix( h, h, 0.0 );
  double** hhn = getmatrix( h, n, 0.0 );
  double** hhp = getmatrix( h, p, 0.0 );
  double** hnp = getmatrix( n, p, 0.0 );
  double* hh = getvector( h, 0.0 );

  // initialization
  double wr = ( double ) ( m );
  double wc = ( double ) ( n );
  for ( size_t i = 1; i <= h; i++ ) {
    for ( size_t j = 1; j <= h; j++ ) {
      double work = 0.0;
      for ( size_t k = 1; k <= m; k++ ) work += q[k][i] * wc * q[k][j];
      hhh[i][j] = qtcq[i][j] = work;
    }
    hhh[i][i] = qtcq[i][i] += rlambda;
  }
  inverse( h, hhh );
  for ( size_t k = 1; k <= h; k++ ) {
    double work = 0.0;
    for ( size_t j = 1; j <= m; j++ ) work += q[j][k];
    for ( size_t i = 1; i <= n; i++ ) hhn[k][i] = work;
  }
  int nfx = 0;
  for ( size_t i = 1; i <= n; i++ ) for ( size_t k = 1; k <= p; k++ ) nfx += fx[i][k];

  // update distances and calculate normalized stress
  dgemm( false, false, m, p, h, 1.0, q, b, 0.0, y );
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

    // compute preliminary updates: xtilde and ytilde
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

    // update x
    for ( size_t k = 1; k <= p; k++ ) {
      double work = 0.0;
      for ( size_t j = 1; j <= m; j++ ) work += y[j][k];
      for ( size_t i = 1; i <= n; i++ ) hnp[i][k] = work;
    }
    for ( size_t i = 1; i <= n; i++ ) {
      for ( size_t j = 1; j <= p; j++ ) if ( fx[i][j] == 0 ) x[i][j] = ( xtilde[i][j] + hnp[i][j] ) / wr;
    }

    // update b
    for ( size_t i = 1; i <= h; i++ ) {
      double work = 0.0;
      for ( size_t j = 1; j <= p; j++ ) work += b[i][j] * b[i][j];
      work = sqrt( work );
      hh[i] = 0.5 * glambda * ( work < TINY ? INVTINY : 1.0 / work );
    }
    dgemm( false, false, h, p, n, 1.0, hhn, x, 0.0, hhp );
    dgemm( true, false, h, p, n, 1.0, q, ytilde, 1.0, hhp );
    for ( size_t k = 1; k <= p; k++ ) {
      dcopy( h * h, &qtcq[1][1], 1, &hhh[1][1], 1 );
      for ( size_t i = 1; i <= h; i++ ) hhh[i][i] += 0.5 * llambda * ( fabs( b[i][k] ) < TINY ? INVTINY : 1.0 / fabs( b[i][k] ) );
      for ( size_t i = 1; i <= h; i++ ) hhh[i][i] += hh[i];
      inverse( h, hhh );
      for ( size_t i = 1; i <= h; i++ ) {
        double work = 0.0;
        for ( size_t j = 1; j <= h; j++ ) work += hhh[i][j] * hhp[j][k];
        b[i][k] = work;
      }
    }

    // update y      
    dgemm( false, false, m, p, h, 1.0, q, b, 0.0, y );

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

  // rotate to principal axes of x
  if ( nfx == 0 ) rotateplus( n, p, x, h, b );

  // de-allocate memory
  freematrix( y );
  freematrix( imb );
  freematrix( xtilde );
  freematrix( ytilde );
  freematrix( qtcq );
  freematrix( hhh );
  freematrix( hhn );
  freematrix( hhp );
  freematrix( hnp );
  freevector( hh );

  return( fnew );
} // pencolresmdu

void Cpencolresmdu( int* rn, int* rm, double* rdelta, int* rp, double* rx, int* rfx, int* rh, double* rq, double* rb, double* rd, double* rrlambda, double* rllambda, double* rglambda, int* rmaxiter, double* rfdif, double* rfvalue, int* recho )
// Function Ccolresmdu() performs column restricted multidimensional unfolding.
{
  // transfer to C
  size_t n = *rn;
  size_t m = *rm;
  size_t h = *rh;
  size_t p = *rp;
  size_t MAXITER = *rmaxiter;
  double** delta = getmatrix( n, m, 0.0 );
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) delta[i][j] = rdelta[k];
  double** x = getmatrix( n, p, 0.0 );
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) x[i][j] = rx[k];
  int** fx = getimatrix( n, p, 0 );
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) fx[i][j] = rfx[k];
  double** q = getmatrix( m, h, 0.0 );
  for ( size_t j = 1, k = 0; j <= h; j++ ) for ( size_t i = 1; i <= m; i++, k++ ) q[i][j] = rq[k];
  double** b = getmatrix( h, p, 0.0 );
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= h; i++, k++ ) b[i][j] = rb[k];
  double** d = getmatrix( n, m, 0.0 );
  const double rlambda = *rrlambda;
  const double llambda = *rllambda;
  const double glambda = *rglambda;
  double FCRIT = *rfdif;
  bool echo = ( *recho ) != 0;

  // run function
  size_t lastiter = 0;
  double lastdif = 0.0;
  double fvalue = pencolresmdu( n, m, delta, p, x, fx, h, q, b, d, rlambda, llambda, glambda, MAXITER, FCRIT, &lastiter, &lastdif, echo );

  // transfer to R
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) rx[k] = x[i][j];
  for ( size_t j = 1, k = 0; j <= h; j++ ) for ( size_t i = 1; i <= m; i++, k++ ) rq[k] = q[i][j];
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= h; i++, k++ ) rb[k] = b[i][j];
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) rd[k] = d[i][j];
  ( *rmaxiter ) = ( int ) ( lastiter );
  ( *rfdif ) = lastdif;
  ( *rfvalue ) = fvalue;

  // de-allocate memory
  freematrix( delta );
  freematrix( x );
  freeimatrix( fx );
  freematrix( q );
  freematrix( b );
  freematrix( d );

} // Cpencolresmdu
