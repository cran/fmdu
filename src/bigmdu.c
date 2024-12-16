//
// Copyright (c) 2020 Frank M.T.A. Busing (e-mail: busing at fsw dot leidenuniv dot nl)
// FreeBSD or 2-Clause BSD or BSD-2 License applies, see Http://www.freebsd.org/copyright/freebsd-license.html
// This is a permissive non-copyleft free software license that is compatible with the GNU GPL. 
//

#include "fmdu.h"

#define IJ2K( n, i, j ) ( j * n + i )

void CRultrafastmdu( int* rn, int* rm, double* rdata, int* rp, double* rx, double* ry, int* rnsteps, double* rminrate, int* rseed )
// function CRultrafastmdu() performs multidimensional unfolding
{
  // transfer to C
  const size_t n = *rn;
  const size_t m = *rm;
  const size_t p = *rp;
  const size_t NSTEPS = *rnsteps;
  const double RCRIT = *rminrate;
  long xseed = ( long )( *rseed );
  randomize( &xseed );

  double* __restrict pdata = &rdata[0];
  double* __restrict px = &rx[0];
  double* __restrict py = &ry[0];

  // set constants
  const double EPS = DBL_EPSILON;                                          // 2.2204460492503131e-16
  const double TOL = sqrt( EPS );                                          // 1.4901161193847656e-08
  const double TINY = pow( 10.0, ( log10( EPS ) + log10( TOL ) ) / 2.0 );  // 1.8189894035458617e-12
  const double MAXRATE = 0.5;
  const size_t NSUBSETS = n + m;
  const double ALPHA = pow( RCRIT / MAXRATE, 1.0 / ( double )( NSTEPS ) );

			
	   

  // start main loop
	  
  double mu = MAXRATE;
  for ( size_t iter = 1; iter <= NSTEPS; iter++ ) {
    const double cmu = 1.0 - mu;

    // start subsets loop
    for( size_t subs = 1; subs <= NSUBSETS; subs++ ) {
		 
		   
		 

      // first and second indices
      const size_t idx = nextsize_t() % n;
      const size_t idy = nextsize_t() % m;
      const size_t idxp = idx * p;
      const size_t idyp = idy * p;

      // update coordinates
      const double d = fdist1( p, &px[idxp], &py[idyp] );
      if ( d < TINY ) continue;
      const double delta = pdata[IJ2K( m, idy, idx )];
      const double b = delta / d;
      for ( size_t k = 0; k < p; k++ ) {
        const double x = px[idxp + k];
        const double y = py[idyp + k];
        const double t = b * ( x - y );
        px[idxp + k] = cmu * x + mu * ( t + y );
        py[idyp + k] = cmu * y + mu * ( x - t );
      }
    }

    // exponentially decrease mu by alpha
    mu *= ALPHA;
  }
} // CRultrafastmdu
		 

void CRultrafastmdufxd( int* rn, int* rm, double* rdata, int* rp, double* rx, int* rfx, double* ry, int* rfy, int* rnsteps, double* rminrate, int* rseed )
// function CRultrafastmdufxd() performs multidimensional unfolding allowing anchors
{
  // transfer to C
  const size_t n = *rn;
  const size_t m = *rm;
  const size_t p = *rp;
  const size_t NSTEPS = *rnsteps;
  const double RCRIT = *rminrate;
  long xseed = ( long )( *rseed );
  randomize( &xseed );

  double* __restrict pdata = &rdata[0];
								 
  double* __restrict px = &rx[0];
  double* __restrict py = &ry[0];
  int* __restrict pfx = &rfx[0];
  int* __restrict pfy = &rfy[0];

  // set constants
  const double EPS = DBL_EPSILON;                                          // 2.2204460492503131e-16
  const double TOL = sqrt( EPS );                                          // 1.4901161193847656e-08
  const double TINY = pow( 10.0, ( log10( EPS ) + log10( TOL ) ) / 2.0 );  // 1.8189894035458617e-12
  const double MAXRATE = 0.5;
  const size_t NSUBSETS = n + m;
  const double ALPHA = pow( RCRIT / MAXRATE, 1.0 / ( double )( NSTEPS ) );

  // start main loop
  double mu = MAXRATE;
  for ( size_t iter = 1; iter <= NSTEPS; iter++ ) {
    const double cmu = 1.0 - mu;

    // start subsets loop
    for( size_t subs = 1; subs <= NSUBSETS; subs++ ) {

      // first and second indices
      const size_t idx = nextsize_t() % n;
      const size_t idy = nextsize_t() % m;
      const size_t idxp = idx * p;
      const size_t idyp = idy * p;

      // update coordinates
      const double d = fdist1( p, &px[idxp], &py[idyp] );
      if ( d < TINY ) continue;
      const double delta = pdata[IJ2K( m, idy, idx )];
      const double b = delta / d;
      for ( size_t k = 0; k < p; k++ ) {
        const double x = px[idxp + k];
        const double y = py[idyp + k];
        const double t = b * ( x - y );
        if ( pfx[idxp + k] == 0 ) px[idxp + k] = cmu * x + mu * ( t + y );
        if ( pfy[idyp + k] == 0 ) py[idyp + k] = cmu * y + mu * ( x - t );
      }
    }

    // exponentially decrease mu by alpha
    mu *= ALPHA;
  }
} // CRultrafastmdufxd

void CRultrafastwgtmdu( int* rn, int* rm, double* rdata, int* rw, int* rp, double* rx, double* ry, int* rnsteps, double* rminrate, int* rseed )
// function CRultrafastwgtmdu() performs weighted multidimensional unfolding
{
  // transfer to C
  const size_t n = *rn;
  const size_t m = *rm;
  const size_t p = *rp;
  const size_t NSTEPS = *rnsteps;
  const double RCRIT = *rminrate;
  long xseed = ( long )( *rseed );
  randomize( &xseed );

  double* __restrict pdata = &rdata[0];
  int* __restrict pw = &rw[0];
  double* __restrict px = &rx[0];
  double* __restrict py = &ry[0];

  // set constants
  const double EPS = DBL_EPSILON;                                          // 2.2204460492503131e-16
  const double TOL = sqrt( EPS );                                          // 1.4901161193847656e-08
  const double TINY = pow( 10.0, ( log10( EPS ) + log10( TOL ) ) / 2.0 );  // 1.8189894035458617e-12
  const double MAXRATE = 0.5;
  const size_t NSUBSETS = n + m;
  const double ALPHA = pow( RCRIT / MAXRATE, 1.0 / ( double )( NSTEPS ) );

  // start main loop
  double mu = MAXRATE;
  for ( size_t iter = 1; iter <= NSTEPS; iter++ ) {
    const double cmu = 1.0 - mu;

    // start subsets loop
    for( size_t subs = 1; subs <= NSUBSETS; subs++ ) {

      // first and second indices
      const size_t idx = nextsize_t() % n;
      const size_t idy = nextsize_t() % m;

      if ( pw[IJ2K( m, idy, idx )] == 0 ) continue;

      const size_t idxp = idx * p;
      const size_t idyp = idy * p;

      // update coordinates
      const double d = fdist1( p, &px[idxp], &py[idyp] );
      if ( d < TINY ) continue;
      const double delta = pdata[IJ2K( m, idy, idx )];
      const double b = delta / d;
      for ( size_t k = 0; k < p; k++ ) {
        const double x = px[idxp + k];
        const double y = py[idyp + k];
        const double t = b * ( x - y );
        px[idxp + k] = cmu * x + mu * ( t + y );
        py[idyp + k] = cmu * y + mu * ( x - t );
      }
    }

    // exponentially decrease mu by alpha
    mu *= ALPHA;
  }
} // CRultrafastwgtmdu

void CRultrafastwgtmdufxd( int* rn, int* rm, double* rdata, int* rw, int* rp, double* rx, int* rfx, double* ry, int* rfy, int* rnsteps, double* rminrate, int* rseed )
// function CRultrafastwgtmdufxd() performs weighted multidimensional unfolding allowing anchors
{
  // transfer to C
  const size_t n = *rn;
  const size_t m = *rm;
  const size_t p = *rp;
  const size_t NSTEPS = *rnsteps;
  const double RCRIT = *rminrate;
  long xseed = ( long )( *rseed );
  randomize( &xseed );

  double* __restrict pdata = &rdata[0];
  int* __restrict pw = &rw[0];
  double* __restrict px = &rx[0];
  double* __restrict py = &ry[0];
  int* __restrict pfx = &rfx[0];
  int* __restrict pfy = &rfy[0];

  // set constants
  const double EPS = DBL_EPSILON;                                          // 2.2204460492503131e-16
  const double TOL = sqrt( EPS );                                          // 1.4901161193847656e-08
  const double TINY = pow( 10.0, ( log10( EPS ) + log10( TOL ) ) / 2.0 );  // 1.8189894035458617e-12
  const double MAXRATE = 0.5;
  const size_t NSUBSETS = n + m;
  const double ALPHA = pow( RCRIT / MAXRATE, 1.0 / ( double )( NSTEPS ) );

  // start main loop
  double mu = MAXRATE;
  for ( size_t iter = 1; iter <= NSTEPS; iter++ ) {
    const double cmu = 1.0 - mu;

    // start subsets loop
    for( size_t subs = 1; subs <= NSUBSETS; subs++ ) {

      // first and second indices
      const size_t idx = nextsize_t() % n;
      const size_t idy = nextsize_t() % m;

      if ( pw[IJ2K( m, idy, idx )] == 0 ) continue;

      const size_t idxp = idx * p;
      const size_t idyp = idy * p;

      // update coordinates
      const double d = fdist1( p, &px[idxp], &py[idyp] );
      if ( d < TINY ) continue;
      const double delta = pdata[IJ2K( m, idy, idx )];
      const double b = delta / d;
      for ( size_t k = 0; k < p; k++ ) {
        const double x = px[idxp + k];
        const double y = py[idyp + k];
        const double t = b * ( x - y );
        if ( pfx[idxp + k] == 0 ) px[idxp + k] = cmu * x + mu * ( t + y );
        if ( pfy[idyp + k] == 0 ) py[idyp + k] = cmu * y + mu * ( x - t );
      }
    }

    // exponentially decrease mu by alpha
    mu *= ALPHA;
  }
} // CRultrafastwgtmdufxd

void CRultrafastmdu2( int* rn, int* rm, double* rdata, int* rp, double* rx, double* ry, int* rnsteps, double* rminrate, int* rseed )
// function CRultrafastmdu() performs multidimensional unfolding
{
  // transfer to C
  const size_t n = *rn;
  const size_t m = *rm;
  const size_t p = *rp;
  const size_t NSTEPS = *rnsteps;
  const double RCRIT = *rminrate;
  long xseed = ( long )( *rseed );
  randomize( &xseed );

  double* __restrict pdata = &rdata[0];
  double* __restrict px = &rx[0];
  double* __restrict py = &ry[0];

  // set constants
  const double EPS = DBL_EPSILON;                                          // 2.2204460492503131e-16
  const double TOL = sqrt( EPS );                                          // 1.4901161193847656e-08
  const double TINY = pow( 10.0, ( log10( EPS ) + log10( TOL ) ) / 2.0 );  // 1.8189894035458617e-12
  const double MAXRATE = 0.5;
  const size_t NSUBSETS = n + m;
  const double ALPHA = pow( RCRIT / MAXRATE, 1.0 / ( double )( NSTEPS ) );

  // start main loop
  double mu = MAXRATE;
  for ( size_t iter = 1; iter <= NSTEPS; iter++ ) {
    const double cmu = 1.0 - mu;

    // start subsets loop
    for( size_t subs = 1; subs <= NSUBSETS; subs++ ) {

      const size_t idx1 = nextsize_t() % n;
      const size_t idx2 = nextsize_t() % n;
      const size_t idy1 = nextsize_t() % m;
      const size_t idy2 = nextsize_t() % m;
      const size_t idxp1 = idx1 * p;
      const size_t idxp2 = idx2 * p;
      const size_t idyp1 = idy1 * p;
      const size_t idyp2 = idy2 * p;

      // update coordinates
      const double d11 = fdist1( p, &px[idxp1], &py[idyp1] );
      const double d12 = fdist1( p, &px[idxp1], &py[idyp2] );
      const double d21 = fdist1( p, &px[idxp2], &py[idyp1] );
      const double d22 = fdist1( p, &px[idxp2], &py[idyp2] );
      const double w11 = 1.0;//pw[IJ2K( m, idy1, idx1 )];
      const double w12 = 1.0;//pw[IJ2K( m, idy2, idx1 )];
      const double w21 = 1.0;//pw[IJ2K( m, idy1, idx2 )];
      const double w22 = 1.0;//pw[IJ2K( m, idy2, idx2 )];
      const double r1 = w11 + w12;
      const double r2 = w21 + w22;
      const double c1 = w11 + w21;
      const double c2 = w12 + w22;
      const double delta11 = pdata[IJ2K( m, idy1, idx1 )];
      const double delta12 = pdata[IJ2K( m, idy2, idx1 )];
      const double delta21 = pdata[IJ2K( m, idy1, idx2 )];
      const double delta22 = pdata[IJ2K( m, idy2, idx2 )];
      const double b11 = ( d11 < TINY ? 0.0 : w11 * delta11 / d11 );
      const double b12 = ( d12 < TINY ? 0.0 : w12 * delta12 / d12 );
      const double b21 = ( d21 < TINY ? 0.0 : w21 * delta21 / d21 );
      const double b22 = ( d22 < TINY ? 0.0 : w22 * delta22 / d22 );
      const double p1 = b11 + b12;
      const double p2 = b21 + b22;
      const double q1 = b11 + b21;
      const double q2 = b12 + b22;
      for ( size_t k = 0; k < p; k++ ) {
        const double x1 = px[idxp1 + k];
        const double x2 = px[idxp2 + k];
        const double y1 = py[idyp1 + k];
        const double y2 = py[idyp2 + k];
        px[idxp1 + k] = cmu * x1 + mu * ( p1 * x1 - b11 * y1 - b12 * y2 + y1 + y2 ) / r1;
        px[idxp2 + k] = cmu * x2 + mu * ( p2 * x2 - b21 * y1 - b22 * y2 + y1 + y2 ) / r2;
        py[idyp1 + k] = cmu * y1 + mu * ( q1 * y1 - b11 * x1 - b21 * x2 + x1 + x2 ) / c1;
        py[idyp2 + k] = cmu * y2 + mu * ( q2 * y2 - b12 * x1 - b22 * x2 + x1 + x2 ) / c2;
      }
    }

    // exponentially decrease mu by alpha
    mu *= ALPHA;
  }
} // CRultrafastmdu2

void CRultrafastrowresmdu( int* rn, int* rm, double* rdata, int* rp, int* rh, double* rq, double* rb, double* ry, int* rnsteps, double* rminrate, int* rseed )
// function CRultrafastrowresmdu() performs multidimensional unfolding
{
  // transfer to C
  const size_t n = *rn;
  const size_t m = *rm;
  const size_t h = *rh;
  const size_t p = *rp;
  const size_t NSTEPS = *rnsteps;
  const double RCRIT = *rminrate;
  long xseed = ( long )( *rseed );
  randomize( &xseed );

  double* __restrict pdata = &rdata[0];
  double* __restrict pq = &rq[0];
  double* __restrict pb = &rb[0];
  double* __restrict py = &ry[0];

  double* __restrict px = ( double* ) calloc( p, sizeof( double ) );
  double* __restrict ps = ( double* ) calloc( p, sizeof( double ) );
  double* __restrict pqtrq = ( double* ) calloc( h, sizeof( double ) );
  for ( size_t j = 0; j < h; j++ ) pqtrq[j] = ( double )( m ) * dssq( p, pq, 1 );

  // set constants
  const double EPS = DBL_EPSILON;                                          // 2.2204460492503131e-16
  const double TOL = sqrt( EPS );                                          // 1.4901161193847656e-08
  const double TINY = pow( 10.0, ( log10( EPS ) + log10( TOL ) ) / 2.0 );  // 1.8189894035458617e-12
  const double MAXRATE = 0.5;
  const size_t NSUBSETS = n + m;
  const double ALPHA = pow( RCRIT / MAXRATE, 1.0 / ( double )( NSTEPS ) );

  // start main loop
  double mu = MAXRATE;
  for ( size_t iter = 1; iter <= NSTEPS; iter++ ) {
    const double cmu = 1.0 - mu;

    const double ssqb = dssq( h * p, pb, 1 );
    dscal( h * p, sqrt( ( double )( h * p ) / ssqb ), pb, 1 ); 

    // start subsets loop
    for( size_t subs = 1; subs <= NSUBSETS; subs++ ) {

      // first and second indices
      const size_t idx = nextsize_t() % n;
      const size_t idy = nextsize_t() % m;
      const size_t idqp = idx * h;
      const size_t idyp = idy * p;

      // update coordinates
      memset( px, 0, p * sizeof( double ) );
      memset( ps, 0, p * sizeof( double ) );
      for ( size_t k = 0; k < p; k++ ) {
        for ( size_t j = 0; j < h; j++ ) px[k] += pq[idqp + j] * pb[k + j * p];
        for ( size_t j = 0; j < h; j++ ) ps[k] += pb[k + j * p];
      }
      const double d = fdist1( p, &px[0], &py[idyp] );
      if ( d < TINY ) continue;
      const double delta = pdata[IJ2K( m, idy, idx )];
      const double b = delta / d;
      for ( size_t k = 0; k < p; k++ ) {
        const double y = py[idyp + k];
        const double t = b * ( px[k] - y );
        py[idyp + k] = cmu * y + mu * ( px[k] - t );
        for ( size_t j = 0; j < h; j++ ) {
          const double qxtilde = pq[idqp + j] * t;
          const double qwy = pq[idqp + j] * y;
          const double smin = ps[j] - pqtrq[j] * pb[k + j * p];
          const double bnew = ( qxtilde + qwy - smin ) / pqtrq[j];
          pb[k + j * p] = ( 1.0 - mu / ( double )( 1000 + iter ) ) * pb[k + j * p] + ( mu / ( double )( 1000 + iter ) ) * bnew;
        }
      }
    }

    // exponentially decrease mu by alpha
    mu *= ALPHA;
  }

  free( px );
  free( ps );
  free( pqtrq );

} // CRultrafastrowresmdu