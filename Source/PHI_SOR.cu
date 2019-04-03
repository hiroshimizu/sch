#pragma once  // ヘッダファイルにはこれを書く

	__constant__ int   Nb_d, Nx2_d, Ny2_d, Nb3__d, Nx_d, Ny_d;
	__constant__ REAL  sgm2_d, pi_d, hx_d, hy_d, hh_d, omg_d, x0_d, y0_d, Bz_d, q_d;
	__constant__ CMPLX beta_d, gamma_d, zeta_d, xi_d;

void SetConst
	(	int Nb, int Nx2, int Ny2, int Nb3_, int Nx, int Ny, REAL sgm2, REAL pi, REAL hx, REAL hy,
		REAL hh, REAL omg, REAL x0, REAL y0, REAL Bz, REAL q, CMPLX beta, CMPLX gamma, CMPLX zeta, CMPLX xi)
{	(cudaMemcpyToSymbol(    Nb_d,    &Nb, sizeof(Nb)   ));
	(cudaMemcpyToSymbol(   Nx2_d,   &Nx2, sizeof(Nx2)  ));
	(cudaMemcpyToSymbol(   Ny2_d,   &Ny2, sizeof(Ny2)  ));
	(cudaMemcpyToSymbol(  Nb3__d,  &Nb3_, sizeof(Nb3_) ));
	(cudaMemcpyToSymbol(    Nx_d,    &Nx, sizeof(Nx)   ));
	(cudaMemcpyToSymbol(    Ny_d,    &Ny, sizeof(Ny)   ));
	(cudaMemcpyToSymbol(  sgm2_d,  &sgm2, sizeof(sgm2) ));
	(cudaMemcpyToSymbol(    pi_d,    &pi, sizeof(pi)   ));
	(cudaMemcpyToSymbol(    hx_d,    &hx, sizeof(hx)   ));
	(cudaMemcpyToSymbol(    hy_d,    &hy, sizeof(hy)   ));
	(cudaMemcpyToSymbol(    hh_d,    &hh, sizeof(hh)   ));
	(cudaMemcpyToSymbol(   omg_d,   &omg, sizeof(omg)  ));
	(cudaMemcpyToSymbol(    x0_d,    &x0, sizeof(x0)   ));
	(cudaMemcpyToSymbol(    y0_d,    &y0, sizeof(y0)   ));
	(cudaMemcpyToSymbol(    Bz_d,    &Bz, sizeof(Bz)   ));
	(cudaMemcpyToSymbol(     q_d,     &q, sizeof(q)    ));
	(cudaMemcpyToSymbol(  beta_d,  &beta, sizeof(beta) ));
	(cudaMemcpyToSymbol( gamma_d, &gamma, sizeof(gamma)));
	(cudaMemcpyToSymbol(  zeta_d,  &zeta, sizeof(zeta) ));
	(cudaMemcpyToSymbol(    xi_d,    &xi, sizeof(xi)   ));
}

#ifdef FieldParticle
__global__ void PHI ( CMPLX *psi, CMPLX *phi, REAL *u, int N_field )
#else
__global__ void PHI ( CMPLX *psi, CMPLX *phi, int N_field )
#endif
{	CMPLX dpsi_dx, dpsi_dy, d2psi_dx2, d2psi_dy2;
	CMPLX BETA_ij, GAMMA_ij, ZETA_ij, XI_ij;
	REAL  qAx, qAy;
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	int ix = idx % Nx2_d;
	int iy = idx / Nx2_d;
	int II = ix  + Nx2_d * iy;
	//phi[II] = 0;
	if ( 0 < ix && ix < Nx2_d - 1 && 0 < iy && iy <= Nb_d )
	{	REAL x  = ( ix - Nx_d ) * hx_d;
		REAL y  = ( iy - Ny_d + N_field * Nb_d ) * hy_d;
		qA ( x, y, Bz_d, q_d, &qAx, &qAy );
		FDM ( II, psi, hx_d, hy_d, Nx2_d, &dpsi_dx, &dpsi_dy, &d2psi_dx2, &d2psi_dy2 );

		#ifdef CUDACOMPLEX_H
			BETA_ij   =  beta_d * ( d2psi_dx2 * d2psi_dy2 );
			GAMMA_ij  = gamma_d * ( qAx * dpsi_dx + qAy * dpsi_dy );
			ZETA_ij   =  zeta_d * ( Square( qAx ) + Square( qAy ) ) * psi[II];
			#ifdef FieldParticle
				XI_ij = xi_d * u[ II + N_field * Nb_d ] * psi[II];
			#else
				XI_ij = xi_d *        potential( x, y, q_d ) * psi[II];
			#endif
			phi[II] = psi[II] + BETA_ij + GAMMA_ij - ZETA_ij - XI_ij;
		#else
			BETA_ij   = cuCmul(  beta_d, cuCadd( d2psi_dx2, d2psi_dy2 ) );
			GAMMA_ij  = cuCmul( gamma_d, cuCadd( cuCmul( mkCMPLX( qAx, 0 ), dpsi_dx ), cuCmul( mkCMPLX( qAy, 0 ), dpsi_dy ) ) );
			ZETA_ij   = cuCmul(  zeta_d, cuCmul( mkCMPLX( Square( qAx ) + Square( qAy ), 0 ), psi[II] ) );
			#ifdef FieldParticle
				XI_ij = cuCmul( xi_d, cuCmul( mkCMPLX( u[ II + N_field * Nb_d ], 0 ), psi[II] ) );
			#else
				XI_ij = cuCmul( xi_d, cuCmul( mkCMPLX( potential( x, y, q_d ), 0 ), psi[II] ) );
			#endif
			phi[II] = cuCsub( cuCsub( cuCadd( cuCadd( psi[II], BETA_ij ), GAMMA_ij ), ZETA_ij ), XI_ij );
		#endif
	}
	__syncthreads();
}


#ifdef FieldParticle
__global__ void SOR ( CMPLX *psi, CMPLX *phi, REAL *res, REAL *u, int N_field )
#else
__global__ void SOR ( CMPLX *psi, CMPLX *phi, REAL *res, int N_field )
#endif
{	CMPLX dpsi_dx, dpsi_dy, d2psi_dx2, d2psi_dy2;
	CMPLX ALPHA_ij, BETA_ij, GAMMA_ij, ZETA_ij, XI_ij;
	CMPLX ResC, tmpC = {0,0};
	REAL  qAx, qAy;
	REAL  tmpR = 0;
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	int ix = idx % Nx2_d;
	int iy = idx / Nx2_d;
	int II = ix + Nx2_d * iy;
	if( 0 < ix && ix < Nx2_d - 1 && 0 < iy && iy <= Nb_d)
	{	REAL x = ( ix - Nx_d ) * hx_d;
		REAL y = ( iy - Ny_d + N_field * Nb_d ) * hy_d;
		qA ( x, y, Bz_d, q_d, &qAx, &qAy );
		FDM ( II, psi, hx_d, hy_d, Nx2_d, &dpsi_dx, &dpsi_dy, &d2psi_dx2, &d2psi_dy2 );

		#ifdef CUDACOMPLEX_H
			BETA_ij  =  beta_d * ( d2psi_dx2 + d2psi_dy2 );
			GAMMA_ij = gamma_d * ( qAx * dpsi_dx + qAy * dpsi_dy );
			ZETA_ij  =  zeta_d * ( Square( qAx ) + Square( qAy ) ) * psi[II];
			#ifdef FieldParticle
				XI_ij = xi_d * u[ II + N_field * Nb_d ] * psi[II];
				ALPHA_ij = 1	+ beta_d * 2 * ( 1 / Square( hx_d ) + 1 / Square( hy_d ) )
					         	+ zeta_d *     (     Square( qAx  ) +     Square( qAy  ) )
					         	+   xi_d * u[ II + N_field * Nb_d ];
			#else
				XI_ij = cuCmul( xi_d, cuCmul( mkCMPLX( potential( x, y, q_d ), 0 ), psi[II] ) );
				ALPHA_ij = cuCadd( cuCadd( cuCadd( mkCMPLX( 1, 0 ),
					cuCmul( beta_d, mkCMPLX( 2 * ( 1 / Square( hx_d ) + 1 / Square( hy_d ) ), 0 ) ) ),
					cuCmul( zeta_d, mkCMPLX(           Square( qAx )  +     Square( qAy  )  , 0 ) ) ),
					cuCmul(   xi_d, mkCMPLX( potential( x, y, q_d )                              , 0 ) ) );
			#endif
			ResC = phi[II] - psi[II] + BETA_ij + GAMMA_ij - ZETA_ij - XI_ij;
			tmpC = omg_d / ALPHA_ij * ResC;
			tmpR = Square( ResC.real() ) + Square( ResC.imag() );
		#else
			BETA_ij  = cuCmul(  beta_d, cuCadd( d2psi_dx2, d2psi_dy2 ) );
			GAMMA_ij = cuCmul( gamma_d, cuCadd( cuCmul( mkCMPLX( qAx, 0 ), dpsi_dx ), cuCmul( mkCMPLX( qAy, 0 ), dpsi_dy ) ) );
			ZETA_ij  = cuCmul(  zeta_d, cuCmul( mkCMPLX( Square( qAx ) + Square( qAy ), 0 ), psi[II] ) );
			#ifdef FieldParticle
				XI_ij = cuCmul( xi_d, cuCmul( mkCMPLX( u[ II + N_field * Nb_d ], 0 ), psi[II] ) );
				ALPHA_ij = cuCadd( cuCadd( cuCadd( mkCMPLX( 1, 0 ),
					cuCmul( beta_d, mkCMPLX( 2 * ( 1 / Square( hx_d ) + 1 / Square( hy_d ) ), 0 ) ) ),
					cuCmul( zeta_d, mkCMPLX( Square( qAx ) + Square( qAy ), 0 ) ) ),
					cuCmul( xi_d, mkCMPLX( u[ II + N_field * Nb_d ], 0 ) ) );
			#else
				XI_ij = cuCmul( xi_d, cuCmul( mkCMPLX( potential( x, y, q_d ), 0 ), psi[II] ) );
				ALPHA_ij = cuCadd( cuCadd( cuCadd( mkCMPLX( 1, 0 ),
					cuCmul( beta_d, mkCMPLX( 2 * ( 1 / Square( hx_d ) + 1 / Square( hy_d ) ), 0 ) ) ),
					cuCmul( zeta_d, mkCMPLX( Square( qAx ) + Square( qAy ), 0 ) ) ),
					cuCmul( xi_d, mkCMPLX( potential( x, y, q_d ), 0 ) ) );
			#endif
			ResC = cuCsub( phi[II], cuCsub( psi[II], cuCadd( BETA_ij, cuCsub( cuCsub( GAMMA_ij, ZETA_ij ), XI_ij ) ) ) ) ;
			tmpC = cuCmul( cuCdiv( mkCMPLX( omg_d, 0 ), ALPHA_ij ), ResC );
			tmpR = Square( cuCabs( ResC ) );
		#endif
	__syncthreads();
	#ifdef CUDACOMPLEX_H
		psi[II] =         psi[II] + tmpC  ;
	#else
		psi[II] = cuCadd( psi[II],  tmpC );
	#endif
	res[II] = tmpR ;
	}
}


#ifdef FieldParticle
__global__ void sor0 ( CMPLX *psi, CMPLX *phi, REAL *res, REAL *u, int N_field )
#else
__global__ void sor0 ( CMPLX *psi, CMPLX *phi, REAL *res, int N_field )
#endif
{	int iy = blockDim.x * blockIdx.x + threadIdx.x;

	__shared__ CMPLX beta_s, gamma_s, zeta_s, xi_s;
	if ( threadIdx.x == 0 )
	{	beta_s  =  beta_d;
		gamma_s = gamma_d;
		zeta_s  =  zeta_d;
		xi_s    =    xi_d;
	}	__syncthreads();

	CMPLX dpsi_dx, dpsi_dy, d2psi_dx2, d2psi_dy2;
	CMPLX ALPHA_ij, BETA_ij, GAMMA_ij, ZETA_ij, XI_ij;
	CMPLX ResC;
	REAL  qAx, qAy;
	if ( iy % 2 == 0 )
	{	if ( 0 < iy && iy <= Nb_d )
		{	REAL y = ( iy - Ny_d + N_field * Nb_d ) * hy_d;
			for( int ix = 1; ix < Nx2_d - 1; ix++ )
			{	int II = ix + Nx2_d * iy;
				REAL x = ( ix - Nx_d ) * hx_d;
				qA ( x, y, Bz_d, q_d, &qAx, &qAy );
				FDM ( II, psi, hx_d, hy_d, Nx2_d, &dpsi_dx, &dpsi_dy, &d2psi_dx2, &d2psi_dy2 );

			#ifdef CUDACOMPLEX_H
				 BETA_ij =  beta_s * ( d2psi_dx2 + d2psi_dy2 );
				GAMMA_ij = gamma_s * ( qAx * dpsi_dx + qAy * dpsi_dy );
				 ZETA_ij =  zeta_s * ( Square( qAx ) + Square( qAy ) ) * psi[II];
				#ifdef FieldParticle
					XI_ij = xi_s * u[ II + N_field * Nb_d ] * psi[II];
					ALPHA_ij = 1	+ beta_s * 2 * ( 1 / Square( hx_d ) + 1 / Square( hy_d ) )
						         	+ zeta_s *     (     Square( qAx  ) +     Square( qAy  ) )
						          	+   xi_s * u[ II + N_field * Nb_d ];
				#else
					XI_ij = xi_s * potential( x, y, q_d ) * psi[II];
					ALPHA_ij = 1 	+ beta_s * 2 * ( 1 / Square( hx_d ) + 1 / Square( hy_d ) )
						          	+ zeta_s *     (     Square( qAx  ) +     Square( qAy  ) )
						          	+   xi_s * potential( x, y, q_d );
				#endif
				ResC = phi[II] - psi[II] + BETA_ij + GAMMA_ij - ZETA_ij - XI_ij ;
				psi[II] = psi[II] + omg_d / ALPHA_ij * ResC;
				res[II] = Square( ResC.real() );
			#else
				BETA_ij  = cuCmul(  beta_s, cuCadd( d2psi_dx2, d2psi_dy2 ) );
				GAMMA_ij = cuCmul( gamma_s, cuCadd( cuCmul( mkCMPLX( qAx, 0 ), dpsi_dx ), cuCmul( mkCMPLX( qAy, 0 ), dpsi_dy ) ) );
				ZETA_ij  = cuCmul(  zeta_s, cuCmul( mkCMPLX( Square( qAx ) + Square( qAy ), 0 ), psi[II] ) );
				#ifdef FieldParticle
					XI_ij = cuCmul( xi_s, cuCmul( mkCMPLX( u[ II + N_field * Nb_d ], 0 ), psi[II] ) );
					ALPHA_ij = cuCadd( cuCadd( cuCadd( mkCMPLX( 1, 0 ),
						cuCmul( beta_s, mkCMPLX( 2 * ( 1 / Square( hx_d ) + 1 / Square( hy_d ) ), 0 ) ) ),
						cuCmul( zeta_s, mkCMPLX( Square( qAx ) + Square( qAy ), 0 ) ) ),
						cuCmul( xi_s, mkCMPLX( u[ II + N_field * Nb_d ], 0 ) ) );
				#else
					XI_ij = cuCmul( xi_s, cuCmul( mkCMPLX( potential( x, y, q_d ), 0 ), psi[II] ) );
					ALPHA_ij = cuCadd( cuCadd( cuCadd( mkCMPLX( 1, 0 ),
						cuCmul( beta_s, mkCMPLX( 2 * ( 1 / Square( hx_d ) + 1 / Square( hy_d ) ), 0 ) ) ),
						cuCmul( zeta_s, mkCMPLX( Square( qAx ) + Square( qAy ), 0 ) ) ),
						cuCmul( xi_s, mkCMPLX( potential( x, y, q_d ), 0 ) ) );
				#endif
				ResC = cuCsub( phi[II], cuCsub( psi[II], cuCadd( BETA_ij, cuCsub( cuCsub( GAMMA_ij, ZETA_ij ), XI_ij ) ) ) ) ;
				psi[II] = cuCadd( psi[II], cuCmul( cuCdiv( mkCMPLX( omg_d, 0 ), ALPHA_ij ), ResC ) );
				res[II] = Square( cuCabs( ResC ) );
			#endif
			}
			res[            Nx2_d * iy] = 0;	//	ix = 0
			res[Nx2_d - 1 + Nx2_d * iy] = 0;	//	ix = Nx2_d - 1
		}
		else if ( iy == 0 || iy == Nb_d + 1 )
		{	for( int ix = 0; ix < Nx2_d; ix++ )
			{	int II  = ix + Nx2_d * iy;
				res[II] = 0;
			}
		}
	}
	__syncthreads();
}

#ifdef FieldParticle
__global__ void sor1 ( CMPLX *psi, CMPLX *phi, REAL *res, REAL *u, int N_field )
#else
__global__ void sor1 ( CMPLX *psi, CMPLX *phi, REAL *res, int N_field )
#endif
{	int iy = blockDim.x * blockIdx.x + threadIdx.x;

	__shared__ CMPLX beta_s, gamma_s, zeta_s, xi_s;
	if ( threadIdx.x == 0 )
	{	beta_s  =  beta_d;
		gamma_s = gamma_d;
		zeta_s  =  zeta_d;
		xi_s    =    xi_d;
	}	__syncthreads();

	CMPLX dpsi_dx, dpsi_dy, d2psi_dx2, d2psi_dy2;
	CMPLX ALPHA_ij, BETA_ij, GAMMA_ij, ZETA_ij, XI_ij;
	CMPLX ResC;
	REAL  qAx, qAy;

	if ( iy % 2 == 1 )
	{	if ( 0 < iy && iy <= Nb_d )
		{	REAL y = ( iy - Ny_d + N_field * Nb_d ) * hy_d;
			for( int ix = 1; ix < Nx2_d - 1; ix++ )
			{	int II = ix + Nx2_d * iy;
				REAL x = ( ix - Nx_d ) * hx_d;
				qA ( x, y, Bz_d, q_d, &qAx, &qAy );
				FDM ( II, psi, hx_d, hy_d, Nx2_d, &dpsi_dx, &dpsi_dy, &d2psi_dx2, &d2psi_dy2 );

			#ifdef CUDACOMPLEX_H
				BETA_ij  =  beta_s * ( d2psi_dx2 + d2psi_dy2 );
				GAMMA_ij = gamma_s * ( qAx * dpsi_dx + qAy * dpsi_dy );
				ZETA_ij  =  zeta_s * ( Square( qAx ) + Square( qAy ) ) * psi[II];
				#ifdef FieldParticle
					XI_ij = xi_s * u[ II + N_field * Nb_d ] * psi[II];
					ALPHA_ij = 1	+ beta_s * 2 * ( 1 / Square( hx_d ) + 1 / Square( hy_d ) )
						          	+ zeta_s *     (     Square( qAx  ) +     Square( qAy  ) )
						          	+   xi_s * u[ II + N_field * Nb_d ];
				#else
					XI_ij = xi_s * potential( x, y, q_d ) * psi[II];
					ALPHA_ij = 1 	+ beta_s * 2 * ( 1 / Square( hx_d ) + 1 / Square( hy_d ) )
						          	+ zeta_s *     (     Square( qAx  ) +     Square( qAy  ) )
						          	+   xi_s * potential( x, y, q_d );
				#endif
				ResC = phi[II] - psi[II] + BETA_ij + GAMMA_ij - ZETA_ij - XI_ij;
				psi[II] = psi[II] + omg_d / ALPHA_ij * ResC;
				res[II] = Square( ResC.real() ) + Square( ResC.imag() );
			#else
				BETA_ij  = cuCmul(  beta_s, cuCadd( d2psi_dx2, d2psi_dy2 ) );
				GAMMA_ij = cuCmul( gamma_s, cuCadd( cuCmul( mkCMPLX( qAx, 0 ), dpsi_dx ), cuCmul( mkCMPLX( qAy, 0 ), dpsi_dy ) ) );
				ZETA_ij  = cuCmul(  zeta_s, cuCmul( mkCMPLX( Square( qAx ) + Square( qAy ), 0 ), psi[II] ) );
				#ifdef FieldParticle
					XI_ij = cuCmul( xi_s, cuCmul( mkCMPLX( u[ II + N_field * Nb_d ], 0 ), psi[II] ) );
					ALPHA_ij = cuCadd( cuCadd( cuCadd( mkCMPLX( 1, 0 ),
						cuCmul( beta_s, mkCMPLX( 2 * ( 1 / Square( hx_d ) + 1 / Square( hy_d ) ), 0 ) ) ),
						cuCmul( zeta_s, mkCMPLX( Square( qAx ) + Square( qAy ), 0 ) ) ),
						cuCmul( xi_s, mkCMPLX( u[ II + N_field * Nb_d ], 0 ) ) );
				#else
					XI_ij = cuCmul( xi_s, cuCmul( mkCMPLX( potential( x, y, q_d ), 0 ), psi[II] ) );
					ALPHA_ij = cuCadd( cuCadd( cuCadd( mkCMPLX( 1, 0 ),
						cuCmul( beta_s, mkCMPLX( 2 * ( 1 / Square( hx_d ) + 1 / Square( hy_d ) ), 0 ) ) ),
						cuCmul( zeta_s, mkCMPLX( Square( qAx ) + Square( qAy ), 0 ) ) ),
						cuCmul( xi_s, mkCMPLX( potential( x, y, q_d ), 0 ) ) );
				#endif
				ResC = cuCsub( phi[II], cuCsub( psi[II], cuCadd( BETA_ij, cuCsub( cuCsub( GAMMA_ij, ZETA_ij ), XI_ij ) ) ) ) ;
				psi[II] = cuCadd( psi[II], cuCmul( cuCdiv( mkCMPLX( omg_d, 0 ), ALPHA_ij ), ResC ) );
				res[II] = Square( cuCabs( ResC ) );
			#endif
			}
			res[            Nx2_d * iy] = 0;	//	ix = 0
			res[Nx2_d - 1 + Nx2_d * iy] = 0;	//	ix = Nx2_d - 1
		}
		else if ( iy == 0 || iy == Nb_d + 1 )
		{	for( int ix = 0; ix < Nx2_d; ix++ )
			{	int II  = ix + Nx2_d * iy;
				res[II] = 0;
			}
		}
	}
	__syncthreads();
}
