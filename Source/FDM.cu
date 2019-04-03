#pragma once  // ヘッダファイルにはこれを書く

__device__ __host__ void FDM (int II, CMPLX *psi, REAL hx, REAL hy, int Nx2,
	CMPLX *dpsi_dx, CMPLX *dpsi_dy, CMPLX *d2psi_dx2, CMPLX *d2psi_dy2 )
{	
	REAL hx2 = Square( hx );
	REAL hy2 = Square( hy );
	#ifdef CUDACOMPLEX_H
	*dpsi_dx   = ( psi[II+1  ] - psi[II-1  ] ) / ( 2 * hx );
	*dpsi_dy   = ( psi[II+Nx2] - psi[II-Nx2] ) / ( 2 * hy );
	*d2psi_dx2 = ( psi[II+1  ] + psi[II-1  ] - 2 * psi[II] ) / hx2;
	*d2psi_dy2 = ( psi[II+Nx2] + psi[II-Nx2] - 2 * psi[II] ) / hy2;
	#else
	*dpsi_dx = cuCdiv( cuCsub( psi[II+1  ], psi[II-1  ] ), mkCMPLX( 2 * hx, 0 ) );
	*dpsi_dy = cuCdiv( cuCsub( psi[II+Nx2], psi[II-Nx2] ), mkCMPLX( 2 * hy, 0 ) );
	*d2psi_dx2 = cuCdiv( cuCsub( cuCadd( psi[II+1  ], psi[II-1  ] ), cuCmul( mkCMPLX( 2, 0 ), psi[II] ) ), mkCMPLX( hx2, 0 ) );
	*d2psi_dy2 = cuCdiv( cuCsub( cuCadd( psi[II+Nx2], psi[II-Nx2] ), cuCmul( mkCMPLX( 2, 0 ), psi[II] ) ), mkCMPLX( hy2, 0 ) );
	#endif
}

