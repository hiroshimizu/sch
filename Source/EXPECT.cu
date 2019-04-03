#pragma once  // ヘッダファイルにはこれを書く

void EXPECT ( CMPLX *psi, REAL *x, REAL *y,
#ifdef FieldParticle
	REAL *u,
#endif
	CMPLX *mvx_av,  CMPLX *mvy_av,  CMPLX *Px,    REAL  *E,      CMPLX *var_mv,
	CMPLX *var_mvx, CMPLX *var_mvy, CMPLX *var_P, CMPLX *var_Px, CMPLX *var_Py,
	REAL *x_av,     REAL *y_av,     REAL *var_r,  REAL *var_x,   REAL *var_y, 
	REAL *err,      REAL *qAx_av,
	REAL Bz, REAL m, REAL q, int Nx2, int Ny2, REAL hx, REAL hy, REAL h_bar )
{	REAL  x_ij, y_ij, rho_ij, qAx, qAy;
	CMPLX dpsi_dx, dpsi_dy, dpsi2dx, dpsi2dy;
	CMPLX psi_cj;
	CMPLX Px_av = {0,0}, Py_av = {0,0}, Px2av = {0,0}, Py2av = {0,0};
	CMPLX mv2av;
	CMPLX mvx2av;
	CMPLX mvy2av;
	CMPLX qAxPx_av = {0,0}; CMPLX qAyPy_av = {0,0};
	REAL  qAy_av = 0, qAx2av = 0, qAy2av = 0;
	REAL  x2av = 0, y2av = 0;
	CMPLX K_av, E_av;
	REAL  U_av = 0;
	REAL  hh = hx * hy;

	/* Calculation */
	/*	Px_av, Px2av, Py_av, Py2av, 
		qAx_av,qAx2av,qAy_av,qAy2av, 
		qAxPx_av, qAyPy_av, U_av	*/
	*qAx_av = 0;
	for(		int iy = 1; iy < Ny2 - 1; iy++) { y_ij = y[iy];
		for(	int ix = 1; ix < Nx2 - 1; ix++) { x_ij = x[ix];
			int II = ix + Nx2 * iy;
			FDM ( II, psi, hx, hy, Nx2, &dpsi_dx, &dpsi_dy, &dpsi2dx, &dpsi2dy );
			#ifdef CUDACOMPLEX_H
				psi_cj = ~psi[II];
				rho_ij = ( psi_cj * psi[II] ).real();

				Px_av = Px_av + psi_cj * dpsi_dx;
				Py_av = Py_av + psi_cj * dpsi_dy;
				Px2av = Px2av + psi_cj * dpsi2dx;
				Py2av = Py2av + psi_cj * dpsi2dy;
				qA ( x_ij, y_ij, Bz, q, &qAx, &qAy );	// qAy == 0が来ているはず
				*qAx_av +=         qAx   * rho_ij;
				 qAy_av +=         qAy   * rho_ij;	// == 0のはず
				qAx2av  += Square( qAx ) * rho_ij;
				qAy2av  += Square( qAy ) * rho_ij;	// == 0のはず
				qAxPx_av = qAxPx_av + psi_cj * ( qAx * dpsi_dx );
				qAyPy_av = qAyPy_av + psi_cj * ( qAy * dpsi_dy );	// == 0のはず
			#else
				rho_ij = Square( cuCabs( psi[II] ) );
				psi_cj = cuConj(         psi[II]   );
				Px_av  = cuCadd( Px_av, cuCmul( psi_cj, dpsi_dx ) );
				Py_av  = cuCadd( Py_av, cuCmul( psi_cj, dpsi_dy ) );
				Px2av  = cuCadd( Px2av, cuCmul( psi_cj, dpsi2dx ) );
				Py2av  = cuCadd( Py2av, cuCmul( psi_cj, dpsi2dy ) );
				qA ( x_ij, y_ij, Bz, q, &qAx, &qAy );	// qAy == 0が来ているはず
				*qAx_av +=         qAx   * rho_ij;
				 qAy_av +=         qAy   * rho_ij;	// == 0のはず
				qAx2av +=  Square( qAx ) * rho_ij;
				qAy2av +=  Square( qAy ) * rho_ij;	// == 0のはず
				qAxPx_av = cuCadd( qAxPx_av, cuCmul( psi_cj, cuCmul( mkCMPLX(qAx,0), dpsi_dx )));
				qAyPy_av = cuCadd( qAyPy_av, cuCmul( psi_cj, cuCmul( mkCMPLX(qAy,0), dpsi_dy )));	// == 0のはず
			#endif
			#ifndef FieldParticle
				U_av += potential( x_ij, y_ij, q ) * rho_ij;
			#else
				U_av +=                   u[II] * rho_ij;
			#endif
		}
	}

	/* Calculation */
	/*	x_av, x2av, y_av, y2av, err */
	*x_av = 0; *y_av = 0; *err = 0;
	for(		int iy = 0; iy < Ny2; iy++) { y_ij = y[iy];
		for(	int ix = 0; ix < Nx2; ix++) { x_ij = x[ix];
			int II = ix + Nx2 * iy;
			#ifdef CUDACOMPLEX_H
				rho_ij = psi[II].abs2();	//rho_ij = ( ~psi[II] * psi[II] ).real();
			#else
				rho_ij = Square( cuCabs( psi[II] ) );
			#endif
			*x_av += 		  x_ij   * rho_ij;
			*y_av += 		  y_ij   * rho_ij;
			x2av  += Square( x_ij ) * rho_ij;
			y2av  += Square( y_ij ) * rho_ij;
			*err	+= 				     rho_ij;
		}
	}
	
	#ifdef CUDACOMPLEX_H
		CMPLX I = {0, 1};
	//	CMPLX I = make_doublecomplex(0, 1);
	//	CMPLX I = make_REALcomplex(0, 1);
		Px_av = - h_bar * hh * I * Px_av;	// 運動の恒量for H = H( y, \vec v)
		Py_av = - h_bar * hh * I * Py_av;
		Px2av = - h_bar * h_bar * hh * Px2av;
		Py2av = - h_bar * h_bar * hh * Py2av;
		qAxPx_av = hh * qAxPx_av;
		qAyPy_av = hh * qAyPy_av;
	#else
		Px_av  = cuCmul( mkCMPLX( 0, - h_bar * hh ), Px_av );	// 運動の恒量for H = H( y, \vec v)
		Py_av  = cuCmul( mkCMPLX( 0, - h_bar * hh ), Py_av );
		Px2av  = cuCmul( mkCMPLX( - h_bar * h_bar * hh, 0 ), Px2av );
		Py2av  = cuCmul( mkCMPLX( - h_bar * h_bar * hh, 0 ), Py2av );
		qAxPx_av = cuCmul( mkCMPLX( hh, 0 ), qAxPx_av );
		qAyPy_av = cuCmul( mkCMPLX( hh, 0 ), qAyPy_av );
	#endif
	*qAx_av *= hh;
	 qAy_av *= hh;
	 qAx2av *= hh;
	 qAy2av *= hh;

	U_av  *= hh;

	*x_av *= hh;
	*y_av *= hh;
	 x2av *= hh;
	 y2av *= hh;
	*err  *= hh;
	
	/* Calculation */
	/*	mvx_av, mvx2av, mvy_av, mvy2_av, mv_av,
		K_av, E_av, and variances */
	#ifdef CUDACOMPLEX_H
		*mvx_av = Px_av - *qAx_av;
		*mvy_av = Py_av -  qAy_av;
		 mvx2av = Px2av + qAx2av + 2 * h_bar * I * qAxPx_av;	//	Axがxの関数でないため、単純に2乗を展開した形
		 mvy2av = Py2av + qAy2av + 2 * h_bar * I * qAyPy_av;	//	同上
		 mv2av  = mvx2av + mvy2av;
		*var_mvx= mvx2av       	- ( *mvx_av * *mvx_av );
		*var_mvy= mvy2av       	- ( *mvy_av * *mvy_av ); 
		*var_mv = var_mvx + var_mvy;
		*var_Px = Px2av - ( Px_av	* Px_av );
		*var_Py = Py2av - ( Py_av	* Py_av );
		*var_P  = *var_Px + *var_Py;
		K_av	  = mv2av / ( 2 * m );
		E_av    = K_av + U_av;
	#else
		*mvx_av = cuCsub( Px_av, mkCMPLX( *qAx_av, 0 ) );
		*mvy_av = cuCsub( Py_av, mkCMPLX(  qAy_av, 0 ) );
		 mvx2av = cuCadd( Px2av, cuCadd( mkCMPLX( qAx2av, 0 ), cuCmul( mkCMPLX( 0, 2 * h_bar ), qAxPx_av ) ) );
		 mvy2av = cuCadd( Py2av, cuCadd( mkCMPLX( qAy2av, 0 ), cuCmul( mkCMPLX( 0, 2 * h_bar ), qAyPy_av ) ) );
		 mv2av  = cuCadd( mvx2av, mvy2av);
		*var_mvx= cuCsub( mvx2av, cuCmul( *mvx_av, *mvx_av ) );
		*var_mvy= cuCsub( mvy2av, cuCmul( *mvy_av, *mvy_av ) );
		*var_mv = cuCadd( *var_mvx, *var_mvy );
		*var_Px = cuCsub( Px2av, cuCmul( Px_av, Px_av ) );
		*var_Py = cuCsub( Py2av, cuCmul( Py_av, Py_av ) );
		*var_P  = cuCadd( *var_Px, *var_Py );
		K_av	= cuCdiv( mv2av, mkCMPLX( 2 * m, 0 ) );
		E_av	= cuCadd( K_av, mkCMPLX( U_av, 0 ) );
	#endif
	*var_x  = x2av - Square( *x_av );
	*var_y  = y2av - Square( *y_av );
	*var_r  = *var_x + *var_y;
	*err   -= 1;
	
	*Px	= Px_av;
	#ifdef CUDACOMPLEX_H
		*E		= E_av.real();
	#else
		*E		= cuCreal( E_av );
	#endif
	/* FOR Check */
//	printf("U_av, K_av, E_av: %25.16e, %25.16e, %25.16e\n", U_av,cuCreal( K_av), cuCreal(E_av));
	/* END Check */
}
