#pragma once  // ヘッダファイルにはこれを書く

__device__ __host__
static inline REAL X ( REAL x ) { return 0; }

__device__ __host__
static inline REAL Y ( REAL y ) {
	//const REAL LB_INV = 1e-6;	// 磁場の勾配長の逆数、0のとき一様磁場
	//const REAL LB_INV = 1e-5;	// 磁場の勾配長の逆数、0のとき一様磁場
	const REAL LB_INV = 0;	// 磁場の勾配長の逆数、0のとき一様磁場
	return y * ( 1 - LB_INV / 2 * y);
}

__device__ __host__ void qA ( REAL x, REAL y, REAL Bz, REAL q, REAL* qAx, REAL* qAy )
{	*qAx = - q * Bz * Y(y);
	*qAy =   q * Bz * X(x);
}


#ifdef FieldParticle
extern "C"
{	void  hiab_();
	void  dehint_(REAL *, REAL *, REAL *, REAL *);
}

void PotEne
	(	REAL *u, const REAL *x, const REAL *y,
		const REAL q, const REAL sgm, const REAL Eps0, const int Nx2, const int Ny2
	)
{	const int  nf = 1;	// number of field particles
	REAL xf[nf], yf[nf], qf[nf];
	REAL zero = 0;		// lower bound of integration (積分区間は[a,b]=[0, infinity])
	REAL EPS  = 1e-12;	// absolute error tolerance
	REAL V;

	xf[0] = 0; yf[0] = 0; qf[0] = 3e-5;
	hiab_();
	for(		int iy = 0; iy < Ny2; iy++) { REAL y_ij = y[iy];
		for(	int ix = 0; ix < Nx2; ix++) { REAL x_ij = x[ix];
			int II = ix + Nx2 * iy;
			REAL phi_ij = 0;
			for( int f = 0; f < nf; f++)
			{	REAL R = sqrt( Square( x_ij - xf[f] ) + Square( y_ij - yf[f] ) );
				REAL eta = R / sgm;
				dehint_ ( &eta, &zero, &EPS, &V );
				phi_ij += qf[f] * V;
			}
			phi_ij =     phi_ij / ( Square( pi ) * Eps0 * sgm );
			u[II]  = q * phi_ij;
		}
	}
}
#else
	__device__ __host__ REAL potential ( REAL x, REAL y, REAL q )
	{
	//	return 0;
		#ifndef _UFM_
		const REAL LE4_INV = INP_NUM_;
		#else
		const REAL LE4_INV = 0;
		#endif
		const REAL Ey = 1;
		REAL y2 = y * y;
		return - q * Ey * y2 * y2 * y / 5 * LE4_INV;
	//	REAL Ey = 0;
	//	REAL Ey = 1;
	//	REAL Ey = 1e-4;
	//	REAL Ey = 1e-5;
	//	REAL Ey = 1e-6;
		//REAL LE_INV = 0;

	//	REAL k_E = 0.1;
	//	REAL k_E = 0.5;
	//	REAL k_E = 1;
	//	REAL k_E = 2;
	//	REAL k_E   = 1e-4;
	//	REAL lmd_E = 1e+4; // lmd_E = 1 / k_E

	/* 1. linear type */
		//return - q * Ey * y *(1 - LE_INV / 2 * y);
	/* 3. sinusoidal type E ~ + cos */
		//return - q * Ey / k_E * sin(k_E * y);
		/* 3+. 1/k_E -> lmd_Eとして小さな数が分母にくることを回避 */
		//return - q * Ey * lmd_E * sin(y / lmd_E);
	/* 3+.sinusoidal type E ~ + sin */
	//	return + q * Ey / k_E * cos(k_E * y);
	}
#endif

