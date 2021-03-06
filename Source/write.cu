#pragma once  // ヘッダファイルにはこれを書く

void Record
	(	int k, REAL tau, REAL E0, REAL E,
		CMPLX var_mvx, CMPLX var_mvy, CMPLX var_mv, CMPLX mvx, CMPLX mvy, CMPLX var_Px, CMPLX var_Py, CMPLX var_P,
		REAL x_av, REAL y_av, REAL var_x, REAL var_y, REAL var_r, CMPLX Px, REAL err)
{	FILE *f;
	if ( k == 0 )
	{
		#ifndef _UFM_
		if ( ( f = fopen("000.dat", "w") ) == NULL ) { printf("file open error!!__000.dat\n"); exit(EXIT_FAILURE); }
		#else
		if ( ( f = fopen("000.ufm", "w") ) == NULL ) { printf("file open error!!__000.ufm\n"); exit(EXIT_FAILURE); }
		#endif
			fprintf(f,"%9s%18s%18s%18s%18s%18s%18s%18s%18s%18s%18s%25s%25s%18s%18s%18s%18s%18s\n",
			"#    time", "<x>", "<y>", "var_x", "var_y", "var_r", "<mv_x>", "<mv_y>", "var_mvx", "var_mvy", "var_mv", "E", "E-E0", "<Px>", "var_Px", "var_Py", "var_P", "err");
	}
	else
	{
		#ifndef _UFM_
		if ( ( f = fopen("000.dat", "a") ) == NULL ) { printf("file open error!!__000.dat\n"); exit(EXIT_FAILURE); }
		#else
		if ( ( f = fopen("000.ufm", "a") ) == NULL ) { printf("file open error!!__000.ufm\n"); exit(EXIT_FAILURE); }
		#endif
	}
/*	fprintf(f,"%9f%18.9e%18.9e%18.9e%18.9e%18.9e%18.9e%18.9e%18.9e%25.16e%25.16e%18.9e%18.9e\n",
	#ifdef CUDACOMPLEX_H
				k*tau, x_av, y_av, var_x, var_y, var_r, mvx.real(), mvy.real(), var_mv.real(), E, E - E0, Px.real(), var_P.real());
	#else
				k*tau, x_av, y_av, var_x, var_y, var_r, mvx.x     , mvy.x     , var_mv.x     , E, E - E0, Px.x     , var_P.x     );
	#endif
*/ fprintf(f,"%9f%18.9e%18.9e%18.9e%18.9e%18.9e%18.9e%18.9e%18.9e%18.9e%18.9e%25.16e%25.16e%18.9e%18.9e%18.9e%18.9e%18.9e\n",
	#ifdef CUDACOMPLEX_H
				k*tau, x_av, y_av, var_x, var_y, var_r, mvx.real(), mvy.real(), var_mvx.real(), var_mvy.real(), var_mv.real(), E, E - E0, Px.real(), var_Px.real(), var_Py.real(), var_P.real(), err);
	#else
				k*tau, x_av, y_av, var_x, var_y, var_r, mvx.x     , mvy.x     , var_mvx.x     , var_mvy.x     , var_mv.x     , E, E - E0, Px.x     , var_Px.x     , var_Py.x     , var_P.x     , err);
	#endif
	fclose(f);
}

void rho ( int k, int Nx2, int Ny2, REAL hh, REAL *x, REAL *y, CMPLX *Psi)
{	char str[20];
	#ifndef _UFM_
	sprintf(str,"OUT/%08d.rho",k);
	#else
	sprintf(str,"out/%08d.rho",k);
	#endif
	FILE *f;
	int rec = 10 * 1;	// hx = hy = 2e-2 のときに rec = 10 が基本
	if ( ( f = fopen(str, "w") ) == NULL ) { printf("file open error!!__%s\n", str); exit(EXIT_FAILURE); }
	for(		int iy = 0; iy < Ny2; iy++)
	{	for(	int ix = 0; ix < Nx2; ix++)
		{	if( ix % rec == 0 && iy % rec == 0 )
			{	int II = ix +  Nx2 * iy;
				#ifdef CUDACOMPLEX_H
				//fprintf(f, "%11.4e %11.4e %11.4e\n", x[ix], y[iy], float( Square(         Psi[II].abs()   ) * hh) );
				fprintf(f, "%25.16e %25.16e %25.16e\n", x[ix], y[iy], float( ( Square( Psi[II].real() ) + Square( Psi[II].imag() ) )  * hh) );
				#else
				//fprintf(f, "%25.16e %25.16e %25.16e %25.16e %25.16e\n", x[ix], y[iy], float( Square( cuCabs( Psi[II] )       ) * hh), float(Psi[II].x) * hh, float(Psi[II].y) * hh );
				fprintf(f, "%11.3e %11.3e %16.7e %16.7e %16.7e\n", x[ix], y[iy], float( Square( cuCabs( Psi[II] )       ) * hh), float(Psi[II].x) * hh, float(Psi[II].y) * hh );
				#endif
			}
		}	if( iy % rec == 0 ) fprintf(f, "\n");
	}	fclose(f);
}

void pot ( int Nx2, int Ny2, REAL *x, REAL *y, REAL *u)
{	FILE *f;
	int rec = 10 * 1;	// hx = hy = 2e-2 のときに rec = 10 が基本
	if ( ( f = fopen("potential.dat", "w") ) == NULL ) { printf("file open error!!__potential.dat\n"); exit(EXIT_FAILURE); }
	for(		int iy = 0; iy < Ny2; iy++)
	{	for(	int ix = 0; ix < Nx2; ix++)
		{	if( ix % rec == 0 && iy % rec == 0 )
			{	int II = ix +  Nx2 * iy;
				fprintf(f, "%11.4e %11.4e %11.4e\n", x[ix], y[iy], float( u[II] ) );
			}
		}	if( iy % rec == 0 ) fprintf(f, "\n");
	}	fclose(f);
}
