#include <stdio.h>
//#include <cutil_inline.h>
//#include"cudacomplex.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>

#define BLAS_V2
#ifdef BLAS_V2
	#include <cublas_v2.h>
#else
	#include <cublas.h>
#endif


//#define FieldParticle
//#define REAL4

#ifndef REAL4
//	#define	eps				5e-21
//	#define	eps				5e-31
	#define	eps				5e-99
	#define	REAL 				double
	#define	REAL2 			double2
	#ifdef CUDACOMPLEX_H
		#define	CMPLX			doublecomplex
		#define	mkCMPLX		make_cudacomplex
	#else
		#define	CMPLX			cuDoubleComplex
		#define	mkCMPLX		make_cuDoubleComplex
	#endif
	#define	SIZEZ				16
	#define	SIZED				 8
#else
	#define	eps				5e-16
	#define	REAL 				float
	#define	REAL2 			float2
	#ifdef CUDACOMPLEX_H
		#define	CMPLX			singlecomplex
		#define	mkCMPLX		make_cudacomplex
	#else
		#define	CMPLX			cuDoubleComplex
		#define	mkCMPLX		make_cuComplex
	#endif
	#define	SIZEZ				8
	#define	SIZED				4
	#define	cuCadd			cuCaddf
	#define	cuCsub			cuCsubf
	#define	cuCmul			cuCmulf
	#define	cuCdiv			cuCdivf
	#define	cuCabs			cuCabsf
	#define	cuConj			cuConjf
	#define	cuCreal			cuCrealf
	#define	cuCimag			cuCimagf
	#define	cublasDasum	cublasSasum
	#define	cublasDznrm2	cublasScnrm2
#endif

#include "constant.cu"
#include "func.cu"
#include "FDM.cu"
#include "potential.cu"
#include "EXPECT.cu"
#include "write.cu"
#include "PHI_SOR.cu"
#include "timer.c"
#include "restart.c"

template < class T > T* host_allocate(size_t size)
{//return static_cast<T*>(malloc(sizeof(T)*size));
	T *ptr; cudaHostAlloc( &ptr, sizeof(T) *size, cudaHostAllocDefault );
	return ptr;
}

template < class T > void host_release(T* ptr) {	/*	free(ptr);*/	cudaFreeHost(ptr);}

//----------------------------------------------------------------------------
//------------------------------main------------------------------------------
//----------------------------------------------------------------------------
int main ( int argc, char** argv) // Program to solve the Schrodinger equation
//{	cutilDeviceInit(argc,argv);		//	cudaSetDevice(1);
{	int devID = findCudaDevice(argc, (const char **)argv); printf("devID =%d\n",devID);
	std::string filename = (argc>=2) ? argv[1] : "sch.inp";

	StopWatchInterface *timer = 0;
	SetTimer(&timer);
	#ifdef BLAS_V2
		cublasStatus_t stat;
		cublasHandle_t handle;
		stat = cublasCreate(&handle);
		if(stat != CUBLAS_STATUS_SUCCESS)
		{	printf("CUBLAS initialization failed\n"); return EXIT_FAILURE; }
	#endif
	int N_rec		;
	int N_rec_rho	;
	int N_step		;
	REAL hx     	;
	REAL hy     	;
	REAL Lx      	;
	REAL Ly      	;
	REAL omg    	;
	REAL v0x    	;
	REAL v0y    	;
	REAL Bz			;
	REAL m			;
	REAL q			;
	REAL x0		   ;
	REAL y0		   ;
	REAL vSIZE	   ;
	REAL v0_	   ;
	REAL Bz_    	;
	REAL mp_	   ;
	REAL e_     	;
	REAL h_bar_ 	;

	if( !read_parameter(filename,
		N_rec, N_rec_rho, N_step, hx, hy, Lx, Ly,
		omg, v0x, v0y, Bz, m, q, x0, y0,
		vSIZE, v0_, Bz_, mp_, e_, h_bar_) )	return -1;

	REAL tau    = 5 * 2 * pi / N_step;
//	REAL tau    =     2 * pi / N_step;
	REAL hh     = hx * hy;
	REAL rl_    = mp_ * v0_ / (e_ * Bz_);
	REAL sgm    = sqrt( h_bar_ / (e_ * Bz_) ) * sqrt( 1 / ( q * Bz ) ) / rl_;
	REAL sgm2   = Square( sgm );
	REAL h_bar  = h_bar_ / (mp_ * v0_ * rl_);
	//REAL k0x    = m * v0x / h_bar;
	REAL k0x    =(m * v0x - q*Bz*y0) / h_bar;
	REAL k0y    = m * v0y / h_bar;
	int Nx 		= int( Lx / hx );
	int Ny 		= int( Ly / hy );
	int Nx2		= 2 * Nx + 1;
	int Ny2		= 2 * Ny + 1;
	int NG		= Nx2 * Ny2;
	//double N_t_s    =     N_step * m / (q * Bz) * 1050 / 100;
	//double N_t_s    =     N_step * m / (q * Bz) * 105 / 100 * (0.5 / 5);	//打ち切る用
	double N_t_s    =     N_step * m / (q * Bz) * 105 / 100 * (5 / 5);	// 5周用
	int    N_t_step = int(N_t_s);
	//N_t_step = 40;	// 打ちきり用
	int k_start = 1;

	CMPLX mvx_av, mvy_av, Px_av;
	CMPLX var_mvx, var_mvy, var_mv, var_Px, var_Py, var_P;
	REAL x_av, y_av, var_x, var_y, var_r;
	REAL qAx_av;
	REAL E0, EN;
	REAL err;
	FILE *interrupt;

	printf("magnetic length = %e\n", sgm);
	printf("Nx2 = %d Ny2 = %d hx = %e hy = %e NG = %d\n", Nx2, Ny2, hx, hy, NG);
	printf("h_bar = %e\n", h_bar);

	#ifndef FieldParticle
		int s3 = sizeof(CMPLX)+sizeof(CMPLX)+sizeof(REAL);	// sizeof Psi + Phi + Res
															// s3=20(REAL4),s3=40(double)
	#else
		int s3 = sizeof(CMPLX)+sizeof(CMPLX)+sizeof(REAL)+sizeof(REAL);	// sizeof Psi + Phi + Res + u
																		// s3=24(REAL4),s3=48(double)
	#endif
	int Nb = int( vSIZE / ( s3 * Nx2 ) - 2 );	// 一度に計算できるだろう行数
	printf("Nb = %d\n", Nb);
	int Nb2 = ( Ny2 - 2 ) / Nb;
	if ( Nb2 == 0 ) Nb = Ny2 - 2;	// Nb2 = 0 ならば一度に全部計算できる
	int Nb_ = Nb + 2;			// デバイスに一度に送る行数（基本）
	int Nb3_= Ny2 - Nb * Nb2;	// デバイスに最後に送る行数（余りの行）
	int Nb3 = Nb3_ - 2;			// 最後に計算する行数（余りの行）

	#ifdef BLAS_V2
		REAL var_tmp = 0;
	#endif

	printf("Ny2= %d Nb= %d Nb2= %d Nb3= %d Nb_= %d\n", Ny2, Nb, Nb2, Nb3, Nb_);
	// Ny2=全行数, Nb=一度に計算する行数

	//----------------------------------------------------------------------------
	REAL *x =(REAL*) host_allocate<REAL>(Nx2);
	REAL *y =(REAL*) host_allocate<REAL>(Ny2);
	for (int ix = 0; ix < Nx2; ix++){ x[ix] = hx * (ix-Nx); }
	for (int iy = 0; iy < Ny2; iy++){ y[iy] = hy * (iy-Ny); }

	#ifdef FieldParticle
		REAL *u = (REAL*) host_allocate< REAL> (NG); if( u == 0 ) { printf("%d Error\n", SIZED * NG); return 0; }
		REAL Eps0 = 8.8542e-12 / Square(e_) * mp_ * Square(v0_) * rl_;	// 真空の誘電率
		printf("Eps0= %e\n", Eps0);
		PotEne( u, x, y, q, sgm, Eps0, Nx2, Ny2 );
		pot( Nx2, Ny2, x, y, u );
	#endif

	CMPLX *Psi = (CMPLX*) host_allocate<CMPLX>(NG); if( Psi == 0 ) { printf("%d Error\n", SIZEZ * NG); return 0; }
	printf("host Psi\n");
	CMPLX *Phi = (CMPLX*) host_allocate<CMPLX>(NG); if( Phi == 0 ) { printf("%d Error\n", SIZEZ * NG); return 0; }
	printf("host Phi\n");

	//---------- LOAD Wavefunction ----------
	if( load_wavefunction(N_t_step, NG, &k_start, Psi, &E0) == 0 ) // 真/偽=1/0
	{	// *.binファイルが見つからない
		// Initialize the wavefunction *********************
		REAL Norm_Psi = 1 / ( sgm * sqrt(pi) );
		REAL Sum_Psi = 0;
		for(	int iy = 0; iy < Ny2; iy++) { REAL dy = y[iy] - y0; REAL dy2 = Square(dy);
			for(int ix = 0; ix < Nx2; ix++) { REAL dx = x[ix] - x0; REAL dx2 = Square(dx);
				int II = ix + Nx2 * iy;
				REAL aaa = Norm_Psi * exp( - ( dx2 + dy2 ) / ( 2 * sgm2 ) );
			//	REAL aaa = Norm_Psi * exp( - (       dy2 ) / ( 2 * sgm2 ) );
				REAL bbb = k0x * x[ix] + k0y * y[iy];
				#ifdef CUDACOMPLEX_H
				Psi[II].real() = aaa * cos( bbb );
				Psi[II].imag() = aaa * sin( bbb );
				Sum_Psi  += Square( Psi[II].real() ) + Square( Psi[II].imag() );
				#else
				Psi[II].x = aaa * cos( bbb );
				Psi[II].y = aaa * sin( bbb );
				Sum_Psi  += Square( cuCabs(Psi[II]) );
				#endif
			}
		}
		Sum_Psi *= hh;	printf("Rho     =%20.16f\n", Sum_Psi );
		// 自然境界条件	begin
		#ifdef CUDACOMPLEX_H
		printf("psi[0]={%11.3e,%11.3e} before BC\n", Psi[0].real(), Psi[0].imag());
		#else
		printf("psi[0]={%11.3e,%11.3e} before BC\n", cuCreal(Psi[0]), cuCimag(Psi[0]));
		#endif
		for( int iy = 0; iy < Ny2  ; iy++)
		{	int II, ix;
			#ifdef CUDACOMPLEX_H
			ix =     0; II = ix + Nx2 * iy; Psi[II] = Psi[II+1] + Psi[II+1] - Psi[II+2];
			ix = Nx2-1; II = ix + Nx2 * iy; Psi[II] = Psi[II-1] + Psi[II-1] - Psi[II-2];
			#else
			ix =     0; II = ix + Nx2 * iy; Psi[II] = cuCsub( cuCadd( Psi[II+1], Psi[II+1] ), Psi[II+2] );
			ix = Nx2-1; II = ix + Nx2 * iy; Psi[II] = cuCsub( cuCadd( Psi[II-1], Psi[II-1] ), Psi[II-2] );
			#endif
		}
		for( int ix = 1; ix < Nx2-1; ix++)
		{	int II, iy;
			#ifdef CUDACOMPLEX_H
			iy =     0; II = ix + Nx2 * iy; Psi[II] = Psi[II+Nx2] + Psi[II+Nx2] - Psi[II+2*Nx2];
			iy = Ny2-1; II = ix + Nx2 * iy; Psi[II] = Psi[II-Nx2] + Psi[II-Nx2] - Psi[II-2*Nx2];
			#else
			iy =     0; II = ix + Nx2 * iy; Psi[II] = cuCsub( cuCadd( Psi[II+Nx2], Psi[II+Nx2] ), Psi[II+2*Nx2] );
			iy = Ny2-1; II = ix + Nx2 * iy; Psi[II] = cuCsub( cuCadd( Psi[II-Nx2], Psi[II-Nx2] ), Psi[II-2*Nx2] );
			#endif
		}
		#ifdef CUDACOMPLEX_H
		printf("psi[0]={%11.3e,%11.3e}  after BC\n", Psi[0].real(), Psi[0].imag());
		#else
		printf("psi[0]={%11.3e,%11.3e}  after BC\n", cuCreal(Psi[0]), cuCimag(Psi[0]));
		#endif
		// 自然境界条件	end
		REAL Sum_Psi_ = 0;
		for(	int iy = 0; iy < Ny2; iy++)
		{	for(int ix = 0; ix < Nx2; ix++)
			{	int II = ix + Nx2 * iy;
				#ifdef CUDACOMPLEX_H
				Psi[II]   = Psi[II] / sqrt(Sum_Psi);
				Sum_Psi_ += Square( Psi[II].real() ) + Square( Psi[II].imag() );
				#else
				Psi[II].x /= sqrt(Sum_Psi);
				Psi[II].y /= sqrt(Sum_Psi);
				Sum_Psi_  += Square( cuCabs(Psi[II]) );
				#endif
			}
		}
		Sum_Psi_ *= hh;	printf("Rho_bar =%20.16f\n", Sum_Psi_);
	/*
		// cublasの使い方例
		{	REAL *s = (REAL*) malloc( SIZED * NG );
			for(	int iy = 0; iy < Ny2; iy++)
			{	for(int ix = 0; ix < Nx2; ix++)
				{	int II = ix + Nx2 * iy;
					s[II] = Square( cuCabs(Psi[II]) ) * hh;
				}
			}
			cublasStatus stat;
			REAL *S;
			cublasInit(); // printf("sizeof(*s) = %d\n", sizeof(*s));
			stat = cublasAlloc ( NG, SIZED, (void**) &S );	// printf ("%s\n", stat);
			if ( stat != CUBLAS_STATUS_SUCCESS) { printf ("device memory allocation failed"); return 1; }
			cublasSetVector( NG, SIZED, s, 1, S, 1);	//printf ("%s\n", cublasGetError());
			printf("blasasum=%20.16f, blas error = %s\n", cublasDasum( NG, S, 1), cublasGetError());
			cublasFree(S);
			cublasShutdown();
			free(s);
		}
	*/
		//-------------init judge---------------------------------------------------
		EXPECT( Psi, x, y,
		#ifdef FieldParticle
		                   u,
		#endif
				&mvx_av, &mvy_av, &Px_av, &E0, &var_mv, &var_mvx, &var_mvy, &var_P, &var_Px, &var_Py,
				&x_av, &y_av, &var_r, &var_x, &var_y, &err, &qAx_av,
				Bz, m, q, Nx2, Ny2, hx, hy, h_bar);

		#ifdef CUDACOMPLEX_H
		printf("<mvx> = %24.16e\n", mvx_av.real());
		printf("<qAx> = %24.16e\n", qAx_av       );
		printf("<P_x> = %24.16e\n",  Px_av.real());
		#else
		printf("<mvx> = %24.16e\n", mvx_av.x     );
		printf("<qAx> = %24.16e\n", qAx_av       );
		printf("<P_x> = %24.16e\n",  Px_av.x     );
		#endif

		EN = E0;
		Record( 0, tau, E0, EN, var_mvx, var_mvy, var_mv, mvx_av, mvy_av, var_Px, var_Py, var_P, x_av, y_av, var_x, var_y, var_r, Px_av, err);
		rho(0, Nx2, Ny2, hh, x, y, Psi);

	//	OUTに一初めの~~~.rhoを書いたところで止めるためのreturn 0;
	//	return 0;
	//
		printf("      k itr     <x>     <y>  var[r]     <u>     <v>  var[p]           <Px>          EN       EN-E0       error\n");
		printf("%7d %3d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %14.6e %11.3e %11.3e %11.3e\n",
		#ifdef CUDACOMPLEX_H
				 0, 0, x_av, y_av, var_r, mvx_av.real(), mvy_av.real(), var_mv.real(), Px_av.real(), E0, EN-E0, err);
		#else
				 0, 0, x_av, y_av, var_r, mvx_av.x     , mvy_av.x     , var_mv.x     , Px_av.x     , E0, EN-E0, err);
		#endif
		/********** SAVE Initial Wavefunction **********/
		printf("Save Initial Wavefunction, k = %d\n", 0);
		save_wavefunction(NG, 0, Psi, E0);
		/***********************************************/
	}

	#ifdef FieldParticle
	REAL     *u_device; cudaMalloc( (void**) &u_device,   SIZED * Nx2 * Nb_ );
	#endif
	CMPLX  *Psi_device; cudaMalloc( (void**) &Psi_device, SIZEZ * Nx2 * Nb_ );
	CMPLX  *Phi_device; cudaMalloc( (void**) &Phi_device, SIZEZ * Nx2 * Nb_ );
	REAL   *Res_device; cudaMalloc( (void**) &Res_device, SIZED * Nx2 * Nb_ );

//	int THREADS =  224;
//	int THREADS =  256;
	int THREADS =  512;
//	int THREADS = 1024;	// for GT-480
	int GRIDS = Nx2 * Nb_ / THREADS; if ( Nx2 * Nb_ % THREADS != 0 ) GRIDS++;
	dim3 grid( GRIDS );
	dim3 threads(THREADS);

	//*************** Set coeffcient ***********************
	#ifdef CUDACOMPLEX_H
	CMPLX beta  =        { 0, tau * h_bar / ( 4 * m ) };
	CMPLX gamma =        {    tau / ( 2 * m ),    0   };
	CMPLX zeta  =        { 0, tau / ( 4 * m * h_bar ) };
	CMPLX xi    =        { 0, tau / ( 2 *     h_bar ) };
	#else
	CMPLX beta  = mkCMPLX( 0, tau * h_bar / ( 4 * m ) );
	CMPLX gamma = mkCMPLX(    tau / ( 2 * m ),    0   );
	CMPLX zeta  = mkCMPLX( 0, tau / ( 4 * m * h_bar ) );
	CMPLX xi    = mkCMPLX( 0, tau / ( 2 *     h_bar ) );
	#endif
	//*************** end Set coeffcient *******************
	SetConst( Nb, Nx2, Ny2, Nb3_, Nx, Ny, sgm2, pi, hx, hy, hh, omg, x0, y0, Bz, q, beta, gamma, zeta, xi);
	//
	if( Nb2 == 0 ) {	// すべてビデオメモリに収まる
		cudaMemcpy( Psi_device, Psi, SIZEZ * NG, cudaMemcpyHostToDevice);
		#ifdef FieldParticle
		cudaMemcpy(   u_device,   u, SIZED * NG, cudaMemcpyHostToDevice);
		#endif
	}
	//
//	cublasInit();
//	******* Total computing loop *************************
	printf("k: %d -> %d\n", k_start, N_t_step);
	int k;
	for(k = k_start; k <= N_t_step; k++)
	{		for(	int k2 = 1; k2 <= Nb2; k2++)	// すべてビデオメモリに収まらない
		{	int N_field = k2 - 1;
			cudaMemcpy ( Psi_device, Psi + Nx2 * Nb * N_field, SIZEZ * Nx2 * Nb_ , cudaMemcpyHostToDevice);
			#ifdef FieldParticle
			cudaMemcpy (   u_device,   u + Nx2 * Nb * N_field, SIZED * Nx2 * Nb_ , cudaMemcpyHostToDevice);
			PHI <<< grid, threads >>> ( Psi_device, Phi_device, u_device, N_field );
			#else
			PHI <<< grid, threads >>> ( Psi_device, Phi_device,           N_field );
			#endif
			cudaMemcpy ( Phi+ Nx2 * Nb * N_field + Nx2, Phi_device + Nx2, SIZEZ * Nx2 * Nb, cudaMemcpyDeviceToHost);
		}
		if( Nb3 > 0 )
		{	int N_field = Nb2;
			if ( Nb2 != 0 ) {
				cudaMemcpy ( Psi_device, Psi + Nx2 * Nb * N_field, SIZEZ * Nx2 * Nb3_, cudaMemcpyHostToDevice);
				#ifdef FieldParticle
				cudaMemcpy (   u_device,   u + Nx2 * Nb * N_field, SIZED * Nx2 * Nb3_, cudaMemcpyHostToDevice);
				#endif
			}
			#ifdef FieldParticle
				PHI <<< grid, threads >>> ( Psi_device, Phi_device, u_device, N_field );
			#else
				PHI <<< grid, threads >>> ( Psi_device, Phi_device,           N_field );
			#endif
cudaThreadSynchronize();
			if ( Nb2 != 0 )
				cudaMemcpy ( Phi + Nx2 * Nb * N_field + Nx2, Phi_device + Nx2, SIZEZ * Nx2 * Nb3, cudaMemcpyDeviceToHost);
		}

		//---------start while--------------------------------------------------------
		REAL var = 2 * eps, var0 = 1e99; int itr = 0;
		while( var > eps )
		{	itr++;	var = 0;
			for(int k2 = 1; k2 <= Nb2; k2++)
			{	// 領域の分割がある場合
				int N_field = k2 - 1;
				cudaMemcpy ( Phi_device, Phi + Nx2 * Nb * N_field, SIZEZ * Nx2 * Nb_, cudaMemcpyHostToDevice);
				cudaMemcpy ( Psi_device, Psi + Nx2 * Nb * N_field, SIZEZ * Nx2 * Nb_, cudaMemcpyHostToDevice);
				#ifdef FieldParticle
				cudaMemcpy ( u_device, u + Nx2 * Nb * N_field, SIZED * Nx2 * Nb_, cudaMemcpyHostToDevice);
				#endif
				if (1)	/* 1 <- 0, on May 18, 2015: oikawa */
				{// 1
					#ifdef FieldParticle
						SOR <<< grid, threads >>> ( Psi_device, Phi_device, Res_device, u_device, N_field );
					#else
						SOR <<< grid, threads >>> ( Psi_device, Phi_device, Res_device, N_field );
					#endif
cudaThreadSynchronize();
				}
				else
				{// 0		/*Does not work on host-b (CUDA 7.0) on May 18, 2015: oikawa */
					int BLOCKS = ( Ny2 % THREADS == 0  ) ? Ny2 / THREADS : Ny2 / THREADS  + 1;
					#ifdef FieldParticle
						sor0<<< BLOCKS, THREADS >>> ( Psi_device, Phi_device, Res_device, u_device, N_field );
cudaThreadSynchronize();
						sor1<<< BLOCKS, THREADS >>> ( Psi_device, Phi_device, Res_device, u_device, N_field );
					#else
						sor0<<< BLOCKS, THREADS >>> ( Psi_device, Phi_device, Res_device,           N_field );
cudaThreadSynchronize();
						sor1<<< BLOCKS, THREADS >>> ( Psi_device, Phi_device, Res_device,           N_field );
					#endif
cudaThreadSynchronize();
				}
				cudaMemcpy ( Psi + Nx2 * Nb * N_field + Nx2, Psi_device + Nx2, SIZEZ * Nx2 * Nb , cudaMemcpyDeviceToHost);
				#ifdef BLAS_V2
					stat = cublasDasum(handle, Nx2 * Nb, Res_device + Nx2, 1, &var_tmp);
					var += var_tmp;
				#else
					var += cublasDasum ( Nx2 * Nb , Res_device + Nx2, 1) ;
				#endif
			}
			if( Nb3 > 0 )
			{	int N_field = Nb2;
				if ( Nb2 != 0 )
				{	( cudaMemcpy ( Phi_device, Phi + Nx2 * Nb * N_field, SIZEZ * Nx2 * Nb3_, cudaMemcpyHostToDevice));
					( cudaMemcpy ( Psi_device, Psi + Nx2 * Nb * N_field, SIZEZ * Nx2 * Nb3_, cudaMemcpyHostToDevice));
					#ifdef FieldParticle
						( cudaMemcpy ( u_device, u + Nx2 * Nb * N_field, SIZED * Nx2 * Nb3_, cudaMemcpyHostToDevice));
					#endif
				}
				if (1)
				{// 1
					#ifdef FieldParticle
						SOR <<< grid, threads >>> ( Psi_device, Phi_device, Res_device, u_device, N_field );
					#else
						SOR <<< grid, threads >>> ( Psi_device, Phi_device, Res_device, N_field );
					#endif
cudaThreadSynchronize();
				}
				else
				{// 0
					int BLOCKS = ( Ny2 % THREADS == 0  ) ? Ny2 / THREADS : Ny2 / THREADS  + 1;
					#ifdef FieldParticle
						sor0<<< BLOCKS, THREADS >>> ( Psi_device, Phi_device, Res_device, u_device, N_field );
cudaThreadSynchronize();
						sor1<<< BLOCKS, THREADS >>> ( Psi_device, Phi_device, Res_device, u_device, N_field );
					#else
						sor0<<< BLOCKS, THREADS >>> ( Psi_device, Phi_device, Res_device,           N_field );
cudaThreadSynchronize();
						sor1<<< BLOCKS, THREADS >>> ( Psi_device, Phi_device, Res_device,           N_field );
cudaThreadSynchronize();
					#endif
				}
				if ( Nb2 != 0 )
					cudaMemcpy ( Psi + Nx2 * Nb * N_field + Nx2, Psi_device + Nx2, SIZEZ * Nx2 * Nb3, cudaMemcpyDeviceToHost );
				// begin BLAS
				#ifdef BLAS_V2
					stat = cublasDasum (handle, Nx2 * Nb3 , Res_device + Nx2, 1, &var_tmp) ;
					var += var_tmp;
				#else
					var += cublasDasum ( Nx2 * Nb3 , Res_device + Nx2, 1) ;
				#endif
				//	end BLAS
				cudaThreadSynchronize();
			}
			var /= NG;
			int NotNumber = isnan(var);
			if ( itr == 99 || NotNumber != 0 )
			{	printf("+++++ Abnormal end +++++\nOVER itr = %d, NAN = %d\n", itr, NotNumber);
				/***** 終了処理 **************************************************************/
				fclose(interrupt);	//	中断チェック用のFILE*変数をクローズする
				cudaFree(Psi_device); cudaFree(Phi_device); cudaFree(Res_device);
				#ifdef FieldParticle
				cudaFree(u_device);
				#endif
				cudaThreadExit();
				// -----------------------------------------------------------------------------
				host_release(x); host_release(y); host_release(Psi); host_release(Phi);
				#ifdef FieldParticle
				host_release(u);
				#endif
				/*****************************************************************************/
				exit(1);
			}
			cudaThreadSynchronize();
			//	***** for test	*******************************
			//	if( NotNumber == 0 ) printf("var = %e\n", var);
			//	******************************************************
			//printf("k=%3d, itr = %2d, var = %12.3e\n", k, itr, var);
			if( var > var0 / 2 )	break;
			var0 = var;
		}
		//----- end while ---------------------------------------------------------------------

		if (  Nb2     == 0 ) ( cudaMemcpy ( Psi, Psi_device, SIZEZ * NG, cudaMemcpyDeviceToHost) );
		// 毎回必ず、デバイスからホストへPsiを戻す
		/********** 期待値らを出力する **********/
		if( k % N_rec == 0 ) {
			EXPECT( Psi, x, y,
			#ifdef FieldParticle
			                   u,
			#endif
					&mvx_av, &mvy_av, &Px_av, &EN, &var_mv, &var_mvx, &var_mvy, &var_P, &var_Px, &var_Py,
					&x_av, &y_av, &var_r, &var_x, &var_y, &err, &qAx_av,
					Bz, m, q, Nx2, Ny2, hx, hy, h_bar);
			Record( k, tau, E0, EN, var_mvx, var_mvy, var_mv, mvx_av, mvy_av, var_Px, var_Py, var_P, x_av, y_av, var_x, var_y, var_r, Px_av, err);
			printf("%7d %3d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %14.6e %11.3e %11.3e %11.3e\n",
			#ifdef CUDACOMPLEX_H
					 k, itr, x_av, y_av, var_r, mvx_av.real(), mvy_av.real(), var_mv.real(), Px_av.real(), EN, EN-E0, err);
			#else
					 k, itr, x_av, y_av, var_r, mvx_av.x, mvy_av.x, var_mv.x, Px_av.x, EN, EN-E0, err);
			#endif
			if ( k == N_t_step )
				printf("      k itr     <x>     <y>  var[r]     <u>     <v>  var[p]           <Px>          EN       EN-E0       error\n");

			// ここにあったらN_rec_rhoはN_recの倍数でないと一切出力されないのでは?
		}	//printf("k=%7d, itr = %2d, var0 = %12.3e, var = %12.3e\n", k, itr, var0, var);
		/********** 粒子保存の確認 **********/
		if( k % 20 == 0 )		// 確認が終わったら\nを\rにする
		{	if( Nb2 == 0 )
			{	//	領域分割なし
				#ifdef BLAS_V2
					REAL Dznrm2_tmp;
					stat = cublasDznrm2(handle, NG,
					#ifdef CUDACOMPLEX_H
						(cuDoubleComplex*)
					#endif
					Psi_device, 1, &Dznrm2_tmp);
					printf("k=%7d, itr=%3d, residue=%11.3e, error using blasnrm2=%11.3e\n",
						     k, itr, var, Square(Dznrm2_tmp) * hh - 1 );
				#else
//					printf("blasasum=%11.3e\n", Square(cublasDasum(2*NG,(REAL*)Psi_device, 1)) * hh - 1 );
					printf("k=%7d, itr=%3d, residue=%11.3e, error using blasnrm2=%11.3e\n",
						     k, itr, var, Square(cublasDznrm2(NG,
					#ifdef CUDACOMPLEX_H
						 		(cuDoubleComplex*)
					#endif
					    		Psi_device, 1))* hh - 1 );
				#endif
				cudaThreadSynchronize();
				/*
				printf("k=%7d, itr=%3d, residue=%11.3e, error using blasdot =%25.16e\n",
					    k, itr, var, cublasDdot ( 2*NG, (REAL*)Psi_device, 1, (REAL*)Psi_device, 1)*hh-1 );
				cudaThreadSynchronize();
				*/
			}
			else
			{	// 領域分割あり
				// ************ cublas使って粒子保存を見る場合はこっち。
				REAL blasnrm2 = 0;	//REAL blasdot  = 0;
				#ifdef BLAS_V2
					REAL blasDznrm2_tmp1 = 0;
					stat = cublasDznrm2( handle, Nx2 * (Nb3 + 1),
					#ifdef CUDACOMPLEX_H
						(cuDoubleComplex*)
					#endif
					(Psi_device + Nx2), 1, &blasDznrm2_tmp1);
					blasnrm2 = blasDznrm2_tmp1 * blasDznrm2_tmp1;
				#else
					blasnrm2 = Square(cublasDznrm2( Nx2 * (Nb3 + 1),
					#ifdef CUDACOMPLEX_H
						(cuDoubleComplex*)
					#endif
					(Psi_device + Nx2), 1) ); // 最後のブロック＋最後の行
					/*
					#ifdef CUDACOMPLEX_H
					blasnrm2 = Square(cublasDznrm2( Nx2 * (Nb3 + 1), (cuDoubleComplex*)(Psi_device + Nx2), 1) ); // 最後のブロック＋最後の行
					#else
					blasnrm2 = Square(cublasDznrm2( Nx2 * (Nb3 + 1),                    Psi_device + Nx2 , 1) ); // 最後のブロック＋最後の行
					#endif
					*/
				#endif
				/*
				cudaThreadSynchronize();
				blasdot  =      cublasDdot( 2 * Nx2 * (Nb3 + 1), (REAL*)(Psi_device + Nx2), 1, (REAL*)(Psi_device + Nx2), 1);
				cudaThreadSynchronize();
				*/
				cudaMemcpy ( Psi_device, Psi, SIZEZ * Nx2 * Nb_, cudaMemcpyHostToDevice);
				#ifdef BLAS_V2
					REAL blasDznrm2_tmp2 = 0;
					stat = cublasDznrm2( handle, Nx2 * (Nb + 1),
					#ifdef CUDACOMPLEX_H
						(cuDoubleComplex*)
					#endif
					(Psi_device + Nx2), 1, &blasDznrm2_tmp2);
					blasnrm2 += blasDznrm2_tmp2 * blasDznrm2_tmp2;
					#else
					blasnrm2 += Square(cublasDznrm2( Nx2 * (Nb + 1),
					#ifdef CUDACOMPLEX_H
						(cuDoubleComplex*)
					#endif
					Psi_device, 1));		 // 最初のブロック＋最初の行
				/*
				#ifdef CUDACOMPLEX_H
				blasnrm2 += Square(cublasDznrm2( Nx2 * (Nb + 1), (cuDoubleComplex*)Psi_device, 1));		 // 最初のブロック＋最初の行
				#else
				blasnrm2 += Square(cublasDznrm2( Nx2 * (Nb + 1),                   Psi_device, 1));		 // 最初のブロック＋最初の行
				#endif
				*/
				#endif
				cudaThreadSynchronize();
				/*
				blasdot  +=        cublasDdot( 2 * Nx2 * (Nb + 1), (REAL*)Psi_device, 1, (REAL*)Psi_device, 1);
				cudaThreadSynchronize();
				*/
				for(int k2 = 2; k2 <= Nb2; k2++)
				{	int N_field = k2 - 1;
					cudaMemcpy ( Psi_device, Psi + Nx2 * Nb * N_field, SIZEZ * Nx2 * Nb_, cudaMemcpyHostToDevice);
					#ifdef BLAS_V2
						REAL blasDznrm2_tmp3 = 0;
						stat = cublasDznrm2( handle,  Nx2 * Nb,
						#ifdef CUDACOMPLEX_H
							(cuDoubleComplex*)
						#endif
						(Psi_device + Nx2), 1, &blasDznrm2_tmp3);	// 中間のブロック
						blasnrm2 += blasDznrm2_tmp3 * blasDznrm2_tmp3;
					#else
						blasnrm2 += Square(cublasDznrm2( Nx2 * Nb,
						#ifdef CUDACOMPLEX_H
							(cuDoubleComplex*)
						#endif
						(Psi_device + Nx2), 1));	// 中間のブロック
						/*
						#ifdef CUDACOMPLEX_H
						blasnrm2 += Square(cublasDznrm2( Nx2 * Nb, (cuDoubleComplex*)(Psi_device + Nx2), 1));	// 中間のブロック
						#else
						blasnrm2 += Square(cublasDznrm2( Nx2 * Nb, Psi_device + Nx2, 1));	// 中間のブロック
						#endif
						*/
					#endif
					cudaThreadSynchronize();
					/*
					blasdot  +=      cublasDdot( 2 * Nx2 * Nb, (REAL*)(Psi_device + Nx2), 1, (REAL*)(Psi_device + Nx2), 1);
					cudaThreadSynchronize();
					*/
				}
				printf("k=%7d, itr=%3d, blasnrm2=%11.3e\n", k, itr, blasnrm2 * hh - 1);
				//printf("k=%7d, itr=%3d, blasdot =%25.16e\n", k, itr, blasdot  * hh - 1);

//				************ cublas使いたくなくて、どうしても粒子保存見たければ使う。
/*				REAL asum = 0;
				for(		int iy = 0; iy < Ny2; iy++)
				{	for(	int ix = 0; ix < Nx2; ix++)
					{	int II = ix + Nx2 * iy;
						#ifdef CUDACOMPLEX_H
						asum += Square( Psi[II].real() ) + Square( Psi[II].imag() );
						#else
						asum += Square( cuCabs( Psi[II] ) );
						#endif
					}
				}
				printf("k=%7d, itr=%3d, asum=%25.16e\r", k, itr, asum * hh - 1);
*/
			}
		}
		/********** ディレクトリOUTに精度を落とした確率密度関数を出力する **********/
		if ( k % N_rec_rho == 0 )
		{	rho ( k, Nx2, Ny2, hh, x, y, Psi);
			/********** Interruption **********/
			if( (interrupt = fopen("SAVE/interruption", "r")) != NULL )
			{	printf(" Interruption, k = %d\n", k); // このタイミングなら今回の波動関数
				fclose(interrupt);	//	中断チェック用のFILE*変数をクローズする
				break;
			}
			/**********************************/
		}
		/********** 10500ごとに定期的に波動関数を出力させる **********/
		// 典型的なパラメータで105000は5周程度、10500は半周くらい
		if(k % 10500 == 0)
		{	printf("Save Wavefunction, k = %d\n", k);
			save_wavefunction(NG, k, Psi, E0);
		}
		/********************************************************/
		fflush(stdout);
	}	printf("Elapsed TIME = %f sec.\n", EndTimer(&timer) / 1000 );
	/*****************************************************************************/
	/***** 終了処理 **************************************************************/
	/*****************************************************************************/
	//---------- SAVE last wavefunction --------------
	//printf("*** Save Last wavefunction ***\n");
	//	普通にループを終了したならばk==N_t_step+1となっているので補正する
	if( k == N_t_step+1 ) k--;
	save_wavefunction(NG, k, Psi, E0);
	//---------- END SAVE last wavefunction ----------
	#ifdef BLAS_V2
		fprintf(stdout, "cublasDestroy(blas4.0)\n"); cublasDestroy(handle);
	#endif
	/* free allocation (GPU) */
	cudaFree(Psi_device); cudaFree(Phi_device); cudaFree(Res_device);
	#ifdef FieldParticle
		cudaFree(u_device);
	#endif
	cudaThreadExit();
	// -----------------------------------------------------------------------------
	/* free allocation (CPU) */
	host_release(x); host_release(y); host_release(Psi); host_release(Phi);
	#ifdef FieldParticle
		host_release(u);
	#endif
	return 0;
}
