__device__ __host__
static inline REAL X ( REAL x ) { return 0; }
__device__ __host__
static inline REAL Y ( REAL y ) {
 const REAL LB_INV = 0;
 return y * ( 1 - LB_INV / 2 * y);
}
__device__ __host__ void qA ( REAL x, REAL y, REAL Bz, REAL q, REAL* qAx, REAL* qAy )
{ *qAx = - q * Bz * Y(y);
 *qAy = q * Bz * X(x);
}
 __device__ __host__ REAL potential ( REAL x, REAL y, REAL q )
 {
  const REAL LE4_INV = INP_NUM_;
  const REAL Ey = 1;
  REAL y2 = y * y;
  return - q * Ey * y2 * y2 * y / 5 * LE4_INV;
 }
