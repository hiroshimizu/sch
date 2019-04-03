#pragma once  // ヘッダファイルにはこれを書く

__device__ __host__ inline REAL Square ( REAL x) { REAL y; y = x * x; return y; }

