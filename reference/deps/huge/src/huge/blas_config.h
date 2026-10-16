// blas_config.h — Thin BLAS abstraction for hot-path linear algebra.
//
// When built as an R package, R provides BLAS via $(BLAS_LIBS).
// When standalone, link against system BLAS (OpenBLAS, Accelerate, etc).
// Fortran-style BLAS symbols use trailing underscore on most platforms.
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// y := alpha * A * x + beta * y  (column-major, no-transpose)
// DGEMV: op(A) is m-by-n, x has length n, y has length m
void dgemv_(const char* trans, const int* m, const int* n,
            const double* alpha, const double* A, const int* lda,
            const double* x, const int* incx,
            const double* beta, double* y, const int* incy);

// C := alpha * op(A) * op(B) + beta * C  (column-major)
// DGEMM: op(A) is m-by-k, op(B) is k-by-n, C is m-by-n
void dgemm_(const char* transa, const char* transb,
            const int* m, const int* n, const int* k,
            const double* alpha, const double* A, const int* lda,
            const double* B, const int* ldb,
            const double* beta, double* C, const int* ldc);

// dot product: sum(x[i]*y[i])
double ddot_(const int* n, const double* x, const int* incx,
             const double* y, const int* incy);

// y := alpha * x + y
void daxpy_(const int* n, const double* alpha,
            const double* x, const int* incx,
            double* y, const int* incy);

// x := alpha * x
void dscal_(const int* n, const double* alpha,
            double* x, const int* incx);

#ifdef __cplusplus
}
#endif
