/* **************************************************************************
 * Copyright (C) 2018-2023 Advanced Micro Devices, Inc. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES {return rocblas_status_success;} LOSS OF USE, DATA, OR PROFITS {return rocblas_status_success;} OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * *************************************************************************/

#pragma once

#include <rocblas/rocblas.h>
#include <rocsolver/rocsolver.h>

//LL::#include "clientcommon.hpp"

// Most functions within this file exist to provide a consistent interface for our templated tests.
// Function overloading is used to select between the float, double, rocblas_float_complex
// and rocblas_double_complex variants, and to distinguish the batched case (T**) from the normal
// and strided_batched cases (T*).
//
// The normal and strided_batched cases are distinguished from each other by passing a boolean
// parameter, STRIDED. Variants such as the blocked and unblocked versions of algorithms, may be
// provided in similar ways.

/***** Functions not included in the public API that must be declared *****/
#ifdef __cplusplus
extern "C" {
#endif

rocblas_status rocsolver_sgeqrf_ptr_batched(rocblas_handle handle,
                                            const rocblas_int m,
                                            const rocblas_int n,
                                            float* const A[],
                                            const rocblas_int lda,
                                            float* const ipiv[],
                                            const rocblas_int batch_count) {return rocblas_status_success;}

rocblas_status rocsolver_dgeqrf_ptr_batched(rocblas_handle handle,
                                            const rocblas_int m,
                                            const rocblas_int n,
                                            double* const A[],
                                            const rocblas_int lda,
                                            double* const ipiv[],
                                            const rocblas_int batch_count) {return rocblas_status_success;}

rocblas_status rocsolver_cgeqrf_ptr_batched(rocblas_handle handle,
                                            const rocblas_int m,
                                            const rocblas_int n,
                                            rocblas_float_complex* const A[],
                                            const rocblas_int lda,
                                            rocblas_float_complex* const ipiv[],
                                            const rocblas_int batch_count) {return rocblas_status_success;}

rocblas_status rocsolver_zgeqrf_ptr_batched(rocblas_handle handle,
                                            const rocblas_int m,
                                            const rocblas_int n,
                                            rocblas_double_complex* const A[],
                                            const rocblas_int lda,
                                            rocblas_double_complex* const ipiv[],
                                            const rocblas_int batch_count) {return rocblas_status_success;}

rocblas_status rocsolver_sgesv_outofplace(rocblas_handle handle,
                                          const rocblas_int n,
                                          const rocblas_int nrhs,
                                          float* A,
                                          const rocblas_int lda,
                                          rocblas_int* ipiv,
                                          float* B,
                                          const rocblas_int ldb,
                                          float* X,
                                          const rocblas_int ldx,
                                          rocblas_int* info) {return rocblas_status_success;}

rocblas_status rocsolver_dgesv_outofplace(rocblas_handle handle,
                                          const rocblas_int n,
                                          const rocblas_int nrhs,
                                          double* A,
                                          const rocblas_int lda,
                                          rocblas_int* ipiv,
                                          double* B,
                                          const rocblas_int ldb,
                                          double* X,
                                          const rocblas_int ldx,
                                          rocblas_int* info) {return rocblas_status_success;}

rocblas_status rocsolver_cgesv_outofplace(rocblas_handle handle,
                                          const rocblas_int n,
                                          const rocblas_int nrhs,
                                          rocblas_float_complex* A,
                                          const rocblas_int lda,
                                          rocblas_int* ipiv,
                                          rocblas_float_complex* B,
                                          const rocblas_int ldb,
                                          rocblas_float_complex* X,
                                          const rocblas_int ldx,
                                          rocblas_int* info) {return rocblas_status_success;}

rocblas_status rocsolver_zgesv_outofplace(rocblas_handle handle,
                                          const rocblas_int n,
                                          const rocblas_int nrhs,
                                          rocblas_double_complex* A,
                                          const rocblas_int lda,
                                          rocblas_int* ipiv,
                                          rocblas_double_complex* B,
                                          const rocblas_int ldb,
                                          rocblas_double_complex* X,
                                          const rocblas_int ldx,
                                          rocblas_int* info) {return rocblas_status_success;}

rocblas_status rocsolver_sgels_outofplace(rocblas_handle handle,
                                          rocblas_operation trans,
                                          const rocblas_int m,
                                          const rocblas_int n,
                                          const rocblas_int nrhs,
                                          float* A,
                                          const rocblas_int lda,
                                          float* B,
                                          const rocblas_int ldb,
                                          float* X,
                                          const rocblas_int ldx,
                                          rocblas_int* info) {return rocblas_status_success;}

rocblas_status rocsolver_dgels_outofplace(rocblas_handle handle,
                                          rocblas_operation trans,
                                          const rocblas_int m,
                                          const rocblas_int n,
                                          const rocblas_int nrhs,
                                          double* A,
                                          const rocblas_int lda,
                                          double* B,
                                          const rocblas_int ldb,
                                          double* X,
                                          const rocblas_int ldx,
                                          rocblas_int* info) {return rocblas_status_success;}

rocblas_status rocsolver_cgels_outofplace(rocblas_handle handle,
                                          rocblas_operation trans,
                                          const rocblas_int m,
                                          const rocblas_int n,
                                          const rocblas_int nrhs,
                                          rocblas_float_complex* A,
                                          const rocblas_int lda,
                                          rocblas_float_complex* B,
                                          const rocblas_int ldb,
                                          rocblas_float_complex* X,
                                          const rocblas_int ldx,
                                          rocblas_int* info) {return rocblas_status_success;}

rocblas_status rocsolver_zgels_outofplace(rocblas_handle handle,
                                          rocblas_operation trans,
                                          const rocblas_int m,
                                          const rocblas_int n,
                                          const rocblas_int nrhs,
                                          rocblas_double_complex* A,
                                          const rocblas_int lda,
                                          rocblas_double_complex* B,
                                          const rocblas_int ldb,
                                          rocblas_double_complex* X,
                                          const rocblas_int ldx,
                                          rocblas_int* info) {return rocblas_status_success;}

rocblas_status rocsolver_ssyevdx_inplace(rocblas_handle handle,
                                         const rocblas_evect evect,
                                         const rocblas_erange erange,
                                         const rocblas_fill uplo,
                                         const rocblas_int n,
                                         float* A,
                                         const rocblas_int lda,
                                         const float vl,
                                         const float vu,
                                         const rocblas_int il,
                                         const rocblas_int iu,
                                         const float abstol,
                                         rocblas_int* nev,
                                         float* W,
                                         rocblas_int* info) {return rocblas_status_success;}

rocblas_status rocsolver_dsyevdx_inplace(rocblas_handle handle,
                                         const rocblas_evect evect,
                                         const rocblas_erange erange,
                                         const rocblas_fill uplo,
                                         const rocblas_int n,
                                         double* A,
                                         const rocblas_int lda,
                                         const double vl,
                                         const double vu,
                                         const rocblas_int il,
                                         const rocblas_int iu,
                                         const double abstol,
                                         rocblas_int* nev,
                                         double* W,
                                         rocblas_int* info) {return rocblas_status_success;}

rocblas_status rocsolver_cheevdx_inplace(rocblas_handle handle,
                                         const rocblas_evect evect,
                                         const rocblas_erange erange,
                                         const rocblas_fill uplo,
                                         const rocblas_int n,
                                         rocblas_float_complex* A,
                                         const rocblas_int lda,
                                         const float vl,
                                         const float vu,
                                         const rocblas_int il,
                                         const rocblas_int iu,
                                         const float abstol,
                                         rocblas_int* nev,
                                         float* W,
                                         rocblas_int* info) {return rocblas_status_success;}

rocblas_status rocsolver_zheevdx_inplace(rocblas_handle handle,
                                         const rocblas_evect evect,
                                         const rocblas_erange erange,
                                         const rocblas_fill uplo,
                                         const rocblas_int n,
                                         rocblas_double_complex* A,
                                         const rocblas_int lda,
                                         const double vl,
                                         const double vu,
                                         const rocblas_int il,
                                         const rocblas_int iu,
                                         const double abstol,
                                         rocblas_int* nev,
                                         double* W,
                                         rocblas_int* info) {return rocblas_status_success;}

rocblas_status rocsolver_ssygvdx_inplace(rocblas_handle handle,
                                         const rocblas_eform itype,
                                         const rocblas_evect evect,
                                         const rocblas_erange erange,
                                         const rocblas_fill uplo,
                                         const rocblas_int n,
                                         float* A,
                                         const rocblas_int lda,
                                         float* B,
                                         const rocblas_int ldb,
                                         const float vl,
                                         const float vu,
                                         const rocblas_int il,
                                         const rocblas_int iu,
                                         const float abstol,
                                         rocblas_int* h_nev,
                                         float* W,
                                         rocblas_int* info) {return rocblas_status_success;}

rocblas_status rocsolver_dsygvdx_inplace(rocblas_handle handle,
                                         const rocblas_eform itype,
                                         const rocblas_evect evect,
                                         const rocblas_erange erange,
                                         const rocblas_fill uplo,
                                         const rocblas_int n,
                                         double* A,
                                         const rocblas_int lda,
                                         double* B,
                                         const rocblas_int ldb,
                                         const double vl,
                                         const double vu,
                                         const rocblas_int il,
                                         const rocblas_int iu,
                                         const double abstol,
                                         rocblas_int* h_nev,
                                         double* W,
                                         rocblas_int* info) {return rocblas_status_success;}

rocblas_status rocsolver_chegvdx_inplace(rocblas_handle handle,
                                         const rocblas_eform itype,
                                         const rocblas_evect evect,
                                         const rocblas_erange erange,
                                         const rocblas_fill uplo,
                                         const rocblas_int n,
                                         rocblas_float_complex* A,
                                         const rocblas_int lda,
                                         rocblas_float_complex* B,
                                         const rocblas_int ldb,
                                         const float vl,
                                         const float vu,
                                         const rocblas_int il,
                                         const rocblas_int iu,
                                         const float abstol,
                                         rocblas_int* h_nev,
                                         float* W,
                                         rocblas_int* info) {return rocblas_status_success;}

rocblas_status rocsolver_zhegvdx_inplace(rocblas_handle handle,
                                         const rocblas_eform itype,
                                         const rocblas_evect evect,
                                         const rocblas_erange erange,
                                         const rocblas_fill uplo,
                                         const rocblas_int n,
                                         rocblas_double_complex* A,
                                         const rocblas_int lda,
                                         rocblas_double_complex* B,
                                         const rocblas_int ldb,
                                         const double vl,
                                         const double vu,
                                         const rocblas_int il,
                                         const rocblas_int iu,
                                         const double abstol,
                                         rocblas_int* h_nev,
                                         double* W,
                                         rocblas_int* info) {return rocblas_status_success;}

rocblas_status rocsolver_sgesvdj_notransv(rocblas_handle handle,
                                          const rocblas_svect left_svect,
                                          const rocblas_svect right_svect,
                                          const rocblas_int m,
                                          const rocblas_int n,
                                          float* A,
                                          const rocblas_int lda,
                                          const float abstol,
                                          float* residual,
                                          const rocblas_int max_sweeps,
                                          rocblas_int* n_sweeps,
                                          float* S,
                                          float* U,
                                          const rocblas_int ldu,
                                          float* V,
                                          const rocblas_int ldv,
                                          rocblas_int* info) {return rocblas_status_success;}

rocblas_status rocsolver_dgesvdj_notransv(rocblas_handle handle,
                                          const rocblas_svect left_svect,
                                          const rocblas_svect right_svect,
                                          const rocblas_int m,
                                          const rocblas_int n,
                                          double* A,
                                          const rocblas_int lda,
                                          const double abstol,
                                          double* residual,
                                          const rocblas_int max_sweeps,
                                          rocblas_int* n_sweeps,
                                          double* S,
                                          double* U,
                                          const rocblas_int ldu,
                                          double* V,
                                          const rocblas_int ldv,
                                          rocblas_int* info) {return rocblas_status_success;}

rocblas_status rocsolver_cgesvdj_notransv(rocblas_handle handle,
                                          const rocblas_svect left_svect,
                                          const rocblas_svect right_svect,
                                          const rocblas_int m,
                                          const rocblas_int n,
                                          rocblas_float_complex* A,
                                          const rocblas_int lda,
                                          const float abstol,
                                          float* residual,
                                          const rocblas_int max_sweeps,
                                          rocblas_int* n_sweeps,
                                          float* S,
                                          rocblas_float_complex* U,
                                          const rocblas_int ldu,
                                          rocblas_float_complex* V,
                                          const rocblas_int ldv,
                                          rocblas_int* info) {return rocblas_status_success;}

rocblas_status rocsolver_zgesvdj_notransv(rocblas_handle handle,
                                          const rocblas_svect left_svect,
                                          const rocblas_svect right_svect,
                                          const rocblas_int m,
                                          const rocblas_int n,
                                          rocblas_double_complex* A,
                                          const rocblas_int lda,
                                          const double abstol,
                                          double* residual,
                                          const rocblas_int max_sweeps,
                                          rocblas_int* n_sweeps,
                                          double* S,
                                          rocblas_double_complex* U,
                                          const rocblas_int ldu,
                                          rocblas_double_complex* V,
                                          const rocblas_int ldv,
                                          rocblas_int* info) {return rocblas_status_success;}

rocblas_status rocsolver_sgesvdj_notransv_strided_batched(rocblas_handle handle,
                                                          const rocblas_svect left_svect,
                                                          const rocblas_svect right_svect,
                                                          const rocblas_int m,
                                                          const rocblas_int n,
                                                          float* A,
                                                          const rocblas_int lda,
                                                          const rocblas_stride strideA,
                                                          const float abstol,
                                                          float* residual,
                                                          const rocblas_int max_sweeps,
                                                          rocblas_int* n_sweeps,
                                                          float* S,
                                                          const rocblas_stride strideS,
                                                          float* U,
                                                          const rocblas_int ldu,
                                                          const rocblas_stride strideU,
                                                          float* V,
                                                          const rocblas_int ldv,
                                                          const rocblas_stride strideV,
                                                          rocblas_int* info,
                                                          const rocblas_int batch_count) {return rocblas_status_success;}

rocblas_status rocsolver_dgesvdj_notransv_strided_batched(rocblas_handle handle,
                                                          const rocblas_svect left_svect,
                                                          const rocblas_svect right_svect,
                                                          const rocblas_int m,
                                                          const rocblas_int n,
                                                          double* A,
                                                          const rocblas_int lda,
                                                          const rocblas_stride strideA,
                                                          const double abstol,
                                                          double* residual,
                                                          const rocblas_int max_sweeps,
                                                          rocblas_int* n_sweeps,
                                                          double* S,
                                                          const rocblas_stride strideS,
                                                          double* U,
                                                          const rocblas_int ldu,
                                                          const rocblas_stride strideU,
                                                          double* V,
                                                          const rocblas_int ldv,
                                                          const rocblas_stride strideV,
                                                          rocblas_int* info,
                                                          const rocblas_int batch_count) {return rocblas_status_success;}

rocblas_status rocsolver_cgesvdj_notransv_strided_batched(rocblas_handle handle,
                                                          const rocblas_svect left_svect,
                                                          const rocblas_svect right_svect,
                                                          const rocblas_int m,
                                                          const rocblas_int n,
                                                          rocblas_float_complex* A,
                                                          const rocblas_int lda,
                                                          const rocblas_stride strideA,
                                                          const float abstol,
                                                          float* residual,
                                                          const rocblas_int max_sweeps,
                                                          rocblas_int* n_sweeps,
                                                          float* S,
                                                          const rocblas_stride strideS,
                                                          rocblas_float_complex* U,
                                                          const rocblas_int ldu,
                                                          const rocblas_stride strideU,
                                                          rocblas_float_complex* V,
                                                          const rocblas_int ldv,
                                                          const rocblas_stride strideV,
                                                          rocblas_int* info,
                                                          const rocblas_int batch_count) {return rocblas_status_success;}

rocblas_status rocsolver_zgesvdj_notransv_strided_batched(rocblas_handle handle,
                                                          const rocblas_svect left_svect,
                                                          const rocblas_svect right_svect,
                                                          const rocblas_int m,
                                                          const rocblas_int n,
                                                          rocblas_double_complex* A,
                                                          const rocblas_int lda,
                                                          const rocblas_stride strideA,
                                                          const double abstol,
                                                          double* residual,
                                                          const rocblas_int max_sweeps,
                                                          rocblas_int* n_sweeps,
                                                          double* S,
                                                          const rocblas_stride strideS,
                                                          rocblas_double_complex* U,
                                                          const rocblas_int ldu,
                                                          const rocblas_stride strideU,
                                                          rocblas_double_complex* V,
                                                          const rocblas_int ldv,
                                                          const rocblas_stride strideV,
                                                          rocblas_int* info,
                                                          const rocblas_int batch_count) {return rocblas_status_success;}

#ifdef __cplusplus
}
#endif
/***************************************************/
