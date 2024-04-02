#include <cuda_runtime.h>
#include <cusolverDn.h>

#include <iostream>
#include <random>
#include <vector>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>

#define NA 512
#define PR_D 0

constexpr int error_exit_code = -1;

inline int report_validation_result(int errors)
{
    if(errors)
    {
        std::cout << "Validation failed. Errors: " << errors << std::endl;
        return error_exit_code;
    }

    std::cout << "Validation passed." << std::endl;
    return 0;
}

template<typename T>
void multiply_matrices(T        alpha,
                       T        beta,
                       int      m,
                       int      n,
                       int      k,
                       const T* A,
                       int      stride1_a,
                       int      stride2_a,
                       const T* B,
                       int      stride1_b,
                       int      stride2_b,
                       T*       C,
                       int      stride_c)
{
    for(int i1 = 0; i1 < m; ++i1)
    {
        for(int i2 = 0; i2 < n; ++i2)
        {
            T t = T(0.0);
            for(int i3 = 0; i3 < k; ++i3)
            {
                t += A[i1 * stride1_a + i3 * stride2_a] * B[i3 * stride1_b + i2 * stride2_b];
            }
            C[i1 + i2 * stride_c] = beta * C[i1 + i2 * stride_c] + alpha * t;
        }
    }
}

template<class BidirectionalIterator>
inline std::string format_range(const BidirectionalIterator begin, const BidirectionalIterator end)
{
    std::stringstream sstream;
    sstream << "[ ";
    for(auto it = begin; it != end; ++it)
    {
        sstream << *it;
        if(it != std::prev(end))
        {
            sstream << ", ";
        }
    }
    sstream << " ]";
    return sstream.str();
}


inline const char* cusolverStatusToString(cusolverStatus_t status)
{
    switch(status)
    {
    case CUSOLVER_STATUS_SUCCESS:
      return "CUSOLVER_STATUS_SUCCESS";
    case CUSOLVER_STATUS_NOT_INITIALIZED:
      return "CUSOLVER_STATUS_NOT_INITIALIZED";
    case CUSOLVER_STATUS_ALLOC_FAILED:
      return "CUSOLVER_STATUS_ALLOC_FAILED";
    case CUSOLVER_STATUS_INVALID_VALUE:
      return "CUSOLVER_STATUS_INVALID_VALUE";
    case CUSOLVER_STATUS_ARCH_MISMATCH:
      return "CUSOLVER_STATUS_ARCH_MISMATCH";
    case CUSOLVER_STATUS_MAPPING_ERROR:
      return "CUSOLVER_STATUS_MAPPING_ERROR";
    case CUSOLVER_STATUS_EXECUTION_FAILED:
      return "CUSOLVER_STATUS_EXECUTION_FAILED";
    case CUSOLVER_STATUS_INTERNAL_ERROR:
      return "CUSOLVER_STATUS_INTERNAL_ERROR";
    case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
      return "CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
    case CUSOLVER_STATUS_NOT_SUPPORTED:
      return "CUSOLVER_STATUS_NOT_SUPPORTED ";
    case CUSOLVER_STATUS_ZERO_PIVOT:
      return "CUSOLVER_STATUS_ZERO_PIVOT";
    case CUSOLVER_STATUS_INVALID_LICENSE:
      return "CUSOLVER_STATUS_INVALID_LICENSE";
    }

    return "<unknown>";
}

#define CUDA_CHECK(condition)                                                                \
    {                                                                                       \
        const cudaError_t error = condition;                                                 \
        if(error != cudaSuccess)                                                             \
        {                                                                                   \
            std::cerr << "An error encountered: \"" << cudaGetErrorString(error) << "\" at " \
                      << __FILE__ << ':' << __LINE__ << std::endl;                          \
            std::exit(error_exit_code);                                                     \
        }                                                                                   \
    }

#define CUSOLVER_CHECK(condition)                                                            \
    {                                                                                         \
        const cusolverStatus_t status = condition;                                           \
        if(status != CUSOLVER_STATUS_SUCCESS)                                                \
        {                                                                                     \
            std::cerr << "cuSOLVER error encountered: \"" << cusolverStatusToString(status) \
                      << "\" at " << __FILE__ << ':' << __LINE__ << std::endl;                \
            std::exit(error_exit_code);                                                       \
        }                                                                                     \
    }

void init_matrix(std::vector<double> &A, int n, int lda)
{
    std::default_random_engine             generator;
    std::uniform_real_distribution<double> distribution(0., 2.);
    auto                                   random_number = std::bind(distribution, generator);

    for(int i = 0; i < n; i++)
    {
        A[(lda + 1) * i] = random_number();
        for(int j = 0; j < i; j++)
        {
            A[i * lda + j] = A[j * lda + i] = random_number();
        }
    }
}

int main(const int argc, char* argv[])
{
    int n = NA;
    if(n <= 0)
    {
        std::cout << "Value of 'n' should be greater than 0" << std::endl;
        return 0;
    }
    const int lda = n;

    // 2. Data vectors
    std::vector<double> A(n * lda); // Input matrix
    std::vector<double> V(n * lda); // Resulting eigenvectors
    std::vector<double> W(n); // resulting eigenvalues

    // 3. Generate a random symmetric matrix
    init_matrix(A, n, lda);

    // 4. Set hipsolver parameters
    const cusolverEigMode_t  jobz = CUSOLVER_EIG_MODE_VECTOR;
    const cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;

    syevjInfo_t  params;
    CUSOLVER_CHECK(cusolverDnCreateSyevjInfo(&params));
    CUSOLVER_CHECK(cusolverDnXsyevjSetMaxSweeps(params, 100));
    CUSOLVER_CHECK(cusolverDnXsyevjSetTolerance(params, 1.e-7));
    CUSOLVER_CHECK(cusolverDnXsyevjSetSortEig(params, 1));

    // 5. Reserve and copy data to device
    double* d_A    = nullptr;
    double* d_W    = nullptr;
    int*    d_info = nullptr;

    CUDA_CHECK(cudaMalloc(&d_A, sizeof(double) * A.size()));
    CUDA_CHECK(cudaMalloc(&d_W, sizeof(double) * W.size()));
    CUDA_CHECK(cudaMalloc(&d_info, sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_A, A.data(), sizeof(double) * A.size(), cudaMemcpyHostToDevice));

    // 6. Initialize hipsolver
    cusolverDnHandle_t cusolver_handle;
    CUSOLVER_CHECK(cusolverDnCreate(&cusolver_handle));

    // 7. Get and reserve the working space on device.
    int     lwork  = 0;
    double* d_work = nullptr;
    CUSOLVER_CHECK(
        cusolverDnDsyevj_bufferSize(cusolver_handle, jobz, uplo, n, d_A, lda, d_W, &lwork, params));


    std::cout<< "LL:: 1 lwork = "<<lwork<<"bytes"<<std::endl;
    lwork += 64;
//    lwork = ((lwork+64-1)/64)*64;

    std::cout<< "LL:: 2 lwork = "<<lwork<<"bytes"<<std::endl;

    CUDA_CHECK(cudaMalloc(&d_work, sizeof(double) * lwork));

    // 8. Compute eigenvectors and eigenvalues
//    HIPSOLVER_CHECK(hipsolverDsyevj(hipsolver_handle,
    CUSOLVER_CHECK(cusolverDnDsyevj(cusolver_handle,
                                    jobz,
                                    uplo,
                                    n,
                                    d_A,
                                    lda,
                                    d_W,
                                    d_work,
                                    lwork,
                                    d_info,
                                    params));
    // 9. Get results from host.
    int info = 0;
    CUDA_CHECK(cudaMemcpy(V.data(), d_A, sizeof(double) * V.size(), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(W.data(), d_W, sizeof(double) * W.size(), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost));

    // 10. Print results
    if(info == 0)
    {
        std::cout << "SYEVJ converges." << std::endl;
    }
    else if(info > 0)
    {
        std::cout << "SYEVJ does not converge (" << info << " elements did not converge)."
                  << std::endl;
    }

    std::cout << "\nGiven the n x n square input matrix A; we computed the linearly independent "
                 "eigenvectors V and the associated eigenvalues W."
              << std::endl;
#if PR_D
    std::cout << "A = " << format_range(A.begin(), A.end()) << std::endl;
    std::cout << "W = " << format_range(W.begin(), W.end()) << std::endl;
    std::cout << "V = " << format_range(V.begin(), V.end()) << std::endl;
#endif

    int    sweeps   = 0;
    double residual = 0;
    CUSOLVER_CHECK(cusolverDnXsyevjGetSweeps(cusolver_handle, params, &sweeps));
    CUSOLVER_CHECK(cusolverDnXsyevjGetResidual(cusolver_handle, params, &residual));

    std::cout << "\nWhich was computed in " << sweeps << " sweeps, with a residual of " << residual
              << std::endl;

    // 11. Validate that 'AV == VD' and that 'AV - VD == 0'.
    std::cout << "\nLet D be the diagonal constructed from W.\n"
              << "The right multiplication of A * V should result in V * D [AV == VD]:"
              << std::endl;

    // Right multiplication of the input matrix with the eigenvectors.
    std::vector<double> AV(n * lda);
    multiply_matrices(1.0, 0.0, n, n, n, A.data(), lda, 1, V.data(), 1, lda, AV.data(), lda);
#if PR_D
    std::cout << "AV = " << format_range(AV.begin(), AV.end()) << std::endl;
#endif
    // Construct the diagonal D from eigenvalues W.
    std::vector<double> D(n * n);
    for(int i = 0; i < n; i++)
    {
        D[(n + 1) * i] = W[i];
    }

    // Scale eigenvectors V with W by multiplying V with D.
    std::vector<double> VD(n * lda);
    multiply_matrices(1.0, 0.0, n, n, n, V.data(), 1, lda, D.data(), lda, 1, VD.data(), lda);
#if PR_D
    std::cout << "VD = " << format_range(VD.begin(), VD.end()) << std::endl;
#endif
    double epsilon = 1.0e5 * std::numeric_limits<double>::epsilon();
    int    errors  = 0;
    double mse     = 0;
    for(int i = 0; i < n * n; i++)
    {
        double diff = (AV[i] - VD[i]);
        diff *= diff;
        mse += diff;

        errors += (diff > epsilon);
    }
    mse /= n * n;
    std::cout << "\nMean Square Error of [AV == VD]:\n  " << mse << std::endl;

    // 12. Clean up device allocations.
    CUSOLVER_CHECK(cusolverDnDestroy(cusolver_handle));
    CUSOLVER_CHECK(cusolverDnDestroySyevjInfo(params));
    CUDA_CHECK(cudaFree(d_A));
    CUDA_CHECK(cudaFree(d_W));
    CUDA_CHECK(cudaFree(d_work));
    CUDA_CHECK(cudaFree(d_info));

    return report_validation_result(errors);
}
