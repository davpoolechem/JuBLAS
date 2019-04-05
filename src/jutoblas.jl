include("dgemm.jl")
import LinearAlgebra.BLAS

function dgemm_check(input::Int64)
    matrix_size::Int64 = input
    matrix_a::Array{Float64,2} = rand(matrix_size,matrix_size)
    matrix_b::Array{Float64,2} = rand(matrix_size,matrix_size)

    BLAS.gemm('N', 'N', 1.0, matrix_a, matrix_b)

    m::Int64 = size(matrix_a)[1]
    n::Int64 = size(matrix_b)[2]
    k::Int64 = size(matrix_a)[2]

    lda::Int64 = max(1,m)
    ldb::Int64 = max(1,k)
    ldc::Int64 = max(1,m)

    matrix_c::Array{Float64,2} = zeros(ldc, n)

    JuBLAS.dgemm('N', 'N', m, n, k,
        1.0, matrix_a, lda,
        matrix_b, ldb, 0.0,
        matrix_c, ldc)
end
