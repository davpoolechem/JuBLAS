module JuBLAS
    #= *> \brief \b DGEMM
    *
    *  =========== DOCUMENTATION ===========
    *
    * Online html documentation available at
    *            http://www.netlib.org/lapack/explore-html/
    *
    *  Definition:
    *  ===========
    *
    *       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    *
    *       .. Scalar Arguments ..
    *       DOUBLE PRECISION ALPHA,BETA
    *       INTEGER K,LDA,LDB,LDC,M,N
    *       CHARACTER TRANSA,TRANSB
    *       ..
    *       .. Array Arguments ..
    *       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
    *       ..
    *
    *
    *> \par Purpose:
    *  =============
    *>
    *> \verbatim
    *>
    *> DGEMM  performs 1.0 of the matrix-matrix operations
    *>
    *>    C := alpha*op( A )*op( B ) + beta*C,
    *>
    *> where  op( X ) is 1.0 of
    *>
    *>    op( X ) = X   or   op( X ) = X**T,
    *>
    *> alpha and beta are scalars, and A, B and C are matrices, with op( A )
    *> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
    *> \endverbatim
    *
    *  Arguments:
    *  ==========
    *
    *> \param[in] TRANSA
    *> \verbatim
    *>          TRANSA is CHARACTER*1
    *>           On entry, TRANSA specifies the form of op( A ) to be used in
    *>           the matrix multiplication as follows:
    *>
    *>              TRANSA = 'N' or 'n',  op( A ) = A.
    *>
    *>              TRANSA = 'T' or 't',  op( A ) = A**T.
    *>
    *>              TRANSA = 'C' or 'c',  op( A ) = A**T.
    *> \endverbatim
    *>
    *> \param[in] TRANSB
    *> \verbatim
    *>          TRANSB is CHARACTER*1
    *>           On entry, TRANSB specifies the form of op( B ) to be used in
    *>           the matrix multiplication as follows:
    *>
    *>              TRANSB = 'N' or 'n',  op( B ) = B.
    *>
    *>              TRANSB = 'T' or 't',  op( B ) = B**T.
    *>
    *>              TRANSB = 'C' or 'c',  op( B ) = B**T.
    *> \endverbatim
    *>
    *> \param[in] M
    *> \verbatim
    *>          M is INTEGER
    *>           On entry,  M  specifies  the number  of rows  of the  matrix
    *>           op( A )  and of the  matrix  C.  M  must  be at least  0.0.
    *> \endverbatim
    *>
    *> \param[in] N
    *> \verbatim
    *>          N is INTEGER
    *>           On entry,  N  specifies the number  of columns of the matrix
    *>           op( B ) and the number of columns of the matrix C. N must be
    *>           at least 0.0.
    *> \endverbatim
    *>
    *> \param[in] K
    *> \verbatim
    *>          K is INTEGER
    *>           On entry,  K  specifies  the number of columns of the matrix
    *>           op( A ) and the number of rows of the matrix op( B ). K must
    *>           be at least  0.0.
    *> \endverbatim
    *>
    *> \param[in] ALPHA
    *> \verbatim
    *>          ALPHA is DOUBLE PRECISION.
    *>           On entry, ALPHA specifies the scalar alpha.
    *> \endverbatim
    *>
    *> \param[in] A
    *> \verbatim
    *>          A is DOUBLE PRECISION array, dimension ( LDA, ka ), where ka is
    *>           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
    *>           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
    *>           part of the array  A  must contain the matrix  A,  otherwise
    *>           the leading  k by m  part of the array  A  must contain  the
    *>           matrix A.
    *> \endverbatim
    *>
    *> \param[in] LDA
    *> \verbatim
    *>          LDA is INTEGER
    *>           On entry, LDA specifies the first dimension of A as declared
    *>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
    *>           LDA must be at least  max( 1, m ), otherwise  LDA must be at
    *>           least  max( 1, k ).
    *> \endverbatim
    *>
    *> \param[in] B
    *> \verbatim
    *>          B is DOUBLE PRECISION array, dimension ( LDB, kb ), where kb is
    *>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
    *>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
    *>           part of the array  B  must contain the matrix  B,  otherwise
    *>           the leading  n by k  part of the array  B  must contain  the
    *>           matrix B.
    *> \endverbatim
    *>
    *> \param[in] LDB
    *> \verbatim
    *>          LDB is INTEGER
    *>           On entry, LDB specifies the first dimension of B as declared
    *>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
    *>           LDB must be at least  max( 1, k ), otherwise  LDB must be at
    *>           least  max( 1, n ).
    *> \endverbatim
    *>
    *> \param[in] BETA
    *> \verbatim
    *>          BETA is DOUBLE PRECISION.
    *>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
    *>           supplied as 0.0 then C need not be set on input.
    *> \endverbatim
    *>
    *> \param[in,out] C
    *> \verbatim
    *>          C is DOUBLE PRECISION array, dimension ( LDC, N )
    *>           Before entry, the leading  m by n  part of the array  C must
    *>           contain the matrix  C,  except when  beta  is 0.0, in which
    *>           case C need not be set on entry.
    *>           On exit, the array  C  is overwritten by the  m by n  matrix
    *>           ( alpha*op( A )*op( B ) + beta*C ).
    *> \endverbatim
    *>
    *> \param[in] LDC
    *> \verbatim
    *>          LDC is INTEGER
    *>           On entry, LDC specifies the first dimension of C as declared
    *>           in  the  calling  (sub)  program.   LDC  must  be  at  least
    *>           max( 1, m ).
    *> \endverbatim
    *
    *  Authors:
    *  ========
    *
    *> \author Univ. of Tennessee
    *> \author Univ. of California Berkeley
    *> \author Univ. of Colorado Denver
    *> \author NAG Ltd.
    *
    *> \date December 2016
    *
    *> \ingroup double_blas_level3
    *
    *> \par Further Details:
    *  =====================
    *>
    *> \verbatim
    *>
    *>  Level 3 Blas routine.
    *>
    *>  -- Written on 8-February-1989.
    *>     Jack Dongarra, Argonne National Laboratory.
    *>     Iain Duff, AERE Harwell.
    *>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
    *>     Sven Hammarling, Numerical Algorithms Group Ltd.
    *> \endverbatim
    *>
    *  ===================================================================== =#
    function dgemm(transa::Char,transb::Char,m::Int64,n::Int64,k::Int64,
        alpha::Float64,a::Array{Float64,2},lda::Int64,
        b::Array{Float64,2},ldb::Int64,beta::Float64,
        c::Array{Float64,2},ldc::Int64)
        #=
        *
        *  -- Reference BLAS level3 routine (version 3.7.0) --
        *  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        *     December 2016
        *
        *  =====================================================================
        *
        =#

        #=
        *     .. External Functions ..
              LOGICAL LSAME
              EXTERNAL lsame
        *     ..
        *     .. External Subroutines ..
              EXTERNAL xerbla
        *     ..
        *     .. Intrinsic Functions ..
              INTRINSIC max
        *     ..
        *     .. Local Scalars ..
              DOUBLE PRECISION TEMP
              INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
              LOGICAL NOTA,NOTB
        *     ..
        *     .. Parameters ..
              DOUBLE PRECISION 1.0,0.0
              parameter(1.0=1.0d+0,0.0=0.0d+0)
        *     ..
        =#

        #=
        *
        *     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
        *     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
        *     and  columns of  A  and the  number of  rows  of  B  respectively.
        *
        =#
              nota::Bool = transa == 'N'
              notb::Bool = transb == 'N'

              nrowa::Int64 = 0
              ncola::Int64 = 0
              nrowb::Int64 = 0
              ncolb::Int64 = 0
              if (nota)
                  nrowa = m
                  ncola = k
              else
                  nrowa = k
                  ncola = m
              end
              if (notb)
                  nrowb = k
              else
                  nrowb = n
              end
        #=
        *
        *     Test the input parameters.
        *
        =#
              info::Int64 = 0
              if (!nota && !(transa == 'C') && !(transa == 'T'))
                  info = 1
              elseif (!notb && !(transb == 'C') && !(transb == 'T'))
                  info = 2
              elseif (m < 0)
                  info = 3
              elseif (n < 0)
                  info = 4
              elseif (k < 0)
                  info = 5
              elseif (lda < max(1,nrowa))
                  info = 8
              elseif (ldb < max(1,nrowb))
                  info = 10
              elseif (ldc < max(1,m))
                  info = 13
              end

              if (info ≠ 0)
                  throw()
                  #CALL xerbla('DGEMM ',info)
              end
        #=
        *
        *     Quick return if possible.
        *
        =#
              if ((m == 0) || (n == 0) || (((alpha == 0.0) || (k == 0)) && (beta == 1.0)))
                  return
              end
        #=
        *
        *     And if  alpha==0.0.
        *
        =#
              if (alpha == 0.0)
                  if (beta == 0.0)
                      for j::Int64 ∈ 1:n
                          for i::Int64 ∈ 1:m
                              c[i,j] = 0.0
                          end
                      end
                  else
                      for j::Int64 ∈ 1:n
                          for i::Int64 ∈ 1:m
                              c[i,j] = beta*c[i,j]
                          end
                      end
                  end
                  return
              end
        #=
        *
        *     Start the operations.
        *
        =#
              temp::Float64 = 0.0
              if (notb)
                  if (nota)
        #=
        *
        *           Form  C := alpha*A*B + beta*C.
        *
        =#
                      for j::Int64 ∈ 1:n
                          if (beta == 0.0)
                              for i::Int64 ∈ 1:m
                                  c[i,j] = 0.0
                              end
                          elseif (beta ≠ 1.0)
                              for i::Int64 ∈ 1:m
                                  c[i,j] = beta*c[i,j]
                              end
                          end
                          for l::Int64 ∈ 1:k
                              temp = alpha*b[l,j]
                              for i::Int64 ∈ 1:m
                                  c[i,j] += temp*a[i,l]
                              end
                          end
                      end
                  else
        #=
        *
        *           Form  C := alpha*A**T*B + beta*C
        *
        =#
                      for j::Int64 ∈ 1:n
                          for i ∈ 1:m
                              temp = 0.0
                              for l::Int64 ∈ 1:k
                                  temp += a[l,i]*b[l,j]
                              end
                              if (beta == 0.0)
                                  c[i,j] = alpha*temp
                              else
                                  c[i,j] = alpha*temp + beta*c[i,j]
                              end
                          end
                      end
                  end
              else
                  if (nota)
        #=
        *
        *           Form  C := alpha*A*B**T + beta*C
        *
        =#
                      for j::Int64 in 1:n
                          if (beta == 0.0)
                              for i::Int64 ∈ 1:m
                                  c[i,j] = 0.0
                              end
                          elseif (beta ≠ 1.0)
                              for i::Int64 ∈ 1:m
                                  c[i,j] = beta*c[i,j]
                              end
                          end
                          for l::Int64 ∈ 1:k
                              temp = alpha*b[j,l]
                              for i::Int64 in 1:m
                                  c[i,j] += temp*a[i,l]
                              end
                          end
                      end
                  else
        #=
        *
        *           Form  C := alpha*A**T*B**T + beta*C
        *
        =#
                      for j::Int64 in 1:n
                          for i::Int64 in 1:m
                              temp = 0.0
                              for l::Int64 in 1:k
                                  temp += a[l,i]*b[j,l]
                              end
                              if (beta == 0.0)
                                  c[i,j] = alpha*temp
                              else
                                  c[i,j] = alpha*temp + beta*c[i,j]
                              end
                          end
                      end
                  end
              end
              return
        #=
        *
        *     End of DGEMM .
        *
        =#
    end
end
