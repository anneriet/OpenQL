/** \file
 * Unitary matrix (decomposition) implementation.
 */

#include "unitary.h"

#include "utils/exception.h"

#ifndef WITHOUT_UNITARY_DECOMPOSITION
#include <Eigen/MatrixFunctions>
#include <complex.h>
#define lapack_complex_float    std::complex<float>
#define lapack_complex_double   std::complex<double>
#include <src/misc/lapacke.h>
#endif

namespace ql {

using namespace utils;

unitary::unitary() :
    name(""),
    is_decomposed(false)
{
}

unitary::unitary(
    const Str &name,
    const Vec<Complex> &array
) :
    name(name),
    array(array),
    is_decomposed(false)
{
}

utils::int unitary::size() const {
    // JvS: Note that the original unitary::size() used
    // Eigen::Matrix::size() if the array is empty. However, if the array
    // is empty, _matrix is never initialized beyond its default ctor,
    // which "allocates" a 0x0 matrix, and is thus size 0, exactly what
    // array.size() would return.
    // Don't get me started about why this returns a Real.
    return array.size();
}

#ifdef WITHOUT_UNITARY_DECOMPOSITION

void unitary::decompose() {
    throw Exception("unitary decomposition was explicitly disabled in this build!");
}

Bool unitary::is_decompose_support_enabled() {
    return false;
}

#else

// JvS: this was originally the class "unitary" itself, but compile times of
// Eigen are so excessive that I moved it into its own compile unit and
// provided a wrapper instead. It doesn't actually NEED to be wrapped like
// this, because the Eigen::Matrix<...> member is actually only used within
// the scope of a single method (calling other methods), but I'm not touching
// this code.
class UnitaryDecomposer {
private:
//    typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> complex_matrix ;
//  complex_matrix _matrix;

public:
    Str name;
    Vec<Complex> array;
    Vec<Complex> SU;
    Bool is_decomposed;
    Vec<Real> instructionlist;

    UnitaryDecomposer() : name(""), is_decomposed(false) {}

    UnitaryDecomposer(
        const Str &name,
        const Vec<Complex> &array
    ) :
        name(name),
        array(array),
        is_decomposed(false)
    {
        QL_DOUT("constructing unitary: " << name
                  << ", containing: " << array.size() << " elements");
    }

    Int size()
    {
		return array.size();
    }

    Eigen::MatrixXcd getMatrix()
    {
    Eigen::MatrixXcd _matrix;
        if (!array.empty())
        {
            Int matrix_size = (int)std::pow(array.size(), 0.5);

            Eigen::Map<Eigen::MatrixXcd> matrix(array.data(), matrix_size, matrix_size);
            _matrix = matrix.transpose();
        }
        return _matrix;
    }

    void decompose() {
        QL_DOUT("decomposing Unitary: " << name);

        getMatrix();
        Int matrix_size = _matrix.rows();
        
        // compute the number of qubits: length of array is collumns*rows, so log2(sqrt(array.size))
        Int numberofbits = uint64_log2(matrix_size);

        Eigen::MatrixXcd identity = Eigen::MatrixXcd::Identity(matrix_size, matrix_size);
        Eigen::MatrixXcd matmatadjoint = (_matrix.adjoint()*_matrix);
        // very little accuracy because of tests using printed-from-matlab code that does not have many digits after the comma
        if (!matmatadjoint.isApprox(identity, 0.001)) {
            //Throw an error
            QL_EOUT("Unitary " << name <<" is not a unitary matrix!");

            throw utils::Exception("Error: Unitary '"+ name+"' is not a unitary matrix. Cannot be decomposed!" + to_string(matmatadjoint), false);
        }
        // initialize the general M^k lookuptable
        genMk(numberofbits);

        decomp_function(_matrix, numberofbits); //needed because the matrix is read in columnmajor

        QL_DOUT("Done decomposing");
        is_decomposed = true;
    }

    Str to_string(
        const complex_matrix &m,
        const Str &vector_prefix = "",
        const Str &elem_sep = ", "
    ) {
        StrStrm ss;
        ss << m << "\n";
        return ss.str();
    }

    void decomp_function(const Eigen::Ref<const complex_matrix>& matrix, Int numberofbits) {
        QL_DOUT("decomp_function: \n" << to_string(matrix));
        if(numberofbits == 1) {
            zyz_decomp(matrix);
        } else {
            Int n = matrix.rows()/2;

            Eigen::MatrixXcd V(n,n);
            Eigen::MatrixXcd W(n,n);
            Eigen::VectorXcd D(n);
            // if q2 is zero, the whole thing is a demultiplexing problem instead of full CSD
            if (matrix.bottomLeftCorner(n,n).isZero(10e-14) && matrix.topRightCorner(n,n).isZero(10e-14)) {
                QL_DOUT("Optimization: q2 is zero, only demultiplexing will be performed.");
                instructionlist.push_back(200.0);
                if (matrix.topLeftCorner(n, n).isApprox(matrix.bottomRightCorner(n,n),10e-4)) {
                    QL_DOUT("Optimization: Unitaries are equal, skip one step in the recursion for unitaries of size: " << n << " They are both: " << matrix.topLeftCorner(n, n));
                    instructionlist.push_back(300.0);
                    decomp_function(matrix.topLeftCorner(n, n), numberofbits-1);
                } else {
                    demultiplexing(matrix.topLeftCorner(n, n), matrix.bottomRightCorner(n,n), V, D, W, numberofbits-1);

                    decomp_function(W, numberofbits-1);
                    multicontrolledZ(D, D.rows());
                    decomp_function(V, numberofbits-1);
                }
            } else if (
                // Check to see if it the kronecker product of a bigger matrix and the identity matrix.
                // By checking if the first row is equal to the second row one over, and if thelast two rows are equal
                // Which means the last qubit is not affected by this gate
                matrix(Eigen::seqN(0, n, 2), Eigen::seqN(1, n, 2)).isZero()
                && matrix(Eigen::seqN(1, n, 2), Eigen::seqN(0, n, 2)).isZero()
                && matrix.block(0,0,1,2*n-1) == matrix.block(1,1,1,2*n-1)
                && matrix.block(2*n-2,0,1,2*n-1) == matrix.block(2*n-1,1,1,2*n-1)
            ) {
                QL_DOUT("Optimization: last qubit is not affected, skip one step in the recursion.");
                // Code for last qubit not affected
                instructionlist.push_back(100.0);
                decomp_function(matrix(Eigen::seqN(0, n, 2), Eigen::seqN(0, n, 2)), numberofbits-1);
            } else {
                complex_matrix ss(n,n);
                complex_matrix L0(n,n);
                complex_matrix L1(n,n);
                complex_matrix R0(n,n);
                complex_matrix R1(n,n);
                CSD(matrix, L0, L1, R0, R1, ss);
                demultiplexing(R0, R1, V, D, W, numberofbits-1);
                decomp_function(W, numberofbits-1);
                multicontrolledZ(D, D.rows());
                decomp_function(V, numberofbits-1);

                multicontrolledY(ss.diagonal(), n);

                demultiplexing(L0, L1, V, D, W, numberofbits-1);
                decomp_function(W, numberofbits-1);
                multicontrolledZ(D, D.rows());
                decomp_function(V, numberofbits-1);
            }
        }
 
    }

    void CSD(
        const Eigen::Ref<const complex_matrix> &U,
        Eigen::Ref<complex_matrix> u1,
        Eigen::Ref<complex_matrix> u2,
        Eigen::Ref<complex_matrix> v1,
        Eigen::Ref<complex_matrix> v2,
        Eigen::Ref<complex_matrix> s
    ) {
        //Cosine sine decomposition
        // U = [q1, U01] = [u1    ][c  s][v1  ]
        //     [q2, U11] = [    u2][-s c][   v2]
        Int n = U.rows();
        // complex_matrix c(n,n); // c matrix is not needed for the higher level
        // complex_matrix q1 = U.topLeftCorner(n/2,m/2);

        Eigen::BDCSVD<Eigen::MatrixXcd> svd(n/2,n/2);
        svd.compute(U.topLeftCorner(n/2,n/2), Eigen::ComputeThinU | Eigen::ComputeThinV); // possible because it's square anyway


        // thinCSD: q1 = u1*c*v1.adjoint()
        //          q2 = u2*s*v1.adjoint()
        Int p = n/2;
        // complex_matrix z = Eigen::MatrixXd::Identity(p, p).colwise().reverse();
        complex_matrix c(svd.singularValues().reverse().asDiagonal());
        u1.noalias() = svd.matrixU().rowwise().reverse();
        v1.noalias() = svd.matrixV().rowwise().reverse(); // Same v as in matlab: u*s*v.adjoint() = q1

        complex_matrix q2 = U.bottomLeftCorner(p,p)*v1;

        Int k = 0;
        for (Int j = 1; j < p; j++) {
            if (c(j,j).real() <= 0.70710678119) {
                k = j;
            }
        }


        Eigen::HouseholderQR<Eigen::MatrixXcd> qr(p,k+1);
        qr.compute(q2.block( 0,0, p, k+1));
        u2 = qr.householderQ();
        s.noalias() = u2.adjoint()*q2;
        if (k < p-1) {
            QL_DOUT("k is smaller than size of q1 = "<< p << ", adjustments will be made, k = " << k);
            k = k+1;
            Eigen::BDCSVD<Eigen::MatrixXcd> svd2(p-k, p-k);
            svd2.compute(s.block(k, k, p-k, p-k), Eigen::ComputeThinU | Eigen::ComputeThinV);
            s.block(k, k, p-k, p-k) = svd2.singularValues().asDiagonal();
            c.block(0,k, p,p-k) = c.block(0,k, p,p-k)*svd2.matrixV();
            u2.block(0,k, p,p-k) = u2.block(0,k, p,p-k)*svd2.matrixU();
            v1.block(0,k, p,p-k) = v1.block(0,k, p,p-k)*svd2.matrixV();

            Eigen::HouseholderQR<complex_matrix> qr2(p-k, p-k);

            qr2.compute(c.block(k,k, p-k,p-k));
            c.block(k,k,p-k,p-k) = qr2.matrixQR().triangularView<Eigen::Upper>();
            u1.block(0,k, p,p-k) = u1.block(0,k, p,p-k)*qr2.householderQ();
        }



        Vec<Int> c_ind;
        Vec<Int> s_ind;
        for (Int j = 0; j < p; j++) {
            if (c(j,j).real() < 0) {
                c_ind.push_back(j);
            }
            if (s(j,j).real() < 0) {
                s_ind.push_back(j);
            }
        }

        c(c_ind,c_ind) = -c(c_ind,c_ind);
        u1(Eigen::all, c_ind) = -u1(Eigen::all, c_ind);

        s(s_ind,s_ind) = -s(s_ind,s_ind);
        u2(Eigen::all, s_ind) = -u2(Eigen::all, s_ind);

        v1.adjointInPlace(); // Use this instead of = v1.adjoint (to avoid aliasing issues)
        s = -s;

        Eigen::MatrixXcd tmp_s = u1.adjoint()*U.topRightCorner(p,p);
        Eigen::MatrixXcd tmp_c = u2.adjoint()*U.bottomRightCorner(p,p);

        // Vec<Int> c_ind_row;
        // Vec<Int> s_ind_row;
        for (Int i = 0; i < p; i++) {
            if (abs(s(i,i)) > abs(c(i,i))) {
                // Vec<Int> s_ind_row;
                v2.row(i).noalias() = tmp_s.row(i)/s(i,i);
            } else {
                // c_ind_row.push_back(i);
                v2.row(i).noalias() = tmp_c.row(i)/c(i,i);
            }
        }

        // U = [q1, U01] = [u1    ][c  s][v1  ]
        //     [q2, U11] = [    u2][-s c][   v2]

        complex_matrix tmp(n,n);
        tmp.topLeftCorner(p,p) = u1*c*v1;
        tmp.bottomLeftCorner(p,p) = -u2*s*v1;
        tmp.topRightCorner(p,p) = u1*s*v2;
        tmp.bottomRightCorner(p,p) = u2*c*v2;
        
        // Just to see if it kinda matches
        if (!tmp.isApprox(U, 10e-2)) {
            throw utils::Exception("CSD of unitary '"+ name+"' is wrong! Failed at matrix: \n"+to_string(tmp) + "\nwhich should be: \n" + to_string(U), false);
        }

    }

    void zyz_decomp(const Eigen::Ref<const complex_matrix> &matrix) {

        Complex det = matrix.determinant();// matrix(0,0)*matrix(1,1)-matrix(1,0)*matrix(0,1);

        Real delta = atan2(det.imag(), det.real())/matrix.rows();
        Complex A = exp(Complex(0,-1)*delta)*matrix(0,0);
        Complex B = exp(Complex(0,-1)*delta)*matrix(0,1); //to comply with the other y-gate definition

        Real sw = sqrt(pow((Real) B.imag(),2) + pow((Real) B.real(),2) + pow((Real) A.imag(),2));
        Real wx = 0;
        Real wy = 0;
        Real wz = 0;

        if (sw > 0) {
            wx = B.imag()/sw;
            wy = B.real()/sw;
            wz = A.imag()/sw;
        }

        Real t1 = atan2(A.imag(),A.real());
        Real t2 = atan2(B.imag(), B.real());
        alpha = t1+t2;
        gamma = t1-t2;
        beta = 2*atan2(sw*sqrt(pow((Real) wx,2)+pow((Real) wy,2)),sqrt(pow((Real) A.real(),2)+pow((wz*sw),2)));
        instructionlist.push_back(-gamma);
        instructionlist.push_back(-beta);
        instructionlist.push_back(-alpha);
    }

    void demultiplexing(
        const Eigen::Ref<const complex_matrix> &U1,
        const Eigen::Ref<const complex_matrix> &U2,
        Eigen::Ref<complex_matrix> V,
        Eigen::Ref<Eigen::VectorXcd> D,
        Eigen::Ref<complex_matrix> W,
        Int numberofcontrolbits
    ) {
        // [U1 0 ]  = [V 0][D 0 ][W 0]
        // [0  U2]    [0 V][0 D*][0 W]
        complex_matrix check = U1*U2.adjoint();
        if (check == check.adjoint()) {
            QL_IOUT("Demultiplexing matrix is self-adjoint()");
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigslv(check);
            D = ((Eigen::MatrixXcd) eigslv.eigenvalues()).cwiseSqrt();
            V.noalias() = eigslv.eigenvectors();
            W.noalias() = D.asDiagonal()*V.adjoint()*U2;
        } else {
            if (numberofcontrolbits < 5) {//schur is faster for small matrices
                Eigen::ComplexSchur<complex_matrix> decomposition(check);
                D.noalias() = decomposition.matrixT().diagonal().cwiseSqrt();
                V.noalias() = decomposition.matrixU();
                W.noalias() = D.asDiagonal() * V.adjoint() * U2;
            } else {
                Eigen::ComplexEigenSolver<complex_matrix> decomposition(check);
                D.noalias() = decomposition.eigenvalues().cwiseSqrt();
                V.noalias() = decomposition.eigenvectors();
                W.noalias() = D.asDiagonal() * V.adjoint() * U2;
            }
        }
    

        if (!(V*V.adjoint()).isApprox(Eigen::MatrixXd::Identity(V.rows(), V.rows()), 10e-3)) {
        {
            QL_DOUT("Eigenvalue decomposition incorrect: V is not unitary, adjustments will be made");
            Eigen::ComplexSchur<Eigen::MatrixXcd> decomposition(check);
            D.noalias() = decomposition.matrixT().diagonal().cwiseSqrt();
            V.noalias() = decomposition.matrixU();
            W.noalias() = D.asDiagonal() * V.adjoint() * U2;             
        }
        
        Eigen::MatrixXcd Dtemp = D.asDiagonal();
        if(!U1.isApprox(V*Dtemp*W, 10e-2) || !U2.isApprox(V*Dtemp.adjoint()*W, 10e-2))
        {
            QL_EOUT("Demultiplexing not correct!");
            throw utils::Exception("Demultiplexing of unitary '"+ name+"' not correct! Failed at matrix U1: \n"+to_string(U1)+ "and matrix U2: \n" +to_string(U2) + "\nwhile they are: \n" + to_string(V*D.asDiagonal()*W) + "\nand \n" + to_string(V*D.conjugate().asDiagonal()*W), false);
        }


		decomp_function(W, numberofcontrolbits);
		multicontrolledZ(D, D.rows());
		decomp_function(V, numberofcontrolbits);
    }

    Vec<Eigen::MatrixXd> genMk_lookuptable;

    // returns M^k = (-1)^(b_(i-1)*g_(i-1)), where * is bitwise inner product, g = binary gray code, b = binary code.
    void genMk(int numberofqubits)
    {
        // int numberqubits = uint64_log2(_matrix.rows());
        for(Int n = 1; n <= numberofqubits; n++)
        {
            Int size=1<<n;
            Eigen::MatrixXd Mk(size,size);
            for (Int i = 0; i < size; i++) {
                for (Int j = 0; j < size ;j++) {
                    Mk(i,j) = pow(-1, bitParity(i&(j^(j>>1))));
                }
            }
            genMk_lookuptable.push_back(Mk);
        }
    }

    // source: https://stackoverflow.com/questions/994593/how-to-do-an-integer-log2-in-c user Todd Lehman
    Int uint64_log2(uint64_t n) {
#define S(k) if (n >= (UINT64_C(1) << k)) { i += k; n >>= k; }
        Int i = -(n == 0); S(32); S(16); S(8); S(4); S(2); S(1); return i;
#undef S
    }

    Int bitParity(Int i) {
        if (i < 2 << 16) {
            i = (i >> 16) ^ i;
            i = (i >> 8) ^ i;
            i = (i >> 4) ^ i;
            i = (i >> 2) ^ i;
            i = (i >> 1) ^ i;
            return i % 2;
        } else {
            throw utils::Exception("Bit parity number too big!", false);
        }
    }

    void multicontrolledY(const Eigen::Ref<const Eigen::VectorXcd> &ss, Int halfthesizeofthematrix) {
        Eigen::VectorXd temp =  2*Eigen::asin(ss.array()).real();
        Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> dec(genMk_lookuptable[uint64_log2(halfthesizeofthematrix)-1]);
        Eigen::VectorXd tr = dec.solve(temp);
        
        // Check is very approximate to account for low-precision input matrices
        if (!temp.isApprox(genMk_lookuptable[uint64_log2(halfthesizeofthematrix)-1]*tr, 10e-2)) {
            QL_EOUT("Multicontrolled Y not correct!");
            throw utils::Exception("Demultiplexing of unitary '"+ name+"' not correct! Failed at demultiplexing of matrix ss: \n"  + to_string(ss), false);
        }

        instructionlist.insert(instructionlist.end(), &tr[0], &tr[halfthesizeofthematrix]);
    }

    void multicontrolledZ(const Eigen::Ref<const Eigen::VectorXcd> &D, Int halfthesizeofthematrix) {

        Eigen::VectorXd temp =  (Complex(0,-2)*Eigen::log(D.array())).real();
        Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> dec(genMk_lookuptable[uint64_log2(halfthesizeofthematrix)-1]);
        Eigen::VectorXd tr = dec.solve(temp);
        
        // Check is very approximate to account for low-precision input matrices
        if (!temp.isApprox(genMk_lookuptable[uint64_log2(halfthesizeofthematrix)-1]*tr, 10e-2)) {
            QL_EOUT("Multicontrolled Z not correct!");
            throw utils::Exception("Demultiplexing of unitary '"+ name+"' not correct! Failed at demultiplexing of matrix D: \n"+ to_string(D), false);
        }

        instructionlist.insert(instructionlist.end(), &tr[0], &tr[halfthesizeofthematrix]);

    }

    ~UnitaryDecomposer() {
        // destroy unitary
        QL_DOUT("destructing unitary: " << name);
    }
};

void unitary::decompose() {
    UnitaryDecomposer decomposer(name, array);
    decomposer.decompose();
    SU = decomposer.SU;
    is_decomposed = decomposer.is_decomposed;
    instructionlist = decomposer.instructionlist;
}

Bool unitary::is_decompose_support_enabled() {
    return true;
}

#endif

} // namespace ql
