/**
 * @file   unitary.cc
 * @date   12/2018
 * @author Imran Ashraf
 * @author Anneriet Krol
 * @brief  unitary matrix (decomposition) implementation
 */

#include <unitary.h>

#ifndef WITHOUT_UNITARY_DECOMPOSITION
#include <Eigen/MatrixFunctions>
#include <complex.h>
#define lapack_complex_float    std::complex<float>
#define lapack_complex_double   std::complex<double>
#include <src/misc/lapacke.h>
#endif

typedef unsigned int uint;

#include <chrono>

namespace ql
{

unitary::unitary() : name(""), is_decomposed(false) {}

unitary::unitary(std::string name, std::vector<std::complex<double>> array) :
name(name), array(array), is_decomposed(false) {}

int unitary::size() {
    return array.size();
}

#ifdef WITHOUT_UNITARY_DECOMPOSITION

void unitary::decompose() {
    throw std::runtime_error("unitary decomposition was explicitly disabled in this build!");
}

bool unitary::is_decompose_support_enabled() {
    return false;
}

#else

// JvS: this was originally the class "unitary" itself, but compile times of
// Eigen are so excessive that I moved it into its own compile unit and
// provided a wrapper instead. It doesn't actually NEED to be wrapped like
// this, because the Eigen::Matrix<...> member is actually only used within
// the scope of a single method (calling other methods), but I'm not touching
// this code.


extern "C"
{
	extern void cuncsd_(char *,		// JOBU1
						char *,		// JOBU2
						char *,		// JOBV1T
						char *,		// JOBV2T
						char *,		// TRANS
						char *,		// SIGNS
						int *,		// M
						int *,		// p
						int *,		// Q
						std::complex<float>  *, // X11
						int *,		// LDX11
						std::complex<float>  *, // X12
						int *,		// LDX12
						std::complex<float>  *, // X21
						int *,		// LDX21
						std::complex<float>  *, // X22
						int *,		// LDX22
						float *,	// THETA
						std::complex<float>  *, // U1
						int *,		// LDU1
						std::complex<float>  *, // U2
						int *,		// LDU2
						std::complex<float>  *, // V1T
						int *,		// LDV1T
						std::complex<float>  *, // V2T
						int *,		// LDV2T
						std::complex<float>  *, // WORK
						int *,		// LWORK
						float *,	// RWORK
						int *,		// LRWORK
						int *,		// IWORK
						int *);		// info
}



class UnitaryDecomposer
{
private:
    std::chrono::duration<float> CSD_time;
    std::chrono::duration<float> CSD_time2;
    std::chrono::duration<float> CSD_time3;
    std::chrono::duration<float> zyz_time;
    std::chrono::duration<float> multiplexing_time;
    std::chrono::duration<float> demultiplexing_time;

public:
    std::string name;
    std::vector<std::complex<double>> array;
    std::vector<std::complex<float>> SU;
    bool is_decomposed;
    std::vector<double> instructionlist;


    UnitaryDecomposer() : name(""), is_decomposed(false) {}

    UnitaryDecomposer(std::string name, std::vector<std::complex<double>> array) : 
            name(name), array(array), is_decomposed(false)
    {
        DOUT("constructing unitary: " << name 
                  << ", containing: " << array.size() << " elements");
    }

    int size()
    {
		return (int) array.size();
    }

    Eigen::MatrixXcf getMatrix()
    {
    Eigen::MatrixXcf _matrix;
        if (!array.empty())
        {
            int matrix_size = (int)std::pow(array.size(), 0.5);

            Eigen::Map<Eigen::MatrixXcd> matrix(array.data(), matrix_size, matrix_size);
            _matrix = matrix.cast<std::complex<float> >().transpose(); // Eigen stores matrices Column-Major
        }
        return _matrix;
    }

    void decompose()
    {
        DOUT("decomposing Unitary: " << name);

        Eigen::MatrixXcf _matrix = getMatrix();
        int matrix_size = _matrix.rows();
        
        // compute the number of qubits: length of array is collumns*rows, so log2(sqrt(array.size))
        int numberofbits = uint64_log2(matrix_size);

        // very little accuracy because of tests using printed-from-matlab code that does not have many digits after the comma    
        if( !_matrix.isUnitary(0.001))
        {
            //Throw an error
            EOUT("Unitary " << name <<" is not a unitary matrix!");

            throw ql::exception("Error: Unitary '"+ name+"' is not a unitary matrix. Cannot be decomposed!", false);
        }
        // initialize the general M^k lookuptable
        genMk(numberofbits);

        decomp_function(_matrix, numberofbits); //needed because the matrix is read in columnmajor
         
		DOUT("CSD_time: 	" << CSD_time.count());
		DOUT("CSD_time2:	" << CSD_time2.count());
		DOUT("CSD_time3:	" << CSD_time3.count());
		DOUT("zyz_time: 	" << zyz_time.count());
		DOUT("multiplexing_time:  	" << multiplexing_time.count());
		DOUT("demultiplexing_time:	" << demultiplexing_time.count());
        DOUT("Done decomposing");
        is_decomposed = true;
    }

    std::string to_string(Eigen::MatrixXcf m, std::string vector_prefix = "",
                            std::string elem_sep = ", ")
    {
        std::ostringstream ss;
        ss << "\n" << m ;
        return ss.str();
    }

    std::string to_string(float  *a, int n)
    {
        std::ostringstream ss;
        ss << "\n";
        for (int i = 0; i < n; i++) {ss << a[i] << ", ";};
        return ss.str();
    }
    std::string to_string(std::complex<float>  *a, int n)
    {
        std::ostringstream ss;
        ss << "\n";
        for (int i = 0; i < n; i++) {ss << a[i] << ", ";};
        return ss.str();
    }

    void decomp_function(const Eigen::Ref<const Eigen::MatrixXcf>& matrix, int numberofbits)
    {                 
        if(numberofbits == 1)
        {
            zyz_decomp(matrix);
        }
        else
        {
            int n = matrix.rows()/2;

            // Eigen::MatrixXcf V(n,n);
            // Eigen::MatrixXcf W(n,n);
            // Eigen::VectorXcf D(n);

            // if q2 is zero, the whole thing is a demultiplexing problem instead of full CSD
            if(matrix.bottomLeftCorner(n,n).isZero(10e-14) && matrix.topRightCorner(n,n).isZero(10e-14))
            {
                DOUT("Optimization: q2 of size " << n << " is zero, only demultiplexing will be performed.");
                instructionlist.push_back(200.0);
                if(matrix.topLeftCorner(n, n).isApprox(matrix.bottomRightCorner(n,n),10e-4))
                {
                    DOUT("Optimization: Unitaries are equal, skip one step in the recursion for unitaries of size: " << n);
                    instructionlist.push_back(300.0);
                    decomp_function(matrix.topLeftCorner(n, n), numberofbits-1);
                }
                else
                {
                    demultiplexing(matrix.topLeftCorner(n, n), matrix.bottomRightCorner(n,n), numberofbits-1);
                }
            }
            // Check to see if it the kronecker product of a bigger matrix and the identity matrix.
            // By checking if the first row is equal to the second row one over, and if thelast two rows are equal 
            // Which means the last qubit is not affected by this gate
            else if (matrix(Eigen::seqN(0, n, 2), Eigen::seqN(1, n, 2)).isZero()  && matrix(Eigen::seqN(1, n, 2), Eigen::seqN(0, n, 2)).isZero()  && matrix.block(0,0,1,2*n-1) == matrix.block(1,1,1,2*n-1) &&  matrix.block(2*n-2,0,1,2*n-1) ==  matrix.block(2*n-1,1,1,2*n-1))
            {
                DOUT("Optimization: last qubit is not affected, skip one step in the recursion for matrix size: " << n);
                // Code for last qubit not affected
                instructionlist.push_back(100.0);
                decomp_function(matrix(Eigen::seqN(0, n, 2), Eigen::seqN(0, n, 2)), numberofbits-1);

            }
            else
            {

            auto start = std::chrono::steady_clock::now();
            // CSD(matrix, L0, L1, R0, R1, ss);
            
           	float *THETA = new float[n];
			std::complex<float>  *U1 = new std::complex<float> [n * n];
			std::complex<float>  *U2 = new std::complex<float> [n * n];
			std::complex<float>  *V1T = new std::complex<float> [n * n];
			std::complex<float>  *V2T = new std::complex<float> [n * n];

			std::complex<float> *Uvec = new std::complex<float> [4*n*n];
			Eigen::MatrixXcf::Map(Uvec, 2*n, 2*n) = matrix;
       	
       		uses_cuncsd(Uvec, n, U1, U2, V1T, V2T, THETA);
       		Eigen::VectorXf theta = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(THETA, n);
            Eigen::MatrixXcf L0 = Eigen::Map<Eigen::MatrixXcf, Eigen::Unaligned>(U1, n, n);
            Eigen::MatrixXcf L1 = Eigen::Map<Eigen::MatrixXcf, Eigen::Unaligned>(U2, n, n);
            Eigen::MatrixXcf R0 = Eigen::Map<Eigen::MatrixXcf, Eigen::Unaligned>(V1T, n, n);
            Eigen::MatrixXcf R1 = Eigen::Map<Eigen::MatrixXcf, Eigen::Unaligned>(V2T, n, n);
       		
            Eigen::MatrixXcf tmp(2*n,2*n);
            Eigen::MatrixXcf c = theta.array().cos().matrix().asDiagonal();
            Eigen::MatrixXcf s = (-1*theta.array().sin()).matrix().asDiagonal();
            tmp.topLeftCorner(n,n) = L0*c*R0;
            tmp.bottomLeftCorner(n,n) = -L1*s*R0;
            tmp.topRightCorner(n,n) = L0*s*R1;
            tmp.bottomRightCorner(n,n) = L1*c*R1;

            if (ql::utils::logger::LOG_LEVEL >= ql::utils::logger::log_level_t::LOG_DEBUG)
            {
                if (!matrix.topLeftCorner(n,n).isApprox(L0 * c * R0, 10e-6) || !matrix.bottomLeftCorner(n,n).isApprox(-L1 * s * R0, 10e-6))
                {
                    if (matrix.topLeftCorner(n,n).isApprox(L0 * c * R0, 10e-6))
                    {
                        DOUT("X11 is correct");
                    }
                    else
                    {
                        DOUT("X11 of size" + std::to_string(n) + " is not correct! (is not usually an issue)");
                        DOUT(to_string(matrix.topLeftCorner(n,n)));
                        DOUT(to_string(tmp.topLeftCorner(n,n)));
                    }
                    if (matrix.bottomLeftCorner(n,n).isApprox(-L1 * s * R0, 10e-6))
                    {
                        DOUT("X21 is correct");
                    }
                    else
                    {
                        DOUT("X21 of size" + std::to_string(n) + " is not correct! (is not usually an issue)");
                        DOUT(to_string(matrix.bottomLeftCorner(n,n)));
                        DOUT(to_string(tmp.bottomLeftCorner(n,n)));
                    }
                }

                if (!matrix.topRightCorner(n,n).isApprox(L0 * s * R1, 10e-6) || !matrix.bottomRightCorner(n,n).isApprox(L1 * c * R1, 10e-6))
                {
                    if (matrix.topRightCorner(n,n).isApprox(L0 * s * R1, 10e-6))
                    {
                        DOUT("X12 is correct");
                    }
                    else
                    {
                        DOUT("X12 of size" + std::to_string(n) + " is not correct! (is not usually an issue)");
                        DOUT(to_string(matrix.topRightCorner(n,n)));
                        DOUT(to_string(tmp.topRightCorner(n,n)));
                    }
                    if (matrix.bottomRightCorner(n,n).isApprox(L1 * c * R1, 10e-6))
                    {
                        DOUT("X22 is correct");
                    }
                    else
                    {
                        DOUT("X22 of size" + std::to_string(n) + " is not correct! (is not usually an issue)");
                        DOUT(to_string(matrix.bottomRightCorner(n,n)));
                        DOUT(to_string(tmp.bottomRightCorner(n,n)));
                    }
                }
            }

            // Just to see if it kinda matches

            if(!tmp.isApprox(matrix, 10e-2))
            {
                throw ql::exception("CSD of unitary '"+ name+"' not correct! Failed at matrix size " + std::to_string(n), false);
            }
            CSD_time += (std::chrono::steady_clock::now() - start);
            
            demultiplexing(R0, R1, numberofbits-1);

            multicontrolledY(s, n);

            demultiplexing(L0, L1, numberofbits-1);
            
            }
        }
 
    }

void uses_cuncsd(std::complex<float>  *U, int n, std::complex<float>  *U1, std::complex<float>  *U2, std::complex<float>  *V1T, std::complex<float>  *V2T, float *THETA  )
{
	DOUT("CSD with cuncsd");

	char Ychar = 'Y';
	char trans = 'C'; // = 'T':      X, U1, U2, V1T, and V2T are stored in row-major     order; otherwise:  X, U1, U2, V1T, and V2T are stored in column-major order.
	char SIGNS = '0'; //= 'O':      The lower-left block is made nonpositive (the "other" convention);     otherwise:  The upper-right block is made nonpositive (the  "default" convention).
	int M, P, Q;	  // rows and columns in X, rows in X11 and X12, columns in X11 and X21;

	M = 2*n; // The numbers of rows and columns in X
	P = n;
	Q = n;
	std::complex<float>  *X11 = new std::complex<float> [n * n]; // complex array; matrix upper left corner
	std::complex<float>  *X12 = new std::complex<float> [n * n]; //matrix upper right corner
	std::complex<float>  *X21 = new std::complex<float> [n * n]; //matrix lower left corner
	std::complex<float>  *X22 = new std::complex<float> [n * n]; //matrix lower right corner

	for (int i = 0; i < n; i++) // row
	{
		for (int j = 0; j < n; j++)  // column
		{											  
			X11[n * i + j] = U[2 * n * i + j];		   // value at n*rows+column from array twice as big 
			X21[n * i + j] = U[2 * n * i + j + n];	   // offset by n from block 11
			X12[n * i + j] = U[2 * n * i + j + 2 * n * n]; // two blocks of n*n before it
			X22[n * i + j] = U[2 * n * i + j + 2 * n * n + n]; // two blocks of n*n before it, offset by n
		}
	}

    DOUT("X11" << to_string(X11, n*n));
    DOUT("X12" << to_string(X12, n*n));
    DOUT("X21" << to_string(X21, n*n));
    DOUT("X22" << to_string(X22, n*n));
	int *iwork = new int[n];
	int info;

	int lwork = -1; // The dimension of the array WORK. If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns this value as the first entry of the work array, and no error message related to LWORK is issued by XERBLA.
	int lrwork = -1;	
	std::complex<float>  wkopt;
	float rworkopt;
	cuncsd_(&Ychar, &Ychar, &Ychar, &Ychar, &trans, &SIGNS, &M, &P, &Q, X11, &n, X12, &n, X21, &n, X22, &n, THETA, U1, &n, U2, &n, V1T, &n, V2T, &n, &wkopt, &lwork, &rworkopt, &lrwork, iwork, &info);

	lwork = (int)wkopt.real();
	lrwork = (int) rworkopt;
	std::complex<float>  *work = new std::complex<float> [lwork];
	float *rwork = new float[lrwork];
	cuncsd_(&Ychar, &Ychar, &Ychar, &Ychar, &trans, &SIGNS, &M, &P, &Q, X11, &n, X12, &n, X21, &n, X22, &n, THETA, U1, &n, U2, &n, V1T, &n, V2T, &n, work, &lwork, rwork, &lrwork, iwork, &info);

	// check for errors
	if (info != 0)
	{
		EOUT("Error: cuncsd returned error code " << info);
	}

	// deallocate
	delete[] X11;
	delete[] X12;
	delete[] X21;
	delete[] X22;

	delete[] work;
	delete[] rwork;
	delete[] iwork;
}

    void CSD(const Eigen::Ref<const Eigen::MatrixXcf>& U, Eigen::Ref<Eigen::MatrixXcf> u1, Eigen::Ref<Eigen::MatrixXcf> u2, Eigen::Ref<Eigen::MatrixXcf> v1, Eigen::Ref<Eigen::MatrixXcf> v2, Eigen::Ref<Eigen::MatrixXcf> s)
    {
        auto start = std::chrono::steady_clock::now();        
        //Cosine sine decomposition
        // U = [q1, U01] = [u1    ][c  s][v1  ]
        //     [q2, U11] = [    u2][-s c][   v2]
        int n = U.rows();

        Eigen::BDCSVD<Eigen::MatrixXcf> svd(n/2,n/2);
        svd.compute(U.topLeftCorner(n/2,n/2), Eigen::ComputeThinU | Eigen::ComputeThinV); // possible because it's square anyway
        

        // thinCSD: q1 = u1*c*v1.adjoint()
        //          q2 = u2*s*v1.adjoint()
        int p = n/2;

        Eigen::MatrixXcf c(svd.singularValues().reverse().asDiagonal());
        u1.noalias() = svd.matrixU().rowwise().reverse();
        v1.noalias() = svd.matrixV().rowwise().reverse(); // Same v as in matlab: u*s*v.adjoint() = q1

        Eigen::MatrixXcf q2 = U.bottomLeftCorner(p,p)*v1;      

        int k = 0;
        for(int j = 1; j < p; j++)
        {
            if(c(j,j).real() <= 0.70710678119)
            {
                k = j;
            }
        }


        Eigen::HouseholderQR<Eigen::MatrixXcf> qr(p,k+1);
        qr.compute(q2.block( 0,0, p, k+1));
        u2 = qr.householderQ();
        s.noalias() = u2.adjoint()*q2;
        if(k < p-1)
        {
            k = k+1;
            Eigen::BDCSVD<Eigen::MatrixXcf> svd2(p-k, p-k);
            svd2.compute(s.block(k, k, p-k, p-k), Eigen::ComputeThinU | Eigen::ComputeThinV);
            s.block(k, k, p-k, p-k) = svd2.singularValues().asDiagonal();
            c.block(0,k, p,p-k) = c.block(0,k, p,p-k)*svd2.matrixV();
            u2.block(0,k, p,p-k) = u2.block(0,k, p,p-k)*svd2.matrixU();
            v1.block(0,k, p,p-k) = v1.block(0,k, p,p-k)*svd2.matrixV();
            
            Eigen::HouseholderQR<Eigen::MatrixXcf> qr2(p-k, p-k);

            qr2.compute(c.block(k,k, p-k,p-k));
            c.block(k,k,p-k,p-k) = qr2.matrixQR().triangularView<Eigen::Upper>();
            u1.block(0,k, p,p-k) = u1.block(0,k, p,p-k)*qr2.householderQ(); 
        }
        CSD_time2 += (std::chrono::steady_clock::now() - start);

        auto start2 = std::chrono::steady_clock::now();



        std::vector<int> c_ind;
        std::vector<int> s_ind;
        for(int j = 0; j < p; j++)
        {
            if(c(j,j).real() < 0)
            {
                c_ind.push_back(j);
            }
            if(s(j,j).real() < 0)
            {
                s_ind.push_back(j);
            } 
        }

        c(c_ind,c_ind) = -c(c_ind,c_ind);
        u1(Eigen::all, c_ind) = -u1(Eigen::all, c_ind);

        s(s_ind,s_ind) = -s(s_ind,s_ind);
        u2(Eigen::all, s_ind) = -u2(Eigen::all, s_ind);

        if ( ql::utils::logger::LOG_LEVEL >= ql::utils::logger::log_level_t::LOG_DEBUG ){
            if (!U.topLeftCorner(p, p).isApprox(u1 * c * v1.adjoint(), 10e-8) || !U.bottomLeftCorner(p, p).isApprox(u2 * s * v1.adjoint(), 10e-8))
            {
                if (U.topLeftCorner(p, p).isApprox(u1 * c * v1.adjoint(), 10e-8))
                {
                    DOUT("q1 is correct");
                }
                else
                {
                    DOUT("q1 of size" + std::to_string(p) + " is not correct! (is not usually an issue)");
                }
                if (U.bottomLeftCorner(p, p).isApprox(u2 * s * v1.adjoint(), 10e-8))
                {
                    DOUT("q2 is correct");
                }
                else
                {
                    DOUT("q2 of size" + std::to_string(p) + " is not correct! (is not usually an issue)");
                }
            }
        }
        v1.adjointInPlace(); // Use this instead of = v1.adjoint (to avoid aliasing issues)
        s = -s;

        Eigen::MatrixXcf tmp_s = u1.adjoint()*U.topRightCorner(p,p);
        Eigen::MatrixXcf tmp_c = u2.adjoint()*U.bottomRightCorner(p,p);

        for(int i = 0; i < p; i++)
        {
            if(std::abs(s(i,i)) > std::abs(c(i,i)))
            {
                v2.row(i).noalias() = tmp_s.row(i)/s(i,i);                
            }
            else
            {
                v2.row(i).noalias() = tmp_c.row(i)/c(i,i);
            }
        }
        // U = [q1, U01] = [u1    ][c  s][v1  ]
        //     [q2, U11] = [    u2][-s c][   v2]
    
        Eigen::MatrixXcf tmp(n,n);
        tmp.topLeftCorner(p,p) = u1*c*v1;
        tmp.bottomLeftCorner(p,p) = -u2*s*v1;
        tmp.topRightCorner(p,p) = u1*s*v2;
        tmp.bottomRightCorner(p,p) = u2*c*v2;
        
        // Just to see if it kinda matches
        if(!tmp.isApprox(U, 10e-2))
        {
            throw ql::exception("CSD of unitary '"+ name+"' not correct! Failed at matrix size " + std::to_string(n), false);
        }
            CSD_time3 += (std::chrono::steady_clock::now() - start2);

    }


    void zyz_decomp(const Eigen::Ref<const Eigen::MatrixXcf>& matrix)
    {
        auto start = std::chrono::steady_clock::now();

        ql::complex_t det = matrix.determinant();

        float delta = atan2(det.imag(), det.real())/matrix.rows();
        std::complex<float> A = exp(std::complex<float>(0,-1)*delta)*matrix(0,0);
        std::complex<float> B = exp(std::complex<float>(0,-1)*delta)*matrix(0,1); // to comply with the other y-gate definition
        
        float sw = sqrt(pow((float) B.imag(),2) + pow((float) B.real(),2) + pow((float) A.imag(),2));
        float wx = 0;
        float wy = 0;
        float wz = 0;
        if(sw > 0)
        {
        wx = B.imag()/sw;
        wy = B.real()/sw;
        wz = A.imag()/sw;
        }
        float t1 = atan2(A.imag(),A.real());
        float t2 = atan2(B.imag(), B.real());
        double alpha = t1+t2;
        double gamma = t1-t2;
        double beta = 2*atan2(sw*sqrt(pow(wx,2)+pow(wy,2)),sqrt(pow( A.real(),2)+pow((wz*sw),2)));
        instructionlist.push_back(-gamma);
        instructionlist.push_back(-beta);
        instructionlist.push_back(-alpha);
        zyz_time += (std::chrono::steady_clock::now() - start);
    }

    void demultiplexing(const Eigen::Ref<const Eigen::MatrixXcf> &U1, const Eigen::Ref<const Eigen::MatrixXcf> &U2, int numberofcontrolbits)
    {
        // [U1 0 ]  = [V 0][D 0 ][W 0]
        // [0  U2]    [0 V][0 D*][0 W]
        
        auto start = std::chrono::steady_clock::now(); 
        Eigen::MatrixXcf check = U1*U2.adjoint();
       	
       	Eigen::MatrixXcf V;
       	Eigen::MatrixXcf W;
		Eigen::VectorXcf D;
		
        if(check == check.adjoint())
        {
            DOUT("Demultiplexing matrix of size " + std::to_string(check.rows())+" is self-adjoint()");
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcf> eigslv(check);
            D = ((Eigen::MatrixXcf) eigslv.eigenvalues()).cwiseSqrt();
            V.noalias() = eigslv.eigenvectors();
            W.noalias() = D.asDiagonal()*V.adjoint()*U2;
        }
        else
        {
            if (numberofcontrolbits < 5) //  Schur is faster for small matrices
            {
                Eigen::ComplexSchur<Eigen::MatrixXcf> decomposition(check);
                D.noalias() = decomposition.matrixT().diagonal().cwiseSqrt();
                V.noalias() = decomposition.matrixU();
                W.noalias() = D.asDiagonal() * V.adjoint() * U2;
            }
            else
            {
                Eigen::ComplexEigenSolver<Eigen::MatrixXcf> decomposition(check);
                D = decomposition.eigenvalues().cwiseSqrt();
                V.noalias() = decomposition.eigenvectors();
                W.noalias() = D.asDiagonal() * V.adjoint() * U2;
            }
        }
    

        if(!V.isUnitary(10e-3))
        {
            DOUT("Eigenvalue decomposition incorrect: V of size " + std::to_string(V.rows())+"is not unitary, adjustments will be made");
                Eigen::ComplexSchur<Eigen::MatrixXcf> decomposition(check);
                D.noalias() = decomposition.matrixT().diagonal().cwiseSqrt();
                V.noalias() = decomposition.matrixU();
                W.noalias() = D.asDiagonal() * V.adjoint() * U2;             
        }

        demultiplexing_time += (std::chrono::steady_clock::now() - start);
        
        Eigen::MatrixXcf Dtemp = D.asDiagonal();
        if(!U1.isApprox(V*Dtemp*W, 10e-2) || !U2.isApprox(V*Dtemp.adjoint()*W, 10e-2))
        {
            EOUT("Demultiplexing not correct!");
            throw ql::exception("Demultiplexing of unitary '"+ name+"' not correct! Failed at matrix size "+ std::to_string(V.rows()), false);
        }


		decomp_function(W, numberofcontrolbits);
		multicontrolledZ(D, D.rows());
		decomp_function(V, numberofcontrolbits);
    }


    std::vector<Eigen::MatrixXf> genMk_lookuptable;

    // returns M^k = (-1)^(b_(i-1)*g_(i-1)), where * is bitwise inner product, g = binary gray code, b = binary code.
    void genMk(int numberofqubits)
    {
        // int numberqubits = uint64_log2(_matrix.rows());
        for(int n = 1; n <= numberofqubits; n++)
        {
            int size=1<<n;
            Eigen::MatrixXf Mk(size,size);
            for(int i = 0; i < size; i++)
            {
                for(int j = 0; j < size ;j++)
                {
                    Mk(i,j) =std::pow(-1, bitParity(i&(j^(j>>1))));
                }
            }
        genMk_lookuptable.push_back(Mk);
        }
    }

    // source: https://stackoverflow.com/questions/994593/how-to-do-an-integer-log2-in-c user Todd Lehman
    int uint64_log2(uint64_t n)
    {
    #define S(k) if (n >= (UINT64_C(1) << k)) { i += k; n >>= k; }

    int i = -(n == 0); S(32); S(16); S(8); S(4); S(2); S(1); return i;

    #undef S
    }

    int bitParity(int i)
    {
        if (i < 2 << 16)
        {
            i = (i >> 16) ^ i;
            i = (i >> 8) ^ i;
            i = (i >> 4) ^ i;
            i = (i >> 2) ^ i;
            i = (i >> 1) ^ i;
            return i % 2;
        }
        else
        {
            throw ql::exception("Bit parity number too big!", false);
        }
    }

    void multicontrolledY(const Eigen::Ref<const Eigen::MatrixXcf> &ss, int halfthesizeofthematrix)
    {
        auto start = std::chrono::steady_clock::now();
        Eigen::VectorXf temp =  2*Eigen::asin(ss.diagonal().array()).real();
        Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXf> dec(genMk_lookuptable[uint64_log2(halfthesizeofthematrix)-1]);
        Eigen::VectorXf tr = dec.solve(temp);

        // #Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXf> dec(genMk_lookuptable[uint64_log2(halfthesizeofthematrix)-1]);
        // Eigen::VectorXf tr = dec.solve(theta);
        
        // Check is very approximate to account for low-precision input matrices
        if(!temp.isApprox(genMk_lookuptable[uint64_log2(halfthesizeofthematrix)-1]*tr, 10e-2))
        {
                EOUT("Multicontrolled Y not correct!");
                throw ql::exception("Demultiplexing of unitary '"+ name+"' not correct! Failed at matrix size " +std::to_string(halfthesizeofthematrix), false);
        }

        instructionlist.insert(instructionlist.end(), &tr[0], &tr[halfthesizeofthematrix]);
        multiplexing_time += std::chrono::steady_clock::now() - start;
    }

    void multicontrolledZ(const Eigen::Ref<const Eigen::VectorXcf> &D, int halfthesizeofthematrix)
    {
        auto start = std::chrono::steady_clock::now();
        
        Eigen::VectorXf temp =  (std::complex<float>(0,-2)*Eigen::log(D.array())).real();
        Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXf> dec(genMk_lookuptable[uint64_log2(halfthesizeofthematrix)-1]);
        Eigen::VectorXf tr = dec.solve(temp);
        
        // Check is very approximate to account for low-precision input matrices
        if(!temp.isApprox(genMk_lookuptable[uint64_log2(halfthesizeofthematrix)-1]*tr, 10e-2))
        {
                EOUT("Multicontrolled Z not correct!");
                throw ql::exception("Demultiplexing of unitary '"+ name+"' not correct! Failed at matrix size "+std::to_string(halfthesizeofthematrix), false);
        }
        

        instructionlist.insert(instructionlist.end(), &tr[0], &tr[halfthesizeofthematrix]);
        multiplexing_time += std::chrono::steady_clock::now() - start;

    }
    ~UnitaryDecomposer()
    {
        // destroy unitary
        DOUT("destructing unitary: " << name);
    }
};

void unitary::decompose() {
    UnitaryDecomposer decomposer(name, array);
    decomposer.decompose();
    SU = decomposer.SU;
    is_decomposed = decomposer.is_decomposed;
    instructionlist = decomposer.instructionlist;
}

bool unitary::is_decompose_support_enabled() {
    return true;
}

#endif

}

