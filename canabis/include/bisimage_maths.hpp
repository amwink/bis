#ifndef BISIMAGE_MATHS_HPP_INCLUDED
#define BISIMAGE_MATHS_HPP_INCLUDED

/* Much of the basic 2/3/4D matrix/vector maths 
 * in this header is inspired by Linas Vepstas'
 * https://fossies.org/linux/gle/src/vvector.h
 *
 * These classes are a middle ground between the 
 * speed of #define directives and the flexibility
 * of class definitions.
 * 
 * 
 *
 * To request more functionality, contact
 * Alle Meije Wink, mail: a.m.wink@gmail.com
 *
 */

#define _USE_MATH_DEFINES
#include<math.h>

#include<vector>
#include<cfloat>
#include<ranges>
#include<climits>
#include<cassert>
#include<numeric>
#include<iostream>
#include<algorithm>



namespace bis {

    

/* templated signum function
 *
 * 
 */
template < typename T >
T signum (T x) {
  if (x > 0) return  1;
  if (x < 0) return -1;
  return x;
}



/* Vectors and matrices in 2 dimensions
 *
 * For manipulating 2-dimensional (2D) coordinates, matrices and 2D
 * vectors, and their operations, are important tools. Examples are
 * operations on digital pictures or slices of a 3D scan. To access
 * and change such data, it is important to optimise the operations
 * that handle 2D images.
 */

/* Vectors and matrices in 3 dimensions
 *
 * For manipulating 3-dimensional (3D) coordinates, 3D matrices and
 * vectors, and their operations, are important tools. Examples are
 * voxel (3D pixel) coordinates in an MRI scan. To access and change
 * data at locations in these scans, it is important to optimise the
 * operations that handle 3D image space.
 */

/* Vectors and matrices in 4 dimensions
 *
 * An important application of 4-dimensional (4D) matrices and vectors
 * is the manipulation of homogeneous co-ordinates of 3D coordinates.
 * An affine transformation of a 3D vector cannot be represented in a
 * single 3D matrix. By extending the vector [x y z] to 4D [x y z 1].
 * Translations and perspective can then be expressed as matrix
 * multiplications.
 */



// forward declarations
template <class, unsigned>
class vecN;
template <class, unsigned>
class matN;

// eigenvalue and -vector
template <typename T, unsigned S>
struct ev {
	vecN<T,S>	v;
	matN<T,S>	m;
};

// store rotations
typedef struct { 
	double c; 
	double s; 
	double t; 
} rotation;

template <typename T, unsigned S>
class vecN {

    template <typename T2, unsigned S2>
    friend std::ostream& operator<<(std::ostream& out, const vecN<T2,S2>& v);

    template <typename T2, unsigned S2>
    friend std::istream& operator>>(std::istream& in, vecN<T2,S2>& v);

    protected:
        std::vector<T> 
			data;
		unsigned
			sz = S;
		

    public:
        vecN (                    ): data ( S )
            {};																			// default constructor

        vecN ( T x, T y ): data ( S )
            { data[0] = x; if (S>1) data[1] = y; }										// constructor from 2 scalars

        vecN ( T x, T y, T z ): data ( S )
            { data[0] = x; if (S>1) data[1] = y; if (S>2) data[2] = z; }				// constructor from 3 scalars

        vecN ( T x, T y, T z, T t ): data ( S )
            { data[0] = x; if (S>1) data[1] = y; 
						   if (S>2) data[2] = z; 
						   if (S>3) data[3] = t; }										// constructor from 4 scalars

		vecN ( T* xyz             ): data ( S )
            { std::copy_n (xyz, S, data.begin() ); }									// constructor from pointer

		vecN ( std::initializer_list<T> l ): data ( S )
			{ auto R = ( S < l.size() ) ? S : l.size();
			  std::copy_n ( l.begin(), R, data.begin() ); }								// constructor from initialiser list

        vecN ( const vecN<T,S>& rhs ): data ( S )
            { std::copy_n ( rhs.data.begin(), S, data.begin() ); }						// copy constructor

		template <unsigned R>
        vecN ( vecN<T,R>& rhs ): data ( S )
            {  if (S>=R) std::copy ( rhs.getdata()->begin(), rhs.getdata()->end(), 
									 data.begin() ); }									// copy from vec2/3 ( for vec3/4 )

		vecN<T,S> operator=( vecN<T,S> rhs ) {
			std::copy( rhs.data.begin(), rhs.data.end(), data.begin() );				// assignment of vector
			return (*this);
		}

		vecN<T,S> operator=( T rhs ) {
			data = rhs;																	// assignment of scalar
			return (*this);
		}

		// passing on the [] operator for accessing elements
			  T& operator[] ( size_t offset)       { return data[offset]; }				// write access
		const T& operator[] ( size_t offset) const { return data[offset]; }				// read  access

		// member function for doing the same inside the class
			  T& at ( size_t offset)       { return data.at ( offset ); }				// write access
		const T& at ( size_t offset) const { return data.at ( offset ); }				// read  access

		// (in) equality
		bool operator== (vecN& v) { return ( data == v.data ); }						// equality
		bool operator!= (vecN& v) { return !( (*this) == v );  }						// inequality

		// some vector computations
		T norm() const { return ( sqrt( (*this) & (*this) ) ); } 						// norm: square root of inner product
		vecN<T,S> reciprocal ( const vecN& v ) const {									// reciprocal: 1/x for all elements x
			auto 
				out = (*this); 
			std::transform ( out.data.begin(), out.data.end(), out.data.begin(), std::bind1st ( std::divides<T>(), 1.0 ) ); 
			return out; 
		}

		// return std::vector with coefficients
		unsigned size() { return S; }													// size of data
		std::vector<T> *getdata() { return &data; }										// data vector

		// left += and + for elements and vectors
		const vecN<T,S>& operator+= ( const T& rhs ) { for ( auto i : std::views::iota(0u, S) ) data[i] += rhs; return (*this); }
		template <typename U>
		const vecN<T,S>& operator+= ( const vecN<U,S>& rhs ) { for ( auto i : std::views::iota(0u, S) ) data[i] += rhs[i]; return (*this); }
		template <typename U>
		const vecN<T,S> operator+   ( const U& rhs ) { auto out = (*this); out += rhs; return out;                   }

		// left -= and - for elements and vectors
		const vecN<T,S>& operator-= ( const T& rhs ) { for ( auto i : std::views::iota(0u, S) ) data[i] -= rhs;    return (*this); }
		template <typename U>
		const vecN<T,S>& operator-= ( const vecN<U,S>& rhs ) { for ( auto i : std::views::iota(0u, S) ) data[i] -= rhs[i]; return (*this); }
		template <typename U>
		const vecN<T,S> operator-   ( const U& rhs ) { auto out = (*this); out -= rhs; return out; }
		const vecN<T,S> operator-   ( void ) { auto out = (*this); out *= -1; return out; }

		// left *= and * for elements and vectors
		const vecN<T,S>& operator*= ( const T& rhs ) { for ( auto i : std::views::iota(0u, S) ) data[i] *= rhs;    return (*this); }
		template <typename U>
		const vecN<T,S>& operator*= ( const vecN<U,S>& rhs ) { for ( auto i : std::views::iota(0u, S) ) data[i] *= rhs[i]; return (*this); }

		// multiplication by scalar
		const vecN<T,S> operator*   ( const T& rhs ) { auto out = (*this); out *= rhs; return out;                   }

		// multiplication by vector
		template <typename U>
		const vecN<T,S> operator*   ( const vecN<U,S>& rhs ) { auto out = (*this); out *= rhs; return out;                   }

		// multiplication by matrix (if v is a row vector)
		template <typename U>
		vecN<T,S> operator*( const matN<U,S>& m ) const { 
			vecN<T,S> out;
			for ( auto i : std::views::iota(0u, S) )		// column in m
				for ( auto j : std::views::iota(0u, S) )	// column in (*this), row in m
					out.data[i] += data[i] * m [i][j];
			return out; 
		}

		// left /= and / for elements and vectors
		const vecN<T,S>& operator/= ( const T& rhs ) { for ( auto i : std::views::iota(0u, S) ) data[i] /= rhs;    return (*this); }
		template <typename U>
		const vecN<T,S>& operator/= ( const vecN<U,S>& rhs ) { for ( auto i : std::views::iota(0u, S) ) data[i] /= rhs[i]; return (*this); }
		template <typename U>
		const vecN<T,S> operator/   ( const U& rhs ) { vecN<T,S> out(*this); out /= rhs; return out;                   }

		// inner product (generalisation of dot product)
		T operator& ( const vecN<T,S>& v ) const { return std::inner_product ( data.begin(), data.end(), v.data.begin(), 0 ); }

		// cross product (exterior/wedge product in other than 3D but mostly used in 3D as cross product)
		vecN<T,S> operator^( const vecN<T,S>& v ) const { 
			switch (S) {
				case ( 2 ): // exterior product
					return vecN( data[0] * v[1] - data[1] * v[0], 0 );   
					break;
				case ( 3 ): // cross product
					return vecN( data[1] * v[2] - data[2] * v[1],
								 data[0] * v[2] - data[2] * v[0],
								 data[0] * v[1] - data[1] * v[0] );   
					break;
			default:
				std::cerr << "exterior / cross product only defined as scalar in 2D / 3D" << std::endl;
				return (*this);
			} // switch
		} // operator

};

// non-members of vecN for vecN
template <typename T, unsigned S>
std::ostream& operator<<(std::ostream& out, const vecN<T,S>& v) { 
	out << "( ";
	for ( auto i : std::views::iota(0u, S) )
		out << v.data[i] << ", "; 
	out << "\b\b )"; // use two ANSI backspace characters '\b' to overwrite final ", "
	return out; 
}

template <typename T, unsigned S>
std::istream& operator>>(std::istream& in , vecN<T,S> &v) { 
	for ( auto i : std::views::iota(0u, S) )
		in >> v.data[i]; 
	return in;
}

// right +, -, * and / operators
template <typename T, unsigned S>
inline vecN <T,S> operator+ ( T x, vecN <T,S> y) {
    return y + x;
}
template <typename T, unsigned S>
inline vecN <T,S> operator* ( T x, vecN <T,S> y) {
    return y * x;
}
template <typename T, unsigned S>
inline vecN <T,S> operator- ( T x, vecN <T,S> y) {
    return -y + x;
}
template <typename T, unsigned S>
inline vecN <T,S> operator/ ( T x, vecN <T,S> y) {
    return reciprocal(y) * x;
}



/* NxN matrix ( for N in 2,3,4 )
 *
 */

template <typename T, unsigned S>
class matN {

    template <typename T2, unsigned S2>
    friend std::ostream& operator<<(std::ostream& out, const matN<T2,S2>& v);

    template <typename T2, unsigned S2>
    friend std::istream& operator>>(std::istream& in, matN<T2,S2>& v);

    protected:
        std::vector<T> 
			data;
		unsigned
			sz = S;

    public:
        matN (                      ): data ( S*S )
            {};																	// default constructor
        matN ( T x0,   T y0,													// constructor from  4 scalars
               T x1,   T y1 ): data( S*S ) 
			{ if ( S>1 ) {
			  data[ 0] = x0; data[ 1] = y0;
			  data[ 2] = x1; data[ 3] = y1;			  
			  }
			}
        matN ( T x0,   T y0,   T z0,
               T x1,   T y1,   T z1,
               T x2,   T y2,   T z2 ): data ( S*S )
            { if ( S>2 ) {														// constructor from 16 scalars
			  data[ 0] = x0; data[ 1] = y0; data[ 2] = z0;
              data[ 3] = x1; data[ 4] = y1; data[ 5] = z1;
              data[ 6] = x2; data[ 7] = y2; data[ 8] = z2;
			  }
			}    
        matN ( T x0,   T y0,   T z0,   T t0,
               T x1,   T y1,   T z1,   T t1,
               T x2,   T y2,   T z2,   T t2,
               T x3,   T y3,   T z3,   T t3 ): data ( S*S )
            { if ( S>3 ) {														// constructor from 16 scalars
			  data[ 0] = x0; data[ 1] = y0; data[ 2] = z0; data[ 3] = t0;
              data[ 4] = x1; data[ 5] = y1; data[ 6] = z1; data[ 7] = t1;
              data[ 8] = x2; data[ 9] = y2; data[10] = z2; data[11] = t2;
              data[12] = x3; data[13] = y3; data[14] = z3; data[15] = t3;
			  }
			}    

        matN ( T* xyz                 ): data ( S*S )
            { std::copy (xyz, xyz+S*S, data.begin() ); }						// constructor from pointer

		matN ( std::initializer_list<T> l ): data ( S*S )
			{ auto R = ( S*S < l.size() ) ? S*S : l.size();
			  std::copy_n ( l.begin(), R, data.begin() ); }						// constructor from initialiser list
			  
        matN ( const matN<T,S>& rhs   ): data ( S*S )
            { std::copy ( rhs.data.begin(), rhs.data.end(), data.begin() ); }	// copy constructor

		matN<T,S> operator=( matN<T,S> rhs ) {
			std::copy( rhs.data.begin(), rhs.data.end(), data.begin() );		// assignment
			return (*this);
		}

		// passing on the [] operator for accessing 2D elements
			  T* operator[] (const size_t offset)       { return &data[S*offset]; }		// write access
		const T* operator[] (const size_t offset) const { return &data[S*offset]; }		// read  access

		// same functionality for inside the class
		      T& at(const size_t r, const size_t c )        { return data.at ( S*r + c ); } 
		const T& at(const size_t r, const size_t c ) const  { return data.at ( S*r + c ); }   

		// (in) equality
		bool operator== (matN& v) { return ( data == v.data ); }
		bool operator!= (matN& v) { return !( (*this) == v );                                   }

		// some matrix computations
		matN<T,S> reciprocal () const {														// reciprocal: 1/x for all elements x
			auto 
				out = (*this); 
			std::transform ( out.data.begin(), out.data.end(), out.data.begin(), std::bind1st ( std::divides<T>(), 1.0 ) ); 
			return out; 
		}
																		 																		
		const matN<T,S> eye ( const T v = 1 ) const {											// identity
			matN<T,S> out; 
			for ( auto i : std::views::iota(0u, S) )
				out[i][i] = v;
			return out;
		}

		// return std::vector with coefficients
		unsigned size() { return S; }														// size of data
		std::vector<T> *getdata() { return &data; }											// data vector

		// left += and + for elements and vectors
		const matN<T,S>& operator+= ( const T& rhs )       { for ( auto i : std::views::iota(0u, S*S) ) data[i] += rhs;    return (*this); }
		template <typename U>
		const matN<T,S>& operator+= ( const matN<U,S>& rhs ) { for ( auto i : std::views::iota(0u, S*S) ) data[i] += rhs.data[i]; return (*this); }

		template <typename U>
		const matN<T,S> operator+   ( const U& rhs )       { auto out = (*this); out += rhs; return out;                   }

		// left -= and - for elements and vectors
		const matN<T,S>& operator-= ( const T& rhs )       { for ( auto i : std::views::iota(0u, S*S) ) data[i] -= rhs; return (*this); }

		template <typename U>
		const matN<T,S>& operator-= ( const matN<U,S>& rhs ) { for ( auto i : std::views::iota(0u, S*S) ) data[i] -= rhs.data[i]; return (*this); }

		template <typename U>
		const matN<T,S> operator-   ( const U& rhs )       { matN<T,S> out(*this); out -= rhs; return out;                   }
		const matN<T,S> operator-   ( void )               { auto out = (*this); out *= -1; return out; }
		
		// matrix-matrix product (only for 2 equal size inputs)
		template <typename U>
		matN<T,S> operator*=( matN<U,S>& m ) { 			
			matN<T,S> product;
			for ( auto i : std::views::iota(0u, S) )
				for ( auto j : std::views::iota(0u, S) )
					for ( auto k : std::views::iota(0u, S) )
						product.data [ S*i + j ] += data [ S*i + k ] * m.data [ S*k + i ];
			return (*this);
		}

		// left *= for scalar
		const matN<T,S>& operator*= ( const T& rhs )       { for ( auto i : std::views::iota(0u, S*S) ) data[i] *= rhs; return (*this);    }
		
		// operator * (matrix pruduct for matrix rhs, element-wise for scalar rhs)
		template <typename U>
		const matN<T,S> operator*   ( matN<U,S>& rhs )       { auto out = (*this); out *= rhs; return out;                        }
		const matN<T,S> operator*   ( const T& rhs )       { auto out = (*this); out *= rhs; return out;                        }

		// hadamard product (element-wise multiplication)
		template <typename U>
		const matN<T,S>& hadamard ( const matN<U,S>& rhs ) { for ( auto i : std::views::iota(0u, S*S) ) data[i] *= rhs.data[i]; return (*this); }

		// left /= and / for elements and matrices
		const matN<T,S>& operator/= ( const T& rhs )       { for ( auto i : std::views::iota(0u, S*S) ) data[i] /= rhs;    return (*this); }
		template <typename U>
		const matN<T,S>& operator/= ( const matN<U,S>& rhs ) { (*this) *= rhs.inverse(); return (*this); }

		template <typename U>
		const matN<T,S> operator/   ( const U& rhs )       { auto out = (*this); out /= rhs; return out;                        }

		// cofactors -- required for adjoint (adjugate) matrix
		const T cofactor_ij ( size_t i, size_t j ) {			
			
			vecN<size_t,S> ii;
			vecN<size_t,S> jj;
			double fac;

			// fill with elements except i
			for (size_t k=0; k<i; k++) ii[k] = k;
			for (size_t k=i; k<3; k++) ii[k] = k+1;

			// fill with elements except j
			for (size_t k=0; k<j; k++) jj[k] = k;
			for (size_t k=j; k<3; k++) jj[k] = k+1;

			(fac)  = (*this)[ii[0]][jj[0]] * (  (*this)[ii[1]][jj[1]] * (*this)[ii[2]][jj[2]]
											 -  (*this)[ii[1]][jj[2]] * (*this)[ii[2]][jj[1]] );
			(fac) -= (*this)[ii[0]][jj[1]] * (  (*this)[ii[1]][jj[0]] * (*this)[ii[2]][jj[2]]
											 -  (*this)[ii[1]][jj[2]] * (*this)[ii[2]][jj[0]] );
			(fac) += (*this)[ii[0]][jj[2]] * (  (*this)[ii[1]][jj[0]] * (*this)[ii[2]][jj[1]]
											 -  (*this)[ii[1]][jj[1]] * (*this)[ii[2]][jj[0]] );

			/* compute sign */
			size_t k = i+j;
			if ( k != (k/2)*2) {
				(fac) = -(fac);
			}

			return T(fac);

		} // cofactor_ij

		// determinant -- 0 for singular matrices
		T determinant () {
			
			switch ( S ) {
				case ( 4 ): 
					return ( cofactor_ij ( 0, 0 ) * (*this)[0][0] +
							 cofactor_ij ( 0, 1 ) * (*this)[0][1] +
							 cofactor_ij ( 0, 2 ) * (*this)[0][2] +
							 cofactor_ij ( 0, 3 ) * (*this)[0][3] );
					break;				
				case ( 3 ):
					return (  data[0] * ( data[4] * data[8] - data[5] * data[7] )
							- data[1] * ( data[3] * data[8] - data[5] * data[6] )
							+ data[2] * ( data[3] * data[7] - data[4] * data[6] ) );
					break;
				default:
					return ( data[0] *  data[3] - data[1] * data[2] );
			} // switch

		} // determinant

		// adjoint (adjugate) matrix -- required for inverse
		matN<T,S> adjoint( const T scale = 1 ) {

			matN<T,S> out;
			switch ( S ) {
				case ( 4 ): 
					for ( auto i : std::views::iota(0u, S) )
						for ( auto j : std::views::iota(0u, S) )
							out[j][i] = scale * cofactor_ij ( i, j );
					return out;
					break;
				case ( 3 ):
					return matN ( scale * (data[4] * data[8] - data[5] * data[7]),
								 -scale * (data[1] * data[8] - data[2] * data[7]),
								  scale * (data[1] * data[5] - data[2] * data[4]),
								 -scale * (data[3] * data[8] - data[6] * data[5]),
								  scale * (data[0] * data[8] - data[2] * data[6]),
								 -scale * (data[0] * data[5] - data[2] * data[3]),
								  scale * (data[3] * data[7] - data[4] * data[6]),
								 -scale * (data[0] * data[7] - data[1] * data[6]),
								  scale * (data[0] * data[4] - data[1] * data[3]) ); 
					break;
				default:
					return matN ( scale * data[3], -scale * data[1],
								 -scale * data[2],  scale * data[0]  );
			} // switch

		}
    
    // inverse
    const matN<T,S> inverse() {
        T det = this->determinant();
        if ( det > FLT_EPSILON || det < -FLT_EPSILON )
            return ( this->adjoint ( 1 / det ) );
        else
            return (*this) * std::numeric_limits<T>::quiet_NaN();
    }

    ////////////////////////////////////////////////////////////////
	// code that omputes eigenvectors of ->only symmetric<- matrices
    ////////////////////////////////////////////////////////////////
	// 
	// the Jacobi algorithm applies successive Givens rotations
	// to the largest above-diagonal element, setting them to 0
	// until the matrix is diagonal. 


	// for a given row, find index of maximum after the diagonal
	unsigned max_col_upper ( unsigned row ) const {
	  unsigned c_max = row+1;
	  for ( unsigned c = row+2; c < S; c++ )
		if ( std::abs ( at ( row, c ) ) > std::abs ( at ( row, c_max ) ) )
		  c_max = c;
	  return c_max;
	}
	
	// find the maximum entry in the matrix M in O(n) time
	void max_pos_upper ( 	unsigned&				i_max, 
							unsigned&				j_max, 
							vecN<unsigned, S-1>&	max_idx_row ) const {
	  i_max = 0;
	  j_max = max_idx_row[i_max];
	  auto 
		max_entry = std::abs ( at ( i_max, j_max ) );
	  unsigned  
		nm1 = S-1;
		
	  for (unsigned i=1; i < nm1; i++) {
		unsigned j = max_idx_row[i];
		if ( std::abs ( at ( i, j ) ) > max_entry) {
		  max_entry = std::abs ( at ( i, j ) );
		  i_max = i;
		  j_max = j;
		} // if abs
	  } // for i

	} // maxEntry

	// Eigenvectors and -values FOR SYMMETRIC MATRICES ONLY
	// This is Jacobi's algorithm, only uses upper diagonal
	//      It diagonalises the matrix, so that the values on the diagonal
	//		are eigenvalues. The required rotations can also be applied to
	//		an identity matrix, to yield the corresponding eigenvectors. 	
	// 
	// based on https://github.com/jewettaij/jacobi_pd
	ev<T,S> diagonalise_sym ( 	bool vectors = true, 			// set to false for eigenvalues only
								unsigned max_iter = 100000 ) {	// decrease to escape slow convergence
		
		matN<T,S> 
			m (*this);							// local copy
		ev<T,S>
			ev;									// output eigenvalues and -vectors
		vecN<unsigned,S-1>
			max_idx_row;						// row 0..S-2: column idx of highest value

		if (vectors)
			ev.m=eye();

		for (unsigned i = 0; i < S-1; i++)			//Initialize the "max_idx_row[]" array 
			max_idx_row[i] = m.max_col_upper(i);	//(which is needed by max_pos_upper())
		
		for ( unsigned iter=0; iter <= max_iter; iter++ ) {
			unsigned i,j;
			m.max_pos_upper( i, j, max_idx_row );	// Find the maximum entry in the matrix. Store in i,j
			
			// If m[i][j] is small compared to m[i][i] and m[j][j], set it to 0 and update max.
			if (  ( m[i][i] + m[i][j] == m[i][i] ) && 
				  ( m[j][j] + m[i][j] == m[j][j] )  ) {
				m[i][j] = 0.0;
				max_idx_row[i] = m.max_col_upper(i);
			}  // if 0
	
			// if the maximum element is 0
			if ( m[i][j] == 0.0 )
				break;

			// Otherwise, apply a rotation to make M[i][j] = 0
			rotation 
				r = m.CalcRot ( i, j );				// Calculate the parameters of the rotation matrix.
			m.ApplyRot ( r, i, j, max_idx_row ); 	// Apply this rotation to the M matrix.
			if ( vectors )							// Optional: the eigenvectors are requested, then
				ev.m.ApplyRotLeft ( r, i, j );		// apply the rotation to the eigenvector matrix
		
		} // for iter 
		
		// copy m's diagonal as the eigenvalues
		for ( size_t r = 0; r<S; r++ ) 
			ev.v[r] = m[r][r]; 
		
		return ev;

	} // eigen_int

// calculate the rotation needed to zero the pivot
// and store the angle's cosine, sine and tangent
rotation CalcRot(	unsigned i,		// row index
					unsigned j) {	// column index
					
	rotation 
		r { 1, 1, 1 };
	double 
		M_jj_ii = at ( j, j ) - at ( i, i );
		
	if (M_jj_ii != 0.0) {
	
		r.t = 0.0;
		double 
			kappa = M_jj_ii,
			M_ij = at ( i, j );
			
		if (M_ij != 0.0) {
			
			kappa /= (2.0*M_ij);
			// t satisfies: t^2 + 2*t*kappa - 1 = 0
			// (choose the root which has the smaller absolute value)
			r.t = 1.0 / (std::sqrt(1 + kappa*kappa) + std::abs(kappa));
			if (kappa < 0.0)
				r.t = -r.t;
				
		} // if Mij
		
	} // if  Mjjii

	r.c = 1.0 / std::sqrt(1 + r.t*r.t);
	r.s = r.c*r.t;
	return r;
	
}

// apply the Givens rotation Q^T * M * Q
void ApplyRot(	rotation r, // angle
				unsigned i,		// row index
				unsigned j,		// column index
				vecN<unsigned,S-1>& max_idx_row ) {

	  // Recall that:
	  // c = cos(θ)
	  // s = sin(θ)
	  // t = tan(θ) (which should be <= 1.0)

	  // Compute the diagonal elements of M which have changed:
	  at(i,i) -= r.t * at(i,j);
	  at(j,j) += r.t * at(i,j);
	  // Note: This is algebraically equivalent to:
	  // M[i][i] = c*c*M[i][i] + s*s*M[j][j] - 2*s*c*M[i][j]
	  // M[j][j] = s*s*M[i][i] + c*c*M[j][j] + 2*s*c*M[i][j]

	  //Update the off-diagonal elements of M which will change (above the diagonal)

	  //assert(i < j);
	  at(i,j) = 0.0;

	  //compute M[w][i] and M[i][w] for all w!=i,considering above-diagonal elements
	  
	for (unsigned w=0; w < i; w++) {		// 0 <= w <  i  <  j < n
		at ( i, w ) = at ( w, i );		// backup the previous value. store below diagonal (i>w)
		at ( w, i ) = r.c * at ( w, i ) - r.s * at ( w, j );	//M[w][i], M[w][j] from previous iteration
		if ( i == max_idx_row[w] ) 
			max_idx_row[w] = max_col_upper(w);
			else if ( std::abs( at ( w, i ) ) > std::abs( at( w, max_idx_row[w]) ) )
				max_idx_row[w]=i;
		//assert(max_idx_row[w] == max_col_upper(M, w));
	} // for w	  
	for (unsigned w=i+1; w < j; w++) {		// 0 <= i <  w  <  j < n
		at ( w, i ) = at ( i, w );		// backup the previous value. store below diagonal (w>i)
		at ( i, w ) = r.c * at ( i, w ) - r.s * at ( w, j ); //M[i][w], M[w][j] from previous iteration
	} // for w
	for (unsigned w=j+1; w < S; w++) {      // 0 <= i < j+1 <= w < n
		at ( w, i ) = at ( i, w ); //backup the previous value. store below diagonal (w>i)
		at ( i, w ) = r.c * at ( i, w ) - r.s * at ( j, w ); //M[i][w], M[j][w] from previous iteration
	} // for w

	// now that we're done modifying row i, we can update max_idx_row[i]
	max_idx_row[i] = max_col_upper(i);

	//compute M[w][j] and M[j][w] for all w!=j,considering above-diagonal elements
	for ( unsigned w=0; w < i; w++ ) {			// 0 <=  w  <  i <  j < n
		at ( w, j ) = r.s * at ( i, w ) + r.c * at ( w, j ); // M[i][w], M[w][j] from previous iteration
		if ( j == max_idx_row[w]) 
			max_idx_row[w] = max_col_upper(w);
			else if ( std::abs( at ( w, j ) ) > std::abs( at ( w, max_idx_row[w] ) ) ) 
					max_idx_row[w]=j;
	} // for w
	for ( unsigned w=i+1; w < j; w++) {      // 0 <= i+1 <= w <  j < n
		at ( w, j ) = r.s*at(w,i) + r.c*at(w,j); //M[w][i], M[w][j] from previous iteration
		if (j == max_idx_row[w]) 
			max_idx_row[w] = max_col_upper(w);
			else if ( std::abs( at ( w, j ) ) > std::abs( at ( w, max_idx_row[w] ) ) )
				max_idx_row[w]=j;
		//assert(max_idx_row[w] == max_col_upper(M, w));
	} // for w
	for (unsigned w=j+1; w < S; w++) {      // 0 <=  i  <  j <  w < n
		at ( j, w ) = r.s * at ( w, i ) + r.c * at ( j, w ); //M[w][i], M[j][w] from previous iteration
	}
	// now that we're done modifying row j, we can update max_idx_row[j]
	max_idx_row[j] = max_col_upper(j);

} //Jacobi::ApplyRot()

// apply one Givens rotation 
void ApplyRotLeft( 	rotation r,	// angle 
					unsigned i,		// row index
					unsigned j) {	// column index

  for (unsigned v = 0; v < S; v++) {
	  
    auto 
		Miv = at ( i, v ); //backup E[i][v]
    at ( i, v ) = r.c * at ( i, v ) - r.s * at ( j, v );
    at ( j, v ) = r.s * Miv     	+ r.c * at ( j, v );
  
  } // for v

} // ApplyRotateLeft




}; // class matN

// non-members of matN for matN
template <typename U, unsigned S>
std::ostream& operator<<(std::ostream& out, const matN<U,S>& m) { 
	out << "( ";
	for ( auto i : std::views::iota(0u, S) ) {
		if ( i ) out << "  ";
		for ( auto j : std::views::iota(0u, S) )
			out << m[i][j] << ", ";
		out << "\b\b\n";
	}
	out << ")"; // use two ANSI backspace characters '\b' to overwrite final ", "
	return out; 
}

template <typename U, unsigned S>
std::istream& operator>>(std::istream& in , matN<U,S> &v) { 
	for ( auto i : std::views::iota(0u, S*S) )
		in >> v.data[i]; 
	return in;
}

// right +, -, * and / operators
template <typename T, unsigned S>
inline matN <T,S> operator+ ( T x, matN <T,S> y) {
    return y             + x;
}
template <typename T, unsigned S>
inline matN <T,S> operator* ( T x, matN <T,S> y) {
    return y             * x;
}
template <typename T, unsigned S>
inline matN <T,S> operator- ( T x, matN <T,S> y) {
    return -y            + x;
}
template <typename T, unsigned S>
inline matN <T,S> operator/ ( T x, matN <T,S> y) {
    return reciprocal(y) * x;
}

// matrix-vector product (v is a column vector)
template <typename T, typename U, unsigned S>
const vecN<T,S> operator* ( const matN<U,S>& m, const vecN<T,S>& v ) { 
	vecN<T,S> out;
	for ( auto i : std::views::iota(0u,S) ) 
		for ( auto j : std::views::iota(0u,S) ) 
			out[i] += m[i][j] * v[j]; 
	return out;
}

// hadamard product (element-wise multiplication)
template <typename T, typename U, unsigned S>
const matN<T,S> hadamard ( const matN<U,S>& m, const matN<T,S>& n ) { return m.hadamard(n); }



// aliases for 2-, 3- and 4-dimensional vectors
template <typename T>
using vec2 = vecN<T,2>;

template <typename T>
using vec3 = vecN<T,3>;

template <typename T>
using vec4 = vecN<T,4>;

template <typename T>
using mat2 = matN<T,2>;

template <typename T>
using mat3 = matN<T,3>;

template <typename T>
using mat4 = matN<T,4>;






























template<class T>
inline T fMaxOf2(T min, T max)
{
    return max > min ? max : min;
}

template<class T>
inline T fMinOf2(T min, T max)
{
    return max > min ? min : max;
}








/** \brief gauss samples the gauss curve
 *
 * given a position x and a width sigma
 */
template <typename T>
inline T gauss(T sigma, T x) {
    T expVal = - .5* pow( x/sigma, 2);
    T divider = sqrt(2 * M_PI * pow(sigma, 2));
    return (1 / divider) * exp(expVal);
}


/** \brief gausskernel in a vector
 *
 * length of the vector is 3.5 * sigma on either side
 * (gauss is .001 at 3.460871782016046838 times sigma)
 * also possible: 5.1 * sigma
 * (gauss is 10^-6 at 5.07869511287291208 times sigma)
 */
template <typename T>
const std::vector<T> gausskernel(T sigma) {

    double
        limit = 3.5;
    std::vector<T>
        out ( 1 + 2 * (unsigned) ceil( sigma * limit ) );

    // update limit as the 1st value on the x axis
    limit = -ceil (sigma * limit);

    // fill the Gaussian vector, whilst updating the x axis
    for (size_t i=0; i<out.size(); i++)
        out[i] = gauss<T> (sigma, limit++);

    return out;
}


/** \brief owfilter() returns an orthogonal wavelet filter given its name
 *
 * Filter is returned as a std::vector of doubles.
 * Only the scaling filter h is given. The wavelet
 * filter g is the QMF of h.
 */
const std::vector <double> owfilter ( const std::string name ) {

    // currently the orthogonal wavelet filters available are:
    // "Haar" and Daubechies2"

    std::vector <double> out;
    const double sq2=sqrt(2);

    if ( name == "Daubechies2" ) { // Daubechies 4-tap filter (2 vanishing moments)
        const double sq3 = sqrt(3);
        const double f2  = 4 * sq2;
        out = { (1 + sq3) / f2, (3 + sq3) / f2, (3 - sq3) / f2, (1 - sq3) / f2 };
    } else {
        out = { sq2, sq2 }; // Haar filter
    }

    return out;

}

} // namespace bis

#endif // BISIMAGE_MATHS_HPP_INCLUDED
