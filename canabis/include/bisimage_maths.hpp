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
#include<climits>
#include<iostream>
#include<cassert>

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



// forward declarations
template <class>
class vec2;
template <class>
class mat2;



/* 2-vector
 *
 * We define this mathematical object because it is often used in 2D scans
 */

template <typename T>
class vec2 {

    template <typename U>
    friend std::ostream& operator<<(std::ostream& out, const vec2<U>& v);

    template <typename U>
    friend std::istream& operator>>(std::istream& in, vec2<U>& v);

    protected:
        std::vector<T> data;

    public:
        vec2 (): data(2)
            {};                                                               // default constructor

        vec2 ( T x, T y ): data(2)
            { data[0] = x; data[1] = y; }                                     // constructor from 2 scalars

        vec2 ( T* xyz             ): data(2)
            { std::copy (xyz, xyz+2, data.begin() ); }                        // constructor from pointer

        vec2 ( const vec2<T>& rhs ): data(2)
            { std::copy ( rhs.data.begin(), rhs.data.end(), data.begin() ); } // copy constructor

    vec2<T> operator=( vec2<T> rhs ) {
        std::copy( rhs.data.begin(), rhs.data.end(), data.begin() );          // assignment of vector
        return (*this);
    }

    vec2<T> operator=( T rhs ) {
        data = rhs;                                                           // assignment of scalar
        return (*this);
    }

    // passing on the [] operator for accessing elements
          T& operator[] ( size_t offset)       { return data[offset]; }       // write access
    const T& operator[] ( size_t offset) const { return data[offset]; }       // read  access

    // (in) equality
    bool operator== (vec2& v) { return ( data[0]==v[0] && data[1]==v[1] ); }
    bool operator!= (vec2& v) { return !( (*this) == v );                  }

    // some vector computations
    T norm() { return ( sqrt( (*this) & (*this) ) ); } // norm: square root of inner product
    vec2<T> reciprocal ( const vec2& v ) const       { return vec2 ( 1 / data[0], 1 / data[1] );          }
    
    // return std::vector with coefficients
    std::vector<T> getdata() { return data; }

    // left += and + for elements and vectors
    const vec2<T>& operator+= ( const T& rhs )       { data[0]+=rhs;    data[1]+=rhs;    return (*this);  }
    template <typename U>
    const vec2<T>& operator+= ( const vec2<U>& rhs ) { data[0]+=rhs[0]; data[1]+=rhs[1]; return (*this);  }
    template <typename U>
    const vec2<T> operator+   ( const U& rhs )       { vec2<T> out(*this); out += rhs; return out;        }

    // left -= and - for elements and vectors
    const vec2<T>& operator-= ( const T& rhs )       { data[0]-=rhs;    data[1]-=rhs;    return (*this);  }
    template <typename U>
    const vec2<T>& operator-= ( const vec2<U>& rhs ) { data[0]-=rhs[0]; data[1]-=rhs[1]; return (*this);  }
    template <typename U>
    const vec2<T> operator-   ( const U& rhs )       { vec2<T> out(*this); out -= rhs; return out;        }
    const vec2<T> operator-   ( void )               { return vec2 ( -data[0], -data[1] );                }

    // left *= and * for elements and vectors
    const vec2<T>& operator*= ( const T& rhs )       { data[0]*=rhs;    data[1]*=rhs;    return (*this);  }
    template <typename U>
    const vec2<T>& operator*= ( const vec2<U>& rhs ) { data[0]*=rhs[0]; data[1]*=rhs[1]; return (*this);  }

    // multiplication by scalar
    const vec2<T> operator*   ( const T& rhs )       { vec2<T> out(*this); out *= rhs; return out;        }

    // multiplication by vector
    template <typename U>
    const vec2<T> operator*   ( const vec2<U>& rhs ) { vec2<T> out(*this); out *= rhs; return out;        }

    // multiplication by matrix (if v is a row vector)
    template <typename U>
    vec2<T> operator*( const mat2<U>& m ) const   { return vec2( data[0] * m[0][0] + data[1] * m[1][0],
                                                                 data[0] * m[0][1] + data[1] * m[1][1] ); }

    // left /= and / for elements and vectors
    const vec2<T>& operator/= ( const T& rhs )       { data[0]/=rhs;    data[1]/=rhs;    return (*this);  }
    template <typename U>
    const vec2<T>& operator/= ( const vec2<U>& rhs ) { data[0]/=rhs[0]; data[1]/=rhs[1]; return (*this);  }
    template <typename U>
    const vec2<T> operator/   ( const U& rhs )       { vec2<T> out(*this); out /= rhs; return out;        }

    // dot product (inner product in the more general case)
    T operator& ( const vec2& v ) const { return v[0] * data[0] + v[1] * data[1]; }

    // cross product (exterior/wedge product in other than 3D but mostly used in 3D as cross product)
    vec2<T> operator^( const vec2& v ) const   { return vec2( data[0] * v[1] - data[1] * v[0] );   }

};

// non-members of vec2 for vec2

template <typename U>
std::ostream& operator<<(std::ostream& out, const vec2<U>& v)
    { return out << '(' << v.data[0] << ',' << v.data[1] << ')'; }

template <typename U>
std::istream& operator>>(std::istream& in , vec2<U> &v)
    { return in >> v.data[0] >> v.data[1]; }

// right +, -, * and / operators
template <typename T>
inline vec2 <T> operator+ ( T x, vec2 <T> y) {
    return y             + x;
}
template <typename T>
inline vec2 <T> operator* ( T x, vec2 <T> y) {
    return y             * x;
}
template <typename T>
inline vec2 <T> operator- ( T x, vec2 <T> y) {
    return -y            + x;
}
template <typename T>
inline vec2 <T> operator/ ( T x, vec2 <T> y) {
    return reciprocal(y) * x;
}



/* 2x2 matrix
 *
 * We define this mathematical object because it is often used in 2D image applications
 */

template <typename T>
class mat2 {

    template <typename U>
    friend std::ostream& operator<<(std::ostream& out, const mat2<U>& v);

    template <typename U>
    friend std::istream& operator>>(std::istream& in, mat2<U>& v);

    protected:
        std::vector<T> data;

    public:
        mat2 (               ): data(4)
            {};                                                               // default constructor

        mat2 ( T x0,   T y0,
               T x1,   T y1  ): data(4)
            { data[0] = x0; data[1] = y0;
              data[2] = x1; data[3] = y1;               }                     // constructor from 4 scalars

        mat2 ( T* xyz               ): data(4)
            { std::copy (xyz, xyz+4, data.begin() ); }                        // constructor from pointer

        mat2 ( const mat2<T>& rhs   ): data(4)
            { std::copy ( rhs.data.begin(), rhs.data.end(), data.begin() ); } // copy constructor

        mat2<T> operator=( mat2<T> rhs ) {
            std::copy( rhs.data.begin(), rhs.data.end(), data.begin() );      // assignment
            return (*this);
        }

    // passing on the [] operator for accessing 2D elements
          T* operator[] (const size_t offset)       { return &data[2*offset]; }       // write access
    const T* operator[] (const size_t offset) const { return &data[2*offset]; }       // read  access

    // (in) equality
    bool operator== (mat2& v) { return ( data[0]==v[0] && data[1]==v[1] && 
                                         data[2]==v[2] && data[3]==v[3]  );                         }
    bool operator!= (mat2& v) { return !( (*this) == v );                                           }

    // some matrix computations
    mat2<T> reciprocal ( const mat2& v ) const       { return mat2 ( 1 / data[0], 1 / data[1], 
                                                                     1 / data[2], 1 / data[3]  );   }
    const mat2<T> eye ( const T v = 1 ) const        { return mat2 ( 1, 0, 0,
                                                                     0, 1, 0,
                                                                     0, 0, 1 );                     }

    // return std::vector with coefficients
    std::vector<T> getdata() { return data; }

    // left += and + for elements and vectors
    const mat2<T>& operator+= ( const T& rhs )       { data[0]+=rhs;    data[1]+=rhs;    
                                                       data[2]+=rhs;    data[3]+=rhs;
                                                       return (*this);                              }
    template <typename U>
    const mat2<T>& operator+= ( const mat2<U>& rhs ) { data[0]+=rhs.data[0]; data[1]+=rhs.data[1]; 
                                                       data[2]+=rhs.data[2]; data[3]+=rhs.data[3];
                                                       return (*this);                              }
    template <typename U>
    const mat2<T> operator+   ( const U& rhs )       { mat2<T> out(*this); out += rhs; return out;  }

    // left -= and - for elements and vectors
    const mat2<T>& operator-= ( const T& rhs )       { data[0]-=rhs;    data[1]-=rhs;    
                                                       data[2]-=rhs;    data[3]-=rhs;
                                                       return (*this);                              }
    template <typename U>
    const mat2<T>& operator-= ( const mat2<U>& rhs ) { data[0]-=rhs[0]; data[1]-=rhs[1]; 
                                                       data[2]-=rhs[2]; data[3]-=rhs[3];
                                                       return (*this);                              }
    template <typename U>
    const mat2<T> operator-   ( const U& rhs )       { mat2<T> out(*this); out -= rhs; return out;  }
    const mat2<T> operator-   ( void )               { return (*this) * -1;                         }

    // matrix-matrix product (only for 2 mat2-sized inputs)
    template <typename U>
    mat2<T> operator*=( mat2<U>& m ) { mat2<T> m2( m.data[0] * data[0] + m.data[2] * data[1],
                                                   m.data[1] * data[0] + m.data[3] * data[1],
                                                   m.data[0] * data[2] + m.data[2] * data[3],
                                                   m.data[1] * data[2] + m.data[3] * data[3]  );
                                       (*this) = m2;
                                       return (*this);                                              }

    // left *= for elements
    const mat2<T>& operator*= ( const T& rhs )       { data[0]*=rhs;    data[1]*=rhs;    
                                                       data[2]*=rhs;    data[3]*=rhs;
                                                       return (*this);                              }

    // operator *
    template <typename U>
    const mat2<T> operator*   ( mat2<U>& rhs )       { mat2<T> out(*this); out *= rhs; return out;  }
    const mat2<T> operator*   ( const T& rhs )       { mat2<T> out(*this); out *= rhs; return out;  }

    // hadamard product (element-wise multiplication)
    template <typename U>
    const mat2<T>& hadamard ( const mat2<U>& rhs ) { data[0]/=rhs.data[0]; data[1]/=rhs.data[1]; 
                                                     data[2]/=rhs.data[2]; data[3]/=rhs.data[3];
                                                     return (*this);                                }

    // left /= and / for elements and vectors
    const mat2<T>& operator/= ( const T& rhs )       { data[0]/=rhs; data[1]/=rhs;
                                                       data[2]/=rhs; data[3]/=rhs;
                                                       return (*this);                              }
    template <typename U>
    const mat2<T>& operator/= ( const mat2<U>& rhs ) { data[0]/=rhs.data[0]; data[1]/=rhs.data[1]; 
                                                       data[2]/=rhs.data[2]; data[3]/=rhs.data[3];
                                                       return (*this);                              }
    template <typename U>
    const mat2<T> operator/   ( const U& rhs )       { mat2<T> out(*this); out /= rhs; return out;  }

    // adjoint (adjugate) matrix -- required for inverse
    mat2<T> adjoint( const T scale = 1 ) const    { return mat2(  scale * data[3],
                                                                 -scale * data[1],
                                                                 -scale * data[2],
                                                                  scale * data[0]  );               }

    // determinant
    T determinant() const { return ( data[0] *  data[3] - data[1] * data[2] );                      }

    // inverse
    const mat2<T> inverse() {
        T det = this->determinant();
        if ( det > FLT_EPSILON || det < -FLT_EPSILON )
            return ( this->adjoint ( 1 / det ) );
        else
            return (*this) * std::numeric_limits<T>::quiet_NaN();
    }

};

// non-members of mat2 for mat2

template <typename U>
std::ostream& operator<<(std::ostream& out, const mat2<U>& m)
    { return out << '(' << m.data[0] << ',' << m.data[1] << ',' << std::endl
                 << ' ' << m.data[2] << ',' << m.data[3] << ')'; }

template <typename U>
std::istream& operator>>(std::istream& in , mat2<U> &m)
    { return in >> m.data[0] >> m.data[1]
                >> m.data[3] >> m.data[4]; }

// right +, -, * and / operators
template <typename T>
inline mat2 <T> operator+ ( T x, mat2 <T> y) {
    return y             + x;
}
template <typename T>
inline mat2 <T> operator* ( T x, mat2 <T> y) {
    return y             * x;
}
template <typename T>
inline mat2 <T> operator- ( T x, mat2 <T> y) {
    return -y            + x;
}
template <typename T>
inline mat2 <T> operator/ ( T x, mat2 <T> y) {
    return reciprocal(y) * x;
}

// matrix-vector product (vector is assumed to be a column vector)
template <typename T, typename U>
const vec2<T> operator* ( const mat2<U>& m, const vec2<T>& v ) { return vec2( v[0] * m[0][0] + v[1] * m[0][1],
                                                                              v[0] * m[1][0] + v[1] * m[1][1]  ); }

// hadamard product (element-wise multiplication)
template <typename T, typename U>
const mat2<T> hadamard ( const mat2<U>& m, const mat2<T>& n ) { return m.hadamard(n); }



/* Vectors and matrices in 3 dimensions
 *
 * For manipulating 3-dimensional (3D) coordinates, 3D matrices and
 * vectors, and their operations, are important tools. Examples are
 * voxel (3D pixel) coordinates in an MRI scan. To access and change
 * data at locations in these scans, it is important to optimise the
 * operations that handle 3D image space.
 */



// forward declarations
template <class>
class vec3;
template <class>
class mat3;



/* 3-vector
 *
 * We define this mathematical object because it is often used in 3D scans
 */

template <typename T>
class vec3 {

    template <typename U>
    friend std::ostream& operator<<(std::ostream& out, const vec3<U>& v);

    template <typename U>
    friend std::istream& operator>>(std::istream& in, vec3<U>& v);

    protected:
        std::vector<T> data;

    public:
        vec3 (                    ): data(3)
            {};                                                               // default constructor

        vec3 ( T x,   T y,   T z  ): data(3)
            { data[0] = x; data[1] = y; data[2] = z; }                        // constructor from 3 scalars

        vec3 ( T* xyz             ): data(3)
            { std::copy (xyz, xyz+3, data.begin() ); }                        // constructor from pointer

        vec3 ( const vec3<T>& rhs ): data(3)
            { std::copy ( rhs.data.begin(), rhs.data.end(), data.begin() ); } // copy constructor
            
        vec3 ( const vec2<T>& rhs ): data(3)
            { data[0] = rhs[0]; data[1]=rhs[1]; }                             // copy from vec2 (add a 0)

    vec3<T> operator=( vec3<T> rhs ) {
        std::copy( rhs.data.begin(), rhs.data.end(), data.begin() );          // assignment of vector
        return (*this);
    }

    vec3<T> operator=( T rhs ) {
        data = rhs;                                                           // assignment of scalar
        return (*this);
    }

    // passing on the [] operator for accessing elements
          T& operator[] ( size_t offset)       { return data[offset]; }       // write access
    const T& operator[] ( size_t offset) const { return data[offset]; }       // read  access

    // (in) equality
    bool operator== (vec3& v) { return ( data[0]==v[0] && data[1]==v[1] && data[2]==v[2] ); }
    bool operator!= (vec3& v) { return !( (*this) == v );                                   }

    // some vector computations
    T norm() { return ( sqrt( (*this) & (*this) ) ); } // norm: square root of inner product
    vec3<T> reciprocal ( const vec3& v ) const       { return vec3 ( 1 / data[0], 1 / data[1], 1 / data[2] );             }

    // return std::vector with coefficients
    std::vector<T> getdata() { return data; }

    // left += and + for elements and vectors
    const vec3<T>& operator+= ( const T& rhs )       { data[0]+=rhs;    data[1]+=rhs;    data[2]+=rhs;    return (*this); }
    template <typename U>
    const vec3<T>& operator+= ( const vec3<U>& rhs ) { data[0]+=rhs[0]; data[1]+=rhs[1]; data[2]+=rhs[2]; return (*this); }
    template <typename U>
    const vec3<T> operator+   ( const U& rhs )       { vec3<T> out(*this); out += rhs; return out;                        }

    // left -= and - for elements and vectors
    const vec3<T>& operator-= ( const T& rhs )       { data[0]-=rhs;    data[1]-=rhs;    data[2]-=rhs;    return (*this); }
    template <typename U>
    const vec3<T>& operator-= ( const vec3<U>& rhs ) { data[0]-=rhs[0]; data[1]-=rhs[1]; data[2]-=rhs[2]; return (*this); }
    template <typename U>
    const vec3<T> operator-   ( const U& rhs )       { vec3<T> out(*this); out -= rhs; return out;                        }
    const vec3<T> operator-   ( void )               { return vec3 ( -data[0], -data[1], -data[2] );                      }

    // left *= and * for elements and vectors
    const vec3<T>& operator*= ( const T& rhs )       { data[0]*=rhs;    data[1]*=rhs;    data[2]*=rhs;    return (*this); }
    template <typename U>
    const vec3<T>& operator*= ( const vec3<U>& rhs ) { data[0]*=rhs[0]; data[1]*=rhs[1]; data[2]*=rhs[2]; return (*this); }

    // multiplication by scalar
    const vec3<T> operator*   ( const T& rhs )       { vec3<T> out(*this); out *= rhs; return out;                        }

    // multiplication by vector
    template <typename U>
    const vec3<T> operator*   ( const vec3<U>& rhs ) { vec3<T> out(*this); out *= rhs; return out;                        }

    // multiplication by matrix (if v is a row vector)
    template <typename U>
    vec3<T> operator*( const mat3<U>& m ) const   { return vec3( data[0] * m[0][0] + data[1] * m[1][0] + data[2] * m[2][0],
                                                                 data[0] * m[0][1] + data[1] * m[1][1] + data[2] * m[2][1],
                                                                 data[0] * m[0][2] + data[1] * m[1][2] + data[2] * m[2][2] ); }

    // left /= and / for elements and vectors
    const vec3<T>& operator/= ( const T& rhs )       { data[0]/=rhs;    data[1]/=rhs;    data[2]/=rhs;    return (*this); }
    template <typename U>
    const vec3<T>& operator/= ( const vec3<U>& rhs ) { data[0]/=rhs[0]; data[1]/=rhs[1]; data[2]/=rhs[2]; return (*this); }
    template <typename U>
    const vec3<T> operator/   ( const U& rhs )       { vec3<T> out(*this); out /= rhs; return out;                        }

    // dot product (inner product in the more general case)
    T operator& ( const vec3& v ) const { return v[0] * data[0] + v[1] * data[1] + v[2] * data[2]; }

    // cross product (wedge product in other than 3D but mostly used in 3D as cross product)
    vec3<T> operator^( const vec3& v ) const   { return vec3( data[1] * v[2] - data[2] * v[1],
                                                              data[0] * v[2] - data[2] * v[0],
                                                              data[0] * v[1] - data[1] * v[0] );   }

};

// non-members of vec3 for vec3

template <typename U>
std::ostream& operator<<(std::ostream& out, const vec3<U>& v)
    { return out << '(' << v.data[0] << ',' << v.data[1] << ',' << v.data[2] <<')'; }

template <typename U>
std::istream& operator>>(std::istream& in , vec3<U> &v)
    { return in >> v.data[0] >> v.data[1] >> v.data[2]; }

// right +, -, * and / operators
template <typename T>
inline vec3 <T> operator+ ( T x, vec3 <T> y) {
    return y             + x;
}
template <typename T>
inline vec3 <T> operator* ( T x, vec3 <T> y) {
    return y             * x;
}
template <typename T>
inline vec3 <T> operator- ( T x, vec3 <T> y) {
    return -y            + x;
}
template <typename T>
inline vec3 <T> operator/ ( T x, vec3 <T> y) {
    return reciprocal(y) * x;
}



/* 3x3 matrix
 *
 * We define this mathematical object because it is often used in 3D scans
 */

template <typename T>
class mat3 {

    template <typename U>
    friend std::ostream& operator<<(std::ostream& out, const mat3<U>& v);

    template <typename U>
    friend std::istream& operator>>(std::istream& in, mat3<U>& v);

    protected:
        std::vector<T> data;

    public:
        mat3 (                      ): data(9)
            {};                                                               // default constructor

        mat3 ( T x0,   T y0,   T z0,
               T x1,   T y1,   T z1,
               T x2,   T y2,   T z2 ): data(9)
            { data[0] = x0; data[1] = y0; data[2] = z0;
              data[3] = x1; data[4] = y1; data[5] = z1;
              data[6] = x2; data[7] = y2; data[8] = z2; }                     // constructor from 9 scalars

        mat3 ( T* xyz               ): data(9)
            { std::copy (xyz, xyz+9, data.begin() ); }                        // constructor from pointer

        mat3 ( const mat3<T>& rhs   ): data(9)
            { std::copy ( rhs.data.begin(), rhs.data.end(), data.begin() ); } // copy constructor

        mat3 ( const mat2<T>& rhs   ): data(9)
            { data[0] = rhs[0][0]; data[1]=rhs[0][1]; 
              data[3] = rhs[1][0]; data[4]=rhs[1][1]; }                             // copy from mat2 (add col + row of 0)

    mat3<T> operator=( mat3<T> rhs ) {
        std::copy( rhs.data.begin(), rhs.data.end(), data.begin() );          // assignment
        return (*this);
    }

    // passing on the [] operator for accessing 2D elements
          T* operator[] (const size_t offset)       { return &data[3*offset]; }       // write access
    const T* operator[] (const size_t offset) const { return &data[3*offset]; }       // read  access

    // (in) equality
    bool operator== (mat3& v) { return ( data[0]==v[0] && data[1]==v[1] && data[2]==v[2] &&
                                         data[3]==v[3] && data[4]==v[4] && data[4]==v[4] &&
                                         data[6]==v[6] && data[7]==v[7] && data[8]==v[8] ); }
    bool operator!= (mat3& v) { return !( (*this) == v );                                   }

    // some matrix computations
    mat3<T> reciprocal ( const mat3& v ) const       { return mat3 ( 1 / data[0], 1 / data[1], 1 / data[2],
                                                                     1 / data[3], 1 / data[4], 1 / data[5],
                                                                     1 / data[6], 1 / data[7], 1 / data[8]  );            }
    const mat3<T> eye ( const T v = 1 ) const        { return mat3 ( 1, 0, 0,
                                                                     0, 1, 0,
                                                                     0, 0, 1 );                                           }

    // return std::vector with coefficients
    std::vector<T> getdata() { return data; }

    // left += and + for elements and vectors
    const mat3<T>& operator+= ( const T& rhs )       { data[0]+=rhs;    data[1]+=rhs;    data[2]+=rhs;
                                                       data[3]+=rhs;    data[4]+=rhs;    data[5]+=rhs;
                                                       data[6]+=rhs;    data[7]+=rhs;    data[8]+=rhs;
                                                       return (*this); }
    template <typename U>
    const mat3<T>& operator+= ( const mat3<U>& rhs ) { data[0]+=rhs.data[0]; data[1]+=rhs.data[1]; data[2]+=rhs.data[2];
                                                       data[3]+=rhs.data[3]; data[4]+=rhs.data[4]; data[5]+=rhs.data[5];
                                                       data[6]+=rhs.data[6]; data[7]+=rhs.data[7]; data[8]+=rhs.data[8];
                                                       return (*this); }
    template <typename U>
    const mat3<T> operator+   ( const U& rhs )       { mat3<T> out(*this); out += rhs; return out;                        }

    // left -= and - for elements and vectors
    const mat3<T>& operator-= ( const T& rhs )       { data[0]-=rhs;    data[1]-=rhs;    data[2]-=rhs;
                                                       data[3]-=rhs;    data[4]-=rhs;    data[5]-=rhs;
                                                       data[6]-=rhs;    data[7]-=rhs;    data[8]-=rhs;
                                                       return (*this); }
    template <typename U>
    const mat3<T>& operator-= ( const mat3<U>& rhs ) { data[0]-=rhs[0]; data[1]-=rhs[1]; data[2]-=rhs[2];
                                                       data[3]-=rhs[3]; data[4]-=rhs[4]; data[5]-=rhs[5];
                                                       data[6]-=rhs[6]; data[7]-=rhs[7]; data[8]-=rhs[8];
                                                       return (*this); }
    template <typename U>
    const mat3<T> operator-   ( const U& rhs )       { mat3<T> out(*this); out -= rhs; return out;                        }
    const mat3<T> operator-   ( void )               { return (*this) * -1;                                               }

    // matrix-matrix product (only for 2 mat3-sized inputs)
    template <typename U>
    mat3<T> operator*=( mat3<U>& m ) { mat3<T> m2( m.data[0] * data[0] + m.data[3] * data[1] + m.data[6] * data[2],
                                                   m.data[1] * data[0] + m.data[4] * data[1] + m.data[7] * data[2],
                                                   m.data[2] * data[0] + m.data[5] * data[1] + m.data[8] * data[2],
                                                   m.data[0] * data[3] + m.data[3] * data[4] + m.data[6] * data[5],
                                                   m.data[1] * data[3] + m.data[4] * data[4] + m.data[7] * data[5],
                                                   m.data[2] * data[3] + m.data[5] * data[4] + m.data[8] * data[5],
                                                   m.data[0] * data[6] + m.data[3] * data[7] + m.data[6] * data[8],
                                                   m.data[1] * data[6] + m.data[4] * data[7] + m.data[7] * data[8],
                                                   m.data[2] * data[6] + m.data[5] * data[7] + m.data[8] * data[8] );
                                       (*this) = m2;
                                       return (*this); }

    // left *= for elements
    const mat3<T>& operator*= ( const T& rhs )       { data[0]*=rhs;    data[1]*=rhs;    data[2]*=rhs;
                                                       data[3]*=rhs;    data[4]*=rhs;    data[5]*=rhs;
                                                       data[6]*=rhs;    data[7]*=rhs;    data[8]*=rhs;
                                                       return (*this); }

    // operator *
    template <typename U>
    const mat3<T> operator*   ( mat3<U>& rhs )       { mat3<T> out(*this); out *= rhs; return out;                        }
    const mat3<T> operator*   ( const T& rhs )       { mat3<T> out(*this); out *= rhs; return out;                        }

    // hadamard product (element-wise multiplication)
    template <typename U>
    const mat3<T>& hadamard ( const mat3<U>& rhs ) { data[0]/=rhs.data[0]; data[1]/=rhs.data[1]; data[2]/=rhs.data[2];
                                                     data[3]/=rhs.data[3]; data[4]/=rhs.data[4]; data[5]/=rhs.data[5];
                                                     data[6]/=rhs.data[6]; data[7]/=rhs.data[7]; data[8]/=rhs.data[8];
                                                     return (*this); }

    // left /= and / for elements and vectors
    const mat3<T>& operator/= ( const T& rhs )       { data[0]/=rhs;    data[1]/=rhs;    data[2]/=rhs;
                                                       data[3]/=rhs;    data[4]/=rhs;    data[5]/=rhs;
                                                       data[6]/=rhs;    data[7]/=rhs;    data[8]/=rhs;
                                                       return (*this); }
    template <typename U>
    const mat3<T>& operator/= ( const mat3<U>& rhs ) { data[0]/=rhs.data[0]; data[1]/=rhs.data[1]; data[2]/=rhs.data[2];
                                                       data[3]/=rhs.data[3]; data[4]/=rhs.data[4]; data[5]/=rhs.data[5];
                                                       data[6]/=rhs.data[6]; data[7]/=rhs.data[7]; data[8]/=rhs.data[8];
                                                       return (*this); }
    template <typename U>
    const mat3<T> operator/   ( const U& rhs )       { mat3<T> out(*this); out /= rhs; return out;                        }

    // adjoint (adjugate) matrix -- required for inverse
    mat3<T> adjoint( const T scale = 1 ) const    { return mat3(  scale * (data[4] * data[8] - data[5] * data[7]),
                                                                 -scale * (data[1] * data[8] - data[2] * data[7]),
                                                                  scale * (data[1] * data[5] - data[2] * data[4]),
                                                                 -scale * (data[3] * data[8] - data[6] * data[5]),
                                                                  scale * (data[0] * data[8] - data[2] * data[6]),
                                                                 -scale * (data[0] * data[5] - data[2] * data[3]),
                                                                  scale * (data[3] * data[7] - data[4] * data[6]),
                                                                 -scale * (data[0] * data[7] - data[1] * data[6]),
                                                                  scale * (data[0] * data[4] - data[1] * data[3]) ); }

    // determinant
    T determinant() const { return   data[0] * ( data[4] * data[8] - data[5] * data[7] )
                                   - data[1] * ( data[3] * data[8] - data[5] * data[6] )
                                   + data[2] * ( data[3] * data[7] - data[4] * data[6] ); }

    // inverse
    const mat3<T> inverse() {
        T det = this->determinant();
        if ( det > FLT_EPSILON || det < -FLT_EPSILON )
            return ( this->adjoint ( 1 / det ) );
        else
            return (*this) * std::numeric_limits<T>::quiet_NaN();
    }

};

// non-members of mat3 for mat3

template <typename U>
std::ostream& operator<<(std::ostream& out, const mat3<U>& m)
    { return out << '(' << m.data[0] << ',' << m.data[1] << ',' << m.data[2] <<        std::endl
                 << ' ' << m.data[3] << ',' << m.data[4] << ',' << m.data[5] <<        std::endl
                 << ' ' << m.data[6] << ',' << m.data[7] << ',' << m.data[8] << ')'; }

template <typename U>
std::istream& operator>>(std::istream& in , mat3<U> &v)
    { return in >> v.data[0] >> v.data[1] >> v.data[2]
                >> v.data[3] >> v.data[4] >> v.data[5]
                >> v.data[6] >> v.data[7] >> v.data[8]; }

// right +, -, * and / operators
template <typename T>
inline mat3 <T> operator+ ( T x, mat3 <T> y) {
    return y             + x;
}
template <typename T>
inline mat3 <T> operator* ( T x, mat3 <T> y) {
    return y             * x;
}
template <typename T>
inline mat3 <T> operator- ( T x, mat3 <T> y) {
    return -y            + x;
}
template <typename T>
inline mat3 <T> operator/ ( T x, mat3 <T> y) {
    return reciprocal(y) * x;
}

// matrix-vector product (vector is assumed to be a column vector)
template <typename T, typename U>
const vec3<T> operator* ( const mat3<U>& m, const vec3<T>& v ) { return vec3( v[0] * m[0][0] + v[1] * m[0][1] + v[2] * m[0][2],
                                                                              v[0] * m[1][0] + v[1] * m[1][1] + v[2] * m[1][2],
                                                                              v[0] * m[2][0] + v[1] * m[2][1] + v[2] * m[2][2] ); }

// hadamard product (element-wise multiplication)
template <typename T, typename U>
const mat3<T> hadamard ( const mat3<U>& m, const mat3<T>& n ) { return m.hadamard(n); }



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
template <class>
class vec4;
template <class>
class mat4;



/* 4-vector
 *
 * We define this mathematical object because it is often used in 4D scans
 */

template <typename T>
class vec4 {

    template <typename U>
    friend std::ostream& operator<<(std::ostream& out, const vec4<U>& v);

    template <typename U>
    friend std::istream& operator>>(std::istream& in, vec4<U>& v);

    protected:
        std::vector<T> data;

    public:
        vec4 (                    ): data(4)
            {};                                                               // default constructor

        vec4 ( T x, T y, T z, T t ): data(4)
            { data[0] = x; data[1] = y; data[2] = z; data[3] = t; }           // constructor from 4 scalars

        vec4 ( T* xyz             ): data(4)
            { std::copy (xyz, xyz+4, data.begin() ); }                        // constructor from pointer

        vec4 ( const vec4<T>& rhs ): data(4)
            { std::copy ( rhs.data.begin(), rhs.data.end(), data.begin() ); } // copy constructor

        vec4 ( const vec2<T>& rhs ): data(4)
            { data[0] = rhs[0]; data[1] = rhs[1]; }                           // copy from vec2 (add 2 0s)

        vec4 ( const vec3<T>& rhs ): data(4)
            { data[0] = rhs[0]; data[1] = rhs[1]; data[2] = rhs[2]; }         // copy from vec3 (add a 0)

    vec4<T> operator=( vec4<T> rhs ) {
        std::copy( rhs.data.begin(), rhs.data.end(), data.begin() );          // assignment of vector
        return (*this);
    }

    vec4<T> operator=( T rhs ) {
        data = rhs;                                                           // assignment of scalar
        return (*this);
    }

    // passing on the [] operator for accessing elements
          T& operator[] ( size_t offset)       { return data[offset]; }       // write access
    const T& operator[] ( size_t offset) const { return data[offset]; }       // read  access

    // (in) equality
    bool operator== (vec4& v) { return ( data[0]==v[0] && data[1]==v[1] && data[2]==v[2] && data[3]==v[3] );         }
    bool operator!= (vec4& v) { return !( (*this) == v );                                                            }

    // some vector computations
    T norm() { return ( sqrt( (*this) & (*this) ) ); } // norm: square root of inner product
    vec4<T> reciprocal ( const vec4& v ) const       { return vec4 ( 1 / data[0], 1 / data[1], 1 / data[2] );        }

    // return std::vector with coefficients
    std::vector<T> getdata() { return data; }

    // left += and + for elements and vectors
    const vec4<T>& operator+= ( const T& rhs )       { for ( auto i : {0,1,2,3} ) data[i] += rhs;    return (*this); }
    template <typename U>
    const vec4<T>& operator+= ( const vec4<U>& rhs ) { for ( auto i : {0,1,2,3} ) data[i] += rhs[i]; return (*this); }
    template <typename U>
    const vec4<T> operator+   ( const U& rhs )       { vec4<T> out(*this); out += rhs; return out;                   }

    // left -= and - for elements and vectors
    const vec4<T>& operator-= ( const T& rhs )       { for ( auto i : {0,1,2,3} ) data[i] -= rhs;    return (*this); }
    template <typename U>
    const vec4<T>& operator-= ( const vec4<U>& rhs ) { for ( auto i : {0,1,2,3} ) data[i] -= rhs[i]; return (*this); }
    template <typename U>
    const vec4<T> operator-   ( const U& rhs )       { vec4<T> out(*this); out -= rhs; return out;                   }
    const vec4<T> operator-   ( void )               { return vec4    ( -data[0], -data[1], -data[2], data[3] );     }

    // left *= and * for elements and vectors
    const vec4<T>& operator*= ( const T& rhs )       { for ( auto i : {0,1,2,3} ) data[i] *= rhs;    return (*this); }
    template <typename U>
    const vec4<T>& operator*= ( const vec4<U>& rhs ) { for ( auto i : {0,1,2,3} ) data[i] *= rhs[i]; return (*this); }

    // multiplication by scalar
    const vec4<T> operator*   ( const T& rhs )       { vec4<T> out(*this); out *= rhs; return out;                   }

    // multiplication by vector
    template <typename U>
    const vec4<T> operator*   ( const vec4<U>& rhs ) { vec4<T> out(*this); out *= rhs; return out;                   }

    // multiplication by matrix (if v is a row vector)
    template <typename U>
    vec4<T> operator*( const mat4<U>& m ) const   { return vec4(   data[0] * m[0][0] + data[1] * m[1][0]
                                                                 + data[2] * m[2][0] + data[2] * m[3][0],
                                                                   data[0] * m[0][1] + data[1] * m[1][1]
                                                                 + data[2] * m[2][1] + data[3] * m[3][1],
                                                                   data[0] * m[0][2] + data[1] * m[1][2]
                                                                 + data[2] * m[2][2] + data[3] * m[3][2],
                                                                   data[0] * m[0][3] + data[1] * m[1][3]
                                                                 + data[2] * m[2][3] + data[3] * m[3][3] ); }

    // left /= and / for elements and vectors
    const vec4<T>& operator/= ( const T& rhs )       { for ( auto i : {0,1,2,3} ) data[i] /= rhs;    return (*this); }
    template <typename U>
    const vec4<T>& operator/= ( const vec4<U>& rhs ) { for ( auto i : {0,1,2,3} ) data[i] /= rhs[i]; return (*this); }
    template <typename U>
    const vec4<T> operator/   ( const U& rhs )       { vec4<T> out(*this); out /= rhs; return out;                   }

    // inner product (generalisation of dot product)
    T operator& ( const vec4& v ) const { return v[0] * data[0] + v[1] * data[1] + v[2] * data[2] + v[3] * data[3];  }

};

// non-members of vec4 for vec4

template <typename U>
std::ostream& operator<<(std::ostream& out, const vec4<U>& v)
    { return out << '(' << v.data[0] << ',' << v.data[1] << ',' << v.data[2] << ',' << v.data[3] <<')'; }

template <typename U>
std::istream& operator>>(std::istream& in , vec4<U> &v)
    { return in >> v.data[0] >> v.data[1] >> v.data[2] >> v.data[3]; }

// right +, -, * and / operators
template <typename T>
inline vec4 <T> operator+ ( T x, vec4 <T> y) {
    return y             + x;
}
template <typename T>
inline vec4 <T> operator* ( T x, vec4 <T> y) {
    return y             * x;
}
template <typename T>
inline vec4 <T> operator- ( T x, vec4 <T> y) {
    return -y            + x;
}
template <typename T>
inline vec4 <T> operator/ ( T x, vec4 <T> y) {
    return reciprocal(y) * x;
}



/* 4x4 matrix
 *
 * We define this mathematical object because it is often used in 4D scans
 */

template <typename T>
class mat4 {

    template <typename U>
    friend std::ostream& operator<<(std::ostream& out, const mat4<U>& v);

    template <typename U>
    friend std::istream& operator>>(std::istream& in, mat4<U>& v);

    protected:
        std::vector<T> data;

    public:
        mat4 (                      ): data(16)
            {};                                                               // default constructor

        mat4 ( T x0,   T y0,   T z0,   T t0,
               T x1,   T y1,   T z1,   T t1,
               T x2,   T y2,   T z2,   T t2,
               T x3,   T y3,   T z3,   T t3 ): data(16)
            { data[ 0] = x0; data[ 1] = y0; data[ 2] = z0; data[ 3] = t0;
              data[ 4] = x1; data[ 5] = y1; data[ 6] = z1; data[ 7] = t1;
              data[ 8] = x2; data[ 9] = y2; data[10] = z2; data[11] = t2;
              data[12] = x3; data[13] = y3; data[14] = z3; data[15] = t3;}    // constructor from 16 scalars

        mat4 ( T* xyz               ): data(16)
            { std::copy (xyz, xyz+16, data.begin() ); }                       // constructor from pointer

        mat4 ( const mat4<T>& rhs   ): data(16)
            { std::copy ( rhs.data.begin(), rhs.data.end(), data.begin() ); } // copy constructor

        mat4 ( const mat2<T>& rhs   ): data(16)
            { data[0] = rhs[0][0]; data[1]=rhs[0][1]; 
              data[4] = rhs[1][0]; data[5]=rhs[1][1]; }                       // copy from mat2 (add cols + rows of 0)

        mat4 ( const mat3<T>& rhs   ): data(16)
            { data[0] = rhs[0][0]; data[1]=rhs[0][1]; data[ 2]=rhs[0][2];
              data[4] = rhs[1][0]; data[5]=rhs[1][1]; data[ 6]=rhs[1][2];
              data[8] = rhs[2][0]; data[9]=rhs[2][1]; data[10]=rhs[1][2]; }   // copy from mat3 (add cols + rows of 0)

    mat4<T> operator=( mat4<T> rhs ) {
        std::copy( rhs.data.begin(), rhs.data.end(), data.begin() );          // assignment
        return (*this);
    }

    // passing on the [] operator for accessing 2D elements
          T* operator[] (const size_t offset)       { return &data[4*offset]; }       // write access
    const T* operator[] (const size_t offset) const { return &data[4*offset]; }       // read  access

    // (in) equality
    bool operator== (mat4& v) { return ( data[ 0]==v[ 0] && data[ 1]==v[ 1] && data[ 2]==v[ 2] && data[ 3]==v[ 3] &&
                                         data[ 4]==v[ 4] && data[ 5]==v[ 5] && data[ 6]==v[ 6] && data[ 7]==v[ 7] &&
                                         data[ 8]==v[ 8] && data[ 9]==v[ 9] && data[10]==v[10] && data[11]==v[11] &&
                                         data[12]==v[12] && data[13]==v[13] && data[14]==v[14] && data[15]==v[15] ); }
    bool operator!= (mat4& v) { return !( (*this) == v );                                   }

    // some matrix computations
    mat4<T> reciprocal ( const mat4& v ) const       { return mat4 ( 1 / data[ 0], 1 / data[ 1], 1 / data[ 2], 1 / data[ 3],
                                                                     1 / data[ 4], 1 / data[ 5], 1 / data[ 6], 1 / data[ 7],
                                                                     1 / data[ 8], 1 / data[ 9], 1 / data[10], 1 / data[11],
                                                                     1 / data[12], 1 / data[13], 1 / data[14], 1 / data[15]  );            }
    const mat4<T> eye ( const T v = 1 ) const        { return mat4 ( 1, 0, 0, 0,
                                                                     0, 1, 0, 0,
                                                                     0, 0, 1, 0,
                                                                     0, 0, 0, 1 );                                           }

    // return std::vector with coefficients
    std::vector<T> getdata() { return data; }

    // left += and + for elements and vectors
    const mat4<T>& operator+= ( const T& rhs )       { for ( int i=0; i<16; i++ ) data[i] += rhs;    return (*this); }
    template <typename U>
    const mat4<T>& operator+= ( const mat4<U>& rhs ) { for ( int i=0; i<16; i++ ) data[i] += rhs.data[i]; return (*this); }

    template <typename U>
    const mat4<T> operator+   ( const U& rhs )       { mat4<T> out(*this); out += rhs; return out;                   }

    // left -= and - for elements and vectors
    const mat4<T>& operator-= ( const T& rhs )       { for ( int i=0; i<16; i++ ) data[i] -= rhs; return (*this); }

    template <typename U>
    const mat4<T>& operator-= ( const mat4<U>& rhs ) { for ( int i=0; i<16; i++ ) data[i] -= rhs.data[i]; return (*this); }

    template <typename U>
    const mat4<T> operator-   ( const U& rhs )       { mat4<T> out(*this); out -= rhs; return out;                   }
    const mat4<T> operator-   ( void )               { return mat4 ( -data[ 0], -data[ 1], -data[ 2], -data[ 3]
                                                                     -data[ 4], -data[ 5], -data[ 6], -data[ 7]
                                                                     -data[ 8], -data[ 9], -data[10], -data[11]
                                                                     -data[12], -data[13], -data[14], -data[15] );        }
    // matrix-matrix product (only for 2 mat4-sized inputs)
    template <typename U>
    mat4<T> operator*=( mat4<U>& m ) { mat4<T> m2( m.data[ 0] * data[ 0] + m.data[ 4] * data[ 1] + m.data[ 8] * data[ 2] + m.data[12] * data[ 3],
                                                   m.data[ 1] * data[ 0] + m.data[ 5] * data[ 1] + m.data[ 9] * data[ 2] + m.data[13] * data[ 3],
                                                   m.data[ 2] * data[ 0] + m.data[ 6] * data[ 1] + m.data[10] * data[ 2] + m.data[14] * data[ 3],
                                                   m.data[ 3] * data[ 0] + m.data[ 7] * data[ 1] + m.data[11] * data[ 2] + m.data[15] * data[ 3],
                                                   m.data[ 0] * data[ 4] + m.data[ 4] * data[ 5] + m.data[ 8] * data[ 6] + m.data[12] * data[ 7],
                                                   m.data[ 1] * data[ 4] + m.data[ 5] * data[ 5] + m.data[ 9] * data[ 6] + m.data[13] * data[ 7],
                                                   m.data[ 2] * data[ 4] + m.data[ 6] * data[ 5] + m.data[10] * data[ 6] + m.data[14] * data[ 7],
                                                   m.data[ 3] * data[ 4] + m.data[ 7] * data[ 5] + m.data[11] * data[ 6] + m.data[15] * data[ 7],
                                                   m.data[ 0] * data[ 8] + m.data[ 4] * data[ 9] + m.data[ 8] * data[10] + m.data[12] * data[11],
                                                   m.data[ 1] * data[ 8] + m.data[ 5] * data[ 9] + m.data[ 9] * data[10] + m.data[13] * data[11],
                                                   m.data[ 2] * data[ 8] + m.data[ 6] * data[ 9] + m.data[10] * data[10] + m.data[14] * data[11],
                                                   m.data[ 3] * data[ 8] + m.data[ 7] * data[ 9] + m.data[11] * data[10] + m.data[15] * data[11],
                                                   m.data[ 0] * data[12] + m.data[ 4] * data[13] + m.data[ 8] * data[14] + m.data[12] * data[15],
                                                   m.data[ 1] * data[12] + m.data[ 5] * data[13] + m.data[ 9] * data[14] + m.data[13] * data[15],
                                                   m.data[ 2] * data[12] + m.data[ 6] * data[13] + m.data[10] * data[14] + m.data[14] * data[15],
                                                   m.data[ 3] * data[12] + m.data[ 7] * data[13] + m.data[11] * data[14] + m.data[15] * data[15]  );
                                       (*this) = m2;
                                       return (*this); }

    // left *= for elements
    const mat4<T>& operator*= ( const T& rhs )       { for ( int i=0; i<16; i++ ) data[i] *= rhs; return (*this);    }
    // operator *
    template <typename U>
    const mat4<T> operator*   ( mat4<U>& rhs )       { mat4<T> out(*this); out *= rhs; return out;                        }
    const mat4<T> operator*   ( const T& rhs )       { mat4<T> out(*this); out *= rhs; return out;                        }

    // hadamard product (element-wise multiplication)
    template <typename U>
    const mat4<T>& hadamard ( const mat4<U>& rhs )   { for ( int i=0; i<16; i++ ) data[i] *= rhs.data[i]; return (*this); }

    // left /= and / for elements and vectors
    const mat4<T>& operator/= ( const T& rhs )       { for ( int i=0; i<16; i++ ) data[i] /= rhs;    return (*this); }
    template <typename U>
    const mat4<T>& operator/= ( const mat4<U>& rhs ) { for ( int i=0; i<16; i++ ) data[i] /= rhs.data[i]; return (*this); }

    template <typename U>
    const mat4<T> operator/   ( const U& rhs )       { mat4<T> out(*this); out /= rhs; return out;                        }

    // cofactors -- required for adjoint (adjugate) matrix
    const T cofactor_ij ( size_t i, size_t j ) {

        vec4<size_t> ii;
        vec4<size_t> jj;
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
    const T determinant () {

    T det =
        cofactor_ij ( 0, 0 ) * (*this)[0][0] +
        cofactor_ij ( 0, 1 ) * (*this)[0][1] +
        cofactor_ij ( 0, 2 ) * (*this)[0][2] +
        cofactor_ij ( 0, 3 ) * (*this)[0][3];

    return det;

    } // determinant

    // adjoint (adjugate) matrix -- required for inverse
    mat4<T> adjoint( const T scale = 1 ) {

        mat4<T> out;

        for ( int i=0; i<4; i++ )
            for ( int j=0; j<4; j++ )
                out[j][i] = scale * cofactor_ij ( i, j );

        return out;

    }

    // inverse
    const mat4<T> inverse() {
        T det = this->determinant();
        if ( det > FLT_EPSILON || det < -FLT_EPSILON )
            return ( this->adjoint ( 1 / det ) );
        else
            return (*this) * std::numeric_limits<T>::quiet_NaN();
    }

};

// non-members of mat4 for mat4

template <typename U>
std::ostream& operator<<(std::ostream& out, const mat4<U>& m)
    { return out << '(' << m.data[ 0] << ',' << m.data[ 1] << ',' << m.data[ 2] << ',' << m.data[ 3] << std::endl
                 << ' ' << m.data[ 4] << ',' << m.data[ 5] << ',' << m.data[ 6] << ',' << m.data[ 7] << std::endl
                 << ' ' << m.data[ 8] << ',' << m.data[ 9] << ',' << m.data[10] << ',' << m.data[11] << std::endl
                 << ' ' << m.data[12] << ',' << m.data[13] << ',' << m.data[14] << ',' << m.data[15] << ')'; }

template <typename U>
std::istream& operator>>(std::istream& in , mat4<U> &v)
    { return in >> v.data[ 0] >> v.data[ 1] >> v.data[ 2] >> v.data[ 3]
                >> v.data[ 4] >> v.data[ 5] >> v.data[ 6] >> v.data[ 7]
                >> v.data[ 8] >> v.data[ 9] >> v.data[10] >> v.data[11]
                >> v.data[12] >> v.data[13] >> v.data[14] >> v.data[15]; }

// right +, -, * and / operators
template <typename T>
inline mat4 <T> operator+ ( T x, mat4 <T> y) {
    return y             + x;
}
template <typename T>
inline mat4 <T> operator* ( T x, mat4 <T> y) {
    return y             * x;
}
template <typename T>
inline mat4 <T> operator- ( T x, mat4 <T> y) {
    return -y            + x;
}
template <typename T>
inline mat4 <T> operator/ ( T x, mat4 <T> y) {
    return reciprocal(y) * x;
}

// matrix-vector product (vector is assumed to be a column vector)
template <typename T, typename U>
const vec4<T> operator* ( const mat4<U>& m, const vec4<T>& v ) { return vec4( v[0] * m[0][0] + v[1] * m[0][1] + v[2] * m[0][2] + v[3] * m[0][3],
                                                                              v[0] * m[1][0] + v[1] * m[1][1] + v[2] * m[1][2] + v[3] * m[1][3],
                                                                              v[0] * m[2][0] + v[1] * m[2][1] + v[2] * m[2][2] + v[3] * m[2][3],
                                                                              v[0] * m[3][0] + v[1] * m[3][1] + v[2] * m[3][2] + v[3] * m[3][3] ); }

// hadamard product (element-wise multiplication)
template <typename T, typename U>
const mat4<T> hadamard ( const mat4<U>& m, const mat4<T>& n ) { return m.hadamard(n); }

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
