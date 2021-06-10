#ifndef BISIMAGE_HPP_INCLUDED
#define BISIMAGE_HPP_INCLUDED



/** \brief Cimg.h: header-only C++ library for handling pictures
 *
 * Downloaded from https://framagit.org/dtschump/CImg/raw/master/CImg.h
 * it can save images as BMP picturs without requiring extra libraries
 */
#include "CImg/CImg.h"



/** \brief bisimage_maths.hpp: maths often used in images
 *
 * Including Gaussian filters/noise and 3D/4D vectors
 */
#include "bisimage_maths.hpp"



/** \brief bisimage_types.hpp: basic input types, enumerations and exceptions
 *
 * Including Gaussian filters/noise and 3D/4D vectors
 */
#include "bisimage_types.hpp"



/** \brief Also using the full benefit of standard library functionality
 *
 */
#include <queue>
#include <random>
#include <limits>
#include <numeric>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <type_traits>
#include <initializer_list>



/** \brief bis: namespace for multidimensional images
 *
 *  This namespace contains the bisimage class for
 *  creating and processing multidimensional images.
 */
namespace bis {



/** Functions outside of the bisimage class.
 *  These can be used by bisimage objects, 
 *  but are not required by the class and do 
 *  not require it themselves.
 *
 */
bool readablefile ( const char *filename );
template <typename T, typename U>
    void filterline ( std::vector<T>& signal, const std::vector<U>& filter );
template <typename T, typename U>
    void fwtline ( std::vector<T>& signal, std::vector<U> filterh, size_t level );
template <typename T, typename U>
    void ifwtline ( std::vector<T>& signal, std::vector<U> filterh, size_t level );
template <typename T, typename U>
    void fwtline ( std::vector<T>& signal, std::vector<U> filterh, size_t level );
template <typename T, typename U>
    void ifwtline ( std::vector<T>& signal, std::vector<U> filterh, size_t level );



/** \brief pixelImport -- memory buffer to typename
 *         
 *  template parameters:
 *  BufElem:    the element type of the input buffer (e.g., float)
 *  Data:       the target container (e.g., a vector<double>)
 * 
 *  function parameters:
 *  pixelbuffer:    pointer to the start of the input buffer
 *  bufsize:        length of the input buffer (in elements)
 *  data:           pointer to the start of the target container's buffer
 * 
 */
template<typename BufElem, typename Data>

void pixelImport ( void* pixelbuffer, unsigned bufsize, Data *data ) {

  typedef typename Data::value_type value_type;

  BufElem*
    buffer = ( BufElem* ) pixelbuffer;
  Data
    &dataref = *data;

  for ( unsigned i = 0; i < bufsize; i++ )
    dataref[i] = ( value_type ) ( buffer[i] );

};



/** \brief class template bisimage for n-dimensional images
 *
 * voxel type T can be any scalar type (char, short, int, float, double)
 */
template <typename value_type>
class bisimage {

    
    
	// self (no superclass)
	using self = bisimage<value_type>;
  

    
	protected:
  
    /** \brief the header, data, sizes and strides should
     *         be usable by subclasses bisnifti and bisdicom
     * 
     * data:    vector with the pixel data
     * sizes:   vector with the dimensions
     * strides: vector with the strides through
     *          the data in each dimension
     *
     */
    std::vector <size_t>
    sizes;

    std::vector <size_t>
    strides;

    std::vector <value_type>
    data;

    /** \brief validate sizes
     *
     * After the sizes vector has changed, update the other vectors
     */   
    void validate_sizes() {         
    
        strides.resize ( sizes.size()+1 );

        strides [ 0 ] = 1;
        for ( size_t i=1; i<=sizes.size(); i++ ) {
            strides [i] = strides [ i - 1 ] * sizes [ i - 1 ];
        }

        data.resize ( *strides.rbegin() );
        
    }



	public:

	/** \brief default constructor
     */
    bisimage() {}

    /** \brief destructor
     *
     * destructor: clears vectors and (if
     * required) frees pointer to header
     */
    ~bisimage() {
        
        data.resize    (0);
        sizes.resize   (0);
        strides.resize (0);

    }

    /** \brief constructor for an empty image
     *                            with sizes
     *
     * constructor for an empty image (no
     * nifti header) with given dimensions
     */
    bisimage ( std::initializer_list <size_t> dims ) {

        newsizes ( dims );
        
    }



    /** \brief (deep) copy constructor
     *
     * copies data, sizes and header from an existing bisimage
     */
    bisimage ( const bisimage& rhs ) {

        data.resize( rhs.data.size() );
        sizes.resize( rhs.sizes.size() );
        strides.resize( rhs.strides.size() );

        std::copy ( rhs.data.begin(),    rhs.data.end(),    data.begin()    );
        std::copy ( rhs.sizes.begin(),   rhs.sizes.end(),   sizes.begin()   );
        std::copy ( rhs.strides.begin(), rhs.strides.end(), strides.begin() );

    }



    /** \brief assignment operator
     *
     * assigns data, sizes and strides of the
     * right-hand side (RHS) image to (*this)
     */
    const bisimage<value_type>& operator= ( bisimage <value_type>& rhs ) {

        if ( this != &rhs ) {

            // just to make sure we don't leave stuff
            data.resize    ( rhs.data.size()    );
            sizes.resize   ( rhs.sizes.size()   );
            strides.resize ( rhs.strides.size() );

            std::copy ( rhs.data.begin(),    rhs.data.end(),    data.begin()    );
            std::copy ( rhs.sizes.begin(),   rhs.sizes.end(),   sizes.begin()   );
            std::copy ( rhs.strides.begin(), rhs.strides.end(), strides.begin() );

        } // if this != rhs

        return *this;

    } // assignment



    /** \brief initialise from a flat array (as a vector)
     *
     * The reason that this is not a constructor is to avoid 
     * confusion with the constructor that initialises size
     * 
     */
    void array_init ( std::vector < value_type > array ) {

            data.resize( array.size() );
            sizes.resize   ( 1 );
            strides.resize ( 2 );
            strides [0] = 1;
            strides [1] = sizes [0] = data.size();

            std::copy ( std::begin(array), std::end(array), data.begin() );

    }



    /** \brief initialise from a flat array (as an initializer_list)
     *
     * The reason that this is not a constructor is to avoid 
     * confusion with the constructor that initialises size
     * 
     */
    void array_init ( std::initializer_list < value_type > array ) {

            data.resize( array.size() );
            sizes.resize   ( 1 );
            strides.resize ( 2 );
            strides [0] = 1;
            strides [1] = sizes [0] = data.size();

            std::copy ( std::begin(array), std::end(array), data.begin() );

    }



    /** \brief initialise from a flat array (c-style)
     *
     * The reason thatthis is not a constructor is to avoid 
     * confusion with the constructor that initialises size
     * 
     */
    template <size_t array_size>
    void array_init ( value_type (&array) [array_size] ) {

            data.resize( array_size );
            sizes.resize   ( 1 );
            strides.resize ( 2 );
            strides [0] = 1;
            strides [1] = sizes [0] = data.size();

            std::copy ( std::begin(array), std::end(array), data.begin() );

    }



    /** \brief newsizes -- changes the image dimensions
     *
     * Resizes data, sizes and strides.
     * WARNING this may lead to undefined data!!!
     * 
     * You probably need to read new data before continuing.
     * (that's why it's newsize not resize)
     * 
     */
    void newsizes ( std::initializer_list<size_t> newdims ) {

        sizes.resize  ( newdims.size() );
        std::copy ( newdims.begin(), newdims.end(), sizes.begin() );        
        validate_sizes();

    }



    /** \brief operator() for positional addressing
     *
     * for an existing image In, this operator can
     * be used with multidimensional co-ordinates to
     * indicate position, so In({x,y,z}) instead of
     * In.data[x + y*sizes[1] + z*sizes[1]*sizes[2]].
     *
     * This operator is for reading only.
     */
    inline value_type const operator() ( std::initializer_list < size_t > const& indices ) const {
        size_t const offset =
            std::inner_product ( indices.begin(), indices.end(),
                                 strides.begin(),
                                 0 );
        return data [ offset ];
    }

    /** \brief operator() for positional addressing
     *
     * This operator is for modifying the data.
     */
    inline value_type& operator() ( std::initializer_list < size_t > const& indices ) {
        size_t const offset =
            std::inner_product ( indices.begin(), indices.end(),
                                 strides.begin(),
                                 0 );
        return data [ offset ] ;
    }

    /** \brief operator[] for positional addressing
     *
     * for an existing image In, this operator can
     * be used with multidimensional co-ordinates to
     * indicate position, so im[{x,y,z}] instead of
     * im.data [x + y * sizes [ 1 ] + z * sizes [ 1 ] * sizes [ 2 ] ] .
     *
     * This operator is for reading only.
     */
    inline value_type const operator[] ( std::initializer_list < size_t > const& indices ) const {
        size_t const offset =
            std::inner_product ( indices.begin(), indices.end(),
                                 strides.begin(),
                                 0 );
        return data [ offset ];
    }

    /** \brief operator[] for positional addressing
     *
     * This operator is for modifying the data.
     */
    inline value_type& operator[] ( std::initializer_list < size_t > const& indices ) {
        size_t const offset =
            std::inner_product ( indices.begin(), indices.end(),
                                 strides.begin(),
                                 0 );
        return data [ offset ] ;
    }

    /** brief compute indices at offset
     *
     * inverse of positional addressing - given a position
     * in the 1D data vector, what are its multidimensional
     * indices?
     */
    inline const std::vector <size_t> indices_at_offset ( size_t pos1d ) {

        auto p = pos1d;
        std::vector <size_t> out ( sizes.size() );

        for ( size_t d = sizes.size()-1; d>0; d-- ) {
            out[d]  =      p / strides[d];
            p      -= out[d] * strides[d];
        }
        out[0] = p;

        return out;

    }

    /** brief compute if two neighbours are valid
     *
     * inverse of positional addressing - given a position
     * in the 1D data vector, what are its multidimensional
     * indices?
     */
    inline bool valid_neighbours ( long offset1, long offset2, long delta = 1 ) {

		long
			off1 = offset1, 
			off2 = offset2,
			pos1, pos2,
			pd = 0;
		bool valid = true;
			
        for ( long d = sizes.size()-1; d>=0; d-- ) {
			
            pos1  = off1 / strides[d];
            pos2  = off2 / strides[d];		
			off1 -= pos1 * strides[d];
			
			if (  ( ( pos1 - pos2 ) >  delta ) ||
				  ( ( pos1 - pos2 ) < -delta ) ) {
					  
				valid = false;
				break;
				
			} else {

				off2  -= pos2 * strides[d];
				
				if ( pos1 != pos2 ) 
					pd++;
				
				if ( pd > delta ) {
					valid = false;
					break;
				}	
			
			} // if pos
					
		} // for d

        return ( valid );

    }
	


    /** \brief getsize() - returns the dimensions
     *
     * return dimensions as a std::vector <size_t>
     */
    std::vector<size_t> getsize (          ) {
        return sizes;
    }

    /** \brief getsize() - returns one dimension
     *
     * return dimension given its index in sizes
     */
    size_t              getsize ( size_t s ) {
        return sizes[s];
    }

    /** \brief getdatasize() - returns the number of intensities
     *
     * the size of the vector 'data'
     */
    size_t              getdatasize ( ) {
        return data.size();
    }

    /** \brief vector_set() - sets the data vector
     *
	 * The input vector should be the same size as the data vector
	 * 
     */
	template <typename T>
    void vector_set ( std::vector <T> newdata ) {
		
		if ( newdata.size() == data.size() )
			
			std::copy ( newdata.begin(), newdata.end(), data.begin() );		
			
    }

    /** \brief vector_mask() - sets the data vector to 0 outside where a mask vector is 0
     *
	 * The mask vector should be the same size as the data vector
	 * 
     */
	template <typename T>
    void vector_mask ( std::vector <T> mask ) {
		
		if ( mask.size() == data.size() )
			
			std::transform ( data.begin(), data.end(), mask.begin(), data.begin(),
								[] ( auto p1, auto p2 ) { return ( ( ! p2 ) ? 0 : p1 ); } );
		
    }

    /** \brief vector_import() - sets the data vector to 0 outside where a mask vector is 0
     *
	 * The mask vector should be the same size as the data vector
	 * 
     */
    void vector_import ( std::vector <value_type> input ) {
		
		if ( input.size() == data.size() )
			std::copy ( input.begin(), input.end(), data.begin() );		
			
    }
	
    /** \brief vector_export() - returns the data vector
     *
	 * The mask vector should be the same size as the data vector
	 * 
     */
    std::vector<value_type> vector_export () {

		std::vector<value_type> 
			output ( data );

		return output;
		
    }

    /** \brief getdata_ptr() - returns the address of the data vector
     *
     * this is a pointer to a vector -- use with care
     */
    std::vector <value_type>*          getdata_ptr ( ) {
        return &data;
    }

    /** \brief reshape() - change dimensions
     *
     * Changes the sizes and strides vectors.
     * Only works if total #elements does not change.
     */
    void reshape ( std::vector < size_t > const& newsizes ) {

        if ( std::accumulate( newsizes.begin(), newsizes.end(), 1, std::multiplies<size_t>() ) ==
             std::accumulate(    sizes.begin(),    sizes.end(), 1, std::multiplies<size_t>() )
           ) {

            sizes = newsizes;
            validate_sizes();

        } else

            std::cout << "reshape impossible because it would change image size";

    } // reshape 



    /** \brief operators += for scalar and bisimage, repectively
     *
     */
    const bisimage<value_type>& operator+= ( const value_type& rhs ) {
        assert ( std::is_arithmetic<value_type>::value );
		for ( auto& lhs: data )
			lhs += rhs;
        return (*this);
    }
	
    template <typename value_type_2>
    const bisimage<value_type>& operator+= ( const bisimage<value_type_2>& rhs ) {
        assert ( std::is_arithmetic<value_type>::value   );
        assert ( std::is_arithmetic<value_type_2>::value );
        assert ( data.size() == rhs.data.size()          );
		std::transform ( data.begin(), data.end(), rhs.data.begin(), data.begin(), std::plus<value_type> () );
        return (*this);
    }

    /** \brief operator + for templated types
     *
     * in this case, types for which += has been defined
     */
    template <typename value_type_2>
    const bisimage<value_type> operator+ ( const value_type_2& rhs ) {
        bisimage out(*this);
        out += rhs;
        return out;
    }

    /** \brief operators *= for scalar and bisimage, repectively
     *
     */
    const bisimage<value_type>& operator*= ( const value_type& rhs ) {
        assert ( std::is_arithmetic<value_type>::value );
		for ( auto& lhs: data )
			lhs *= rhs;
        return (*this);
    }
    template <typename value_type_2>
    const bisimage<value_type>& operator*= ( const bisimage<value_type_2>& rhs ) {
        assert ( std::is_arithmetic<value_type>::value   );
        assert ( std::is_arithmetic<value_type_2>::value );
        assert ( data.size() == rhs.data.size()          );
		std::transform ( data.begin(), data.end(), rhs.data.begin(), data.begin(), std::multiplies<value_type> () );
        return (*this);
    }

    /** \brief operator * for templated types
     *
     * in this case, types for which *= has been defined
     */
    template <typename value_type_2>
    const bisimage<value_type> operator* ( const value_type_2& rhs ) {
        assert ( std::is_arithmetic<value_type>::value );
        bisimage out(*this);
        out *= rhs;
        return out;
    }

    /** \brief operators -= for scalar and bisimage, repectively
     *
     */
    const bisimage<value_type>& operator-= ( const value_type rhs ) {
        assert ( std::is_arithmetic<value_type>::value );
		for ( auto& lhs: data )
			lhs -= rhs;
        return (*this);
    }
    template <typename value_type_2>
    const bisimage<value_type>& operator-= ( const bisimage<value_type_2>& rhs ) {
        assert ( std::is_arithmetic<value_type>::value   );
        assert ( std::is_arithmetic<value_type_2>::value );
        assert ( data.size() == rhs.data.size()          );
		std::transform ( data.begin(), data.end(), rhs.data.begin(), data.begin(), std::minus<value_type> () );
        return (*this);
    }

    /** \brief operator - for templated types
     *
     * in this case, types for which -= has been defined
     */
    template <typename value_type_2>
    const bisimage<value_type> operator- ( const value_type_2& rhs ) {
        bisimage out(*this);
        out -= rhs;
        return out;
    }

    /* \brief operator - without operand
     *
     * negate yourself
     */
    const bisimage<value_type> operator- ( void ) {
        bisimage out(*this);
        out *= -1;
        return out;
    }

    /** \brief operators /= for scalar and bisimage, repectively
     *
     */
    const bisimage<value_type>& operator/= ( const value_type rhs ) {
        assert ( std::is_arithmetic<value_type>::value );
		for ( auto& lhs: data )
			lhs /= rhs;
        return (*this);
    }
    template <typename value_type_2>
    const bisimage<value_type>& operator/= ( const bisimage<value_type_2>& rhs ) {
        assert ( std::is_arithmetic<value_type>::value   );
        assert ( std::is_arithmetic<value_type_2>::value );
        assert ( data.size() == rhs.data.size()          );
		std::transform ( data.begin(), data.end(), rhs.data.begin(), data.begin(), std::divides<value_type> () );
        return (*this);
    }

    /** \brief operator / for templated types
     *
     * in this case, types for which /= has been defined
     */
    template <typename value_type_2>
    const bisimage<value_type> operator/ ( const value_type_2& rhs ) {
        bisimage<value_type> out(*this);
        out /= rhs;
        return out;
    }

    /** \brief operator ^, raise all elements to a power
     *
     */
    template <typename value_type_2>
    const bisimage<value_type> operator^ ( const value_type_2& rhs ) {
        bisimage<value_type> out (*this);
        std::transform ( data.begin(), data.end(), data.begin(), 
                         [ &rhs ] ( auto in ) { return ( static_cast <value_type> ( std::pow ( static_cast<float> ( in ), rhs ) ) ); } );
        return out;
    }

    /** \brief reciprocal function, returns 1/c
     *
     * Return an image with 1/c for all coefficients c.
     * This function is used for dividing by an image.
     */
    const bisimage<value_type>& reciprocal () {
        bisimage<value_type> out(*this);
        std::transform ( data.begin(), data.end(), data.begin(), 
                         [] ( auto in ) { return ( 1. / in ); } );
        return out;
    }
    
    /** \brief sqrt function, returns square root of each element
     *
     */
    const bisimage<value_type>& sqrt () {
        bisimage<value_type> out(*this);
        return ( out ^ .5 );
    }
    
    /** \brief sqr function, returns square of each element
     *
     */
    const bisimage<value_type>& sqr () {
        bisimage<value_type> out(*this);
        return ( out ^ 2. );
    }

    

/** Functions outside of the bisimage class.
 *  These can be used by bisimage objects, but
 *  are not required by the class and do not
 *  require it themselves.
 *
 */



// Reducers are used to walk through the data as if they were 3D,
// squashing the dimensions before the 'at' dimension together
// (inner) and the ones after the 'at' dimension too (outer).
void set_reducers ( size_t dimension, 
                    reducer& dims, 
                    reducer& limits,
                    reducer& steps ) {
    
    auto 
        at_dim  = dimension,
        oldsize = getsize();
    
    // for first dimension, insert a 1 before
    if ( ! at_dim ) {
        oldsize.insert ( oldsize.begin(), 1 );
        at_dim++;                            
    }                                        
        
    // for last dimension, insert a 1 after
    if ( at_dim == ( oldsize.size() - 1 ) )
        oldsize.push_back ( 1 );
    
    // compute the dimensions during the 3D stage
    dims.atr = oldsize [ at_dim ];
    for ( size_t d = 0; d < at_dim; d++ ) 
        dims.inr *= oldsize [ d ];             // inr: product of dimensions up to dim
    for ( size_t d = at_dim + 1; d < oldsize.size(); d++ ) 
        dims.out *= oldsize [ d ];        // outer: product of dimensions dim+1 and up

    // compute the steps (the 'inr' step is always 1)
    steps.inr = 1;
    steps.atr    = dims.inr;
    steps.out = dims.atr * dims.inr; 

    // only the limit of the outer strides is fixed: data size
    limits.out = strides.back();
                        
}



/** \brief sum -- sum the image [ along a dimension ]
 * 
 * parameters:
 *      size_t dimension: the dimension that is reduced to its sum
 *                        default: 0
 * 
 * the sum of the whole image (no parameters, output one value) is a one-liner
 * 
 */
 
const value_type 
    sum () { return ( std::accumulate ( data.begin(), data.end(), 0 ) ); }

const bisimage < value_type > 
    sum ( size_t dimension = 0 ) {

        auto             // new size is 1 dimension smaller
            newsize = sizes;
        newsize.erase ( newsize.begin() + dimension );    
                   
        reducer          // Reducers are used to walk through the data as is it were 3D,
            dims, lims, offs, steps;
        
        // get the values for this particular combination
        set_reducers ( dimension, dims, lims, steps );
        
        bisimage < value_type > 
            output ( { dims.inr, dims.out } );
        size_t 
            index_2d;
  
        for ( offs.out = 0, 
			  index_2d = 0;   
			   offs.out < lims.out;   
				offs.out += steps.out )
			for ( offs.inr = offs.out, 
				  lims.inr = offs.out + dims.inr;   
				   offs.inr < lims.inr;
					offs.inr ++, 
					index_2d ++ ) 
				for ( offs.atr = offs.inr, 
					  lims.atr = offs.inr + dims.atr * steps.atr;   
					   offs.atr < lims.atr;   
					    offs.atr += steps.atr )                    
					output.data [ index_2d ] += data [ offs.atr ]; // formula for summation

        output.reshape ( newsize );
        return ( output );
            
    }



/** \brief mean -- mean of the image [ along a dimension ]
 * 
 * parameters:
 *      size_t dimension: the dimension along whhich variance is computed
 *                        default: 0
 *      bool sample: use sample variance (normalise ny n-1 not n)
 *                   default: false
 * 
 * the mean of the whole image (no parameters, output one value) is a one-liner
 * 
 */
 
const value_type 
    mean () { return ( sum() / getdatasize() ); }

const bisimage < value_type > 
    mean ( size_t dimension = 0 ) {
        auto 
            output = sum ( dimension ) / sizes [ dimension ];
        return ( output );           
    }

/** \brief var -- variance in the image along a dimension
 * 
 * parameters:
 *      size_t dimension: the dimension along whhich variance is computed
 *                        default: 0
 *      bool sample: use sample variance (normalise ny n-1 not n)
 *                   default: false
 * 
 * the variance of the whole image (only correction parameter) 
 *                  is simpler but bit of bookkeeping required
 * 
 */

const value_type
	var( std::string group = "population" ) {
		
		float_t 
			bessel_correction = ( group == "sample" ) ? 1 : getdatasize() / ( getdatasize() - 1 );
		value_type
			mymeansq = std::inner_product ( data.begin(), data.end(), data.begin(), 0 ),
			mymean   = std::accumulate    ( data.begin(), data.end(),               0 ) / getdatasize();
		return 
			( bessel_correction * mymean * mymean / mymeansq );
		
	}
		
const bisimage < value_type > 
    var ( size_t dimension = 0, std::string group = "population" ) {		

		float_t
			bessel_correction = ( group == "sample" ) ? 1 : sizes [ dimension ] / ( sizes [ dimension -1 ] );
		auto 
			mymeansq   = (*this) ^ 2,
			mymeand    = mean ( dimension ) * bessel_correction;
		return
			( mymeand.sqr() - mymeansq.mean ( dimension ) );
            
    }
    


/** \brief get3Dline() returns a 3D line (along a dimension) from a volume
 *
 * Result is a std::vector of value_type
 *
 * The line is sampled along dimension dim ( 0, 1 or 2 respectively)
 * at position pos1, pos2 in the other dimensions
 * ( 1 and 2, 0 and 2, or 0 and 1, respectively).
 */
const std::vector<value_type> get3Dline ( size_t dim, size_t pos1, size_t pos2,
                                              size_t linestart, size_t lineend ) {

    std::vector <value_type> out;

    if (sizes.size() != 3) {

        std::cout << "3D lines must be selected from a 3D image\n";

    } else {

        size_t
        step    = strides[dim],
        slicex  = (dim>0) ? 0 : 1,
        slicey  = (dim>1) ? 1 : 2,
        line0   = std::max <size_t> ( linestart, 0          ),
        line1   = std::min <size_t> ( lineend,   sizes[dim] );

        value_type*
        dptr    = data.data() + pos1 * strides[slicex] + pos2 * strides[slicey];

        out.resize ( line1-line0 );

        for ( size_t i = line0; i < line1; i++, dptr+=step )
            out[i] = *dptr;

    } // if sizes

    return out;

} // get3Dline



/** \brief set3Dline() puts a 3D line (along a dimension) in a volume
 *
 * Input is a std::vector of value_type
 *
 * The line is inserted along dimension dim ( 0, 1 or 2 respectively)
 * at position pos1, pos2 in the other dimensions
 * ( 1 and 2, 0 and 2, or 0 and 1, respectively).
 */
void set3Dline ( std::vector <value_type>& in,
                              size_t dim, size_t pos1, size_t pos2,
                              size_t linestart, size_t lineend ) {

    if (sizes.size() != 3) {

        std::cout << "3D lines must be selected from a 3D image\n";

    } else {

        size_t
        step = strides[dim],
        slicex = (dim>0) ? 0 : 1,
        slicey = (dim>1) ? 1 : 2,
        line0 = std::max <size_t> ( linestart, 0          ),
        line1 = std::min <size_t> ( line0 + lineend, sizes[dim] );

        value_type* dptr = data.data() + pos1 * strides[slicex] + pos2 * strides[slicey];

        for ( size_t i = line0; i < line1; i++, dptr+=step )
            *dptr = in[i];

    } // if sizes

} // set3Dline



/** \brief getSlice() returns a 2D slice from a volume
 *
 * Result is an bisimage of value_type
 * which is a copy of slice no. <sli>
 * along dimension <dim> (0, 1, or 2)
 * and at position <sli>.
 */
const bisimage<value_type> getSlice( size_t dim, size_t sli, std::string filename ) {
    // get a slice from a volume
    // optional: write it out as a .bmp file
    //

    bisimage<value_type> out;

    if (sizes.size() != 3) {

        std::cout << "slices can only be selected from a 3D image\n";

    } else {

        // slice sizes are called slicex (lowest 3D dim index) and slicey (highest)
        size_t slicex = (dim>0) ? 0 : 1;
        size_t slicey = (dim>1) ? 1 : 2;

        // set sizes for sizes, strides and data start with 3D
        out.sizes     = {    sizes[slicex], sizes[slicey], 1                                             };
        out.strides   = { 1, sizes[slicex], sizes[slicex] * sizes[slicey], sizes[slicex] * sizes[slicey] };
        out.data.resize (    *out.strides.rbegin()                                                       );

        // fill the slice by calling get3Dline (from folume) for each y line in the slice
        // loop over highest (outer == slower with largest strides) dimension first
        //
        // dim x -> yz slice, slicex 1, slicey 2 -> lines in y, loop over z -> line pos [ sli z ] = [  sli ypos ]
        // dim y -> xz slice, slicex 0, slicey 2 -> lines in x, loop over z -> line pos [ sli z ] = [  sli ypos ]
        // dim z -> xy slice, slicex 0, slicey 1 -> lines in x, loop over y -> line pos [ y sli ] = [ ypos  sli ]
        for ( size_t ypos=0; ypos<sizes[slicey]; ypos++ ) {

            // position where the line is taken:
            //
            // an x line is taken from an y,z position
            size_t linx = ( slicey>1 ) ?  sli : ypos;
            size_t liny = ( slicey>1 ) ? ypos :  sli;

            std::vector<value_type> sli_line = get3Dline( slicex, linx, liny );
            out.set3Dline ( sli_line, 0, ypos, 0); // x line (0), put in y position liny, 'z' position 0

        } // for ypos

    } // if sizes

    if ( !filename.empty() ) {

        cimg_library::CImg<value_type>*
        my_bitmap = new cimg_library::CImg<value_type> (out.data.data(),
                out.sizes[0],
                out.sizes[1],
                1, 1, true);
        my_bitmap->rotate(180);
        my_bitmap->save_bmp( filename.c_str() );
        delete ( my_bitmap );

    } // if filename

    out.reshape( { out.sizes[0], out.sizes[1] } ); // remove dimension 3 (which is 1) of the output slice
    return out;

} // getSlice

/** \brief setSlice() insert a 2D slice into a volume
 *
 * Slice no. <sli> along
 * dimension <dim> (0, 1, or 2)
 * of the current object is copied from the input
 */
void setSlice ( bisimage<value_type>& input, size_t dim, size_t sli ) {

    if (sizes.size() != 3) {

        std::cout << "slices can only be inserted into a 3D image\n";

    } else {

        // slice sizes are called slicex (lowest 3D dim index) and slicey (highest)
        size_t slicex = (dim>0) ? 0 : 1;
        size_t slicey = (dim>1) ? 1 : 2;

        // check if input sizes match the current imagNd dimensions
        if ( ( input.sizes[0] != sizes[slicex] ) | ( input.sizes[1] != sizes[slicey] ) ) {

            std::cout << "input slice deimensions do not match volume slice size \n";

        }

        // briefly make our slice 3D for using get3Dline()
        input.reshape ( input.sizes[0], input.size[1], 1 );

        for ( size_t ypos=0; ypos<sizes[slicey]; ypos++ ) {

            size_t linx = ( slicey>1 ) ?  sli : ypos;
            size_t liny = ( slicey>1 ) ? ypos :  sli;

            std::vector<value_type> sli_line = input.get3Dline ( sli_line, 0, ypos, 0);
            set3Dline( slicex, linx, liny );

        } // for ypos

        // make our slice 2D again
        input.reshape ( input.sizes[0], input.size[1] );

    } // if sizes

} // getSlice



/** \brief getsubvolume() returns a 3D subvolume from a 3D image
 *
 * Ranges are given as xmin, xmax, ymin, ymax, zmin, zmax
 * Input and output are both of type bisimage, and this
 * routine works exclusively with 3D data
 */
const bisimage<value_type> getsubvolume ( size_t startx, size_t endx,
                                              size_t starty, size_t endy,
                                              size_t startz, size_t endz) {

    bisimage <value_type> out;

    if (sizes.size() != 3) {

        std::cout << "subvolumes can only be selected from a 3D image\n";

    } else {

        size_t  x0 = std::max<size_t> ( startx, 0 ),
                y0 = std::max<size_t> ( starty, 0 ),
                z0 = std::max<size_t> ( startz, 0 ),
                x1 = std::min<size_t> ( endx, sizes[0] ),
                y1 = std::min<size_t> ( endy, sizes[1] ),
                z1 = std::min<size_t> ( endz, sizes[2] );

        out.sizes   = { std::max<size_t> (x1 - x0, 1),
                        std::max<size_t> (y1 - y0, 1),
                        std::max<size_t> (z1 - z0, 1)
                      };
        out.strides = { 1, out.sizes[0], out.sizes[1] * out.sizes[0], out.sizes[2] * out.sizes [1] * out.sizes[0] };
        out.data.resize( *out.strides.rbegin() );

        value_type *dptr = out.data.data();

        for ( size_t z=z0; z<z1; z++ )
            for ( size_t y=y0; y<y1; y++ )
                for ( size_t x=x0; x<x1; x++ )
                    *dptr++ = operator[] ( { x, y, z } );

    } // if sizes

    return out;

} // getsubvolume

/** \brief setsubvolume() inserts a 3D subvolume into a 3D image
 *
 * Ranges are given as xmin, ymin, zmin for where to insert
 * Source and destination are both of type bisimage, and this
 * routine works exclusively with 3D data
 */
void setsubvolume ( bisimage <value_type>& in,
                                 size_t startx,
                                 size_t starty,
                                 size_t startz ) {

    if ( (sizes.size() != 3) | (in.sizes.size() != 3) ) {

        std::cout << "only 3D can be put in only 3D images\n";

    } else {

        size_t  x0 = std::max<size_t> ( startx, 0 ),
                y0 = std::max<size_t> ( starty, 0 ),
                z0 = std::max<size_t> ( startz, 0 ),
                x1 = std::min<size_t> ( startx + in.sizes[0], sizes[0] ),
                y1 = std::min<size_t> ( starty + in.sizes[1], sizes[1] ),
                z1 = std::min<size_t> ( startz + in.sizes[2], sizes[2] );

        value_type *dptr = &in.data[0];

        for ( size_t z=z0; z<z1; z++ )
            for ( size_t y=y0; y<y1; y++ )
                for ( size_t x=x0; x<x1; x++ )
                    operator[] ({ x, y, z }) = *dptr++;

    } // if sizes

} // getsubvolume



/** \brief filter() filter along dimension {0, 1 or 2} in a 3D image
 *
 * The filter is given as a numrical std::vector
 * this method uses the function filterline
 * (outside this class)
 */
void filter ( std::vector<double> filt, size_t dim ) {

    if (sizes.size() != 3) {

        std::cout << "currently filter works only for 3D images\n";

    } else {

        size_t slicex = (dim>0) ? 0 : 1;
        size_t slicey = (dim>1) ? 1 : 2;

        for ( size_t posx = 0; posx < sizes[slicex]; posx++ )

            for ( size_t posy = 0; posy < sizes[slicey]; posy++ ) {

                std::vector <value_type> sign = get3Dline( dim, posx, posy );
                filterline ( sign, filt );
                set3Dline( sign, dim, posx, posy);

            } // for posy


    } // if sizes

} //filter

/** \brief fwt_slice() apply the fwt to the 2 dimensions of each slice
 *
 * Slice direction is taken as the 1st dimension (this is not standard!).
 * Dimensions 0 and 1 in each slice are dimension 1 and 2 of the input.
 * This methos uses fwtline (outside this class), wavelet basis is given
 * by its name ("Daubecchies4", "Haar").
 */
void fwt_slices ( std::string wname, size_t level ) {

    // get filter as a vector not a name
    std::vector <double> wvec = owfilter( wname );

    // set sizes
    size_t xmax = sizes[0];
    size_t ymax = sizes[1];
    size_t zmax = sizes[2];

    // loop over levels
    for ( size_t l=0; l<level; l++ ) {

        // loop over x slices, do fwt in z then y
        for ( size_t x=0; x<xmax; x++ ) {

            for ( size_t y=0; y<ymax; y++ ) {
                std::vector <value_type> dwt1 = get3Dline( 2, x, y, 0, zmax );
                fwtline( dwt1, wvec, 1 );
                set3Dline( dwt1, 2, x, y, 0, zmax );
            } // for y

            for ( size_t z=0; z<zmax; z++ ) {
                std::vector <value_type> dwt1 = get3Dline( 1, x, z, 0, ymax );
                fwtline( dwt1, wvec, 1 );
                set3Dline( dwt1, 1, x, z, 0, ymax );
            } // for y

        } // for x

        ymax /= 2;
        zmax /= 2;

    } // for l

} // fwtslices

/** \brief ifwt_slice() apply the fwt to the 2 dimensions of each slice
 *
 * Slice direction is taken as the 1st dimension (this is not standard!).
 * Dimensions 0 and 1 in each slice are dimension 1 and 2 of the input.
 * This methos uses ifwtline (outside this class), wavelet basis is given
 * by its name ("Daubecchies4", "Haar").
 */
void ifwt_slices ( std::string wname, size_t level ) {

    // get filter as a vector not a name
    std::vector <double> wvec = owfilter( wname );

    size_t fac=1;
    for (size_t i=1; i<level; i++)
        fac *= 2;

    // set sizes
    size_t xmax = sizes[0];
    size_t ymax = sizes[1]/fac;
    size_t zmax = sizes[2]/fac;

    // loop over levels
    for ( size_t l=0; l<level; l++ ) {

        // loop over x slices, do fwt in z then y
        for ( size_t x=0; x<xmax; x++ ) {

            for ( size_t y=0; y<ymax; y++ ) {
                std::vector <value_type> dwt1 = get3Dline( 2, x, y, 0, zmax );
                ifwtline( dwt1, wvec, 1 );
                set3Dline( dwt1, 2, x, y, 0, zmax );
            } // for y

            for ( size_t z=0; z<zmax; z++ ) {
                std::vector <value_type> dwt1 = get3Dline( 1, x, z, 0, ymax );
                ifwtline( dwt1, wvec, 1 );
                set3Dline( dwt1, 1, x, z, 0, ymax );
            } // for y

        } // for x

        ymax *= 2;
        zmax *= 2;

    } // for l

} // ifwtslices



/** \brief waveletVisuThresh() apply hard VisuThresh shrinkage
 *         to 2D wavelet channels of each sllice in a 3D image
 *
 * This method is only useful on images after applying
 * the fwt_slices() method, see above, at the same level
 */
void waveletVisuThresh ( size_t level, value_type booster ) {

    // for 2D wavelet channels in slices, apply the VisuThresh
    // denoising threshold to the detail coefficients. This is
    // defined as sqrt(2 * log(|c|)) wher c is the size of the
    // wavelet channel. These channels are the same size in all
    // slices, so we can select them together.

    size_t xmax = sizes[0];
    size_t ymax = sizes[1];
    size_t zmax = sizes[2];

    for ( size_t l = 0; l< level; l++ ) {

        ymax /= 2;
        zmax /= 2;

        double threshold = booster * sqrt ( 2 * log (zmax) / log (2) );

        {
            // y only
            bisimage <value_type> workspace = getsubvolume(  0, xmax,
                                             ymax, 2*ymax,
                                             0, zmax );
            for ( size_t i=0; i<workspace.data.size(); i++)
                workspace.data[i] *= (workspace.data[i]>threshold);
            setsubvolume( workspace, 0, ymax, 0 );
        }

        {
            // z only
            bisimage <value_type> workspace = getsubvolume(  0, xmax,
                                             0, ymax,
                                             zmax, 2*zmax );
            for ( size_t i=0; i<workspace.data.size(); i++)
                workspace.data[i] *= (workspace.data[i]>threshold);
            setsubvolume( workspace, 0, 0, zmax );
        }

        {
            // both y and z
            bisimage <value_type> workspace = getsubvolume(  0, xmax,
                                             ymax, 2*ymax,
                                             zmax, 2*zmax );
            for ( size_t i=0; i<workspace.data.size(); i++)
                workspace.data[i] *= (workspace.data[i]>threshold);
            setsubvolume( workspace, 0, ymax, zmax );
        }

    } // for l

} // waveletVisuThresh



int GetNNeigbours(int ip, int* NeighbourIndices, int ndim, size_t* dimensions) {
    if(ndim<=0 || ndim>3) {
        std::cout<<"ERROR: MImage::GetNNeigbours(). ndim (="<< ndim <<") out of range. \n";
        return 0;
    }
    if(NeighbourIndices==NULL) {
        std::cout<<"ERROR: MImage::GetNNeigbours(). Invalid NULL argument. \n";
        return 0;
    }

    int dimx  = dimensions[0];
    int dimy  = dimensions[1];
    int dimz  = dimensions[2];
    int dimxy = dimx*dimy;

// Test for out of range
    int ix = ndim>0 ? (  ip    % dimx       ) : 0;
    int iy = ndim>1 ? (((ip-ix)/ dimx)%dimy ) : 0;
    int iz = ndim>2 ? (  ip    / dimxy      ) : 0;

    if((ndim>0 && (ix<0 || ix>=dimx))  ||
            (ndim>1 && (iy<0 || iy>=dimy))  ||
            (ndim>2 && (iz<0 || iz>=dimz))) {
        std::cout<<"ERROR: MImage::GetNNeigbours(). point index out of range (ix, iy, iz) = ("<< ix <<", "<< iy <<" "<< iz << ")\n";
        return 0;
    }

    int NNeig = 0;
    if(ndim>0 && dimx>1) {
        if(ix>0     )
            NeighbourIndices[NNeig++]=ip-1;
        if(ix<dimx-1)
            NeighbourIndices[NNeig++]=ip+1;
    }
    if(ndim>1 && dimy>1) {
        if(iy>0     )
            NeighbourIndices[NNeig++]=ip-dimx;
        if(iy<dimy-1)
            NeighbourIndices[NNeig++]=ip+dimx;
    }
    if(ndim>2 && dimz>1) {
        if(iz>0     )
            NeighbourIndices[NNeig++]=ip-dimxy;
        if(iz<dimz-1)
            NeighbourIndices[NNeig++]=ip+dimxy;
    }
    return NNeig;
}



bool GetWatershedImage() {

    if ( sizes.size() > 3 ) {
        std::cout<<"ERROR: MImage::GetWatershedImage(). Invalid dimensionality ("<< sizes.size() <<"). \n";
        return false;
    }

    value_type max_value = *std::max_element(data.begin(),data.end());
    value_type min_value = *std::min_element(data.begin(),data.end());
    size_t ndim=sizes.size();
    size_t NP=getdatasize();

    std::vector<int>
    Index   (NP),
            Dist    (NP),
            Label   (NP),
            Hist    (max_value + 2),
            CHist   (max_value + 2),
            NeigArr (200);

    // check if image needs inverting and do so if yes
    // (if pixel [0] has lower value than maximum/2 ?)
    if ( data[0] < (max_value/2) ) {
        std::cout << "inverting ... \n";
        for ( size_t i=0; i< NP; i++ )
            data [i] = max_value - min_value - data[i];
    }

    // build the histogram
    for (unsigned n=0; n < NP; n++)
        Hist[data[n]]++;

    // build the cumulative histogram (differs from histogram after index 0)
    for (unsigned k=1; k < max_value+1; k++)
        CHist[k] = CHist[k-1] + Hist[k-1];

    // label point based on value in cumulative histogram -- increasing index to number within intensity
    for (unsigned n=0; n < NP; n++)
        Index[CHist[data[n]]++] = n;

    // subtract histogram from cumulative after labelling
    for (unsigned k=0; k< max_value+1; k++)
        CHist[k] -= Hist[k]; // restore cumulative histogram

    CHist[max_value+1] = NP; // this was still 0

    const int LABELINIT  =   -1;
    const int MASK       =   -2;
    const int WSHED      =    0;
    const int FICTITIOUS =   -3;

    // initialise labels
    for ( unsigned n=0; n< NP; n++)
        Label[n] = LABELINIT;

    std::queue<int> fifoQueue;
    int curlab = 0;

    // Geodesic SKIZ of level h-1 inside level h. INCLUDE LAST LEVEL!
    for( value_type h = min_value; h<=max_value; h++) {
        for( int pixelIndex = CHist[h]; pixelIndex < CHist[h+1]; pixelIndex++) { //mask all pixels at level h
            int   ip  = Index[pixelIndex];
            Label[ip] = MASK;

            int NNEig = GetNNeigbours(ip, NeigArr.data(), ndim, sizes.data());

            for(int i=0; i<NNEig; i++) {
                if(Label[NeigArr[i]] < 0 && Label[NeigArr[i]] != WSHED)
                    continue;

                Dist[ip] = 1;  //Initialise queue with neighbours at level h of current basins or watersheds
                fifoQueue.push(ip);
                break;
            }
        }

        int curdist = 1;
        fifoQueue.push(FICTITIOUS);

        while(true) { // extend basins
            int voxelIndex = fifoQueue.front();
            fifoQueue.pop();

            if(voxelIndex == FICTITIOUS) {
                if(fifoQueue.empty())
                    break;

                fifoQueue.push(FICTITIOUS);
                curdist++;
                voxelIndex = fifoQueue.front();
                fifoQueue.pop();
            }

            int NNEig = GetNNeigbours(voxelIndex, NeigArr.data(), ndim, sizes.data());
            for(int i=0; i<NNEig; i++) { // Labelling p by inspecting neighbours
                if(Dist[NeigArr[i]] < curdist && (Label[NeigArr[i]] > 0 || Label[NeigArr[i]]==WSHED)) {
                    if(Label[NeigArr[i]] > 0) { // q belongs to an existing basin or to a watershed
                        if(Label[voxelIndex] == MASK || Label[voxelIndex] ==WSHED)
                            Label[voxelIndex] = Label[NeigArr[i]]; // Removed from original algorithm || p.isLabelWSHED() )
                        else if(Label[voxelIndex] != Label[NeigArr[i]])
                            Label[voxelIndex] = WSHED;

                    } // end if lab>0
                    else if (Label[voxelIndex]==MASK)
                        Label[voxelIndex] = WSHED;
                } else if(Label[NeigArr[i]]==MASK && Dist[NeigArr[i]]==0) {
                    Dist[NeigArr[i]] = curdist + 1;   //q is plateau pixel
                    fifoQueue.push(NeigArr[i]);
                }
            } // end for, end processing neighbours
        } // end while (loop)

        // Detect and process new minima at level h
        for(int pixelIndex = CHist[h]; pixelIndex < CHist[h+1]; pixelIndex++) { //mask all pixels at level h
            int ip   = Index[pixelIndex];
            Dist[ip] = 0;       // Reset distance to zero

            if(Label[ip]!=MASK)
                continue;
            curlab++;       // The pixel is inside a new minimum , create new label
            fifoQueue.push(ip);
            Label[ip] = curlab;

            while(fifoQueue.size()) {
                int voxelIndex = fifoQueue.front();
                fifoQueue.pop();

                int NNEig = GetNNeigbours(voxelIndex, NeigArr.data(), ndim, sizes.data());  // replaced ip by voxelIndex
                for(int i=0; i<NNEig; i++) { // inspect neighbours of q
                    if(Label[NeigArr[i]]!=MASK)
                        continue;

                    fifoQueue.push(NeigArr[i]);
                    Label[NeigArr[i]] = curlab;
                }
            } // end while
        } // end for
    } // loop over h

    int MINS = (1<<15) -1;
    for ( unsigned i=0; i<NP; i++)
        data[i] = short(fMinOf2<int>(MINS, Label[i]));

    return true;
}  



/*
 * make std::cout so that it prints respecting the dimensions
 * this requires the friend declaration for std::cout
 * 
 */
friend std::ostream& operator<< ( std::ostream& sout, const bisimage output) {
    sout << std::setw ( 6 ) << std::setfill ( ' ' ) << "{";          
    switch ( output.sizes.size() ) {
        case (4): 
        for ( size_t t = 0; t < output.sizes [ 3 ]; t++ ) {
            sout << std::endl << "     ";
            for ( size_t y = 0; y < output.sizes [ 1 ]; y++ ) {
                for ( size_t z = 0; z < output.sizes [ 2 ]; z++ ) {
                    for ( size_t x = 0; x < output.sizes [ 0 ]; x++ )
                        sout << output ( { x, y, z, t } ) << " ";
                    sout << "\t";
                }
                sout << std::endl << "     ";
            }
        }
        break;
        case (3): 
        for ( size_t z = 0; z < output.sizes [ 2 ]; z++ ) {
            sout << std::endl << "     ";
            for ( size_t y = 0; y < output.sizes [ 1 ]; y++ ) {
                for ( size_t x = 0; x < output.sizes [ 0 ]; x++ )
                    sout << output ( { x, y, z } ) << " ";
                sout << std::endl << "     ";
            }
        }
        break;
        case (2): 
        for ( size_t y = 0; y < output.sizes [ 1 ]; y++ ) {
            sout << std::endl << "     ";
            for ( size_t x = 0; x < output.sizes [ 0 ]; x++ )
                sout << output ( { x, y } ) << " ";
        }
        break;
        default: 
        for ( auto elem: output.data )
            sout << elem << " ";
        break;
    }
    sout << "}" << std::endl;
    return ( sout );
}

 
 
/** \brief addNormalNoise() add normally distributed noise
 *
 * Uses standard library functions to produce random numbers
 * Parameters mu and sigma are doubles, used in the expected way
 */
void addNormalNoise ( double mu, double sigma ) {

    // random device class instance, source of 'true' randomness for initializing random seed
    std::random_device randomdevice{};

    // Mersenne twister PRNG, initialized with seed from previous random device instance
    std::mt19937 engine { randomdevice() };

    // set up the distribution
    std::normal_distribution <double> normaldist ( mu, sigma );

    // add noise to the data
    // N (mu, sigma) all with <double> values
    // (non-float leads to unexpected behaviour)
    for ( size_t i = 0; i < data.size(); i++ )
        data[i] += normaldist ( engine );

}

/** \brief addRicianNoise() add Rician distributed noise
 *
 * Uses standard library functions to produce random numbers
 * Parameters mu and sigma are doubles, used in the expected way
 * for two normal distributions.
 */
void addRicianNoise ( double mu, double sigma ) {

    // random device class instance, source of 'true' randomness for initializing random seed
    std::random_device randomdevice{};

    // Mersenne twister PRNG, initialized with seed from previous random device instance
    std::mt19937 engine { randomdevice() };

    // set up the distribution
    std::normal_distribution <double> normaldist ( mu, sigma );

    // add Rician noise to the data using 2 normally distributed noise values
    for ( size_t i = 0; i < data.size(); i++ ) {

        double n1 = data[i] + normaldist ( engine );
        double n2           = normaldist ( engine );
        data [i] = sqrt ( ( n1 * n1 ) + ( n2 * n2 ) );

    }

}


 
}; // class



/** \brief overloaded operators +, *, - and / for basic numerical types
  *
  * This makes it possible to not only do bisimage +  <type>
  * but also                               <type> + bisimage, etc.
  */
template <typename T, typename U>
inline bisimage <T> operator+ ( U x, bisimage <T> y) {
    return y             + x;
}
template <typename T, typename U>
inline bisimage <T> operator* ( U x, bisimage <T> y) {
    return y             * x;
}
template <typename T, typename U>
inline bisimage <T> operator- ( U x, bisimage <T> y) {
    return -y            + x;
}
template <typename T, typename U>
inline bisimage <T> operator/ ( U x, bisimage <T> y) {
    return reciprocal(y) * x;
}



}; // namespace



#endif // BISIMAGE_HPP_INCLUDED
