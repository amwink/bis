#ifndef BISMAXTREE_HPP_INCLUDED
#define BISMAXTREE_HPP_INCLUDED

/** \brief support parallel sorting
 *
 */
#include <execution>

/** \brief bismaxtree is a subclass of bisimage
 */
#include "bisimage.hpp"

/** \brief more goodies from the standard library
 *
 */
#include <queue>
#include <stack>
#include <type_traits>

/** \brief bis: namespace for multidimensional images
 *
 *  This namespace contains the bisimage class for
 *  creating and processing multidimensional images.
 */
namespace bis {

template <typename value_type> class bismaxtree : public bisimage<value_type> {

	// self and superclass
	using self       = bismaxtree<value_type>;
	using superclass =   bisimage<value_type>;

	using superclass::data;
	using superclass::sizes;
	using superclass::strides;

	typedef int level_t;
	typedef struct component {
		level_t		value;  			// quantised value
		size_t 		uniq;				// number of points with exactly this value
		size_t		size;    			// number of points with at most this value
		size_t		root;    			// offset of 1st point		
		size_t		parent;				// parent component number
		std::vector <size_t> children;	// components on top of this
		std::vector <size_t> points;	// 1D points (only unique - with own value)
	} component;
	
protected:
  
    /** \brief the header, data, sizes and strides should
     *         be usable by subclasses
     * 
     * components: vector with the max-tree components, including the points in the tree structure
     *
     */
	std::vector <component>
		components;

	static constexpr size_t
		undefined = std::numeric_limits<size_t>::max();

public:

	/** \brief default constructor
	 *
	 * (leves members empty)
	 *
	 */
	bismaxtree() {
	}

	/** \brief destructor
	 *
	 * no allocations with new are used, so no 
	 * delete required - members go out of scope
	 * 
	 */
	~bismaxtree() {
	}

	/** \brief constructor from a bisimage
	 *
	 * This constructor takes an n-dimensional image
	 * and bulds its maxtree representation
	 */
	bismaxtree ( const bisimage<value_type>& rhs,
	             const level_t levels = UINT16_MAX,
	             const char connectivity = 6,
	             const std::string method = "Berger" )
		: bisimage<value_type> ( rhs ) {

		////////////////////////////////////////////////////////////////////////////////
		//
		// offsets in the image for connectivities
		// assuming strides in 2d / 3d images
		//			 6: only horizontal, vertical or sideways neighbours
		//			26: all horizontal, vertical sideways and combinations
		//

		std::vector<long long> neighbours;
		switch ( connectivity ) {
		case 4:
			assert ( sizes.size() >= 2 );
			neighbours = { -strides[1], -1, 1, strides[1] };
			break;
		case 8:
			assert ( sizes.size() >= 2 );
			neighbours = { -strides[1] - 1, -strides[1], -strides[1] + 1, 
									   - 1,							   1,
							strides[1] - 1,	 strides[1],  strides[1] + 1
			             };
			break;
		case 6:
			assert ( sizes.size() >= 3 );
			neighbours = { -strides[2], -strides[1], -1, 1, strides[1], strides[2] };
			break;
		case 26:
			assert ( sizes.size() >= 3 );
			neighbours = { -strides[2] - strides[1] - 1, - strides[2] - strides[1], -strides[2] - strides[1] + 1,
			               -strides[2] - 1,              - strides[2],							- strides[2] + 1, 
						   -strides[2] + strides[1] - 1, - strides[2] + strides[1], -strides[2] + strides[1] + 1, 
						                -strides[1] - 1,              - strides[1],              -strides[1] + 1, 
						                             -1,													   1, 
										 strides[1] - 1,			    strides[1], 			  strides[1] + 1, 
						    strides[2] - strides[1] - 1,   strides[2] - strides[1],	 strides[2] - strides[1] + 1, 
						    strides[2] - 1, 			   strides[2], 				 strides[2]				 + 1, 
						    strides[2] + strides[1] - 1,   strides[2] + strides[1],	 strides[2] + strides[1] + 1
			             };
			break;
		}

		////////////////////////////////////////////////////////////////////////////////
		//
		// sort sorted intensities in the data vector "sorted"
		// and also keep an array "indices" to see where they were in the image
		//

		std::vector<size_t> indices ( data.size() );

		std::iota ( indices.begin(), indices.end(), 0 );      // fill with 0, 1, 2, ..
		//if ( method == "Berger" )                           // to label last index as
		//   std::reverse ( indices.begin(), indices.end() ); // root (ICIP 2007 paper)
		std::stable_sort ( indices.begin(), indices.end(),    // 'data' sorted for determining 'indices'
						   [&] ( size_t i, size_t j ) { return ( data[i] < data[j] ); } );

		auto 
			mn = data[indices.front()], 
			mx = data[indices.back()];

		////////////////////////////////////////////////////////////////////////////////
		//
		// quantise the image, store in "quant"
		//

		std::vector <level_t> 
			quant ( data.size(), 0 );
		std::vector <size_t> 
			cdata ( data.size(), 0 );

		if ( ( std::is_same<value_type, float>::value ) || ( std::is_same<value_type, float>::value ) )
			// proper quantisation for floating point data types -- otherwise too many checks needed
			std::transform ( std::execution::seq, data.begin(), data.end(), quant.begin(),
			[&mn, &mx, &levels] ( auto input ) {
			return ( input - mn ) * levels / ( mx - mn );
		} );
		else
			// simpler algorithm for integers: first check case if our image fits in "levels" - in which case just copy
			if ( ( mx - mn ) < levels )
				std::transform ( std::execution::seq, data.begin(), data.end(), quant.begin(),
				[&mn] ( auto input ) {
				return ( input - mn );
			} );
		else {
			// if the image contains more unique intensity values than "levels", divide by the proper scalar
			auto factor = static_cast<float> ( mx - mn ) / levels;
			std::transform ( std::execution::seq, data.begin(), data.end(), quant.begin(),
			[&mn, &factor] ( auto input ) {
				return ( ( input - mn ) / factor );
			} );
		}

		// the parent vector
		std::vector<size_t> parent ( quant.size(), undefined );

		////////////////////////////////////////////////////////////////////////////////
		//
		// build the max-tree from the quantised image
		// using the method from Berger's ICIP 2007 paper
		// http://dx.doi.org/10.1109/ICIP.2007.4379949
		//

		if ( method == "Berger" ) {

			std::vector<size_t>	zpar ( quant.size(), undefined );
			std::vector<size_t>	root ( quant.size(),		 0 );
			std::vector<size_t>	rank ( quant.size(),		 0 );

			// std::string letters = "CDHAFBIGEJ"; // from Berger's 2007 ICIP paper
			for ( size_t i = 0; i < indices.size(); i++ ) { 

				size_t 
					p     = indices[ indices.size() - i -1 ]; // point at index, step from top - bottom
				parent[p] = p;          // pixel at this (higher level) starts as parent
				zpar  [p] = p;          //							as union-find parent
				root  [p] = p;

				auto x    = p;          // keep this as zpar

				// this is valid for Bergers 2007 ICIP example
				// std::cout << "level " << quant[p] << ", procesing node " << p << " (" << letters[p] << ")" <<
				// std::endl;			

				for ( unsigned k = 0; k < neighbours.size(); k++ ) {	// k: neighbour offset
					long long n = p + neighbours[k];					// n: storage position of neighbour

					if ( ( n > -1 )											&& 
					     ( ( static_cast<size_t> ( n ) ) < indices.size() ) && 
						 this->valid_neighbours ( n, p, 1 ) ) {

						if ( zpar[n] !=undefined ) {					// if n has been visited it has a zpar

							// this is valid for Bergers 2007 ICIP example
							// std::cout << "looking at neighbour " << n << " (" << letters[n] << ")";

							size_t r = static_cast<size_t> ( n );		// r = root: index of neighbour q
							while ( r != zpar[r] ) 						//     whose zpar points to itself
								r = zpar[r];							//     (short version of 'findroot')

							if ( r != x ) {
								
								parent[ root [ r ] ] = p;
								
								if ( rank [ x ] < rank [ r ] ) {
									
									zpar [ x ] = r;
									root [ r ] = p;
									         x = r;
											 
								} else {
									
									zpar [ r ] = p;

									if ( rank [ r ] >= rank [ p ] ) 
										rank [ p ] += 1;										
									
								}
								
								// this is valid for Bergers 2007 ICIP example
								// std::cout << " whose new root is now " << p << " (" << letters[p] << ")";

							} // if r and p need joining

							// std::cout << std::endl;

						} // if n has a root r

					} // if neighbour n in image

				} // for n neighbours n of p

			} // for p

			////////////////////////////////////////////////////////////////////////////////
			//
			// link to level roots
			//
			
			for ( auto p : indices ) {
				auto q = parent[p];
				if ( quant[parent[q]] == quant[q] )
					parent[p] = parent[q];
			} // for pi

			////////////////////////////////////////////////////////////////////////////////
			//
			// identify components ( ≥1 per level )
			//
			
			level_t
				ccount  = 0;
			cdata [ indices [ 0 ] ] = 0;
			for ( auto p : indices ) {
				if ( quant [ parent [ p ] ] == quant [ p ] )
					cdata [ p ] = cdata [ parent [ p ] ];
				else
					cdata [ p ] = ++ccount;
			} // for pi
			
			////////////////////////////////////////////////////////////////////////////////
			//
			// store components: { number, intensity, area, pos. root, number parent }
			//
			components.resize ( ccount + 1 );
			for ( size_t p = 0; p < cdata.size(); p++ ) {
				auto  c = cdata [ p ];
				if ( ! 	components [ c ].size ) {
						components [ c ].root   = p;
						components [ c ].value  = quant                    [ p ];
						components [ c ].parent = cdata  [ parent [ parent [ p ] ] ];			
				}
				components [ c ].size++;
				components [ c ].points.push_back ( p );
			} 
			
			////////////////////////////////////////////////////////////////////////////////
			//
			// add higher components to lower
			//
			for ( size_t c = 0; c < components.size(); c++ ) 
				components [ c ].uniq = components [ c ].size;	
			for ( size_t c = components.size() - 1; c > 0; c-- ) {
				components [ components [ c ].parent ].size += components [ c ].size;
				components [ components [ c ].parent ].children.push_back ( c );
			}
			for ( size_t c = 0; c < components.size(); c++ ) 
				std::reverse ( components [ c ].children.begin(), components [ c ].children.end() );

			//std::copy ( quant.begin(), quant.end(), data.begin() );
			//std::copy ( parent.begin(), parent.end(), data.begin() );
			std::copy ( parent.begin(), parent.end(), data.begin() );

			
		} // if Berger method used

	} // constructor from image

	/** \brief (deep) copy constructor
	 *
	 * copies from an existing bismaxtree
	 */
	bismaxtree ( const bismaxtree& rhs ) {
	}

	/** \brief assignment operator
	 *
	 * assigns contents of the right-hand side (RHS) maxtree to (*this)
	 */
	const bismaxtree<value_type>& operator= ( bismaxtree<value_type>& rhs ) {
	} // assignment

	/** \brief get a component's points
	 *
	 * assigns contents of the right-hand side (RHS) maxtree to (*this)
	 */
	const std::vector<size_t> getpoints ( size_t comp_start, size_t comp_end = 0, bool sorted = false ) {
		
		std::vector<size_t> 
			mypoints;
		auto chigh = ( ! comp_end ) ? components.size() : comp_end;
			
		mypoints.insert ( mypoints.end(), components [ comp_start ].points.begin(), components [ comp_start ].points.end() );
		for ( auto c: components [ comp_start ].children ) {
			// if ( c <= chigh ) {
				auto vec = getpoints( c, chigh, sorted );
				mypoints.insert ( mypoints.end(), vec.begin(), vec.end() );
			}
		
		if (sorted)
			std::sort( mypoints.begin(), mypoints.end() );
		
		return ( mypoints );
		
	} // assignment

	/** \brief return a bisimage with the intensities of selected components
	 *
	 * sets only the points in the image that belong to certain components
	 * 
	 */
	bisimage<value_type> setpoints ( size_t comp_start, size_t comp_end = 0 ) {
		
		std::vector<size_t> 
			mypoints = getpoints ( comp_start, comp_end, true );
		std::vector<char> 
			found ( superclass::data.size(), 0 );
		bisimage<value_type>
			output = (*this);
		
		for ( auto p: mypoints )
			found [ p ] = 1;
		output.vector_mask ( found );

		return ( output );
		
	} // assignment

	/*
	* make std::cout so that it summarises the tree
	* 
	*/
	friend std::ostream& operator<< ( std::ostream& sout, const bismaxtree& output) {
		
		for ( size_t c = 0; c < output.components.size(); c++ ) {
			
			sout	<< " component "   << std::setfill ( ' ' ) << std::setw ( 3 ) << c 
					<<  ": { size: "   << std::setfill ( ' ' ) << std::setw ( 3 ) << output.components[c].size    
					<<  ", (unique: "  << std::setfill ( ' ' ) << std::setw ( 3 ) << output.components[c].uniq
					<< "), root: "     << std::setfill ( ' ' ) << std::setw ( 3 ) << output.components[c].root   
					<<  ", value: "    << std::setfill ( ' ' ) << std::setw ( 3 ) << output.components[c].value   
					<<  ", parent: "   << std::setfill ( ' ' ) << std::setw ( 3 ) << output.components[c].parent
					<<  ", children: " << std::setfill ( ' ' ) << std::setw ( 3 ) << output.components[c].children << " }" 
					<< std::endl;
					
			sout	<< "points: " << output.components [ c ].points << std::endl;
			
		}

		return ( sout );
		
	}



}; // class bismaxtree

} // namespace bis

#endif // BISNIFTIIMAGE_HPP_INCLUDED
