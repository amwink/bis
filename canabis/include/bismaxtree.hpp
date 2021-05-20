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
	using self = bismaxtree<value_type>;
	using superclass = bisimage<value_type>;

	using superclass::data;
	using superclass::sizes;
	using superclass::strides;

public:
	/** \brief default constructor
	 *
	 * (used by subclass bisdicom)
	 *
	 */
	bismaxtree() {
	}

	/** \brief destructor
	 *
	 * destructor: clears vectors and (if
	 * required) frees pointer to header
	 */
	~bismaxtree() {
	}

	/** \brief constructor from a bisimage
	 *
	 * This constructor takes an n-dimensional image
	 * and bulds its maxtree representation
	 */
	typedef int level_t;
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

		std::vector<long> neighbours;
		switch ( connectivity ) {
		case 4:
			assert ( sizes.size() == 2 );
			neighbours = { -strides[1], -1, 1, strides[1] };
			break;
		case 8:
			assert ( sizes.size() == 2 );
			neighbours = { -strides[1] - 1, -strides[1], -strides[1] + 1, 
									   - 1,							   1,
							strides[1] - 1,	 strides[1],  strides[1] + 1
			             };
			break;
		case 6:
			assert ( sizes.size() == 3 );
			neighbours = { -strides[2], -strides[1], -1, 1, strides[1], strides[2] };
			break;
		case 26:
			assert ( sizes.size() == 3 );
			neighbours = { -strides[2] - strides[1] - 1, -strides[2] - strides[1], -strides[2] - strides[1] + 1,
			               -strides[2] - 1,              -strides[2],                           -strides[2] + 1, 
						   -strides[2] + strides[1] - 1, -strides[2] + strides[1], -strides[2] + strides[1] + 1, 
						                -strides[1] - 1,              -strides[1],              -strides[1] + 1, 
						                             -1,                                                      1, 
										 strides[1] - 1,			   strides[1], 				 strides[1] + 1, 
						   strides[2] - strides[1] - 1,   strides[2] - strides[1],	strides[2] - strides[1] + 1, 
						   strides[2] - 1, 				  strides[2], 				strides[2]				+ 1, 
						   strides[2] + strides[1] - 1,	  strides[2] + strides[1],	strides[2] + strides[1] + 1
			             };
			break;
		}

		////////////////////////////////////////////////////////////////////////////////
		//
		// sort sorted intensities in the data vector "sorted"
		// and also keep an array "indices" to see where they were in the image
		//

		std::vector<size_t> indices ( data.size() );

		std::iota ( indices.begin(), indices.end(), 0 );  // fill with 0, 1, 2, ..
		if ( method == "Berger" )                         // to label last index as
			std::reverse ( indices.begin(), indices.end() ); // root (ICIP 2007 paper)
		std::stable_sort ( indices.begin(), indices.end(), // 'data' sorted for determining 'indices'
						   [&] ( size_t i, size_t j ) { return ( data[i] < data[j] ); } );

		auto mn = data[indices.front()], mx = data[indices.back()];

		////////////////////////////////////////////////////////////////////////////////
		//
		// quantise the image, store in "quant"
		//

		std::vector<level_t> quant ( data.size() );

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
		std::vector<long> parent ( quant.size(), -1 );

		////////////////////////////////////////////////////////////////////////////////
		//
		// build the max-tree from the quantised image
		// using the method from Berger's ICIP 2007 paper
		// http://dx.doi.org/10.1109/ICIP.2007.4379949
		//

		if ( method == "Berger" ) {

			std::vector<long> zpar ( quant.size(), -1 );

			// std::string letters = "CDHAFBIGEJ"; // from Berger's 2007 ICIP paper
			for ( long i = indices.size() - 1; i >= 0; i-- ) { // i: index ( sorted from high to low )

				long p = indices[i]; // point at index
				parent[p] = p;       // pixel at this (higher level) starts as parent
				zpar[p] = p;         //							as union-find parent

				// this is valid for Bergers 2007 ICIP example
				// std::cout << "level " << quant[p] << ", procesing node " << p << " (" << letters[p] << ")" <<
				// std::endl;

				for ( int k = 0; k < neighbours.size(); k++ ) { // k: neighbour offset
					long n = p + neighbours[k];              // q: neighbour position

					if ( ( n > -1 ) && ( n < indices.size() ) && this->valid_neighbours ( n, p ) ) {

						if ( zpar[n] > -1 ) { // if q has been visited it has a zpar

							// this is valid for Bergers 2007 ICIP example
							// std::cout << "looking at neighbour " << n << " (" << letters[n] << ")";

							long r = n;         // r = root: index of neighbour q
							while ( r != zpar[r] ) //     whose zpar points to itself
								r = zpar[r];

							if ( r != p ) {
								parent[r] = p;
								zpar[r] = p;

								// this is valid for Bergers 2007 ICIP example
								// std::cout << " whose new root is now " << p << " (" << letters[p] << ")";

							} // if r and p need joining

							std::cout << std::endl;

						} // if n has a root r

					} // if neighbour n in image

				} // for n neighbours n of p

			} // for p

			////////////////////////////////////////////////////////////////////////////////
			//
			// Canonisation (platslaan)
			//

			for ( auto p : indices ) {
				auto q = parent[p];
				if ( quant[parent[q]] == quant[q] )
					parent[p] = parent[q];
			} // for pi
			
			
			
		} else { // if Wilkinson method used



			////////////////////////////////////////////////////////////////////////////////
			//
			// build the max-tree from the quantised image
			// using the method from Wilkinson's paper
			// 
			//

			// stack for flooding levels
			std::vector<long> 
				level_vec ( levels );
			std::stack<long, std::vector<long>> 
				level_stack ( std::move ( level_vec ) );

			// priority queue for filtered pixel values
			std::priority_queue<long, std::vector<long>, std::less<size_t>> 
				node_queue ( quant.begin(), quant.end() );

			// apart from a node's parent, we also store the area...
			std::vector <long> 
				area ( quant.size(), 1 );
			// ... and a flag that the node has been prcessed
			std::vector <bool>
				processed ( quant.size(), false );

			// init tree: push the index of a minimal intensity onto queue and stack
			// indices [ 0 ] is the first index in the image with the lowest intensity-
			// lowest because the indices are of the sorted intensities, and first be-
			// cause a stable sort has been used.
			
			// node_queue.push  ( indices[0] ); // queue already has all
			level_stack.push ( indices[0] );
			parent [ indices [ 0] ] = 0;

			while ( !node_queue.empty() ) {

				bool ready_to_pop = true;

				// start flooding
				auto p = node_queue.top();       // pixel that represents a level   (in the filter queue)
				auto l = level_stack.top();      // level as represented by a pixel (in the levels stack)
				assert ( quant[p] == quant[l] ); // so they should be the same at the start of flooding



				for ( int k = 0; k < neighbours.size(); k++ ) { // k: neighbour offset
					long n = p + neighbours[k];                 // q: neighbour position

					// only continue if the neighbour position is inside the image
					if ( ( n > -1 ) && ( n < indices.size() ) && this->valid_neighbours ( n, p ) ) {

						if ( parent[n] == -1 ) {    // if neighbour not yet visited -- no curent parent

							node_queue.push ( n );  // then place neighbour in the queue
							parent[n] = 0;          // and change parent to intermediate 

							if ( quant[p] < quant[n] ) {  // if one of the unprocessed neighbours is higher
								level_stack.push ( n );   // leave 'for', return to the start of 'while'
								ready_to_pop = false;
								break;							
							} // if quant[n] < quant[p]

						} // if processed

					} // if valid neighbours
					
				} // for neighbours



				if ( ready_to_pop ) {

					node_queue.pop();
					parent[p] = l;
					
					
					
				}



			} // while node_queue not empty
			
		} // if wilkinson method

		std::cout << "data" << std::endl;
		std::cout << ( *this ) << std::endl;

		std::cout << "quant" << std::endl;
		std::copy ( quant.begin(), quant.end(), data.begin() );
		std::cout << ( *this ) << std::endl;

		std::cout << "indices" << std::endl;
		std::copy ( indices.begin(), indices.end(), data.begin() );
		std::cout << ( *this ) << std::endl;

		std::cout << "parent" << std::endl;
		std::copy ( parent.begin(), parent.end(), data.begin() );
		std::cout << ( *this ) << std::endl;

		std::cout << "done." << std::endl;

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

}; // class bismaxtree

} // namespace bis

#endif // BISNIFTIIMAGE_HPP_INCLUDED
