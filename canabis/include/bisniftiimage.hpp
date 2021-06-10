#ifndef BISNIFTIIMAGE_HPP_INCLUDED
#define BISNIFTIIMAGE_HPP_INCLUDED

/** \brief bisniftiimage is a subclass of bisimage
 */
#include "bisimage.hpp"

/** \brief bisniftiimage.hpp uses nifti2.h for defining
 *         the latest version of the NIfTI header.
 */
#include "nifti2_io.h"

/** \brief bis: namespace for multidimensional images
 *
 *  This namespace contains the bisimage class for
 *  creating and processing multidimensional images.
 */
namespace bis {

    /** \brief niftiExport -- typename to memory buffer (in nifti record)
     *
     *  template parameters:
     *  BufElem:    the element type of the input buffer (e.g., float)
     *  Data:       the target container (e.g., a vector<double>)
     *
     *  function parameters:
     *  nifti_image:    a nifti_image record with a pointer to an output data buffer
     *  data:           pointer to the start of the target container's buffer
     *
     */

    template <typename BufElem, typename Data>

    void niftiExport ( nifti_image* nim, Data* data ) {

        const unsigned bufsize = data->size();
        BufElem* buffer = new BufElem[bufsize];
        Data& dataref = *data;

        for ( unsigned i = 0; i < bufsize; i++ )
            buffer[i] = ( BufElem ) dataref[i];

        nim->data = ( void* ) buffer;
    };

    /** \brief getNiftiBricks: import bricks from nifti file
     *
     *  template parameters:
     *  DataVec:        the target container (e.g., a vector<double>)
     *
     *  function parameters:
     *  nim:            pointer to a nifti_image
     *  nifti_blob:     pointer to the data buffer belonging to the nifti_image
     *  bufsize:        number of pixels
     *  vec:            target container
     *
     */
    template <typename DataVec>

    void getNiftiBricks ( nifti_image* nim, void* nifti_blob, unsigned bufsize, DataVec* vec ) {

        switch ( nim->datatype ) {

            case ( NIFTI_TYPE_UINT8 ) :
                pixelImport<unsigned char> ( nifti_blob, bufsize, vec );
                break;
            case ( NIFTI_TYPE_INT16 ) :
                pixelImport<signed short> ( nifti_blob, bufsize, vec );
                break;
            case ( NIFTI_TYPE_INT32 ) :
                pixelImport<signed int> ( nifti_blob, bufsize, vec );
                break;
            case ( NIFTI_TYPE_FLOAT32 ) :
                pixelImport<float> ( nifti_blob, bufsize, vec );
                break;
            /*case(NIFTI_TYPE_COMPLEX64):
                pixelImport<complex_float> ( nifti_blob, bufsize, vec);
                break;*/
            case ( NIFTI_TYPE_FLOAT64 ) :
                pixelImport<double> ( nifti_blob, bufsize, vec );
                break;
            /*case(NIFTI_TYPE_RGB24):
                pixelImport<rgb_byte> ( nifti_blob, bufsize, vec);
                break;*/
            case ( NIFTI_TYPE_INT8 ) :
                pixelImport<signed char> ( nifti_blob, bufsize, vec );
                break;
            case ( NIFTI_TYPE_UINT16 ) :
                pixelImport<unsigned short> ( nifti_blob, bufsize, vec );
                break;
            case ( NIFTI_TYPE_UINT32 ) :
                pixelImport<unsigned int> ( nifti_blob, bufsize, vec );
                break;
            case ( NIFTI_TYPE_INT64 ) :
                pixelImport<signed long long> ( nifti_blob, bufsize, vec );
                break;
            case ( NIFTI_TYPE_UINT64 ) :
                pixelImport<unsigned long long> ( nifti_blob, bufsize, vec );
                break;
            case ( NIFTI_TYPE_FLOAT128 ) :
                pixelImport<long double> ( nifti_blob, bufsize, vec );
                break;
            /*case(NIFTI_TYPE_COMPLEX128):
              // pixelImport<complex_double> ( nifti_blob, bufsize, vec);
              // break;*/
            /*case(NIFTI_TYPE_COMPLEX256):
              // pixelImport<complex_longdouble> ( nifti_blob, bufsize, vec);
              // break;*/
            /*case(NIFTI_TYPE_RGBA32):
              // pixelImport<> ( nifti_blob, bufsize, vec);
              // break;*/
            default:
                std::cout << nim->fname << " has unsupported data type " << nim->datatype << std::endl;
                break;
        } /* switch ( nim->datatype )                                          */

    } /* getNiftiBricks()                                                          */

    /** \brief setNiftiBricks: export bricks to nifti file
     *  template parameters:
     *  DataVec:        the input container (e.g., a vector<double>)
     *
     *  function parameters:
     *  nim:            pointer to a nifti_image
     *  vec:            input container
     *
     */
    template <typename DataVec>

    void setNiftiBricks ( nifti_image* nim, DataVec* vec ) {

        switch ( nim->datatype ) {

            case ( NIFTI_TYPE_UINT8 ) :
                niftiExport<unsigned char> ( nim, vec );
                break;
            case ( NIFTI_TYPE_INT16 ) :
                niftiExport<signed short> ( nim, vec );
                break;
            case ( NIFTI_TYPE_INT32 ) :
                niftiExport<signed int> ( nim, vec );
                break;
            case ( NIFTI_TYPE_FLOAT32 ) :
                niftiExport<float> ( nim, vec );
                break;
            /*case(NIFTI_TYPE_COMPLEX64):
              niftiExport<complex_float> ( nim, vec);
              break;*/
            case ( NIFTI_TYPE_FLOAT64 ) :
                niftiExport<double> ( nim, vec );
                break;
            /*case(NIFTI_TYPE_RGB24):
              niftiExport<rgb_byte> ( nim, vec);
              break;*/
            case ( NIFTI_TYPE_INT8 ) :
                niftiExport<signed char> ( nim, vec );
                break;
            case ( NIFTI_TYPE_UINT16 ) :
                niftiExport<unsigned short> ( nim, vec );
                break;
            case ( NIFTI_TYPE_UINT32 ) :
                niftiExport<unsigned int> ( nim, vec );
                break;
            case ( NIFTI_TYPE_INT64 ) :
                niftiExport<signed long long> ( nim, vec );
                break;
            case ( NIFTI_TYPE_UINT64 ) :
                niftiExport<unsigned long long> ( nim, vec );
                break;
            case ( NIFTI_TYPE_FLOAT128 ) :
                niftiExport<long double> ( nim, vec );
                break;
            /*case(NIFTI_TYPE_COMPLEX128):
              niftiExport<complex_double> ( nim, vec);
              break;*/
            /*case(NIFTI_TYPE_COMPLEX256):
              niftiExport<complex_longdouble> ( nim, vec);
              break;*/
            /*case(NIFTI_TYPE_RGBA32):
              niftiExport<> ( nim, vec);
              break;*/
            default:
                std::cout << nim->fname << " has unsupported data type " << nim->datatype << std::endl;
                break;
        } /* switch ( nim->datatype )                                          */

    } /* setNiftiBricks()                                                          */

    /** \brief get the correct nifti datatype for T the header
     *
     */
    template <typename T> unsigned getdatatype() {

        unsigned data_out = DT_UNKNOWN;

        if ( std::is_same<T, unsigned char>::value )
            data_out = NIFTI_TYPE_UINT8;
        if ( std::is_same<T, char>::value )
            data_out = NIFTI_TYPE_INT8;
        if ( std::is_same<T, unsigned short>::value )
            data_out = NIFTI_TYPE_UINT16;
        if ( std::is_same<T, short>::value )
            data_out = NIFTI_TYPE_INT16;
        if ( std::is_same<T, unsigned int>::value )
            data_out = NIFTI_TYPE_UINT32;
        if ( std::is_same<T, int>::value )
            data_out = NIFTI_TYPE_INT32;
        if ( std::is_same<T, unsigned long int>::value )
            data_out = NIFTI_TYPE_UINT64;
        if ( std::is_same<T, long int>::value )
            data_out = NIFTI_TYPE_INT64;
        if ( std::is_same<T, float>::value )
            data_out = NIFTI_TYPE_FLOAT32;
        if ( std::is_same<T, double>::value )
            data_out = NIFTI_TYPE_FLOAT64;
        if ( std::is_same<T, long double>::value )
            data_out = NIFTI_TYPE_FLOAT128;

        return ( data_out );

    } // getdatatype

    template <typename value_type>

    class bisnifti : public bisimage<value_type> {

            // self and superclass
            using self = bisnifti<value_type>;
            using superclass = bisimage<value_type>;

        protected:
            /** \brief the header should be usable by subclas bisdicom
             *
             * header:  pointer to a NIfTI header
             *          for reading / writing files
             */
            nifti_image* header = nullptr;

        public:
            /** \brief default constructor
             *
             * (used by subclass bisdicom)
             *
             */
            bisnifti() {
            }

            /** \brief destructor
             *
             * destructor: clears vectors and (if
             * required) frees pointer to header
             */
            ~bisnifti() {

                superclass::data.resize ( 0 );
                superclass::sizes.resize ( 0 );
                superclass::strides.resize ( 0 );
                if ( header != NULL )
                    free ( header );
            }

            /** \brief constructor from a NIfTI image
             *
             * This constructor reads an n-dimensional image
             * ( for NIfTI 1 <= n <= 7 ) from a NIfTI file and
             * assigns its data, sizes and header to (*this).
             */
            bisnifti ( std::string filename, bis::readdata readornot = bis::DO_READ_DATA ) {

                header = nifti_image_read ( filename.c_str(), readornot );

                superclass::sizes.resize ( header->dim[0] );
                superclass::strides.resize ( header->dim[0] + 1 );

                // make the array 'strides' so that it uses the last dimension as well
                superclass::strides[0] = 1;
                for ( size_t i = 1; i <= superclass::sizes.size(); i++ ) {
                    superclass::sizes[i - 1] = header->dim[i];
                    superclass::strides[i] = superclass::strides[i - 1] * superclass::sizes[i - 1];
                }

                if ( readornot == bis::DO_READ_DATA ) {

                    superclass::data.resize ( * ( superclass::strides.rbegin() ) ); // the end of strides holds the image's size

                    bis::getNiftiBricks ( header, header->data, superclass::data.size(), superclass::getdata_ptr() );

                } // if readornot
            }

            /** \brief (deep) copy constructor
             *
             * copies data, sizes and header from an existing bisnifti
             */
            bisnifti ( const bisnifti& rhs ) {

                if ( rhs.header != NULL )
                    header = nifti_copy_nim_info ( rhs.header );

                superclass::data.resize ( rhs.data.size() );
                superclass::sizes.resize ( rhs.sizes.size() );
                superclass::strides.resize ( rhs.strides.size() );

                std::copy ( rhs.data.begin(), rhs.data.end(), superclass::data.begin() );
                std::copy ( rhs.sizes.begin(), rhs.sizes.end(), superclass::sizes.begin() );
                std::copy ( rhs.strides.begin(), rhs.strides.end(), superclass::strides.begin() );
            }

            /** \brief constructor from a bisImage
             *
             * copies data and sizes from an existing bisimage, create header
             */
            bisnifti ( bisimage<value_type>& rhs ) {

				int64_t 
					dims [ 8 ] = { 0, 0, 0, 0, 0, 0, 0, 0 };
				
                superclass::data.resize ( rhs.getdatasize() );
                superclass::sizes = rhs.getsize();
                superclass::validate_sizes();
                std::copy ( rhs.getdata_ptr()->begin(), rhs.getdata_ptr()->end(), superclass::data.begin() );

                if ( header != NULL )
                    free ( header );
					
				dims [ 0 ] = superclass::sizes.size();	
				for ( unsigned d = 0; d < dims [ 0 ]; d++ )
					dims [ d + 1 ] = superclass::sizes [ d ];
				
				header  = nifti_make_new_nim ( dims, getdatatype<value_type>(), 0 ); // the 0 is for not creating the intensities array
				header -> data = superclass::getdata_ptr();

            }
			
            /** \brief assignment operator
             *
             * assigns data, sizes and header of the
             * right-hand side (RHS) image to (*this)
             */
            const bisnifti<value_type>& operator= ( bisnifti<value_type>& rhs ) {

                if ( this != &rhs ) {

                    // just to make sure we don't leave stuff
                    if ( this->header != NULL ) {
                        free ( header );
                        header = NULL;
                    }

                    if ( rhs.header != NULL )
                        header = nifti_copy_nim_info ( rhs.header );

                    // just to make sure we don't leave stuff
                    superclass::data.resize ( rhs.data.size() );
                    superclass::sizes.resize ( rhs.sizes.size() );
                    superclass::strides.resize ( rhs.strides.size() );

                    std::copy ( rhs.data.begin(), rhs.data.end(), superclass::data.begin() );
                    std::copy ( rhs.sizes.begin(), rhs.sizes.end(), superclass::sizes.begin() );
                    std::copy ( rhs.strides.begin(), rhs.strides.end(), superclass::strides.begin() );

                } // if this != rhs

                return *this;

            } // assignment

            /** \brief change the data type for writing a NIfTI file
             *
             * The data type is changed to e.g. float for
             * more flexibility (at the cost of file size).
             */
            void setNIIdatatype ( unsigned dtype ) {

                header->datatype = dtype;
                nifti_datatype_sizes ( header->datatype, &header->nbyper, &header->swapsize );
            }

            /** \brief write to a NIfTI image
             *
             * Writes the contents of (*this) to a NIfTI
             * file using the data and header information.
             * The data type is changed to float for more
             * flexibility (at the cost of file size).
             *
             * NIfTI routines are run on a temporary copy
             *       as they seem to be memory-unsafe
             */
            void saveNII ( std::string filename = "" ) {

                if ( ( filename != "" ) || ( strlen ( header->fname ) > 0 ) ) {

                    if ( filename == "" )
                        filename = header->fname;

                    // for now, make sure we write nifti-1
                    // (not many programs support i/o for nifti-2)
                    header->nifti_type = NIFTI_FTYPE_NIFTI1_1;
                    // header -> nifti_type = NIFTI_FTYPE_NIFTI1_2;

                    // also set header offset, as this is not automatically updated
                    nifti_set_iname_offset ( header, 1 );

                    nifti_set_filenames ( header, filename.c_str(), 0, 0 );
                    setNiftiBricks ( header, & ( superclass::data ) );

                    auto checkit = true;
                    if ( checkit ) {
                        std::cout << "data size: " << superclass::data.size() << std::endl;
                        std::cout << "sizes:     " << superclass::sizes << std::endl;
                        std::cout << "strides:   " << superclass::strides << std::endl;
                        //auto checkit2 = false;
                        //if ( checkit2 ) {
                        //    auto pp = header->data;
                        //    for ( auto p = superclass::data.begin(); p < superclass::data.end(); p += 2000, pp += 2000 ) {
                        //        std::cout << "position " << p - superclass::data.begin() << ", bis value " << *p;
                        //        std::cout << ", nii value " << * ( unsigned short* ) ( pp ) << std::endl;
                        //    }
                        //}
                        infodump();
                        auto checkit3 = false;
                        if ( checkit3 )
                            std::cin.ignore();
                    }

                    nifti_image_write ( header );

                } else

                    std::cerr << "saveNII: no file name for image given";
            }

            /** \brief getinfo() - returns the address of the
             *
             * this is a pointer not the nifti record itself -- use with care
             */
            nifti_image* getinfo() {
                return header;
            }

            /** \brief infodump() - returns the address of the
             *
             * do an infodump to stderr, just using the basic nifti functionality
             */
            void infodump() {
                nifti_image_infodump ( getinfo() );
            }

    }; // class bisnifti

} // namespace bis

#endif // BISNIFTIIMAGE_HPP_INCLUDED
