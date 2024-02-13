#ifndef BISBIDSIMAGE_HPP_INCLUDED
#define BISBIDSIMAGE_HPP_INCLUDED

#include <regex>



/** \brief bisbidsimage is a subclass of bisniftiimage
 */
#include "bisniftiimage.hpp"



/** \brief DICOM headers may contain BIDS data fields that
 *         cannot be stored in NIfTI headers; these are
 *         stored in a JSON 'sidecar' with the same name
 */
#include <nlohmann/json.hpp>
using json = nlohmann::json;



/** \brief bis: namespace for multidimensional images
 *
 *  This namespace contains the bisimage class for
 *  creating and processing multidimensional images.
 */



namespace bis {
		
	std::string sidecarname ( std::string niftifile ) {
	
		std::string 
			jsonfname = std::regex_replace ( niftifile, std::regex ( ".nii.gz" ), ".json" );  // replace .nii.gz for compressed files
		jsonfname = std::regex_replace ( jsonfname,  std::regex ( ".nii" ),    ".json" );     // in case it wasn't compressed
	
		return ( jsonfname );
	
	}
	
	
	
    template <typename value_type>

    class bisbids : public bisnifti<value_type> {

		// self (own class) and superclasses
		using self       = bisbids <value_type>;
		using superclass = bisnifti<value_type>;
		using supersuper = bisimage<value_type>;

        protected:

            json 
                sidecar; // JSON record that holds BIDS information not stored in the NIfTI header

        public:
        
			bisbids () : bisnifti <value_type> () {}; // default constructor (required for bisdicom subclass)
		
            bisbids ( std::string fname ) : bisnifti <value_type> ( fname ) { // constructor that reads from a file

				std::string	jsonfile = sidecarname ( fname );
				getjson ( jsonfile );
			
            } // constructor that reads from a file

			void getjson ( std::string json_fname ) {
				
				std::ifstream 
					j_in ( json_fname );
				sidecar = json::parse ( j_in );					

			}

			// overloaded function from bisnifti
			void read ( std::string fname ) {
				
				superclass::read( fname );
				std::string	jsonfile = sidecarname ( fname );
				getjson ( jsonfile );
				
			}

			// write the sidecar
			void write_sidecar ( std::string niiname ) {
				
				std::ofstream jfile ( sidecarname ( niiname ) );
				jfile << std::setw ( 4 ) << sidecar << std::endl;
				jfile.close();
				
			}

			// overloaded function from bisnifti
			void write ( std::string filename = "" ) {
				
				superclass::write ( filename );
				write_sidecar( filename );
				
			}



    }; // class bisbids

} // namespace bis

#endif // BISBIDSIMAGE_HPP_INCLUDED