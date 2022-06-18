#ifndef BISDICOMIMAGE_HPP_INCLUDED
#define BISDICOMIMAGE_HPP_INCLUDED

/** \brief use the standard library for trawling dicom directories
 */
#include <filesystem>
#include <stdlib.h>		// strtof
#include <iomanip>      // std::setprecision

/** \brief bisniftiimage is a subclass of bisimage
 */
#include "bisbidsimage.hpp"

/** \brief bisdicomimage.hpp uses dcmtk for DICOM file I/O
 *         (basically input), in a single function
 */
#include "dcmtk/config/osconfig.h"
#include "dcmtk/dcmdata/dcistrmf.h"
#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmimgle/dcmimage.h"

/** \brief bis: namespace for multidimensional images
 *
 *  This namespace contains the bisimage class for
 *  creating and processing multidimensional images.
 */

namespace bis {

// for gathering the slices, which may be in different files
typedef struct slicestats { 	// make a struct:
    std::string filename;		//	string	filename, for single-slice dicom every slice has a unique file name
    unsigned innumber;			//	int	instance number: should be from 1 to #files
    unsigned acqnumber;			//	int	acquisition number
    float aqtime;				//	float	acquisition time
    float slipos;				//	float	slice position, computed from ImagePOsitionPatient and ImageOrientationPatient
	float sliloc;				//	float 	slice location, backup method sometimes (eg Philips DTI) using SliceLocation
	unsigned index;				//	int	index, for multi-slice dicoms (e.g. mosaic), this number can be higher than 0
} slicestats;

/** \brief getItemString: get a string value from a given Dicom tag
 *
 * parameters:
 *
 * item:      pointer to the DICOM item that is queried
 * tagKey:    requested DICOM tag as named in DCMTK
 * position:  if the field has more than 1 values, get this position
 * searchsub: if OFTrue (default) then search into header subsequences
 *
 * return value
 * std::string with the contents of the DICOM field
 *
 */
std::string
    getItemString(DcmItem* item, const DcmTagKey& tagKey, const size_t position = 0, const OFBool searchsub = OFFalse) {

    OFString tmpofstring;

    if(item->findAndGetOFString(tagKey, tmpofstring, position, searchsub).good())
	return std::string(tmpofstring.c_str());
    else
	return "-1";
}

/** \brief getItemStringArray: get a string value from a given Dicom tag
 *
 * parameters:
 *
 * item:      pointer to the DICOM item that is queried
 * tagKey:    requested DICOM tag as named in DCMTK
 *
 * return value
 * std::string with the contents of the DICOM field
 *
 */
std::string getItemStringArray(DcmItem* item, const DcmTagKey& tagKey) {

    OFString tmpofstring;

    if(item->findAndGetOFStringArray(tagKey, tmpofstring).good())
	return std::string(tmpofstring.c_str());
    else
	return "-1";
}

/** \brief getItemDouble: get a float64 value from a given Dicom tag
 *
 * parameters:
 *
 * item:      pointer to the DICOM item that is queried
 * tagKey:    requested DICOM tag as named in DCMTK
 * position:  if the field has more than 1 values, get this position
 * searchsub: if OFTrue (default) then search into header subsequences
 *
 * return value
 * double with the contents of the DICOM field
 *
 */
double
    getItemDouble(DcmItem* item, const DcmTagKey& tagKey, const size_t position = 0, const OFBool searchsub = OFTrue) {

    double tmpdouble;

    if(item->findAndGetFloat64(tagKey, tmpdouble, position, searchsub).good())
	return (tmpdouble);
    else
	return -1.;
}

/** \brief getDicomData: import DICOM slices/bricks/timeseries from file
 *
 *  template parameters:
 *  DataVec:        the target container (e.g., a vector<double>)
 *
 *  function parameters:
 *  nim:            pointer to a nifti_image
 *  nifti_blob:     pointer to the data buffer belonging to the nifti_image
 *
 *  filelist:       list of input dicom images
 *  bufsize:        number of pixels
 *  vec:            target container
 *
 */
template <typename DataVec>

void getDicomData(std::vector<std::string> file_list, unsigned bufsize, DataVec* vec) {

    std::vector<unsigned char> buffer;      // for reading binary data
    size_t slicesize = 0, bufferoffset = 0; // for reading pixels
    int bitsperpixel = 8;                   // usually changed to 16 (intensities to about 2000)
    unsigned bytesperpixel = 1;             // usually changed to 2

    for(size_t im = 0; im < file_list.size(); im++) {

	DicomImage* image = new DicomImage(file_list[im].c_str());

	// determine some stats in the first image
	if(!im) {

	    auto bpp = image->getDepth(), // assume this is the same for all
	        bitsperpixel = (bpp > 32) ? 64 : (bpp > 16) ? 32 : (bpp > 8) ? 16 : 8;
	    bytesperpixel = bitsperpixel / 8;

	    slicesize = image->getOutputDataSize(); // likewise here

	    buffer.resize(slicesize * file_list.size());
	}

	// use frame 0 (last parameter)
	if(!image->getOutputData(&buffer[bufferoffset], slicesize, 8, 0))
	    std::cerr << "no bytes read at buffer position " << bufferoffset;
	else
	    bufferoffset += (slicesize / bytesperpixel);

	delete(image);

    } // for im

    switch(bitsperpixel) {
	case 64:
	    pixelImport<unsigned long>	(buffer.data(), bufsize, vec);
	    break;
	case 32:
	    pixelImport<unsigned>		(buffer.data(), bufsize, vec);
	    break;
	case 16:
	    pixelImport<unsigned short>	(buffer.data(), bufsize, vec);
	    break;
	case 8:
	    pixelImport<unsigned char>	(buffer.data(), bufsize, vec);
	    break;
    }

    buffer.resize(0);
}

void parseProtocol ( const char* protbuffer , unsigned long protlength, json& sidecar, std::string csafield, std::string protfield ) {
	
	std::string 
		instr, lhs, rhs,
		stringbuffer ( (char*) protbuffer, protlength );
	std::stringstream 
		sbuffer;
	sbuffer.str ( stringbuffer );
	
	while ( instr != "ASCCONV" ) // "### ASCCONV BEGIN ###": start of information block
		sbuffer >> instr;
	sbuffer >> instr; // BEGIN
	sbuffer >> instr; // ###
		
	json protocol;
	while ( instr != "###" ) { // "### ASCCONV END ###": end of information block
	
		std::getline ( sbuffer, lhs, '\n' );
		// std::cout << lhs << "\n";
		
		if ( lhs.find("__attribute__") == std::string::npos ) {			// this means: not found, leaving this in makes structs from arrays
		
			instr = lhs.substr ( 0, lhs.find ( " " ) );
			rhs = trim ( lhs.substr ( lhs.find ( "=" ) +1 ) );
			lhs = trim ( lhs.substr ( 0, lhs.find ( "=" ) ) );
			
			if ( ( instr != "###" ) && ( instr != "" ) ) {
				
				// deal with LHS
				lhs ="." + lhs;											// leading struct
				static const std::regex rgx_pattern("(\\.)[a-z]*([A-Z])");
				lhs = std::regex_replace(lhs, rgx_pattern, "$1$2");		// remove lowercase letters before field
				std::replace(lhs.begin(), lhs.end(), '.', '/');			// turn dot to path -> nlohmann makes struct
				std::replace(lhs.begin(), lhs.end(), '[', '/');			// turn bracket to path -> nlohmann makes array
				std::replace(lhs.begin(), lhs.end(), ']', '/');			// and numbers -> indices, accounting for gaps
				lhs = std::regex_replace(lhs, std::regex("//"), "/");	// empty index/struct -> core dump
				lhs.erase ( lhs.find_last_not_of ( "/" ) +1 );			// remove leading lowercase			
				json::json_pointer key { lhs };
				
				// deal with RHS
				errno = 0; 
				char *test;
				auto ftest = strtof ( rhs.c_str(), &test );

				if ( !errno )
					if ( ftest == floor ( ftest ) )
						protocol [ key ] = static_cast<int> ( ftest );
					else
						protocol [ key ] = static_cast<int> ( roundf ( ftest * 1000000. ) ) / 1000000.;
				else
						protocol [ key ] = trim ( rhs );

				//if ( lhs.find ( "SliceAcqOrder" ) != std::string::npos )
				//	std::cout << std::setw(4) << protocol [ "SliceArray" ] << "lhs: " << lhs << ", instr: " << instr << ", rhs: " << rhs << ", ftest: " << ftest << "\n";
				
			} // if instr not empty
			
		} // if line not about array attributes
		
	} // while instr not ###
	
	if ( instr == "###" ) // that should be on the last line
		sidecar [ csafield ] [ protfield ] = protocol;
		//std::cout << "Siemens CSA parsed successfully\n";
	else
		std::cout << "Siemens CSA: left with instr " << instr << "\n";
	
	return;

}

void parseCSA ( const unsigned char* csabuffer , unsigned long csalength, json& sidecar, std::string jsonfield ) {
	
	// based on https://github.com/neurodebian/spm12/blob/master/spm_dicom_headers.m#L383
	std::string 
		stringbuffer ( (char*) csabuffer, csalength );
	std::stringstream 
		sbuffer;
	sbuffer.str ( stringbuffer );
	
	unsigned
		csa_ver = 1;
	char
		csav [ 4+1 ] = { 0, 0, 0, 0 };
	sbuffer.read ( csav, 4 );
	if ( std::string ( csav ) == "SV10" ) { // csa2 -> read even further
		csa_ver = 2;		
		sbuffer.read ( csav, 4 );			// read 4 more bytes (should be 4b,3b,2b,1b)
	} else 									// csa1 -> go back to start
		sbuffer.seekg ( 0, sbuffer.beg );		
	sidecar [ jsonfield ] [ "CSA_headertype" ] = csa_ver;	
	// std::cout << "0. determined version, position " << sbuffer.tellg() << "\n";
	
	uint32_t 
		num_tags;
	sbuffer.read ( (char*) &num_tags, 4 ); 
	// std::cout << "1. determined #tags(" << num_tags << "), position " << sbuffer.tellg() << "\n";
	
	uint32_t
		them;
	if ( ( num_tags < 1 ) || ( num_tags > 128 ) ) {
		std::cerr << "number of CSA elements below 1 or above 128\n";
		return;
	} else {
		sbuffer.read ( (char*) &them, 4 ); 
		if ( them != 'M' ) { // if next CSA character is 'M'
			std::cerr << "CSA did not dontain the expected 'M'\n";
			return;
		} // if 'M'
	} // if num_tags
	// std::cout << "2. read 'M', position " << sbuffer.tellg() << "\n";
	
	char 
		tagname[64+1], // character array to read each tag's name
		vr[4+1],       // character array to read tag representation
		*val;          // pointer to store (variable-length) values
	uint32_t 
		ni0 = 0;       // number of items of the 1st tag
		
	// process each tag 0 -- the name of the CSA header field
	for ( uint32_t t = 0; t < num_tags; t++ ) {
		
		// we need to read at least 100 bytes ( 64 + 20 + 16 ) to the next value
		if ( ( t < ( num_tags - 1 ) ) && ( 100 + ( long unsigned ) ( sbuffer.tellg() ) ) >= csalength ) {
			std::cout << "remaining CSA block too small";
			break;
		}
				
		// read 64 bytes -> tag name
		tagname [ 64 ] = 0;
		vr      [  4 ] = 0;
		sbuffer.read ( (char*) &tagname, 64 );
		std::string tagstr ( tagname );
		// std::cout << "\t   tag " << t << " name: " << tagstr << "\n";
		// std::cout << "\t3. tag name (64) read, position " << sbuffer.tellg() << "\n";
		
		// read 20 bytes
		uint32_t 
			vm, dt, ni, mm;
		sbuffer.read ( (char*) &vm, 4 ); // value multiplicity
		sbuffer.read ( (char*) &vr, 4 ); // value representation
		std::string tagrep ( vr, 3 );
		sbuffer.read ( (char*) &dt, 4 ); // value representation in syngo
		sbuffer.read ( (char*) &ni, 4 ); // number of items
		sbuffer.read ( (char*) &mm, 4 ); // check: should be 77 (ASCI for 'M')
		// std::cout << "\t   vm " << vm << ", vr " << vr;
		// std::cout <<    ", dt " << dt << ", ni " << ni << ", M " << mm << "\n";
		// std::cout << "\t4. 5 more 4-byte words (20) read, position " << sbuffer.tellg() << "\n";
			
		// dt 3:   // DS -- decimal string
		// dt 4:   // FD -- floating point double
		// dt 5:   // FL -- floating point (single)
		// dt 6:   // IS -- integer string
		// dt 7:   // SL -- signed long
		// dt 8:   // SS -- signed short
		// dt 9:   // UL -- unsigned long
		// dt 10:  // US -- unsigned short
		bool numvalue = ( ( dt >= 3 ) && ( dt <= 10 ) ) ? true : false;
		
		// used in CSA1 (see below)
		if ( !t )
			ni0 = ni;
		
		// process each item -- the values 
		uint32_t itemlen;
		std::vector<std::string> values;
		std::vector<double> numvalues;

		for ( uint32_t i = 0; i < ni; i++ ) {
			
			// std::cout << "\t\t   ";
			// read 4 int32 (numbers of values in item), store in itemlen
			for ( uint32_t k = 0; k < 4; k++ ) { // apparently first {P itemlen 77 Q} for some P,Q
				sbuffer.read ( (char*) &mm, 4 ); // check: should be 77 (ASCI for 'M') 
				if ( ( k == 0 ) && ( csa_ver == 1 ) )
					itemlen = mm - ni0; // see https://github.com/neurodebian/spm12/blob/master/spm_dicom_headers.m#L438
				if ( ( k == 1 ) && ( csa_ver == 2 ) )
					itemlen = mm;
				// std::cout << "mm " << mm << ", ";
			} //for k			
			// std::cout << "\n\t\t5. 4 item-specific 4-byte words (16) read, itemlen " << itemlen << ", position " << sbuffer.tellg() << "\n";
			
			// if the item size takes us past the end of the CSA block, adjust item length
			if ( sbuffer.tellg() > long ( csalength - itemlen ) ) {
				std::cerr << "tried to read past end of CSA block, changing item size from " << itemlen;
				itemlen = csalength - sbuffer.tellg();
				std::cerr << " to " << itemlen << "\n";
			}
			
			// read the value as a char array then string
			val = new char [ itemlen + 1 ];
			val [ itemlen ] = 0;
			sbuffer.read ( val, itemlen );
			if ( itemlen > 0 ) {
				if ( ( tagstr == "MrProtocol" ) || (tagstr == "MrPhoenixProtocol" ) ) 
					parseProtocol ( val, itemlen, sidecar, jsonfield, tagstr );
				else {
					std::string itemval ( val );
					std::stringstream parsed_item;
					if ( numvalue )
						numvalues.push_back ( std::stod ( trim ( itemval ) ) );
					else
						values.push_back( trim ( itemval ) );
				}
			} // if itemlen
			delete[] val;
			
			// put value in sidecar
			switch ( values.size() ) {
			case 0:
				break;
			case 1:
				if ( numvalue )
					sidecar [ jsonfield ] [ tagname ] = numvalues [ 0 ];
				else
					sidecar [ jsonfield ] [ tagname ] = values [ 0 ];
				break;
			default:
				if ( numvalue )
					sidecar [ jsonfield ] [ tagname ] = numvalues;
				else
					sidecar [ jsonfield ] [ tagname ] = values;
				break;
			} // switch itemlen

			// read some bytes into array 'vr' to align with 4-byte address
			uint32_t remain = ( ( sizeof ( uint32_t ) - ( itemlen % sizeof ( uint32_t) ) ) % sizeof ( uint32_t ) ); 
			sbuffer.read ( vr, remain );
			// std::cout << "now read up until position " << sbuffer.tellg() << "\n";;
						
		} // for i (tag items)
				
	} // for t (tags)
	
	// check
	auto check = true;
	if ( check ) {
		std::ofstream ofile ( "/tmp/" + jsonfield + ".json" );
		ofile << std::setw(4) << sidecar [ jsonfield ] << std::endl;
		ofile.close();
	}
	
	return;
	
}



template <typename value_type>

class bisdicom : public bisbids<value_type> {

    // self (own class) and superclasses
    using self 		 = bisdicom	<value_type>;
    using superclass = bisbids	<value_type>;
    using supersuper = bisnifti	<value_type>;
    using superhyper = bisimage	<value_type>;

  private:
    vec3<float> 
		voxsize, 
		slice_orx, 
		slice_ory, 
		slice_pos,
		slice_norm;
    std::vector<size_t> 
		im_size;
	bool	
		mosaic = false;

    // necessary forward declarations -- function definitions are after class definition for general tidiness
    std::vector<slicestats>		get_seriesslicedata	( std::string dirname						);
    std::vector<std::string> 	get_intensitydata	( std::vector<slicestats>& series_slicedata	);
    void 						get_metadata		( std::vector<std::string>& series_files	);

  public:
    bisdicom(std::string fname) { // constructor that reads from a file

	// increase warning threshold
	// see https://support.dcmtk.org/redmine/projects/dcmtk/wiki/howto_logprogram
	OFLog::configure(OFLogger::ERROR_LOG_LEVEL);

	// this is the file format used by DCMtk
	DcmFileFormat dcmfile;

	if(dcmfile.loadFile(fname.c_str()).good()) {

	    // load the header into memory
	    DcmDataset* const 
			dcmdata = dcmfile.getDataset();

	    // for storing DCMTK strings
	    OFString 
			tmpofstring;

	    if(dcmdata->findAndGetOFString(DCM_SeriesInstanceUID, tmpofstring).good())
		superclass::sidecar["SeriesInstanceUID"] = std::string(tmpofstring.c_str());

	    if(!superclass::sidecar["SeriesInstanceUID"].empty()) {

		// scan directory for files with that series UID
		// and store their data in a vector of slicestats
		const std::string dirname = std::filesystem::path(fname).remove_filename();
		auto seriesslicedata = get_seriesslicedata(dirname); // data for DICOM files that match

		// use the file list to load the intensities in all the slices,
		// according to the scan dimensions and directions 
		// (this order is given via a sorted file list)
		auto seriesfiles = get_intensitydata(seriesslicedata);

		// The voxel data are in now. For BIDS, all that is left are the JSON and the NIfTI header
		get_metadata(seriesfiles);

	    } // series UID not empty

	} // if status.good

    } // constructor that reads from a file

    // We don't want to write DICOM, just BIDS
    // overloaded function from bisbids
    void write(std::string filename = "") {

	superclass::write(filename);
    }

}; // class bisdicom

/** \brief getSeriesslicesdata:	given a directory name and series instance UID
 *								return the file data where all slices are found
 * 								in a std::vector of slicestats (see above)
 *
 */
template <class value_type> std::vector<slicestats> bisdicom<value_type>::get_seriesslicedata(std::string dirname) {

	// store image file names / slices info
	std::vector<std::string> valid_files;
    std::vector<slicestats>  seriesslicedata;

    // for storing DCMTK strings
    OFString tmpofstring;

    // loop over all files in the directory and add one with matching tag
    for ( const auto& seriesfile : std::filesystem::directory_iterator ( dirname) ) {
		DcmFileFormat tmpfile;
		if ( tmpfile.loadFile ( seriesfile.path().c_str() ).good() )
			if ( tmpfile.getDataset() -> findAndGetOFString ( DCM_SeriesInstanceUID, tmpofstring).good() )
				if ( std::string ( tmpofstring.c_str() ) == superclass::sidecar [ "SeriesInstanceUID" ] )
					valid_files.push_back ( seriesfile.path() );
    }

    unsigned 
		added_slices = 0,
		added_files  = 0;					// count added files, use 1st file for dimensions
    auto 
		file_iter = valid_files.begin();	// iterator for files in the list with correct UID

    // loop over the list of dicom files with the right series UID
    while(file_iter != valid_files.end()) {
		
	    seriesslicedata.push_back ( { .filename = ( *file_iter ) } );
		
		DcmFileFormat tmpfile;
		tmpfile.loadFile ( ( *file_iter ).c_str() );
		DcmDataset* const tmpdata = tmpfile.getDataset();

		if(!added_slices) { // use first file to populate the info record

			// general info, image size
			superclass::sidecar["Modality"]   = getItemString (                      tmpdata, DCM_Modality             );
			superclass::sidecar["BitsStored"] = atoi          ( getItemStringArray ( tmpdata, DCM_BitsStored ).c_str() );
			int 
				number_fread = atoi(getItemStringArray(tmpdata, DCM_NumberOfTemporalPositions).c_str());
			if (number_fread < 1 ) // sometimes this is omitted if there is just 1 timepoint?
				number_fread = 1;
			superclass::sidecar["NumberOfTemporalPositions"] = number_fread;

			// size of a / the slice in X (rows) and Y (columns)
			im_size.push_back( atoi ( getItemStringArray ( tmpdata, DCM_Rows).c_str()    ) );
			im_size.push_back( atoi ( getItemStringArray ( tmpdata, DCM_Columns).c_str() ) );

			// voxel size
			voxsize[0] =
				static_cast<float>(std::stof(std::string(getItemString(tmpdata, DCM_PixelSpacing, 0).c_str())));
			voxsize[1] =
				static_cast<float>(std::stof(std::string(getItemString(tmpdata, DCM_PixelSpacing, 1).c_str())));
			// Siemens / Philips voxel size: SliceThickness
			voxsize[2] =
				static_cast<float>(std::stof(std::string(getItemString(tmpdata, DCM_SliceThickness).c_str())));
			if ( voxsize[2] < 0 )
				// Philips sometimes
				voxsize[2] =
					static_cast<float>(std::stof(std::string(getItemString(tmpdata, DCM_SpacingBetweenSlices).c_str())));

			if (	( voxsize[0]	< 0 ) ||
					( voxsize[1]	< 0 ) ||
					( voxsize[2]	< 0 )	) {				

				std::cerr << "missing data for file \n " << ( *file_iter ) << ",\ncannot read required info\n";
				break;
				
			}

			// store the voxel size in the sidecar
			superclass::sidecar["VoxelSize"] = {voxsize[0], voxsize[1], voxsize[2]};

			// orientation of the slice plane (orx, ory: perpendicular sides of 1 slice in the stack)
			slice_orx[0] = static_cast<float>(
				std::stof(std::string(getItemString(tmpdata, DCM_ImageOrientationPatient, 0).c_str())));
			slice_orx[1] = static_cast<float>(
				std::stof(std::string(getItemString(tmpdata, DCM_ImageOrientationPatient, 1).c_str())));
			slice_orx[2] = static_cast<float>(
				std::stof(std::string(getItemString(tmpdata, DCM_ImageOrientationPatient, 2).c_str())));
			slice_ory[0] = static_cast<float>(
				std::stof(std::string(getItemString(tmpdata, DCM_ImageOrientationPatient, 3).c_str())));
			slice_ory[1] = static_cast<float>(
				std::stof(std::string(getItemString(tmpdata, DCM_ImageOrientationPatient, 4).c_str())));
			slice_ory[2] = static_cast<float>(
				std::stof(std::string(getItemString(tmpdata, DCM_ImageOrientationPatient, 5).c_str())));
			slice_norm = slice_orx ^ slice_ory; // slice normal ( ^ gives cross product -- see bisimage_maths.hpp )

			// store the orientations in the sidecar
			superclass::sidecar["ImageOrientationPatientX"] = {  slice_orx[0],  slice_orx[1],  slice_orx[2] };
			superclass::sidecar["ImageOrientationPatientY"] = {  slice_ory[0],  slice_ory[1],  slice_ory[2] };
			superclass::sidecar["ImageOrientationPatientZ"] = { slice_norm[0], slice_norm[1], slice_norm[2] };

			// determine slice direction -- see https://stackoverflow.com/questions/34782409
			if(fabs(slice_orx[0]) > fabs(slice_orx[1])) {
				if(fabs(slice_ory[1]) > fabs(slice_ory[2]))
					superclass::sidecar["SlicesDimension"] = 2;	// transverse
				else
					superclass::sidecar["SlicesDimension"] = 1;	// coronal
			} else
				superclass::sidecar["SlicesDimension"] = 0; 	// sagittal

			// determine phase encoding direction
			// assume column direction ( 1 in slice ), i.e. 1 if slice direction == 2, 2 otherwise
			// this covers the "COL" and "OTHER" case for InPlanePhaseEncodingDirection
			superclass::sidecar["PhaseEncodingDimension"] =
				(superclass::sidecar["SlicesDimension"] > 1) ? 1 : 2;
			tmpdata->findAndGetOFString(DCM_InPlanePhaseEncodingDirection, tmpofstring);
			if ( tmpofstring[0] == 'R') // if "ROW" direction ( 0 in slice ), i.e. 0 if slice direction != 0, 1 otherwise
				superclass::sidecar["PhaseEncodingDimension"] =
					(superclass::sidecar["SlicesDimension"] > 0) ? 0 : 1;

			// for now just setting frequency encoding dimension
			// as the one that is left after slice and phase
			superclass::sidecar["FrequencyEncodingDimension"] =
				(superclass::sidecar["SlicesDimension"] == 2) ?
					((superclass::sidecar["PhaseEncodingDimension"] == 0) ? 1 : 0) : // transverse: either 0 or 1
					(superclass::sidecar["SlicesDimension"] == 1) ?
					((superclass::sidecar["PhaseEncodingDimension"] == 0) ? 2 : 0) :  //    coronal: either 0 or 2
						(superclass::sidecar["PhaseEncodingDimension"] == 1) ? 2 : 1; //   sagittal: either 1 or 2
						
			// if this is a multislice image (mosaic, advanced dicom or other) 
			// store that in the JSON for later use: geometry adjustments etc.
			DcmElement 
				*tmpelem;
			tmpdata->findAndGetElement(DCM_ImageType, tmpelem, OFTrue, OFTrue);
			tmpelem->getOFStringArray(tmpofstring);
			std::string
				tmpstring;
			std::stringstream
				sstr ( tmpofstring.c_str() );
			std::vector<std::string> 
				image_types;
			while ( std::getline ( sstr, tmpstring, '\\' ) )
				image_types.push_back ( tmpstring );
			for ( auto s: image_types ) 
				if ( s == "MOSAIC" )
					mosaic = true;

		} // ! added_slices


		// If we have a Siemens mosaic file, read multiple slices from 1 file
		// See https://github.com/InsightSoftwareConsortium/DCMTK/blob/master/dcmdata/data/private.dic
		//     (this may have the number of slices in the string DCM_SourceImageSequence)
		//  or https://dicom-parser.readthedocs.io/en/latest/siemens/private_tags.html
		//     https://github.com/malaterre/GDCM/blob/master/Source/DataDictionary/privatedicts.xml
		size_t
			slicesinfile = 1;
		std::vector<double> 
			slicetimes;

		// Siemens mosaic requires special treatment
		if ( mosaic ) {

			const uint8_t
				*tmpbuffer;
			const double
				*tmpbuf_fd;
			unsigned long 
				csalength;

			slicesinfile = atoi ( getItemStringArray ( tmpdata, DcmTagKey ( 0x0019, 0x100a ) ).c_str() );
			// std::cerr << "found mosaic with " << slicesinfile << " slices\n";
			tmpdata->findAndGetUint8Array  ( DcmTagKey ( 0x0029, 0x1010 ), tmpbuffer, &csalength, OFTrue ); // image data
			parseCSA ( tmpbuffer, csalength, superclass::sidecar, "imageCSA" );			
			tmpdata->findAndGetUint8Array  ( DcmTagKey ( 0x0029, 0x1020 ), tmpbuffer, &csalength, OFTrue ); // series data
			parseCSA ( tmpbuffer, csalength, superclass::sidecar, "seriesCSA" );
			tmpdata->findAndGetFloat64Array( DcmTagKey ( 0x0019, 0x1029 ), tmpbuf_fd, &csalength, OFTrue ); // MosaicRefAcqTimes
			if ( csalength >= slicesinfile )
				for ( unsigned i = 0; i < csalength; i++ )
					slicetimes.push_back ( tmpbuf_fd [ i ] );      // put the times in the vector
					
			auto
				mos_rows = std::ceil ( std::sqrt ( slicesinfile ) ),
				mos_cols = std::ceil ( slicesinfile / mos_rows );
			im_size[0] /= mos_rows;		
			im_size[1] /= mos_cols;

		} // if mosaic

		// now we know the X size (rows) and Y size (columns) of a slice, even for mosaic
		superclass::sidecar [ "Rows"    ] = im_size[0];
		superclass::sidecar [ "Columns" ] = im_size[1];

		for ( size_t s = 0; s<slicesinfile; s++ ) {
			
			slicestats *currentslice = ( &(seriesslicedata[added_slices]) );
			(*currentslice).index = s; 
			
			// get instance number of single file / slice 0
			(*currentslice).innumber =
				static_cast<unsigned>(std::stof(getItemString(tmpdata, DCM_InstanceNumber, 0, 0).c_str()));

														/////////////////////////////////////////////////////////////
			if ( !mosaic ) {							// if no mosaic, fill from header
				
				// get acquisition number of single file / slice 0
				(*currentslice).acqnumber =
					static_cast<unsigned>(std::stof(getItemString(tmpdata, DCM_AcquisitionNumber, 0, 0).c_str()));
				
				// get acquisition time of single file / slice 0
				(*currentslice).aqtime = 
					static_cast<float>(std::stof(getItemString(tmpdata, DCM_AcquisitionTime).c_str()));
				// get acquisition time of single file / slice 0
				(*currentslice).sliloc = 
					static_cast<float>(std::stof(getItemString(tmpdata, DCM_SliceLocation).c_str()));				

				// position / offset from start of the slice in the current volume
				slice_pos[0] =
				static_cast<float>(std::stof(std::string(getItemString(tmpdata, DCM_ImagePositionPatient, 0).c_str())));
				slice_pos[1] =
				static_cast<float>(std::stof(std::string(getItemString(tmpdata, DCM_ImagePositionPatient, 1).c_str())));
				slice_pos[2] =
				static_cast<float>(std::stof(std::string(getItemString(tmpdata, DCM_ImagePositionPatient, 2).c_str())));
															
			} else { 
				
				(*currentslice).aqtime = static_cast<float>(std::stof(getItemString(tmpdata, DCM_AcquisitionTime).c_str()));
														/////////////////////////////////////////////////////////
				if (!s) {								// if mosaic -> set values for slice 0 (not in header)		
					(*currentslice).acqnumber = 0;
					(*currentslice).sliloc    = 0;					
				} else {								/////////////////////////////////////////////////////////
														// for slices 1 and up these are in the header
					(*currentslice).acqnumber = superclass::sidecar [ "seriesCSA"  ] [ "MrPhoenixProtocol" ] [ "SliceArray" ] [ "SliceAcqOrder"     ] [ s ];				
					(*currentslice).sliloc    = superclass::sidecar [ "seriesCSA"  ] [ "MrPhoenixProtocol" ] [ "SliceArray" ] [ "Pos"               ] [ s ];					
					(*currentslice).aqtime   += slicetimes [ s-1 ];
				}
				
				// we assume thet innumber[0] is actually the mosaic number (copy for all slices)
				// and that acqnumber[0] is as well -> set to 0 // see the s==0 section
				//(*currentslice).innumber= superclass::sidecar [ "seriesCSA"  ] [ "MrPhoenixProtocol" ] [ "SliceArray" ] [ "Asc"               ] [ s ];			
				slice_pos[0] = superclass::sidecar [ "seriesCSA"  ] [ "MrPhoenixProtocol" ] [ "SliceArray" ] [ "Slice" ] [ s ] [ "Position" ] [ "Sag" ]; 
				slice_pos[1] = superclass::sidecar [ "seriesCSA"  ] [ "MrPhoenixProtocol" ] [ "SliceArray" ] [ "Slice" ] [ s ] [ "Position" ] [ "Cor" ]; 
				slice_pos[2] = superclass::sidecar [ "seriesCSA"  ] [ "MrPhoenixProtocol" ] [ "SliceArray" ] [ "Slice" ] [ s ] [ "Position" ] [ "Tra" ]; 
			
			} // if mosaic
			
			// compute the z position as projected along the normal
			(*currentslice).slipos =
				(slice_pos & slice_norm); // for type vec3 (see bisimage_maths.hpp) the & operator gives dot product

			// add copy of current slice for next - 
			// we only do this if our file has more than one slices 
			// (i.e. inside the 'for s' loop) to have a new 
			if ( s < ( slicesinfile - 1 ) ) {
				slicestats 
					last_copy {
						seriesslicedata.back().filename,		//	string	filename, for single-slice dicom every slice has a unique file name
						seriesslicedata.back().innumber,			//	int	instance number: should be from 1 to #files
						seriesslicedata.back().acqnumber,			//	int	acquisition number
						seriesslicedata.back().aqtime,				//	float	acquisition time
						seriesslicedata.back().slipos,				//	float	slice position, computed from ImagePOsitionPatient and ImageOrientationPatient
						seriesslicedata.back().sliloc,				//	float 	slice location, backup method sometimes (eg Philips DTI) using SliceLocation
						seriesslicedata.back().index				//	int	index, for multi-slice dicoms (e.g. mosaic), this number can be higher than 0
					}; // last_copy 
				seriesslicedata.push_back ( last_copy );				
				added_slices++; // max value: 1 less than #slice in file
			} // slicesinfile


			
		} // for s slicesinfile

		// advance file list iterator
		added_files++;
		file_iter++;


    } // while file_iter

	if ( ! mosaic ) { // mosaic orders slices in its own special way

		// sort by instance number into subsequent frames caution
		// sometimes (e.g. in ASL) one acquisition is multiple frames
		std::sort(begin(seriesslicedata), end(seriesslicedata),
				  [](auto const& t0, auto const& t1) { return ((t0).innumber < (t1).innumber); });
				  
		// Check if slice 0 and 1 have the same slice position. This is true in e.g. Philips fMRI
		// (see https://github.com/rordenlab/dcm2niix/tree/master/Philips), stable sort by aq time.
		// Stable sort respects previous ordering where possible.
		if ( seriesslicedata[0].slipos == seriesslicedata[1].slipos )
			std::stable_sort(begin(seriesslicedata), end(seriesslicedata),
							 [](auto const& t0, auto const& t1) { return ((t0).aqtime < (t1).aqtime); });

		// Check if slice 0 and 1 still have the same slice position. This happens in e.g. Philips DWI
		// (see https://github.com/rordenlab/dcm2niix/tree/master/Philips), stable sort by SliceLocation.
		// Stable sort respects previous ordering for ties.
		if ( seriesslicedata[0].slipos == seriesslicedata[1].slipos ) {
			char fac = ( seriesslicedata[0].slipos < seriesslicedata[added_slices - 1].slipos ) ? 1 : -1;
			std::stable_sort( begin(seriesslicedata), end(seriesslicedata),
							  [&](auto const& t0, auto const& t1) { return ( fac*(t0).slipos < fac*(t1).slipos ) ; } );
			auto lowest_innumber = seriesslicedata[0].innumber;
			auto first_slipos = std::count_if ( begin(seriesslicedata), end(seriesslicedata), 
												[&](auto const& fd0) { return ( (fd0).sliloc == seriesslicedata[0].sliloc ); } ) ;
			std::stable_sort( begin(seriesslicedata), end(seriesslicedata),
							  [&](auto const& t0, auto const& t1) { return ( ((t0).innumber-lowest_innumber)%first_slipos < ((t1).innumber-lowest_innumber)%first_slipos ); } );
		}
	
	} else { // mosaic -> sort by innumber then sliloc
	
			std::stable_sort( begin(seriesslicedata), end(seriesslicedata),
							  [&](auto const& t0, auto const& t1) { return ( ( (t0).innumber < (t1).innumber ) || ( ( (t0).innumber == (t1).innumber ) && ( (t0).acqnumber < (t1).acqnumber ) ) ); } );

	}
	
    auto diag = false;
    if (diag) {
		
		// print diagnostics -- for weird slice orderings like
		//							Philips fMRI   -> slice 0 of every volume, slice 1 of every volume, etc
		//                          Siemens Mosaic -> multiple slices turned into patchwork

		for(auto slicedata : seriesslicedata)
			std::cout << std::setprecision(6) 
					<< slicedata.filename << "\t" << slicedata.innumber << "\t" << slicedata.acqnumber << "\t"
					<< slicedata.aqtime   << "\t" << slicedata.slipos   << "\t" << slicedata.sliloc   << std::endl;
				  
    }

    return (seriesslicedata);
	
}

/** \brief get_intensitydata: 	given a std::vector of slicestats (see above)
 * 								process this list to determin the intensity
 * 								data's dimensions and read them in properly.
 *
 */
template <class value_type>
std::vector<std::string> bisdicom<value_type>::get_intensitydata(std::vector<slicestats>& series_slicedata) {

    auto 
		num_slices  = series_slicedata.size(), // number of correct slices with correct UID
		num_files   = num_slices;
	if ( mosaic ) 
		num_files /= static_cast<int> ( superclass::sidecar [ "seriesCSA"  ] [ "MrPhoenixProtocol" ] [ "SliceArray" ] [ "Size" ] );

    // Determine the #slices in a volume by checking
    // when the slice position is back at its initial value.
    unsigned sli = 1;
    for(auto startpos = (series_slicedata[0]).sliloc; sli < series_slicedata.size(); sli++)
	if((series_slicedata[sli]).sliloc == startpos)
	    break;
    superclass::sidecar["Slices"] = sli;

    std::cout << "counted slices (before returning to pos 0): " << sli << std::endl;

    // double check if this is still true at the end
    sli = series_slicedata.size() - 2;
    for(auto startpos = (series_slicedata[num_slices - 1]).sliloc; sli > 0; sli--)
	if((series_slicedata[sli]).sliloc == startpos) { 
		sli++;
	    break;
	}
    sli = num_slices - sli; // counting from n-1 down to 0

    std::cout << "counted slices (counting backwards from end): " << sli << std::endl;

    // In the first volume, find the sign of the slice positioning: if the
    // position of slice 1 is lower than slice 0 then -1, otherwise +1.
    // THIS NEEDS TO CORRESPOND TO THE SLICE NORMAL for determining rotate 180 or flip
	float 
		slice_delta = series_slicedata[1].slipos - series_slicedata[0].slipos;
    superclass::sidecar["Slicedirection"] = static_cast<int>(bis::signum<float>( slice_delta ));
	
    // in the case of 2D files, #volumes = #files / (slices per volume)
    // TR of a volume = start time of volume 1 - start time of volume 0
    // (we know sometimes images are acquired in pairs - cater for that)
    superclass::sidecar["Volumes"] = num_slices / int(superclass::sidecar["Slices"]);
    superclass::sidecar["RepetitionTime"] = ((series_slicedata[2 * sli]).aqtime - (series_slicedata[0]).aqtime) / 2.;

    // std::cout << superclass::sidecar << std::endl;

    // continue if the first and last volume have the same #slices
    // and if the #files is the product of #slices/vol and #volumes
    auto vol = unsigned(superclass::sidecar["Volumes"]);
    if( ( unsigned ( superclass::sidecar [ "Slices" ] ) == sli ) && ( num_slices == ( sli * vol ) ) ) {

	// output size at least 2d, check if 3d or 4d and if yes add those dimensions
	if ( sli > 1 ) {
	    im_size.push_back(sli);

	    // if one acquisition has more than one volume, change acquisition index to volume index
	    // (if acquisition 0 is 4 slices in 2 volumes, change from { 0,0,0,0,.. } to { 0,0,1,1,.. } )
	    if(vol > 1) {

		// check if sizes based on header data and pixel data are the same
		if(vol != (num_slices / sli)) {
		    std::cerr << "DICOM tag for #vol different than #sli / #vol";
		    vol = num_slices / sli; // only workable formula for pixel data
		}
		im_size.push_back(vol);

	    } // if vol > 1

	} // if sli > 1

    } // #files is the product of #vol * *#(slices/vol)

    // re-order slice aquisition sequence
    for(unsigned i = 0; i < num_slices; i++)
		( series_slicedata[i] ).acqnumber = i / int ( superclass::sidecar [ "Slices" ] );

    // resize the data array in bisimage for putting the voxels in
    superhyper::newsizes( { im_size[0], im_size[1], sli, vol } ); // first 2 always exist, 3 and 4 don't

    // sort by aqcuisition number and if that is the same (multiple vol/timepoint), slice
    // position: sorting by instance number is not working in the case of interleaved stacks
    std::sort(begin(series_slicedata), end(series_slicedata), [](auto const& t0, auto const& t1) {
	return (((t0).acqnumber < (t1).acqnumber) || (((t0).acqnumber == (t1).acqnumber) && ((t0).slipos < (t1).slipos)));
    });

    // get a list of only the files, e.g. extracting from the structs
    // see https://www.fluentcpp.com/2018/11/09/retrieve-firsts-collection-pairs for how this works for tuples
    std::vector<std::string> seriesfiles;
    std::transform(begin(series_slicedata), end(series_slicedata), std::back_inserter(seriesfiles),
                   [](auto const& fdata) { return (fdata).filename; });

    // check bisimage size
    std::cout << "data sizes: " << superhyper::getsize() << std::endl;

    // read in the slices -- don't go and look, dirty bit of code...
    getDicomData(seriesfiles, superhyper::data.size(), superhyper::getdata_ptr());

    return (seriesfiles);
}

/** \brief get_metadata: 	given a std::vector of filenames (see above)
 *							use the values in the headers of these files to
 * 							complete the JSON sidecar and the NIfTI header.
 */
template <class value_type> void bisdicom<value_type>::get_metadata(std::vector<std::string>& series_files) {

    // the locals
    DcmFileFormat tmpfile;
    vec3<float> sli_0_pos, sli_n_pos;
    auto sli = unsigned(superclass::sidecar["Slices"]) - 1;

    // get the slice position from volume 0 slice 0
    tmpfile.loadFile(series_files[0].c_str());
    DcmDataset* tmpdata = tmpfile.getDataset();
    sli_0_pos[0] =
        static_cast<float>(std::stof(std::string(getItemString(tmpdata, DCM_ImagePositionPatient, 0).c_str())));
    sli_0_pos[1] =
        static_cast<float>(std::stof(std::string(getItemString(tmpdata, DCM_ImagePositionPatient, 1).c_str())));
    sli_0_pos[2] =
        static_cast<float>(std::stof(std::string(getItemString(tmpdata, DCM_ImagePositionPatient, 2).c_str())));

    // get the slice position from volume 0 slice N
    tmpfile.loadFile(series_files[sli].c_str());
    tmpdata = tmpfile.getDataset();
    sli_n_pos[0] =
        static_cast<float>(std::stof(std::string(getItemString(tmpdata, DCM_ImagePositionPatient, 0).c_str())));
    sli_n_pos[1] =
        static_cast<float>(std::stof(std::string(getItemString(tmpdata, DCM_ImagePositionPatient, 1).c_str())));
    sli_n_pos[2] =
        static_cast<float>(std::stof(std::string(getItemString(tmpdata, DCM_ImagePositionPatient, 2).c_str())));

    // Dicom voxel to mm transform
    mat4<float> Rdcm(slice_orx[0] * voxsize[0], slice_ory[0] * voxsize[1], (sli_n_pos[0] - sli_0_pos[0]) / sli, sli_0_pos[0], 
					 slice_orx[1] * voxsize[0], slice_ory[1] * voxsize[1], (sli_n_pos[1] - sli_0_pos[1]) / sli, sli_0_pos[1], 
					 slice_orx[2] * voxsize[0], slice_ory[2] * voxsize[1], (sli_n_pos[2] - sli_0_pos[2]) / sli, sli_0_pos[2], 
					 0, 0, 0, 1);

    // NIfTI voxel to mm transform
    mat4<float> Rnii(	-Rdcm[0][0], -Rdcm[0][1], -Rdcm[0][2], -Rdcm[0][3], 
						-Rdcm[1][0], -Rdcm[1][1], -Rdcm[1][2], -Rdcm[1][3], 
						 Rdcm[2][0],  Rdcm[2][1],  Rdcm[2][2],  Rdcm[2][3], 
						 0, 0, 0, 1);

    // usually (if slices go top -> bottom) we want to flip the Y axis
    //                  (dicom: positive is RL, nifti: positive is LR)
    // see https://github.com/rordenlab/dcm2niix/blob/master/console/nii_dicom.cpp#L2659
    if(superclass::sidecar["Slicedirection"] < 0) {

	Rnii[0][1] *= -1; //
	Rnii[1][1] *= -1; // flip the 2nd (y) column of Rnii
	Rnii[2][1] *= -1; //

	vec4<float> yflip_ori(0, im_size[1] - 1, 0, 0); // yflip_ori: vector of the new origin
	yflip_ori = (Rnii * yflip_ori);

	Rnii[0][3] -= yflip_ori[0]; //
	Rnii[1][3] -= yflip_ori[1]; // restore origin after flipping
	Rnii[2][3] -= yflip_ori[2]; //

	// and flip the voxels!
	size_t st = (superhyper::sizes.size() > 3) ? superhyper::sizes[3] : 1,
	       sz = (superhyper::sizes.size() > 2) ? superhyper::sizes[2] : 1,
	       sy = (superhyper::sizes.size() > 1) ? superhyper::sizes[1] : 1,
	       sx = (superhyper::sizes.size() > 0) ? superhyper::sizes[0] : 1, sy1 = sy - 1, sy2 = sy / 2;
	for(size_t t = 0; t < st; t++)
	    for(size_t k = 0; k < sz; k++)
		for(size_t j = 0; j < sy2; j++)
		    for(size_t i = 0; i < sx; i++)
			std::swap<value_type>( superhyper::at ( { i, j, k , t } ), superhyper::at ( { i, sy1 - j, k, t } ) );
			// std::swap<value_type>( (*this)[{i, j, k, t}], (*this)[{i, sy1 - j, k, t}]);
			

    }; // if slicedirection

    //
    // populate the NIfTI header
    //

    // set the own (nifti type) header
    supersuper::header = nifti_simple_init_nim();
    supersuper::header->datatype = getniitype<value_type>();
    nifti_datatype_sizes(supersuper::header->datatype, &(supersuper::header->nbyper), &(supersuper::header->swapsize));
    auto dicompath = std::filesystem::path(series_files[0]);
    nifti_set_filenames(supersuper::header, dicompath.replace_extension("nii.gz").c_str(), 0, 0);

    // set dimensions and spacings
    supersuper::header->dim[0] = superhyper::sizes.size();
    supersuper::header->pixdim[0] = superhyper::sizes.size();
    if(supersuper::header->dim[0]) {
	supersuper::header->dim[1] = superhyper::sizes[0];
	supersuper::header->pixdim[1] = voxsize[0];
	if(supersuper::header->dim[0] > 1) {
	    supersuper::header->dim[2] = superhyper::sizes[1];
	    supersuper::header->pixdim[2] = voxsize[1];
	    if(supersuper::header->dim[0] > 2) {
		supersuper::header->dim[3] = superhyper::sizes[2];
		supersuper::header->pixdim[3] = voxsize[2];
		if(supersuper::header->dim[0] > 3) {
		    supersuper::header->dim[4] = superhyper::sizes[3];
		    supersuper::header->pixdim[4] = float(superclass::sidecar["RepetitionTime"]);
		    // never encountered 5D but you can go on of course
		} // if >3
	    }     // if >2
	}         // if >1
    }             // if >0

    // get rescale slope and intercept - first try DCM_RescaleSlope/Intercept
    //        and if that does not work try DCM_RealWorldValueSlope/Intercept
    supersuper::header->scl_slope = getItemDouble(tmpdata, DCM_RescaleSlope);
    if(supersuper::header->scl_slope)
	supersuper::header->scl_slope = getItemDouble(tmpdata, DCM_RescaleIntercept);
    else {
	supersuper::header->scl_slope = getItemDouble(tmpdata, DCM_RealWorldValueSlope);
	supersuper::header->scl_inter = getItemDouble(tmpdata, DCM_RealWorldValueIntercept);
    } // if scl_slope

    // set cal_min and cal_max as highest and lowest nonzero values
    auto ordered_nonzero = std::vector<value_type>(0);
    for(auto intensity = superhyper::data.begin(); intensity != superhyper::data.end(); intensity++)
	if(*intensity != 0)
	    ordered_nonzero.push_back(*intensity);
    std::sort(ordered_nonzero.begin(), ordered_nonzero.end());
    supersuper::header->cal_min = float(ordered_nonzero.front());
    supersuper::header->cal_max = float(ordered_nonzero.back());
    ordered_nonzero.resize(0); // does that clear the memory?

    // assume mm as vox units, s as time units
    supersuper::header->xyz_units = NIFTI_UNITS_MM;
    supersuper::header->time_units =
        (float(superclass::sidecar["RepetitionTime"])) < 3600 ? NIFTI_UNITS_SEC : NIFTI_UNITS_MSEC;

    // set slice, phase and frequency encoding directions
    // supersuper::header -> dim_info   = FPS_INTO_DIM_INFO ( sidecar [ "FrequencyEncodingDimension" ],
    // sidecar [ "PhaseEncodingDimension" ], sidecar [ "SlicesDimension" ] );
    supersuper::header->freq_dim = int(superclass::sidecar["FrequencyEncodingDimension"]) + 1;
    supersuper::header->phase_dim = int(superclass::sidecar["PhaseEncodingDimension"]) + 1;
    supersuper::header->slice_dim = int(superclass::sidecar["SlicesDimension"]) + 1;

    // The S-form can handle more transforms than the Qform,
    // so we copy the Q-form to the S-form (not vice versa)
    //
    // we need some data from vol 0 slice 0, and vol 0 slice n-1, to get qform
    // see https://discovery.ucl.ac.uk/id/eprint/1495621

    // set Qform in NIfTI header ( mat44 is an internal nifti data type )
    mat44 m;
    m.m[0][0] = Rnii[0][0];    m.m[0][1] = Rnii[0][1];    m.m[0][2] = Rnii[0][2];    m.m[0][3] = Rnii[0][3];
    m.m[1][0] = Rnii[1][0];    m.m[1][1] = Rnii[1][1];    m.m[1][2] = Rnii[1][2];    m.m[1][3] = Rnii[1][3];
    m.m[2][0] = Rnii[2][0];    m.m[2][1] = Rnii[2][1];    m.m[2][2] = Rnii[2][2];    m.m[2][3] = Rnii[2][3];
    m.m[3][0] = Rnii[3][0];    m.m[3][1] = Rnii[3][1];    m.m[3][2] = Rnii[3][2];    m.m[3][3] = Rnii[3][3];

    std::cout << "q-form matrix:" << std::endl;
    std::cout << m.m[0][0] << " " << m.m[0][1] << " " << m.m[0][2] << " " << m.m[0][3] << " " << std::endl;
    std::cout << m.m[1][0] << " " << m.m[1][1] << " " << m.m[1][2] << " " << m.m[1][3] << " " << std::endl;
    std::cout << m.m[2][0] << " " << m.m[2][1] << " " << m.m[2][2] << " " << m.m[2][3] << " " << std::endl;
    std::cout << m.m[3][0] << " " << m.m[3][1] << " " << m.m[3][2] << " " << m.m[3][3] << " " << std::endl;

    float qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac = 1;
    nifti_mat44_to_quatern(m, &qb, &qc, &qd, &qx, &qy, &qz, &dx, &dy, &dz, &qfac);
    supersuper::header->qform_code = NIFTI_XFORM_SCANNER_ANAT;
    supersuper::header->quatern_b = qb;
    supersuper::header->quatern_c = qc;
    supersuper::header->quatern_d = qd;
    supersuper::header->qoffset_x = qx;
    supersuper::header->qoffset_y = qy;
    supersuper::header->qoffset_z = qz;
    supersuper::header->pixdim[0] = supersuper::header->qfac = superclass::sidecar["Slicedirection"]; // see nifti1.h

    supersuper::header->sform_code = NIFTI_XFORM_SCANNER_ANAT;
    supersuper::header->sto_xyz.m[0][0] = Rnii[0][0];
    supersuper::header->sto_xyz.m[0][1] = Rnii[0][1];
    supersuper::header->sto_xyz.m[0][2] = Rnii[0][2];
    supersuper::header->sto_xyz.m[0][3] = Rnii[0][3];
    supersuper::header->sto_xyz.m[1][0] = Rnii[1][0];
    supersuper::header->sto_xyz.m[1][1] = Rnii[1][1];
    supersuper::header->sto_xyz.m[1][2] = Rnii[1][2];
    supersuper::header->sto_xyz.m[1][3] = Rnii[1][3];
    supersuper::header->sto_xyz.m[2][0] = Rnii[2][0];
    supersuper::header->sto_xyz.m[2][1] = Rnii[2][1];
    supersuper::header->sto_xyz.m[2][2] = Rnii[2][2];
    supersuper::header->sto_xyz.m[2][3] = Rnii[2][3];
    supersuper::header->sto_xyz.m[3][0] = Rnii[3][0];
    supersuper::header->sto_xyz.m[3][1] = Rnii[3][1];
    supersuper::header->sto_xyz.m[3][2] = Rnii[3][2];
    supersuper::header->sto_xyz.m[3][3] = Rnii[3][3];

    // check the header
    nifti_update_dims_from_array(supersuper::header);

    //
    // complete the JSON sidecar and write
    //
    superclass::sidecar["TotalAcquiredPairs"] = im_size[3] / int(superclass::sidecar["NumberOfTemporalPositions"]);

    return;
}

} // namespace bis

#endif // BISDICOMIMAGE_HPP_INCLUDED