#ifndef BISDICOMIMAGE_HPP_INCLUDED
#define BISDICOMIMAGE_HPP_INCLUDED



/** \brief use the standard library for trawling dicom directories
 */
#include <filesystem>

/** \brief bisniftiimage is a subclass of bisimage
 */
#include "bisniftiimage.hpp"

/** \brief bisdicomimage.hpp uses dcmtk for DICOM file I/O
 *         (basically input), in a single function
 */
#include "dcmtk/config/osconfig.h"
#include "dcmtk/dcmdata/dcistrmf.h"
#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmimgle/dcmimage.h"

/** \brief DICOM headers may contain BIDS data fields that
 *         cannot be stored in NIfTI headers; these are
 *         stored in a JSON 'sidecar' with the same name
 */
#include "nlohmann/json.hpp"
using json = nlohmann::json;



/** \brief bis: namespace for multidimensional images
 *
 *  This namespace contains the bisimage class for
 *  creating and processing multidimensional images.
 */

namespace bis {

	
	
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
    getItemString ( DcmItem* item, const DcmTagKey& tagKey, const size_t position = 0, const OFBool searchsub = OFFalse ) {

        OFString tmpstring;

        if ( item->findAndGetOFString ( tagKey, tmpstring, position, searchsub ).good() )
            return std::string ( tmpstring.c_str() );
        else
            return "";
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
    std::string getItemStringArray ( DcmItem* item, const DcmTagKey& tagKey ) {

        OFString tmpstring;

        if ( item->findAndGetOFStringArray ( tagKey, tmpstring ).good() )
            return std::string ( tmpstring.c_str() );
        else
            return "";
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
    double getItemDouble ( DcmItem* item, const DcmTagKey& tagKey, const size_t position = 0, const OFBool searchsub = OFTrue ) {

        double tmpdouble;

        if ( item->findAndGetFloat64 ( tagKey, tmpdouble, position, searchsub ).good() )
            return ( tmpdouble );
        else
            return 0.;
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

    void getDicomData ( std::vector<std::string> file_list, unsigned bufsize, DataVec* vec ) {

        std::vector<unsigned char> buffer;      // for reading binary data
        size_t slicesize = 0, bufferoffset = 0; // for reading pixels	
        int       bitsperpixel = 8;             // usually changed to 16 (intensities to about 2000)
        unsigned bytesperpixel = 1;             // usually changed to 2  

        for ( size_t im = 0; im < file_list.size(); im++ ) {

            DicomImage* image = new DicomImage ( file_list[im].c_str() );

            // determine some stats in the first image
            if ( !im ) {

                auto bpp = image->getDepth(), // assume this is the same for all
                     bitsperpixel = ( bpp > 32 ) ? 64 : ( bpp > 16 ) ? 32 : ( bpp > 8 ) ? 16 : 8;
                bytesperpixel = bitsperpixel / 8;

                slicesize = image->getOutputDataSize(); // likewise here

                buffer.resize ( slicesize * file_list.size() );
            }

            // use frame 0 (last parameter)
            if ( !image->getOutputData ( &buffer[bufferoffset], slicesize, 8, 0 ) )
                std::cerr << "no bytes read at buffer position " << bufferoffset;
            else
                bufferoffset += ( slicesize / bytesperpixel );

            delete ( image );

        } // for im

        switch ( bitsperpixel ) {
            case 64:
                pixelImport<unsigned long> ( buffer.data(), bufsize, vec );
                break;
            case 32:
                pixelImport<unsigned> ( buffer.data(), bufsize, vec );
                break;
            case 16:
                pixelImport<unsigned short> ( buffer.data(), bufsize, vec );
                break;
            case 8:
                pixelImport<unsigned char> ( buffer.data(), bufsize, vec );
                break;
        }

        buffer.resize ( 0 );
    }



    template <typename value_type>

    class bisdicom : public bisnifti<value_type> {

		// self (own class) and superclasses
		using self       = bisdicom<value_type>;
		using superclass = bisnifti<value_type>;
		using supersuper = bisimage<value_type>;

        private:

            json 
                sidecar;
				
            vec3 <float> 
                voxsize, 
                slice_orx, 
                slice_ory, 
                slice_norm;
				
            std::vector <size_t>
                im_size;

        public:
        
            bisdicom ( std::string fname ) {                    // constructor that reads from a file

                // increase warning threshold
                // see https://support.dcmtk.org/redmine/projects/dcmtk/wiki/howto_logprogram
                OFLog::configure ( OFLogger::ERROR_LOG_LEVEL );

				// this is the file format used by DCMtk
                DcmFileFormat 
                    dcmfile;

                if ( dcmfile.loadFile ( fname.c_str() ).good() ) {

                    // load the header into memory
                    DcmDataset* const 
                        dcmdata = dcmfile.getDataset();

                    // get the series UID to find the files needed for this 2D/3D/4D data set
                    OFString 
                        tmpstring;
                        
                    if ( dcmdata->findAndGetOFString ( DCM_SeriesInstanceUID, tmpstring ).good() )
                        sidecar["SeriesInstanceUID"] = std::string ( tmpstring.c_str() );

                    if ( !sidecar["SeriesInstanceUID"].empty() ) {

                        //
                        // for gathering the slices, which may be in different files
                        //
                        // make a vector of tuples:
                        // < std::string; # filename
                        //   int          # instance number: should be from 1 to #files
                        //   int          # volume number (can be different in 1 acq)
                        //   float        # acquisition time
                        //   float;       # slice position
                        // >                               -- tuples are easy to sort!
                        
						typedef struct filestats {
							std::string filename;
							unsigned    innumber;
							unsigned    aqnumber;
							float		aqtime;
							float		slipos;
						} filestats;
						
                        std::vector < filestats > 
                            seriesdata;						

                        // scan directory for files with that series UID
                        const std::string 
                            dirname = std::filesystem::path ( fname ).remove_filename();

                        // loop over all files in the directory and add one with matching tag
                        for ( const auto& seriesfile : std::filesystem::directory_iterator ( dirname ) ) {                            
                            DcmFileFormat 
                                tmpfile;
                            if ( tmpfile.loadFile ( seriesfile.path().c_str() ).good() )
                                if ( tmpfile.getDataset()->findAndGetOFString ( DCM_SeriesInstanceUID, tmpstring ).good() )
                                    if ( std::string ( tmpstring.c_str() ) == sidecar["SeriesInstanceUID"] )
                                        seriesdata.push_back ( { .filename = seriesfile.path() } );
                        }

                        // concatenate files together
                        unsigned added_files = 0;            // count added files, use 1st file for dimensions
                        auto file_iter = seriesdata.begin(); // iterator for files in the list with correct UID

                        // loop over the list of dicom files with the right series UID
                        while ( file_iter != seriesdata.end() ) {

                            DcmFileFormat tmpfile;
                            tmpfile.loadFile ( ( *file_iter ).filename.c_str() );
                            DcmDataset* const tmpdata = tmpfile.getDataset();

                            if ( !added_files ) { // use first file to populate the info record

                                // general info, image size
                                sidecar [ "Modality" ] = 
                                    getItemString ( tmpdata, DCM_Modality );
                                sidecar [ "BitsStored"                ] = 
                                    atoi ( getItemStringArray ( tmpdata, DCM_BitsStored                ).c_str() );
                                sidecar [ "NumberOfTemporalPositions" ] =
                                    atoi ( getItemStringArray ( tmpdata, DCM_NumberOfTemporalPositions ).c_str() );
                                // size of a / the slice
                                im_size.push_back ( atoi ( getItemStringArray ( tmpdata, DCM_Rows    ).c_str() ) );
                                sidecar [ "Rows" ] = 
                                    im_size.back();
                                im_size.push_back ( atoi ( getItemStringArray ( tmpdata, DCM_Columns ).c_str() ) );
                                sidecar [ "Columns"                   ] = 
                                    im_size.back();

                                // voxel size
                                voxsize[0] = static_cast<float> (
                                                 std::stof ( std::string ( getItemString ( tmpdata, DCM_PixelSpacing, 0 ).c_str() ) ) );
                                voxsize[1] = static_cast<float> (
                                                 std::stof ( std::string ( getItemString ( tmpdata, DCM_PixelSpacing, 1 ).c_str() ) ) );
                                voxsize[2] = static_cast<float> (
                                                 std::stof ( std::string ( getItemString ( tmpdata, DCM_SpacingBetweenSlices ).c_str() ) ) );

                                // store the voxel size in the sidecar
                                sidecar["VoxelSize"] = { voxsize[0], voxsize[1], voxsize[2] };

                                // orientation of the slice plane (orx, ory: perpendicular sides of 1 slice in the stack)
                                slice_orx[0] = static_cast<float> (
                                                   std::stof ( std::string ( getItemString ( tmpdata, DCM_ImageOrientationPatient, 0 ).c_str() ) ) );
                                slice_orx[1] = static_cast<float> (
                                                   std::stof ( std::string ( getItemString ( tmpdata, DCM_ImageOrientationPatient, 1 ).c_str() ) ) );
                                slice_orx[2] = static_cast<float> (
                                                   std::stof ( std::string ( getItemString ( tmpdata, DCM_ImageOrientationPatient, 2 ).c_str() ) ) );
                                slice_ory[0] = static_cast<float> (
                                                   std::stof ( std::string ( getItemString ( tmpdata, DCM_ImageOrientationPatient, 3 ).c_str() ) ) );
                                slice_ory[1] = static_cast<float> (
                                                   std::stof ( std::string ( getItemString ( tmpdata, DCM_ImageOrientationPatient, 4 ).c_str() ) ) );
                                slice_ory[2] = static_cast<float> (
                                                   std::stof ( std::string ( getItemString ( tmpdata, DCM_ImageOrientationPatient, 5 ).c_str() ) ) );
                                slice_norm =
                                    slice_orx ^ slice_ory; // slice normal ( ^ gives cross product -- see bisimage_maths.hpp )

                                // store the orientations in the sidecar
                                sidecar["ImageOrientationPatientX"] = { slice_orx[0], slice_orx[1], slice_orx[2] };
                                sidecar["ImageOrientationPatientY"] = { slice_ory[0], slice_ory[1], slice_ory[2] };
                                sidecar["ImageOrientationPatientZ"] = { slice_norm[0], slice_norm[1], slice_norm[2] };

                                // determine slice direction -- see https://stackoverflow.com/questions/34782409
                                if ( fabs ( slice_orx[0] ) > fabs ( slice_orx[1] ) )
                                    if ( fabs ( slice_ory[1] ) > fabs ( slice_ory[2] ) )
                                        sidecar["SlicesDimension"] = 2; // transverse
                                    else
                                        sidecar["SlicesDimension"] = 1; // coronal
                                else
                                    sidecar["SlicesDimension"] = 0; // sagittal

                                // determine phase encoding direction
                                tmpdata->findAndGetOFString ( DCM_InPlanePhaseEncodingDirection, tmpstring );
                                switch ( tmpstring[0] ) {
                                    case 'R': // row direction    ( 0 in slice ), i.e. 0 if slice direction != 0, 1 otherwise
                                        sidecar["PhaseEncodingDimension"] = ( sidecar["SlicesDimension"] > 0 ) ? 0 : 1;
                                    case 'C': // column direction ( 1 in slice ), i.e. 1 if slice direction == 2, 2 otherwise
                                        sidecar["PhaseEncodingDimension"] = ( sidecar["SlicesDimension"] > 1 ) ? 1 : 2;
                                }

                                // for now just setting frequency encoding dimension
                                // as the one that is left after slice and phase
                                sidecar["FrequencyEncodingDimension"] = ( sidecar["SlicesDimension"] == 2 ) ?
                                                                        ( ( sidecar["PhaseEncodingDimension"] == 0 ) ? 1 : 0 ) : // transverse: either 0 or 1
                                                                        ( sidecar["SlicesDimension"] == 1 ) ?
                                                                        ( ( sidecar["PhaseEncodingDimension"] == 0 ) ? 2 : 0 ) : //    coronal: either 0 or 2
                                                                        ( sidecar["PhaseEncodingDimension"] == 1 ) ? 2 : 1; //   sagittal: either 1 or 2

                            } // ! added_files

                            // get instance number of single file
                            ( *file_iter ).innumber =
                                static_cast<unsigned> ( std::stof ( getItemString ( tmpdata, DCM_InstanceNumber, 0, 0 ).c_str() ) );
                            // get acquisition number of single file
                            ( *file_iter ).aqnumber =
                                static_cast<unsigned> ( std::stof ( getItemString ( tmpdata, DCM_AcquisitionNumber, 0, 0 ).c_str() ) );
                            // get acquisition time of single file
                            ( *file_iter ).aqtime =
                                static_cast<float> ( std::stof ( getItemString ( tmpdata, DCM_AcquisitionTime ).c_str() ) );

                            // position / offset from start of the slice in the current volume
                            vec3<float> slicepos;
                            slicepos[0] = static_cast<float> (
                                              std::stof ( std::string ( getItemString ( tmpdata, DCM_ImagePositionPatient, 0 ).c_str() ) ) );
                            slicepos[1] = static_cast<float> (
                                              std::stof ( std::string ( getItemString ( tmpdata, DCM_ImagePositionPatient, 1 ).c_str() ) ) );
                            slicepos[2] = static_cast<float> (
                                              std::stof ( std::string ( getItemString ( tmpdata, DCM_ImagePositionPatient, 2 ).c_str() ) ) );

                            // get position along slice normal of single file and put in tuple <4> for this file
                            ( *file_iter ).slipos =
                                ( slicepos & slice_norm ); // for type vec3 (see bisimage_maths) the & operator gives dot product

                            // increase number of files added
                            added_files++;

                            // advance file list iterator
                            file_iter++;

                        } // while file_iter

                        // sort by tuple value <1> (instance number) into subsequent frames
						// caution: sometimes (e.g. in ASL) one acquisition is multiple frames
                        std::sort ( begin ( seriesdata ), end ( seriesdata ),
                        [] ( auto const & t0, auto const & t1 ) {
                            return ( ( t0 ).innumber < ( t1 ).innumber );
                        } );

                        // In the first volume, find the sign of the slice positioning: if the
                        // position of slice 1 is lower than slice 0 then -1, otherwise +1.
                        sidecar["Slicedirection"] =
                            static_cast<int> ( bis::signum<float> ( seriesdata[1].slipos - seriesdata[0].slipos ) );

                        // Now determine the #slices in a volume by checking
                        // when the slice position is back at it's initial value.
                        unsigned sli = 1;
                        for ( auto startpos = ( seriesdata[0] ).slipos; sli < added_files; sli++ )
                            if ( ( seriesdata[sli] ).slipos == startpos )
                                break;
                        sidecar["Slices"] = sli;

                        // double check if this is still true at the end
                        sli = added_files - 2;
                        for ( auto startpos = ( seriesdata[added_files - 1] ).slipos; sli >= 0; sli-- )
                            if ( ( seriesdata[sli] ).slipos == startpos )
                                break;
                        sli = added_files - sli - 1; // counting from n-1 down to 0, ya gotta love it

                        // in the case of 2D files, #volumes = #files / (slices per volume)
                        // TR of a volume = start time of volume 1 - start time of volume 0
                        // (we know sometimes images are acquired in pairs - cater for that)
                        sidecar["Volumes"] = added_files / int ( sidecar["Slices"] );
                        sidecar["RepetitionTime"] = ( ( seriesdata[2 * sli] ).aqtime - ( seriesdata[0] ).aqtime ) / 2.;

                        // continue if the first and last volume have the same #slices
                        // and if the #files is the product of #slices/vol and #volumes
                        if ( ( unsigned ( sidecar["Slices"] ) == sli ) && ( added_files == ( sli * unsigned ( sidecar["Volumes"] ) ) ) ) {

                            // output size at least 2d, check if 3d or 4d and if yes add those dimensions
                            if ( sli > 1 ) {
                                im_size.push_back ( sli );
                                if ( int ( sidecar["Volumes"] > 1 ) ) {
                                    im_size.push_back ( sidecar["Volumes"] );
                                }
                            }

                            // if one acquisition has more than one volume, change acquisition index to volume index
                            // (if acquisition 0 is 4 slices in 2 volumes, change from { 0,0,0,0,.. } to { 0,0,1,1,.. } )
                            sidecar["NumberOfTemporalPositions"] =
                                im_size[3] / unsigned ( sidecar["NumberOfTemporalPositions"] );
                            for ( unsigned i = 0; i < added_files; i++ )
                                ( seriesdata[i] ).aqnumber = i / int ( sidecar["Slices"] );

                            // resize the data array in bisimage for putting the voxels in
                            supersuper::newsizes ( { im_size[0], im_size[1], im_size[2], im_size[3] } );

                            // sort by tuple value <2> (volume) and if that is the same, by <4> (slice position)
                            // sorting by instance number is not working in the case of interleaved stacks
                            std::sort ( begin ( seriesdata ), end ( seriesdata ), [] ( auto const & t0, auto const & t1 ) {
                                return ( ( ( t0 ).aqnumber < ( t1 ).aqnumber ) ||
                                         ( ( ( t0 ).aqnumber == ( t1 ).aqnumber ) && ( ( t0 ).slipos < ( t1 ).slipos ) ) );
                            } );

                            // get a list of only the files, e.g. extracting tuple<0>
                            // see https://www.fluentcpp.com/2018/11/09/retrieve-firsts-collection-pairs/
                            std::vector<std::string> seriesfiles;
                            std::transform ( begin ( seriesdata ), end ( seriesdata ), std::back_inserter ( seriesfiles ),
                            [] ( auto const & tup ) {
                                return ( tup ).filename;
                            } );

                            // read in the slices -- don't go and look, dirty bit of code...
                            getDicomData ( seriesfiles, superclass::data.size(), superclass::getdata_ptr() );

                            // we need some data from vol 0 slice 0, and vol 0 slice n-1, to get qform
                            // see https://discovery.ucl.ac.uk/id/eprint/1495621
                            DcmFileFormat tmpfile;
                            vec3<float> sli_0_pos, sli_n_pos;

                            tmpfile.loadFile ( seriesfiles[0].c_str() );
                            DcmDataset* tmpdata = tmpfile.getDataset();
                            sli_0_pos[0] = static_cast<float> (
                                               std::stof ( std::string ( getItemString ( tmpdata, DCM_ImagePositionPatient, 0 ).c_str() ) ) );
                            sli_0_pos[1] = static_cast<float> (
                                               std::stof ( std::string ( getItemString ( tmpdata, DCM_ImagePositionPatient, 1 ).c_str() ) ) );
                            sli_0_pos[2] = static_cast<float> (
                                               std::stof ( std::string ( getItemString ( tmpdata, DCM_ImagePositionPatient, 2 ).c_str() ) ) );

                            tmpfile.loadFile ( seriesfiles[sli - 1].c_str() );
                            tmpdata = tmpfile.getDataset();
                            sli_n_pos[0] = static_cast<float> (
                                               std::stof ( std::string ( getItemString ( tmpdata, DCM_ImagePositionPatient, 0 ).c_str() ) ) );
                            sli_n_pos[1] = static_cast<float> (
                                               std::stof ( std::string ( getItemString ( tmpdata, DCM_ImagePositionPatient, 1 ).c_str() ) ) );
                            sli_n_pos[2] = static_cast<float> (
                                               std::stof ( std::string ( getItemString ( tmpdata, DCM_ImagePositionPatient, 2 ).c_str() ) ) );

                            sli--; // we need to use sidecar [ "Slices" ] - 1

                            // Dicom voxel to mm transform
                            mat4<float> Rdcm ( slice_orx[0] * voxsize[0], slice_ory[0] * voxsize[1], ( sli_n_pos[0] - sli_0_pos[0] ) / sli, sli_0_pos[0], 
                                               slice_orx[1] * voxsize[0], slice_ory[1] * voxsize[1], ( sli_n_pos[1] - sli_0_pos[1] ) / sli, sli_0_pos[1],
                                               slice_orx[2] * voxsize[0], slice_ory[2] * voxsize[1], ( sli_n_pos[2] - sli_0_pos[2] ) / sli, sli_0_pos[2], 
                                                                       0,                         0,                                     0,            1  );

                            // NIfTI voxel to mm transform
                            mat4<float> Rnii ( -Rdcm[0][0], -Rdcm[0][1], -Rdcm[0][2], -Rdcm[0][3], 
                                               -Rdcm[1][0], -Rdcm[1][1], -Rdcm[1][2], -Rdcm[1][3], 
                                                Rdcm[2][0],  Rdcm[2][1],  Rdcm[2][2],  Rdcm[2][3], 
                                                         0,           0,           0,           1  );

                            // usually (if slices go top -> bottom) we want to flip the Y axis
							//                  (dicom: positive is RL, nifti: positive is LR)
                            // see https://github.com/rordenlab/dcm2niix/blob/master/console/nii_dicom.cpp#L2659
							if ( sidecar["Slicedirection"] < 0 ) {
								
								Rnii[0][1] *= -1;                                 //
								Rnii[1][1] *= -1;                                 // flip the 2nd (y) column of Rnii
								Rnii[2][1] *= -1;                                 //
								
								vec4<float> yflip_ori ( 0, im_size[1] - 1, 0, 0 ); // yflip_ori: vector of the new origin								
								yflip_ori = ( Rnii * yflip_ori );
								
								Rnii[0][3] -= yflip_ori[0]; //
								Rnii[1][3] -= yflip_ori[1]; // restore origin after flipping
								Rnii[2][3] -= yflip_ori[2]; //

								// and flip the voxels!
								size_t 
									st  = superclass::sizes [ 3 ],
									sz  = superclass::sizes [ 2 ],
									sy  = superclass::sizes [ 1 ],
									sx  = superclass::sizes [ 0 ],
									sy1 = sy-1,
									sy2 = sy/2;
								for ( size_t t = 0; t < st; t++ )
									for ( size_t k = 0; k < sz; k++ )
										for ( size_t j = 0; j < sy2; j++ )
											for ( size_t i = 0; i < sx; i++ )
												std::swap<value_type> ( (*this)[ { i, j, k, t } ], (*this)[ { i, sy1-j, k, t } ] );

							}


                            //
                            // populate the NIfTI header
                            //

                            // set the own (nifti type) header
                            superclass::header = nifti_simple_init_nim();
                            superclass::header->datatype = getdatatype<value_type>();
                            nifti_datatype_sizes (
                                superclass::header->datatype, &( superclass::header->nbyper ), &( superclass::header->swapsize ) );
                            auto dicompath = std::filesystem::path ( seriesfiles[0] );
                            nifti_set_filenames ( superclass::header, dicompath.replace_extension ( "nii.gz" ).c_str(), 0, 0 );

                            // set dimensions and spacings
                            superclass::header->dim[0] = superclass::sizes.size();
                            superclass::header->pixdim[0] = superclass::sizes.size();
                            if ( superclass::header->dim[0] ) {
                                 superclass::header->dim[1] = superclass::sizes[0];
                                 superclass::header->pixdim[1] = voxsize[0];
                                 if ( superclass::header->dim[0] > 1 ) {
                                      superclass::header->dim[2] = superclass::sizes[1];
                                      superclass::header->pixdim[2] = voxsize[1];
                                      if ( superclass::header->dim[0] > 2 ) {
                                           superclass::header->dim[3] = superclass::sizes[2];
                                           superclass::header->pixdim[3] = voxsize[2];
                                           if ( superclass::header->dim[0] > 3 ) {
                                                superclass::header->dim[4] = superclass::sizes[3];
                                                superclass::header->pixdim[4] = float ( sidecar["RepetitionTime"] );
                                                // never encountered 5D but you can go on of course
                                           }    // if >3
                                      }         // if >2
                                 }              // if >1
                            }                   // if >0

                            // get rescale slope and intercept - first try DCM_RescaleSlope/Intercept
                            //        and if that does not work try DCM_RealWorldValueSlope/Intercept
                            superclass::header->scl_slope = getItemDouble ( tmpdata, DCM_RescaleSlope );
                            if ( superclass::header->scl_slope )
                                 superclass::header->scl_slope = getItemDouble ( tmpdata, DCM_RescaleIntercept );
                            else {
                                 superclass::header->scl_slope = getItemDouble ( tmpdata, DCM_RealWorldValueSlope );
                                 superclass::header->scl_inter = getItemDouble ( tmpdata, DCM_RealWorldValueIntercept );
                            } // if scl_slope

                            // set cal_min and cal_max as highest and lowest nonzero values
                            auto ordered_nonzero = std::vector<value_type> ( 0 );
                            for ( auto intensity = supersuper::data.begin(); intensity != supersuper::data.end(); intensity++ )
                                if ( *intensity != 0 )
                                     ordered_nonzero.push_back ( *intensity );
                            std::sort ( ordered_nonzero.begin(),  ordered_nonzero.end() );
                            superclass::header->cal_min = float ( ordered_nonzero.front() );
                            superclass::header->cal_max = float ( ordered_nonzero.back()  );
                            ordered_nonzero.resize(0); // does that clear the memory?

                            // assume mm as vox units, s as time units
                            superclass::header->xyz_units = NIFTI_UNITS_MM;
                            superclass::header->time_units =
                                ( float ( sidecar["RepetitionTime"] ) ) < 3600 ? NIFTI_UNITS_SEC : NIFTI_UNITS_MSEC;

                            // set slice, phase and frequency encoding directions
                            // superclass::header -> dim_info   = FPS_INTO_DIM_INFO ( sidecar [ "FrequencyEncodingDimension" ],
                            // sidecar [ "PhaseEncodingDimension" ], sidecar [ "SlicesDimension" ] );
                            superclass::header->freq_dim  = int ( sidecar["FrequencyEncodingDimension"] ) + 1;
                            superclass::header->phase_dim = int ( sidecar["PhaseEncodingDimension"]     ) + 1;
                            superclass::header->slice_dim = int ( sidecar["SlicesDimension"]            ) + 1;

                            // set Qform in NIfTI header ( mat44 is an internal nifti data type )
                            mat44 m;
                            m.m[0][0] = Rnii[0][0];     m.m[0][1] = Rnii[0][1];     m.m[0][2] = Rnii[0][2];     m.m[0][3] = Rnii[0][3];
                            m.m[1][0] = Rnii[1][0];     m.m[1][1] = Rnii[1][1];     m.m[1][2] = Rnii[1][2];     m.m[1][3] = Rnii[1][3];
                            m.m[2][0] = Rnii[2][0];     m.m[2][1] = Rnii[2][1];     m.m[2][2] = Rnii[2][2];     m.m[2][3] = Rnii[2][3];
                            m.m[3][0] = Rnii[3][0];     m.m[3][1] = Rnii[3][1];     m.m[3][2] = Rnii[3][2];     m.m[3][3] = Rnii[3][3];

                            std::cout << "q-form matrix:" << std::endl;
                            std::cout << m.m[0][0] << " " << m.m[0][1] << " " << m.m[0][2] << " " << m.m[0][3] << " " << std::endl;
                            std::cout << m.m[1][0] << " " << m.m[1][1] << " " << m.m[1][2] << " " << m.m[1][3] << " " << std::endl;
                            std::cout << m.m[2][0] << " " << m.m[2][1] << " " << m.m[2][2] << " " << m.m[2][3] << " " << std::endl;
                            std::cout << m.m[3][0] << " " << m.m[3][1] << " " << m.m[3][2] << " " << m.m[3][3] << " " << std::endl;

                            float qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac = 1;
                            nifti_mat44_to_quatern ( m, &qb, &qc, &qd, &qx, &qy, &qz, &dx, &dy, &dz, &qfac );
                            superclass::header->qform_code = NIFTI_XFORM_SCANNER_ANAT;
                            superclass::header->quatern_b = qb;
                            superclass::header->quatern_c = qc;
                            superclass::header->quatern_d = qd;
                            superclass::header->qoffset_x = qx;
                            superclass::header->qoffset_y = qy;
                            superclass::header->qoffset_z = qz;
                            superclass::header->pixdim[0] = superclass::header->qfac = sidecar["Slicedirection"]; // see nifti1.h

							superclass::header->sform_code = NIFTI_XFORM_SCANNER_ANAT;
							superclass::header->sto_xyz.m[0][0] = Rnii[0][0];
							superclass::header->sto_xyz.m[0][1] = Rnii[0][1];
							superclass::header->sto_xyz.m[0][2] = Rnii[0][2];
							superclass::header->sto_xyz.m[0][3] = Rnii[0][3];
							superclass::header->sto_xyz.m[1][0] = Rnii[1][0];
							superclass::header->sto_xyz.m[1][1] = Rnii[1][1];
							superclass::header->sto_xyz.m[1][2] = Rnii[1][2];
							superclass::header->sto_xyz.m[1][3] = Rnii[1][3];
							superclass::header->sto_xyz.m[2][0] = Rnii[2][0];
							superclass::header->sto_xyz.m[2][1] = Rnii[2][1];
							superclass::header->sto_xyz.m[2][2] = Rnii[2][2];
							superclass::header->sto_xyz.m[2][3] = Rnii[2][3];
							superclass::header->sto_xyz.m[3][0] = Rnii[3][0];
							superclass::header->sto_xyz.m[3][1] = Rnii[3][1];
							superclass::header->sto_xyz.m[3][2] = Rnii[3][2];
							superclass::header->sto_xyz.m[3][3] = Rnii[3][3];

                            // check the header
                            nifti_update_dims_from_array ( superclass::header );

                            //
                            // complete the JSON sidecar and write
                            //
                            sidecar["TotalAcquiredPairs"] = im_size[3] / int ( sidecar["NumberOfTemporalPositions"] );
                            std::ofstream jfile ( dirname + "sidecar_test.json" );
                            jfile << std::setw ( 4 ) << sidecar << std::endl;

                        } // if sli == sz

                    } // series UID not empty

                } // if status.good

            } // constructor that reads from a file

    }; // class dicom

} // namespace bis

#endif // BISDICOMIMAGE_HPP_INCLUDED
