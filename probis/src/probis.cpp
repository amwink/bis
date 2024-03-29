#include "bisdicomimage.hpp"
#include "bismaxtree.hpp"
#include <cstdlib>

// using the zip library
#include <zip.h>

using namespace bis;



// forward declaration for a helper function
std::vector < std::filesystem::path > fullpath_extract ( std::filesystem::path dicomsdir, 
														 std::filesystem::path dicomzip );



// running the tests
int main ( int argc, char* argv[] )
{
    
    
    
    // if fsldir found
    std::string fsldir;
    if ( std::getenv( "FSLDIR" ) )
        fsldir = ( std::getenv( "FSLDIR" ) );
    else
        fsldir = "/usr/local/fsl";
    std::string 
        dcmdir;
        
    std::cout << "hello from canabis! \n";

    auto 
        do_2d      = true,	// test 2D vectors & matrices
        do_3d      = true,	// test 3D vectors & matrices ( only tested if do_2d )
        do_4d      = true,	// test 3D vectors & matrices ( only tested if do_3d )    
        do_nifti   = true,	// test nifti file interactions
        do_dicom   = false,	// test dicom import
        do_reduce  = true,	// test dimensionality reducer 
		do_maxtree = true;	// test maxtree

    ////////////////////////////////////////
    // 1.
    //
    // Test vectors and matrices
    //


if ( do_2d || do_3d || do_4d ) {

	// vectors        
	vec2<float>
		v2_1 ( { 1, 1 } ),
		v2_2 ( v2_1 );

	
		if ( do_2d ) {	

			std::cout << "2D\nvector v2_1: \t";
			std::cout << v2_1 << "\n";

			std::cout << "vector v2_2: \t" << v2_2 << " \n";
			std::cout << "v2_2 + 2 v2_1: \t" << v2_2 + v2_1 * float(2.) << " \n";
			
			// matrices		
			mat2<float>
				m2_1 ( { 1, 2, 2, 1 } ),
				m2_2 ( - m2_1 );
			m2_2[0][0] = m2_1[0][0];
			m2_2[1][1] += 4;
				
			std::cout << "matrix m2_1: \n" << m2_1 << " \n";
			std::cout << "inverse m2_1: \n" << m2_1.inverse() << " \n";
			
			auto e2_1 = m2_1.diagonalise_sym();
			std::cout << "eigenvalues and -vectors\n";
			std::cout << "values\n"  << e2_1.v << "\n";
			std::cout << "vectors\n" << e2_1.m << "\n";
			
			std::cout << "matrix m2_2: \n" << m2_2 << " \n";
			std::cout << "inverse m2_2: \n" << m2_2.inverse() << " \n";			
			auto e2_2 = m2_2.diagonalise_sym();
			std::cout << "eigenvalues and -vectors\n";
			std::cout << "values\n"  << e2_2.v << "\n";
			std::cout << "vectors\n" << e2_2.m << "\n";
									
		}
		
        if( do_3d || do_4d ) {

            // vectors
            vec3<float> v3_1(v2_1);
            v3_1[2] = 30;

			if (do_3d) {

				std::cout << "3D\nvector v3_1: \t";
				std::cout << v3_1 << " \n";

				auto v3_2(v3_1);
				v3_2 += 1.;
				std::cout << "vector v3_2: \t" << v3_2 << " \n";
				std::cout << "v3_2 + 2 v3_1: \t" << v3_2 + v3_1 * float(2.) << " \n";
		
				// matrices
				mat3<double> 
					m3_1 ( { 1, 1, 0, 1, 1, 0, 0, 0, 1 } );
		 
				std::cout << "matrix m3_1: \n" << m3_1 << " \n";
				std::cout << "inverse m3_1: \n" << m3_1.inverse() << " \n";
				
				auto e3_1 = m3_1.diagonalise_sym();
				std::cout << "eigenvalues and -vectors\n";
				std::cout << "vectors\n" << e3_1.v << "\n";
				std::cout << "vectors\n" << e3_1.m << "\n";
			
				auto m3_2 ( m3_1 + m3_1 );
				m3_2 [0][2] = -1;

				std::cout << "matrix m3_2: \n" << m3_2 << " \n";
				std::cout << "inverse m3_2: \n" << m3_2.inverse() << " \n";

			}
        
            if ( do_4d ) {

                // vectors
                vec4<float> v4_1(v3_1);
                v4_1[3] = 40;

                std::cout << "4D\nvector v4_1: \t";
                std::cout << v4_1 << " \n";

                auto v4_2(v4_1);
                v4_2 += 1.;
                std::cout << "vector v4_2: \t" << v4_2 << " \n";
                std::cout << "v4_2 + 2 v4_1: \t" << v4_2 + v4_1 * float(2.) << " \n";

                // matrices
                mat4<double> 
                    m4_1 ( { 1,  0,  0, -1, 
                             0,  1,  0,  0, 
                             0,  0,  1,  0, 
                             0,  0,  0,  1  } );
     
                std::cout << "matrix m4_1: \n" << m4_1 << " \n";
                std::cout << "inverse m4_1: \n" << m4_1.inverse() << " \n";
            
				// auto m4_2 ( m4_1 );
                mat4<double> 
                    m4_2 ( { 1,  0,  2,  1, 
                             0,  2,  0,  2, 
                             2,  0,  3,  0, 
                             1,  2,  0,  4  } );

                std::cout << "matrix m4_2: \n" << m4_2 << " \n";
                std::cout << "inverse m4_2: \n" << m4_2.inverse() << " \n";
				
				auto e4_2 = m4_2.diagonalise_sym();
				std::cout << "eigenvalues and -vectors\n";
				std::cout << "values\n"  << e4_2.v << "\n";
				std::cout << "vectors\n" << e4_2.m << "\n";

            } // if do_4d

        } // if do_3d

    } // if do_2d



    if ( do_nifti ) {
        
        ////////////////////////////////////////
        // 1.
        //
        // Test NIfTI I/O
        //

        std::string test_file ( fsldir + "/data/standard/MNI152_T1_2mm_brain.nii.gz" );
        std::cout << "loading " << test_file << " ..." << std::endl;
        bisnifti<float> test_image ( test_file, bis::DO_READ_DATA );        
        test_image.write ( "/tmp/MNI152_T1_2mm_brain.nii.gz" );   
    
    }
      


    if ( do_dicom ) {

        ////////////////////////////////////////
        // 2.
        //
        // Test DICOM I/O
        //
		// Find zip-files of dicoms in <bis directory>/data/dicom
		// and extract and convert each of them. Pixel array sizes,
		// coordinate spaces, bvecs and bvals can then be checked.
		
		// first find the bis directory from which we were called
		auto 
			mypath = std::filesystem::path(argv[0]).parent_path();

		// given that we're in */bis/cmake-build-Debug/output" (or Release)
		// we can now find bis/data/dicoms
		auto 
			dicomsdir = mypath.parent_path().parent_path();
		dicomsdir += "/data/dicom";
		
		std::vector < std::filesystem::path > 
			dicomzips;
		for ( auto& possible_zip: std::filesystem::recursive_directory_iterator ( dicomsdir ) )
			if ( ( possible_zip.path().extension() == ".zip" ) && 
				 ( possible_zip.path().filename() != "philipsDicom.zip" ) && 
				 ( possible_zip.path().filename() != "siemensDicom.zip" ) )
				dicomzips.push_back ( possible_zip.path() );		

		for ( auto f: dicomzips ) std::cout << f << "\n";

		for ( auto dicomzip: dicomzips ) {
		//for ( auto dicomzip: dicomzips | std::views::reverse ) {
			
			// put in the faulty ones first 
			auto 
				dicomtest = "/home/amwink/software/cpp/bis/data/dicom/siemensDicom/dwi_dwi/MRI_DTI_MR-EPAD-SingleShell-DTI48-19_.zip";
			dicomzip  = dicomtest;
			std::cout << "unpacking " << dicomzip << "... " << std::endl;
			
			auto 
				dicomlist = fullpath_extract ( dicomsdir, dicomzip );
			auto 
				counter = 0;
				
			// list files (optional)
			bool listfiles = true;
			if ( listfiles )
				for ( auto v: dicomlist )
					printf ( "dicom %04d in %s: %s\r", counter++, dicomsdir.c_str(), v.c_str() );

			// find the first dicom (assuming there's one series in the ZIP file) 
			DcmFileFormat 
				dcmfile;
			for ( counter = 0; auto v: dicomlist ) {
				if ( dcmfile.loadFile ( v.c_str() ).good() )
					break;
				counter++;
			}
			
			std::string 
				dcmslice = dicomlist[counter].string();
			std::cout << "converting: " << dcmslice << "... " << std::endl;	
			
			// make a dicom object of this DICOM's series and write as a bids image
			bis::bisdicom<unsigned short> test_dicom ( dcmslice );

			std::string 
				newname = dcmslice.substr ( 0, dcmslice.find_last_of('.') ) + ".nii.gz";
			std::cout << "writing: " << newname << "... " << std::endl;	
			test_dicom.write ( newname );   
			
			if ( dicomzip == dicomtest ) {
				std::cout << "Using the test file " << dicomtest << ", 1 DICOM conversion only.";
				break;
			}
			
		} // for dicomzip
	
	} // if do_dicom
    
    
    
    if ( do_reduce ) {

        ////////////////////////////////////////
        // 2.
        //
        // Test reducers -- functions applied on one dimension
        //

        bisimage <int> image4d ( { 2, 3, 4, 5 } );
        std::vector <int> values ( 120 );
        std::iota ( values.begin(), values.end(), 0 );
        image4d.vector_import( values );
        std::cout << image4d.vector_export() << std::endl;
		std::cout << "element { 0, 1, 2, 3 }: " << image4d ( { 0, 1, 2, 3 } ) << std::endl;
		//std::cout << image4d << std::endl;
        for ( auto dim: {3,2,1,0} ) {
            auto sum_image = image4d.sum(dim);
            std::cout << sum_image << std::endl;
        }
		std::string mri_name ( fsldir + "/data/standard/MNI152_T1_2mm.nii.gz" );
        std::cout << "loading " << mri_name << " ..." << std::endl;
        bisnifti<unsigned short> 
			mri_data ( mri_name, bis::DO_READ_DATA );
		mri_data = mri_data.mean( 2 );
		mri_data.write ( "/tmp/2dmean.nii.gz" );
    }



    if ( do_maxtree ) {

        ////////////////////////////////////////
        // 2.
        //
        // Test maxtree -- for doing mathematical morphology
        //

		std::cout << "FSL directory: " << bis::getenv( "FSLDIR" ) << std::endl;

		std::string mri_name ( fsldir + "/data/standard/MNI152_T1_2mm.nii.gz" );
        std::cout << "loading " << mri_name << " ..." << std::endl;
        bisnifti<unsigned short> 
			mri_data ( mri_name, bis::DO_READ_DATA );						// load as nifti image
        std::cout << "building maxtree ..." << std::endl;
        bismaxtree<unsigned short> mri_mt (	mri_data,						// image of which to build the maxtree
											8,								// number of levels ( test image: 4 )
											2 * mri_data.getsize().size(),	// type of connectivity (4|8 for 2D, 6|26 for 3D)
											"Berger" );						// method ( "Berger" works, hopefully will add others)

		bisimage < unsigned short >
			mri_labels ( mri_mt );											// first cast maxtree to bisnifti: labels -> voxels
		bisnifti < unsigned short > 
			nifti_labels ( mri_labels );									// then cast bisimage to bisnifti (minimal header info!)
		nifti_labels.write ( "/tmp/mri_mt.nii.gz" );

		mri_mt.addattr ( "mass" );

		auto cnum = 0;
        std::cout << "writing mask of component " << cnum << " and kids ..." << std::endl;	
		auto 
			mri_set = mri_mt.setpoints ( 1, 0, DO_SORT, DO_LEVEL );			// auto: output type setpoints() fixed
		bisnifti < unsigned short >											// not auto: cast bisimage to bisnifti
			nifti_set ( mri_set );		
		nifti_set.write ( "/tmp/mri_mask.nii.gz" );

    }



    return 0;
    
}



std::vector < std::filesystem::path > fullpath_extract ( std::filesystem::path dicomsdir, 
														 std::filesystem::path dicomzip ) {
															 
	std::vector < std::filesystem::path > 
		file_list;

	int  
		zip_err = 0;
	auto 
		zip_archive = zip_open ( dicomzip.c_str(), 0, &zip_err );
	auto  
		files_total = zip_get_num_files ( zip_archive );
	struct zip_stat
		zst;
	
	// init the zip archive stats, including its list of files
	zip_stat_init ( &zst );
	
	// extract all files from the list
	for ( int i = 0; i < files_total; i++ ) {
		
		auto 
			file_in_zip = zip_fopen_index ( zip_archive, i, 0 );
		
		if ( file_in_zip ) {

			// init the archived file's stats,, including its name and uncompressed size
			zip_stat_index ( zip_archive, i, 0, &zst );
			std::string
				fname ( zst.name );								
			auto 
				*contents = new char [ zst.size ];
			auto 
				inpath = dicomsdir;			

			// get the current name -- assume files in directories occur after them in the list
			inpath += ( "/tmp/" + fname );
			
			// if directory (trailing '/') create that, if file extract it
			if ( *fname.rbegin() == '/' )
				std::filesystem::create_directories ( inpath );				
			else {
				zip_fread ( file_in_zip, contents, zst.size );
				if ( std::ofstream ( inpath , std::ofstream::binary ).write( contents, zst.size ) )
					file_list.push_back( inpath );
			} // if file not dir
				
			// free the space used for extraction
			delete [] contents;
				
		} // if file_in_zip
		
	} // for i	

	return ( file_list );

}



