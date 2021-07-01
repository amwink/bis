#include "bisdicomimage.hpp"
#include "bismaxtree.hpp"
#include <cstdlib>

using namespace bis;

int main()
{
    
    
    
    // if fsldir found
    std::string fsldir;
    if ( std::getenv( "FSLDIR" ) )
        fsldir = ( std::getenv( "FSLDIR" ) );
    else
        fsldir = "/usr/local/fsl";
    std::string 
        dcmdir = "/home/amwink/work/projects/epad/example_data/data/epad-download/040-00003/Screening/ASL/040-00003_Screening_MRI_ASL_MR-EPAD_ASL_2025_version1";
        
    std::cout << "hello from canabis! \n";
    std::cout << "fsl directory: " << fsldir << " \n";

    auto 
        do_2d    = false,  // test 2D vectors & matrices
        do_3d    = false,  // test 3D vectors & matrices ( only tested if do_2d )
        do_4d    = false,  // test 3D vectors & matrices ( only tested if do_3d )    
        do_nifti = false,   // test nifti file interactions
        do_dicom = false,  // test nifti file interactions
        do_reduce = false, // test dimensionality reducer 
		do_maxtree = true; // test maxtree

    ////////////////////////////////////////
    // 1.
    //
    // Test vectors and matrices
    //



    float  
        v[2] = { 10, 20 };
    double 
        m[4] = { 1, -1, 0, 1 };
    vec2<float>  
        v2_1 ( v );
    mat2<double> 
        m2_1 ( m );

    if(do_2d) {

        // vectors
        
        std::cout << "2D\nvector v2_1: \t";
        std::cout << v2_1 << "\n";

        auto v2_2 ( v2_1 );
        v2_2 += 1.;

        std::cout << "vector v2_2: \t" << v2_2 << " \n";
        std::cout << "v2_2 + 2 v2_1: \t" << v2_2 + v2_1 * float(2.) << " \n";
        
        // matrices
        
        std::cout << "matrix m2_1: \n" << m2_1 << " \n";
        std::cout << "inverse m2_1: \n" << m2_1.inverse() << " \n";
        
        auto m2_2 ( m2_1 + m2_1 );
        m2_2 [0][1] = -1;

        std::cout << "matrix m2_2: \n" << m2_2 << " \n";
        std::cout << "inverse m2_2: \n" << m2_2.inverse() << " \n";
            
        if(do_3d) {

            // vectors
            vec3<float> v3_1(v2_1);
            v3_1[2] = 30;

            std::cout << "3D\nvector v3_1: \t";
            std::cout << v3_1 << " \n";

            auto v3_2(v3_1);
            v3_2 += 1.;
            std::cout << "vector v3_2: \t" << v3_2 << " \n";
            std::cout << "v3_2 + 2 v3_1: \t" << v3_2 + v3_1 * float(2.) << " \n";
    
            // matrices
            mat3<double> 
                m3_1 ( { 1, 0, -1, 0, 1, 0, 0, 0, 1 } );
     
            std::cout << "matrix m3_1: \n" << m3_1 << " \n";
            std::cout << "inverse m3_1: \n" << m3_1.inverse() << " \n";
            
            auto m3_2 ( m3_1 + m3_1 );
            m3_2 [0][2] = -1;

            std::cout << "matrix m3_2: \n" << m3_2 << " \n";
            std::cout << "inverse m3_2: \n" << m3_2.inverse() << " \n";
        
            if(do_4d) {

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
            
                auto m4_2 ( 2.  * m4_1 );
                m4_2 [0][3] = -1;

                std::cout << "matrix m4_2: \n" << m4_2 << " \n";
                std::cout << "inverse m4_2: \n" << m4_2.inverse() << " \n";

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
        test_image.saveNII( dcmdir + "/MNI152_T1_2mm_brain.nii.gz" );   
    
    }
      


    if ( do_dicom ) {

        ////////////////////////////////////////
        // 2.
        //
        // Test DICOM I/O
        //

        std::string test_file2 ( dcmdir + "/EPAD_040-00003_Screening_002159.dcm" );
        bis::bisdicom<unsigned short> test_image2 ( test_file2 );
        std::cerr << "file successfully read\n";
        test_image2.saveNII( dcmdir + "/EPAD_040-00003_Screening_002159.nii.gz" );   

    }
    
    
    
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
        for ( auto dim: {0,1,2,3} ) {
            auto sum_image = image4d.sum(dim);
            std::cout << sum_image << std::endl;
        }
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
		nifti_labels.saveNII ( "/tmp/mri_mt.nii.gz" );

		mri_mt.getattr ( "mass" );

		auto cnum = 240;
        std::cout << "writing mask of component " << cnum << " and kids ..." << std::endl;
		auto 
			mri_set = mri_mt.setpoints ( cnum, 0 );							// auto: output type setpoints() fixed
		bisnifti < unsigned short >											// not auto: cast bisimage to bisnifti
			nifti_set ( mri_set );		
		nifti_set.saveNII ( "/tmp/mri_mask.nii.gz" );

    }



    return 0;
    
}
