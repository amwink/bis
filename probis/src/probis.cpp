#include "bisdicomimage.hpp"
#include "bismaxtree.hpp"
// #include <cstdlib>

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
        std::copy ( values.begin(), values.end(), image4d.getdata_ptr()->begin() );
        std::cout << *(image4d.getdata_ptr()) << std::endl;
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

        bisnifti<unsigned short> test_image;  
      
		/* 3d test 
		test_image.array_init( { 2, 2, 1, 1, 1,
								 1, 1, 1, 2, 2,      // slice 1
								 
								 3, 2, 1, 1, 1,
								 1, 1, 1, 2, 3 } );  // slice 2
		test_image.reshape ( { 5, 2, 2 } );          // see the shape above */

		/* 2d test */ 
		test_image.array_init ( { 3, 3, 1, 4, 2,     // see 10.1109/ICIP.2007.4379949
								  4, 1, 2, 3, 1 } ); // only 1 slice 
		test_image.reshape    ( { 5, 2 }          ); // see the shape above */
		
        bismaxtree<unsigned short> test_mt ( test_image,                      // image of which to build the maxtree
											 4,                               // number of levels ( test image: 4 )
											 4, // type of connectivity (4|8 for 2D, 6|26 for 3D)
											 "Berger" );                      // method ( "Berger" works, "Wilkinson" will be added)

    }



    return 0;
    
}
