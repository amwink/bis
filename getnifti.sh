rm -rf nifticlib*
git clone https://github.com/NIFTI-Imaging/nifti_clib.git
mkdir -p nifti/include
mkdir -p nifti/src
cd nifti/include
ln -sf ../../nifti_clib/niftilib/nifti1.h .
ln -sf ../../nifti_clib/niftilib/nifti1_io.h
ln -sf ../../nifti_clib/niftilib/nifti1_io_version.h
ln -sf ../../nifti_clib/nifti2/nifti2.h .
ln -sf ../../nifti_clib/nifti2/nifti2_io.h
ln -sf ../../nifti_clib/nifti2/nifti2_io_version.h
cd ../src
ln -sf ../../nifti_clib/niftilib/nifti1_io.c
ln -sf ../../nifti_clib/nifti2/nifti2_io.c
cd ../..
mkdir -p znz/include
mkdir -p znz/src
cd znz/include 
ln -sf ../../nifti_clib/znzlib/znzlib.h
ln -sf ../../nifti_clib/znzlib/znzlib_version.h
cd ../src
ln -sf ../../nifti_clib/znzlib/znzlib.c
cd ../..
