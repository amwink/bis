
/*********************************************************************************/
/*                                                                               */
/*     Command-line tool to access support vector machine                        */
/*                          classification functionality.                        */
/*                                                                               */
/*     Alle Meije Wink                                                           */
/*                                                                               */
/*********************************************************************************/

/*
  Update history

  Who    When       What
  AMW    13-10-12   creation

*/

#define DLIB_PNG_SUPPORT

#include <dlib/image_io.h>                              //        -- be able to write png
#include <dlib/svm.h>                                   // n-dimensional vectors (n!=3) should use matrix in dlib
#include <dlib/matrix.h>                                //        -- see http://dlib.net/linear_algebra.html#vector
#include <dlib/svm/svm_c_linear_dcd_trainer.h>          // use svm trainer that supports "warm starting"
                                                        // see http://dlib.net/dlib/svm/active_learning.h.html

#include "combisDesign.hpp"

using namespace dlib;

int combisSVM(bis::bisnifti<value_type> *currentImage, std::string designfile)
{
    // test combis with
    // -i ~/work/documents/memorabel/PRNI2016/vumc/ECM/allmask_fmri.nii.gz --svm ~/work/documents/memorabel/PRNI2016/vumc/ECM/image_matrix.json -o ~/work/documents/memorabel/PRNI2016/vumc/ECM/weights.nii.gz

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    // main data types

    // classification table
    std::vector<value_type>                   y;
    std::vector<dlib::matrix<value_type,0,1>> x;

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    // data input from file

    // read the design
    auto design_isbinary=combisSVMdesign(designfile, &x, &y, currentImage);   // read text / binary design

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    // kernel design


    // type of sample (1 image) and kernel
    typedef dlib::matrix<value_type,0,1> sample_type;                              // one sample == one row of x
    typedef dlib::linear_kernel<sample_type> kernel_type;                         // kernel data type

    std::vector<sample_type> samples;                                             // sample set container
    std::vector<value_type> labels;                                                // label set container

    samples.assign(x.begin(), x.end());
    labels.assign(y.begin(), y.end());

    for (auto l: labels)
        std::cout << l << std::endl;
    if (!design_isbinary)
        for (auto s: samples)
            std::cout << s << std::endl;

    // trainer for this type of kernel
    // This trainer solves the "C" formulation of the SVM.  See the documentation for
    // details.
    dlib::svm_c_linear_dcd_trainer<kernel_type> linear_dcd_trainer;
    linear_dcd_trainer.set_c(1000);

    // normalise samples of x -- see http://dlib.net/svm_ex.cpp.html
    //dlib::vector_normalizer<sample_type> normalizer;
    //normalizer.train(samples);
    //for (auto sx:samples)
    //    sx = normalizer(sx);

    // preserve the state of the classifier for warm-starting (see active_learing.h)
    typedef typename dlib::svm_c_linear_dcd_trainer<kernel_type>::optimizer_state optimizer_state;
    optimizer_state state;

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    // construct decision / projection vector
    typedef decision_function< kernel_type > dftype;
    dftype decision_function = linear_dcd_trainer.train(samples, labels, state);

    sample_type m(2);
    m(0)=m(1)=0;
    for (long i = 0; i < decision_function.alpha.nr(); ++i) {
        std::cout << i << std::endl;
        std::cout << decision_function.alpha(i) << std::endl;
        std::cout << decision_function.basis_vectors(i) << std::endl;
        m += decision_function.alpha(i) * decision_function.basis_vectors(i);
    }
    std::cout << "b:" << std::endl;
    std::cout << decision_function.b << std::endl;

    #define RED     "\033[31m"      /* Red */
    #define GREEN   "\033[32m"      /* Green */
    #define YELLOW  "\033[33m"      /* Yellow */
    #define RESET   "\033[0m"

    dlib::array2d<float> map_image;
    float top=5,step=.49;
    map_image.set_size(top/step+1,top/step+1);
    for (float j=top; j>=0.; j-=step) {
        for (float i=0.1; i<=top; i+=step) {
            m(0)=i; m(1)=j;
            auto f = decision_function(m);
            map_image[i][j]=f;
            if (f<-1.)
                printf("%1.01f,%1.01f -> %s%5.02f%s  ",i,j,GREEN,-1.,RESET);
            else if (f>1.)
                printf("%1.01f,%1.01f -> %s%5.02f%s  ",i,j,RED,1.,RESET);
            else
                printf("%1.01f,%1.01f -> %s%5.02f%s  ",i,j,YELLOW,0.,RESET);
            }
        std::cout << std::endl; }
    //dlib::save_png(map_image, "/tmp/map_image.png");

    #undef RED
    #undef GREEN
    #undef YELLOW
    #undef RESET

    /*

    if (!my_machine.margin_set.empty())
        for (size_t k=0; k<my_machine.margin_set.size(); k++) {
            if (my_machine.output(my_machine.margin_key[k]))
                projection += (my_machine.weight[k] * x[my_machine.margin_key[k]]);
            else
                projection -= (my_machine.weight[k] * x[my_machine.margin_key[k]]);
        }
    if (!my_machine.error_set.empty())
        for (size_t i=0; i<my_machine.error_set.size(); i++) {
                projection -= (C * x[my_machine.every_key[my_machine.error_set[i]]]);
        }
    if (!my_machine.error_star_set.empty())
        for (size_t i=0; i<my_machine.error_star_set.size(); i++) {
                projection -= (C * x[my_machine.every_key[my_machine.error_star_set[i]]]);
        }
    auto projectionbias=my_machine.bias;

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    // demonstrate output

    if (y.size() == 4)
    {

        // test data with 0..5 x 0..5 grid and print classes
        double stp=1.0;//                                                             // sample the input space and show the classification
        for ( double y=5.; y>=0.; y-=stp )                                            // at all sampled points
        {
            for ( double x=0.; x<=5.; x+=stp )
            {
                ublas::vector <value_type> testxy(2);
                testxy[0]=x;
                testxy[1]=y;
                std::cout << "("
                        << std::fixed << std::setprecision(1) << x << ","
                        << std::fixed << std::setprecision(1) << y << ") -> "
                        << my_machine (testxy) << ", " ;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;

        // sample the input space and show the classification
        for ( double y=5.; y>=0.; y-=stp )                                            // at all sampled points
        {
            for ( double x=0.; x<=5.; x+=stp )
            {
                ublas::vector <value_type> testxy(2);
                testxy[0]=x;
                testxy[1]=y;
                std::cout << "("
                        << std::fixed << std::setprecision(1) << x << ","
                        << std::fixed << std::setprecision(1) << y << ") -> "
                        << projectionbias + blas::dot ( projection, testxy ) << ", " ;
            }
            std::cout << std::endl;
        }

        std::cout << "weights:\n";
        auto alpha=my_machine.weight;
        for (auto a: alpha) std::cout << a << " ";
        std::cout << std::endl;
        std::cout << "margin vectors:\n";
        {auto vecset=my_machine.margin_set;
        for (auto v: vecset) std::cout << v << " ";}
        std::cout << std::endl;
        std::cout << "bias:\n" << my_machine.bias << std::endl;
        std::cout << "C:\n" << my_machine.C << std::endl;

        std::cout << "error vectors:\n";
        {auto vecset=my_machine.error_set;
        for (auto v: vecset) std::cout << v << " ";}
        std::cout << std::endl;
        std::cout << "remaining vectors:\n";
        {auto vecset=my_machine.remaining_set;
        for (auto v: vecset) std::cout << v << " ";}
        std::cout << std::endl;

    }

    // if images -> make a projection image (weights map)
    // that shows the voting rights of each brain region
    if (design_isbinary) { // images design

        currentImage->bisArray::operator*=(0);
        {size_t i=0;
            for (auto m: mask)
                currentImage->my_data[m]=projection[i++]/mstd;
        }
        std::cout << "mask size: " << mask.size() << std::endl;

        std::cout << "projections of training images on weight map: " << std::endl;
        for (auto k:keys)
            std::cout << "specimen " << k << ", class " << y[k] << " SVM output " << my_machine(x[k]) << " projection " << projectionbias + blas::dot ( projection, x[k] ) << std::endl;

        auto iminfo=nifti_copy_nim_info(currentImage->getHeader());

        // std::cout << "storing bias in toffset " << std::endl;
        iminfo->toffset=float(projectionbias);
        //iminfo->scl_inter=float(projectionbias);

        // std::cout << "map in current image " << std::endl;
        currentImage->setHeader(iminfo);

        //nifti_image_infodump(iminfo);

    }

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    // cross-validation

    auto cv=1;

    if (cv) {

    std::vector<double> Cs = {1./16000., 1./8000, 1./4000., 1./2000., 1./1000.,};

    for (auto C2: Cs) {

        std::cout << "C value = " << std::setw(6) << C2 << std::endl;

        unsigned n=0,p=0,tn=0,tp=0,sv=0;
        for (auto k2:keys) {
            auto keys2 = keys;
            keys2.erase(keys2.begin()+k2);
            onlinesvm_machine_type my_machine2( C2, my_kernel, training_data );        // C = 1.0
            my_machine2.learn( keys2.begin(), keys2.end() );                           // start the learning
            auto in2=y[k2];
            auto out2=2*my_machine2(x[k2])-1;
            std::cout << "crossval " << k2 << ", class " << in2 << " SVM output " << out2
                      << ", #SV " << my_machine2.margin_set.size()+my_machine2.margin_set.size()
                      << ", bias " << my_machine2.bias << std::endl;
            if (in2<0)
                {n++; tn+=(out2<0);}
            else
                {p++; tp+=(out2>0);}
            sv=sv+my_machine2.margin_set.size();
        }

        sv/=(keys.size()-1);
        auto tnr=double(tn)/double(n);
        auto tpr=double(tp)/double(p);
        std::cout << "  true negative ratio = " << tnr
                  << ", true positive ratio = " << tpr
                  << ", balanced accuracy   = " << (tnr+tpr)/2.
                  << ", average #SV = " << sv << std::endl;


    }

    }
    */
	
	return 0;
	
}


