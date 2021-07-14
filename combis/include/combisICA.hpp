
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

#include "combisDesign.hpp"
#include <dlib/svm.h>

using namespace dlib;

namespace cb = canabis;                       // computer analysis of brain image statistics

int combisICA(cb::bisNiftiImage<intensity> *currentImage, std::string designfile)
{
    // classification table
    std::vector<size_t>               mask(1);
    std::vector<intensity>               y(1);
    std::vector<std::vector <intensity>> x(1);

    // read the design
    auto design_isbinary=combisDesign(designfile, &x, &y, &mask, currentImage);                    // read text / binary design

    double noiseStdDev = 0.25; //0.175;
    size_t nReplicates = 10;   //30;
    size_t nAngles     = 10;   //150;
    size_t nSweeps     = 0;

    auto xs  = x.size();
    auto x0s = x[0].size();
    arma::mat matX ( x0s, xs ); // column major

    for ( size_t i=0; i<xs; i++ )
        for (size_t j=0; j<x0s; j++)
            matX(j,i)=x[i][j];

    // initialise radical
    nSweeps = matX.n_rows - 1;
    //mlpack::radical::Radical rad(noiseStdDev, nReplicates, nAngles, nSweeps);

    // run radical
    arma::mat matY;
    arma::mat matW;
    //rad.DoRadical(matX.t(), matY, matW);

    std::cout << "completed";
}


