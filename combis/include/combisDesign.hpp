#ifndef COMBISDESIGN_HPP_INCLUDED
#define COMBISDESIGN_HPP_INCLUDED

#include <fstream>
#include "bisniftiimage.hpp"

#include <dlib/matrix.h>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

typedef float_t value_type;

bool combisSVMdesign(std::string                               designfile,
                     std::vector <dlib::matrix<value_type,0,1>> *x,
                     std::vector <value_type>                   *y,
                     bis::bisnifti<value_type>              *curr) {

    bool isbinary = true;
    json jdesign;

    try {

        // read JSON design file
        std::ifstream infile(designfile);
        infile >> jdesign;

        // load labels
        auto    ytemp ( jdesign[ "labels" ].get<std::vector<value_type>>() );
        y->assign(ytemp.begin(),ytemp.end());

        // determine whether we we have a text design (ascii matrix) or binary inputs (images)
        std::string designtype=jdesign["designtype"].get<std::string>();

        // text vectors: read from ascii matrix
        if (!designtype.compare("text")) {

            isbinary = false;

            std::vector<std::vector<value_type>> xtesti(jdesign["vectors"].get<std::vector<std::vector<value_type>>>());

            // load a vector as a sequence of floats (converted already by json)
            for (auto v: xtesti) {

                dlib::matrix<value_type,0,1> xi = dlib::mat(v);
                x->push_back(xi);

            }

        } else { // binary vectors: read from images

            // test for a mask in json
            std::string maskFile = jdesign["mask"].get<std::string>();

            bis::bisnifti<value_type> *maskImage = new bis::bisnifti<value_type>(maskFile);
            auto 
				nzvec = maskImage->nonzeros();

            std::vector<std::string> xtests(jdesign["vectors"].get<std::vector<std::string>>());
            for (auto v: xtests) {
				
                bis::bisnifti<value_type> 
					oneImage = bis::bisnifti<value_type>(v);
                auto 
					selectedvoxels = oneImage.vector_get ( nzvec );

                dlib::matrix<value_type,0,1> xs = dlib::mat(selectedvoxels);
                x->push_back(xs);
            }

       }

    } catch(const json::parse_error &) {
        throw (bis::bisException("error parsing design json"));
    }

    return ( isbinary );

} // readdesign

#endif // COMBISDESIGN_HPP_INCLUDED
