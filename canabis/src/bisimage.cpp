#include "bisimage.hpp"



namespace bis {

    
    
/** Functions outside of the bisimage class.
 *  These can be used by bisimage objects, 
 *  but are not required by the class and do 
 *  not require it themselves.
 *
 */



/** \brief readablefile -- returns true or thows exception
 *
 * the filename is and old-fashioned character array
 * 
 */
bool readablefile ( const char *filename ) {

    std::ifstream file;
    try {
        file.open ( filename );        
        file.get();
        file.close();   
    } catch (std::ifstream::failure & e) {
        std::cerr << "Exception opening/reading/closing file " << filename << "\n";
        return false;
    }
    
    return true;
    
} // readablefile



/** \brief filterline for 1-dimensional convolution
 *
 * signal and filter are both of type std::vector
 */
template <typename T, typename U>
void filterline ( std::vector<T>& signal, const std::vector<U>& filter ) {

    // do circular for now -> after signal[size -1] comes signal[0]
    std::vector <T> workspace ( signal.size() );
    size_t flen2 = filter.size() / 2;

    for ( size_t si = 0; si < signal.size(); si++ )
        for ( size_t fi = 0; fi < filter.size(); fi++ )
            workspace [ (si + flen2) % signal.size() ] += signal [ (si + fi) % signal.size() ] * filter [ fi ];

    std::copy ( workspace.begin(), workspace.end(), signal.begin() );

}


/** \brief fwtline for 1-dimensional fast wavelet transform
 *
 * signal and filterh are both of type std::vector
 * level is of type size_t
 */
template <typename T, typename U>
void fwtline ( std::vector<T>& signal, std::vector<U> filterh, size_t level ) {

    // do circular for now -> after signal[size -1] comes signal[0]
    size_t
    newlength = signal.size();
    std::vector <U>
    filterg = ( filterh );

    // filter g not only the reverse of h, every second coefficient must be - as well
    std::reverse ( filterg.begin(), filterg.end() );
    for ( size_t i=1; i < filterg.size(); i+=2 )
        filterg[i] *= -1;

    for ( size_t l=0; l<level; l++ ) {

        std::vector <T>
        tmp ( newlength, 0 );

        newlength /= 2;

        for ( size_t wi=0; wi<newlength; wi++ )
            for ( size_t fi=0; fi<filterh.size(); fi++ ) {

                tmp[wi            ] += signal [ (2*wi + fi) % signal.size() ] * filterh [fi];
                tmp[wi + newlength] += signal [ (2*wi + fi) % signal.size() ] * filterg [fi];

            } // for wi

        std::copy ( tmp.begin(), tmp.end(), signal.begin() );

    } // for l

} // fwtline



/** \brief ifwtline for 1-dimensional fast inverse wavelet transform (reconstruction)
 *
 * signal and filterh are both of type std::vector
 * level is of type size_t
 */
template <typename T, typename U>
void ifwtline ( std::vector<T>& signal, std::vector<U> filterh, size_t level ) {

    // do circular for now -> after signal[size -1] comes signal[0]
    size_t
    newlength = signal.size();
    std::vector <U>
    filterg = ( filterh );

    // filter g not only the reverse of h, every second coefficient must be 0 as well
    std::reverse ( filterg.begin(), filterg.end() );
    for ( size_t i=1; i < filterg.size(); i+=2 )
        filterg[i] *= -1;

    // length of approximation is 2^level times smaller
    for ( size_t i=0; i<level; i++ )
        newlength /= 2;

    for ( size_t l=0; l<level; l++ ) {

        std::vector <T>
        tmpc ( newlength*2, 0 );
        std::vector <T>
        tmpd ( newlength*2, 0 );

        // upsample c and d from current level
        for ( size_t i=0; i<newlength; i++ ) {
            tmpc [ 2 * i ] = signal [ i             ];
            tmpd [ 2 * i ] = signal [ newlength + i ];
        }

        std::fill( signal.begin(), signal.end(), 0 );
        newlength *=2;

        for ( size_t wi=0; wi<newlength; wi++ )
            for ( size_t fi=0; fi<filterh.size(); fi++ ) {

                signal[wi] += tmpc [ (wi - fi) % tmpc.size() ] * filterh[fi];
                signal[wi] += tmpd [ (wi - fi) % tmpd.size() ] * filterg[fi];

            } // for wi

    } // for l

} // ifwtline



} // namespace