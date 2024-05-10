/*=========================================================================
*
*  Copyright NumFOCUS
*
*  Licensed under the Apache License, Version 2.0 (the "License");
*  you may not use this file except in compliance with the License.
*  You may obtain a copy of the License at
*
*         http://www.apache.org/licenses/LICENSE-2.0.txt
*
*  Unless required by applicable law or agreed to in writing, software
*  distributed under the License is distributed on an "AS IS" BASIS,
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*  See the License for the specific language governing permissions and
*  limitations under the License.
*
*=========================================================================*/
#ifndef sitkRealToHalfHermitianForwardFFTImageFilter_h
#define sitkRealToHalfHermitianForwardFFTImageFilter_h

/*
 * WARNING: DO NOT EDIT THIS FILE!
 * THIS FILE IS AUTOMATICALLY GENERATED BY THE SIMPLEITK BUILD PROCESS.
 * Please look at sitkImageFilterTemplate.h.in to make changes.
 */

#include <memory>

#include "sitkBasicFilters.h"
#include "sitkImageFilter.h"

namespace itk {
  namespace simple {

    /**\class RealToHalfHermitianForwardFFTImageFilter
\brief Base class for specialized real-to-complex forward Fast Fourier Transform .

This is a base class for the "forward" or "direct" discrete Fourier Transform . This is an abstract base class: the actual implementation is provided by the best child class available on the system when the object is created via the object factory system.

This class transforms a real input image into its complex Fourier transform. The Fourier transform of a real input image has Hermitian symmetry: \f$ f(\mathbf{x}) = f^*(-\mathbf{x}) \f$ . That is, when the result of the transform is split in half along the X-dimension, the values in the second half of the transform are the complex conjugates of values in the first half reflected about the center of the image in each dimension. This filter takes advantage of the Hermitian symmetry property and reduces the size of the output in the first dimension to N/2+1, where N is the size of the input image in that dimension and the division by 2 is rounded down.

\see HalfHermitianToRealInverseFFTImageFilter 


\see ForwardFFTImageFilter
\sa itk::simple::RealToHalfHermitianForwardFFT for the procedural interface
\sa itk::RealToHalfHermitianForwardFFTImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT RealToHalfHermitianForwardFFTImageFilter : public ImageFilter {
    public:
      using Self = RealToHalfHermitianForwardFFTImageFilter;

      /** Destructor */
      virtual ~RealToHalfHermitianForwardFFTImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      RealToHalfHermitianForwardFFTImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = RealPixelIDTypeList;


      /** Name of this class */
      std::string GetName() const { return std::string ("RealToHalfHermitianForwardFFTImageFilter"); }

      /** Print ourselves out */
      std::string ToString() const;


      /** Execute the filter on the input image */

      Image Execute ( const Image& image1 );

    private:

      /** Setup for member function dispatching */

      using MemberFunctionType = Image (Self::*)( const Image& image1 );
      template <class TImageType> Image ExecuteInternal ( const Image& image1 );


      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;

      std::unique_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;



    };

    /**\
     * \brief Base class for specialized real-to-complex forward Fast Fourier Transform .
     *
     * This function directly calls the execute method of RealToHalfHermitianForwardFFTImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::RealToHalfHermitianForwardFFTImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image RealToHalfHermitianForwardFFT ( const Image& image1 );

     /** @} */
  }
}
#endif