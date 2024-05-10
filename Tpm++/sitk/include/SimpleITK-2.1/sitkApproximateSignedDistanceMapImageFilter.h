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
#ifndef sitkApproximateSignedDistanceMapImageFilter_h
#define sitkApproximateSignedDistanceMapImageFilter_h

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

    /**\class ApproximateSignedDistanceMapImageFilter
\brief Create a map of the approximate signed distance from the boundaries of a binary image.

The ApproximateSignedDistanceMapImageFilter takes as input a binary image and produces a signed distance map. Each pixel value in the output contains the approximate distance from that pixel to the nearest "object" in the binary image. This filter differs from the DanielssonDistanceMapImageFilter in that it calculates the distance to the "object edge" for pixels within the object.

Negative values in the output indicate that the pixel at that position is within an object in the input image. The absolute value of a negative pixel represents the approximate distance to the nearest object boundary pixel.

WARNING: This filter requires that the output type be floating-point. Otherwise internal calculations will not be performed to the appropriate precision, resulting in completely incorrect (read: zero-valued) output.

The distances computed by this filter are Chamfer distances, which are only an approximation to Euclidian distances, and are not as exact approximations as those calculated by the DanielssonDistanceMapImageFilter . On the other hand, this filter is faster.

This filter requires that an "inside value" and "outside value" be set as parameters. The "inside value" is the intensity value of the binary image which corresponds to objects, and the "outside value" is the intensity of the background. (A typical binary image often represents objects as black (0) and background as white (usually 255), or vice-versa.) Note that this filter is slightly faster if the inside value is less than the outside value. Otherwise an extra iteration through the image is required.

This filter uses the FastChamferDistanceImageFilter and the IsoContourDistanceImageFilter internally to perform the distance calculations.

\see DanielssonDistanceMapImageFilter 


\see SignedDanielssonDistanceMapImageFilter 


\see SignedMaurerDistanceMapImageFilter 


\see FastChamferDistanceImageFilter 


\see IsoContourDistanceImageFilter 


\author Zach Pincus
\sa itk::simple::ApproximateSignedDistanceMap for the procedural interface
\sa itk::ApproximateSignedDistanceMapImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT ApproximateSignedDistanceMapImageFilter : public ImageFilter {
    public:
      using Self = ApproximateSignedDistanceMapImageFilter;

      /** Destructor */
      virtual ~ApproximateSignedDistanceMapImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      ApproximateSignedDistanceMapImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = IntegerPixelIDTypeList;
\

      /**
       * Set/Get intensity value representing the interior of objects in the mask.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetInsideValue ( double InsideValue ) { this->m_InsideValue = InsideValue; return *this; }

      /**
       * Set/Get intensity value representing the interior of objects in the mask.
       */
      double GetInsideValue() const { return this->m_InsideValue; }\

      /**
       * Set/Get intensity value representing non-objects in the mask.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetOutsideValue ( double OutsideValue ) { this->m_OutsideValue = OutsideValue; return *this; }

      /**
       * Set/Get intensity value representing the interior of objects in the mask.
       */
      double GetOutsideValue() const { return this->m_OutsideValue; }

      /** Name of this class */
      std::string GetName() const { return std::string ("ApproximateSignedDistanceMapImageFilter"); }

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


      double  m_InsideValue{1u};

      double  m_OutsideValue{0u};


    };

    /**\
     * \brief Create a map of the approximate signed distance from the boundaries of a binary image.
     *
     * This function directly calls the execute method of ApproximateSignedDistanceMapImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::ApproximateSignedDistanceMapImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image ApproximateSignedDistanceMap ( const Image& image1, double insideValue = 1u, double outsideValue = 0u );

     /** @} */
  }
}
#endif
