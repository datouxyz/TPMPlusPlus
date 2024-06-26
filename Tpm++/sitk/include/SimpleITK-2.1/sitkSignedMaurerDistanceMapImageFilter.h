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
#ifndef sitkSignedMaurerDistanceMapImageFilter_h
#define sitkSignedMaurerDistanceMapImageFilter_h

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

    /**\class SignedMaurerDistanceMapImageFilter
\brief This filter calculates the Euclidean distance transform of a binary image in linear time for arbitrary dimensions.

\par Inputs and Outputs
This is an image-to-image filter. The dimensionality is arbitrary. The only dimensionality constraint is that the input and output images be of the same dimensions and size. To maintain integer arithmetic within the filter, the default output is the signed squared distance. This implies that the input image should be of type "unsigned int" or "int" whereas the output image is of type "int". Obviously, if the user wishes to utilize the image spacing or to have a filter with the Euclidean distance (as opposed to the squared distance), output image types of float or double should be used.


The inside is considered as having negative distances. Outside is treated as having positive distances. To change the convention, use the InsideIsPositive(bool) function.

\par Parameters
Set/GetBackgroundValue specifies the background of the value of the input binary image. Normally this is zero and, as such, zero is the default value. Other than that, the usage is completely analogous to the itk::DanielssonDistanceImageFilter class except it does not return the Voronoi map.


Reference: C. R. Maurer, Jr., R. Qi, and V. Raghavan, "A Linear Time Algorithm for Computing Exact Euclidean Distance Transforms of Binary Images in Arbitrary Dimensions", IEEE - Transactions on Pattern Analysis and Machine Intelligence, 25(2): 265-270, 2003.
\sa itk::simple::SignedMaurerDistanceMap for the procedural interface
\sa itk::SignedMaurerDistanceMapImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT SignedMaurerDistanceMapImageFilter : public ImageFilter {
    public:
      using Self = SignedMaurerDistanceMapImageFilter;

      /** Destructor */
      virtual ~SignedMaurerDistanceMapImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      SignedMaurerDistanceMapImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = IntegerPixelIDTypeList;
\

      /**
       * Set if the inside represents positive values in the signed distance map. By convention ON pixels are treated as inside pixels.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetInsideIsPositive ( bool InsideIsPositive ) { this->m_InsideIsPositive = InsideIsPositive; return *this; }

      /** Set the value of InsideIsPositive to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER InsideIsPositiveOn() { return this->SetInsideIsPositive(true); }
      SITK_RETURN_SELF_TYPE_HEADER InsideIsPositiveOff() { return this->SetInsideIsPositive(false); }

      /**
       * Get if the inside represents positive values in the signed distance map. \see GetInsideIsPositive()
       */
      bool GetInsideIsPositive() const { return this->m_InsideIsPositive; }\

      /**
       * Set if the distance should be squared.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetSquaredDistance ( bool SquaredDistance ) { this->m_SquaredDistance = SquaredDistance; return *this; }

      /** Set the value of SquaredDistance to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER SquaredDistanceOn() { return this->SetSquaredDistance(true); }
      SITK_RETURN_SELF_TYPE_HEADER SquaredDistanceOff() { return this->SetSquaredDistance(false); }

      /**
       * Get the distance squared.
       */
      bool GetSquaredDistance() const { return this->m_SquaredDistance; }\

      /**
       * Set if image spacing should be used in computing distances.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetUseImageSpacing ( bool UseImageSpacing ) { this->m_UseImageSpacing = UseImageSpacing; return *this; }

      /** Set the value of UseImageSpacing to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER UseImageSpacingOn() { return this->SetUseImageSpacing(true); }
      SITK_RETURN_SELF_TYPE_HEADER UseImageSpacingOff() { return this->SetUseImageSpacing(false); }

      /**
       * Get whether spacing is used.
       */
      bool GetUseImageSpacing() const { return this->m_UseImageSpacing; }\

      /**
       * Set the background value which defines the object. Usually this value is = 0.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetBackgroundValue ( double BackgroundValue ) { this->m_BackgroundValue = BackgroundValue; return *this; }

      /**
       * Set the background value which defines the object. Usually this value is = 0.
       */
      double GetBackgroundValue() const { return this->m_BackgroundValue; }

      /** Name of this class */
      std::string GetName() const { return std::string ("SignedMaurerDistanceMapImageFilter"); }

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


      bool  m_InsideIsPositive{false};

      bool  m_SquaredDistance{true};

      bool  m_UseImageSpacing{false};

      double  m_BackgroundValue{0.0};


    };

    /**\
     * \brief This filter calculates the Euclidean distance transform of a binary image in linear time for arbitrary dimensions.
     *
     * This function directly calls the execute method of SignedMaurerDistanceMapImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::SignedMaurerDistanceMapImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image SignedMaurerDistanceMap ( const Image& image1, bool insideIsPositive = false, bool squaredDistance = true, bool useImageSpacing = false, double backgroundValue = 0.0 );

     /** @} */
  }
}
#endif
