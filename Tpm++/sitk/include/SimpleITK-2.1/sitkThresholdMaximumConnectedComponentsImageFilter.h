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
#ifndef sitkThresholdMaximumConnectedComponentsImageFilter_h
#define sitkThresholdMaximumConnectedComponentsImageFilter_h

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

    /**\class ThresholdMaximumConnectedComponentsImageFilter
\brief Finds the threshold value of an image based on maximizing the number of objects in the image that are larger than a given minimal size.

\par 
This method is based on Topological Stable State Thresholding to calculate the threshold set point. This method is particularly effective when there are a large number of objects in a microscopy image. Compiling in Debug mode and enable the debug flag for this filter to print debug information to see how the filter focuses in on a threshold value. Please see the Insight Journal's MICCAI 2005 workshop for a complete description. References are below.


\par Parameters
The MinimumObjectSizeInPixels parameter is controlled through the class Get/SetMinimumObjectSizeInPixels() method. Similar to the standard itk::BinaryThresholdImageFilter the Get/SetInside and Get/SetOutside values of the threshold can be set. The GetNumberOfObjects() and GetThresholdValue() methods return the number of objects above the minimum pixel size and the calculated threshold value.


\par Automatic Thresholding in ITK
There are multiple methods to automatically calculate the threshold intensity value of an image. As of version 4.0, ITK has a Thresholding ( ITKThresholding ) module which contains numerous automatic thresholding methods.implements two of these. Topological Stable State Thresholding works well on images with a large number of objects to be counted.


\par References:
1) Urish KL, August J, Huard J. "Unsupervised segmentation for myofiber
counting in immunoflourescent images". Insight Journal. ISC/NA-MIC/MICCAI Workshop on Open-Source Software (2005) Dspace handle: https://hdl.handle.net/1926/48 2) Pikaz A, Averbuch, A. "Digital image thresholding based on topological
stable-state". Pattern Recognition, 29(5): 829-843, 1996.


\par 
Questions: email Ken Urish at ken.urish(at)gmail.com Please cc the itk list serve for archival purposes.
\sa itk::simple::ThresholdMaximumConnectedComponents for the procedural interface
\sa itk::ThresholdMaximumConnectedComponentsImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT ThresholdMaximumConnectedComponentsImageFilter : public ImageFilter {
    public:
      using Self = ThresholdMaximumConnectedComponentsImageFilter;

      /** Destructor */
      virtual ~ThresholdMaximumConnectedComponentsImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      ThresholdMaximumConnectedComponentsImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = ScalarPixelIDTypeList;
\

      /**
       * The pixel type must support comparison operators. Set the minimum pixel area used to count objects on the image. Thus, only objects that have a pixel area greater than the minimum pixel area will be counted as an object in the optimization portion of this filter. Essentially, it eliminates noise from being counted as an object. The default value is zero.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetMinimumObjectSizeInPixels ( uint32_t MinimumObjectSizeInPixels ) { this->m_MinimumObjectSizeInPixels = MinimumObjectSizeInPixels; return *this; }

      /**
       * The pixel type must support comparison operators. Set the minimum pixel area used to count objects on the image. Thus, only objects that have a pixel area greater than the minimum pixel area will be counted as an object in the optimization portion of this filter. Essentially, it eliminates noise from being counted as an object. The default value is zero.
       */
      uint32_t GetMinimumObjectSizeInPixels() const { return this->m_MinimumObjectSizeInPixels; }\

      /**
       * The following Set/Get methods are for the binary threshold function. This class automatically calculates the lower threshold boundary. The upper threshold boundary, inside value, and outside value can be defined by the user, however the standard values are used as default if not set by the user. The default value of the: Inside value is the maximum pixel type intensity. Outside value is the minimum pixel type intensity. Upper threshold boundary is the maximum pixel type intensity.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetUpperBoundary ( double UpperBoundary ) { this->m_UpperBoundary = UpperBoundary; return *this; }

      /**
       * The following Set/Get methods are for the binary threshold function. This class automatically calculates the lower threshold boundary. The upper threshold boundary, inside value, and outside value can be defined by the user, however the standard values are used as default if not set by the user. The default value of the: Inside value is the maximum pixel type intensity. Outside value is the minimum pixel type intensity. Upper threshold boundary is the maximum pixel type intensity.
       */
      double GetUpperBoundary() const { return this->m_UpperBoundary; }\

      /**
       * The following Set/Get methods are for the binary threshold function. This class automatically calculates the lower threshold boundary. The upper threshold boundary, inside value, and outside value can be defined by the user, however the standard values are used as default if not set by the user. The default value of the: Inside value is the maximum pixel type intensity. Outside value is the minimum pixel type intensity. Upper threshold boundary is the maximum pixel type intensity.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetInsideValue ( uint8_t InsideValue ) { this->m_InsideValue = InsideValue; return *this; }

      /**
       * The following Set/Get methods are for the binary threshold function. This class automatically calculates the lower threshold boundary. The upper threshold boundary, inside value, and outside value can be defined by the user, however the standard values are used as default if not set by the user. The default value of the: Inside value is the maximum pixel type intensity. Outside value is the minimum pixel type intensity. Upper threshold boundary is the maximum pixel type intensity.
       */
      uint8_t GetInsideValue() const { return this->m_InsideValue; }\

      /**
       * The following Set/Get methods are for the binary threshold function. This class automatically calculates the lower threshold boundary. The upper threshold boundary, inside value, and outside value can be defined by the user, however the standard values are used as default if not set by the user. The default value of the: Inside value is the maximum pixel type intensity. Outside value is the minimum pixel type intensity. Upper threshold boundary is the maximum pixel type intensity.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetOutsideValue ( uint8_t OutsideValue ) { this->m_OutsideValue = OutsideValue; return *this; }

      /**
       * The following Set/Get methods are for the binary threshold function. This class automatically calculates the lower threshold boundary. The upper threshold boundary, inside value, and outside value can be defined by the user, however the standard values are used as default if not set by the user. The default value of the: Inside value is the maximum pixel type intensity. Outside value is the minimum pixel type intensity. Upper threshold boundary is the maximum pixel type intensity.
       */
      uint8_t GetOutsideValue() const { return this->m_OutsideValue; }

      /** Name of this class */
      std::string GetName() const { return std::string ("ThresholdMaximumConnectedComponentsImageFilter"); }

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


      uint32_t  m_MinimumObjectSizeInPixels{0u};

      double  m_UpperBoundary{std::numeric_limits<double>::max()};

      uint8_t  m_InsideValue{1u};

      uint8_t  m_OutsideValue{0u};


    };

    /**\
     * \brief Finds the threshold value of an image based on maximizing the number of objects in the image that are larger than a given minimal size.
     *
     * This function directly calls the execute method of ThresholdMaximumConnectedComponentsImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::ThresholdMaximumConnectedComponentsImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image ThresholdMaximumConnectedComponents ( const Image& image1, uint32_t minimumObjectSizeInPixels = 0u, double upperBoundary = std::numeric_limits<double>::max(), uint8_t insideValue = 1u, uint8_t outsideValue = 0u );

     /** @} */
  }
}
#endif