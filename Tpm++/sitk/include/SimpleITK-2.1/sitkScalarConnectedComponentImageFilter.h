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
#ifndef sitkScalarConnectedComponentImageFilter_h
#define sitkScalarConnectedComponentImageFilter_h

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

    /**\class ScalarConnectedComponentImageFilter
\brief A connected components filter that labels the objects in an arbitrary image. Two pixels are similar if they are within threshold of each other. Uses ConnectedComponentFunctorImageFilter .


\sa itk::simple::ScalarConnectedComponent for the procedural interface
\sa itk::ScalarConnectedComponentImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT ScalarConnectedComponentImageFilter : public ImageFilter {
    public:
      using Self = ScalarConnectedComponentImageFilter;

      /** Destructor */
      virtual ~ScalarConnectedComponentImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      ScalarConnectedComponentImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;
\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetDistanceThreshold ( double DistanceThreshold ) { this->m_DistanceThreshold = DistanceThreshold; return *this; }

      /**
       */
      double GetDistanceThreshold() const { return this->m_DistanceThreshold; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetFullyConnected( bool FullyConnected ) { this->m_FullyConnected = FullyConnected; return *this; }

      /** Set the value of FullyConnected to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER FullyConnectedOn() { return this->SetFullyConnected(true); }
      SITK_RETURN_SELF_TYPE_HEADER FullyConnectedOff() { return this->SetFullyConnected(false); }

      /**
       */
      bool GetFullyConnected() const { return this->m_FullyConnected; }

      /** Name of this class */
      std::string GetName() const { return std::string ("ScalarConnectedComponentImageFilter"); }

      /** Print ourselves out */
      std::string ToString() const;


      /** Execute the filter on the input image */

      Image Execute ( const Image & image, const Image & maskImage );
      Image Execute ( const Image & image );

    private:

      /** Setup for member function dispatching */

      using MemberFunctionType = Image (Self::*)( const Image * image, const Image * maskImage );
      template <class TImageType> Image ExecuteInternal ( const Image * image, const Image * maskImage );


      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;

      std::unique_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;


      double  m_DistanceThreshold{0.0};

      /*  */
      bool  m_FullyConnected{false};


    };

    /**\
     * \brief A connected components filter that labels the objects in an arbitrary image. Two pixels are similar if they are within threshold of each other. Uses ConnectedComponentFunctorImageFilter .
     *
     * This function directly calls the execute method of ScalarConnectedComponentImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::ScalarConnectedComponentImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image ScalarConnectedComponent ( const Image & image, const Image & maskImage, double distanceThreshold = 0.0, bool fullyConnected = false );
     SITKBasicFilters_EXPORT Image ScalarConnectedComponent ( const Image & image, double distanceThreshold = 0.0, bool fullyConnected = false );
     /** @} */
  }
}
#endif
