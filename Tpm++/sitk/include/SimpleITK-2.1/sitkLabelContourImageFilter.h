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
#ifndef sitkLabelContourImageFilter_h
#define sitkLabelContourImageFilter_h

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

    /**\class LabelContourImageFilter
\brief Labels the pixels on the border of the objects in a labeled image.

LabelContourImageFilter takes a labeled image as input, where the pixels in the objects are the pixels with a value different of the BackgroundValue. Only the pixels on the contours of the objects are kept. The pixels not on the border are changed to BackgroundValue. The labels of the object are the same in the input and in the output image.

The connectivity can be changed to minimum or maximum connectivity with SetFullyConnected() . Full connectivity produces thicker contours.

https://hdl.handle.net/1926/1352 

\author Gaetan Lehmann. Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France.


\see BinaryContourImageFilter
\sa itk::simple::LabelContour for the procedural interface
\sa itk::LabelContourImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT LabelContourImageFilter : public ImageFilter {
    public:
      using Self = LabelContourImageFilter;

      /** Destructor */
      virtual ~LabelContourImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      LabelContourImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = IntegerPixelIDTypeList;
\

      /**
       * Set/Get whether the connected components are defined strictly by face connectivity or by face+edge+vertex connectivity. Default is FullyConnectedOff. \note For objects that are 1 pixel wide, use FullyConnectedOn.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetFullyConnected ( bool FullyConnected ) { this->m_FullyConnected = FullyConnected; return *this; }

      /** Set the value of FullyConnected to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER FullyConnectedOn() { return this->SetFullyConnected(true); }
      SITK_RETURN_SELF_TYPE_HEADER FullyConnectedOff() { return this->SetFullyConnected(false); }

      /**
       * Set/Get whether the connected components are defined strictly by face connectivity or by face+edge+vertex connectivity. Default is FullyConnectedOff. \note For objects that are 1 pixel wide, use FullyConnectedOn.
       */
      bool GetFullyConnected() const { return this->m_FullyConnected; }\

      /**
       * Set/Get the background value used to identify the objects and mark the pixels not on the border of the objects.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetBackgroundValue ( double BackgroundValue ) { this->m_BackgroundValue = BackgroundValue; return *this; }

      /**
       * Set/Get the background value used to identify the objects and mark the pixels not on the border of the objects.
       */
      double GetBackgroundValue() const { return this->m_BackgroundValue; }

      /** Name of this class */
      std::string GetName() const { return std::string ("LabelContourImageFilter"); }

      /** Print ourselves out */
      std::string ToString() const;


      /** Execute the filter on the input image */
#ifndef SWIG
      Image Execute ( Image&& image1 );
#endif
      Image Execute ( const Image& image1 );

    private:

      /** Setup for member function dispatching */

      using MemberFunctionType = Image (Self::*)( const Image& image1 );
      template <class TImageType> Image ExecuteInternal ( const Image& image1 );


      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;

      std::unique_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;


      /*  */
      bool  m_FullyConnected{false};

      double  m_BackgroundValue{0};


      bool m_InPlace{false};
    };

    /**\
     * \brief Labels the pixels on the border of the objects in a labeled image.
     *
     * This function directly calls the execute method of LabelContourImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::LabelContourImageFilter for the object oriented interface
     * @{
     */
#ifndef SWIG
     SITKBasicFilters_EXPORT Image LabelContour ( Image&& image1, bool fullyConnected = false, double backgroundValue = 0 );
#endif
     SITKBasicFilters_EXPORT Image LabelContour ( const Image& image1, bool fullyConnected = false, double backgroundValue = 0 );

     /** @} */
  }
}
#endif
