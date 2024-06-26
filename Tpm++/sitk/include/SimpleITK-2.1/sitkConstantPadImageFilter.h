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
#ifndef sitkConstantPadImageFilter_h
#define sitkConstantPadImageFilter_h

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

    /**\class ConstantPadImageFilter
\brief Increase the image size by padding with a constant value.

ConstantPadImageFilter changes the output image region. If the output image region is larger than the input image region, the extra pixels are filled in by a constant value. The output image region must be specified.

Visual explanation of padding regions.


This filter is implemented as a multithreaded filter. It provides a DynamicThreadedGenerateData() method for its implementation.

\see WrapPadImageFilter , MirrorPadImageFilter
\sa itk::simple::ConstantPad for the procedural interface
\sa itk::ConstantPadImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT ConstantPadImageFilter : public ImageFilter {
    public:
      using Self = ConstantPadImageFilter;

      /** Destructor */
      virtual ~ConstantPadImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      ConstantPadImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;
\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetPadLowerBound ( std::vector<unsigned int> PadLowerBound ) { this->m_PadLowerBound = std::move(PadLowerBound); return *this; }

      /**
       */
      std::vector<unsigned int> GetPadLowerBound() const { return this->m_PadLowerBound; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetPadUpperBound ( std::vector<unsigned int> PadUpperBound ) { this->m_PadUpperBound = std::move(PadUpperBound); return *this; }

      /**
       */
      std::vector<unsigned int> GetPadUpperBound() const { return this->m_PadUpperBound; }\

      /**
       * Set/Get the pad value. Default is Zero.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetConstant ( double Constant ) { this->m_Constant = Constant; return *this; }

      /**
       * Set/Get the pad value. Default is Zero.
       */
      double GetConstant() const { return this->m_Constant; }

      /** Name of this class */
      std::string GetName() const { return std::string ("ConstantPadImageFilter"); }

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


      /* 	odo what does this itk_type attribute do? */
      std::vector<unsigned int>  m_PadLowerBound{std::vector<unsigned int>(3, 0)};

      /* 	odo what does this itk_type attribute do? */
      std::vector<unsigned int>  m_PadUpperBound{std::vector<unsigned int>(3, 0)};

      /*  */
      double  m_Constant{0.0};


    };

    /**\
     * \brief Increase the image size by padding with a constant value.
     *
     * This function directly calls the execute method of ConstantPadImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::ConstantPadImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image ConstantPad ( const Image& image1, std::vector<unsigned int> padLowerBound = std::vector<unsigned int>(3, 0), std::vector<unsigned int> padUpperBound = std::vector<unsigned int>(3, 0), double constant = 0.0 );

     /** @} */
  }
}
#endif
