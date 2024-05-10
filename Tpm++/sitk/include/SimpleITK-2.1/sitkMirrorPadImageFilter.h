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
#ifndef sitkMirrorPadImageFilter_h
#define sitkMirrorPadImageFilter_h

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

    /**\class MirrorPadImageFilter
\brief Increase the image size by padding with replicants of the input image value.

MirrorPadImageFilter changes the image bounds of an image. Any added pixels are filled in with a mirrored replica of the input image. For instance, if the output image needs a pixel that is two pixels to the left of the LargestPossibleRegion of the input image, the value assigned will be from the pixel two pixels inside the left boundary of the LargestPossibleRegion. The image bounds of the output must be specified.

Visual explanation of padding regions.


This filter is implemented as a multithreaded filter. It provides a DynamicThreadedGenerateData() method for its implementation.

Exponential decay in the bounds is enabled when DecayBase has to be in the range (0.0, 1.0]. When it is 1.0 it is disabled. The decay rate is based on the Manhattan distance.

\see WrapPadImageFilter , ConstantPadImageFilter
\sa itk::simple::MirrorPad for the procedural interface
\sa itk::MirrorPadImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT MirrorPadImageFilter : public ImageFilter {
    public:
      using Self = MirrorPadImageFilter;

      /** Destructor */
      virtual ~MirrorPadImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      MirrorPadImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = NonLabelPixelIDTypeList;
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
       * Get/Set the base for exponential decay in mirrored region.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetDecayBase ( double DecayBase ) { this->m_DecayBase = DecayBase; return *this; }

      /**
       * Get/Set the base for exponential decay in mirrored region.
       */
      double GetDecayBase() const { return this->m_DecayBase; }

      /** Name of this class */
      std::string GetName() const { return std::string ("MirrorPadImageFilter"); }

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

      double  m_DecayBase{1.0};


    };

    /**\
     * \brief Increase the image size by padding with replicants of the input image value.
     *
     * This function directly calls the execute method of MirrorPadImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::MirrorPadImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image MirrorPad ( const Image& image1, std::vector<unsigned int> padLowerBound = std::vector<unsigned int>(3, 0), std::vector<unsigned int> padUpperBound = std::vector<unsigned int>(3, 0), double decayBase = 1.0 );

     /** @} */
  }
}
#endif