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
#ifndef sitkBinShrinkImageFilter_h
#define sitkBinShrinkImageFilter_h

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

    /**\class BinShrinkImageFilter
\brief Reduce the size of an image by an integer factor in each dimension while performing averaging of an input neighborhood.

The output image size in each dimension is given by:

outputSize[j] = max( std::floor(inputSize[j]/shrinkFactor[j]), 1 );

The algorithm implemented can be describe with the following equation for 2D: \f[ \mathsf{I}_{out}(x_o,x_1) = \frac{\sum_{i=0}^{f_0}\sum_{j=0}^{f_1}\mathsf{I}_{in}(f_0 x_o+i,f_1 x_1+j)}{f_0 f_1} \f] 

This filter is implemented so that the starting extent of the first pixel of the output matches that of the input.

The change in image geometry from a 5x5 image binned by a factor of 2x2.


This code was contributed in the Insight Journal paper: "BinShrink: A multi-resolution filter with cache efficient averaging" by Lowekamp B., Chen D. https://hdl.handle.net/10380/3450
\sa itk::simple::BinShrink for the procedural interface
\sa itk::BinShrinkImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT BinShrinkImageFilter : public ImageFilter {
    public:
      using Self = BinShrinkImageFilter;

      /** Destructor */
      virtual ~BinShrinkImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      BinShrinkImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = NonLabelPixelIDTypeList;
      /** Custom public declarations */
      SITK_RETURN_SELF_TYPE_HEADER SetShrinkFactor( unsigned int s ) { this->m_ShrinkFactors = std::vector<unsigned int>(3, s ); return *this; }

\

      /**
       * Set the shrink factors. Values are clamped to a minimum value of 1. Default is 1 for all dimensions.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetShrinkFactors ( std::vector<unsigned int> ShrinkFactors ) { this->m_ShrinkFactors = std::move(ShrinkFactors); return *this; }

      /**
       * Get the shrink factors.
       */
      std::vector<unsigned int> GetShrinkFactors() const { return this->m_ShrinkFactors; }

      /** Name of this class */
      std::string GetName() const { return std::string ("BinShrinkImageFilter"); }

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


      /*  */
      std::vector<unsigned int>  m_ShrinkFactors{std::vector<unsigned int>(3, 1)};


    };

    /**\
     * \brief Reduce the size of an image by an integer factor in each dimension while performing averaging of an input neighborhood.
     *
     * This function directly calls the execute method of BinShrinkImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::BinShrinkImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image BinShrink ( const Image& image1, std::vector<unsigned int> shrinkFactors = std::vector<unsigned int>(3, 1) );

     /** @} */
  }
}
#endif
