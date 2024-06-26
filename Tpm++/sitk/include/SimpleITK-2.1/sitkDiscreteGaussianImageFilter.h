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
#ifndef sitkDiscreteGaussianImageFilter_h
#define sitkDiscreteGaussianImageFilter_h

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

    /**\class DiscreteGaussianImageFilter
\brief Blurs an image by separable convolution with discrete gaussian kernels. This filter performs Gaussian blurring by separable convolution of an image and a discrete Gaussian operator (kernel).

The Gaussian operator used here was described by Tony Lindeberg (Discrete Scale-Space Theory and the Scale-Space Primal Sketch. Dissertation. Royal Institute of Technology, Stockholm, Sweden. May 1991.) The Gaussian kernel used here was designed so that smoothing and derivative operations commute after discretization.

The variance or standard deviation (sigma) will be evaluated as pixel units if SetUseImageSpacing is off (false) or as physical units if SetUseImageSpacing is on (true, default). The variance can be set independently in each dimension.

When the Gaussian kernel is small, this filter tends to run faster than itk::RecursiveGaussianImageFilter .

\see GaussianOperator 


\see Image 


\see Neighborhood 


\see NeighborhoodOperator 


\see RecursiveGaussianImageFilter
\sa itk::simple::DiscreteGaussian for the procedural interface
\sa itk::DiscreteGaussianImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT DiscreteGaussianImageFilter : public ImageFilter {
    public:
      using Self = DiscreteGaussianImageFilter;

      /** Destructor */
      virtual ~DiscreteGaussianImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      DiscreteGaussianImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;
\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetVariance ( std::vector<double> Variance ) { this->m_Variance = std::move(Variance); return *this; }

      /** Set the values of the Variance vector all to value */
      SITK_RETURN_SELF_TYPE_HEADER SetVariance( double value ) { this->m_Variance = std::vector<double>(3, value); return *this; }

      /**
       * The variance for the discrete Gaussian kernel. Sets the variance independently for each dimension, but see also SetVariance(const double v) . The default is 0.0 in each dimension. If UseImageSpacing is true, the units are the physical units of your image. If UseImageSpacing is false then the units are pixels.
       */
      std::vector<double> GetVariance() const { return this->m_Variance; }\

      /**
       * Set the kernel to be no wider than MaximumKernelWidth pixels, even if MaximumError demands it. The default is 32 pixels.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetMaximumKernelWidth ( unsigned int MaximumKernelWidth ) { this->m_MaximumKernelWidth = MaximumKernelWidth; return *this; }

      /**
       * Set the kernel to be no wider than MaximumKernelWidth pixels, even if MaximumError demands it. The default is 32 pixels.
       */
      unsigned int GetMaximumKernelWidth() const { return this->m_MaximumKernelWidth; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetMaximumError ( std::vector<double> MaximumError ) { this->m_MaximumError = std::move(MaximumError); return *this; }

      /** Set the values of the MaximumError vector all to value */
      SITK_RETURN_SELF_TYPE_HEADER SetMaximumError( double value ) { this->m_MaximumError = std::vector<double>(3, value); return *this; }

      /**
       * The algorithm will size the discrete kernel so that the error resulting from truncation of the kernel is no greater than MaximumError. The default is 0.01 in each dimension.
       */
      std::vector<double> GetMaximumError() const { return this->m_MaximumError; }\

      /**
       * Set/Get whether or not the filter will use the spacing of the input image in its calculations
       */
      SITK_RETURN_SELF_TYPE_HEADER SetUseImageSpacing ( bool UseImageSpacing ) { this->m_UseImageSpacing = UseImageSpacing; return *this; }

      /** Set the value of UseImageSpacing to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER UseImageSpacingOn() { return this->SetUseImageSpacing(true); }
      SITK_RETURN_SELF_TYPE_HEADER UseImageSpacingOff() { return this->SetUseImageSpacing(false); }

      /**
       * Set/Get whether or not the filter will use the spacing of the input image in its calculations
       */
      bool GetUseImageSpacing() const { return this->m_UseImageSpacing; }

      /** Name of this class */
      std::string GetName() const { return std::string ("DiscreteGaussianImageFilter"); }

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


      std::vector<double>  m_Variance{std::vector<double>(3,1.0)};

      unsigned int  m_MaximumKernelWidth{32u};

      std::vector<double>  m_MaximumError{std::vector<double>(3, 0.01)};

      bool  m_UseImageSpacing{true};


    };

    /**\
     * \brief Blurs an image by separable convolution with discrete gaussian kernels. This filter performs Gaussian blurring by separable convolution of an image and a discrete Gaussian operator (kernel).
     *
     * This function directly calls the execute method of DiscreteGaussianImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::DiscreteGaussianImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image DiscreteGaussian ( const Image& image1, std::vector<double> variance = std::vector<double>(3,1.0), unsigned int maximumKernelWidth = 32u, std::vector<double> maximumError = std::vector<double>(3, 0.01), bool useImageSpacing = true );

     /** @} */
  }
}
#endif
