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
#ifndef sitkDiscreteGaussianDerivativeImageFilter_h
#define sitkDiscreteGaussianDerivativeImageFilter_h

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

    /**\class DiscreteGaussianDerivativeImageFilter
\brief Calculates image derivatives using discrete derivative gaussian kernels. This filter calculates Gaussian derivative by separable convolution of an image and a discrete Gaussian derivative operator (kernel).

The Gaussian operators used here were described by Tony Lindeberg (Discrete Scale-Space Theory and the Scale-Space Primal Sketch. Dissertation. Royal Institute of Technology, Stockholm, Sweden. May 1991.)

The variance or standard deviation (sigma) will be evaluated as pixel units if SetUseImageSpacing is off (false) or as physical units if SetUseImageSpacing is on (true, default). The variance can be set independently in each dimension.

When the Gaussian kernel is small, this filter tends to run faster than itk::RecursiveGaussianImageFilter .

\author Ivan Macia, VICOMTech, Spain, http://www.vicomtech.es 


This implementation was taken from the Insight Journal paper: https://hdl.handle.net/1926/1290 

\see GaussianDerivativeOperator 


\see Image 


\see Neighborhood 


\see NeighborhoodOperator
\sa itk::simple::DiscreteGaussianDerivative for the procedural interface
\sa itk::DiscreteGaussianDerivativeImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT DiscreteGaussianDerivativeImageFilter : public ImageFilter {
    public:
      using Self = DiscreteGaussianDerivativeImageFilter;

      /** Destructor */
      virtual ~DiscreteGaussianDerivativeImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      DiscreteGaussianDerivativeImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;
\

      /**
       * Convenience Set methods for setting all dimensional parameters to the same values.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetVariance ( std::vector<double> Variance ) { this->m_Variance = std::move(Variance); return *this; }

      /** Set the values of the Variance vector all to value */
      SITK_RETURN_SELF_TYPE_HEADER SetVariance( double value ) { this->m_Variance = std::vector<double>(3, value); return *this; }

      /**
       * The variance for the discrete Gaussian kernel. Sets the variance independently for each dimension, but see also SetVariance(const double v) . The default is 0.0 in each dimension. If UseImageSpacing is true, the units are the physical units of your image. If UseImageSpacing is false then the units are pixels.
       */
      std::vector<double> GetVariance() const { return this->m_Variance; }\

      /**
       * Convenience Set methods for setting all dimensional parameters to the same values.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetOrder ( std::vector<unsigned int> Order ) { this->m_Order = std::move(Order); return *this; }

      /** Set the values of the Order vector all to value */
      SITK_RETURN_SELF_TYPE_HEADER SetOrder( unsigned int value ) { this->m_Order = std::vector<unsigned int>(3, value); return *this; }

      /**
       * Order of derivatives in each dimension. Sets the derivative order independently for each dimension, but see also SetOrder(const unsigned int v) . The default is 1 in each dimension.
       */
      std::vector<unsigned int> GetOrder() const { return this->m_Order; }\

      /**
       * Set the kernel to be no wider than MaximumKernelWidth pixels, even if MaximumError demands it. The default is 32 pixels.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetMaximumKernelWidth ( unsigned int MaximumKernelWidth ) { this->m_MaximumKernelWidth = MaximumKernelWidth; return *this; }

      /**
       * Set the kernel to be no wider than MaximumKernelWidth pixels, even if MaximumError demands it. The default is 32 pixels.
       */
      unsigned int GetMaximumKernelWidth() const { return this->m_MaximumKernelWidth; }\

      /**
       * Convenience Set methods for setting all dimensional parameters to the same values.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetMaximumError ( double MaximumError ) { this->m_MaximumError = MaximumError; return *this; }

      /**
       * The algorithm will size the discrete kernel so that the error resulting from truncation of the kernel is no greater than MaximumError. The default is 0.01 in each dimension.
       */
      double GetMaximumError() const { return this->m_MaximumError; }\

      /**
       * Set/Get whether or not the filter will use the spacing of the input image in its calculations. Default is ImageSpacingOn.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetUseImageSpacing ( bool UseImageSpacing ) { this->m_UseImageSpacing = UseImageSpacing; return *this; }

      /** Set the value of UseImageSpacing to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER UseImageSpacingOn() { return this->SetUseImageSpacing(true); }
      SITK_RETURN_SELF_TYPE_HEADER UseImageSpacingOff() { return this->SetUseImageSpacing(false); }

      /**
       * Set/Get whether or not the filter will use the spacing of the input image in its calculations. Default is ImageSpacingOn.
       */
      bool GetUseImageSpacing() const { return this->m_UseImageSpacing; }\

      /**
       * Set/Get the flag for calculating scale-space normalized derivatives. Normalized derivatives are obtained multiplying by the scale parameter t.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetNormalizeAcrossScale ( bool NormalizeAcrossScale ) { this->m_NormalizeAcrossScale = NormalizeAcrossScale; return *this; }

      /** Set the value of NormalizeAcrossScale to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER NormalizeAcrossScaleOn() { return this->SetNormalizeAcrossScale(true); }
      SITK_RETURN_SELF_TYPE_HEADER NormalizeAcrossScaleOff() { return this->SetNormalizeAcrossScale(false); }

      /**
       * Set/Get the flag for calculating scale-space normalized derivatives. Normalized derivatives are obtained multiplying by the scale parameter t.
       */
      bool GetNormalizeAcrossScale() const { return this->m_NormalizeAcrossScale; }

      /** Name of this class */
      std::string GetName() const { return std::string ("DiscreteGaussianDerivativeImageFilter"); }

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


      std::vector<double>  m_Variance{std::vector<double>(3, 0.0)};

      std::vector<unsigned int>  m_Order{std::vector<unsigned int>(3, 1)};

      unsigned int  m_MaximumKernelWidth{32u};

      double  m_MaximumError{0.01};

      bool  m_UseImageSpacing{true};

      bool  m_NormalizeAcrossScale{false};


    };

    /**\
     * \brief Calculates image derivatives using discrete derivative gaussian kernels. This filter calculates Gaussian derivative by separable convolution of an image and a discrete Gaussian derivative operator (kernel).
     *
     * This function directly calls the execute method of DiscreteGaussianDerivativeImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::DiscreteGaussianDerivativeImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image DiscreteGaussianDerivative ( const Image& image1, std::vector<double> variance = std::vector<double>(3, 0.0), std::vector<unsigned int> order = std::vector<unsigned int>(3, 1), unsigned int maximumKernelWidth = 32u, double maximumError = 0.01, bool useImageSpacing = true, bool normalizeAcrossScale = false );

     /** @} */
  }
}
#endif
