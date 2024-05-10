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
#ifndef sitkTikhonovDeconvolutionImageFilter_h
#define sitkTikhonovDeconvolutionImageFilter_h

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

    /**\class TikhonovDeconvolutionImageFilter
\brief An inverse deconvolution filter regularized in the Tikhonov sense.

The Tikhonov deconvolution filter is the inverse deconvolution filter with a regularization term added to the denominator. The filter minimizes the equation \f[ ||\hat{f} \otimes h - g||_{L_2}^2 + \mu||\hat{f}||^2 \f] where \f$\hat{f}\f$ is the estimate of the unblurred image, \f$h\f$ is the blurring kernel, \f$g\f$ is the blurred image, and \f$\mu\f$ is a non-negative real regularization function.

The filter applies a kernel described in the Fourier domain as \f$H^*(\omega) / (|H(\omega)|^2 + \mu)\f$ where \f$H(\omega)\f$ is the Fourier transform of \f$h\f$ . The term \f$\mu\f$ is called RegularizationConstant in this filter. If \f$\mu\f$ is set to zero, this filter is equivalent to the InverseDeconvolutionImageFilter .

\author Gaetan Lehmann, Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France 


\author Cory Quammen, The University of North Carolina at Chapel Hill
\sa itk::simple::TikhonovDeconvolution for the procedural interface
\sa itk::TikhonovDeconvolutionImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT TikhonovDeconvolutionImageFilter : public ImageFilter {
    public:
      using Self = TikhonovDeconvolutionImageFilter;

      /** Destructor */
      virtual ~TikhonovDeconvolutionImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      TikhonovDeconvolutionImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;
\

      /**
       * The regularization factor. Larger values reduce the dominance of noise in the solution, but results in higher approximation error in the deblurred image. Default value is 0.0, yielding the same results as the InverseDeconvolutionImageFilter .
       */
      SITK_RETURN_SELF_TYPE_HEADER SetRegularizationConstant ( double RegularizationConstant ) { this->m_RegularizationConstant = RegularizationConstant; return *this; }

      /**
       * The regularization factor. Larger values reduce the dominance of noise in the solution, but results in higher approximation error in the deblurred image. Default value is 0.0, yielding the same results as the InverseDeconvolutionImageFilter .
       */
      double GetRegularizationConstant() const { return this->m_RegularizationConstant; }\

      /**
       * Normalize the output image by the sum of the kernel components

       */
      SITK_RETURN_SELF_TYPE_HEADER SetNormalize ( bool Normalize ) { this->m_Normalize = Normalize; return *this; }

      /** Set the value of Normalize to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER NormalizeOn() { return this->SetNormalize(true); }
      SITK_RETURN_SELF_TYPE_HEADER NormalizeOff() { return this->SetNormalize(false); }

      /**
       */
      bool GetNormalize() const { return this->m_Normalize; }

      typedef enum {ZERO_PAD,ZERO_FLUX_NEUMANN_PAD,PERIODIC_PAD} BoundaryConditionType;\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetBoundaryCondition ( BoundaryConditionType BoundaryCondition ) { this->m_BoundaryCondition = BoundaryCondition; return *this; }

      /**
       */
      BoundaryConditionType GetBoundaryCondition() const { return this->m_BoundaryCondition; }

      typedef enum {SAME,VALID} OutputRegionModeType;\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetOutputRegionMode ( OutputRegionModeType OutputRegionMode ) { this->m_OutputRegionMode = OutputRegionMode; return *this; }

      /**
       */
      OutputRegionModeType GetOutputRegionMode() const { return this->m_OutputRegionMode; }

      /** Name of this class */
      std::string GetName() const { return std::string ("TikhonovDeconvolutionImageFilter"); }

      /** Print ourselves out */
      std::string ToString() const;


      /** Execute the filter on the input images */

      Image Execute ( const Image& image1, const Image& image2 );

    private:

      /** Setup for member function dispatching */

      using MemberFunctionType = Image (Self::*)( const Image& image1, const Image& image2 );
      template <class TImageType> Image ExecuteInternal ( const Image& image1, const Image& image2 );


      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;

      std::unique_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;


      double  m_RegularizationConstant{0.0};

      /*  */
      bool  m_Normalize{false};

      BoundaryConditionType  m_BoundaryCondition{itk::simple::TikhonovDeconvolutionImageFilter::ZERO_FLUX_NEUMANN_PAD};

      OutputRegionModeType  m_OutputRegionMode{itk::simple::TikhonovDeconvolutionImageFilter::SAME};


    };

    /**\
     * \brief An inverse deconvolution filter regularized in the Tikhonov sense.
     *
     * This function directly calls the execute method of TikhonovDeconvolutionImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::TikhonovDeconvolutionImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image TikhonovDeconvolution ( const Image& image1, const Image& image2, double regularizationConstant = 0.0, bool normalize = false, TikhonovDeconvolutionImageFilter::BoundaryConditionType boundaryCondition = itk::simple::TikhonovDeconvolutionImageFilter::ZERO_FLUX_NEUMANN_PAD, TikhonovDeconvolutionImageFilter::OutputRegionModeType outputRegionMode = itk::simple::TikhonovDeconvolutionImageFilter::SAME );

     /** @} */
  }
}
#endif