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
#ifndef sitkPatchBasedDenoisingImageFilter_h
#define sitkPatchBasedDenoisingImageFilter_h

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

    /**\class PatchBasedDenoisingImageFilter
\brief Derived class implementing a specific patch-based denoising algorithm, as detailed below.

This class is derived from the base class PatchBasedDenoisingBaseImageFilter; please refer to the documentation of the base class first. This class implements a denoising filter that uses iterative non-local, or semi-local, weighted averaging of image patches for image denoising. The intensity at each pixel 'p' gets updated as a weighted average of intensities of a chosen subset of pixels from the image.

This class implements the denoising algorithm using a Gaussian kernel function for nonparametric density estimation. The class implements a scheme to automatically estimated the kernel bandwidth parameter (namely, sigma) using leave-one-out cross validation. It implements schemes for random sampling of patches non-locally (from the entire image) as well as semi-locally (from the spatial proximity of the pixel being denoised at the specific point in time). It implements a specific scheme for defining patch weights (mask) as described in Awate and Whitaker 2005 IEEE CVPR and 2006 IEEE TPAMI.

\see PatchBasedDenoisingBaseImageFilter

\sa itk::PatchBasedDenoisingImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT PatchBasedDenoisingImageFilter : public ImageFilter {
    public:
      using Self = PatchBasedDenoisingImageFilter;

      /** Destructor */
      virtual ~PatchBasedDenoisingImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      PatchBasedDenoisingImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;
\

      /**
       * Set/Get initial kernel bandwidth estimate. To prevent the class from automatically modifying this estimate, set KernelBandwidthEstimation to false in the base class.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetKernelBandwidthSigma ( double KernelBandwidthSigma ) { this->m_KernelBandwidthSigma = KernelBandwidthSigma; return *this; }

      /**
       * Set/Get initial kernel bandwidth estimate. To prevent the class from automatically modifying this estimate, set KernelBandwidthEstimation to false in the base class.
       */
      double GetKernelBandwidthSigma() const { return this->m_KernelBandwidthSigma; }\

      /**
       * Set/Get the patch radius specified in physical coordinates. Patch radius is preferably set to an even number. Currently, only isotropic patches in physical space are allowed; patches can be anisotropic in voxel space.

       */
      SITK_RETURN_SELF_TYPE_HEADER SetPatchRadius ( uint32_t PatchRadius ) { this->m_PatchRadius = PatchRadius; return *this; }

      /**
       * Set/Get the patch radius specified in physical coordinates. Patch radius is preferably set to an even number. Currently, only isotropic patches in physical space are allowed; patches can be anisotropic in voxel space.

       */
      uint32_t GetPatchRadius() const { return this->m_PatchRadius; }\

      /**
       * Set/Get the number of denoising iterations to perform. Must be a positive integer. Defaults to 1.

       */
      SITK_RETURN_SELF_TYPE_HEADER SetNumberOfIterations ( uint32_t NumberOfIterations ) { this->m_NumberOfIterations = NumberOfIterations; return *this; }

      /**
       * Set/Get the number of denoising iterations to perform. Must be a positive integer. Defaults to 1.

       */
      uint32_t GetNumberOfIterations() const { return this->m_NumberOfIterations; }\

      /**
       * Set/Get the number of patches to sample for each pixel.

       */
      SITK_RETURN_SELF_TYPE_HEADER SetNumberOfSamplePatches ( uint32_t NumberOfSamplePatches ) { this->m_NumberOfSamplePatches = NumberOfSamplePatches; return *this; }

      /**
       */
      uint32_t GetNumberOfSamplePatches() const { return this->m_NumberOfSamplePatches; }\

      /**
       * Set/Get the variance of the domain where patches are sampled.

       */
      SITK_RETURN_SELF_TYPE_HEADER SetSampleVariance ( double SampleVariance ) { this->m_SampleVariance = SampleVariance; return *this; }

      /**
       * Set/Get the variance of the domain where patches are sampled.

       */
      double GetSampleVariance() const { return this->m_SampleVariance; }

      typedef enum {NOMODEL,GAUSSIAN,RICIAN,POISSON} NoiseModelType;\

      /**
       * Set/Get the noise model type. Defaults to GAUSSIAN. To use the noise model during denoising, FidelityWeight must be positive.

       */
      SITK_RETURN_SELF_TYPE_HEADER SetNoiseModel ( NoiseModelType NoiseModel ) { this->m_NoiseModel = NoiseModel; return *this; }

      /**
       * Set/Get the noise model type. Defaults to GAUSSIAN. To use the noise model during denoising, FidelityWeight must be positive.

       */
      NoiseModelType GetNoiseModel() const { return this->m_NoiseModel; }\

      /**
       * Set/Get the noise sigma. Used by the noise model where appropriate, defaults to 5% of the image intensity range
       */
      SITK_RETURN_SELF_TYPE_HEADER SetNoiseSigma ( double NoiseSigma ) { this->m_NoiseSigma = NoiseSigma; return *this; }

      /**
       */
      double GetNoiseSigma() const { return this->m_NoiseSigma; }\

      /**
       * Set/Get the weight on the fidelity term (penalizes deviations from the noisy data). This option is used when a noise model is specified. This weight controls the balance between the smoothing and the closeness to the noisy data.

       */
      SITK_RETURN_SELF_TYPE_HEADER SetNoiseModelFidelityWeight ( double NoiseModelFidelityWeight ) { this->m_NoiseModelFidelityWeight = NoiseModelFidelityWeight; return *this; }

      /**
       * Set/Get the weight on the fidelity term (penalizes deviations from the noisy data). This option is used when a noise model is specified. This weight controls the balance between the smoothing and the closeness to the noisy data.

       */
      double GetNoiseModelFidelityWeight() const { return this->m_NoiseModelFidelityWeight; }\

      /**
       * Set/Get flag indicating whether all components should always be treated as if they are in euclidean space regardless of pixel type. Defaults to false.

       */
      SITK_RETURN_SELF_TYPE_HEADER SetAlwaysTreatComponentsAsEuclidean ( bool AlwaysTreatComponentsAsEuclidean ) { this->m_AlwaysTreatComponentsAsEuclidean = AlwaysTreatComponentsAsEuclidean; return *this; }

      /** Set the value of AlwaysTreatComponentsAsEuclidean to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER AlwaysTreatComponentsAsEuclideanOn() { return this->SetAlwaysTreatComponentsAsEuclidean(true); }
      SITK_RETURN_SELF_TYPE_HEADER AlwaysTreatComponentsAsEuclideanOff() { return this->SetAlwaysTreatComponentsAsEuclidean(false); }

      /**
       * Set/Get flag indicating whether all components should always be treated as if they are in euclidean space regardless of pixel type. Defaults to false.

       */
      bool GetAlwaysTreatComponentsAsEuclidean() const { return this->m_AlwaysTreatComponentsAsEuclidean; }\

      /**
       * Set/Get flag indicating whether kernel-bandwidth should be estimated automatically from the image data. Defaults to true.

       */
      SITK_RETURN_SELF_TYPE_HEADER SetKernelBandwidthEstimation ( bool KernelBandwidthEstimation ) { this->m_KernelBandwidthEstimation = KernelBandwidthEstimation; return *this; }

      /** Set the value of KernelBandwidthEstimation to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER KernelBandwidthEstimationOn() { return this->SetKernelBandwidthEstimation(true); }
      SITK_RETURN_SELF_TYPE_HEADER KernelBandwidthEstimationOff() { return this->SetKernelBandwidthEstimation(false); }

      /**
       * Set/Get flag indicating whether kernel-bandwidth should be estimated automatically from the image data. Defaults to true.

       */
      bool GetKernelBandwidthEstimation() const { return this->m_KernelBandwidthEstimation; }\

      /**
       * Set/Get the kernel bandwidth sigma multiplication factor used to modify the automatically-estimated kernel bandwidth sigma. At times, it may be desirable to modify the value of the automatically-estimated sigma. Typically, this number isn't very far from 1. Note: This is used only when KernelBandwidthEstimation is True/On.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetKernelBandwidthMultiplicationFactor ( double KernelBandwidthMultiplicationFactor ) { this->m_KernelBandwidthMultiplicationFactor = KernelBandwidthMultiplicationFactor; return *this; }

      /**
       * Set/Get the kernel bandwidth sigma multiplication factor used to modify the automatically-estimated kernel bandwidth sigma. At times, it may be desirable to modify the value of the automatically-estimated sigma. Typically, this number isn't very far from 1. Note: This is used only when KernelBandwidthEstimation is True/On.
       */
      double GetKernelBandwidthMultiplicationFactor() const { return this->m_KernelBandwidthMultiplicationFactor; }\

      /**
       * Set/Get the update frequency for the kernel bandwidth estimation. An optimal bandwidth will be re-estimated based on the denoised image after every 'n' iterations. Must be a positive integer. Defaults to 3, i.e. bandwidth updated after every 3 denoising iteration.

       */
      SITK_RETURN_SELF_TYPE_HEADER SetKernelBandwidthUpdateFrequency ( uint32_t KernelBandwidthUpdateFrequency ) { this->m_KernelBandwidthUpdateFrequency = KernelBandwidthUpdateFrequency; return *this; }

      /**
       * Set/Get the update frequency for the kernel bandwidth estimation. An optimal bandwidth will be re-estimated based on the denoised image after every 'n' iterations. Must be a positive integer. Defaults to 3, i.e. bandwidth updated after every 3 denoising iteration.

       */
      uint32_t GetKernelBandwidthUpdateFrequency() const { return this->m_KernelBandwidthUpdateFrequency; }\

      /**
       * Set/Get the fraction of the image to use for kernel bandwidth sigma estimation. To reduce the computational burden for computing sigma, a small random fraction of the image pixels can be used.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetKernelBandwidthFractionPixelsForEstimation ( double KernelBandwidthFractionPixelsForEstimation ) { this->m_KernelBandwidthFractionPixelsForEstimation = KernelBandwidthFractionPixelsForEstimation; return *this; }

      /**
       * Set/Get the fraction of the image to use for kernel bandwidth sigma estimation. To reduce the computational burden for computing sigma, a small random fraction of the image pixels can be used.
       */
      double GetKernelBandwidthFractionPixelsForEstimation() const { return this->m_KernelBandwidthFractionPixelsForEstimation; }

      /** Name of this class */
      std::string GetName() const { return std::string ("PatchBasedDenoisingImageFilter"); }

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


      /* KernelBandwidthSigma */
      double  m_KernelBandwidthSigma{400.0};

      uint32_t  m_PatchRadius{4u};

      /* Number of iterations to run */
      uint32_t  m_NumberOfIterations{1u};

      uint32_t  m_NumberOfSamplePatches{200u};

      double  m_SampleVariance{400.0};

      NoiseModelType  m_NoiseModel{itk::simple::PatchBasedDenoisingImageFilter::NOMODEL};

      double  m_NoiseSigma{0.0};

      double  m_NoiseModelFidelityWeight{0.0};

      bool  m_AlwaysTreatComponentsAsEuclidean{false};

      bool  m_KernelBandwidthEstimation{false};

      double  m_KernelBandwidthMultiplicationFactor{1.0};

      uint32_t  m_KernelBandwidthUpdateFrequency{3u};

      double  m_KernelBandwidthFractionPixelsForEstimation{0.2};


    };


  }
}
#endif
