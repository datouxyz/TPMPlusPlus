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
#ifndef sitkFastSymmetricForcesDemonsRegistrationFilter_h
#define sitkFastSymmetricForcesDemonsRegistrationFilter_h

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

    /**\class FastSymmetricForcesDemonsRegistrationFilter
\brief Deformably register two images using a symmetric forces demons algorithm.

This class was contributed by Tom Vercauteren, INRIA & Mauna Kea Technologies based on a variation of the DemonsRegistrationFilter .

FastSymmetricForcesDemonsRegistrationFilter implements the demons deformable algorithm that register two images by computing the deformation field which will map a moving image onto a fixed image.

A deformation field is represented as a image whose pixel type is some vector type with at least N elements, where N is the dimension of the fixed image. The vector type must support element access via operator []. It is assumed that the vector elements behave like floating point scalars.

This class is templated over the fixed image type, moving image type and the deformation field type.

The input fixed and moving images are set via methods SetFixedImage and SetMovingImage respectively. An initial deformation field maybe set via SetInitialDisplacementField or SetInput. If no initial field is set, a zero field is used as the initial condition.

The output deformation field can be obtained via methods GetOutput or GetDisplacementField.

This class make use of the finite difference solver hierarchy. Update for each iteration is computed in DemonsRegistrationFunction .

\author Tom Vercauteren, INRIA & Mauna Kea Technologies


This implementation was taken from the Insight Journal paper: https://hdl.handle.net/1926/510 

\warning This filter assumes that the fixed image type, moving image type and deformation field type all have the same number of dimensions.


\see DemonsRegistrationFilter 


\see DemonsRegistrationFunction

\sa itk::FastSymmetricForcesDemonsRegistrationFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT FastSymmetricForcesDemonsRegistrationFilter : public ImageFilter {
    public:
      using Self = FastSymmetricForcesDemonsRegistrationFilter;

      /** Destructor */
      virtual ~FastSymmetricForcesDemonsRegistrationFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      FastSymmetricForcesDemonsRegistrationFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;
\

      /**
       * Set/Get the Gaussian smoothing standard deviations for the displacement field. The values are set with respect to pixel coordinates.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetStandardDeviations ( std::vector<double> StandardDeviations ) { this->m_StandardDeviations = std::move(StandardDeviations); return *this; }

      /** Set the values of the StandardDeviations vector all to value */
      SITK_RETURN_SELF_TYPE_HEADER SetStandardDeviations( double value ) { this->m_StandardDeviations = std::vector<double>(3, value); return *this; }

      /**
       * Set/Get the Gaussian smoothing standard deviations for the displacement field. The values are set with respect to pixel coordinates.
       */
      std::vector<double> GetStandardDeviations() const { return this->m_StandardDeviations; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetNumberOfIterations ( uint32_t NumberOfIterations ) { this->m_NumberOfIterations = NumberOfIterations; return *this; }

      /**
       */
      uint32_t GetNumberOfIterations() const { return this->m_NumberOfIterations; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetMaximumRMSError ( double MaximumRMSError ) { this->m_MaximumRMSError = MaximumRMSError; return *this; }

      /**
       */
      double GetMaximumRMSError() const { return this->m_MaximumRMSError; }

      typedef enum {Symmetric,Fixed,WarpedMoving,MappedMoving} UseGradientTypeType;\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetUseGradientType ( UseGradientTypeType UseGradientType ) { this->m_UseGradientType = UseGradientType; return *this; }

      /**
       */
      UseGradientTypeType GetUseGradientType() const { return this->m_UseGradientType; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetMaximumUpdateStepLength ( double MaximumUpdateStepLength ) { this->m_MaximumUpdateStepLength = MaximumUpdateStepLength; return *this; }

      /**
       */
      double GetMaximumUpdateStepLength() const { return this->m_MaximumUpdateStepLength; }\

      /**
       * Set/Get whether the displacement field is smoothed (regularized). Smoothing the displacement yields a solution elastic in nature. If SmoothDisplacementField is on, then the displacement field is smoothed with a Gaussian whose standard deviations are specified with SetStandardDeviations()
       */
      SITK_RETURN_SELF_TYPE_HEADER SetSmoothDisplacementField ( bool SmoothDisplacementField ) { this->m_SmoothDisplacementField = SmoothDisplacementField; return *this; }

      /** Set the value of SmoothDisplacementField to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER SmoothDisplacementFieldOn() { return this->SetSmoothDisplacementField(true); }
      SITK_RETURN_SELF_TYPE_HEADER SmoothDisplacementFieldOff() { return this->SetSmoothDisplacementField(false); }

      /**
       * Set/Get whether the displacement field is smoothed (regularized). Smoothing the displacement yields a solution elastic in nature. If SmoothDisplacementField is on, then the displacement field is smoothed with a Gaussian whose standard deviations are specified with SetStandardDeviations()
       */
      bool GetSmoothDisplacementField() const { return this->m_SmoothDisplacementField; }\

      /**
       * Set/Get whether the update field is smoothed (regularized). Smoothing the update field yields a solution viscous in nature. If SmoothUpdateField is on, then the update field is smoothed with a Gaussian whose standard deviations are specified with SetUpdateFieldStandardDeviations()
       */
      SITK_RETURN_SELF_TYPE_HEADER SetSmoothUpdateField ( bool SmoothUpdateField ) { this->m_SmoothUpdateField = SmoothUpdateField; return *this; }

      /** Set the value of SmoothUpdateField to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER SmoothUpdateFieldOn() { return this->SetSmoothUpdateField(true); }
      SITK_RETURN_SELF_TYPE_HEADER SmoothUpdateFieldOff() { return this->SetSmoothUpdateField(false); }

      /**
       * Set/Get whether the update field is smoothed (regularized). Smoothing the update field yields a solution viscous in nature. If SmoothUpdateField is on, then the update field is smoothed with a Gaussian whose standard deviations are specified with SetUpdateFieldStandardDeviations()
       */
      bool GetSmoothUpdateField() const { return this->m_SmoothUpdateField; }\

      /**
       * Set the Gaussian smoothing standard deviations for the update field. The values are set with respect to pixel coordinates.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetUpdateFieldStandardDeviations ( std::vector<double> UpdateFieldStandardDeviations ) { this->m_UpdateFieldStandardDeviations = std::move(UpdateFieldStandardDeviations); return *this; }

      /** Set the values of the UpdateFieldStandardDeviations vector all to value */
      SITK_RETURN_SELF_TYPE_HEADER SetUpdateFieldStandardDeviations( double value ) { this->m_UpdateFieldStandardDeviations = std::vector<double>(3, value); return *this; }

      /**
       * Set the Gaussian smoothing standard deviations for the update field. The values are set with respect to pixel coordinates.
       */
      std::vector<double> GetUpdateFieldStandardDeviations() const { return this->m_UpdateFieldStandardDeviations; }\

      /**
       * Set/Get the desired limits of the Gaussian kernel width.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetMaximumKernelWidth ( unsigned int MaximumKernelWidth ) { this->m_MaximumKernelWidth = MaximumKernelWidth; return *this; }

      /**
       * Set/Get the desired limits of the Gaussian kernel width.
       */
      unsigned int GetMaximumKernelWidth() const { return this->m_MaximumKernelWidth; }\

      /**
       * Set/Get the desired maximum error of the Guassian kernel approximate.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetMaximumError ( double MaximumError ) { this->m_MaximumError = MaximumError; return *this; }

      /**
       * Set/Get the desired maximum error of the Guassian kernel approximate.
       */
      double GetMaximumError() const { return this->m_MaximumError; }\

      /**
       * Set/Get the threshold below which the absolute difference of intensity yields a match. When the intensities match between a moving and fixed image pixel, the update vector (for that iteration) will be the zero vector. Default is 0.001.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetIntensityDifferenceThreshold ( double IntensityDifferenceThreshold ) { this->m_IntensityDifferenceThreshold = IntensityDifferenceThreshold; return *this; }

      /**
       */
      double GetIntensityDifferenceThreshold() const { return this->m_IntensityDifferenceThreshold; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetUseImageSpacing ( bool UseImageSpacing ) { this->m_UseImageSpacing = UseImageSpacing; return *this; }

      /** Set the value of UseImageSpacing to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER UseImageSpacingOn() { return this->SetUseImageSpacing(true); }
      SITK_RETURN_SELF_TYPE_HEADER UseImageSpacingOff() { return this->SetUseImageSpacing(false); }

      /**
       */
      bool GetUseImageSpacing() const { return this->m_UseImageSpacing; }
     /** \brief Number of iterations run.
      *
      *
      * This is an active measurement. It may be accessed while the
      * filter is being executing in command call-backs and can be
      * accessed after execution.
      */
     uint32_t GetElapsedIterations() const { return this->m_pfGetElapsedIterations(); };

     /**
      * Set/Get the root mean squared change of the previous iteration. May not be used by all solvers.
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     double GetRMSChange() const { return this->m_RMSChange; };

     /**
      * Get the metric value. The metric value is the mean square difference in intensity between the fixed image and transforming moving image computed over the the overlapping region between the two images. This value is calculated for the current iteration
      *
      * This is an active measurement. It may be accessed while the
      * filter is being executing in command call-backs and can be
      * accessed after execution.
      */
     double GetMetric() const { return this->m_pfGetMetric(); };


      /** Name of this class */
      std::string GetName() const { return std::string ("FastSymmetricForcesDemonsRegistrationFilter"); }

      /** Print ourselves out */
      std::string ToString() const;


      /** Execute the filter on the input image */

      Image Execute ( const Image & fixedImage, const Image & movingImage, const Image & initialDisplacementField );
      Image Execute ( const Image & fixedImage, const Image & movingImage );

    private:

      /** Setup for member function dispatching */

      using MemberFunctionType = Image (Self::*)( const Image * fixedImage, const Image * movingImage, const Image * initialDisplacementField );
      template <class TImageType> Image ExecuteInternal ( const Image * fixedImage, const Image * movingImage, const Image * initialDisplacementField );


      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;

      std::unique_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;


      std::vector<double>  m_StandardDeviations{std::vector<double>(3, 1.0)};

      /* Number of iterations to run */
      uint32_t  m_NumberOfIterations{10u};

      /* Value of RMS change below which the filter should stop. This is a convergence criterion. */
      double  m_MaximumRMSError{0.02};

      UseGradientTypeType  m_UseGradientType{itk::simple::FastSymmetricForcesDemonsRegistrationFilter::Symmetric};

      double  m_MaximumUpdateStepLength{0.5};

      bool  m_SmoothDisplacementField{true};

      bool  m_SmoothUpdateField{false};

      std::vector<double>  m_UpdateFieldStandardDeviations{std::vector<double>(3, 1.0)};

      unsigned int  m_MaximumKernelWidth{30u};

      double  m_MaximumError{0.1};

      double  m_IntensityDifferenceThreshold{0.001};

      bool  m_UseImageSpacing{true};


      std::function<uint32_t()> m_pfGetElapsedIterations;

      double m_RMSChange{0.0};

      std::function<double()> m_pfGetMetric;

      // Holder of process object for active measurements
      itk::ProcessObject *m_Filter{nullptr};

    };


  }
}
#endif
