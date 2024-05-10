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
#ifndef sitkMinMaxCurvatureFlowImageFilter_h
#define sitkMinMaxCurvatureFlowImageFilter_h

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

    /**\class MinMaxCurvatureFlowImageFilter
\brief Denoise an image using min/max curvature flow.

MinMaxCurvatureFlowImageFilter implements a curvature driven image denoising algorithm. Iso-brightness contours in the grayscale input image are viewed as a level set. The level set is then evolved using a curvature-based speed function:

 \f[ I_t = F_{\mbox{minmax}} |\nabla I| \f] 

where \f$ F_{\mbox{minmax}} = \max(\kappa,0) \f$ if \f$ \mbox{Avg}_{\mbox{stencil}}(x) \f$ is less than or equal to \f$ T_{threshold} \f$ and \f$ \min(\kappa,0) \f$ , otherwise. \f$ \kappa \f$ is the mean curvature of the iso-brightness contour at point \f$ x \f$ .

In min/max curvature flow, movement is turned on or off depending on the scale of the noise one wants to remove. Switching depends on the average image value of a region of radius \f$ R \f$ around each point. The choice of \f$ R \f$ , the stencil radius, governs the scale of the noise to be removed.

The threshold value \f$ T_{threshold} \f$ is the average intensity obtained in the direction perpendicular to the gradient at point \f$ x \f$ at the extrema of the local neighborhood.

This filter make use of the multi-threaded finite difference solver hierarchy. Updates are computed using a MinMaxCurvatureFlowFunction object. A zero flux Neumann boundary condition is used when computing derivatives near the data boundary.

\warning This filter assumes that the input and output types have the same dimensions. This filter also requires that the output image pixels are of a real type. This filter works for any dimensional images, however for dimensions greater than 3D, an expensive brute-force search is used to compute the local threshold.


Reference: "Level Set Methods and Fast Marching Methods", J.A. Sethian, Cambridge Press, Chapter 16, Second edition, 1999.

\see MinMaxCurvatureFlowFunction 


\see CurvatureFlowImageFilter 


\see BinaryMinMaxCurvatureFlowImageFilter
\sa itk::simple::MinMaxCurvatureFlow for the procedural interface
\sa itk::MinMaxCurvatureFlowImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT MinMaxCurvatureFlowImageFilter : public ImageFilter {
    public:
      using Self = MinMaxCurvatureFlowImageFilter;

      /** Destructor */
      virtual ~MinMaxCurvatureFlowImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      MinMaxCurvatureFlowImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = RealPixelIDTypeList;
\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetTimeStep ( double TimeStep ) { this->m_TimeStep = TimeStep; return *this; }

      /**
       */
      double GetTimeStep() const { return this->m_TimeStep; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetNumberOfIterations ( uint32_t NumberOfIterations ) { this->m_NumberOfIterations = NumberOfIterations; return *this; }

      /**
       */
      uint32_t GetNumberOfIterations() const { return this->m_NumberOfIterations; }\

      /**
       * Set/Get the stencil radius.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetStencilRadius ( int StencilRadius ) { this->m_StencilRadius = StencilRadius; return *this; }

      /**
       * Set/Get the stencil radius.
       */
      int GetStencilRadius() const { return this->m_StencilRadius; }

      /** Name of this class */
      std::string GetName() const { return std::string ("MinMaxCurvatureFlowImageFilter"); }

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


      /* Time step for PDE solver */
      double  m_TimeStep{0.05};

      /* Number of iterations to run */
      uint32_t  m_NumberOfIterations{5u};

      int  m_StencilRadius{2};


      bool m_InPlace{false};
    };

    /**\
     * \brief Denoise an image using min/max curvature flow.
     *
     * This function directly calls the execute method of MinMaxCurvatureFlowImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::MinMaxCurvatureFlowImageFilter for the object oriented interface
     * @{
     */
#ifndef SWIG
     SITKBasicFilters_EXPORT Image MinMaxCurvatureFlow ( Image&& image1, double timeStep = 0.05, uint32_t numberOfIterations = 5u, int stencilRadius = 2 );
#endif
     SITKBasicFilters_EXPORT Image MinMaxCurvatureFlow ( const Image& image1, double timeStep = 0.05, uint32_t numberOfIterations = 5u, int stencilRadius = 2 );

     /** @} */
  }
}
#endif
