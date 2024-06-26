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
#ifndef sitkFastMarchingUpwindGradientImageFilter_h
#define sitkFastMarchingUpwindGradientImageFilter_h

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

    /**\class FastMarchingUpwindGradientImageFilter
\brief Generates the upwind gradient field of fast marching arrival times.

This filter adds some extra functionality to its base class. While the solution T(x) of the Eikonal equation is being generated by the base class with the fast marching method, the filter generates the upwind gradient vectors of T(x), storing them in an image.

Since the Eikonal equation generates the arrival times of a wave traveling at a given speed, the generated gradient vectors can be interpreted as the slowness (1/velocity) vectors of the front (the quantity inside the modulus operator in the Eikonal equation).

Gradient vectors are computed using upwind finite differences, that is, information only propagates from points where the wavefront has already passed. This is consistent with how the fast marching method works.

One more extra feature is the possibility to define a set of Target points where the propagation stops. This can be used to avoid computing the Eikonal solution for the whole domain. The front can be stopped either when one Target point is reached or all Target points are reached. The propagation can stop after a time TargetOffset has passed since the stop condition is met. This way the solution is computed a bit downstream the Target points, so that the level sets of T(x) corresponding to the Target are smooth.

For an alternative implementation, see itk::FastMarchingUpwindGradientImageFilterBase .

\author Luca Antiga Ph.D. Biomedical Technologies Laboratory, Bioengineering Department, Mario Negri Institute, Italy.
\sa itk::simple::FastMarchingUpwindGradient for the procedural interface
\sa itk::FastMarchingUpwindGradientImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT FastMarchingUpwindGradientImageFilter : public ImageFilter {
    public:
      using Self = FastMarchingUpwindGradientImageFilter;

      /** Destructor */
      virtual ~FastMarchingUpwindGradientImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      FastMarchingUpwindGradientImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;
\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetTrialPoints ( std::vector< std::vector<unsigned int> > TrialPoints ) { this->m_TrialPoints = std::move(TrialPoints); return *this; }

      /**
       */
      std::vector< std::vector< unsigned int > > GetTrialPoints() const { return this->m_TrialPoints; }
      /** \brief Add TrialPoints point */
      SITK_RETURN_SELF_TYPE_HEADER AddTrialPoint( std::vector< unsigned int > point ) { this->m_TrialPoints.push_back(std::move(point)); return *this;}
      /** \brief Remove all TrialPoints points */
      SITK_RETURN_SELF_TYPE_HEADER ClearTrialPoints( ) { this->m_TrialPoints.clear(); return *this;}
\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetNumberOfTargets ( unsigned int NumberOfTargets ) { this->m_NumberOfTargets = NumberOfTargets; return *this; }

      /**
       * Get the number of targets.
       */
      unsigned int GetNumberOfTargets() const { return this->m_NumberOfTargets; }\

      /**
       * Set the container of Target Points. If a target point is reached, the propagation stops. Trial points are represented as a VectorContainer of LevelSetNodes.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetTargetPoints ( std::vector< std::vector<unsigned int> > TargetPoints ) { this->m_TargetPoints = std::move(TargetPoints); return *this; }

      /**
       * Get the container of Target Points.
       */
      std::vector< std::vector< unsigned int > > GetTargetPoints() const { return this->m_TargetPoints; }
      /** \brief Add TargetPoints point */
      SITK_RETURN_SELF_TYPE_HEADER AddTargetPoint( std::vector< unsigned int > point ) { this->m_TargetPoints.push_back(std::move(point)); return *this;}
      /** \brief Remove all TargetPoints points */
      SITK_RETURN_SELF_TYPE_HEADER ClearTargetPoints( ) { this->m_TargetPoints.clear(); return *this;}
\

      /**
       * Set how long (in terms of arrival times) after targets are reached the front must stop. This is useful to ensure that the level set of target arrival time is smooth.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetTargetOffset ( double TargetOffset ) { this->m_TargetOffset = TargetOffset; return *this; }

      /**
       * Get the TargetOffset ivar.
       */
      double GetTargetOffset() const { return this->m_TargetOffset; }\

      /**
       * Set/Get the Normalization Factor for the Speed Image . The values in the Speed Image is divided by this factor. This allows the use of images with integer pixel types to represent the speed.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetNormalizationFactor ( double NormalizationFactor ) { this->m_NormalizationFactor = NormalizationFactor; return *this; }

      /**
       * Set/Get the Normalization Factor for the Speed Image . The values in the Speed Image is divided by this factor. This allows the use of images with integer pixel types to represent the speed.
       */
      double GetNormalizationFactor() const { return this->m_NormalizationFactor; }
     /**
      * Get the gradient image.
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     Image GetGradientImage() const { return this->m_GradientImage; };

     /**
      * Get the arrival time corresponding to the last reached target. If TargetReachedMode is set to NoTargets, TargetValue contains the last (aka largest) Eikonal solution value generated.
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     double GetTargetValue() const { return this->m_TargetValue; };


      /** Name of this class */
      std::string GetName() const { return std::string ("FastMarchingUpwindGradientImageFilter"); }

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


      std::vector< std::vector<unsigned int> >  m_TrialPoints{std::vector< std::vector<unsigned int > >()};

      unsigned int  m_NumberOfTargets{0u};

      std::vector< std::vector<unsigned int> >  m_TargetPoints{std::vector< std::vector<unsigned int > >()};

      double  m_TargetOffset{1.0};

      double  m_NormalizationFactor{1.0};

      /* Docs */
      Image m_GradientImage{Image()};
      /* Docs */
      double m_TargetValue{0.0};


    };

    /**\
     * \brief Generates the upwind gradient field of fast marching arrival times.
     *
     * This function directly calls the execute method of FastMarchingUpwindGradientImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::FastMarchingUpwindGradientImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image FastMarchingUpwindGradient ( const Image& image1, std::vector< std::vector<unsigned int> > trialPoints = std::vector< std::vector<unsigned int > >(), unsigned int numberOfTargets = 0u, std::vector< std::vector<unsigned int> > targetPoints = std::vector< std::vector<unsigned int > >(), double targetOffset = 1.0, double normalizationFactor = 1.0 );

     /** @} */
  }
}
#endif
