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
#ifndef sitkCollidingFrontsImageFilter_h
#define sitkCollidingFrontsImageFilter_h

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

    /**\class CollidingFrontsImageFilter
\brief Selects a region of space where two independent fronts run towards each other.

The filter can be used to quickly segment anatomical structures (e.g. for level set initialization).

The filter uses two instances of FastMarchingUpwindGradientImageFilter to compute the gradients of arrival times of two wavefronts propagating from two sets of seeds. The input of the filter is used as the speed of the two wavefronts. The output is the dot product between the two gradient vector fields.

The filter works on the following basic idea. In the regions where the dot product between the two gradient fields is negative, the two fronts propagate in opposite directions. In the regions where the dot product is positive, the two fronts propagate in the same direction. This can be used to extract the region of space between two sets of points.

If StopOnTargets is On, then each front will stop as soon as all seeds of the other front have been reached. This can markedly speed up the execution of the filter, since wave propagation does not take place on the complete image.

Optionally, a connectivity criterion can be applied to the resulting dot product image. In this case, the only negative region in the output image is the one connected to the seeds.

\author Luca Antiga Ph.D. Biomedical Technologies Laboratory, Bioengineering Department, Mario Negri Institute, Italy.
\sa itk::simple::CollidingFronts for the procedural interface
\sa itk::CollidingFrontsImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT CollidingFrontsImageFilter : public ImageFilter {
    public:
      using Self = CollidingFrontsImageFilter;

      /** Destructor */
      virtual ~CollidingFrontsImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      CollidingFrontsImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;
\

      /**
       * Set the container of Seed Points representing the first initial front. Seed points are represented as a VectorContainer of LevelSetNodes.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetSeedPoints1 ( std::vector< std::vector<unsigned int> > SeedPoints1 ) { this->m_SeedPoints1 = std::move(SeedPoints1); return *this; }

      /**
       * Get the container of Seed Points representing the first initial front.
       */
      std::vector< std::vector< unsigned int > > GetSeedPoints1() const { return this->m_SeedPoints1; }
      /** \brief Add SeedPoints1 point */
      SITK_RETURN_SELF_TYPE_HEADER AddSeedPoint1( std::vector< unsigned int > point ) { this->m_SeedPoints1.push_back(std::move(point)); return *this;}
      /** \brief Remove all SeedPoints1 points */
      SITK_RETURN_SELF_TYPE_HEADER ClearSeedPoints1( ) { this->m_SeedPoints1.clear(); return *this;}
\

      /**
       * Set the container of Seed Points representing the second initial front. Seed points are represented as a VectorContainer of LevelSetNodes.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetSeedPoints2 ( std::vector< std::vector<unsigned int> > SeedPoints2 ) { this->m_SeedPoints2 = std::move(SeedPoints2); return *this; }

      /**
       * Get the container of Seed Points representing the second initial front.
       */
      std::vector< std::vector< unsigned int > > GetSeedPoints2() const { return this->m_SeedPoints2; }
      /** \brief Add SeedPoints2 point */
      SITK_RETURN_SELF_TYPE_HEADER AddSeedPoint2( std::vector< unsigned int > point ) { this->m_SeedPoints2.push_back(std::move(point)); return *this;}
      /** \brief Remove all SeedPoints2 points */
      SITK_RETURN_SELF_TYPE_HEADER ClearSeedPoints2( ) { this->m_SeedPoints2.clear(); return *this;}
\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetApplyConnectivity ( bool ApplyConnectivity ) { this->m_ApplyConnectivity = ApplyConnectivity; return *this; }

      /** Set the value of ApplyConnectivity to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER ApplyConnectivityOn() { return this->SetApplyConnectivity(true); }
      SITK_RETURN_SELF_TYPE_HEADER ApplyConnectivityOff() { return this->SetApplyConnectivity(false); }

      /**
       */
      bool GetApplyConnectivity() const { return this->m_ApplyConnectivity; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetNegativeEpsilon ( double NegativeEpsilon ) { this->m_NegativeEpsilon = NegativeEpsilon; return *this; }

      /**
       */
      double GetNegativeEpsilon() const { return this->m_NegativeEpsilon; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetStopOnTargets ( bool StopOnTargets ) { this->m_StopOnTargets = StopOnTargets; return *this; }

      /** Set the value of StopOnTargets to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER StopOnTargetsOn() { return this->SetStopOnTargets(true); }
      SITK_RETURN_SELF_TYPE_HEADER StopOnTargetsOff() { return this->SetStopOnTargets(false); }

      /**
       */
      bool GetStopOnTargets() const { return this->m_StopOnTargets; }

      /** Name of this class */
      std::string GetName() const { return std::string ("CollidingFrontsImageFilter"); }

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


      std::vector< std::vector<unsigned int> >  m_SeedPoints1{std::vector< std::vector<unsigned int > >()};

      std::vector< std::vector<unsigned int> >  m_SeedPoints2{std::vector< std::vector<unsigned int > >()};

      bool  m_ApplyConnectivity{true};

      double  m_NegativeEpsilon{-1e-6};

      bool  m_StopOnTargets{false};


    };

    /**\
     * \brief Selects a region of space where two independent fronts run towards each other.
     *
     * This function directly calls the execute method of CollidingFrontsImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::CollidingFrontsImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image CollidingFronts ( const Image& image1, std::vector< std::vector<unsigned int> > seedPoints1 = std::vector< std::vector<unsigned int > >(), std::vector< std::vector<unsigned int> > seedPoints2 = std::vector< std::vector<unsigned int > >(), bool applyConnectivity = true, double negativeEpsilon = -1e-6, bool stopOnTargets = false );

     /** @} */
  }
}
#endif
