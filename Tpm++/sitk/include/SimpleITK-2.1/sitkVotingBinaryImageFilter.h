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
#ifndef sitkVotingBinaryImageFilter_h
#define sitkVotingBinaryImageFilter_h

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

    /**\class VotingBinaryImageFilter
\brief Applies a voting operation in a neighborhood of each pixel.

\note Pixels which are not Foreground or Background will remain unchanged.


\see Image 


\see Neighborhood 


\see NeighborhoodOperator 


\see NeighborhoodIterator
\sa itk::simple::VotingBinary for the procedural interface
\sa itk::VotingBinaryImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT VotingBinaryImageFilter : public ImageFilter {
    public:
      using Self = VotingBinaryImageFilter;

      /** Destructor */
      virtual ~VotingBinaryImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      VotingBinaryImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = IntegerPixelIDTypeList;
\

      /**
       * Set the radius of the neighborhood used to compute the median.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetRadius ( std::vector<unsigned int> Radius ) { this->m_Radius = std::move(Radius); return *this; }

      /** Set the values of the Radius vector all to value */
      SITK_RETURN_SELF_TYPE_HEADER SetRadius( unsigned int value ) { this->m_Radius = std::vector<unsigned int>(3, value); return *this; }

      /**
       * Get the radius of the neighborhood used to compute the median
       */
      std::vector<unsigned int> GetRadius() const { return this->m_Radius; }\

      /**
       * Birth threshold. Pixels that are OFF will turn ON when the number of neighbors ON is larger than the value defined in this threshold.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetBirthThreshold ( unsigned int BirthThreshold ) { this->m_BirthThreshold = BirthThreshold; return *this; }

      /**
       * Birth threshold. Pixels that are OFF will turn ON when the number of neighbors ON is larger than the value defined in this threshold.
       */
      unsigned int GetBirthThreshold() const { return this->m_BirthThreshold; }\

      /**
       * Survival threshold. Pixels that are ON will turn OFF when the number of neighbors ON is smaller than the value defined in this survival threshold.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetSurvivalThreshold ( unsigned int SurvivalThreshold ) { this->m_SurvivalThreshold = SurvivalThreshold; return *this; }

      /**
       * Survival threshold. Pixels that are ON will turn OFF when the number of neighbors ON is smaller than the value defined in this survival threshold.
       */
      unsigned int GetSurvivalThreshold() const { return this->m_SurvivalThreshold; }\

      /**
       * Set the value associated with the Foreground (or the object) on the binary input image and the Background .
       */
      SITK_RETURN_SELF_TYPE_HEADER SetForegroundValue ( double ForegroundValue ) { this->m_ForegroundValue = ForegroundValue; return *this; }

      /**
       * Get the value associated with the Foreground (or the object) on the binary input image and the Background .
       */
      double GetForegroundValue() const { return this->m_ForegroundValue; }\

      /**
       * Set the value associated with the Foreground (or the object) on the binary input image and the Background .
       */
      SITK_RETURN_SELF_TYPE_HEADER SetBackgroundValue ( double BackgroundValue ) { this->m_BackgroundValue = BackgroundValue; return *this; }

      /**
       * Get the value associated with the Foreground (or the object) on the binary input image and the Background .
       */
      double GetBackgroundValue() const { return this->m_BackgroundValue; }

      /** Name of this class */
      std::string GetName() const { return std::string ("VotingBinaryImageFilter"); }

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
      std::vector<unsigned int>  m_Radius{std::vector<unsigned int>(3, 1)};

      unsigned int  m_BirthThreshold{1u};

      unsigned int  m_SurvivalThreshold{1u};

      double  m_ForegroundValue{1.0};

      double  m_BackgroundValue{0.0};


    };

    /**\
     * \brief Applies a voting operation in a neighborhood of each pixel.
     *
     * This function directly calls the execute method of VotingBinaryImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::VotingBinaryImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image VotingBinary ( const Image& image1, std::vector<unsigned int> radius = std::vector<unsigned int>(3, 1), unsigned int birthThreshold = 1u, unsigned int survivalThreshold = 1u, double foregroundValue = 1.0, double backgroundValue = 0.0 );

     /** @} */
  }
}
#endif
