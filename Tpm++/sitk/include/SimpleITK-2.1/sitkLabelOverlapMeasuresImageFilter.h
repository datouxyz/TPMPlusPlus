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
#ifndef sitkLabelOverlapMeasuresImageFilter_h
#define sitkLabelOverlapMeasuresImageFilter_h

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

    /**\class LabelOverlapMeasuresImageFilter
\brief Computes overlap measures between the set same set of labels of pixels of two images. Background is assumed to be 0.

This code was contributed in the Insight Journal paper: "Introducing Dice, Jaccard, and Other Label Overlap Measures To ITK" by Nicholas J. Tustison, James C. Gee https://hdl.handle.net/10380/3141 http://www.insight-journal.org/browse/publication/707 

\author Nicholas J. Tustison 


\see LabelOverlapMeasuresImageFilter

\sa itk::LabelOverlapMeasuresImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT LabelOverlapMeasuresImageFilter : public ImageFilter {
    public:
      using Self = LabelOverlapMeasuresImageFilter;

      /** Destructor */
      virtual ~LabelOverlapMeasuresImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      LabelOverlapMeasuresImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = IntegerPixelIDTypeList;

     /**
      * Get the false negative error for the specified individual label.
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     double GetFalseNegativeError() const { return this->m_FalseNegativeError; };

     /**
      * Get the false negative error over all labels.
      *
      * This is an active measurement. It may be accessed while the
      * filter is being executing in command call-backs and can be
      * accessed after execution.
      */
     double GetFalseNegativeError(int64_t label) const { return this->m_pfGetFalseNegativeError(label); };

     /**
      * Get the false positive error for the specified individual label.
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     double GetFalsePositiveError() const { return this->m_FalsePositiveError; };

     /**
      * Get the false positive error over all labels.
      *
      * This is an active measurement. It may be accessed while the
      * filter is being executing in command call-backs and can be
      * accessed after execution.
      */
     double GetFalsePositiveError(int64_t label) const { return this->m_pfGetFalsePositiveError(label); };

     /**
      * Get the mean overlap (Dice coefficient) for the specified individual label.
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     double GetMeanOverlap() const { return this->m_MeanOverlap; };

     /**
      * Get the mean overlap (Dice coefficient) over all labels.
      *
      * This is an active measurement. It may be accessed while the
      * filter is being executing in command call-backs and can be
      * accessed after execution.
      */
     double GetMeanOverlap(int64_t label) const { return this->m_pfGetMeanOverlap(label); };

     /**
      * Get the union overlap (Jaccard coefficient) for the specified individual label.
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     double GetUnionOverlap() const { return this->m_UnionOverlap; };

     /**
      * Get the union overlap (Jaccard coefficient) over all labels.
      *
      * This is an active measurement. It may be accessed while the
      * filter is being executing in command call-backs and can be
      * accessed after execution.
      */
     double GetUnionOverlap(int64_t label) const { return this->m_pfGetUnionOverlap(label); };

     /**
      * Get the volume similarity for the specified individual label.
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     double GetVolumeSimilarity() const { return this->m_VolumeSimilarity; };

     /**
      * Get the volume similarity over all labels.
      *
      * This is an active measurement. It may be accessed while the
      * filter is being executing in command call-backs and can be
      * accessed after execution.
      */
     double GetVolumeSimilarity(int64_t label) const { return this->m_pfGetVolumeSimilarity(label); };

     /**
      * Get the union overlap (Jaccard coefficient) for the specified individual label.
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     double GetJaccardCoefficient() const { return this->m_JaccardCoefficient; };

     /**
      * Get the union overlap (Jaccard coefficient) over all labels.
      *
      * This is an active measurement. It may be accessed while the
      * filter is being executing in command call-backs and can be
      * accessed after execution.
      */
     double GetJaccardCoefficient(int64_t label) const { return this->m_pfGetJaccardCoefficient(label); };

     /**
      * Get the mean overlap (Dice coefficient) for the specified individual label.
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     double GetDiceCoefficient() const { return this->m_DiceCoefficient; };

     /**
      * Get the mean overlap (Dice coefficient) over all labels.
      *
      * This is an active measurement. It may be accessed while the
      * filter is being executing in command call-backs and can be
      * accessed after execution.
      */
     double GetDiceCoefficient(int64_t label) const { return this->m_pfGetDiceCoefficient(label); };


      /** Name of this class */
      std::string GetName() const { return std::string ("LabelOverlapMeasuresImageFilter"); }

      /** Print ourselves out */
      std::string ToString() const;


      /** Execute the filter on the input images */

      void Execute ( const Image& image1, const Image& image2 );

    private:

      /** Setup for member function dispatching */

      using MemberFunctionType = void (Self::*)( const Image& image1, const Image& image2 );
      template <class TImageType> void ExecuteInternal ( const Image& image1, const Image& image2 );


      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;

      std::unique_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;



      double m_FalseNegativeError{0.0};

      std::function<double(int64_t)> m_pfGetFalseNegativeError;

      double m_FalsePositiveError{0.0};

      std::function<double(int64_t)> m_pfGetFalsePositiveError;

      double m_MeanOverlap{0.0};

      std::function<double(int64_t)> m_pfGetMeanOverlap;

      double m_UnionOverlap{0.0};

      std::function<double(int64_t)> m_pfGetUnionOverlap;

      double m_VolumeSimilarity{0.0};

      std::function<double(int64_t)> m_pfGetVolumeSimilarity;

      double m_JaccardCoefficient{0.0};

      std::function<double(int64_t)> m_pfGetJaccardCoefficient;

      double m_DiceCoefficient{0.0};

      std::function<double(int64_t)> m_pfGetDiceCoefficient;

      // Holder of process object for active measurements
      itk::ProcessObject *m_Filter{nullptr};

    };


  }
}
#endif