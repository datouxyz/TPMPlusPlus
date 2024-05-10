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
#ifndef sitkPhysicalPointImageSource_h
#define sitkPhysicalPointImageSource_h

/*
 * WARNING: DO NOT EDIT THIS FILE!
 * THIS FILE IS AUTOMATICALLY GENERATED BY THE SIMPLEITK BUILD PROCESS.
 * Please look at sitkImageSourceTemplate.h.in to make changes.
 */

#include <memory>

#include "sitkBasicFilters.h"
#include "sitkImageFilter.h"

namespace itk {
  namespace simple {

    /**\class PhysicalPointImageSource
\brief Generate an image of the physical locations of each pixel.

This image source supports image which have a multi-component pixel equal to the image dimension, and variable length VectorImages. It is recommended that the component type be a real valued type.
\sa itk::simple::PhysicalPointSource for the procedural interface
\sa itk::PhysicalPointImageSource for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT PhysicalPointImageSource : public ImageFilter {
    public:
      using Self = PhysicalPointImageSource;

      /** Destructor */
      virtual ~PhysicalPointImageSource();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      PhysicalPointImageSource();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = VectorPixelIDTypeList;



\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetOutputPixelType ( PixelIDValueEnum OutputPixelType ) { this->m_OutputPixelType = OutputPixelType; return *this; }

      /**
       */
      PixelIDValueEnum GetOutputPixelType() const { return this->m_OutputPixelType; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetSize ( std::vector<unsigned int> Size ) { this->m_Size = std::move(Size); return *this; }

      /**
       */
      std::vector<unsigned int> GetSize() const { return this->m_Size; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetOrigin ( std::vector<double> Origin ) { this->m_Origin = std::move(Origin); return *this; }

      /**
       */
      std::vector<double> GetOrigin() const { return this->m_Origin; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetSpacing ( std::vector<double> Spacing ) { this->m_Spacing = std::move(Spacing); return *this; }

      /**
       */
      std::vector<double> GetSpacing() const { return this->m_Spacing; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetDirection ( std::vector<double> Direction ) { this->m_Direction = Direction; return *this; }

      /**
       */
      std::vector<double> GetDirection() const { return this->m_Direction; }

      /** Name of this class */
      std::string GetName() const { return std::string ("PhysicalPointImageSource"); }

      /** Print ourselves out */
      std::string ToString() const;


      /** Execute the filter on the input image */

      Image Execute (  );


      /** This methods sets the size, origin, spacing and direction to that of the provided image */
      void SetReferenceImage(const Image & refImage );


    private:

      /** Setup for member function dispatching */

      using MemberFunctionType = Image (Self::*)(  );
      template <class TImageType> Image ExecuteInternal (  );


      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;

      std::unique_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;


      PixelIDValueEnum  m_OutputPixelType{itk::simple::sitkVectorFloat32};

      std::vector<unsigned int>  m_Size{std::vector<unsigned int>(3, 64)};

      std::vector<double>  m_Origin{std::vector<double>(3, 0.0)};

      std::vector<double>  m_Spacing{std::vector<double>(3, 1.0)};

      /* Passing a zero sized array, defaults to identiy matrix. The size of the array must exactly match the direction matrix for the dimension of the image. */
      std::vector<double>  m_Direction{std::vector<double>()};





    };



   /**
     * \brief Generate an image of the physical locations of each pixel.
     *
     * This function directly calls the execute method of PhysicalPointImageSource
     * in order to support a procedural API
     *
     * \sa itk::simple::PhysicalPointImageSource for the object oriented interface
     */
SITKBasicFilters_EXPORT Image PhysicalPointSource ( PixelIDValueEnum outputPixelType = itk::simple::sitkVectorFloat32, std::vector<unsigned int> size = std::vector<unsigned int>(3, 64), std::vector<double> origin = std::vector<double>(3, 0.0), std::vector<double> spacing = std::vector<double>(3, 1.0), std::vector<double> direction = std::vector<double>() );
  }
}
#endif
