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
#ifndef sitkConfigure_h
#define sitkConfigure_h


/* #undef SITK_BUILD_SHARED_LIBS */
#ifdef SITK_BUILD_SHARED_LIBS
#define SITKDLL
#else
#define SITKSTATIC
#endif
/* #undef SITK_SimpleITKExplit_STATIC */

#define SITK_MAX_DIMENSION 5

// defined if compiler supports using template keyword to disambiguate
// dependent names
#define SITK_HAS_TEMPLATE_DISAMBIGUATOR_DEPENDENT_NAME

#define SITK_INT64_PIXELIDS

/* #undef SITK_EXPLICIT_INSTANTIATION */

// Include ITK version reported in CMake with SITK prefix, so that
// SimpleITK doesn't need ITK header in our headers.
#define SITK_ITK_VERSION_MAJOR 5
#define SITK_ITK_VERSION_MINOR 2
#define SITK_ITK_VERSION_PATCH 0

#endif // sitkConfigure_h
