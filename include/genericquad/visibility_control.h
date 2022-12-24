#ifndef GENERICQUAD__VISIBILITY_CONTROL_H_
#define GENERICQUAD__VISIBILITY_CONTROL_H_
#if defined _WIN32 || defined __CYGWIN__
  #ifdef __GNUC__
    #define GENERICQUAD_EXPORT __attribute__ ((dllexport))
    #define GENERICQUAD_IMPORT __attribute__ ((dllimport))
  #else
    #define GENERICQUAD_EXPORT __declspec(dllexport)
    #define GENERICQUAD_IMPORT __declspec(dllimport)
  #endif
  #ifdef GENERICQUAD_BUILDING_LIBRARY
    #define GENERICQUAD_PUBLIC GENERICQUAD_EXPORT
  #else
    #define GENERICQUAD_PUBLIC GENERICQUAD_IMPORT
  #endif
  #define GENERICQUAD_PUBLIC_TYPE GENERICQUAD_PUBLIC
  #define GENERICQUAD_LOCAL
#else
  #define GENERICQUAD_EXPORT __attribute__ ((visibility("default")))
  #define GENERICQUAD_IMPORT
  #if __GNUC__ >= 4
    #define GENERICQUAD_PUBLIC __attribute__ ((visibility("default")))
    #define GENERICQUAD_LOCAL  __attribute__ ((visibility("hidden")))
  #else
    #define GENERICQUAD_PUBLIC
    #define GENERICQUAD_LOCAL
  #endif
  #define GENERICQUAD_PUBLIC_TYPE
#endif
#endif  // GENERICQUAD__VISIBILITY_CONTROL_H_
// Generated 24-Dec-2022 10:32:42
// Copyright 2019-2020 The MathWorks, Inc.
