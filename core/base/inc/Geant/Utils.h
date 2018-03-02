//===------------------ Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantError.h
 * @brief Error handling routines.
 * @author Philippe Canal
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_UTILS_H
#define GEANT_UTILS_H

#include "Geant/Config.h"
#include <string>
#include <libgen.h>

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

inline std::string GetDataFileLocation(int argc, char **argv, const char *dataFilename)
{
  const char *dir = argc > 0 ? dirname(argv[0]) : ".";

  if (dir == nullptr) dir = ".";
  return std::string(dir) + "/" + dataFilename;
}

inline std::string GetDataFileLocation(int argc, char **argv, std::string dataFilename)
{
  return GetDataFileLocation(argc, argv, dataFilename.c_str());
}

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif // GEANT_UTILS_H
