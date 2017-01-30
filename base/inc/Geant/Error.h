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

#ifndef GEANT_ERROR_H
#define GEANT_ERROR_H

#include "Geant/Config.h"

#include <cstdio>
#include <mutex>

namespace Geant {
  enum class EMsgLevel {
    kUnset    =  -1,
    kPrint    =   0,
    kInfo     =   1000,
    kWarning  =   2000,
    kError    =   3000,
    kBreak    =   4000,
    kSysError =   5000,
    kFatal    =   6000
  };
   
  inline namespace GEANT_IMPL_NAMESPACE {

#ifndef VECCORE_CUDA
    void ErrorHandlerImpl(EMsgLevel level, const char *location, const char *msgfmt, ...);
#endif

    template <typename... ArgsTypes>
    VECCORE_ATT_HOST_DEVICE
    void MessageHandler(EMsgLevel level, const char *location, const char *msgfmt, ArgsTypes... params) {
#ifdef VECCORE_CUDA
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
      static std::mutex prntMutex;
#endif
      const char *type = nullptr;
      switch(level) {
        case EMsgLevel::kPrint: type = "Print"; break;
        case EMsgLevel::kInfo: type = "Info"; break;
        case EMsgLevel::kWarning: type = "Warning"; break;
        case EMsgLevel::kError: type = "Error"; break;
        case EMsgLevel::kBreak: type = "Break"; break;
        case EMsgLevel::kSysError: type = "SysError"; break;
        case EMsgLevel::kFatal: type = "Fatal"; break;
        default: type = "Unknown Level"; break;
      }
      { // print mutex scope
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
        std::lock_guard<std::mutex> lock(prntMutex);
#endif
        if (level == EMsgLevel::kPrint)
          printf("%s:",location);
        else
          printf("%s in <%s>:",type, location);
        printf(msgfmt,params...);
        printf("\n");
      }
      if (level >= EMsgLevel::kFatal) {
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
        // Did not find a way to halt a kernel from within yet.
        //cudaDeviceReset();
        //cudaThreadExit();
        //throw("Fatal error in CUDA kernel");
#else
        exit( EXIT_FAILURE );
#endif
      }
#else
      // Currently we use the ROOT message handler on the host/gcc code.
      ErrorHandlerImpl(level,location,msgfmt,params...);
#endif
    }

    template <typename... ArgsTypes>
    VECCORE_ATT_HOST_DEVICE
    void Printf(const char *msgfmt, ArgsTypes... params)
    {
      MessageHandler(EMsgLevel::kPrint,"",msgfmt, params...);
    }

    template <typename... ArgsTypes>
    VECCORE_ATT_HOST_DEVICE
    void Print(const char *location, const char *msgfmt, ArgsTypes... params)
    {
      MessageHandler(EMsgLevel::kPrint,location,msgfmt, params...);
    }

    template <typename... ArgsTypes>
    VECCORE_ATT_HOST_DEVICE
    void Info(const char *location, const char *msgfmt, ArgsTypes... params)
    {
      MessageHandler(EMsgLevel::kInfo,location,msgfmt, params...);
    }

    template <typename... ArgsTypes>
    VECCORE_ATT_HOST_DEVICE
    void Warning(const char *location, const char *msgfmt, ArgsTypes... params)
    {
      MessageHandler(EMsgLevel::kWarning,location,msgfmt, params...);
    }

    template <typename... ArgsTypes>
    VECCORE_ATT_HOST_DEVICE
    void Error(const char *location, const char *msgfmt, ArgsTypes... params)
    {
      MessageHandler(EMsgLevel::kError,location,msgfmt, params...);
    }

    template <typename... ArgsTypes>
    VECCORE_ATT_HOST_DEVICE
    void Fatal(const char *location, const char *msgfmt, ArgsTypes... params)
    {
      MessageHandler(EMsgLevel::kFatal,location,msgfmt, params...);
    }

  } // GEANT_IMPL_NAMESPACE
} // Geant

#endif // GEANT_ERROR_H
