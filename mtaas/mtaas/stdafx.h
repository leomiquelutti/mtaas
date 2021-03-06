// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
// Windows Header Files:
#include <windows.h>
//#include <memory.h>
//#include <stdio.h>
//#include <tchar.h>
#include <string>
#include <complex>
#include <vector>

#define N_CHANNEL_INPUTS 2
typedef std::complex<double> Complex;

#include "matCUDA.h"
#include "MT.h"
//#include "jni.h"