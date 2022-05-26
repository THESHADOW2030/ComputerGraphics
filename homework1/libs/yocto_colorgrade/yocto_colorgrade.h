//
// Yocto/Grade: Tiny library for color grading.
//

//
// LICENSE:
//
// Copyright (c) 2020 -- 2020 Fabio Pellacini
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//

#ifndef _YOCTO_COLORGRADE_
#define _YOCTO_COLORGRADE_

#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>

// -----------------------------------------------------------------------------
// COLOR GRADING FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Color grading parameters
struct grade_params {
  float exposure   = 0.0f;
  bool  filmic     = false;
  bool  srgb       = true;
  vec3f tint       = vec3f{1, 1, 1};
  float saturation = 0.5f;
  float contrast   = 0.5f;
  float vignette   = 0.0f;
  float grain      = 0.0f;
  int   mosaic     = 0;
  int   grid       = 0;
  int   r          = 3;
 // vec3f   segregation = vec3f{0, 0, 1};
  //vec3f sepia       = {0, 0, 0};
  float StrongdgeDetectionMatrix[3][3]     = {{1, -2, 1}, {-2, 4, -2}, {1, -2, 1}};
  float MediumedgeDetectionMatrix[3][3] = {{-1, -1, -1}, {-1, 8, -1}, {-1, -1, -1}};
  float LowdgeDetectionMatrix[3][3] = {{0, 1, 0}, {1, -4, 1}, {0, 1, 0}};
  float embossMatrix[3][3] = {{-0.5, 0, 0}, {0, 0, 0}, {0, 0, 0.5}};
  float enhancedDetailsMatrix[3][3] = {{0, -1, 0}, {-1, 10, -1}, {0, -1, 0}};
  float HeavySoftenMatrix[3][3] = {{11, 11, 11}, {11, 11, 11}, {11, 11, 11}};
  float MediumSoftenMatrix[3][3] = {{10, 10, 10}, {10, 25, 10}, {10, 10, 10}};
  float LightSoftenMatrix[3][3]  = {{6, 12, 6}, {12, 25, 12}, {6, 12, 6}};

  //float blurLights[3][3] = {{1, 2, 1}, {2, 2, 2}, {1, 2, 1}};
  //bool  channelSegregation = false;
  bool  inverti            = false;
  bool  blur               = false;
  int   blurIntensity      = 10;
  bool  cut                = false;
  bool  edgeDetection      = false;
  int   livelloEdgeDetection            = 1;
  bool  emboss             = false;
  bool  extraCredits                    = false;
  bool  soften                          = false;
  bool  enhanceDetails                  = false;
  bool  oilPainting                     = false;
  int   intensityOil                    = 20;
  int   rang                            = 5;
  const int   kernelW            = 3;
  

  bool  gauss              = false;
  float sigma = 6;
  int   gaussBlurIntensity = 20;
  int   anaglyphOffSet     = 30;
  bool  trueAnaglyph       = false;
  bool  grayAnaglyph       = false;
  bool  colorAnaglyph       = false;
  bool  halfColorAnaglyph       = false;
  bool  optimizedColorAnaglyph = false;
  int  livelloSoften           = 1;
  
  bool sketch = false;
  int  boh    = 9;
  bool CRTTV       = false;
};

// Grading functions
color_image grade_image(const color_image& image, const grade_params& params);

};  // namespace yocto

#endif
