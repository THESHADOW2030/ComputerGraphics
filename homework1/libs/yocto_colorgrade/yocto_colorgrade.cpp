//
// Implementation for Yocto/Grade.
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

#include "yocto_colorgrade.h"
#include <list>
#include <yocto/yocto_color.h>
#include <yocto/yocto_sampling.h>
#include <map>

// -----------------------------------------------------------------------------
// COLOR GRADING FUNCTIONS
// 
// http://www.3dtv.at/knowhow/anaglyphcomparison_en.aspx
// https://old.cescg.org/CESCG97/boros/
// -----------------------------------------------------------------------------
namespace yocto {

color_image homework24(const color_image& image, const grade_params& params) {
  auto graded = image;
  auto rng    = make_rng(11234);

  for (int j : range(image.height)) {
    for (int i : range(image.width)) {
      auto& c = graded[{i, j}];
      
      // 24punti

      // mosaic Effect
      if (params.mosaic != 0) {
        c = image[{i - i % params.mosaic, j - j % params.mosaic}];
      }
      // grid
      if (!params.grid == 0 && (0 == i % params.grid || 0 == j % params.grid)) {
        c = 0.5 * c;
      }

      // tone mapping
      c = c * pow(2, params.exposure);
      if (params.filmic) {
        c *= 0.6;
        c = (pow(c, 2) * 2.51 + c * 0.03) /
            (pow(c, 2) * 2.43 + c * 0.59 + 0.14);
      }
      if (params.srgb && !params.extraCredits) {
        c = pow(c, 1 / 2.2);
      }
      c = clamp(c, 0, 1);

      // color tint
      auto  boh  = c.w;
      vec3f junk = {c.x, c.y, c.z};
      junk       = junk * params.tint;
      c          = {junk.x, junk.y, junk.z, boh};
      // saturation
      auto g = (c.x + c.y + c.z) / 3;
      c      = g + (c - g) * (params.saturation * 2);
      // contrast
      c = gain(c, 1 - params.contrast);
      // vignette
      auto  vr   = 1 - params.vignette;
      vec2f ji   = {j, i};
      vec2f size = {image.height, image.width};
      auto  r    = length(ji - size / 2) / length(size / 2);
      c          = c * (1 - smoothstep(vr, 2 * vr, r));
      // film grain
      c = c + (rand1f(rng) - 0.5) * params.grain;
      c.w = 1;
    }  // 24punti
    
  }
  return graded;
}


color_image inverti(  // funziona
    const color_image& image, const grade_params& params) {
  auto  graded = image;
  auto& pixels = graded.pixels;
  for (auto&(pixel) : pixels) {
    // printf("\nprima\n%f, %f, %f",pixel.x, pixel.y, pixel.z);
    pixel.x = 1 - pixel.x;
    pixel.y = 1 - pixel.y;
    pixel.z = 1 - pixel.z;
    //  printf("\n dopo\n%f, %f, %f", pixel.x, pixel.y, pixel.z);
  }
  return graded;
}

color_image blur(const color_image& image, const grade_params& params) {

  auto graded = image;

  float mediaR = 0;
  float mediaG = 0;
  float mediaB = 0;
  int   count  = 0;
  vec4f pixel  = image[{0, 0}];
  int   w      = graded.width;
  int   h      = graded.height;
  int   k      = params.blurIntensity;
  int   boh    = 0;
  auto& pixels = graded.pixels;

  for (auto j = 0; j < h; j++) {
    for (auto i = 0; i < w; i++) {
      //  printf("////////////////\ni: %d\n////////////////////", i);
      //  auto c = image[{i, j}];

      mediaR = 0;
      mediaG = 0;
      mediaB = 0;
      count  = 0;

      for (auto y = j - k; y <= j + k && y < image.height; y++) {
        if (y < 0) {
          continue;
        } else {
          for (auto x = i - k; x <= i + k && x < image.width; x++) {
            if (x < 0) {             
              continue;
            } else {
              auto rabbia = image[{x, y}];
              mediaR += rabbia.x;
              mediaG += rabbia.y;
              mediaB += rabbia.z;
              count = count + 1;
            }
          }
        }
      }
      
      pixels[boh].x = mediaR / count;
      pixels[boh].y = mediaG / count;
      pixels[boh].z = mediaB / count;
      boh           = boh + 1;
      
    }
  }
  return graded;

}

float seconda(float ciao) { return ciao * ciao; }





color_image cut(const color_image& image, const grade_params& params) { //funziona

      auto  graded = image;
      auto& pixels = graded.pixels;
      int   X      = graded.width;
      int   Y      = graded.height;
      int   boh    = 0;

      for (auto y = 0; y < Y; y++) {
        for (auto x = 0; x < X; x++) {
          auto& pixel = pixels[boh]; //pixel attuale
          if (((seconda(x-X/2)/seconda(X/2))  + seconda(y-Y/2) /seconda(Y/2))  < 1 ) {

              
          
          } else {

              pixel.x = 0;
              pixel.y = 0;
              pixel.z = 0;
          
          
          }
          boh = boh + 1;
          
        }
      }
      return graded;
     

}


enum Intensity{
    LOW = 0, MEDIUM = 1, STRONG = 2

};

//funziona
color_image edgeDetection(const color_image& image, const grade_params& params, Intensity en ) { //funziona
  auto graded = image;

  float mediaR = 0;
  float mediaG = 0;
  float mediaB = 0;
  int   count  = 0;
  vec4f pixel  = image[{0, 0}];
  int   w      = graded.width;
  int   h      = graded.height;
  int   k      = 10;
  int   boh    = 0;
  auto& pixels = graded.pixels;
  auto  M      = params.StrongdgeDetectionMatrix;
  if (en == Intensity::STRONG) {
     M = params.StrongdgeDetectionMatrix;
  } else if (en == Intensity::MEDIUM) {
     M = params.MediumedgeDetectionMatrix;
  } else if (en == Intensity::LOW) {
    M = params.LowdgeDetectionMatrix;
  }
  

  for (auto j = 0; j < h; j++) {
    for (auto i = 0; i < w; i++) {
      //  printf("////////////////\ni: %d\n////////////////////", i);
      //  auto c = image[{i, j}];

      mediaR = 0;
      mediaG = 0;
      mediaB = 0;
      count  = 0;
      auto dy = 0;
      auto dx = 0;
      for (auto y = j - 1; y <= j + 1 && y < image.height; y++) {
        dx = 0; 
        if (y < 0) {
          continue;
        } else {
          for (auto x = i - 1; x <= i + 1 && x < image.width; x++) {
            
          
            if (x < 0) {
              
              continue;
            } else {
              auto rabbia = image[{x, y}];
              mediaR += rabbia.x * M[dy][dx];
              mediaG += rabbia.y * M[dy][dx];
              mediaB += rabbia.z * M[dy][dx];
              count = count +  1;
              
            }
            dx = dx + 1;
          }
          

        }
        dy = dy + 1;
      }

      pixels[boh].x = mediaR ;
      pixels[boh].y = mediaG ;
      pixels[boh].z = mediaB;
      boh           = boh + 1;
      
    }
  }
  return graded;
  
}

color_image emboss(
    const color_image& image, const grade_params& params) {  // funziona
  auto graded = image;

  float mediaR = 0;
  float mediaG = 0;
  float mediaB = 0;
  int   count  = 0;
  vec4f pixel  = image[{0, 0}];
  int   w      = graded.width;
  int   h      = graded.height;
  int   k      = 10;
  int   boh    = 0;
  auto& pixels = graded.pixels;
  auto  M      = params.embossMatrix;

  for (auto j = 0; j < h; j++) {
    for (auto i = 0; i < w; i++) {
      //  printf("////////////////\ni: %d\n////////////////////", i);
      //  auto c = image[{i, j}];

      mediaR  = 0;
      mediaG  = 0;
      mediaB  = 0;
      count   = 0;
      auto dy = 0;
      auto dx = 0;
      for (auto y = j - 1; y <= j + 1 && y < image.height; y++) {
        dx = 0;
        if (y < 0) {
          continue;
        } else {
          for (auto x = i - 1; x <= i + 1 && x < image.width; x++) {
            // printf("////////////////\nj: %d i: %d\n////////////////////", j,
            // i);
            //  printf("////////////////\ny: %d x: %d\n////////////////////", y,
            //  x);
            if (x < 0) {
              continue;
            } else {
              auto rabbia = image[{x, y}];
              mediaR += rabbia.x * M[dy][dx];
              mediaG += rabbia.y * M[dy][dx];
              mediaB += rabbia.z * M[dy][dx];
              count = count + 1;
            }
            dx = dx + 1;
          }
        }
        dy = dy + 1;
      }

      //  printf("////////////////\ny: %d x: %d\n////////////////////", j, i);

      pixels[boh].x = mediaR + 0.5;
      pixels[boh].y = mediaG + 0.5;
      pixels[boh].z = mediaB + 0.5;
      boh           = boh + 1;
    }
  }
  return graded;
}


color_image EnhancedDetails(
    const color_image& image, const grade_params& params) {  // funziona
  auto graded = image;

  float mediaR = 0;
  float mediaG = 0;
  float mediaB = 0;
  int   count  = 0;
  vec4f pixel  = image[{0, 0}];
  int   w      = graded.width;
  int   h      = graded.height;
  int   k      = 10;
  int   boh    = 0;
  auto& pixels = graded.pixels;
  auto  M      = params.enhancedDetailsMatrix;

  for (auto j = 0; j < h; j++) {
    for (auto i = 0; i < w; i++) {
      //  printf("////////////////\ni: %d\n////////////////////", i);
      //  auto c = image[{i, j}];

      mediaR  = 0;
      mediaG  = 0;
      mediaB  = 0;
      count   = 0;
      auto dy = 0;
      auto dx = 0;
      for (auto y = j - 1; y <= j + 1 && y < image.height; y++) {
        dx = 0;
        if (y < 0) {
          continue;
        } else {
          for (auto x = i - 1; x <= i + 1 && x < image.width; x++) {
            // printf("////////////////\nj: %d i: %d\n////////////////////", j,
            // i);
            //  printf("////////////////\ny: %d x: %d\n////////////////////", y,
            //  x);
            if (x < 0) {
              continue;
            } else {
              auto rabbia = image[{x, y}];
              mediaR += rabbia.x * M[dy][dx];
              mediaG += rabbia.y * M[dy][dx];
              mediaB += rabbia.z * M[dy][dx];
              count = count + 1;
            }
            dx = dx + 1;
          }
        }
        dy = dy + 1;
      }

      //  printf("////////////////\ny: %d x: %d\n////////////////////", j, i);

      pixels[boh].x = mediaR/6;
      pixels[boh].y = mediaG / 6;
      pixels[boh].z = mediaB / 6;
      boh           = boh + 1;
    }
  }
  return graded;
}



color_image softenImage(const color_image& image, const grade_params& params,
    Intensity en) {  // funziona
  auto graded = image;

  float mediaR = 0;
  float mediaG = 0;
  float mediaB = 0;
  int   count  = 0;
  vec4f pixel  = image[{0, 0}];
  int   koeff  = 1;
  int   w      = graded.width;
  int   h      = graded.height;
  int   k      = 10;
  int   boh    = 0;
  auto& pixels = graded.pixels;
  auto  M      = params.StrongdgeDetectionMatrix;
  if (en == Intensity::STRONG) {
    M = params.HeavySoftenMatrix;
    koeff = 99;
  } else if (en == Intensity::MEDIUM) {
    M = params.MediumSoftenMatrix;
    koeff = 100;
  } else if (en == Intensity::LOW) {
    M = params.LightSoftenMatrix;
    koeff = 97;

  }
  
 
  for (auto j = 0; j < h; j++) {
    for (auto i = 0; i < w; i++) {
      //  printf("////////////////\ni: %d\n////////////////////", i);
      //  auto c = image[{i, j}];

      mediaR  = 0;
      mediaG  = 0;
      mediaB  = 0;
      count   = 0;
      auto dy = 0;
      auto dx = 0;
      for (auto y = j - 1; y <= j + 1 && y < image.height; y++) {
        dx = 0;
        if (y < 0) {
          continue;
        } else {
          for (auto x = i - 1; x <= i + 1 && x < image.width; x++) {
            // printf("////////////////\nj: %d i: %d\n////////////////////", j,
            // i);
            //  printf("////////////////\ny: %d x: %d\n////////////////////", y,
            //  x);
            if (x < 0) {
              continue;
            } else {
              auto rabbia = image[{x, y}];
              mediaR += rabbia.x * M[dy][dx];
              mediaG += rabbia.y * M[dy][dx];
              mediaB += rabbia.z * M[dy][dx];
              count = count + 1;
            }
            dx = dx + 1;
          }
        }
        dy = dy + 1;
      }

      //  printf("////////////////\ny: %d x: %d\n////////////////////", j, i);

      pixels[boh].x = mediaR / koeff;
      pixels[boh].y = mediaG / koeff;
      pixels[boh].z = mediaB / koeff;
      boh           = boh + 1;
    }
  }
  return graded;
}



color_image oil( 
    const color_image& image, const grade_params& params) {  // funziona
  auto graded = image;

  
  vec4f pixel  = image[{0, 0}];
  int   w      = graded.width;
  int   h      = graded.height;
  
  int   boh    = 0;
  auto& pixels = graded.pixels;

  //  printf("sto dentro al for di OIL");
  int              min = 0;
  std::list<vec3f> ciao;
  std::map<int, int> dizCount;
  std::map<int, int> dizR;
  std::map<int, int> dizG;
  std::map<int, int> dizB;
  int              max = 0;
  int                intensityLevels = params.intensityOil;
  int                rang            = params.rang;
  for (auto j = 0; j < h; j++) {
    for (auto i = 0; i < w; i++) {
      
      auto& c = graded[{i, j}];
     
      int              min = 0;
      ciao.clear();
      int              max = 0;
      int maxIntens = 0;
      dizCount.clear();
      dizR.clear();
      dizG.clear();
      dizB.clear();
      
      // printf("//////////////////////////////////////");
      for (int y = j - rang; y <= j + rang && y < image.height; y++) {
        for (int x = i - rang; x <= i + rang && x < image.width; x++) {
          if (x < 0 || y < 0) {
            continue;
          }
          auto c            = image[{x, y}];
          auto& pixel = graded[{i, j}];
          //f:1 = x
          int r            = c.x * 255;
          int g            = c.y * 255;
          int b            = c.z * 255;
          
          int curIntensity = (int)(double) ((((r+g+b)/3)*intensityLevels)/255);
          if (curIntensity > maxIntens) maxIntens = curIntensity;
          if (dizCount.find(curIntensity) == dizCount.end()) {
          
              dizCount.insert(std::pair<int, int>(curIntensity, 1));
          
          } else {

              dizCount[curIntensity] = dizCount[curIntensity] + 1; 
          }

          if (dizR.find(curIntensity) == dizR.end()) {
            dizR.insert(std::pair<int, int>(curIntensity, r));

          } else {
              dizR[curIntensity] = dizR[curIntensity] + r;
          }
          if (dizG.find(curIntensity) == dizG.end()) {
            dizG.insert(std::pair<int, int>(curIntensity, g));

          } else {
            dizG[curIntensity] = dizG[curIntensity] + g;
          }
          if (dizB.find(curIntensity) == dizB.end()) {
            dizB.insert(std::pair<int, int>(curIntensity, b));

          } else {
            dizB[curIntensity] = dizB[curIntensity] + b;
          }
          
        }
      }
      int ValueMax = 0;
      int keyMax = 0;
      for (auto const&  [ key, val ] : dizCount) {
        
          if (ValueMax < val) {
          ValueMax = val;
          keyMax   = key;
        }
      }
      float finalR = dizR[keyMax] / ValueMax;
      float finalG = dizG[keyMax] / ValueMax;
      float finalB = dizB[keyMax] / ValueMax;
      //int : 255 = float : 1
     // printf("Prima: %f %f %f\n", pixel.x, pixel.y, pixel.z);
      auto& pixel = graded[{i, j}];
      pixel.x = finalR / 255;
      pixel.y = finalG / 255;
      pixel.z = finalB / 255;
   //   printf("Dopo: %f %f %f\n", pixel.x, pixel.y, pixel.z);
      
    }
  }
  return graded;
}



color_image gaussianBlur(color_image& image, const grade_params& params) {

    float sigma = params.sigma;
    int   nVolte     = params.gaussBlurIntensity;
    int   esec   = 0;

    float values[11];
    float Sr             = 0;
    float Sg             = 0;
    float Sb             = 0;
    float normalizeValue = 0;
    for (int i = -5; i < 6; ++i) {
      values[i + 5] = exp(-(pow((float)i, 2) / (2 * pow(sigma, 2)))) /
                      sqrt(2 * pi * pow(sigma, 2));

      normalizeValue += values[i + 5];
    }
    for (int i = 0; i < 11; ++i) values[i] = values[i] / normalizeValue;
  
  
  while (esec < nVolte){
  for (int y = 0; y < image.height; ++y) {
      for (int x = 0; x < image.width; ++x) {
        Sr = 0;
        Sg = 0;
        Sb = 0;
        int cont = 0;
        for (int i = -5; i < 6; ++i) {
          auto val = values[cont];
          if (x + i >= 0 && x + i < image.width) {

            Sr += image[{x + i, y}].x * val;
            Sg += image[{x + i, y}].y * val;
            Sb += image[{x + i, y}].z * val;
          } else {
            
            Sr += image[{x, y}].x * val;
            Sg += image[{x, y}].y * val;
            Sb += image[{x, y}].z * val;
            
          }
          cont = cont + 1;
        }
        image[{x, y}].x = Sr;
        image[{x, y}].y = Sg;
        image[{x, y}].z = Sb;
      }
    }

  for (int y = 0; y < image.height; ++y) {
      for (int x = 0; x < image.width; ++x) {
        Sr = 0;
        Sg = 0;
        Sb = 0;
        int cont = 0;
        for (int i = -5; i < 6; ++i) {
          auto val = values[cont];
          if (y + i >= 0 && y + i < image.height) {
            Sr += image[{x, y + i}].x * val;
            Sg += image[{x, y + i}].y * val;
            Sb += image[{x, y + i}].z * val;
          } else {
            Sr += image[{x, y}].x * val;
            Sg += image[{x, y}].y * val;
            Sb += image[{x, y}].z * val;
          }
          cont = cont + 1;
        }
        image[{x, y}].x = Sr;
        image[{x, y}].y = Sg;
        image[{x, y}].z = Sb;
      }
    }
    esec  = esec + 1;
  }
  return image;
}



color_image trueAnaglyph(const color_image& image, const grade_params& params){
  

  auto graded = image;

  int   w      = graded.width;
  int   h      = graded.height;
  int   k      = params.anaglyphOffSet;
  

  for (auto j = 0; j < h; j++) {
    for (auto i = 0; i < w; i++) {
      auto pixel = xyz(image[{i, j}]);
      auto& c     = graded[{i, j}];
      int   boh   = k + i;
      if (boh >= image.width) {
        
        c.x = pixel.x * 0.299 + pixel.y * 0.587 + pixel.z * 0.113;
        c.y = 0;    
        c.z = pixel.x * 0.299 + pixel.y * 0.587 + pixel.z * 0.113;
      
      } else {
        auto pix = image[{boh, j}];
      
          c.x         = pixel.x * 0.299 + pixel.y * 0.587 + pixel.z * 0.113;
          c.y         = 0;
          c.z         = pix.x * 0.299 + pix.y * 0.587 + pix.z * 0.113;
      }
    }
  }
  return graded;
  


}



color_image grayAnaglyph(const color_image& image, const grade_params& params) {
  auto graded = image;

  int w = graded.width;
  int h = graded.height;
  int k = params.anaglyphOffSet;

  for (auto j = 0; j < h; j++) {
    for (auto i = 0; i < w; i++) {
      auto  pixel = xyz(image[{i, j}]);
      auto& c     = graded[{i, j}];
      int   boh   = k + i;
      if (boh >= image.width) {
        c.x = pixel.x * 0.299 + pixel.y * 0.587 + pixel.z * 0.113;
        c.y = pixel.x * 0.299 + pixel.y * 0.587 + pixel.z * 0.113;
        c.z = pixel.x * 0.299 + pixel.y * 0.587 + pixel.z * 0.113;

      } else {
        auto pix = image[{boh, j}];

        c.x = pixel.x * 0.299 + pixel.y * 0.587 + pixel.z * 0.113;
        c.y = pix.x * 0.299 + pix.y * 0.587 + pix.z * 0.113;
        c.z = pix.x * 0.299 + pix.y * 0.587 + pix.z * 0.113;
      }
    }
  }
  return graded;
}


color_image colorAnaglyph(const color_image& image, const grade_params& params) {
  auto graded = image;

  int w = graded.width;
  int h = graded.height;
  int k = params.anaglyphOffSet;

  for (auto j = 0; j < h; j++) {
    for (auto i = 0; i < w; i++) {
      auto  pixel = xyz(image[{i, j}]);
      auto& c     = graded[{i, j}];
      int   boh   = k + i;
      if (boh >= image.width) {
        c.x = pixel.x * 1 + pixel.y * 0 + pixel.z * 0;
        c.y = pixel.x * 0 + pixel.y * 1 + pixel.z * 0;
        c.z = pixel.x * 0 + pixel.y * 0 + pixel.z * 1;

      } else {
        auto pix = image[{boh, j}];

        c.x = pixel.x * 1 + pixel.y * 0 + pixel.z * 0;
        c.y = pix.x * 0 + pix.y * 1 + pix.z * 0;
        c.z = pix.x * 0 + pix.y * 0 + pix.z * 1;
      }
    }
  }
  return graded;
}



color_image halfColorAnaglyph(
    const color_image& image, const grade_params& params) {
  auto graded = image;

  int w = graded.width;
  int h = graded.height;
  int k = params.anaglyphOffSet;

  for (auto j = 0; j < h; j++) {
    for (auto i = 0; i < w; i++) {
      auto  pixel = xyz(image[{i, j}]);
      auto& c     = graded[{i, j}];
      int   boh   = k + i;
      if (boh >= image.width) {
        c.x = pixel.x * 0.299 + pixel.y * 0.587 + pixel.z * 0.113;
        c.y = pixel.x * 0 + pixel.y * 1 + pixel.z * 0;
        c.z = pixel.x * 0 + pixel.y * 0 + pixel.z * 1;

      } else {
        auto pix = image[{boh, j}];

        c.x = c.x = pixel.x * 0.299 + pixel.y * 0.587 + pixel.z * 0.113;
        c.y = pix.x * 0 + pix.y * 1 + pix.z * 0;
        c.z = pix.x * 0 + pix.y * 0 + pix.z * 1;
      }
    }
  }
  return graded;
}





color_image optimizedColorAnaglyph(
    const color_image& image, const grade_params& params) {
  auto graded = image;

  int w = graded.width;
  int h = graded.height;
  int k = params.anaglyphOffSet;

  for (auto j = 0; j < h; j++) {
    for (auto i = 0; i < w; i++) {
      auto  pixel = xyz(image[{i, j}]) * (1 / 1.5);
      auto& c     = graded[{i, j}];
      int   boh   = k + i;
      if (boh >= image.width) {
        c.x = pixel.x * 0 + pixel.y * 0.7 + pixel.z * 0.3;
        c.y = pixel.x * 0 + pixel.y * 1 + pixel.z * 0;
        c.z = pixel.x * 0 + pixel.y * 0 + pixel.z * 1;

      } else {
        auto pix = xyz(image[{boh, j}]) * (1/1.5);

        c.x = c.x = pixel.x * 0 + pixel.y * 0.7 + pixel.z * 0.7;
        c.y       = pix.x * 0 + pix.y * 1 + pix.z * 0;
        c.z       = pix.x * 0 + pix.y * 0 + pix.z * 1;
      }
    }
  }
  return graded;
}





color_image sketch(const color_image& image, const grade_params& params) {
  //https://www.shadertoy.com/view/7scSRX
    auto graded = image;
  for (int j : range(image.height)) {
    for (int i : range(image.width)) {
      auto& c = graded[{i, j}];
      

      vec2f size = {image.width, image.height};
      vec2f ij   = {i, j};
      vec2f uv   = ij / size;
      float  r    = sin(c.x * 6.28 * .9);
      auto  silhouette = vec4f{1, 1, 1} * (smoothstep(0., .8, 1. - pow(r, 6)));
      auto   col = lerp(c * silhouette, silhouette , float(sin(params.boh * .8) * .5 + .5));
      c.x      = col.x;
      c.y      = col.y;
      c.z      = col.z;

    }
  }
   


  return graded;
}













const float light = 1.15;
const float dark  = 0.85;

bool tri(const vec2f p1, const vec2f p2, const vec2f p3, const vec2f p) {
  float alpha = ((p2.y - p3.y) * (p.x - p3.x) + (p3.x - p2.x) * (p.y - p3.y)) /
                ((p2.y - p3.y) * (p1.x - p3.x) + (p3.x - p2.x) * (p1.y - p3.y));

  float beta = ((p3.y - p1.y) * (p.x - p3.x) + (p1.x - p3.x) * (p.y - p3.y)) /
               ((p2.y - p3.y) * (p1.x - p3.x) + (p3.x - p2.x) * (p1.y - p3.y));

  float gamma = 1.0 - alpha - beta;

  return (alpha >= 0.0 && beta >= 0.0 && gamma >= 0.0);
}



color_image triangoli(const color_image& image, const grade_params& params) {
  // https://www.shadertoy.com/view/XscSWH
    
  return image;
  
}

color_image CRTTV(const color_image& image, const grade_params& params) {

  auto graded = image;
  int w = graded.width;
  int h = graded.height;


  for (auto j = 0; j < h; j++) {
    for (auto i = 0; i < w; i++) {
      auto& c = graded[{i, j}];
      
      vec2f f = {i, j};
      float iTime = 2.0f;
      vec2f r = {image.width, image.height} , u = f / r;
      u.x += sin(iTime + u.y * r.y) * tan(iTime * u.y * 100.) *
             sin(u.y + iTime) * 0.007;
      float uu = u.y;
      float sinsin = sin(2.0f + uu * 200.0f);
      float rr = eval_image(image,
          
          u + vec2f{sinsin * length(u - 0.5) * 0.07f, 0.0f})
                     .x;
      float gg = eval_image(image, u).y;
      float bb = eval_image(image,
          u + vec2f{sin(iTime + u.y * 200.0f) * length(u - 0.5) * 0.05f, 0.0f})
                     .z;


      c.x = rr;
      c.y = gg;
      c.z = bb;

    }
  }

  return graded;

}





color_image grade_image(const color_image& image, const grade_params& params) {


  

 
  Intensity  potenza = Intensity::MEDIUM ;
  auto      ciao    = image;
  
  
  ciao = homework24(image, params);
  
  //if (channelSegregation)  ciao = channelSegregation(image, params);
 if (params.inverti) ciao = inverti(ciao, params);
  if (params.blur) ciao = blur(ciao, params);
  if (params.cut) ciao = cut(ciao, params);
  if (params.edgeDetection) {
    if (params.livelloEdgeDetection == 0) {
      potenza = Intensity::LOW;
    } else if (params.livelloEdgeDetection == 1) {
      potenza = Intensity::MEDIUM;
    } else if (params.livelloEdgeDetection == 2) {
      potenza = Intensity::STRONG;
    }
    ciao = edgeDetection(ciao, params, potenza);
  }

  if (params.soften) {
    if (params.livelloSoften == 0) {
      potenza = Intensity::LOW;
    } else if (params.livelloSoften == 1) {
      potenza = Intensity::MEDIUM;
    } else if (params.livelloSoften == 2) {
      potenza = Intensity::STRONG;
    }
    ciao = softenImage(ciao, params, potenza);
  }

  if (params.enhanceDetails) ciao = EnhancedDetails(ciao, params);

  if (params.emboss) ciao = emboss(ciao, params);
  
   if(params.oilPainting)  ciao = oil(ciao, params);

  
  

  if(params.gauss) ciao =gaussianBlur(ciao, params);
  if (params.trueAnaglyph) ciao = trueAnaglyph(ciao, params);
  if (params.grayAnaglyph) ciao = grayAnaglyph(ciao, params);
  if (params.colorAnaglyph) ciao = colorAnaglyph(ciao, params);
  if (params.halfColorAnaglyph) ciao = halfColorAnaglyph(ciao, params);
  if (params.optimizedColorAnaglyph) ciao = optimizedColorAnaglyph(ciao, params);
  if (params.sketch) ciao = sketch(ciao, params);
  if (params.CRTTV) ciao = CRTTV(ciao, params);


  
  return ciao;

 
}
}


  // namespace yocto-