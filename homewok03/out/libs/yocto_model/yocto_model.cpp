//
// Implementation for Yocto/Model
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_model.h"

#include <yocto/yocto_sampling.h>

#include "ext/perlin-noise/noise1234.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::string;
using std::vector;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR EXAMPLE OF PROCEDURAL MODELING
// -----------------------------------------------------------------------------
namespace yocto {

 vec3f hash3(vec2f p);
vec2f fract(vec2f x);

 vec2f hash2(vec2f p) {

     auto value = vec2f{dot(p, vec2f{127.1, 311.7}), dot(p, vec2f{269.5, 183.3})};
     value.x    = sin(value.x) * 43758.5453;
     value.y    = sin(value.y) * 43758.5453;

     return fract(value);

 }

 vec3f cellnoise(vec3f x) {

     vec2f n = vec2f{floor(x.x), floor(x.y)};
   vec2f f = vec2f{x.x - n.x, x.y - n.y};

    vec2f mg, mr;

   float md = 8.0;
    for (int j = -1; j <= 1; j++)
     for (int i = -1; i <= 1; i++) {

         auto g = vec2f{float(i), float(j)};
         auto o = hash2(n + g);
         auto r = g + o - f;
         float d = dot(r, r);
         if (d < md) {
           md = d;
           mr = r;
           mg = g;
         }
     }

    
    md = 8.0;
    for (int j = -2; j <= 2; j++)
     for (int i = -2; i <= 2; i++) {
       auto  g = mg + vec2f{float(i), float(j)};
       auto  o = hash2(n + g);
       auto  r = g + o - f;
       float d = dot(r, r);
       if (dot(mr - r, mr - r) > 0.00001) {
         md = min(md, dot(0.5 * (mr + r), normalize(r - mr)));
       }
     }


    return vec3f{md, mr.x, mr.y};
 }


































vec3f fract(vec3f x) {
  auto y = vec3f{floor(x.x), floor(x.y), floor(x.z)};
  return x - y;
}
vec2f fract(vec2f x) {
  auto y = vec2f{floor(x.x), floor(x.y)};
  return x - y;
}







float extraSmoothVoronoi(vec2f x, float falloff) {


    auto  p = vec2f{floor(x.x), floor(x.y)};
    vec2f f = fract(x);

    float res = 0.0;
    for (int j = -1; j <= 1; j++)
      for (int i = -1; i <= 1; i++) {
        vec2f b = vec2f{float(i), float(j)};
        auto  ciao = hash3(p + b);
        vec2f r    = b - f + vec2f{ciao.x, ciao.y};
        float d = length(r);

        res += exp2(-falloff * d);
      }
    return -(1.0 / falloff) * log2(res);


}



















    vec3f hash3(vec2f p) {
      auto q = vec3f{dot(p, vec2f{127.1f, 311.7f}),
          dot(p, vec2f{269.5f, 183.3f}), dot(p, vec2f{419.2f, 371.9f})};
      return fract(vec3f{(sin(q.x) * 43758.5453f),
          (sin(q.y) * 43758.5453f),
          (sin(q.z) * 43758.5453f)});
    }

    float extraVoronoise(vec2f x, float u, float v) {
      auto  p = vec2f{floor(x.x), floor(x.y)};
      vec2f f = fract(x);

      float k  = 1.0 + 63.0 * pow(1.0 - v, 4.0);
      float va = 0.0;
      float wt = 0.0;
      for (int j = -2; j <= 2; j++)
        for (int i = -2; i <= 2; i++) {
            auto g = vec2f{float(i), float(j)};
            auto o = hash3(p + g) * vec3f{u, u, 1.0};
            auto r = g - f + vec2f{o.x, o.y};
            float d = dot(r, r);
            float w = pow(1.0f - smoothstep(0.0f, 1.414f, sqrt(d)), k);
            va += w * o.z;
            wt += w;

        }
    
      return va/wt;

 }

float noise(const vec3f& p) { return ::noise3(p.x, p.y, p.z); }
vec2f noise2(const vec3f& p) {
  return {noise(p + vec3f{0, 0, 0}), noise(p + vec3f{3, 7, 11})};
}
vec3f noise3(const vec3f& p) {
  return {noise(p + vec3f{0, 0, 0}), noise(p + vec3f{3, 7, 11}),
      noise(p + vec3f{13, 17, 19})};
}
float fbm(const vec3f& p, int octaves) {
  auto sum    = 0.0f;
  auto weight = 1.0f;
  auto scale  = 1.0f;
  for (auto octave = 0; octave < octaves; octave++) {
    sum += weight * fabs(noise(p * scale));
    weight /= 2;
    scale *= 2;
  }
  return sum;
}
float turbulence(const vec3f& p, int octaves) {
  auto sum    = 0.0f;
  auto weight = 1.0f;
  auto scale  = 1.0f;
  for (auto octave = 0; octave < octaves; octave++) {
    sum += weight * fabs(noise(p * scale));
    weight /= 2;
    scale *= 2;
  }
  return sum;
}
float ridge(const vec3f& p, int octaves) {
  auto sum    = 0.0f;
  auto weight = 0.5f;
  auto scale  = 1.0f;
  for (auto octave = 0; octave < octaves; octave++) {
    sum += weight * (1 - fabs(noise(p * scale))) * (1 - fabs(noise(p * scale)));
    weight /= 2;
    scale *= 2;
  }
  return sum;
}

void add_polyline(shape_data& shape, const vector<vec3f>& positions,
    const vector<vec4f>& colors, float thickness = 0.0001f) {
  auto offset = (int)shape.positions.size();
  shape.positions.insert(
      shape.positions.end(), positions.begin(), positions.end());
  shape.colors.insert(shape.colors.end(), colors.begin(), colors.end());
  shape.radius.insert(shape.radius.end(), positions.size(), thickness);
  for (auto idx = 0; idx < positions.size() - 1; idx++) {
    shape.lines.push_back({offset + idx, offset + idx + 1});
  }
}

void sample_shape(vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, const shape_data& shape, int num) {
  auto triangles  = shape.triangles;
  auto qtriangles = quads_to_triangles(shape.quads);
  triangles.insert(triangles.end(), qtriangles.begin(), qtriangles.end());
  auto cdf = sample_triangles_cdf(triangles, shape.positions);
  auto rng = make_rng(19873991);
  for (auto idx = 0; idx < num; idx++) {
    auto [elem, uv] = sample_triangles(cdf, rand1f(rng), rand2f(rng));
    auto q          = triangles[elem];
    positions.push_back(interpolate_triangle(
        shape.positions[q.x], shape.positions[q.y], shape.positions[q.z], uv));
    normals.push_back(normalize(interpolate_triangle(
        shape.normals[q.x], shape.normals[q.y], shape.normals[q.z], uv)));
    if (!texcoords.empty()) {
      texcoords.push_back(interpolate_triangle(shape.texcoords[q.x],
          shape.texcoords[q.y], shape.texcoords[q.z], uv));
    } else {
      texcoords.push_back(uv);
    }
  }
}

void make_terrain(shape_data& shape, const terrain_params& params) { 

    for (int i = 0; i < shape.positions.size(); i++) {
    
        auto position = shape.positions[i];

      position += shape.normals[i] *
                  ridge(position * params.scale, params.octaves) *
                  params.height *
                  (1.f - length(position - params.center) / params.size);


        float perc = position.y / params.height;
      if (perc <= 0.3f)
        shape.colors.push_back(params.bottom);
      else if (perc > 0.6)
        shape.colors.push_back(params.top);
      else
        shape.colors.push_back(params.middle);
      shape.positions[i] = position;
    
    }
    compute_normals(shape.normals, shape);


}



void make_displacementVoronoise(
    shape_data& shape, const displacement_params& params);

void make_displacementSmoothVoronoi(
    shape_data& shape, const displacement_params& params);

void make_displacementCellNoise(
    shape_data& shape, const displacement_params& params);

void make_displacement(
    shape_data& shape, const displacement_params& params) {
  // YOUR CODE GOES HERE
  if (params.smoothVoronoi) return make_displacementSmoothVoronoi(shape, params);
  if (params.cellnoise) return make_displacementCellNoise(shape, params);
  if (params.voronoise) return make_displacementVoronoise(shape, params);
  for (int i = 0; i < shape.positions.size(); i++) {
    auto position = shape.positions[i];
    auto last_p   = position;
    position +=
        shape.normals[i] *
        (turbulence(position * params.scale, params.octaves) * params.height);
    shape.colors.push_back(interpolate_line(
        params.bottom, params.top, distance(position, last_p) / params.height));
    shape.positions[i] = position;
  }
  compute_normals(shape.normals, shape);
}



void make_displacementSmoothVoronoi(shape_data& shape, const displacement_params& params) {
  // YOUR CODE GOES HERE
  for (int i = 0; i < shape.positions.size(); i++) {
    auto position = shape.positions[i];
    auto last_p   = position;
    position +=
        shape.normals[i] *
        (extraSmoothVoronoi(vec2f{position.x, position.y} * params.scale, params.falloff) *
            params.height);
    shape.colors.push_back(interpolate_line(
        params.bottom, params.top, distance(position, last_p) / params.height));
    shape.positions[i] = position;
  }
  compute_normals(shape.normals, shape);
}



void make_displacementCellNoise(
    shape_data& shape, const displacement_params& params) {
  // YOUR CODE GOES HERE
  for (int i = 0; i < shape.positions.size(); i++) {
    auto position = shape.positions[i];
    auto last_p   = position;
    auto cn       = cellnoise(position * params.scale);
    auto col      = cn.x * (0.5 + 0.5 * sin(64 * cn.x)) * vec3f{1, 1, 1};
    col           = interpolate_line(
        vec3f{1, 0.6, 0.0}, col, smoothstep(0.04f, 0.07f, cn.x));
    float dd = length(vec2f{cn.y, cn.z});

    col = interpolate_line(vec3f{1, 0.6, 0.1}, col, smoothstep(0.0f, 0.12f, dd));
    col += vec3f{1, 0.6, 0.1} * (1.0 - smoothstep(0.0f, 0.04f, dd));
    shape.colors.push_back(rgb_to_rgba(col));
    shape.positions[i] = position;
  }
  compute_normals(shape.normals, shape);
}



void make_displacementVoronoise(shape_data& shape, const displacement_params& params) {
  // YOUR CODE GOES HERE
  for (int i = 0; i < shape.positions.size(); i++) {
    auto position = shape.positions[i];
    auto last_p   = position;
    position +=
        shape.normals[i] *
        (extraVoronoise(vec2f{position.x, position.y} * params.scale, 1, 1) *
            params.height);
    shape.colors.push_back(interpolate_line(
        params.bottom, params.top, distance(position, last_p) / params.height));
    shape.positions[i] = position;
  }
  compute_normals(shape.normals, shape);
}


void make_hair(
    shape_data& hair, const shape_data& shapee, const hair_params& params) {
  // YOUR CODE GOES HERE
  auto rngDensity = make_rng(10123);
  auto shape = shapee;
    int size = shape.positions.size();
  int  debug      = 1;
  sample_shape(
      shape.positions, shape.normals, shape.texcoords, shape, params.num);
    
    float t = params.lenght / params.steps;
  for (int i = size; i < shape.positions.size(); i++) {
      if (rand1f(rngDensity) > params.density) {
      
      continue;
    }
      std::vector<vec3f> positions;
      std::vector<vec4f> colors;
      positions.push_back(shape.positions[i]);
      colors.push_back(params.bottom);
      vec3f normal = shape.normals[i];
      for (int k = 0; k < params.steps; k++) {
        vec3f next = positions[k] + t * normal +
                     noise3(positions[k] * params.scale) * params.strength;
        next.y -= params.gravity;
        
        normal = normalize(next - positions[k]);
        positions.push_back(next);
        colors.push_back(
            interpolate_line(params.bottom, params.top,
                distance(next, positions[0]) / params.lenght));
      }
      colors[params.steps] = params.top;
      add_polyline(hair, positions, colors);
  }
  
   
  auto tang = lines_tangents(hair.lines, hair.positions);
  for (int i = 0; i < tang.size(); i++)
   hair.tangents.push_back(vec4f{tang[i].x, tang[i].y, tang[i].z, 0.f});
}


void make_grass(scene_data& scene, const instance_data& object,
    const vector<instance_data>& grasses, const grass_params& params) {
  // YOUR CODE GOES HERE

  auto rng = make_rng(2030);
  auto rngDensity = make_rng(10123);
  auto shape = scene.shapes[object.shape];  
  int  size  = shape.positions.size();
  int  debug      = 1;
  sample_shape(
      shape.positions, shape.normals, shape.texcoords, shape, params.num);

  for (int i = size; i < shape.positions.size(); i++) {
    if (rand1f(rngDensity) > params.density) {

      
        continue;
    }
    vec3f  position = shape.positions[i];
    auto grass = grasses[rand1i(rng, grasses.size())];
    grass.frame.y = shape.normals[i];
    grass.frame.x = normalize(
        vec3f{1, 0, 0} - dot(vec3f{1, 0, 0}, grass.frame.y) * grass.frame.y);
    grass.frame.z = cross(grass.frame.x, grass.frame.y);
    grass.frame.o = position;
    float rand = 0.9f + rand1f(rng) * 0.1f;
    grass.frame *= scaling_frame(vec3f{rand, rand, rand});
    rand = rand1f(rng) * 2 * pif;
    grass.frame *= rotation_frame(grass.frame.y, rand);
    rand = 0.1f + rand1f(rng) * 0.1f;
    grass.frame *= rotation_frame(grass.frame.z, rand);
    scene.instances.push_back(grass);
  } 
}



}  // namespace yocto
