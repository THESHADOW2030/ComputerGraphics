//
// Implementation for Yocto/RayTrace.
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

#include "yocto_raytrace.h"
#include <yocto/yocto_cli.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_parallel.h>
#include <yocto/yocto_sampling.h>
#include <yocto/yocto_shading.h>
#include <yocto/yocto_shape.h>

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SCENE EVALUATION
// -----------------------------------------------------------------------------
namespace yocto {

// Generates a ray from a camera for yimg::image plane coordinate uv and
// the lens coordinates luv.
static ray3f eval_camera(const camera_data& camera, const vec2f& uv) {
  auto film = camera.aspect >= 1
                  ? vec2f{camera.film, camera.film / camera.aspect}
                  : vec2f{camera.film * camera.aspect, camera.film};
  auto q    = transform_point(camera.frame,
      {film.x * (0.5f - uv.x), film.y * (uv.y - 0.5f), camera.lens});
  auto e    = transform_point(camera.frame, {0, 0, 0});
  return {e, normalize(e - q)};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING
// -----------------------------------------------------------------------------
namespace yocto {

double reflectance(double cosine, double ref_idx) {
  // Use Schlick's approximation for reflectance.
  auto r0 = (1 - ref_idx) / (1 + ref_idx);
  r0      = r0 * r0;
  return r0 + (1 - r0) * pow((1 - cosine), 5);
}

// Raytrace renderer.
static vec4f shade_raytrace(const scene_data& scene, const bvh_scene& bvh,
    const ray3f& ray, int bounce, rng_state& rng,
    const raytrace_params& params) {
  // YOUR CODE GOES HERE ----
  auto isec = intersect_bvh(bvh, scene, ray);
  if (!isec.hit) return rgb_to_rgba(eval_environment(scene, ray.d));
  auto& instance = scene.instances[isec.instance];
  auto& shape    = scene.shapes[instance.shape];
  auto& material = eval_material(scene, instance, isec.element, isec.uv);
  auto  normal   = transform_direction(
      instance.frame, eval_normal(shape, isec.element, isec.uv));
  auto position = transform_point(
      instance.frame, eval_position(shape, isec.element, isec.uv));
  auto radiance = material.emission;
  auto outgoing = -ray.d;

  /*                     */
  if (!shape.points.empty()) {
    normal = -ray.d;
  } else if (!shape.lines.empty()) {
    normal = orthonormalize(-ray.d, normal);
  } else if (!shape.triangles.empty()) {
    if (dot(-ray.d, normal) < 0) {
      normal = -normal;
    }
  }
  /*                     */

  if (rand1f(rng) < 1 - material.opacity)
    return shade_raytrace(scene, bvh, ray3f{position, ray.d}, bounce + 1, rng,
        params);  // opacita
  if (bounce >= params.bounces) return rgb_to_rgba(radiance);
  // non ho scritto il look up delle texture
  auto color = material.color;
  switch (material.type) {
    case (material_type::matte): {
      // handle diffuse
      auto incoming = sample_hemisphere_cos(normal, rand2f(rng));
      radiance += color *
                  rgba_to_rgb(shade_raytrace(scene, bvh,
                      ray3f{position, incoming}, bounce + 1, rng, params));
    } break;
    case (material_type::reflective): {
      if (material.roughness>0.0f)  // Se è ruvido
      {
        //<handle rough metals>
        auto halfway = sample_hemisphere_cospower(
            2 / (material.roughness * material.roughness), normal, rand2f(rng));
        vec3f incoming = reflect(outgoing, halfway);
        radiance += fresnel_schlick(color, halfway, outgoing) *
                    rgba_to_rgb(shade_raytrace(scene, bvh,
                        ray3f{position, incoming}, bounce + 1, rng, params));
      } else {
        //<handle polished metals>
        vec3f incoming = reflect(outgoing, normal);
        radiance += fresnel_schlick(color, normal, outgoing) *
                    rgba_to_rgb(shade_raytrace(scene, bvh,
                        ray3f{position, incoming}, bounce + 1, rng, params));
      }
    } break;
    case (material_type::glossy): {  // rough plastic

      auto halfway = sample_hemisphere_cospower(
          2 / (material.roughness * material.roughness), normal, rand2f(rng));
      auto fs = fresnel_schlick(vec3f{0.04, 0.04, 0.04}, halfway, outgoing);
      if (rand1f(rng) < fs.x) {
        auto incoming = reflect(outgoing, halfway);
        radiance += rgba_to_rgb(shade_raytrace(
            scene, bvh, ray3f{position, incoming}, bounce + 1, rng, params));
      } else {
        auto incoming = sample_hemisphere_cos(normal, rand2f(rng));
        radiance += color *
                    rgba_to_rgb(shade_raytrace(scene, bvh,
                        ray3f{position, incoming}, bounce + 1, rng, params));
      }
    } break;
    case (material_type::transparent): {  // polished dialetrics

    
      if (rand1f(rng) <
          fresnel_schlick(vec3f{0.04, 0.04, 0.04}, normal, outgoing)
              .x) {  // da correggere
        auto incoming = reflect(outgoing, normal);
        radiance += rgba_to_rgb(shade_raytrace(
            scene, bvh, ray3f{position, incoming}, bounce + 1, rng, params));
      } else {
        auto incoming = -outgoing;
        radiance += material.color *
                    rgba_to_rgb(shade_raytrace(scene, bvh,
                        ray3f{position, incoming}, bounce + 1, rng, params));
      }
    } break;
    case material_type::refractive: {
      auto  ior       = material.ior;
      auto  normal2 = normal;
      auto  ior2   = ior;
      if (dot(outgoing, normal) < 0) {
        normal2 = -normal;
        ior2   = (1 / ior);
      }
      if (rand1f(rng) < fresnel_schlick({0.04, 0.04, 0.04}, normal, outgoing).x) {
        auto incoming = reflect(outgoing, normal2);
        radiance += rgba_to_rgb(shade_raytrace(
            scene, bvh, ray3f{position, incoming}, bounce + 1, rng, params));
      } else {
        auto incoming = refract(outgoing, normal2, ior2);
        radiance += color *
                    rgba_to_rgb(shade_raytrace(scene, bvh,
                        ray3f{position, incoming}, bounce + 1, rng, params));
       
      }
    }
  }
      return rgb_to_rgba(radiance);
}




// Matte renderer.;
static vec4f shade_matte(const scene_data& scene, const bvh_scene& bvh,
    const ray3f& ray, int bounce, rng_state& rng,
    const raytrace_params& params) {
  // YOUR CODE GOES HERE ----
  return {0, 0, 0, 0};
}
static vec4f cel_shading(const scene_data& scene, const bvh_scene& bvh,
    const ray3f& ray, int bounce, rng_state& rng,
    const raytrace_params& params) {
  // YOUR CODE GOES HERE ----
  auto isec = intersect_bvh(bvh, scene, ray);
  if (!isec.hit) return rgb_to_rgba(eval_environment(scene, ray.d));
  auto& instance = scene.instances[isec.instance];
  auto& shape    = scene.shapes[instance.shape];

  auto material = eval_material(scene, instance, isec.element, isec.uv);
  auto normal   = transform_direction(
      instance.frame, eval_normal(shape, isec.element, isec.uv));
  auto position = transform_point(
      instance.frame, eval_position(shape, isec.element, isec.uv));
  auto radiance = material.emission;
  auto outgoing = -ray.d;

  vec3f light_pos = {0.6043607234954834, 0.800000011920929, 0.800000011920929};
  float nDotl     = dot(light_pos, normal);
  auto  color     = material.color;

  float lightIntensity = smoothstep(0.0f, 0.01f, nDotl);

  // VA BENE FINO A QUI

  vec4f specularColor     = vec4f{1, 1, 1, 1} * 0.9;
  specularColor.w         = 1;
  float glossiness        = 32;
  auto  halfVector        = normalize(-ray.d + light_pos);
  float NdotH             = dot(normal, halfVector);
  float specularIntensity = pow(
      NdotH * lightIntensity, glossiness * glossiness);
  float specularIntensitySmooth = smoothstep(0.005f, 0.01f, specularIntensity);
  vec4f specular                = specularIntensitySmooth * specularColor;

  float rimDot = 1 - dot(-ray.d, normal);
  auto  rimAmount = 0.716f;  
  float  rimIntensity       = rimDot * pow(nDotl, 0.1f);
  rimIntensity = smoothstep(rimAmount - 0.01f, rimAmount + 0.01f, rimIntensity);
  auto  rim          = rimIntensity * vec3f{1, 1, 1};
  if (!instance.material == 0) {
    color = color * (lightIntensity * 0.3 + eval_environment(scene, light_pos) +
                        rgba_to_rgb(specular) + rim);
  
  } else {
    
  color = color * (eval_environment(scene, light_pos));
  }

  
  auto isec2 = intersect_bvh(bvh, scene, ray3f{position, light_pos});
  if (isec2.hit &&
      scene.materials[scene.instances[isec2.instance].material].emission.x ==
          0 &&
      scene.materials[scene.instances[isec2.instance].material].emission.y ==
          0 &&
      scene.materials[scene.instances[isec2.instance].material].emission.z ==
          0) {
    color = color * 0.8;
    
  }

  return rgb_to_rgba(color);
}



static vec4f matcap(const scene_data& scene, const bvh_scene& bvh,
    const ray3f& ray, int bounce, rng_state& rng,
    const raytrace_params& params) {
  // YOUR CODE GOES HERE ----
  auto isec = intersect_bvh(bvh, scene, ray);
  if (!isec.hit) return rgb_to_rgba(eval_environment(scene, ray.d));
  auto& instance = scene.instances[isec.instance];
  auto& shape    = scene.shapes[instance.shape];
  // auto& material = eval_material(scene, instance, isec.element, isec.uv);
  auto normal = transform_direction(
      instance.frame, eval_normal(shape, isec.element, isec.uv));
  auto position = transform_point(
      instance.frame, eval_position(shape, isec.element, isec.uv));
  if (instance.material == 0) {
    auto ciao = eval_material(scene, instance, isec.element, isec.uv);

    return rgb_to_rgba(ciao.color);
  }

  auto  reflected = reflect(ray.o, normal);
  float m         = 2.8284271247461903 * sqrt(reflected.z + 1.0);
  auto  uv        = vec2f{reflected.x, reflected.y} / m + 0.5;

  auto  text      = eval_texture(scene, isec.instance, uv);
  auto material = rgb_to_rgba(scene.materials[instance.material].color) * text;
  return material;
}



// Eyelight renderer.
static vec4f shade_eyelight(const scene_data& scene, const bvh_scene& bvh,
    const ray3f& ray, int bounce, rng_state& rng,
    const raytrace_params& params) {
  // YOUR CODE GOES HERE ----
  auto isec = intersect_bvh(bvh, scene, ray);
  if (!isec.hit) {
    return zero4f;
  }
  auto& instance = scene.instances[isec.instance];
  auto& material = scene.materials[instance.material];
  auto& shape    = scene.shapes[instance.shape];
  auto  normal   = transform_direction(
      instance.frame, eval_normal(shape, isec.element, isec.uv));
  auto ciao = material.color * dot(normal, -ray.d);
  auto idk  = vec4f{ciao.x, ciao.y, ciao.z};
  return idk;
}

static vec4f shade_normal(const scene_data& scene, const bvh_scene& bvh,
    const ray3f& ray, int bounce, rng_state& rng,
    const raytrace_params& params) {
  // YOUR CODE GOES HERE ----
  auto intersezione = intersect_bvh(bvh, scene, ray);
  if (!intersezione.hit) {
    return zero4f;
  }
  auto  isec     = intersect_bvh(bvh, scene, ray);
  auto& instance = scene.instances[isec.instance];
  auto& shape    = scene.shapes[instance.shape];
  auto  normal   = transform_direction(
      instance.frame, eval_normal(shape, isec.element, isec.uv));

  auto boh = vec4f{normal.x, normal.y, normal.z, 1} * 0.5 + 0.5;
  return boh;
}

static vec4f shade_texcoord(const scene_data& scene, const bvh_scene& bvh,
    const ray3f& ray, int bounce, rng_state& rng,  // fatto
    const raytrace_params& params) {
  // YOUR CODE GOES HERE ----

  auto intersezione = intersect_bvh(bvh, scene, ray);

  auto  isec     = intersect_bvh(bvh, scene, ray);
  auto& instance = scene.instances[isec.instance];
  auto& shape    = scene.shapes[instance.shape];
  auto  normal   = transform_direction(
      instance.frame, eval_normal(shape, isec.element, isec.uv));

  auto ciao = eval_texcoord(shape, isec.element, isec.uv);
  auto bohX = fmod(ciao.x, 1);
  auto bohY = fmod(ciao.y, 1);
  return vec4f{bohX, bohY, 0, 1};
}

static vec4f shade_color(const scene_data& scene, const bvh_scene& bvh,
    const ray3f& ray, int bounce, rng_state& rng,  // fatto
    const raytrace_params& params) {
  // YOUR CODE GOES HERE ----

  auto intersezione = intersect_bvh(bvh, scene, ray);
  auto ciao         = intersezione.instance;
  auto boh          = scene.instances.at(ciao);
  auto hola         = scene.materials.at(boh.material);
  auto grad         = hola.color;
  return vec4f{grad.x, grad.y, grad.z, 1};
}

// Trace a single ray from the camera using the given algorithm.
using raytrace_shader_func = vec4f (*)(const scene_data& scene,
    const bvh_scene& bvh, const ray3f& ray, int bounce, rng_state& rng,
    const raytrace_params& params);
static raytrace_shader_func get_shader(const raytrace_params& params) {
  
   
   if (params.toon) return cel_shading;
  if (params.matcap) return matcap;
  switch (params.shader) {
    
    case raytrace_shader_type::raytrace: return shade_raytrace;
    case raytrace_shader_type::matte: return shade_matte;
    case raytrace_shader_type::eyelight: return shade_eyelight;
    case raytrace_shader_type::normal: return shade_normal;
    case raytrace_shader_type::texcoord: return shade_texcoord;
    case raytrace_shader_type::color: return shade_color;
    
    default: {
      throw std::runtime_error("sampler unknown");
      return nullptr;
    }
  }
}

// Build the bvh acceleration structure.
bvh_scene make_bvh(const scene_data& scene, const raytrace_params& params) {
  return make_bvh(scene, false, false, params.noparallel);
}

// Init a sequence of random number generators.
raytrace_state make_state(
    const scene_data& scene, const raytrace_params& params) {
  auto& camera = scene.cameras[params.camera];
  auto  state  = raytrace_state{};
  if (camera.aspect >= 1) {
    state.width  = params.resolution;
    state.height = (int)round(params.resolution / camera.aspect);
  } else {
    state.height = params.resolution;
    state.width  = (int)round(params.resolution * camera.aspect);
  }
  state.samples = 0;
  state.image.assign(state.width * state.height, {0, 0, 0, 0});
  state.hits.assign(state.width * state.height, 0);
  state.rngs.assign(state.width * state.height, {});
  auto rng_ = make_rng(1301081);
  for (auto& rng : state.rngs) {
    rng = make_rng(961748941ull, rand1i(rng_, 1 << 31) / 2 + 1);
  }
  return state;
}

// Progressively compute an image by calling trace_samples multiple times.
void raytrace_samples(raytrace_state& state, const scene_data& scene,
    const bvh_scene& bvh, const raytrace_params& params) {
  if (state.samples >= params.samples) return;
  auto& camera = scene.cameras[params.camera];
  auto  shader = get_shader(params);
  state.samples += 1;
  if (params.samples == 1) {
    for (auto idx = 0; idx < state.width * state.height; idx++) {
      auto i = idx % state.width, j = idx / state.width;
      auto u = (i + 0.5f) / state.width, v = (j + 0.5f) / state.height;
      auto ray      = eval_camera(camera, {u, v});
      auto radiance = shader(scene, bvh, ray, 0, state.rngs[idx], params);
      if (!isfinite(radiance)) radiance = {0, 0, 0};
      state.image[idx] += radiance;
      state.hits[idx] += 1;
    }
  } else if (params.noparallel) {
    for (auto idx = 0; idx < state.width * state.height; idx++) {
      auto i = idx % state.width, j = idx / state.width;
      auto u        = (i + rand1f(state.rngs[idx])) / state.width,
           v        = (j + rand1f(state.rngs[idx])) / state.height;
      auto ray      = eval_camera(camera, {u, v});
      auto radiance = shader(scene, bvh, ray, 0, state.rngs[idx], params);
      if (!isfinite(radiance)) radiance = {0, 0, 0};
      state.image[idx] += radiance;
      state.hits[idx] += 1;
    }
  } else {
    parallel_for(state.width * state.height, [&](int idx) {
      auto i = idx % state.width, j = idx / state.width;
      auto u        = (i + rand1f(state.rngs[idx])) / state.width,
           v        = (j + rand1f(state.rngs[idx])) / state.height;
      auto ray      = eval_camera(camera, {u, v});
      auto radiance = shader(scene, bvh, ray, 0, state.rngs[idx], params);
      if (!isfinite(radiance)) radiance = {0, 0, 0};
      state.image[idx] += radiance;
      state.hits[idx] += 1;
    });
  }
}

// Check image type
static void check_image(
    const color_image& image, int width, int height, bool linear) {
  if (image.width != width || image.height != height)
    throw std::invalid_argument{"image should have the same size"};
  if (image.linear != linear)
    throw std::invalid_argument{
        linear ? "expected linear image" : "expected srgb image"};
}

// Get resulting render
color_image get_render(const raytrace_state& state) {
  auto image = make_image(state.width, state.height, true);
  get_render(image, state);
  return image;
}
void get_render(color_image& image, const raytrace_state& state) {
  check_image(image, state.width, state.height, true);
  auto scale = 1.0f / (float)state.samples;
  for (auto idx = 0; idx < state.width * state.height; idx++) {
    image.pixels[idx] = state.image[idx] * scale;
  }
}

}  // namespace yocto
