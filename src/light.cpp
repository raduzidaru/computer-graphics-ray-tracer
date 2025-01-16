#include "light.h"
#include "bvh_interface.h"
#include "config.h"
#include "draw.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "shading.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()


// TODO: Standard feature
// Given a single segment light, transform a uniformly distributed 1d sample in [0, 1),
// into a uniformly sampled position and an interpolated color on the segment light,
// and write these into the reference return values.
// - sample;    a uniformly distributed 1d sample in [0, 1)
// - light;     the SegmentLight object, see `common.h`
// - position;  reference return value of the sampled position on the light
// - color;     reference return value of the color emitted by the light at the sampled position
// This method is unit-tested, so do not change the function signature.
void sampleSegmentLight(const float& sample, const SegmentLight& light, glm::vec3& position, glm::vec3& color)
{
    // TODO: implement this function.
    position = glm::mix(light.endpoint0, light.endpoint1, sample);
    color = glm::mix(light.color0, light.color1, sample);
}

// TODO: Standard feature
// Given a single paralellogram light, transform a uniformly distributed 2d sample in [0, 1),
// into a uniformly sampled position and interpolated color on the paralellogram light,
// and write these into the reference return values.
// - sample;   a uniformly distributed 2d sample in [0, 1)
// - light;    the ParallelogramLight object, see `common.h`
// - position; reference return value of the sampled position on the light
// - color;    reference return value of the color emitted by the light at the sampled position
// This method is unit-tested, so do not change the function signature.
void sampleParallelogramLight(const glm::vec2& sample, const ParallelogramLight& light, glm::vec3& position, glm::vec3& color)
{
    // TODO: implement this function.
    position = light.v0 + sample.x * light.edge01 + sample.y * light.edge02;

    color = light.color0 + sample.x * (light.color1 - light.color0) + sample.y * (light.color2 - light.color0);
}

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return whether
// or not the light is visible from the provided ray/intersection.
// For a description of the method's arguments, refer to 'light.cpp'
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        whether the light is visible (true) or not (false)
// This method is unit-tested, so do not change the function signature.
bool visibilityOfLightSampleBinary(RenderState& state, const glm::vec3& lightPosition, const glm::vec3 &lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return true;
    } else {
        // Shadows are enabled in the renderer
        // TODO: implement this function; currently, the light simply passes through
        Ray lightRay;
        lightRay.origin = ray.origin + ray.t * ray.direction + 1e-5f * glm::normalize(hitInfo.normal);
        lightRay.direction = glm::normalize(lightPosition - lightRay.origin);
        lightRay.t = glm::length(lightPosition - lightRay.origin);

        HitInfo info;
        bool intersection = state.bvh.intersect(state, lightRay, info);
        

        if (!intersection) {
            drawRay(lightRay, lightColor);
        } else {
            drawRay(lightRay, glm::vec3(1.0f, 0.0f, 0.0f));
        }
        return !intersection;
    }
}

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// Use the following blending operation: lightColor = lightColor * kd * (1 - alpha)
// Please reflect within 50 words in your report on why this is incorrect, and illustrate
// two examples of what is incorrect.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSampleTransparency(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
        glm::vec3 res = lightColor;

        glm::vec3 p = ray.origin + ray.direction * ray.t;
        Ray lightRay;
        lightRay.direction = glm::normalize(p-lightPosition);
        lightRay.origin = lightPosition + 0.001f * lightRay.direction;
        lightRay.t = glm::length(lightRay.origin - p) - 0.001f;
        if (state.features.enableDebugDraw) {
            drawRay(lightRay, res);
        }
        HitInfo lightRayInfo;
        while (state.bvh.intersect(state, lightRay, lightRayInfo)) {
            
            res *= lightRayInfo.material.kd * (1.0f - lightRayInfo.material.transparency);

            lightRay.origin += lightRay.direction * (lightRay.t + 0.001f);
            lightRay.t = glm::length(lightRay.origin - p) - 0.001f;
            if (state.features.enableDebugDraw) {
                drawRay(lightRay, res);
            }
        }

       /* if (state.features.enableDebugDraw) {
            drawRay(lightRay,res);
        }*/
        return res;
    
}

// TODO: Standard feature
// Given a single point light, compute its contribution towards an incident ray at an intersection point.
//
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the light is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;   the active scene, feature config, bvh, and a thread-safe sampler
// - light;   the PointLight object, see `common.h`
// - ray;     the incident ray to the current intersection
// - hitInfo; information about the current intersection
// - return;  reflected light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionPointLight(RenderState& state, const PointLight& light, const Ray& ray, const HitInfo& hitInfo)
{
    // TODO: modify this function to incorporate visibility corerctly
    glm::vec3 p = ray.origin + ray.t * ray.direction;
    glm::vec3 l = glm::normalize(light.position - p);
    glm::vec3 v = -ray.direction;
    return computeShading(state, v, l, light.color, hitInfo) * visibilityOfLightSample(state, light.position, light.color, ray, hitInfo);
}

// TODO: Standard feature
// Given a single segment light, compute its contribution towards an incident ray at an intersection point
// by integrating over the segment, taking `numSamples` samples from the light source.
//
// Hint: you can sample the light by using `sampleSegmentLight(state.sampler.next_1d(), ...);`, which
//       you should implement first.
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the sample is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;      the active scene, feature config, bvh, and a thread-safe sampler
// - light;      the SegmentLight object, see `common.h`
// - ray;        the incident ray to the current intersection
// - hitInfo;    information about the current intersection
// - numSamples; the number of samples you need to take
// - return;     accumulated light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionSegmentLight(RenderState& state, const SegmentLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    // TODO: implement this function; repeat numSamples times:
    // - sample the segment light
    // - test the sample's visibility
    // - then evaluate the phong model
    glm::vec3 res(0);

    for (uint32_t i = 1; i <= numSamples; i++) {
        float samplePosition = state.sampler.next_1d();
        glm::vec3 lightPosition, lightColor;
        sampleSegmentLight(samplePosition, light, lightPosition, lightColor);
        
        glm::vec3 p = ray.origin + ray.t * ray.direction;
        glm::vec3 l = glm::normalize(lightPosition - p);
        glm::vec3 v = -ray.direction;
        // Check the visibility of the light sample
        glm::vec3 col = visibilityOfLightSample(state, lightPosition, lightColor, ray, hitInfo);
        if (glm::length(col) > 0.0f) {
           

            glm::vec3 shadingContribution = computeShading(state, v, l, lightColor, hitInfo);
            res += shadingContribution;
        }

    }

    return res / static_cast<float>(numSamples);
}

// TODO: Standard feature
// Given a single parralelogram light, compute its contribution towards an incident ray at an intersection point
// by integrating over the parralelogram, taking `numSamples` samples from the light source, and applying
// shading.
//
// Hint: you can sample the light by using `sampleParallelogramLight(state.sampler.next_1d(), ...);`, which
//       you should implement first.
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the sample is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;      the active scene, feature config, bvh, and a thread-safe sampler
// - light;      the ParallelogramLight object, see `common.h`
// - ray;        the incident ray to the current intersection
// - hitInfo;    information about the current intersection
// - numSamples; the number of samples you need to take
// - return;     accumulated light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionParallelogramLight(RenderState& state, const ParallelogramLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    // TODO: implement this function; repeat numSamples times:
    // - sample the parallellogram light
    // - test the sample's visibility
    // - then evaluate the phong model
    glm::vec3 res(0.0f); // Initialize the accumulated light contribution

    for (uint32_t i = 0; i < numSamples; ++i) {
        glm::vec3 lightPosition, lightColor;
        /*float sampleU = state.sampler.next_1d();
        float sampleV = state.sampler.next_1d();
        glm::vec2 smp = glm::vec2(sampleU, sampleV);*/
        sampleParallelogramLight(state.sampler.next_2d(), light, lightPosition, lightColor);

        drawSphere(lightPosition, 0.005f, glm::vec3(1.0f,0.0f,0.0f));

        

        glm::vec3 col = visibilityOfLightSample(state, lightPosition, lightColor, ray, hitInfo);
        if (glm::length(col) > 0.0f) {

            glm::vec3 p = ray.origin + ray.t * ray.direction;
            glm::vec3 l = glm::normalize(lightPosition - p);
            glm::vec3 v = -ray.direction;
            glm::vec3 shadingContribution = computeShading(state, v, l, lightColor, hitInfo);
            res += shadingContribution;
        }
    }

   

    // Average the accumulated light contribution over the number of samples
        return res / static_cast<float>(numSamples);

}

// This function is provided as-is. You do not have to implement it.
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// This forowards to `visibilityOfLightSampleBinary`/`visibilityOfLightSampleTransparency` based on settings.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSample(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return lightColor;
    } else if (!state.features.enableTransparency) {
        // Shadows are enabled but transparency is disabled
        return visibilityOfLightSampleBinary(state, lightPosition, lightColor, ray, hitInfo) ? lightColor : glm::vec3(0);
    } else {
        // Shadows and transparency are enabled
        return visibilityOfLightSampleTransparency(state, lightPosition, lightColor, ray, hitInfo);
    }
}

// This function is provided as-is. You do not have to implement it.
glm::vec3 computeLightContribution(RenderState& state, const Ray& ray, const HitInfo& hitInfo)
{
    // Iterate over all lights
    glm::vec3 Lo { 0.0f };
    for (const auto& light : state.scene.lights) {
        if (std::holds_alternative<PointLight>(light)) {
            Lo += computeContributionPointLight(state, std::get<PointLight>(light), ray, hitInfo);
        } else if (std::holds_alternative<SegmentLight>(light)) {
            Lo += computeContributionSegmentLight(state, std::get<SegmentLight>(light), ray, hitInfo, state.features.numShadowSamples);
        } else if (std::holds_alternative<ParallelogramLight>(light)) {
            Lo += computeContributionParallelogramLight(state, std::get<ParallelogramLight>(light), ray, hitInfo, state.features.numShadowSamples);
        }
    }
    return Lo;
}