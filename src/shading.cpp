#include "render.h"
#include "texture.h"
#include <algorithm>
#include <cmath>
#include <fmt/core.h>
#include <glm/geometric.hpp>
#include <glm/gtx/string_cast.hpp>
#include <shading.h>
#include <ranges>

// This function is provided as-is. You do not have to implement it (unless
// you need to for some extra feature).
// Given render state and an intersection, based on render settings, sample
// the underlying material data in the expected manner.
glm::vec3 sampleMaterialKd(RenderState& state, const HitInfo& hitInfo)
{
    if (state.features.enableTextureMapping && hitInfo.material.kdTexture) {
        if (state.features.enableBilinearTextureFiltering) {
            return sampleTextureBilinear(*hitInfo.material.kdTexture, hitInfo.texCoord);
        } else {
            return sampleTextureNearest(*hitInfo.material.kdTexture, hitInfo.texCoord);
        }
    } else {
        return hitInfo.material.kd;
    }
}

// This function is provided as-is. You do not have to implement it.
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the scene-selected shading model, returning the reflected light towards the target.
glm::vec3 computeShading(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    // Hardcoded color gradient. Feel free to modify this
    static LinearGradient gradient = {
        .components = {
            { -1.0f, glm::vec3(1.0f, 0.0f, 0.0f) },
            { -0.7f, glm::vec3(1.0f, 0.5f, 0.0f) },
            { -0.4f, glm::vec3(1.0f, 1.0f, 0.0f) },
            { 0.0f, glm::vec3(0.0f, 1.0f, 0.0f) },
            { 0.4f, glm::vec3(0.0f, 0.0f, 1.0f) },
            { 0.7f, glm::vec3(0.29f, 0.0f, 0.51f) },
            { 1.0f, glm::vec3(0.56f, 0.0f, 1.0f) }

        }
    };

    if (state.features.enableShading) {
        const glm::vec3 kd = sampleMaterialKd(state, hitInfo);
        switch (state.features.shadingModel) {
        case ShadingModel::Lambertian:
            return computeLambertianModel(state, cameraDirection, lightDirection, lightColor, hitInfo, kd);
        case ShadingModel::Phong:
            return computePhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo, kd);
        case ShadingModel::BlinnPhong:
            return computeBlinnPhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo, kd);
        case ShadingModel::LinearGradient:
            return computeLinearGradientModel(state, cameraDirection, lightDirection, lightColor, hitInfo, gradient);
        case ShadingModel::LinearGradientComparison:
            return computeLinearGradientModelComparison(state, cameraDirection, lightDirection, lightColor, hitInfo, gradient);
        };
    }

    return lightColor * sampleMaterialKd(state, hitInfo);
}

// TODO: Standard feature
// Given a number ti between [-1, 1], sample from the gradient's components and return the
// linearly interpolated color, for which ti lies in the interval between the t-values of two
// components, or on a boundary. If ti falls outside the gradient's smallest/largest components,
// the nearest component must be sampled.
// - ti; a number between [-1, 1]
// This method is unit-tested, so do not change the function signature.
glm::vec3 LinearGradient::sample(float ti) const
{
    if (components.empty())
        return { 0.f, 0.f, 0.f };

    ti = glm::clamp(ti, -1.0f, 1.0f);
    std::vector<Component> sorted = components;
    std::sort(sorted.begin(), sorted.end(),
        [](const Component& a, const Component& b) { return a.t < b.t; });

    for (size_t i = 0; i < sorted.size() - 1; ++i) {
        const auto& c1 = sorted[i];
        const auto& c2 = sorted[i + 1];

        if (ti >= c1.t && ti <= c2.t) {
            float alpha = (ti - c1.t) / (c2.t - c1.t);
            return glm::mix(c1.color, c2.color, alpha);
        }
    }

    if (ti < sorted.front().t) {
        return sorted.front().color;
    } else {
        return sorted.back().color;
    }
}

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate a diffuse shading model, such that the diffuse component is sampled not
// from the intersected material, but a provided linear gradient, based on the cosine of theta
// as defined in the diffuse shading part of the Phong model.
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - gradient;        the linear gradient object
// - return;          the result of shading
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeLinearGradientModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo, const LinearGradient& gradient)
{
    float cos_theta = glm::dot(lightDirection, hitInfo.normal);
    cos_theta = glm::clamp(cos_theta, -1.0f, 1.0f);
    /*if (cos_theta <= 1e-6f && !state.features.enableTransparency)
        return { 0.f, 0.f, 0.f };*/

    glm::vec3 gradientColor = gradient.sample(cos_theta);

    return lightColor * gradientColor * cos_theta;
}

glm::vec3 computeLinearGradientModelComparison(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo, const LinearGradient& gradient)
{
    float cos_theta = glm::dot(lightDirection, hitInfo.normal);
    glm::vec3 kd = gradient.sample(cos_theta);
    glm::vec3 phong = computePhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo, kd);
    glm::vec3 blinnPhong = computeBlinnPhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo, kd);
    float d = glm::length(phong - blinnPhong);
    static LinearGradient diffgrad = {
        .components = {
            { -1.0f, glm::vec3(0.0f, 0.0f, 1.0f) }, // blue for having big diff
            { 0.0f, glm::vec3(0.0f, 1.0f, 0.0f) }, // green for middle, signaling no difference between models
            { 1.0f, glm::vec3(1.0f, 0.0f, 0.0f) } } // red for phong having much more light than blinn

    };

    // cap the diff
    float colfac = d * 10;
    if (colfac > 1)
        colfac = 1;
    if (colfac < -1)
        colfac = -1;

    return diffgrad.sample(colfac); // we multiply by 10 to see the diff more pronounced
}
