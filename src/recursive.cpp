#include "recursive.h"
#include "draw.h"
#include "bvh_interface.h"
#include "intersect.h"
#include "extra.h"
#include "light.h"

// This function is provided as-is. You do not have to implement it.
// Given a range of rays, render out all rays and average the result
glm::vec3 renderRays(RenderState& state, std::span<const Ray> rays, int rayDepth)
{
    glm::vec3 L { 0.f };
    for (const auto& ray : rays) {
        L += renderRay(state, ray, rayDepth);
    }
    return L / static_cast<float>(rays.size());
}

// This method is provided as-is. You do not have to implement it.
// Given a camera ray (or secondary ray), tests for a scene intersection, and
// dependent on the results, evaluates the following functions which you must
// implement yourself:
// - `computeLightContribution()` and its submethods
// - `renderRaySpecularComponent()`, `renderRayTransparentComponent()`, `renderRayGlossyComponent()`
glm::vec3 renderRay(RenderState& state, Ray ray, int rayDepth)
{
    // Trace the ray into the scene. If nothing was hit, return early
    HitInfo hitInfo;
    if (!state.bvh.intersect(state, ray, hitInfo)) {
        if (state.features.enableDebugDraw) {
            drawRay(ray, glm::vec3(1, 0, 0));
        }
        return sampleEnvironmentMap(state, ray);
    }

    // Return value: the light along the ray
    // Given an intersection, estimate the contribution of scene lights at this intersection
    glm::vec3 Lo = computeLightContribution(state, ray, hitInfo);

    // Draw an example debug ray for the incident ray (feel free to modify this for yourself)
    drawRay(ray, glm::vec3(1.0f));

    // Given that recursive components are enabled, and we have not exceeded maximum depth,
    // estimate the contribution along these components
    if (rayDepth < 6) {
        bool isReflective = glm::any(glm::notEqual(hitInfo.material.ks, glm::vec3(0.0f)));
        bool isTransparent = hitInfo.material.transparency != 1.f;

        // Default, specular reflections
        if (state.features.enableReflections && !state.features.extra.enableGlossyReflection && isReflective) {
            renderRaySpecularComponent(state, ray, hitInfo, Lo, rayDepth);
        }

        // Alternative, glossy reflections
        if (state.features.enableReflections && state.features.extra.enableGlossyReflection && isReflective) {
            renderRayGlossyComponent(state, ray, hitInfo, Lo, rayDepth);
        }

        // Transparency passthrough
        if (state.features.enableTransparency && isTransparent) {
            renderRayTransparentComponent(state, ray, hitInfo, Lo, rayDepth);
        }
    }

    return Lo;
}

// TODO: Standard feature
// Given an incident ray and a intersection point, generate a mirrored ray
// - Ray;     the indicent ray
// - HitInfo; hit struct for the intersection point
// - return;  a reflected ray
// This method is unit-tested, so do not change the function signature.
Ray generateReflectionRay(Ray ray, HitInfo hitInfo)
{
    // TODO: generate a mirrored ray
    //       if you use glm::reflect, you will not get points for this method!
    Ray result;
    glm::vec3 intersection = ray.origin + ray.direction * ray.t;
    glm::vec3 reflected = ray.direction - 2.0f * glm::dot(ray.direction, hitInfo.normal) * hitInfo.normal;
    result.direction = reflected;
    result.origin = intersection - 0.0001f * ray.direction;
    return result;
}

// TODO: Standard feature
// Given an incident ray and a intersection point, generate a passthrough ray for transparency,
// starting at the intersection point and continuing in the same direction.
// - Ray;     the indicent ray
// - HitInfo; hit struct for the intersection point
// - return;  a passthrough ray for transparency
// This method is unit-tested, so do not change the function signature.
Ray generatePassthroughRay(Ray ray, HitInfo hitInfo)
{
    // TODO: generate a passthrough ray
    Ray result;
    glm::vec3 intersection = ray.origin + ray.direction * ray.t;
    result.origin = intersection + 0.0001f * ray.direction;
    result.direction = ray.direction;
    return result;
}

// TODO: standard feature
// Given a camera ray (or secondary ray) and an intersection, evaluates the contribution
// of a mirrored ray, recursively evaluating renderRay(..., depth + 1) along this ray,
// and adding the result times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is unit-tested, so do not change the function signature.
void renderRaySpecularComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // TODO; you should first implement generateReflectionRay()
    Ray r = generateReflectionRay(ray, hitInfo);
    hitColor += renderRay(state, r, rayDepth + 1) * hitInfo.material.ks;
    if (state.features.enableDebugDraw) {
        Ray n;
        n.origin = ray.origin + ray.direction * ray.t;
        n.t = 1.0f;
        n.direction = hitInfo.normal;
        drawRay(n, glm::vec3(1.0f, 1.0f, 0.0f));

        drawRay(r, glm::vec3(0.0f, 1.0f, 0.0f));

        for (int i = 0; i < state.scene.lights.size(); i++) {
            Ray light;
            light.origin = ray.origin + ray.direction * ray.t;
            if (std::holds_alternative<PointLight>(state.scene.lights[i])) {
                light.direction = std::get<PointLight>(state.scene.lights[i]).position - light.origin;
                light.t = glm::length(light.direction);
                light.direction = glm::normalize(light.direction);
                drawRay(light, glm::vec3(1.0f, 0.0f, 0.0f));
            } else if (std::holds_alternative<SegmentLight>(state.scene.lights[i])) {
                light.direction = std::get<SegmentLight>(state.scene.lights[i]).endpoint0 - light.origin;
                light.t = glm::length(light.direction);
                light.direction = glm::normalize(light.direction);
                drawRay(light, glm::vec3(1.0f, 0.0f, 0.0f));

                Ray light2;
                light2.origin = light.origin;
                light2.direction = std::get<SegmentLight>(state.scene.lights[i]).endpoint1 - light2.origin;
                light2.t = glm::length(light2.direction);
                light2.direction = glm::normalize(light2.direction);
                drawRay(light2, glm::vec3(1.0f, 0.0f, 0.0f));
            } else if (std::holds_alternative<ParallelogramLight>(state.scene.lights[i])) {
                light.direction = std::get<ParallelogramLight>(state.scene.lights[i]).v0 - light.origin;
                light.t = glm::length(light.direction);
                light.direction = glm::normalize(light.direction);
                drawRay(light, glm::vec3(1.0f, 0.0f, 0.0f));
            }
        }
    }
}

// TODO: standard feature
// Given a camera ray (or secondary ray) and an intersection, evaluates the contribution
// of a passthrough transparent ray, recursively evaluating renderRay(..., depth + 1) along this ray,
// and correctly alpha blending the result with the current intersection's hit color
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is unit-tested, so do not change the function signature.
void renderRayTransparentComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // TODO; you should first implement generatePassthroughRay()
    Ray r = generatePassthroughRay(ray, hitInfo);
    hitColor = renderRay(state, r, rayDepth + 1) * hitInfo.material.transparency + (1.0f - hitInfo.material.transparency) * hitColor;
    if (state.features.enableDebugDraw) {
        drawRay(r, hitColor);
        drawSphere(ray.origin + ray.direction * ray.t, 0.01f, hitInfo.material.transparency * glm::vec3(hitInfo.material.kd.r, hitInfo.material.kd.g, hitInfo.material.kd.b));
    }
}