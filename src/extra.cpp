#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include <framework/trackball.h>
#include "draw.h"
#include "texture.h"
#include "interpolate.h"

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of Depth of Field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithDepthOfField(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableDepthOfField) {
        return;
    }

    // ...
}

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of motion blur. Here, you integrate over a time domain, and not just the pixel's image domain,
// to give objects the appearance of "fast movement".
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithMotionBlur(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableMotionBlur) {
        return;
    }

}

// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
    }

    // ...
}


// TODO; Extra feature
// Given a camera ray (or reflected camera ray) and an intersection, evaluates the contribution of a set of
// glossy reflective rays, recursively evaluating renderRay(..., depth + 1) along each ray, and adding the
// results times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderRayGlossyComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    auto numSamples = state.features.extra.numGlossySamples;
    Ray reflectedRay = generateReflectionRay(ray, hitInfo);
    glm::vec3 basis1 = glm::normalize(glm::cross(hitInfo.normal, reflectedRay.direction));
    glm::vec3 basis2 = glm::normalize(glm::cross(basis1, reflectedRay.direction));
    glm::vec3 center = reflectedRay.origin;

    float diskRadius = 1.0f / (hitInfo.material.shininess + 1.0f);

    for (int i = 0; i < numSamples; ++i) {
        glm::vec2 sample = state.sampler.next_2d();

        float r = std::sqrt(sample.x);
        float theta = 2.0f * glm::pi<float>() * sample.y;
        float x = r * std::cos(theta);
        float y = r * std::sin(theta);

        glm::vec3 offset = diskRadius * (x * basis1 + y * basis2);
        Ray secondRay;
        secondRay.origin = center;
        secondRay.direction = reflectedRay.direction + offset;
        hitColor += renderRay(state, secondRay, rayDepth + 1) * hitInfo.material.ks / static_cast<float>(numSamples);
    }
}

bool isPointInTriangle(const glm::vec3& point, const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    glm::vec3 e0 = v1 - v0;
    glm::vec3 e1 = v2 - v1;
    glm::vec3 e2 = v0 - v2;
    glm::vec3 C0 = point - v0;
    glm::vec3 C1 = point - v1;
    glm::vec3 C2 = point - v2;
    if (glm::dot(glm::cross(e0, C0), glm::cross(e0, e1)) < 0)
        return false;
    if (glm::dot(glm::cross(e1, C1), glm::cross(e1, e2)) < 0)
        return false;
    if (glm::dot(glm::cross(e2, C2), glm::cross(e2, e0)) < 0)
        return false;

    return true;
}

glm::vec3 getTextureColor(const RenderState& state, const Ray& ray)
{
    Mesh envMapMesh = state.scene.envMapCube[0];
    float closestT = std::numeric_limits<float>::max();
    glm::vec2 closestTexCoord;
    bool hit = false;

    for (const auto& triangle : envMapMesh.triangles) {
        const glm::vec3& vertex0 = envMapMesh.vertices[triangle.x].position;
        const glm::vec3& vertex1 = envMapMesh.vertices[triangle.y].position;
        const glm::vec3& vertex2 = envMapMesh.vertices[triangle.z].position;
        const glm::vec2& tex0 = envMapMesh.vertices[triangle.x].texCoord;
        const glm::vec2& tex1 = envMapMesh.vertices[triangle.y].texCoord;
        const glm::vec2& tex2 = envMapMesh.vertices[triangle.z].texCoord;
        glm::vec3 normal = glm::normalize(glm::cross(vertex1 - vertex0, vertex2 - vertex0));
        float denom = glm::dot(ray.direction, normal);
        if (std::abs(denom) < 1e-6f)
            continue;
        float t = glm::dot(vertex0 - ray.origin, normal) / denom;
        if (t < 1e-6f || t >= closestT)
            continue;
        glm::vec3 intersection = ray.origin + t * ray.direction;
        if (!isPointInTriangle(intersection, vertex0, vertex1, vertex2))
            continue;
        closestT = t;
        hit = true;
        glm::vec3 barycentricCoords = computeBarycentricCoord(vertex0, vertex1, vertex2, intersection);
        closestTexCoord = barycentricCoords.x * tex0 + barycentricCoords.y * tex1 + barycentricCoords.z * tex2;
    }

    if (hit && state.scene.envMap != nullptr) {
        closestTexCoord = glm::clamp(closestTexCoord, glm::vec2(0.0f), glm::vec2(1.0f));
        if (state.features.enableBilinearTextureFiltering) {
            return sampleTextureBilinear(*state.scene.envMap, closestTexCoord);
        }
        return sampleTextureNearest(*state.scene.envMap, closestTexCoord);
    }

    return glm::vec3(0.0f);
}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You will have to add support for environment textures
// to the Scene object, and provide a scene with the right data to supply this.
// - state; the active scene, feature config, bvh, and sampler
// - ray;   ray object
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
glm::vec3 sampleEnvironmentMap(RenderState& state, Ray ray)
{
    if (state.features.extra.enableEnvironmentMap && state.scene.envMap != nullptr) {
        return getTextureColor(state, ray);
    } else {
        return glm::vec3(0.f);
    }
}
