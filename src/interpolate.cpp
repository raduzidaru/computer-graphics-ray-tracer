#include "interpolate.h"
#include <glm/geometric.hpp>
#include <iostream>

// TODO Standard feature
// Given three triangle vertices and a point on the triangle, compute the corresponding barycentric coordinates of the point.
// and return a vec3 with the barycentric coordinates (alpha, beta, gamma).
// - v0;     Triangle vertex 0
// - v1;     Triangle vertex 1
// - v2;     Triangle vertex 2
// - p;      Point on triangle
// - return; Corresponding barycentric coordinates for point p.
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeBarycentricCoord(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    glm::vec3 v0v1 = v1 - v0;
    glm::vec3 v0v2 = v2 - v0;
    glm::vec3 v0p = p - v0;

    float d00 = glm::dot(v0v2, v0v2);
    float d01 = glm::dot(v0v2, v0v1);
    float d11 = glm::dot(v0v1, v0v1);
    float d20 = glm::dot(v0p, v0v2);
    float d21 = glm::dot(v0p, v0v1);

    float denom = d00 * d11 - d01 * d01;

    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;

    if (u >= 0.0f && v >= 0.0f && w >= 0.0f && u <= 1.0f && v <= 1.0f && w <= 1.0f) {
        return glm::vec3(u, w, v);
    }

    glm::vec3 closestPoint;

    glm::vec3 edge0Closest = v0 + glm::clamp(glm::dot(v0p, v0v1) / glm::dot(v0v1, v0v1), 0.0f, 1.0f) * v0v1;
    float dist0 = glm::length(p - edge0Closest);

    glm::vec3 edge1Closest = v1 + glm::clamp(glm::dot(p - v1, v2 - v1) / glm::dot(v2 - v1, v2 - v1), 0.0f, 1.0f) * (v2 - v1);
    float dist1 = glm::length(p - edge1Closest);

    glm::vec3 edge2Closest = v2 + glm::clamp(glm::dot(p - v2, v0 - v2) / glm::dot(v0 - v2, v0 - v2), 0.0f, 1.0f) * (v0 - v2);
    float dist2 = glm::length(p - edge2Closest);

    if (dist0 < dist1 && dist0 < dist2) {
        closestPoint = edge0Closest;
    } else if (dist1 < dist2) {
        closestPoint = edge1Closest;
    } else {
        closestPoint = edge2Closest;
    }

    glm::vec3 closestVec = closestPoint - v0;
    d20 = glm::dot(closestVec, v0v2);
    d21 = glm::dot(closestVec, v0v1);

    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.0f - v - w;

    return glm::vec3(u, w, v);
}

// TODO Standard feature
// Linearly interpolate three normals using barycentric coordinates.
// - n0;     Triangle normal 0
// - n1;     Triangle normal 1
// - n2;     Triangle normal 2
// - bc;     Barycentric coordinate
// - return; The smoothly interpolated normal.
// This method is unit-tested, so do not change the function signature.
glm::vec3 interpolateNormal(const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 bc)
{
    glm::vec3 interpolatedNormal = bc.x * n0 + bc.y * n1 + bc.z * n2;
    return glm::normalize(interpolatedNormal);
}

// TODO Standard feature
// Linearly interpolate three texture coordinates using barycentric coordinates.
// - n0;     Triangle texture coordinate 0
// - n1;     Triangle texture coordinate 1
// - n2;     Triangle texture coordinate 2
// - bc;     Barycentric coordinate
// - return; The smoothly interpolated texturre coordinate.
// This method is unit-tested, so do not change the function signature.
glm::vec2 interpolateTexCoord(const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 bc)
{
    glm::vec2 interpolatedTexCoord = bc.x * t0 + bc.y * t1 + bc.z * t2;
    return interpolatedTexCoord;
}
