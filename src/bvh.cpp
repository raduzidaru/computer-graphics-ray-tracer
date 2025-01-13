#include "bvh.h"
#include <stack>
#include "scene.h"
#include "render.h"
#include "intersect.h"
#include "interpolate.h"
#include "draw.h"

bool intersectRayWithAABB(const Ray& ray, const AxisAlignedBox& aabb)
{
    float t_min = 0.0f;
    float t_max = ray.t;
    for (int axis = 0; axis < 3; ++axis) {
        float lower = aabb.lower[axis];
        float upper = aabb.upper[axis];
        float origin = ray.origin[axis];
        float direction = ray.direction[axis];
        if (direction == 0.0f) {
            if (origin < lower || origin > upper) {
                return false;
            }
            continue;
        }
        float t1 = (lower - origin) / direction;
        float t2 = (upper - origin) / direction;
        if (t1 > t2) {
            std::swap(t1, t2);
        }
        t_min = std::max(t_min, t1);
        t_max = std::min(t_max, t2);
        if (t_min > t_max) {
            return false;
        }
    }
    return true;
}

glm::vec3 intersectionPointRayWithAABB(const Ray& ray, const AxisAlignedBox& aabb)
{
    float t_min = 0.0f;
    float t_max = ray.t;
    glm::vec3 intersectionPoint = glm::vec3(std::numeric_limits<float>::max());
    for (int axis = 0; axis < 3; ++axis) {
        float lower = aabb.lower[axis];
        float upper = aabb.upper[axis];
        float origin = ray.origin[axis];
        float direction = ray.direction[axis];
        if (direction == 0.0f) {
            if (origin < lower || origin > upper) {
                return glm::vec3(std::numeric_limits<float>::max());
            }
            continue;
        }
        float t1 = (lower - origin) / direction;
        float t2 = (upper - origin) / direction;
        if (t1 > t2) {
            std::swap(t1, t2);
        }
        t_min = std::max(t_min, t1);
        t_max = std::min(t_max, t2);
        if (t_min > t_max) {
            return glm::vec3(std::numeric_limits<float>::max());
        }
    }
    intersectionPoint = ray.origin + ray.direction * t_min;
    return intersectionPoint;
}


void updateHitInfo(RenderState& state, const Scene& scene, const BVHInterface::Primitive& primitive, const Ray& ray, HitInfo& hitInfo)
{
    const auto& [v0, v1, v2] = std::tie(primitive.v0, primitive.v1, primitive.v2);
    const auto& mesh = scene.meshes[primitive.meshID];
    glm::vec3 intersect = ray.origin + ray.direction * ray.t;

    hitInfo.material = mesh.material;
    if (state.features.enableNormalInterp) {
        hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, intersect);
        hitInfo.normal = interpolateNormal(primitive.v0.normal, primitive.v1.normal, primitive.v2.normal, hitInfo.barycentricCoord);
        hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);
    } else {
        const auto n = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));
        hitInfo.normal = n;
    }

    // Catch flipped normals
    if (glm::dot(ray.direction, hitInfo.normal) > 0.0f) {
        hitInfo.normal = -hitInfo.normal;
    }
}

int countIntersectedAABBs(const BVHInterface& bvh, const Ray& ray)
{
    std::span<const BVHInterface::Node> nodes = bvh.nodes();
    int intersectedCount = 0;

    std::stack<uint32_t> stack;
    stack.push(0);

    while (!stack.empty()) {
        uint32_t currentNodeIdx = stack.top();
        stack.pop();

        const auto& currentNode = nodes[currentNodeIdx];

        if (intersectRayWithAABB(ray, currentNode.aabb)) {
            intersectedCount++;

            if (!currentNode.isLeaf()) {
                stack.push(currentNode.leftChild());
                stack.push(currentNode.rightChild());
            }
        }
    }

    return intersectedCount;
}


bool intersectRayWithBVHImplemented(RenderState& state, const BVHInterface& bvh, Ray& ray, HitInfo& hitInfo)
{
    std::span<const BVHInterface::Node> nodes = bvh.nodes();
    std::span<const BVHInterface::Primitive> primitives = bvh.primitives();

    bool is_hit = false;
    std::stack<std::pair<uint32_t, int>> stack;
    stack.push({ 0, 0 });
    std::vector<int> intersectedNodes;


    while (!stack.empty()) {
        auto [currentNodeIdx, level] = stack.top();
        stack.pop();

        BVHInterface::Node currentNode = nodes[currentNodeIdx];

        if (!intersectRayWithAABB(ray, currentNode.aabb)) {
            continue;
        }
        intersectedNodes.push_back(currentNodeIdx);

        if (state.features.enableDebugDraw) {
            glm::vec3 intersectionPoint = intersectionPointRayWithAABB(ray, currentNode.aabb);
            if (intersectionPoint.x != std::numeric_limits<float>::max()) {
                float normalizedOrder = float(intersectedNodes.size()) / float(countIntersectedAABBs(bvh, ray));
                glm::vec3 gradientColor = glm::mix(glm::vec3(1.0f, 1.0f, 0.0f), glm::vec3(1.0f, 0.0f, 0.0f), normalizedOrder);
                drawSphere(intersectionPoint, 0.005f, gradientColor);
                drawAABB(currentNode.aabb, DrawMode::Wireframe, gradientColor);
            }
        }

        if (currentNode.isLeaf()) {
            size_t offset = currentNode.primitiveOffset();
            size_t count = currentNode.primitiveCount();

            for (size_t i = 0; i < count; ++i) {
                auto& prim = primitives[offset + i];
                const auto& [v0, v1, v2] = std::tie(prim.v0, prim.v1, prim.v2);
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                    is_hit = true;
                    updateHitInfo(state, state.scene, prim, ray, hitInfo);
                    if (state.features.enableDebugDraw) {
                        drawTriangle(v0, v1, v2);
                    }
                }
            }
        } else {
            stack.push({ currentNode.rightChild(), level + 1 });
            stack.push({ currentNode.leftChild(), level + 1 });
        }
    }

    return is_hit;
}



// TODO Standard Feature
// Hierarchy traversal routine; you must implement this method and implement it carefully!
//
// The default implementation uses precompiled methods and only performs rudimentary updates to hitInfo.
// For correct normal interpolation, barycentric coordinates, etc., you have to implement this traversal yourself.
// If you have implemented the traversal method for assignment 4B, you might be able to reuse most of it.
//
// This method returns `true` if geometry was hit, and `false` otherwise. On first/closest hit, the
// distance `t` in the `ray` object is updated, and information is updated in the `hitInfo` object.
// - state;    current render state (containing scene, features, ...)
// - bvh;      the actual bvh which should be traversed for faster intersection
// - ray;      the ray intersecting the scene's geometry
// - hitInfo;  the return object, with info regarding the hit geometry
// - return;   boolean, if geometry was hit or not
bool BVH::intersect(RenderState& state, Ray& ray, HitInfo& hitInfo) const
{
    if (state.features.enableAccelStructure) {
        return intersectRayWithBVHImplemented(state, *this, ray, hitInfo);
    }
    return intersectRayWithBVH(state, *this, ray, hitInfo);

    // TODO: implement here your (probably stack-based) BVH traversal.
    // Some hints (refer to bvh_interface.h either way). BVH nodes are packed, so the
    // data is not easily extracted. Helper methods are available, however:
    // - For a given node, you can test if the node is a leaf with `node.isLeaf()`.
    // - If the node is not a leaf, you can obtain the left/right children with `node.leftChild()` etc.
    // - If the node is a leaf, you can obtain the offset to and nr. of primitives in the bvh's list
    //   of underlying primitives with `node.primitiveOffset()` and `node.primitiveCount()`
    //
    // In short, you will have to step down the bvh, node by node, and intersect your ray
    // with the node's AABB. If this intersection passes, you should:
    // - if the node is a leaf, intersect with the leaf's primitives
    // - if the node is not a leaf, test the left and right children as well!
    //
    // Note that it is entirely possible for a ray to hit a leaf node, but not its primitives,
    // and it is likewise possible for a ray to hit both children of a node.
    //
    // Make sure to update the hitInfo in case of a hit
}
