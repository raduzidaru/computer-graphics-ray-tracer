#include "bvh.h"

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
