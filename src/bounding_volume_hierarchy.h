#pragma once
#include "ray_tracing.h"
#include "scene.h"

class BoundingVolumeHierarchy {
public:
    BoundingVolumeHierarchy(Scene* pScene);

    // Use this function to visualize your BVH. This can be useful for debugging.
    void debugDraw(int level);
    int numLevels() const;
    void splitBVH();
    void SetFirstTwoBox(Scene* pScene);
    void SetFirstBox(Scene* pScene);
    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo) const;

    struct Node {
        //interior means its an parent
        bool interior = false; //interior or leaf
        std::vector<float> indices; //indices of children OR triangles
        //std::vector<float> triIndices;
        AxisAlignedBox aabb;
        int level = 0;
        int meshNumber = 0;
    };

    Node getAabb(Node root, Ray& ray) const;

    void SetThem(Node &hello);
    std::vector<Node> nodes; //stores all nodes
    float testIntersect(Node test, Ray& ray) const;
    bool finalIntersect(AxisAlignedBox box, Ray& ray)const ;

private:
    Scene* m_pScene;
};
