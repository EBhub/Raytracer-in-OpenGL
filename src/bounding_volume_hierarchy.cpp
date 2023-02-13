#include <glm/geometric.hpp>
#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

float minX = 0;
float minY = 0;
float minZ = 0;
float maxX = 0;
float maxY = 0;
float maxZ = 0;

int h = 0;

int first = 0;

int totalLevels = 0;

int axisDevide = 0;

std::vector<Mesh> scene;
Scene* currentScene;

std::vector<Triangle> allTriangles;
std::vector<Vertex> allVertices;
std::vector<int> intersectedNodes;
std::vector<BoundingVolumeHierarchy::Node> allLeafs;

void SetMax(float x, float y, float z) {
    if (x > maxX) {
        maxX = x;
    }
    if (y > maxY) {
        maxY= y;
    }
    if (z > maxZ) {
        maxZ = z;
    }
}

void SetMin(float x, float y, float z) {
    if (x < minX) {
        minX = x;
    }
    if (y < minY) {
        minY = y;
    }
    if (z < minZ) {
        minZ = z;
    }
}

void ResetValues() {
    maxX = 0.0f;
    maxY = 0.0f;
    maxZ = 0.0f;

    minX = 0.0f;
    minY = 0.0f;
    minZ = 0.0f;
}

void BoundingVolumeHierarchy::SetThem(Node &hello) {

    float meanX = 0.0f;
    float meanY = 0.0f;
    float meanZ = 0.0f;
    float size = 1.0f;

    bool first = true;

    //method that loops over all triangles of node and sets the minMaxVals accordingly
    //for (const auto& mesh : m_pScene->meshes) {
        //b for the triangle index to check 
        //for every index in nodes indices (the triangle indexes)
        for (int c = 0; c < hello.indices.size(); c++) {
            //the triangle index <- this triangle has the indexes to the points in space
            if (hello.indices[c] < allTriangles.size()) {
                int triIndex = hello.indices[c];
                //TRIANGLE IS FOUND!!

                /*
                * loop over the 3 verices of the triangle
                * and set the min and maxes
                * and set the MEAN!
                */
                for (int x = 0; x < 3; x++) {
                    const auto& test = allVertices[allTriangles[triIndex][x]];
                    //if first loop, set the min/max/avg to the first vertex. Use this as reference now
                    if (first) {
                        minX = test.p.x;
                        minY = test.p.y;
                        minZ = test.p.z;
                        maxX = test.p.x;
                        maxY = test.p.y;
                        maxZ = test.p.z;
                        first = false;
                    }
                    //std::cout << "-------    X: " << test.p.x << std::endl;
                    meanX += test.p.x;
                    meanY += test.p.y;
                    meanZ += test.p.z;
                    size++;
                    SetMin(test.p.x, test.p.y, test.p.z);
                    SetMax(test.p.x, test.p.y, test.p.z);
                }
                
            }
        }
    //}

    //set the min/max/avg values
        /*
    std::cout << minX << std::endl;
    std::cout << minY << std::endl;
    std::cout << minZ << std::endl;
    std::cout << maxX << std::endl;
    std::cout << maxY << std::endl;
    std::cout << maxZ << std::endl;
    std::cout << std::endl;
    */

    hello.aabb.lower.x = minX;
    hello.aabb.lower.y = minY;
    hello.aabb.lower.z = minZ;
    hello.aabb.upper.x = maxX;
    hello.aabb.upper.y = maxY;
    hello.aabb.upper.z = maxZ;
    ResetValues();
}

/* Here the very first AABB gets made.
* This method goes over the vector with all the triangles and finds the 
* minima and maxima.
 */
void BoundingVolumeHierarchy::SetFirstBox(Scene* pScene)
{
    Node root;
    // 1. Get all the maxes and minimums of the mesh
    //for (const auto& mesh : pScene->meshes) {
        //for (const auto& tri : mesh.triangles) {
        for (const auto& tri : allTriangles) {
            for (int x = 0; x < 3; x++) {
                //loop over the 3 verices
                const auto test = allVertices[tri[x]];

                SetMax(test.p.x, test.p.y, test.p.z);
                SetMin(test.p.x, test.p.y, test.p.z);
            }
        }
    //}

    int triangleIndex = 0;

    //add all the triangles in the scene to this first box
    //for (const auto& mesh : pScene->meshes) {
        //for every triangle in this mesh do this
        for (const auto& tri : allTriangles) {
            root.indices.push_back(triangleIndex);
            triangleIndex++;
        }
    //}

    SetThem(root);
    nodes.push_back(root);

}

void BoundingVolumeHierarchy::SetFirstTwoBox(Scene* pScene) {
    nodes[0].interior = true;
    //check when avg = 0?
    //float averageX = (maxX + minX) / 2;

    //Node root;
    Node a1;
    Node b1;

    //own + 1 and 2
    //root.interior = true;
    //root.indices.push_back(1);
    //root.indices.push_back(2);


    //std::cout << "Average on FIRST SPLIT:" << averageX << std::endl << std::endl;

    //the index of the triangle IN THIS MESH
    int triangleIndex = 0;

    //now fill the nodes with all the triangle indices
    // If one point is contained in the box, INSERT THE VERTICES of that triangle
    // Here we construct the first 2 sub AABB's
    //for every mesh do this
    //for (const auto& mesh : pScene->meshes) {
        //for every triangle in this mesh do this
        for (const auto& tri : allTriangles) {
            //loop over the 3 verices of the triangle
            //for (int x = 0; x < 3; x++) {
                //const auto test = mesh.vertices[tri[0]];
                //test if these vertices are contained in the AABB, test on avgX
                //is maller, put in 1a
                //insert the trangles index
                float centroidX = (allVertices[tri[0]].p.x + allVertices[tri[1]].p.x + allVertices[tri[2]].p.x) / 3;
                //std::cout << "CENTROID: " << centroidX << std::endl;

                //if (test.p.x > 0) {
                if (centroidX < 0) { 
                    a1.indices.push_back(triangleIndex);
                    //std::cout << "AABB A1 -> Index: " << triangleIndex << std::endl;
                    
                } else 
                //if (test.p.x < 0) {
                if (centroidX > 0) {
                    b1.indices.push_back(triangleIndex);
                    //std::cout << "  AABB B1 -> Index: " << triangleIndex << std::endl;
                    
                }
            triangleIndex++;
        }
    //}

    //Set the min and max values
    SetThem(a1);
    SetThem(b1);


    //set the level of the second and third box
    // the level of the root box defaults to 0 (all boxes do)
    a1.level = 1;
    b1.level = 1;

    //nodes.push_back(root);
    nodes.push_back(a1);
    nodes.push_back(b1);
    nodes[0].indices.clear();
    nodes[0].indices.push_back(1);
    nodes[0].indices.push_back(2);
}

void BoundingVolumeHierarchy::splitBVH() {
    // 1. traverse the tree to the leaf nodes. 
    // if leaf, split it!

    //TODO: SPLIT UNTIL LEVEL OR SIZE OF NODE?

    Node left;
    Node right;
    int t{};

    float average = 0;

    //traverse
    int size = nodes.size();
    for (int a = 0; a < size; a++) {
        //if leaf do this: Thus FOR ALL LEAFS
        if (nodes[a].interior == false) {
            //how big the indices need to be to split
            if (nodes[a].indices.size() > 1) {
                //node at indice a is a leaf node. This leaf node has the traingle indices that are in its space

                //check what axis we devide now
                float average = 0;
                if (axisDevide == 0) {
                    average = (nodes[a].aabb.lower.x + nodes[a].aabb.upper.x) / 2;
                    //average = nodes[a].minMaxVals[6];
                    std::cout << "      We split on the X axis " << average << std::endl;
                }
                else if (axisDevide == 1) {
                    average = (nodes[a].aabb.lower.y + nodes[a].aabb.upper.y) / 2;
                    //average = nodes[a].minMaxVals[7];
                    std::cout << "      We split on the Y axis " << average << std::endl;
                }
                else if (axisDevide == 2) {
                    average = (nodes[a].aabb.lower.z + nodes[a].aabb.upper.z) / 2;
                    //average = nodes[a].minMaxVals[8];
                    std::cout << "      We split on the Z axis " << average << std::endl;
                }
                //now use this average to split

                //reset the vectors
                left.indices.clear();
                right.indices.clear();

                //for every mesh
                //for (const auto& mesh : m_pScene->meshes) {
                //for (const auto& tri : allTriangles) {
                    //b for the triangle index to check 
                    //for every index in nodes indices (the triangle indexes)

                for (int c = 0; c < nodes[a].indices.size(); c++) {
                    //the triangle index <- this triangle has the indexes to the points in space
                    //if (nodes[a].indices[c] < mesh.triangles.size()) {
                    if (nodes[a].indices[c] < allTriangles.size()) {
                        int triIndex = nodes[a].indices[c];
                        //TRIANGLE IS FOUND!!

                        //test if the middle of these triangles is on what side of the AABB
                        //is maller, put in left

                        float centroid = (
                            allVertices[allTriangles[triIndex][0]].p[axisDevide] +
                            allVertices[allTriangles[triIndex][1]].p[axisDevide] +
                            allVertices[allTriangles[triIndex][2]].p[axisDevide]) / 3;

                        if (centroid < average) {
                            left.indices.push_back(triIndex);
                            //std::cout << "AABB A1 -> Index: " << triangleIndex << std::endl;
                            //std::cout << "push back global triangle: " << triIndex << " for LEFT node on place: " << left.indices.size() << std::endl;
                        }
                        else if (centroid > average) {
                            right.indices.push_back(triIndex);
                            //std::cout << "AABB B1 -> Index: " << triangleIndex << std::endl;
                            //std::cout << "push back global triangle: " << triIndex << " for RIGHT node on place: "<< right.indices.size() << std::endl;
                        }

                    }
                }

                //std::cout << "NODES size before PUSH: " << nodes.size() << std::endl << std::endl;

                SetThem(left);
                SetThem(right);

                left.level = totalLevels + 2;
                right.level = totalLevels + 2;

                nodes.push_back(left);
                nodes.push_back(right);

                //clear indices of parent node and set the indices to the child AND set to iterior node
                nodes[a].interior = true;
                nodes[a].indices.clear();
                nodes[a].indices.push_back(nodes.size() - 1);
                nodes[a].indices.push_back(nodes.size() - 2);

                std::cout << "SPLITTED: left size: " << left.indices.size() << ", Right size: " << right.indices.size() << std::endl;
                //}
            }
            //size is not bigger then 1 (look at the if statement of this else)
            else {

                Node push;
                push.level = nodes[a].level + 1;

                for (int d = 0; d < nodes[a].indices.size(); d++) {
                    push.indices.push_back(nodes[a].indices[d]);
                }

                push.aabb.lower.x = nodes[a].aabb.lower.x;
                push.aabb.lower.y = nodes[a].aabb.lower.y;
                push.aabb.lower.z = nodes[a].aabb.lower.z;
                push.aabb.upper.x = nodes[a].aabb.upper.x;
                push.aabb.upper.y = nodes[a].aabb.upper.y;
                push.aabb.upper.z = nodes[a].aabb.upper.z;
                push.interior = false;
                nodes.push_back(push);
            }
            t = nodes[a].indices.size();
        }
    }

    axisDevide = (axisDevide + 1) % 3;

    if (t > 1) {
        totalLevels++;
        splitBVH();
    }
}

/* Method to create one big vector with all the triangles from
* all mmeshes in the scene. 
*/
void createVector(Scene* pScene) {
    allTriangles.clear();
    allVertices.clear();
    
    std::cout << "CALLED the vector creation" << std::endl << std::endl;

    int totalMeshes = 0;
    int totalTriangles = 0;
    int vectorOffset = 0;
    int indicesPerTriangle = 0;
    Triangle input;
    int j = 0;
    // insert the triangle indices of the mesh
    // For every mesh in the scene 
    for (auto& mesh : pScene->meshes) {
        totalMeshes++;
        vectorOffset = allTriangles.size();
        /*
        std::cout << "MESH FOUND" << std::endl;
        std::cout << std::endl << "TRIANGLE OFFSET: " << vectorOffset << std::endl;
        std::cout << "MESH INDICES SIZE: " << mesh.vertices.size() << std::endl;
        std::cout << "MESH TRIANGLES SIZE: " << mesh.triangles.size() << std::endl << std::endl;
        */

        indicesPerTriangle = mesh.vertices.size() / mesh.triangles.size();

        // For every triangle in this mesh
        for (auto& tri : mesh.triangles) {
            //Add triangle to vector 
            //std::cout << "TRIANGLE FOUND" << std::endl;


            // Times 3 because thats how many indices there are
            input.x = tri.x + vectorOffset * indicesPerTriangle;
            input.y = tri.y + vectorOffset * indicesPerTriangle;
            input.z = tri.z + vectorOffset * indicesPerTriangle;
            

            /*
            std::cout << "Triangle " << j << " x: " << input.x << std::endl;
            std::cout << "Triangle " << j << " y: " << input.y << std::endl;
            std::cout << "Triangle " << j << " z: " << input.z << std::endl;
            j++;
            */

            allTriangles.push_back(input);
            totalTriangles++;
        }
    }

    
    for (auto& mesh : pScene->meshes) {
        //std::cout << "MESH FOUND" << std::endl;

        //vectorOffset = allVertices.size();
        //std::cout << "Vertex FOUND " << mesh.vertices.size() << std::endl;
        // For every triangle in this mesh
        for (auto& ver : mesh.vertices) {
            //Add triangle to vector 
            
            allVertices.push_back(ver);


        }
    }
  

    //for (int i = 0; i < allTriangles.size(); i++) {
    //    std::cout << "Triangle " << i << " x: " << allTriangles[i].x << std::endl;
    //    std::cout << "Triangle " << i << " y: " << allTriangles[i].y << std::endl;
    //    std::cout << "Triangle " << i << " z: " << allTriangles[i].z << std::endl;
    //}

    std::cout << std::endl << "             TOTAL MESHES: " << totalMeshes << std::endl;
    std::cout << std::endl << "             TOTAL TRIANGLES: " << totalTriangles << std::endl;
    std::cout << std::endl << "             TRIANGLE VECTOR SIZE: " << allTriangles.size() << std::endl;
    std::cout << std::endl << "             VERTICES VECTOR SIZE: " << allVertices.size() << std::endl;
}

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    // as an example of how to iterate over all meshes in the scene, look at the intersect method below
    
    //for (int a = 0; a < pScene->meshes.size(); a++) {
    //    std::cout << "MESHES COUNTER, " << a << std::endl;
    //}

    //MAKE VECTOR WIT HALL THGE TRIANGLES
    createVector(m_pScene);


    ResetValues();
    SetFirstBox(m_pScene);

    SetFirstTwoBox(m_pScene);
   // std::cout << "First two are created." << std::endl;
    //first 2 are created, now keep dividing them

    totalLevels = 0;
    //next devidde is on y axis
    axisDevide = 1;
    //was 1

    //std::cout << "NUMBER OF MESHEES:  " << m_pScene->meshes.size() << std::endl;

    //std::cout << "NUMBER OF TRIANGLES:  " << m_pScene->meshes[0].vertices.size() << std::endl;

    splitBVH();
}

// Use this function to visualize your BVH. This can be useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
/////// loop over all the boxes
void BoundingVolumeHierarchy::debugDraw(int level)
{
    /*
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);
    // Draw the AABB as a (white) wireframe box.
    ////AxisAlignedBox aabb { glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe);
    ////drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1);
    //iterate over all the AABB's that are created. On each, call drawAABB()!
    //AxisAlignedBox aabb{ glm::vec3(0.0f), glm::vec3(0.5f, 0.5f, 0.5f) };
    //drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1);
    //AxisAlignedBox aabb{ glm::vec3(minX, minY, minZ), glm::vec3(maxX, maxY, maxZ) };
    //AxisAlignedBox aabb{ glm::vec3(0.0f), glm::vec3(0.5f, 0.5f, 0.5f) };
    //drawAABB(aabb, DrawMode::Filled, glm::vec3(0.5f, 1.0f, 0.5f), .3);
    std::cout << "CALLED " << std::endl;
    std::cout << "min--- x:" << minX << ", y: " << minY << ", z: " << minZ << std::endl;
    std::cout << "max--- x:" << maxX << ", y: " << maxY << ", z " << maxZ << std::endl;
    //float averageX = (std::abs(maxX) + std::abs(minX)) / 2;
    //std::cout << "MAXX " << maxX << std::endl;
    //std::cout << "MINX " << minX << std::endl;
    //std::cout << "AVGX " << averageX << std::endl;
    //subdivide it creating two new child AABBs
    //split on the X axis first
    //AxisAlignedBox split1a = { glm::vec3(minX, minY, minZ), glm::vec3(maxX - averageX, maxY, maxZ) };
    //AxisAlignedBox split1b = { glm::vec3(minX + averageX, minY, minZ), glm::vec3(maxX, maxY, maxZ) };
    //drawAABB(split1a, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 1.0f), .3);
    //drawAABB(split1b, DrawMode::Filled, glm::vec3(0.5f, 1.0f, 0.5f), .3);
    */
    
    //check all the nodes if they are leaf nodes
    //std::cout << "NODES size: " << nodes.size() << std::endl;
    for (float a = 0.0; a < nodes.size(); a++) {
        //if leaf
        //if (nodes[a].interior == false) {
        if(nodes[a].level == level){


            //std::cout << "CHILDREN: "  << nodes[a].indices.size() << std::endl;


            //if (nodes[a].indices.size() == 1) {
                AxisAlignedBox draw = { 
                    nodes[a].aabb.lower,
                    nodes[a].aabb.upper };
                //AxisAlignedBox draw = { glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(1.0f, 1.0f, 1.0f) };
                drawAABB(draw, DrawMode::Wireframe, glm::vec3(1.0f, 1.0f, 0.5f), 0.3f);
            //}

        }


    }

    
    //std::cout << "RETURN " << std::endl;
}

//The amount the slider displays 
int BoundingVolumeHierarchy::numLevels() const
{
    return totalLevels;
}

// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h .
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo) const
{


    intersectedNodes.clear();
    Material res;
    res.kd = {1.0f, 1.0f, 1.0f};
    //res.ks = { 0.5f };

    
    bool hit = false;
    

    /*
    // Intersect with all triangles of all meshes.
    for (const auto& mesh : m_pScene->meshes) {
        //every mesh has triangles and for that triangle do:
        for (const auto& tri : mesh.triangles) {
            const auto v0 = mesh.vertices[tri[0]];
            const auto v1 = mesh.vertices[tri[1]];
            const auto v2 = mesh.vertices[tri[2]];
            if (intersectRayWithTriangle(v0.p, v1.p, v2.p, ray, hitInfo)) {
                hitInfo.material = mesh.material;
<<<<<<< HEAD

                glm::vec3 P = ray.origin + ray.t * ray.direction;

                float A_tot = glm::length(glm::cross(v0.p - v1.p, v0.p - v2.p));
                float A0 = glm::length(glm::cross(v1.p - P, v2.p - P)) / A_tot;
                float A1 = glm::length(glm::cross(v0.p - P, v2.p - P)) / A_tot;
                float A2 = glm::length(glm::cross(v0.p - P, v1.p - P)) / A_tot;

                glm::vec3 N0 = A0 * v0.n;
                glm::vec3 N1 = A1 * v1.n;
                glm::vec3 N2 = A2 * v2.n;

                glm::vec3 N = N0 + N1 + N2;

                hitInfo.normal = N;

                if (glm::dot(hitInfo.normal, ray.direction) < 0) {
                    hitInfo.normal = glm::normalize(hitInfo.normal);
                } else {
                    hitInfo.normal = -glm::normalize(hitInfo.normal);
                }

=======
                //hitInfo.material = res;
>>>>>>> 19b045af3b0dee6f20167929bb5b3a2cc68770b0
                hit = true;
            }
        }
    }
    return hit;
    */
    

    //Traverse the tree
    //start at the root
    
    //base case
    //if (firstCall) {

    //}
    //std::cout << nodes.size() << std::endl;
    //std::cout << nodes[0].indices.size() << std::endl;
    auto root = nodes[0];
    /*
    std::cout << "WE ARE INSUDE TJE FIRST: " << std::endl;
    std::cout << "MinX: " << root.aabb.lower.x << std::endl;
    std::cout << "MinY: " << root.aabb.lower.y << std::endl;
    std::cout << "MinZ: " << root.aabb.lower.z << std::endl;
    std::cout << "MaxX: " << root.aabb.upper.x << std::endl;
    std::cout << "MaxY: " << root.aabb.upper.y << std::endl;
    std::cout << "MaxZ: " << root.aabb.upper.z << std::endl;
    */

    //std::cout << " Index left child: " << root.indices[0] << std::endl;
    //std::cout << " Index right child: " << root.indices[1] << std::endl;
    //now check the intersection over this aabb
    Node theNode = getAabb(root, ray);

    /*
    std::cout << "    WWHEN RETURNDED THIS IS THE VAL" << std::endl;
    std::cout << "MinX: " << theAabb.lower.x << std::endl;
    std::cout << "MinY: " << theAabb.lower.y << std::endl;
    std::cout << "MinZ: " << theAabb.lower.z << std::endl;
    std::cout << "MaxX: " << theAabb.upper.x << std::endl;
    std::cout << "MaxY: " << theAabb.upper.y << std::endl;
    std::cout << "MaxZ: " << theAabb.upper.z << std::endl;
    */

    //use finalIntersec
    //bool intersect = finalIntersect(theNode.aabb, ray);
    AxisAlignedBox draw = {
        theNode.aabb.lower,
        theNode.aabb.upper };
    drawAABB(draw, DrawMode::Wireframe, glm::vec3(1.0f, 1.0f, 1.0f), 0.8f);

    //std::cout << " Intersect at t: " << ray.t << std::endl;
    //std::cout << std::endl << std::endl << std::endl << std::endl;

    /*
    std::cout << "Node size: " << theNode.indices.size() << std::endl;
    std::cout << "Node First Index: " << theNode.indices[0] << std::endl;
    std::cout << "Ray.t: " << ray.t << std::endl;
    std::cout << "Ray.dir. x: " << ray.direction.x << ", y: " << ray.direction.y << ", z: " << ray.direction.z << std::endl;
    std::cout << "Ray.origin. x: " << ray.origin.x << ", y: " << ray.origin.y << ", z: " << ray.origin.z << std::endl << std::endl;
    */

    
    //loop over all the triangles in this node and return intersection
    for (int j = 0; j < theNode.indices.size(); j++) {

        //loop over every vertex
        //for (int k = 0; k < 3; k++) {
        //allTriangles
        //allVertices[allTriangles[triIndex][x]]

        glm::vec3 first = allVertices[allTriangles[theNode.indices[j]][0]].p;
        glm::vec3 second = allVertices[allTriangles[theNode.indices[j]][1]].p;
        glm::vec3 third = allVertices[allTriangles[theNode.indices[j]][2]].p;

        /*
        m_pScene->meshes[0].triangles[];

        allVertices[allTriangles[j][0]].p.x = 10;
        allVertices[allTriangles[j][0]].p.y = 10;
        allVertices[allTriangles[j][0]].p.z = 10;
        allVertices[allTriangles[j][1]].p.x = 10;
        allVertices[allTriangles[j][1]].p.y = 10;
        allVertices[allTriangles[j][1]].p.z = 10;
        allVertices[allTriangles[j][2]].p.x = 10;
        allVertices[allTriangles[j][2]].p.y = 10;
        allVertices[allTriangles[j][2]].p.z = 10;
        */
        /*
        std::cout << "First vertex. x: " << first.x << ", y: " << first.y << ", z: " << first.z << std::endl;
        std::cout << "Second vertex. x: " << second.x << ", y: " << second.y << ", z: " << second.z << std::endl;
        std::cout << "Third vertex. x: " << third.x << ", y: " << third.y << ", z: " << third.z << std::endl << std::endl;
        */

        if (!hit) {
            hit = intersectRayWithTriangle(first, second, third, ray, hitInfo);
        }

        if (hit) {
            if (intersectRayWithTriangle(first, second, third, ray, hitInfo)) {
                hit = intersectRayWithTriangle(first, second, third, ray, hitInfo);
            }
        }

        //hit = intersectRayWithTriangle(first, second, third, ray, hitInfo);
    }
    
    /*
    // IDK why it doesnt work. Instead itterate over ALL leafs
    //if false, check all the intersected boxes
    //no intersection. Get all the leaf nodes the ray goes through
    if (!hit) {
        //traverse whole nodes[a] and get the leafs
        for (int m = 0; m < nodes.size(); m++) {
            //if leaf
            if (!nodes[m].interior) {
                //check intersection with this aabb
                if (testIntersect(nodes[m], ray) < std::numeric_limits<float>::max()) {
                    //it intersects so add this boi AND create

                    //Node add;

                    //for (int q = 0; q < nodes[m].indices.size(); q++) {
                    //    add.indices[q] = nodes[m].indices[q];
                    //}


                    allLeafs.push_back(nodes[m]);
                    //draw them
                    //AxisAlignedBox draw = {
                    //nodes[m].aabb.lower,
                    //nodes[m].aabb.upper };
                   // drawAABB(draw, DrawMode::Wireframe, glm::vec3(1.0f, 1.0f, 1.0f), 0.8f);
                }
            }
        }

        //now traverse this vector with the nodes
        //for every leaf
        for (int n = 0; n < allLeafs.size(); n++) {
            //traverse the triangles to look for hit
            for (int p = 0; p < allLeafs[n].indices.size(); p++) {

                Node test = allLeafs[n];

                glm::vec3 first = allVertices[allTriangles[test.indices[p]][0]].p;
                glm::vec3 second = allVertices[allTriangles[test.indices[p]][1]].p;
                glm::vec3 third = allVertices[allTriangles[test.indices[p]][2]].p;
                hit = intersectRayWithTriangle(first, second, third, ray, hitInfo);
            }
        }
    }
    */
    

    // Intersect with spheres.
    //for (const auto& sphere : m_pScene->spheres) {
    //    hit |= intersectRayWithShape(sphere, ray, hitInfo);
    //}

    const auto& mesh = m_pScene->meshes[0];
    hitInfo.material = mesh.material;

    return hit;
}

BoundingVolumeHierarchy::Node BoundingVolumeHierarchy::getAabb(Node root, Ray& ray) const {
    Node res;
    res.aabb.lower.x = 0;
    res.aabb.lower.y = 0;
    res.aabb.lower.z = 0;
    res.aabb.upper.x = 0;
    res.aabb.upper.y = 0;
    res.aabb.upper.z = 0;
    //base case is lowest level, aka leaf
    if (!root.interior) {

        /*
        std::cout << "MinX: " << root.aabb.lower.x << std::endl;
        std::cout << "MinY: " << root.aabb.lower.y << std::endl;
        std::cout << "MinZ: " << root.aabb.lower.z << std::endl;
        std::cout << "MaxX: " << root.aabb.upper.x << std::endl;
        std::cout << "MaxY: " << root.aabb.upper.y << std::endl;
        std::cout << "MaxZ: " << root.aabb.upper.z << std::endl;
       
        */

        //std::cout << "      -------------------    BOTTOM" << std::endl;

        //fill the vector with all the triangles in this 

        /*
        AxisAlignedBox res;
        res.lower.x = root.aabb.lower.x;
        res.lower.y = root.aabb.lower.y;
        res.lower.z = root.aabb.lower.z;
        res.upper.x = root.aabb.upper.x;
        res.upper.y = root.aabb.upper.y;
        res.upper.z = root.aabb.upper.z;
        */
        
        return root;
    }

    /*
    std::cout << "MinX: " << root.aabb.lower.x << std::endl;
    std::cout << "MinY: " << root.aabb.lower.y << std::endl;
    std::cout << "MinZ: " << root.aabb.lower.z << std::endl;
    std::cout << "MaxX: " << root.aabb.upper.x << std::endl;
    std::cout << "MaxY: " << root.aabb.upper.y << std::endl;
    std::cout << "MaxZ: " << root.aabb.upper.z << std::endl;
    */

    //check if ray hits the given box
    if (testIntersect(root, ray) < std::numeric_limits<float>::max()) {
        float leftT = testIntersect(nodes[root.indices[0]], ray);
        float rightT = testIntersect(nodes[root.indices[1]], ray);

        if (leftT < rightT) {
            //std::cout << "---- Intersect Left: " << root.indices[0] << std::endl;
            //intersectedNodes.push_back(root.indices[0]);
            /*
            AxisAlignedBox draw = {
            nodes[root.indices[0]].aabb.lower,
            nodes[root.indices[0]].aabb.upper };
            drawAABB(draw, DrawMode::Wireframe, glm::vec3(1.0f, 1.0f, 1.0f), 0.8f);
            */
            return res = getAabb(nodes[root.indices[0]], ray);
        } 
        //right
        else if (leftT > rightT) {
            //std::cout << "---- Intersect Right: " << root.indices[1] << std::endl;
            //intersectedNodes.push_back(root.indices[1]);
            /*
            AxisAlignedBox draw = {
            nodes[root.indices[1]].aabb.lower,
            nodes[root.indices[1]].aabb.upper };
            drawAABB(draw, DrawMode::Wireframe, glm::vec3(1.0f, 1.0f, 1.0f), 0.8f);
            */
            return res = getAabb(nodes[root.indices[1]], ray);
        }
        //both t values are the MAX. This means that the parent aabb gets intersected but not one of the child aabb
        else if (leftT == rightT && leftT == std::numeric_limits<float>::max() && rightT == std::numeric_limits<float>::max()) {
            //std::cout << "---- Intersection with parent box, but not children. Return Parent box" << std::endl;

            /*
            float leftT1 = std::numeric_limits<float>::max();
            float leftT2 = std::numeric_limits<float>::max();
            float rightT1 = std::numeric_limits<float>::max();
            float rightT2 = std::numeric_limits<float>::max();
            */

            /*
            std::cout << "---- Interior?: " << nodes[root.indices[0]].interior << std::endl;
            leftT1 = testIntersect(nodes[nodes[root.indices[0]].indices[0]], ray);
            if (nodes[root.indices[0]].interior) {
                leftT2 = testIntersect(nodes[nodes[root.indices[0]].indices[1]], ray);
            }

            std::cout << "---- Interior?: " << nodes[root.indices[1]].interior << std::endl;
            rightT1 = testIntersect(nodes[nodes[root.indices[1]].indices[0]], ray);
            if (nodes[root.indices[1]].interior) {
                rightT2 = testIntersect(nodes[nodes[root.indices[1]].indices[1]], ray);
            }

            if (leftT1 < leftT2, rightT1, rightT2) {
                std::cout << "---- left1 best" << std::endl;
                return res = getAabb(nodes[nodes[root.indices[0]].indices[0]], ray);
            } 
            else if (leftT2 < leftT1, rightT1, rightT2) {
                std::cout << "---- left2 best" << std::endl;
                return res = getAabb(nodes[nodes[root.indices[0]].indices[0]], ray);
            }
            else if (rightT1 < leftT1, leftT2, rightT2) {
                std::cout << "---- right1 best" << std::endl;
                return res = getAabb(nodes[nodes[root.indices[0]].indices[0]], ray);
            }
            else if (rightT2 < leftT1, leftT2, rightT1) {
                std::cout << "---- right2 best" << std::endl;
                return res = getAabb(nodes[nodes[root.indices[0]].indices[0]], ray);
            }
            */
            
            //float the t to
            float tToCheck = std::numeric_limits<float>::max();
            float noIntersect = 0;
            int theNode = 0;
            // check all the aabb's because you dont know
            for (int i = 0; i < nodes.size(); i++) {
                //if leaf check intersect
                if (!nodes[i].interior) {
                    
                    //tToCheck = ray.t;

                    //first find all intersection points with the axis
                    float tXmin = (nodes[i].aabb.lower.x - ray.origin.x) / ray.direction.x;
                    float tXmax = (nodes[i].aabb.upper.x - ray.origin.x) / ray.direction.x;
                    float tYmin = (nodes[i].aabb.lower.y - ray.origin.y) / ray.direction.y;
                    float tYmax = (nodes[i].aabb.upper.y - ray.origin.y) / ray.direction.y;
                    float tZmin = (nodes[i].aabb.lower.z - ray.origin.z) / ray.direction.z;
                    float tZmax = (nodes[i].aabb.upper.z - ray.origin.z) / ray.direction.z;

                    //find all the in and out points
                    float tInX = std::min(tXmin, tXmax);
                    float tOutX = std::max(tXmin, tXmax);
                    float tInY = std::min(tYmin, tYmax);
                    float tOutY = std::max(tYmin, tYmax);
                    float tInZ = std::min(tZmin, tZmax);
                    float tOutZ = std::max(tZmin, tZmax);

                    //In the box after crossing ALL the in points
                    //Out the box after crossing ONE out point

                    //Find the global T of the ray
                    //Global tIn is the smallest t where we crossed ALL the in points
                    //Use lists to compare multiple values
                    float tIn = std::max({ tInX, tInY, tInZ });

                    //Global tOut is the smallest t where we crossed ATLEAST ONE of the out points
                    float tOut = std::min({ tOutX, tOutY, tOutZ });

                    //Now we know the in and out T-values, we need to check if we mis the box or not
                    if (tIn > tOut || tOut < 0) {
                        //std::cout << "We miss the box:" << std::endl;
                        //tToCheck = std::numeric_limits<float>::max();
                        noIntersect = std::numeric_limits<float>::max();
                    }

                    /*
                    std::cout << "tInX " << tInX << std::endl;
                    std::cout << "tInY " << tInY << std::endl;
                    std::cout << "tInZ " << tInZ << std::endl;
                    std::cout << "tOutX " << tOutX << std::endl;
                    std::cout << "tOutY " << tOutY << std::endl;
                    std::cout << "tOutZ " << tOutZ << std::endl;
                    */

                    //IF one of the ins is negative, its inside the box
                    // thus tout = t
                    //std::cout << "THE T NOW =  " << tNow << std::endl;
                    //std::cout << "tIn =  " << tIn << std::endl;
                    if (tIn < 0) {
                        //only assign t if its smaller


                        //std::cout << "Assigned T When INSIDE, OUT: " << tOut << std::endl;
                        //ray.t = tOut;
                        //tToCheck = std::numeric_limits<float>::max();
                        noIntersect = std::numeric_limits<float>::max();
                    }


                    //intersection and smaller t
                    //std::cout << "Assigned T: " << tIn << std::endl;
                    //std::cout << "TRUE:" << tIn << std::endl;
                    //ray.t = tIn;
                    //if (tToCheck > tIn) {
                    //    tToCheck = tIn;
                    //    theNode = i;
                    //}
                    
                    //means it intersects. Add to vector and find closest one
                    //if (tToCheck != std::numeric_limits<float>::max()) {
                    //    intersectedNodes.push_back(i);
                    //}

                    //float checkit = testIntersect(nodes[i], ray);

                    //meaning there was a intersection
                    if (noIntersect == 0) {
                        if (tToCheck > tIn) {
                            theNode = i;
                            tToCheck = tIn;
                        }
                        //intersectedNodes.push_back(i);
                    }

                    //if (tToCheck > checkit) {
                    //    tToCheck = checkit;
                    //    theNode = i;
                    //}

                }
                noIntersect = 0;
            }
            //std::cout << "---- WE CHECKED AND THE INDEX OF THE NODE = " << theNode << std::endl;
            if (theNode == 0) {
                return res;
            }
            return res = getAabb(nodes[theNode], ray);
        }
        
        else if (leftT == rightT) {
            //std::cout << "Left and right same T" << std::endl;
            //intersectedNodes.push_back(root.indices[0]);

            return res = getAabb(nodes[root.indices[0]], ray);
            //return res = root;
        }
    }
    else {
        //std::cout << "NO RESULT" << std::endl;
        return res;
    }
    return res;


}

//returns the t of the intersection wit hthe aabb from this node
float BoundingVolumeHierarchy::testIntersect(Node root, Ray& ray) const {
    //this is the T value to revert to 
    float tNow = ray.t;
    /*
    std::cout << "NOW WE CHECK FOR THE INTERSECTION: "  << std::endl;
    std::cout << "MinX: " << root.aabb.lower.x << std::endl;
    std::cout << "MinY: " << root.aabb.lower.y << std::endl;
    std::cout << "MinZ: " << root.aabb.lower.z << std::endl;
    std::cout << "MaxX: " << root.aabb.upper.x << std::endl;
    std::cout << "MaxY: " << root.aabb.upper.y << std::endl;
    std::cout << "MaxZ: " << root.aabb.upper.z << std::endl;
    */

    //first find all intersection points with the axis
    float tXmin = (root.aabb.lower.x - ray.origin.x) / ray.direction.x;
    float tXmax = (root.aabb.upper.x - ray.origin.x) / ray.direction.x;
    float tYmin = (root.aabb.lower.y - ray.origin.y) / ray.direction.y;
    float tYmax = (root.aabb.upper.y - ray.origin.y) / ray.direction.y;
    float tZmin = (root.aabb.lower.z - ray.origin.z) / ray.direction.z;
    float tZmax = (root.aabb.upper.z - ray.origin.z) / ray.direction.z;

    //find all the in and out points
    float tInX = std::min(tXmin, tXmax);
    float tOutX = std::max(tXmin, tXmax);
    float tInY = std::min(tYmin, tYmax);
    float tOutY = std::max(tYmin, tYmax);
    float tInZ = std::min(tZmin, tZmax);
    float tOutZ = std::max(tZmin, tZmax);

    //In the box after crossing ALL the in points
    //Out the box after crossing ONE out point

    //Find the global T of the ray
    //Global tIn is the smallest t where we crossed ALL the in points
    //Use lists to compare multiple values
    float tIn = std::max({ tInX, tInY, tInZ });

    //Global tOut is the smallest t where we crossed ATLEAST ONE of the out points
    float tOut = std::min({ tOutX, tOutY, tOutZ });

    //Now we know the in and out T-values, we need to check if we mis the box or not
    if (tIn > tOut || tOut < 0) {
        //std::cout << "We miss the box:" << std::endl;
        return std::numeric_limits<float>::max();
    }

    /*
    std::cout << "tInX " << tInX << std::endl;
    std::cout << "tInY " << tInY << std::endl;
    std::cout << "tInZ " << tInZ << std::endl;
    std::cout << "tOutX " << tOutX << std::endl;
    std::cout << "tOutY " << tOutY << std::endl;
    std::cout << "tOutZ " << tOutZ << std::endl;
    */

    //IF one of the ins is negative, its inside the box
    // thus tout = t
    //std::cout << "THE T NOW =  " << tNow << std::endl;
    //std::cout << "tIn =  " << tIn << std::endl;
    if (tIn < 0) {
        //only assign t if its smaller
        if (tOut > tNow) {
            //it NOT closer, dont update
            //std::cout << "NOT BIGGER, asssign nothing, tOut:" << tOut << std::endl;
            return std::numeric_limits<float>::max();
        }

        //std::cout << "Assigned T When INSIDE, OUT: " << tOut << std::endl;
        //ray.t = tOut;
        //return tOut;
    }

    //only assign t if its smaller
    if (tIn > tNow) {
        //it NOT closer, revert update
        //ray.t = tNow;
        //std::cout << "NOT BIGGER, asssign nothing, tIn:" << tIn << std::endl;
        return std::numeric_limits<float>::max();
    }

    //intersection and smaller t
    //std::cout << "Assigned T: " << tIn << std::endl;
    //std::cout << "TRUE:" << tIn << std::endl;
    //ray.t = tIn;
    return tIn;
}

bool BoundingVolumeHierarchy::finalIntersect(AxisAlignedBox box, Ray& ray) const {

    AxisAlignedBox draw = {
    box.lower,
    box.upper };
    //AxisAlignedBox draw = { glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(1.0f, 1.0f, 1.0f) };
    drawAABB(draw, DrawMode::Wireframe, glm::vec3(1.0f, 1.0f, 1.0f), 0.8f);


    //this is the T value to revert to 
    float tNow = ray.t;

    //first find all intersection points with the axis
    float tXmin = (box.lower.x - ray.origin.x) / ray.direction.x;
    float tXmax = (box.upper.x - ray.origin.x) / ray.direction.x;
    float tYmin = (box.lower.y - ray.origin.y) / ray.direction.y;
    float tYmax = (box.upper.y - ray.origin.y) / ray.direction.y;
    float tZmin = (box.lower.z - ray.origin.z) / ray.direction.z;
    float tZmax = (box.upper.z - ray.origin.z) / ray.direction.z;

    //find all the in and out points
    float tInX = std::min(tXmin, tXmax);
    float tOutX = std::max(tXmin, tXmax);
    float tInY = std::min(tYmin, tYmax);
    float tOutY = std::max(tYmin, tYmax);
    float tInZ = std::min(tZmin, tZmax);
    float tOutZ = std::max(tZmin, tZmax);

    //In the box after crossing ALL the in points
    //Out the box after crossing ONE out point

    //Find the global T of the ray
    //Global tIn is the smallest t where we crossed ALL the in points
    //Use lists to compare multiple values
    float tIn = std::max({ tInX, tInY, tInZ });

    //Global tOut is the smallest t where we crossed ATLEAST ONE of the out points
    float tOut = std::min({ tOutX, tOutY, tOutZ });

    //Now we know the in and out T-values, we need to check if we mis the box or not
    if (tIn > tOut || tOut < 0) {
        return false;
    }

    /*
    std::cout << "tInX " << tInX << std::endl;
    std::cout << "tInY " << tInY << std::endl;
    std::cout << "tInZ " << tInZ << std::endl;
    std::cout << "tOutX " << tOutX << std::endl;
    std::cout << "tOutY " << tOutY << std::endl;
    std::cout << "tOutZ " << tOutZ << std::endl;
    */

    //IF one of the ins is negative, its inside the box
    // thus tout = t
    //std::cout << "THE T NOW =  " << tNow << std::endl;
    //std::cout << "tIn =  " << tIn << std::endl;
    if (tIn < 0) {
        //only assign t if its smaller
        if (tOut > tNow) {
            //it NOT closer, dont update
            //std::cout << "NOT BIGGER, asssign nothing, tOut:" << tOut << std::endl;
            return false;
        }

        //std::cout << "Assigned T When INSIDE, OUT: " << tOut << std::endl;
        ray.t = tOut;
        return true;
    }

    //only assign t if its smaller
    if (tIn > tNow) {
        //it NOT closer, revert update
        //ray.t = tNow;
        //std::cout << "NOT BIGGER, asssign nothing, tIn:" << tIn << std::endl;
        return false;
    }

    //std::cout << "Assigned T: " << tIn << std::endl;

    ray.t = tIn;
    return true;
}