#ifndef ACCELERATION_STRUCTURES_H
#define ACCELERATION_STRUCTURES_H

#include <vector>

class AccelerationStructure
{
public:
    AccelerationStructure(std::vector<std::unique_ptr<const Mesh>>& m) : meshes(std::move(m)) {}
    virtual ~AccelerationStructure() {}
    virtual bool intersect(const Vec3f& orig, const Vec3f& dir, const uint32_t& rayId, float& tHit) const
    {
        const Mesh* intersectedMesh = nullptr;
        float t = kInfinity;
        for (const auto& mesh: meshes) {
            if (mesh->intersect(orig, dir, t) && t < tHit) {
                intersectedMesh = mesh.get();
                tHit = t;
            }
        }

        return (intersectedMesh != nullptr);
    }
protected:
    const std::vector<std::unique_ptr<const Mesh>> meshes;
};


class BBoxAcceleration : public AccelerationStructure
{
public:
    BBoxAcceleration(std::vector<std::unique_ptr<const Mesh>>& m) : AccelerationStructure(m) {}
    virtual bool intersect(const Vec3f& orig, const Vec3f& dir, const uint32_t& rayId, float& tHit) const
    {
        const Mesh* intersectedMesh = nullptr;
        const Vec3f invDir = 1 / dir;
        const Vec3b sign(dir.x < 0, dir.y < 0, dir.z < 0);
        float t = kInfinity;
        for (const auto& mesh : meshes) {
            // If you intersect the box
            if (mesh->bbox.intersect(orig, invDir, sign, t)) {
                // Then test if the ray intersects the mesh and if does then first check
                // if the intersection distance is the nearest and if we pass that test as well
                // then update tNear variable with t and keep a pointer to the intersected mesh
                if (mesh->intersect(orig, dir, t) && t < tHit) {
                    tHit = t;
                    intersectedMesh = mesh.get();
                }
            }
        }

        // Return true if the variable intersectedMesh is not null, false otherwise
        return (intersectedMesh != nullptr);
    }
};


class BVH : public AccelerationStructure
{
    static const uint8_t kNumPlaneSetNormals = 7;
    static const Vec3f planeSetNormals[kNumPlaneSetNormals];
    struct Extents
    {
        Extents()
        {
            for (uint8_t i = 0;  i < kNumPlaneSetNormals; ++i)
                d[i][0] = kInfinity, d[i][1] = -kInfinity;
        }
        void extendBy(const Extents& e)
        {

            for (uint8_t i = 0;  i < kNumPlaneSetNormals; ++i) {
                if (e.d[i][0] < d[i][0]) d[i][0] = e.d[i][0];
                if (e.d[i][1] > d[i][1]) d[i][1] = e.d[i][1];
            }
        }
        /* inline */
        Vec3f centroid() const
        {
            return Vec3f(
                d[0][0] + d[0][1] * 0.5,
                d[1][0] + d[1][1] * 0.5,
                d[2][0] + d[2][1] * 0.5);
        }
        bool intersect(const float*, const float*, float&, float&, uint8_t&) const;
        float d[kNumPlaneSetNormals][2];
        const Mesh* mesh;
    };

    struct Octree
    {
        Octree(const Extents& sceneExtents)
        {
            float xDiff = sceneExtents.d[0][1] - sceneExtents.d[0][0];
            float yDiff = sceneExtents.d[1][1] - sceneExtents.d[1][0];
            float zDiff = sceneExtents.d[2][1] - sceneExtents.d[2][0];
            float maxDiff = std::max(xDiff, std::max(yDiff, zDiff));
            Vec3f minPlusMax(
                sceneExtents.d[0][0] + sceneExtents.d[0][1],
                sceneExtents.d[1][0] + sceneExtents.d[1][1],
                sceneExtents.d[2][0] + sceneExtents.d[2][1]);
            bbox[0] = (minPlusMax - maxDiff) * 0.5;
            bbox[1] = (minPlusMax + maxDiff) * 0.5;
            root = new OctreeNode;
        }

        ~Octree() { deleteOctreeNode(root); }

        void insert(const Extents* extents) { insert(root, extents, bbox, 0); }
        void build() { build(root, bbox); }

        struct OctreeNode
        {
            OctreeNode* child[8] = { nullptr };
            std::vector<const Extents *> nodeExtentsList; // pointer to the objects extents
            Extents nodeExtents; // extents of the octree node itself
            bool isLeaf = true;
        };

        struct QueueElement
        {
            const OctreeNode *node; // octree node held by this element in the queue
            float t; // distance from the ray origin to the extents of the node
            QueueElement(const OctreeNode *n, float tn) : node(n), t(tn) {}
            // priority_queue behaves like a min-heap
            friend bool operator < (const QueueElement &a, const QueueElement &b) { return a.t > b.t; }
        };

        OctreeNode* root = nullptr; // make unique son don't have to manage deallocation
        BBox<> bbox;

    private:

        void deleteOctreeNode(OctreeNode*& node)
        {
            for (uint8_t i = 0; i < 8; i++) {
                if (node->child[i] != nullptr) {
                    deleteOctreeNode(node->child[i]);
                }
            }
            delete node;
        }

        void insert(OctreeNode*& node, const Extents* extents, const BBox<>& bbox, uint32_t depth)
        {
            if (node->isLeaf) {
                if (node->nodeExtentsList.size() == 0 || depth == 16) {
                    node->nodeExtentsList.push_back(extents);
                }
                else {
                    node->isLeaf = false;
                    // Re-insert extents held by this node
                    while (node->nodeExtentsList.size()) {
                        insert(node, node->nodeExtentsList.back(), bbox, depth);
                        node->nodeExtentsList.pop_back();
                    }
                    // Insert new extent
                    insert(node, extents, bbox, depth);
                }
            }
            else {
                // Need to compute in which child of the current node this extents should
                // be inserted into
                Vec3f extentsCentroid = extents->centroid();
                Vec3f nodeCentroid = (bbox[0] + bbox[1]) * 0.5;
                BBox<> childBBox;
                uint8_t childIndex = 0;
                // x-axis
                if (extentsCentroid.x > nodeCentroid.x) {
                    childIndex = 4;
                    childBBox[0].x = nodeCentroid.x;
                    childBBox[1].x = bbox[1].x;
                }
                else {
                    childBBox[0].x = bbox[0].x;
                    childBBox[1].x = nodeCentroid.x;
                }
                // y-axis
                if (extentsCentroid.y > nodeCentroid.y) {
                    childIndex += 2;
                    childBBox[0].y = nodeCentroid.y;
                    childBBox[1].y = bbox[1].y;
                }
                else {
                    childBBox[0].y = bbox[0].y;
                    childBBox[1].y = nodeCentroid.y;
                }
                // z-axis
                if (extentsCentroid.z > nodeCentroid.z) {
                    childIndex += 1;
                    childBBox[0].z = nodeCentroid.z;
                    childBBox[1].z = bbox[1].z;
                }
                else {
                    childBBox[0].z = bbox[0].z;
                    childBBox[1].z = nodeCentroid.z;
                }

                // Create the child node if it doesn't exsit yet and then insert the extents in it
                if (node->child[childIndex] == nullptr)
                    node->child[childIndex] = new OctreeNode;
                insert(node->child[childIndex], extents, childBBox, depth + 1);
            }
        }

        void build(OctreeNode*& node, const BBox<>& bbox)
        {
            if (node->isLeaf) {
                for (const auto& e: node->nodeExtentsList) {
                    node->nodeExtents.extendBy(*e);
                }
            }
            else {
                for (uint8_t i = 0; i < 8; ++i) {
                        if (node->child[i]) {
                        BBox<> childBBox;
                        Vec3f centroid = bbox.centroid();
                        // x-axis
                        childBBox[0].x = (i & 4) ? centroid.x : bbox[0].x;
                        childBBox[1].x = (i & 4) ? bbox[1].x : centroid.x;
                        // y-axis
                        childBBox[0].y = (i & 2) ? centroid.y : bbox[0].y;
                        childBBox[1].y = (i & 2) ? bbox[1].y : centroid.y;
                        // z-axis
                        childBBox[0].z = (i & 1) ? centroid.z : bbox[0].z;
                        childBBox[1].z = (i & 1) ? bbox[1].z : centroid.z;

                        // Inspect child
                        build(node->child[i], childBBox);

                        // Expand extents with extents of child
                        node->nodeExtents.extendBy(node->child[i]->nodeExtents);
                    }
                }
            }
        }
    };

    std::vector<Extents> extentsList;
    Octree* octree = nullptr;
public:
    BVH(std::vector<std::unique_ptr<const Mesh>>& m);
    bool intersect(const Vec3f&, const Vec3f&, const uint32_t&, float&) const;
    ~BVH() { delete octree; }
};

const Vec3f BVH::planeSetNormals[BVH::kNumPlaneSetNormals] = {
    Vec3f(1, 0, 0),
    Vec3f(0, 1, 0),
    Vec3f(0, 0, 1),
    Vec3f( sqrtf(3) / 3.f,  sqrtf(3) / 3.f, sqrtf(3) / 3.f),
    Vec3f(-sqrtf(3) / 3.f,  sqrtf(3) / 3.f, sqrtf(3) / 3.f),
    Vec3f(-sqrtf(3) / 3.f, -sqrtf(3) / 3.f, sqrtf(3) / 3.f),
    Vec3f( sqrtf(3) / 3.f, -sqrtf(3) / 3.f, sqrtf(3) / 3.f)
};

BVH::BVH(std::vector<std::unique_ptr<const Mesh>>& m) : AccelerationStructure(m)
{
    Extents sceneExtents; // that's the extent of the entire scene which we need to compute for the octree
    extentsList.reserve(meshes.size());
    for (uint32_t i = 0; i < meshes.size(); ++i) {
        for (uint8_t j = 0; j < kNumPlaneSetNormals; ++j) {
            for (const auto vtx : meshes[i]->vertexPool) {
                float d = dot(planeSetNormals[j], vtx);
                // set dNEar and dFar
                if (d < extentsList[i].d[j][0]) extentsList[i].d[j][0] = d;
                if (d > extentsList[i].d[j][1]) extentsList[i].d[j][1] = d;
            }
        }
        sceneExtents.extendBy(extentsList[i]); // expand the scene extent of this object's extent
        extentsList[i].mesh = meshes[i].get(); // the extent itself needs to keep a pointer to the object its holds
    }

    // Now that we have the extent of the scene we can start building our octree
    // Using C++ make_unique function here but you don't need to, just to learn something...
    octree = new Octree(sceneExtents);

    for (uint32_t i = 0; i < meshes.size(); ++i) {
        octree->insert(&extentsList[i]);
    }

    // Build from bottom up
    octree->build();
}

bool BVH::Extents::intersect(
    const float* precomputedNumerator,
    const float* precomputedDenominator,
    float& tNear,   // tn and tf in this method need to be contained
    float& tFar,    // within the range [tNear:tFar]
    uint8_t& planeIndex) const
{
    numRayBoundingVolumeTests++;
    for (uint8_t i = 0; i < kNumPlaneSetNormals; ++i) {
        float tNearExtents = (d[i][0] - precomputedNumerator[i]) / precomputedDenominator[i];
        float tFarExtents = (d[i][1] - precomputedNumerator[i]) / precomputedDenominator[i];
        if (precomputedDenominator[i] < 0) std::swap(tNearExtents, tFarExtents);
        if (tNearExtents > tNear) tNear = tNearExtents, planeIndex = i;
        if (tFarExtents < tFar) tFar = tFarExtents;
        if (tNear > tFar) return false;
    }

    return true;
}

bool BVH::intersect(const Vec3f& orig, const Vec3f& dir, const uint32_t& rayId, float& tHit) const
{
    tHit = kInfinity;
    const Mesh* intersectedMesh = nullptr;
    float precomputedNumerator[BVH::kNumPlaneSetNormals];
    float precomputedDenominator[BVH::kNumPlaneSetNormals];
    for (uint8_t i = 0; i < kNumPlaneSetNormals; ++i) {
        precomputedNumerator[i] = dot(planeSetNormals[i], orig);
        precomputedDenominator[i] = dot(planeSetNormals[i], dir);
    }

    /*
    tNear = kInfinity; // set
    for (uint32_t i = 0; i < meshes.size(); ++i) {
        numRayVolumeTests++;
        float tn = -kInfinity, tf = kInfinity;
        uint8_t planeIndex;
        if (extents[i].intersect(precomputedNumerator, precomputedDenominator, tn, tf, planeIndex)) {
            if (tn < tNear) {
                intersectedMesh = meshes[i].get();
                tNear = tn;
                // normal = planeSetNormals[planeIndex];
            }
        }
    }
    */

    uint8_t planeIndex;
    float tNear = 0, tFar = kInfinity; // tNear, tFar for the intersected extents
    if (!octree->root->nodeExtents.intersect(precomputedNumerator, precomputedDenominator, tNear, tFar, planeIndex) || tFar < 0)
        return false;
    tHit = tFar;
    std::priority_queue<BVH::Octree::QueueElement> queue;
    queue.push(BVH::Octree::QueueElement(octree->root, 0));
    while (!queue.empty() && queue.top().t < tHit) {
        const Octree::OctreeNode *node = queue.top().node;
        queue.pop();
        if (node->isLeaf) {
            for (const auto& e: node->nodeExtentsList) {
                float t = kInfinity;
                if (e->mesh->intersect(orig, dir, t) && t < tHit) {
                    tHit = t;
                    intersectedMesh = e->mesh;
                }
            }
        }
        else {
            for (uint8_t i = 0; i < 8; ++i) {
                if (node->child[i] != nullptr) {
                    float tNearChild = 0, tFarChild = tFar;
                    if (node->child[i]->nodeExtents.intersect(precomputedNumerator, precomputedDenominator, tNearChild, tFarChild, planeIndex)) {
                        float t = (tNearChild < 0 && tFarChild >= 0) ? tFarChild : tNearChild;
                        queue.push(BVH::Octree::QueueElement(node->child[i], t));
                    }
                }
            }
        }
    }

    return (intersectedMesh != nullptr);
}


class Grid : public AccelerationStructure
{
    struct Cell
    {
        Cell() {}
        struct TriangleDesc
        {
            TriangleDesc(const Mesh* m, const uint32_t &t) : mesh(m), tri(t) {}
            const Mesh* mesh;
            uint32_t tri;
        };

        void insert(const Mesh* mesh, uint32_t t)
        { triangles.push_back(Grid::Cell::TriangleDesc(mesh, t)); }

        bool intersect(const Vec3f&, const Vec3f&, const uint32_t&, float&, const Mesh*&) const;

        std::vector<TriangleDesc> triangles;
    };
public:
    Grid(std::vector<std::unique_ptr<const Mesh>>& m);
    ~Grid()
    {
        for (uint32_t i = 0; i < resolution[0] * resolution[1] * resolution[2]; ++i)
            if (cells[i] != NULL) delete cells[i];
        delete [] cells;
    }
    bool intersect(const Vec3f&, const Vec3f&, const uint32_t&, float&) const;
    Cell **cells;
    BBox<> bbox;
    Vec3<uint32_t> resolution;
    Vec3f cellDimension;
};

Grid::Grid(std::vector<std::unique_ptr<const Mesh>>& m) : AccelerationStructure(m)
{
    uint32_t totalNumTriangles = 0;
    for (const auto& m : meshes) {
        bbox.extendBy(m->bbox[0]);
        bbox.extendBy(m->bbox[1]);
        totalNumTriangles += m->numTriangles;
    }
    // Create the grid
    Vec3f size = bbox[1] - bbox[0];
    float cubeRoot = std::powf(totalNumTriangles / (size.x * size.y * size.z), 1. / 3.f);
    for (uint8_t i = 0; i < 3; ++i) {
        resolution[i] = std::floor(size[i] * cubeRoot);
        if (resolution[i] < 1) resolution[i] = 1;
        if (resolution[i] > 128) resolution[i] = 128;
    }
    cellDimension = size / resolution;
    uint32_t numCells = resolution.x * resolution.y * resolution.z;
    cells = new Grid::Cell* [numCells];
    memset(cells, 0x0, sizeof(Grid::Grid*) * numCells);

    for (const auto& m : meshes) {
        for (uint32_t i = 0, off = 0; i < m->numTriangles; ++i, off += 3) {
            Vec3f min(kInfinity), max(-kInfinity);
            const Vec3f& v0 = m->vertexPool[m->triangleIndicesInVertexPool[off]];
            const Vec3f& v1 = m->vertexPool[m->triangleIndicesInVertexPool[off + 1]];
            const Vec3f& v2 = m->vertexPool[m->triangleIndicesInVertexPool[off + 2]];
            for (uint8_t j = 0; j < 3; ++j) {
                if (v0[j] < min[j]) min[j] = v0[j];
                if (v1[j] < min[j]) min[j] = v1[j];
                if (v2[j] < min[j]) min[j] = v2[j];
                if (v0[j] > max[j]) max[j] = v0[j];
                if (v1[j] > max[j]) max[j] = v1[j];
                if (v2[j] > max[j]) max[j] = v2[j];
            }
            // Convert to cell coordinates
            min = (min - bbox[0]) / cellDimension;
            max = (max - bbox[0]) / cellDimension;
            uint32_t zmin = clamp<uint32_t>(std::floor(min[2]), 0, resolution[2] - 1);
            uint32_t zmax = clamp<uint32_t>(std::floor(max[2]), 0, resolution[2] - 1);
            uint32_t ymin = clamp<uint32_t>(std::floor(min[1]), 0, resolution[1] - 1);
            uint32_t ymax = clamp<uint32_t>(std::floor(max[1]), 0, resolution[1] - 1);
            uint32_t xmin = clamp<uint32_t>(std::floor(min[0]), 0, resolution[0] - 1);
            uint32_t xmax = clamp<uint32_t>(std::floor(max[0]), 0, resolution[0] - 1);
            // Loop over the cells the triangle overlaps and insert
            for (uint32_t z = zmin; z <= zmax; ++z) {
                for (uint32_t y = ymin; y <= ymax; ++y) {
                    for (uint32_t x = xmin; x <= xmax; ++x) {
                        uint32_t index = z * resolution[0] * resolution[1] + y * resolution[0] + x;
                        if (cells[index] == NULL) cells[index] = new Grid::Cell;
                        cells[index]->insert(m.get(), i);
                    }
                }
            }
        }
    }
}

bool Grid::Cell::intersect(
    const Vec3f& orig, const Vec3f& dir, const uint32_t& rayId,
    float& tHit, const Mesh*& intersectedMesh) const
{
    float uhit, vhit;
    for (uint32_t i = 0; i < triangles.size(); ++i) {
        if (rayId != triangles[i].mesh->mailbox[triangles[i].tri]) {
            triangles[i].mesh->mailbox[triangles[i].tri] = rayId;
            const Mesh *mesh = triangles[i].mesh;
            uint32_t j = triangles[i].tri * 3;
            const Vec3f &v0 = mesh->vertexPool[mesh->triangleIndicesInVertexPool[j    ]];
            const Vec3f &v1 = mesh->vertexPool[mesh->triangleIndicesInVertexPool[j + 1]];
            const Vec3f &v2 = mesh->vertexPool[mesh->triangleIndicesInVertexPool[j + 2]];
            float t, u, v;
            if (rayTriangleIntersect(orig, dir, v0, v1, v2, t, u, v)) {
                if (t < tHit) {
                    tHit = t;
                    uhit = u;
                    vhit = v;
                    intersectedMesh = triangles[i].mesh;
                }
            }
        }
    }
        return (intersectedMesh != nullptr);
}


bool Grid::intersect(const Vec3f& orig, const Vec3f& dir, const uint32_t& rayId, float& tHit) const
{
    const Vec3f invDir = 1 / dir;
    const Vec3b sign(dir.x < 0, dir.y < 0, dir.z < 0);
    float tHitBox;
    if (!bbox.intersect(orig, invDir, sign, tHitBox)) return false;

    // initialization step
    Vec3i exit, step, cell;
    Vec3f deltaT, nextCrossingT;
    for (uint8_t i = 0; i < 3; ++i) {
        // convert ray starting point to cell coordinates
        float rayOrigCell = ((orig[i] + dir[i] * tHitBox) -  bbox[0][i]);
        cell[i] = clamp<uint32_t>(std::floor(rayOrigCell / cellDimension[i]), 0, resolution[i] - 1);
        if (dir[i] < 0) {
            deltaT[i] = -cellDimension[i] * invDir[i];
            nextCrossingT[i] = tHitBox + (cell[i] * cellDimension[i] - rayOrigCell) * invDir[i];
            exit[i] = -1;
            step[i] = -1;
        }
        else {
            deltaT[i] = cellDimension[i] * invDir[i];
            nextCrossingT[i] = tHitBox + ((cell[i] + 1)  * cellDimension[i] - rayOrigCell) * invDir[i];
            exit[i] = resolution[i];
            step[i] = 1;
        }
    }

    // Walk through each cell of the grid and test for an intersection if
    // current cell contains geometry
    const Mesh* intersectedMesh = nullptr;
    while (1) {
        uint32_t o = cell[2] * resolution[0] * resolution[1] + cell[1] * resolution[0] + cell[0];
        if (cells[o] != nullptr) {
            cells[o]->intersect(orig, dir, rayId, tHit, intersectedMesh);
            //if (intersectedMesh != nullptr) { ray.color = cells[o]->color; }
        }
        uint8_t k =
            ((nextCrossingT[0] < nextCrossingT[1]) << 2) +
            ((nextCrossingT[0] < nextCrossingT[2]) << 1) +
            ((nextCrossingT[1] < nextCrossingT[2]));
        static const uint8_t map[8] = {2, 1, 2, 1, 2, 2, 0, 0};
        uint8_t axis = map[k];

        if (tHit < nextCrossingT[axis]) break;
        cell[axis] += step[axis];
        if (cell[axis] == exit[axis]) break;
        nextCrossingT[axis] += deltaT[axis];
    }

    return (intersectedMesh != nullptr);
}


#endif // ACCELERATION_STRUCTURES_H
