#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>
#include <utility>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <chrono>

#include <random>

#include "geometry.h"
#include "lights.h"
#include "objects.h"

static const Vec3f kDefaultBackgroundColor = Vec3f(0.235294, 0.67451, 0.843137);

struct IsectInfo
{
    const Object *hitObject = nullptr;
    float tNear = kInfinity;
    Vec2f uv;
    uint32_t index = 0;
};

struct Options
{
    uint32_t width = 640;
    uint32_t height = 480;
    float fov = 90;
    Vec3f backgroundColor = kDefaultBackgroundColor;
    Matrix44f cameraToWorld;
    float bias = 0.0001;
    uint32_t maxDepth = 2;
};


inline
float clamp(const float &lo, const float &hi, const float &v)
{ return std::max(lo, std::min(hi, v)); }

inline
float deg2rad(const float &deg)
{ return deg * M_PI / 180; }

inline
Vec3f mix(const Vec3f &a, const Vec3f& b, const float &mixValue)
{ return a * (1 - mixValue) + b * mixValue; }


TriangleMesh* loadPolyMeshFromFile(const char *file, const Matrix44f &o2w)
{
    std::ifstream ifs;
    try {
        ifs.open(file);
        if (ifs.fail()) throw;
        std::stringstream ss;
        ss << ifs.rdbuf();
        uint32_t numFaces;
        ss >> numFaces;
        std::unique_ptr<uint32_t []> faceIndex(new uint32_t[numFaces]);
        uint32_t vertsIndexArraySize = 0;
        // reading face index array
        for (uint32_t i = 0; i < numFaces; ++i) {
            ss >> faceIndex[i];
            vertsIndexArraySize += faceIndex[i];
        }
        std::unique_ptr<uint32_t []> vertsIndex(new uint32_t[vertsIndexArraySize]);
        uint32_t vertsArraySize = 0;
        // reading vertex index array
        for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
            ss >> vertsIndex[i];
            if (vertsIndex[i] > vertsArraySize) vertsArraySize = vertsIndex[i];
        }
        vertsArraySize += 1;
        // reading vertices
        std::unique_ptr<Vec3f []> verts(new Vec3f[vertsArraySize]);
        for (uint32_t i = 0; i < vertsArraySize; ++i) {
            ss >> verts[i].x >> verts[i].y >> verts[i].z;
        }
        // reading normals
        std::unique_ptr<Vec3f []> normals(new Vec3f[vertsIndexArraySize]);
        for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
            ss >> normals[i].x >> normals[i].y >> normals[i].z;
        }
        // reading st coordinates
        std::unique_ptr<Vec2f []> st(new Vec2f[vertsIndexArraySize]);
        for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
            ss >> st[i].x >> st[i].y;
        }

        return new TriangleMesh(o2w, numFaces, faceIndex, vertsIndex, verts, normals, st);
    }
    catch (...) {
        ifs.close();
        std::cout << "open file error!!!";
    }
    ifs.close();

    return nullptr;
}

bool trace(
    const Vec3f &orig, const Vec3f &dir,
    const std::vector<std::unique_ptr<Object>> &objects,
    IsectInfo &isect,
    RayType rayType = kPrimaryRay)
{
    isect.hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNear = kInfinity;
        uint32_t index = 0;
        Vec2f uv;
        if (objects[k]->intersect(orig, dir, tNear, index, uv) && tNear < isect.tNear) {
            isect.hitObject = objects[k].get();
            isect.tNear = tNear;
            isect.index = index;
            isect.uv = uv;
        }
    }

    return (isect.hitObject != nullptr);
}


void creatCoordinateSystem(const Vec3f &N, Vec3f &Nt, Vec3f &Nb)
{
    if(std::fabs(N.x) > std::fabs(N.y))
        Nt = Vec3f(N.z, 0, -N.x) / sqrtf(N.x * N.x + N.z * N.z);
    else
        Nt = Vec3f(0, -N.z, N.y) / sqrtf(N.y * N.y + N.z * N.z);
    Nb = N.crossProduct(Nt);
}

Vec3f uniformSampleHemisphere(const float &r1, const float &r2)
{
    float sinTheta = sqrtf(1 - r1 * r1);
    float phi = 2 * M_PI * r2;
    float x = sinTheta * cosf(phi);
    float z = sinTheta * sinf(phi);
    return Vec3f(x, r1, z);
}

std::default_random_engine generator;
std::uniform_real_distribution<float> distribution(0, 1);

Vec3f castRay(
    const Vec3f &orig, const Vec3f &dir,
    const std::vector<std::unique_ptr<Object>> &objects,
    const std::vector<std::unique_ptr<Light>> &lights,
    const Options &options,
    const uint32_t & depth = 0)
{
    if (depth > options.maxDepth) return options.backgroundColor;
    Vec3f hitColor = 0;
    IsectInfo isect;
    if (trace(orig, dir, objects, isect)) {
        Vec3f hitPoint = orig + dir * isect.tNear;
        Vec3f hitNormal;
        Vec2f hitTexCoordinates;
        isect.hitObject->getSurfaceProperties(hitPoint, dir, isect.index, isect.uv, hitNormal, hitTexCoordinates);
        switch (isect.hitObject->type) {
            case kDiffuse:
                {
                Vec3f directLighting = 0;
                for (uint32_t i = 0; i < lights.size(); ++i) {
                    Vec3f lightDir, lightIntensity;
                    IsectInfo isectShad;
                    lights[i]->illuminate(hitPoint, lightDir, lightIntensity, isectShad.tNear);
                    bool vis = !trace(hitPoint + hitNormal * options.bias, -lightDir, objects, isectShad, kShadowRay);
                    directLighting = vis * lightIntensity * std::max(0.f, hitNormal.dotProduct(-lightDir));
                }
            Vec3f indirectLigthing = 0;
#ifdef GI
            uint32_t N = 128;// / (depth + 1);
            Vec3f Nt, Nb;
            createCoordinateSystem(hitNormal, Nt, Nb);
            float pdf = 1 / (2 * M_PI);
            for (uint32_t n = 0; n < N; ++n) {
                float r1 = distribution(generator);
                float r2 = distribution(generator);
                Vec3f sample = uniformSampleHemisphere(r1, r2);
                Vec3f sampleWorld(
                    sample.x * Nb.x + sample.y * hitNormal.x + sample.z * Nt.x,
                    sample.x * Nb.y + sample.y * hitNormal.y + sample.z * Nt.y,
                    sample.x * Nb.z + sample.y * hitNormal.z + sample.z * Nt.z);
                // don't forget to divide by PDF and multiply by cos(theta)
                indirectLigthing += r1 * castRay(hitPoint + sampleWorld * options.bias,
                    sampleWorld, objects, lights, options, depth + 1) / pdf;
            }
            // divide by N
            indirectLigthing /= (float)N;
#endif

            hitColor = (directLighting / M_PI + 2 * indirectLigthing) * isect.hitObject->albedo;
            break;
        }
        default:
            break;
    }
}
else {
    hitColor = 1;
}

return hitColor;
}

//主渲染功能。
//我们迭代图像中的所有像素，生成主光线并将这些光线投射到场景中。
//帧缓冲区的内容保存到文件中。

void render(
    const Options &options,
    const std::vector<std::unique_ptr<Object>> &objects,
    const std::vector<std::unique_ptr<Light>> &lights)
{
    std::unique_ptr<Vec3f []> framebuffer(new Vec3f[options.width * options.height]);
    Vec3f *pix = framebuffer.get();
    float scale = tan(deg2rad(options.fov * 0.5));
    float imageAspectRatio = options.width / (float)options.height;
    Vec3f orig;
    options.cameraToWorld.multVecMatrix(Vec3f(0), orig);
    auto timeStart = std::chrono::high_resolution_clock::now();
    for (uint32_t j = 0; j < options.height; ++j) {
        for (uint32_t i = 0; i < options.width; ++i) {
            // generate primary ray direction
            float x = (2 * (i + 0.5) / (float)options.width - 1) * imageAspectRatio * scale;
            float y = (1 - 2 * (j + 0.5) / (float)options.height) * scale;

            if(j == 50 && i == 20)
                std::cout << "x:" << x << ";" << "y:" << y << ";" << std::endl;
            if(j == 301 && i == 20)
                std::cout << "x:" << x << ";" << "y:" << y << ";" << std::endl;
            Vec3f dir;
            options.cameraToWorld.multDirMatrix(Vec3f(x, y, -1), dir);
            dir.normalize();
            *(pix++) = castRay(orig, dir, objects, lights, options);
        }
        fprintf(stderr, "\r%3d%c", uint32_t(j / (float)options.height * 100), '%');
    }
    auto timeEnd = std::chrono::high_resolution_clock::now();
    auto passedTime = std::chrono::duration<double, std::milli>(timeEnd - timeStart).count();
    fprintf(stderr, "\rDone: %.2f (sec)\n", passedTime / 1000);

    // save framebuffer to file
    float gamma = 1;
    std::ofstream ofs;
    ofs.open("out.ppm");
    ofs << "P6\n" << options.width << " " << options.height << "\n255\n";
    for (uint32_t i = 0; i < options.height * options.width; ++i) {
        char r = (char)(255 * clamp(0, 1, powf(framebuffer[i].x, 1/gamma)));
        char g = (char)(255 * clamp(0, 1, powf(framebuffer[i].y, 1/gamma)));
        char b = (char)(255 * clamp(0, 1, powf(framebuffer[i].z, 1/gamma)));
        ofs << r << g << b;
    }
    ofs.close();
}

int main(int argc, char **argv)
{
    // loading gemetry
    std::vector<std::unique_ptr<Object>> objects;
    // lights
    std::vector<std::unique_ptr<Light>> lights;
    Options options;

    // aliasing example
    options.fov = 39.89;
    options.width = 512;
    options.height = 512;
    options.cameraToWorld = Matrix44f(0.965926, 0, -0.258819, 0, 0.0066019, 0.999675, 0.0246386, 0, 0.258735, -0.0255078, 0.965612, 0, 0.764985, 0.791882, 5.868275, 1);

    TriangleMesh *plane = loadPolyMeshFromFile("./planegi.geo", Matrix44f::kIdentity);
    if (plane != nullptr) {
        plane->albedo = Vec3f(0.225, 0.144, 0.144);
        objects.push_back(std::unique_ptr<Object>(plane));
    }

    TriangleMesh *cube = loadPolyMeshFromFile("./cubegi.geo", Matrix44f::kIdentity);
    if (cube != nullptr) {
        cube->albedo = Vec3f(0.188559, 0.287, 0.200726);
        objects.push_back(std::unique_ptr<Object>(cube));
    }

    Matrix44f xformSphere;
    xformSphere[3][1] = 1;
    Sphere *sph = new Sphere(xformSphere, 1);
    objects.push_back(std::unique_ptr<Object>(sph));

    Matrix44f l2w(0.916445, -0.218118, 0.335488, 0, 0.204618, -0.465058, -0.861309, 0, 0.343889, 0.857989, -0.381569, 0, 0, 0, 0, 1);
    lights.push_back(std::unique_ptr<Light>(new DistantLight(l2w, 1, 16)));

    // finally, render
    render(options, objects, lights);

    return 0;
}
