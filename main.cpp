#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <chrono>
#include <cstdint>

#include <random>

#include "geometry.h"
#include "lights.h"
#include "objects.h"
#include "bezier_data.h"

uint32_t row = 50;
uint32_t column = 300;

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


Vec3f evalBezierCurve(const Vec3f *P, const float u)
{
    float K0 = (1 - u) * (1 - u) * (1 - u);
    float K1 = 3 * u * (1 -u) * (1 - u);
    float K2 = 3 * u * u * (1 - u);
    float K3 = u * u * u;
    return Vec3f(P[0] * K0 + P[1] * K1 + P[2] * K2 + P[3] *K3);
}

Vec3f evalBezierSurface(const Vec3f *P, const float u, const float v)
{
    Vec3f uCurve[4];
    for(int i = 0; i < 4; ++i)
        uCurve[i] = evalBezierCurve(P + i * 4, u);
    return evalBezierCurve(uCurve, v);
}

void evalBezierCurveFFD(const uint32_t &divs, const Vec3f &P0, const Vec3f &P1, const Vec3f &P2, const Vec3f &P3, Vec3f *B)
{
#if 1
    float h = 1 / (float)divs;
    Vec3f b0 = P0;
    Vec3f fph = 3 * (P1 - P0) * h;
    Vec3f fpphh = (6 * P0 - 12 * P1 + 6 * P2) * h * h;
    Vec3f fppphhh = (-6 * P0 + 18 * P1 - 18 * P2 + 6 * P3) * h * h * h;
    B[0] = b0;
    for(uint32_t i = 1; i <= divs; ++i) {
        B[i] = b0 + fph + fpphh + fppphhh / 6;
        fph = fph + fpphh + fppphhh / 2;
        fpphh = fpphh + fppphhh;
    }
#else
    Vec3f b0 = P0;
    Vec3f bd0 = 3 * (P1 - P0);
    Vec3f bdd0 = (6 * P0 - 12 * P1 + 6 * P2);
    Vec3f bddd0 = (-6 * P0 + 18 * P1 - 18 * P2 + 6 * P3);
    for(uint32_t i = 0; i < divs; ++i) {
        float x = i / (float)divs;
        B[i] = b0 + bd0 * x + bdd0 * x / 2 + bddd0 * x / 6;
    }
#endif
}

void evalBezierPatchFFD(const uint32_t &divs, const Vec3f *controlPoints, Vec3f *&P)
{
    Vec3f controlPointsV[4][divs + 1];
    for(uint32_t i = 0; i < 4; ++i) {
        evalBezierCurveFFD(divs, controlPoints[i], controlPoints[i + 4],
                controlPoints[i + 8], controlPoints[i + 12], controlPointsV[i]);
    }
    for(uint32_t i = 0; i <= divs; ++i) {
        evalBezierCurveFFD(divs, controlPointsV[0][i], controlPointsV[1][i],
                controlPointsV[2][i], controlPointsV[3][i], P + i * (divs + 1));
    }
}


Vec3f derivBezier(const Vec3f *P, const float &t)
{
    return -3 * (1 - t) * (1 - t) * P[0] +
        (3 * (1 - t) * (1 - t) - 6 * t * (1 - t)) * P[1] +
        (6 * t * (1 - t) - 3 * t * t) * P[2] +
        3 * t * t * P[3];
}


Vec3f dUBezier(const Vec3f *controlPoints, const float &u, const float &v)
{
    Vec3f P[4];
    Vec3f vCurve[4];
    for(int i = 0; i < 4; ++i) {
        P[0] = controlPoints[i];
        P[1] = controlPoints[4 + i];
        P[2] = controlPoints[8 + i];
        P[3] = controlPoints[12 + i];
        vCurve[i] = evalBezierCurve(P, v);
    }

    return derivBezier(vCurve, u);
}

Vec3f dVBezier(const Vec3f *controlPoints, const float &u, const float &v)
{
    Vec3f uCurve[4];
    for(int i = 0; i < 4; ++i) {
        uCurve[i] = evalBezierCurve(controlPoints + 4 * i, u);
    }
    return derivBezier(uCurve, v);
}

void createCurveGeometry(std::vector<std::unique_ptr<Object>> &object)
{
    uint32_t ndivs = 16;
    uint32_t ncurves = 1 + (curveNumPts - 4) / 3;
    Vec3f pts[4];
    std::unique_ptr<Vec3f []> P(new Vec3f[(ndivs + 1) * ndivs * ncurves + 1]);//1904
    std::unique_ptr<Vec3f[]> N(new Vec3f[(ndivs + 1) * ndivs * ncurves + 1]);
    std::unique_ptr<Vec2f[]> st(new Vec2f[(ndivs + 1) * ndivs * ncurves + 1]);
    for(uint32_t i = 0; i < ncurves; ++i) {
        for(uint32_t j = 0; j < ndivs; ++j) {
            pts[0] = curveData[i * 3];
            pts[1] = curveData[i * 3 + 1];
            pts[2] = curveData[i * 3 + 2];
            pts[3] = curveData[i * 3 + 3];
            float s = j / (float)ndivs;
            Vec3f pt = evalBezierCurve(pts, s);
            Vec3f tangent = derivBezier(pts, s).normalize();
            bool swap = false;

            uint8_t maxAxis;
            if(std::abs(tangent.x) > std::abs(tangent.y))
                if(std::abs(tangent.x) > std::abs(tangent.z))
                    maxAxis = 0;
                else
                    maxAxis = 2;
            else if(std::abs(tangent.y) > std::abs(tangent.z))
                maxAxis = 1;
            else
                maxAxis = 2;

            Vec3f up, forward, right;

            switch(maxAxis) {
            case 0:
            case 1:
                up = tangent;
                forward = Vec3f(0, 0, 1);
                right = up.crossProduct(forward);
                forward = right.crossProduct(up);
                break;
            case 2:
                up = tangent;
                right = Vec3f(0, 0, 1);
                forward = right.crossProduct(up);
                right = up.crossProduct(forward);
            default:
                break;
            };

            float sNormalized = (i * ndivs + j) / float(ndivs * ncurves);
            float rad = 0.1 *(1 - sNormalized);
            for(uint32_t k = 0; k <= ndivs; ++k) {
                float t = k / (float)ndivs;
                float theta = t * 2 * M_PI;
                Vec3f pc(cos(theta) * rad, 0, sin(theta) * rad);
                float x = pc.x * right.x + pc.y * up.x + pc.z * forward.x;
                float y = pc.x * right.y + pc.y * up.y + pc.z * forward.y;
                float z = pc.x * right.z + pc.y * up.z + pc.z * forward.z;
                P[i * (ndivs + 1) * ndivs + j * (ndivs + 1) + k] = Vec3f(pt.x + x, pt.y + y, pt.z + z);
                N[i * (ndivs + 1) * ndivs + j * (ndivs + 1) + k] = Vec3f(x, y, z).normalize();
                st[i * (ndivs + 1) * ndivs + j * (ndivs + 1) + k] = Vec2f(sNormalized, t);
            }
        }
    }
    P[(ndivs + 1) * ndivs * ncurves] = curveData[curveNumPts - 1];
    N[(ndivs + 1) * ndivs * ncurves] = (curveData[curveNumPts - 2] - curveData[curveNumPts - 1]).normalize();
    st[(ndivs + 1) * ndivs * ncurves] = Vec2f(1, 0.5);
    uint32_t numFace = ndivs * ndivs * ncurves;
    std::unique_ptr<uint32_t []> verts(new uint32_t[numFace]);
    for(uint32_t i = 0; i < numFace; ++i)
        verts[i] = (i < (numFace - ndivs)) ? 4 : 3;
    std::unique_ptr<uint32_t []> vertIndices(new uint32_t[ndivs * ndivs * ncurves * 4 + ndivs * 3]);
    uint32_t nf = 0, ix = 0;
    for(uint32_t k = 0; k < ncurves; ++k) {
        for(uint32_t j = 0; j < ndivs; ++j) {
            if(k == (ncurves - 1) && j == (ndivs - 1)) {break; }
            for(uint32_t i = 0; i < ndivs; ++i) {
                vertIndices[ix] = nf;
                vertIndices[ix + 1] = nf +(ndivs + 1);
                vertIndices[ix + 2] = nf +(ndivs + 1) + 1;
                vertIndices[ix + 3] = nf + 1;
                ix += 4;
                ++nf;
            }
            nf++;
        }
    }

    for(uint32_t i = 0; i < ndivs; ++i) {
        vertIndices[ix] = nf;
        vertIndices[ix + 1] = (ndivs + 1) * ndivs * ncurves;
        vertIndices[ix + 2] = nf + 1;
        ix += 3;
        nf++;
    }

    object.push_back(std::unique_ptr<TriangleMesh>(new TriangleMesh(Matrix44f::kIdentity, numFace, verts, vertIndices, P, N, st)));
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
                    //std::cout << "call;" << std::endl;
                    lights[i]->illuminate(hitPoint, lightDir, lightIntensity, isectShad.tNear);
                    bool vis = !trace(hitPoint + hitNormal * options.bias, -lightDir, objects, isectShad, kShadowRay);
                    directLighting = vis * lightIntensity * std::max(0.f, hitNormal.dotProduct(-lightDir));
                }
            Vec3f indirectLigthing = 0;
#ifdef GI
            uint32_t N = 128;// / (depth + 1);
            Vec3f Nt, Nb;
            creatCoordinateSystem(hitNormal, Nt, Nb);
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
    ofs << "P6\n" << options.height << " " << options.width << "\n255\n";
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
    std::cout << "call" << std::endl;

    createCurveGeometry(objects);

    // lights
    std::vector<std::unique_ptr<Light>> lights;
    Options options;
    // aliasing example
    options.fov = 39.89;
    options.width = 512;
    options.height = 512;
    options.maxDepth = 1;
    // to render the teapot
    //options.cameraToWorld = Matrix44f(0.897258, 0, -0.441506, 0, -0.288129, 0.757698, -0.585556, 0, 0.334528, 0.652606, 0.679851, 0, 5.439442, 11.080794, 10.381341, 1);

    // to render the curve as geometry
    options.cameraToWorld = Matrix44f(0.707107, 0, -0.707107, 0, -0.369866, 0.85229, -0.369866, 0, 0.60266, 0.523069, 0.60266, 0, 2.634, 3.178036, 2.262122, 1);

    Matrix44f l2w(0.916445, -0.218118, 0.335488, 0, 0.204618, -0.465058, -0.861309, 0, 0.343889, 0.857989, -0.381569, 0, 0, 0, 0, 1);
    lights.push_back(std::unique_ptr<Light>(new DistantLight(l2w, 1, 16)));

    // finally, render
    render(options, objects, lights);

    return 0;
}
