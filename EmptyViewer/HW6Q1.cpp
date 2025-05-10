#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#define NOMINMAX
#include <Windows.h>
#include <iostream>
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/freeglut.h>
#define GLFW_INCLUDE_GLU
#define GLFW_DLL
#include <GLFW/glfw3.h>

#include <vector>
#include <cmath>
#include <cstring>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int Width = 512, Height = 512;
constexpr int IMG_W = 512, IMG_H = 512;
std::vector<float> OutputImage;
static unsigned char framebuffer[IMG_H][IMG_W][3];
static float zbuffer[IMG_H][IMG_W];
int gNumVertices = 0, gNumTriangles = 0;
int* gIndexBuffer = nullptr;
float* gVertexBuffer = nullptr;


// ===== Vec3 Utilities =====
struct Vec3 {
    float x, y, z;
    Vec3() = default;
    Vec3(float X, float Y, float Z) : x(X), y(Y), z(Z) {}
    Vec3 operator+(const Vec3& o) const { return { x + o.x, y + o.y, z + o.z }; }
    Vec3 operator-(const Vec3& o) const { return { x - o.x, y - o.y, z - o.z }; }
    Vec3 operator*(float s) const { return { x * s, y * s, z * s }; }
};

inline float dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(Vec3 a, Vec3 b) {
    return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
}
inline Vec3 normalize(Vec3 v) {
    float L = std::sqrtf(dot(v, v));
    return (L > 0) ? v * (1.0f / L) : v;
}
inline float clamp01(float c) { return std::min(1.0f, std::max(0.0f, c)); }
inline unsigned char toSRGB(float c) {
    float gammaInv = 1.0f / 2.2f;
    return static_cast<unsigned char>(std::pow(clamp01(c), gammaInv) * 255.0f + 0.5f);
}


// ===== Fragment Processing (Flat shading per-triangle) =====
Vec3 shade(Vec3 Pcam, Vec3 Ncam) {
    Vec3 ka(0, 1, 0), kd(0, 0.5f, 0), ks(0.5f, 0.5f, 0.5f);
    float shininess = 32.f, Ia = 0.2f;
    Vec3 lightPos(-4, 4, -3);
    Vec3 L = normalize(lightPos - Pcam);
    Vec3 V = normalize(Pcam * -1.0f);
    Vec3 H = normalize(L + V); // Halfway vector

    float diff = std::max(0.f, dot(Ncam, L));
    float spec = std::pow(std::max(0.f, dot(Ncam, H)), shininess);

    return Vec3{
        ka.x * Ia + kd.x * diff + ks.x * spec,
        ka.y * Ia + kd.y * diff + ks.y * spec,
        ka.z * Ia + kd.z * diff + ks.z * spec
    };
}


struct VScreen {
    float x, y, depth;
    Vec3 cam;
};

// ===== Framebuffer Initialization =====
void clear_buffers() {
    std::memset(framebuffer, 0, sizeof(framebuffer));
    for (int y = 0; y < IMG_H; ++y)
        for (int x = 0; x < IMG_W; ++x)
            zbuffer[y][x] = 1e9f;
}

// ===== Fragment output =====
void put_pixel(int x, int y, float z, Vec3 color) {
    if ((unsigned)x >= IMG_W || (unsigned)y >= IMG_H) return;
    if (z < zbuffer[y][x]) {
        zbuffer[y][x] = z;
        framebuffer[y][x][0] = toSRGB(color.x);
        framebuffer[y][x][1] = toSRGB(color.y);
        framebuffer[y][x][2] = toSRGB(color.z);
    }
}

// ===== Rasterization: Pass-through color =====
void raster_triangle_flat(const VScreen& a, const VScreen& b, const VScreen& c, Vec3 col) {
    int x0 = int(a.x), y0 = int(a.y);
    int x1 = int(b.x), y1 = int(b.y);
    int x2 = int(c.x), y2 = int(c.y);
    int minX = std::max(0, std::min({ x0, x1, x2 }));
    int maxX = std::min(IMG_W - 1, std::max({ x0, x1, x2 }));
    int minY = std::max(0, std::min({ y0, y1, y2 }));
    int maxY = std::min(IMG_H - 1, std::max({ y0, y1, y2 }));
    float area = float((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
    if (area == 0) return;
    float invA = 1.0f / area;
    for (int y = minY; y <= maxY; ++y)
        for (int x = minX; x <= maxX; ++x) {
            float w0 = float((x1 - x0) * (y - y0) - (x - x0) * (y1 - y0));
            float w1 = float((x2 - x1) * (y - y1) - (x - x1) * (y2 - y1));
            float w2 = float((x0 - x2) * (y - y2) - (x - x2) * (y0 - y2));
            if (!((w0 >= 0 && w1 >= 0 && w2 >= 0) || (w0 <= 0 && w1 <= 0 && w2 <= 0))) continue;
            float aB = w1 * invA, bB = w2 * invA, cB = 1 - aB - bB;
            float depth = aB * a.depth + bB * b.depth + cB * c.depth;
            put_pixel(x, y, depth, col);
        }
}

// ===== Input: Geometry generation =====
void create_scene() {
        int width = 32, height = 16;
        gNumVertices = (height - 2) * width + 2;
        gNumTriangles = (height - 2) * (width - 1) * 2;
        gVertexBuffer = new float[3 * gNumVertices];
        gIndexBuffer = new int[3 * gNumTriangles];
        int t = 0;
        for (int j = 1; j < height - 1; ++j)
            for (int i = 0; i < width; ++i) {
                float theta = float(j) / (height - 1) * M_PI;
                float phi = float(i) / (width - 1) * M_PI * 2.f;
                gVertexBuffer[3 * t + 0] = sinf(theta) * cosf(phi);
                gVertexBuffer[3 * t + 1] = cosf(theta);
                gVertexBuffer[3 * t + 2] = -sinf(theta) * sinf(phi);
                ++t;
            }
        gVertexBuffer[3 * t + 0] = 0; gVertexBuffer[3 * t + 1] = 1; gVertexBuffer[3 * t + 2] = 0; ++t;
        gVertexBuffer[3 * t + 0] = 0; gVertexBuffer[3 * t + 1] = -1; gVertexBuffer[3 * t + 2] = 0; ++t;
        t = 0;
        for (int j = 0; j < height - 3; ++j)
            for (int i = 0; i < width - 1; ++i) {
                gIndexBuffer[t++] = j * width + i;
                gIndexBuffer[t++] = (j + 1) * width + i + 1;
                gIndexBuffer[t++] = j * width + i + 1;
                gIndexBuffer[t++] = j * width + i;
                gIndexBuffer[t++] = (j + 1) * width + i;
                gIndexBuffer[t++] = (j + 1) * width + i + 1;
            }
        for (int i = 0; i < width - 1; ++i) {
            gIndexBuffer[t++] = (height - 2) * width;
            gIndexBuffer[t++] = i;
            gIndexBuffer[t++] = i + 1;
            gIndexBuffer[t++] = (height - 2) * width + 1;
            gIndexBuffer[t++] = (height - 3) * width + i + 1;
            gIndexBuffer[t++] = (height - 3) * width + i;
        }
    }

// ===== Vertex Processing + Rasterization =====
void render() {
    clear_buffers();
    create_scene();
    std::vector<VScreen> vs(gNumVertices);
    float l = -0.1f, r = 0.1f, b = -0.1f, t = 0.1f, n = -0.1f, f = -1000.f;

    float persp_mat[4][4] = {
        { 2 * n / (r - l), 0, (r + l) / (r - l), 0 },
        { 0, 2 * n / (t - b), (t + b) / (t - b), 0 },
        { 0, 0, -(f + n) / (f - n), -2 * f * n / (f - n) },
        { 0, 0, -1, 0 }
    };

    for (int i = 0; i < gNumVertices; ++i) {
        Vec3 obj{ gVertexBuffer[3 * i + 0], gVertexBuffer[3 * i + 1], gVertexBuffer[3 * i + 2] };
        Vec3 cam = obj * 2.f + Vec3{ 0, 0, -7.f }; // modeling transform

        float x = cam.x, y = cam.y, z = cam.z;
        float px = persp_mat[0][0] * x + persp_mat[0][2] * z;
        float py = persp_mat[1][1] * y + persp_mat[1][2] * z;
        float pw = -z;

        float sx = (1 - px / pw) * 0.5f * IMG_W;
        float sy = (1 - py / pw) * 0.5f * IMG_H;


        vs[i] = { sx, sy, pw, cam };
    }
    for (int t = 0; t < gNumTriangles; ++t) {
        VScreen& A = vs[gIndexBuffer[3 * t + 0]];
        VScreen& B = vs[gIndexBuffer[3 * t + 1]];
        VScreen& C = vs[gIndexBuffer[3 * t + 2]];
        Vec3 P0 = A.cam, P1 = B.cam, P2 = C.cam;
        Vec3 N = normalize(cross(P1 - P0, P2 - P0));
        Vec3 Pc = (P0 + P1 + P2) * (1.f / 3.f);
        if (dot(N, Pc * -1.0f) < 0) N = N * -1.f;
        Vec3 col = shade(Pc, N);
        raster_triangle_flat(A, B, C, col);
    }
    OutputImage.resize(IMG_W * IMG_H * 3);
    for (int j = 0; j < IMG_H; ++j)
        for (int i = 0; i < IMG_W; ++i) {
            int idx = (j * IMG_W + i) * 3;
            OutputImage[idx + 0] = framebuffer[j][i][0] / 255.0f;
            OutputImage[idx + 1] = framebuffer[j][i][1] / 255.0f;
            OutputImage[idx + 2] = framebuffer[j][i][2] / 255.0f;
        }
}

void resize_callback(GLFWwindow*, int nw, int nh) {
    Width = nw;
    Height = nh;
    glViewport(0, 0, nw, nh);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, Width, 0.0, Height, 1.0, -1.0);
    OutputImage.reserve(Width * Height * 3);
    render();
}

int main() {
    if (!glfwInit()) return -1;
    GLFWwindow* window = glfwCreateWindow(Width, Height, "Flat Shaded Sphere", NULL, NULL);
    if (!window) { glfwTerminate(); return -1; }
    glfwMakeContextCurrent(window);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glfwSetFramebufferSizeCallback(window, resize_callback);
    resize_callback(NULL, Width, Height);
    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);
        glRasterPos2i((Width - IMG_W) / 2, (Height - IMG_H) / 2);
        glDrawPixels(IMG_W, IMG_H, GL_RGB, GL_FLOAT, OutputImage.data());
        glfwSwapBuffers(window);
        glfwPollEvents();
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS ||
            glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
            glfwSetWindowShouldClose(window, true);
    }
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
