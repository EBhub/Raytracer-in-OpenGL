#include "bounding_volume_hierarchy.h"
#include "disable_all_warnings.h"
#include "draw.h"
#include "image.h"
#include "ray_tracing.h"
#include "screen.h"
#include "trackball.h"
#include "window.h"
// Disable compiler warnings in third-party code (which we cannot change).
DISABLE_WARNINGS_PUSH()
#include <glm/gtc/type_ptr.hpp>
#include <glm/vec2.hpp>
#include <glm/vec4.hpp>
#include <imgui.h>
DISABLE_WARNINGS_POP()
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <optional>
#include <random>
#include <string>
#include <type_traits>
#ifdef USE_OPENMP
#include <omp.h>
#endif
#define MAX_LEVEL 3

// This is the main application. The code in here does not need to be modified.
constexpr glm::ivec2 windowResolution{ 800, 800 };
const std::filesystem::path dataPath{ DATA_DIR };
const std::filesystem::path outputPath{ OUTPUT_DIR };

enum class ViewMode {
    Rasterization = 0,
    RayTracing = 1
};

static glm::vec3 shade(int level, BoundingVolumeHierarchy bvh, Ray ray, HitInfo info, const Scene& scene);
static glm::vec3 getDiffuse(HitInfo hitInfo, float nL, glm::vec3 light);
static glm::vec3 getSpecular(HitInfo hitInfo, glm::vec3 L, glm::vec3 V);
static glm::vec3 getReflected(int level, Scene scene, BoundingVolumeHierarchy bvh, glm::vec3 hit, Ray ray, HitInfo hitInfo);
static glm::vec3 getSoftShadows(SphericalLight light, glm::vec3 hit, BoundingVolumeHierarchy bvh, HitInfo hitInfo);
static glm::vec3 getRefracted(Ray ray, HitInfo info);

// NOTE(Mathijs): separate function to make recursion easier (could also be done with lambda + std::function).
static glm::vec3 getFinalColor(int level, const Scene& scene, const BoundingVolumeHierarchy& bvh, Ray ray) {
    HitInfo hitInfo;

    if (bvh.intersect(ray, hitInfo) && level < MAX_LEVEL) {
        level++;

        glm::vec3 color = shade(level, bvh, ray, hitInfo, scene);

        glm::vec3 reflected = {0.0f, 0.0f, 0.0f};

        glm::vec3 hit = ray.origin + ray.t * ray.direction;
        if (hitInfo.material.ks != glm::vec3{ 0.f, 0.f, 0.f } && level < MAX_LEVEL) {
            reflected = getReflected(level, scene, bvh, hit, ray, hitInfo);

        }
    /**
        // Set default to black
        glm::vec3 refractionColor = {0.0f, 0.0f, 0.0f};
        if (hitInfo.material.indexRefraction > 1.0f && level < MAX_LEVEL) {
            glm::vec3 N = hitInfo.normal;
            float ni = hitInfo.material.indexRefraction;
            float ratio = 1.0f;
            glm::vec3 refracted = getRefracted(ray, hitInfo);
            //glm::vec3 reflected = getReflected(level, scene, bvh, hit, ray, hitInfo);

            Ray refract = Ray{hit + refracted * 0.00001f, glm::normalize(refracted), ray.t};
            glm::vec3 refractionColor = getFinalColor(level, scene, bvh, refract);

        }

        glm::vec3 transparency = {0.0f, 0.0f, 0.0f};
        if (hitInfo.material.transparency < 1.0f) {
            // Continue the ray to the next object
            // level-- since we continued, not reflect/refract
            // return transparency * color
            float t = hitInfo.material.transparency;
            Ray onwards = Ray {hit + 0.0001f * ray.direction, ray.direction, ray.t};

            if (bvh.intersect(onwards, hitInfo)) {
                transparency =  (1.0f - t) * getFinalColor((level - 1), scene, bvh, onwards);
            }
        }*/

        color = color + reflected;

        // Draw the ray to be the color at the hit point
        drawRay(ray, color);
        // Set the color of the pixel to white if the ray hits.
        return color;
    } else {
        // Draw a red debug ray if the ray missed.
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        // Set the color of the pixel to black if the ray misses.
        return glm::vec3(0.0f);
    }
}

static glm::vec3 shade(int level, BoundingVolumeHierarchy bvh, Ray ray, HitInfo info, const Scene& scene) {
    glm::vec3 hit = ray.origin + ray.t * ray.direction;
    glm::vec3 color = glm::vec3(0);
    glm::vec3 N = glm::normalize(info.normal);

    // Go through all of the point lights in the scene
    for (PointLight lights : scene.pointLights) {
        glm::vec3 L = glm::normalize(lights.position - hit);
        glm::vec3 V = glm::normalize(ray.origin - hit);
        // Check if it's actually lit
        if (glm::dot(N, L) * glm::dot(N, V) <= 0.f) {
            continue;
        }

        // Check for hard shadows from point lights
        glm::vec3 shadowDirection = glm::normalize(lights.position - hit);
        float distToLight = glm::distance(lights.position, hit);
        Ray shadowRay = { hit + shadowDirection * 1e-5f, shadowDirection, distToLight };

        // If the shadow ray hits something before the light, draw red ray
        if (bvh.intersect(shadowRay, info) && shadowRay.t < distToLight) {
            drawRay(shadowRay, glm::vec3(1.0f, 0.0f, 0.0f));
            continue;
        } else {
            // Draw a white ray
            drawRay(shadowRay, glm::vec3(1.0f, 1.0f, 1.0f));
        }

        // Diffuse shading
        glm::vec3 diffuse = getDiffuse(info, glm::max(glm::dot(N, L), 0.0f), lights.color);

        // Specular shading
        glm::vec3 specular = getSpecular(info, L, V);

        color += (diffuse + specular);
    }

    // Go through spherical lights, check for soft shadows, and calculate shading appropriately
    for (SphericalLight light : scene.sphericalLight) {
        glm::vec3 L = glm::normalize(light.position - hit);
        glm::vec3 V = glm::normalize(ray.origin - hit);

        // Check if it's actually lit
        if (glm::dot(N, L) * glm::dot(N, V) <= 0.f) {
            continue;
        }

        // Get diffuse shading
        glm::vec3 diffuse = getDiffuse(info, glm::max(glm::dot(N, L), 0.0f), light.color);

        // Get specular shading
        glm::vec3 specular = getSpecular(info, L, V);

        // Get intensity of shadow
        glm::vec3 intensity = getSoftShadows(light, hit, bvh, info);

        // Put it all together
        color += intensity * (diffuse + specular);
    }

    return color;
}

// Using Snell's law = ni * sin(thetai) = nt * sin(thetat)
// https://en.wikipedia.org/wiki/Snell%27s_law
// Found t using lecture Ray Tracing, slide on Refraction
static glm::vec3 getRefracted(Ray ray, HitInfo info) {
    // Incidence vector
    glm::vec3 I = ray.direction;
    glm::vec3 N = info.normal;
    float theta = glm::dot(I, N);
    float ni, nt = 1.0f;

    if (theta > 0.0001f) {
        nt = info.material.indexRefraction;
    } else {
        ni = info.material.indexRefraction;
        N = -N;
    }

    float n = ni / nt;
    float cosi = glm::cos(glm::dot(I, N));
    float cosphi = 1.0f - ((n * n) * (1.0f - cosi * cosi));


    float num = 1.0f - ((ni * ni) * (1.0f - glm::pow(glm::dot(ray.direction, info.normal), 2.0f) / (nt * nt)));

    // No refraction, but total reflection
    if (num < 0.00001f) {
        return glm::reflect(I, N);
    }

    glm::vec3 t = n * (ray.direction - (glm::dot(ray.direction, info.normal) * info.normal)) - info.normal * glm::sqrt(num);

    return t;
}

static glm::vec3 getDiffuse(HitInfo hitInfo, float nL, glm::vec3 light) {
    glm::vec3 kd = hitInfo.material.kd;

    return kd * nL * light;
}

static glm::vec3 getSpecular(HitInfo hitInfo, glm::vec3 L, glm::vec3 V) {
    glm::vec3 N = glm::normalize(hitInfo.normal);
    glm::vec3 ks = hitInfo.material.ks;
    float shininess = hitInfo.material.shininess;

    glm::vec3 reflected = glm::reflect(-L, N);
    float spec = glm::pow(glm::max(glm::dot(V, reflected), 0.0f), shininess);

    return ks * spec;
}

static glm::vec3 getSoftShadows(SphericalLight light, glm::vec3 hit, BoundingVolumeHierarchy bvh, HitInfo hitInfo) {
    // Variable amount of sampling
    int samples = 100;
    int total = 0;
    glm::vec3 base = glm::vec3(0);
    float r = light.radius;

    // After research, using the fibonacci lattice to project points seems to give a more equal distribution
    // Further supported by various stackoverflow posts and research paper listed below
    // Ex: https://math.stackexchange.com/questions/1358046/is-the-fibonacci-lattice-the-very-best-way-to-evenly-distribute-n-points-on-a-sp
    float phi = (1.0f + glm::sqrt(5.0f)) / 2.0f; //https://en.wikipedia.org/wiki/Golden_ratio

    // Calculate sample number of shadow rays to cast towards spherical light
    for (int i = 0; i < samples; i++) {
        // Latitude and longitude are taken from ideas here: https://www.researchgate.net/publication/45891871_Measurement_of_Areas_on_a_Sphere_Using_Fibonacci_and_Latitude-Longitude_Lattices
        float latitude = glm::asin(2.0f * float(i) / (2.0f * samples + 1.0f));
        float longitude = glm::pow(phi, -2.0f) * 2.0f * glm::pi<float>() * i;

        // Using
        float x = r * std::cos(longitude) * std::cos(latitude);
        float y = r * std::sin(longitude) * std::cos(latitude);
        float z = r * std::sin(latitude);

        base = {x, y, z};
        // Translate to light position
        base = glm::normalize((base + light.position) - hit);
        float distToLight = glm::distance(light.position, hit);
        Ray shadowRay = {hit + base * 1e-5f, base, distToLight};

        if (bvh.intersect(shadowRay, hitInfo)) {
            // Draw red line
            drawRay(shadowRay, glm::vec3(1.0f, 0, 0));
            total++;
        } else {
            // If nothing is in the way, draw white ray
            // Convert sphere light to "sphere"
            Sphere l = Sphere {light.position, light.radius + (0.05f * light.radius), light.color};
            bool intersect = intersectRayWithShape(l, shadowRay, hitInfo);
            drawRay(shadowRay, glm::vec3(1.0f, 1.0f, 1.0f));
        }
    }

    float obstructed = (float)total / (float)samples;
    float intensity = 1.0f - obstructed;

    return intensity * light.color;
}

static glm::vec3 getReflected(int level, Scene scene, BoundingVolumeHierarchy bvh, glm::vec3 hit, Ray ray, HitInfo hitInfo) {
    glm::vec3 R = glm::normalize(glm::reflect(ray.direction, hitInfo.normal));
    Ray reflected = {hit + 1e-4f * R, R, ray.t};

    return getFinalColor(level, scene, bvh, reflected);
}


static void setOpenGLMatrices(const Trackball& camera);
static void renderOpenGL(const Scene& scene, const Trackball& camera, int selectedLight);

// This is the main rendering function. You are free to change this function in any way (including the function signature).
static void renderRayTracing(const Scene& scene, const Trackball& camera, const BoundingVolumeHierarchy& bvh, Screen& screen) {
    int level = 0;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos{
                    float(x) / windowResolution.x * 2.0f - 1.0f,
                    float(y) / windowResolution.y * 2.0f - 1.0f
            };
            const Ray cameraRay = camera.generateRay(normalizedPixelPos);
            screen.setPixel(x, y, getFinalColor(level, scene, bvh, cameraRay));
        }
    }
}

int main(int argc, char** argv) {
    Trackball::printHelp();
    std::cout << "\n Press the [R] key on your keyboard to create a ray towards the mouse cursor" << std::endl
              << std::endl;

    Window window{ "Final Project - Part 2", windowResolution, OpenGLVersion::GL2 };
    Screen screen{ windowResolution };
    Trackball camera{ &window, glm::radians(50.0f), 3.0f };
    camera.setCamera(glm::vec3(0.0f, 0.0f, 0.0f), glm::radians(glm::vec3(20.0f, 20.0f, 0.0f)), 3.0f);

    SceneType sceneType{ SceneType::SingleTriangle };
    std::optional<Ray> optDebugRay;
    Scene scene = loadScene(sceneType, dataPath);
    BoundingVolumeHierarchy bvh{ &scene };

    int bvhDebugLevel = 0;
    bool debugBVH{ false };
    ViewMode viewMode{ ViewMode::Rasterization };

    window.registerKeyCallback([&](int key, int /* scancode */, int action, int /* mods */) {
        if (action == GLFW_PRESS) {
            switch (key) {
                case GLFW_KEY_R: {
                    // Shoot a ray. Produce a ray from camera to the far plane.
                    const auto tmp = window.getNormalizedCursorPos();
                    optDebugRay = camera.generateRay(tmp * 2.0f - 1.0f);
                    viewMode = ViewMode::Rasterization;
                } break;
                case GLFW_KEY_ESCAPE: {
                    window.close();
                } break;
            };
        }
    });

    int selectedLight{ 0 };
    while (!window.shouldClose()) {
        window.updateInput();

        // === Setup the UI ===
        ImGui::Begin("Final Project - Part 2");
        {
            constexpr std::array items{ "SingleTriangle", "Cube", "Cornell Box (with mirror)", "Cornell Box (spherical light and mirror)", "Monkey", "Dragon", /* "AABBs",*/ "Spheres", /*"Mixed",*/ "Custom" };
            if (ImGui::Combo("Scenes", reinterpret_cast<int*>(&sceneType), items.data(), int(items.size()))) {
                optDebugRay.reset();
                scene = loadScene(sceneType, dataPath);
                bvh = BoundingVolumeHierarchy(&scene);
                if (optDebugRay) {
                    HitInfo dummy{};
                    bvh.intersect(*optDebugRay, dummy);
                }
            }
        }
        {
            constexpr std::array items{ "Rasterization", "Ray Traced" };
            ImGui::Combo("View mode", reinterpret_cast<int*>(&viewMode), items.data(), int(items.size()));
        }
        if (ImGui::Button("Render to file")) {
            {
                using clock = std::chrono::high_resolution_clock;
                const auto start = clock::now();
                renderRayTracing(scene, camera, bvh, screen);
                const auto end = clock::now();
                std::cout << "Time to render image: " << std::chrono::duration<float, std::milli>(end - start).count() << " milliseconds" << std::endl;
            }
            screen.writeBitmapToFile(outputPath / "render.bmp");
        }
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Text("Debugging");
        if (viewMode == ViewMode::Rasterization) {
            ImGui::Checkbox("Draw BVH", &debugBVH);
            if (debugBVH)
                ImGui::SliderInt("BVH Level", &bvhDebugLevel, 0, bvh.numLevels() - 1);
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Text("Lights");
        if (!scene.pointLights.empty() || !scene.sphericalLight.empty()) {
            {
                std::vector<std::string> options;
                for (size_t i = 0; i < scene.pointLights.size(); i++) {
                    options.push_back("Point Light " + std::to_string(i + 1));
                }
                for (size_t i = 0; i < scene.sphericalLight.size(); i++) {
                    options.push_back("Spherical Light " + std::to_string(i + 1));
                }

                std::vector<const char*> optionsPointers;
                std::transform(std::begin(options), std::end(options), std::back_inserter(optionsPointers),
                               [](const auto& str) { return str.c_str(); });

                ImGui::Combo("Selected light", &selectedLight, optionsPointers.data(), static_cast<int>(optionsPointers.size()));
            }

            {
                const auto showLightOptions = [](auto& light) {
                    ImGui::DragFloat3("Light position", glm::value_ptr(light.position), 0.01f, -3.0f, 3.0f);
                    ImGui::ColorEdit3("Light color", glm::value_ptr(light.color));
                    if constexpr (std::is_same_v<std::decay_t<decltype(light)>, SphericalLight>) {
                        ImGui::DragFloat("Light radius", &light.radius, 0.01f, 0.01f, 0.5f);
                    }
                };
                if (selectedLight < static_cast<int>(scene.pointLights.size())) {
                    // Draw a big yellow sphere and then the small light sphere on top.
                    showLightOptions(scene.pointLights[selectedLight]);
                }
                else {
                    // Draw a big yellow sphere and then the smaller light sphere on top.
                    showLightOptions(scene.sphericalLight[selectedLight - scene.pointLights.size()]);
                }
            }
        }

        if (ImGui::Button("Add point light")) {
            scene.pointLights.push_back(PointLight{ glm::vec3(0.0f), glm::vec3(1.0f) });
            selectedLight = int(scene.pointLights.size() - 1);
        }
        if (ImGui::Button("Add spherical light")) {
            scene.sphericalLight.push_back(SphericalLight{ glm::vec3(0.0f), 0.1f, glm::vec3(1.0f) });
            selectedLight = int(scene.pointLights.size() + scene.sphericalLight.size() - 1);
        }
        if (ImGui::Button("Remove selected light")) {
            if (selectedLight < static_cast<int>(scene.pointLights.size())) {
                scene.pointLights.erase(std::begin(scene.pointLights) + selectedLight);
            }
            else {
                scene.sphericalLight.erase(std::begin(scene.sphericalLight) + (selectedLight - scene.pointLights.size()));
            }
            selectedLight = 0;
        }

        // Clear screen.
        glClearDepth(1.0f);
        glClearColor(0.0, 0.0, 0.0, 0.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Draw either using OpenGL (rasterization) or the ray tracing function.
        switch (viewMode) {
            case ViewMode::Rasterization: {
                glPushAttrib(GL_ALL_ATTRIB_BITS);
                renderOpenGL(scene, camera, selectedLight);
                if (optDebugRay) {
                    // Call getFinalColor for the debug ray. Ignore the result but tell the function that it should
                    // draw the rays instead.
                    enableDrawRay = true;
                    int level = 0;
                    (void)getFinalColor(level, scene, bvh, *optDebugRay);
                    enableDrawRay = false;
                }
                glPopAttrib();
            } break;
            case ViewMode::RayTracing: {
                screen.clear(glm::vec3(0.0f));
                renderRayTracing(scene, camera, bvh, screen);
                screen.setPixel(0, 0, glm::vec3(1.0f));
                screen.draw(); // Takes the image generated using ray tracing and outputs it to the screen using OpenGL.
            } break;
            default:
                break;
        };

        if (debugBVH) {
            glPushAttrib(GL_ALL_ATTRIB_BITS);
            setOpenGLMatrices(camera);
            glDisable(GL_LIGHTING);
            glEnable(GL_DEPTH_TEST);

            // Enable alpha blending. More info at:
            // https://learnopengl.com/Advanced-OpenGL/Blending
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            bvh.debugDraw(bvhDebugLevel);
            glPopAttrib();
        }

        ImGui::End();
        window.swapBuffers();
    }

    return 0; // execution never reaches this point
}

static void setOpenGLMatrices(const Trackball& camera) {
    // Load view matrix.
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    const glm::mat4 viewMatrix = camera.viewMatrix();
    glMultMatrixf(glm::value_ptr(viewMatrix));

    // Load projection matrix.
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    const glm::mat4 projectionMatrix = camera.projectionMatrix();
    glMultMatrixf(glm::value_ptr(projectionMatrix));
}

static void renderOpenGL(const Scene& scene, const Trackball& camera, int selectedLight) {
    // Normals will be normalized in the graphics pipeline.
    glEnable(GL_NORMALIZE);
    // Activate rendering modes.
    glEnable(GL_DEPTH_TEST);
    // Draw front and back facing triangles filled.
    glPolygonMode(GL_FRONT, GL_FILL);
    glPolygonMode(GL_BACK, GL_FILL);
    // Interpolate vertex colors over the triangles.
    glShadeModel(GL_SMOOTH);
    setOpenGLMatrices(camera);

    glDisable(GL_LIGHTING);
    // Render point lights as very small dots
    for (const auto& light : scene.pointLights)
        drawSphere(light.position, 0.01f, light.color);
    for (const auto& light : scene.sphericalLight)
        drawSphere(light.position, light.radius, light.color);

    if (!scene.pointLights.empty() || !scene.sphericalLight.empty()) {
        if (selectedLight < static_cast<int>(scene.pointLights.size())) {
            // Draw a big yellow sphere and then the small light sphere on top.
            const auto& light = scene.pointLights[selectedLight];
            drawSphere(light.position, 0.05f, glm::vec3(1, 1, 0));
            glDisable(GL_DEPTH_TEST);
            drawSphere(light.position, 0.01f, light.color);
            glEnable(GL_DEPTH_TEST);
        }
        else {
            // Draw a big yellow sphere and then the smaller light sphere on top.
            const auto& light = scene.sphericalLight[selectedLight - scene.pointLights.size()];
            drawSphere(light.position, light.radius + 0.01f, glm::vec3(1, 1, 0));
            glDisable(GL_DEPTH_TEST);
            drawSphere(light.position, light.radius, light.color);
            glEnable(GL_DEPTH_TEST);
        }
    }

    // Activate the light in the legacy OpenGL mode.
    glEnable(GL_LIGHTING);

    int i = 0;
    const auto enableLight = [&](const auto& light) {
        glEnable(GL_LIGHT0 + i);
        const glm::vec4 position4{ light.position, 1 };
        glLightfv(GL_LIGHT0 + i, GL_POSITION, glm::value_ptr(position4));
        const glm::vec4 color4{ glm::clamp(light.color, 0.0f, 1.0f), 1.0f };
        const glm::vec4 zero4{ 0.0f, 0.0f, 0.0f, 1.0f };
        glLightfv(GL_LIGHT0 + i, GL_AMBIENT, glm::value_ptr(zero4));
        glLightfv(GL_LIGHT0 + i, GL_DIFFUSE, glm::value_ptr(color4));
        glLightfv(GL_LIGHT0 + i, GL_SPECULAR, glm::value_ptr(zero4));
        // NOTE: quadratic attenuation doesn't work like you think it would in legacy OpenGL.
        // The distance is not in world space but in NDC space!
        glLightf(GL_LIGHT0 + i, GL_CONSTANT_ATTENUATION, 1.0f);
        glLightf(GL_LIGHT0 + i, GL_LINEAR_ATTENUATION, 0.0f);
        glLightf(GL_LIGHT0 + i, GL_QUADRATIC_ATTENUATION, 0.0f);
        i++;
    };
    for (const auto& light : scene.pointLights)
        enableLight(light);
    for (const auto& light : scene.sphericalLight)
        enableLight(light);

    // Draw the scene and the ray (if any).
    drawScene(scene);

    // Draw a colored sphere at the location at which the trackball is looking/rotating around.
    glDisable(GL_LIGHTING);
    drawSphere(camera.lookAt(), 0.01f, glm::vec3(0.2f, 0.2f, 1.0f));
}
