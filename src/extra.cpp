#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include <framework/trackball.h>

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of Depth of Field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithDepthOfField(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableDepthOfField) {
        return;
    }

    // ...
}

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of motion blur. Here, you integrate over a time domain, and not just the pixel's image domain,
// to give objects the appearance of "fast movement".
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithMotionBlur(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableMotionBlur) {
        return;
    }

}


// Given two number n and k, this computes the binomial coefficient.
int binomialCoefficients(int n, int k)
{
    if (k == 0 || k == n)
        return 1;
    return binomialCoefficients(n - 1, k - 1) + binomialCoefficients(n - 1, k);
}


// This method generates the weights for the filter using the binomial coefficient.
void calculateWeightsBinomial(std::vector<float>& weights, int n)
{
    float total = 0.0f;

    for (int k = 0; k <= n; k++) {
        weights.push_back(binomialCoefficients(n, k));
        total += weights[k];
    }

    for (int k = 0; k <= n; k++) {
        weights[k] = (float)(weights[k]) / (float)(total);
    }
}

// This method generates the weights for the filter using the binomial coefficient.
void calculateWeightsGrayscale(const Screen& screen, std::vector<float>& weights)
{
    weights.clear(); // Ensure the vector is empty before filling it

    float totalValue = 0.0f;

    // First pass: Compute the total intensity of all pixels
    for (const auto& pixel : screen.pixels()) {
        float intensity = pixel.r + pixel.g + pixel.b; // Sum RGB components
        totalValue += intensity;
    }

    // Prevent division by zero
    if (totalValue == 0.0f) {
        totalValue = 1.0f;
    }

    // Second pass: Compute normalized weights
    for (const auto& pixel : screen.pixels()) {
        float intensity = pixel.r + pixel.g + pixel.b;
        weights.push_back(intensity / totalValue);
    }
}


void translateImageToScreen(const Image& image, Screen& screen, int n)
{
    n = std::clamp(n, 2, 20);

    int regionSize = 2 * n;
    int usedWidth = std::min(regionSize, image.width);
    int usedHeight = std::min(regionSize, image.height);

    screen = Screen(glm::ivec2(usedWidth, usedHeight));

    for (int y = 0; y < usedHeight; y++) {
        for (int x = 0; x < usedWidth; x++) {
            // Convert (x, y) coordinates to a 1D index in the image
            int imageIndex = y * image.width + x;

            glm::vec3 pixelColor = image.get_pixel<3>(imageIndex);
            screen.setPixel(x, y, pixelColor);
        }
    }
}


glm::vec3 binaryMapping(glm::vec3 pixelColor, float t) {
    if (pixelColor.x > t || pixelColor.y > t || pixelColor.z > t) {
        return glm::vec3(1.0f);
    } else{
        return glm::vec3(0.0f);
    }
}

glm::vec3 linearMapping(glm::vec3 pixelColor, float t)
{
    if (pixelColor.x > t || pixelColor.y > t || pixelColor.z > t) {
        return (pixelColor - t) / (1.0f - t);
    } else {
        return glm::vec3(0.0f);
    }
}

glm::vec3 truncateMapping(glm::vec3 pixelColor, float t)
{
    if (pixelColor.x > t || pixelColor.y > t || pixelColor.z > t) {
        return pixelColor;
    } else {
        return glm::vec3(0.0f);
    }
}



void applyFilter(const Scene& scene, const Features& features, const Trackball& camera, Screen& image, std::vector<float> weights)
{
    Screen lights(image.resolution());
    float t = features.extra.treshold;

    if (features.extra.linear || features.extra.truncate) {
        if (features.extra.linear) {
            
            //LINEAR MAPPING
            for (int i = 0; i < image.resolution().x; i++) {
                for (int j = 0; j < image.resolution().y; j++) {
                    glm::vec3 pixelColor = image.pixels()[image.indexAt(i, j)];
                    lights.setPixel(i, j, linearMapping(pixelColor, t));
                }
            }
        } else {
            //TRUNCATE MAPPING
            for (int i = 0; i < image.resolution().x; i++) {
                for (int j = 0; j < image.resolution().y; j++) {
                    glm::vec3 pixelColor = image.pixels()[image.indexAt(i, j)];
                    lights.setPixel(i, j, truncateMapping(pixelColor, t));
                }
            }
        }
    } else {
        //BINARY MAPPING
        for (int i = 0; i < image.resolution().x; i++) {
            for (int j = 0; j < image.resolution().y; j++) {
                glm::vec3 pixelColor = image.pixels()[image.indexAt(i, j)];
                lights.setPixel(i, j, binaryMapping(pixelColor, t));
            }
        }
    }

    int size = weights.size();

    Screen horizontal(image.resolution());
    for (int j = 0; j < image.resolution().y; j++) {
        for (int i = 0; i < image.resolution().x; i++) {
            glm::vec3 blurredPixel { 0.0f };
            for (int k = -size / 2; k <= size / 2; k++) {
                int idx = i + k;
                if (idx >= 0 && idx < image.resolution().x) {
                    blurredPixel += lights.pixels()[lights.indexAt(idx, j)] * weights[k + size / 2];
                }
            }
            horizontal.setPixel(i, j, blurredPixel);
        }
    }

    Screen vertical(image.resolution());
    for (int i = 0; i < image.resolution().x; i++) {
        for (int j = 0; j < image.resolution().y; j++) {
            glm::vec3 blurredPixel { 0.0f };
            for (int k = -size / 2; k <= size / 2; k++) {
                int idx = j + k;
                if (idx >= 0 && idx < image.resolution().y) {
                    blurredPixel += horizontal.pixels()[horizontal.indexAt(i, idx)] * weights[k + size / 2];
                }
            }
            vertical.setPixel(i, j, blurredPixel);
        }
    }

    // Combine the bloom effect with the original image.
    for (int i = 0; i < image.resolution().x; i++) {
        for (int j = 0; j < image.resolution().y; j++) {
            glm::vec3 mask = lights.pixels()[lights.indexAt(i, j)];
            glm::vec3 originalColor = image.pixels()[image.indexAt(i, j)];
            glm::vec3 bloomColor = vertical.pixels()[vertical.indexAt(i, j)];
            glm::vec3 finalColor = originalColor + bloomColor;
            // show the bloom effect on the image
            image.setPixel(i, j, finalColor);
            
            //show only the MASK
            //image.setPixel(i, j, mask);
            
            // show the mask with the applied Bloomfiltrer
            //image.setPixel(i, j, bloomColor);
        }
    }
}

// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
    }

    int n = features.extra.n;
    if (n % 2 == 1) {
        n++;
    }
    std::vector<float> weights;


    //If a valid path is given
    if (!features.extra.filterImagePath.empty()) {
        Image a(features.extra.filterImagePath);
        Screen upload(glm::ivec2(1, 1));
        translateImageToScreen(a,upload, n);
        calculateWeightsGrayscale(upload, weights);

    } else {
        calculateWeightsBinomial(weights, n);
    }
            
    applyFilter(scene, features, camera, image, weights);
    

}


// TODO; Extra feature
// Given a camera ray (or reflected camera ray) and an intersection, evaluates the contribution of a set of
// glossy reflective rays, recursively evaluating renderRay(..., depth + 1) along each ray, and adding the
// results times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderRayGlossyComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // Generate an initial specular ray, and base secondary glossies on this ray
    // auto numSamples = state.features.extra.numGlossySamples;
    // ...
}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You will have to add support for environment textures
// to the Scene object, and provide a scene with the right data to supply this.
// - state; the active scene, feature config, bvh, and sampler
// - ray;   ray object
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
glm::vec3 sampleEnvironmentMap(RenderState& state, Ray ray)
{
    if (state.features.extra.enableEnvironmentMap) {
        // Part of your implementation should go here
        return glm::vec3(0.f);
    } else {
        return glm::vec3(0.f);
    }
}
