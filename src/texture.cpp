#include "texture.h"
#include "render.h"
#include <framework/image.h>
#include <fmt/core.h>
#include <glm/gtx/string_cast.hpp>

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// the nearest texel to the coordinates is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the nearest corresponding texel
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureNearest(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)

    //Get image dimensions
    int height = image.height;
    int width = image.width;

    //Adjust texCoord for pixel center (0.5, 0.5)
    float u = texCoord.x * width - 0.5f;
    float v = texCoord.y * height - 0.5f;

    //Round to nearest texel
    int x = std::round(u);
    int y = std::round(v);

    //Clamp the indices within bounds
    x = std::clamp(x, 0, width - 1);
    y = std::clamp(y, 0, height - 1);

    //Calculate the image index
    int index = y * width + x;

    //Return the calculated index
    return image.get_pixel(index);
}

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// a bilinearly interpolated texel is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the filter of the corresponding texels
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureBilinear(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)

    //Get image dimensions
    int height = image.height;
    int width = image.width;

    //Adjust texCoord for pixel center (0.5, 0.5)
    float u = texCoord.x * width - 0.5f;
    float v = texCoord.y * height - 0.5f;

    //Calculate texel indices
    int x1 = std::clamp(static_cast<int>(std::floor(u)), 0, width - 1);
    int x2 = std::clamp(x1 + 1, 0, width - 1);
    int y1 = std::clamp(static_cast<int>(std::floor(v)), 0, height - 1);
    int y2 = std::clamp(y1 + 1, 0, height - 1);

    //Compute weight of the pixels
    float xWeight = u - x1;
    float yWeight = v - y1;

    //Get the image pixels
    glm::vec3 pixel1 = image.get_pixel(y1 * width + x1);
    glm::vec3 pixel2 = image.get_pixel(y1 * width + x2);
    glm::vec3 pixel3 = image.get_pixel(y2 * width + x1);
    glm::vec3 pixel4 = image.get_pixel(y2 * width + x2);

    //Bilinear interpolation
    glm::vec3 mix1 = glm::mix(pixel1, pixel2, xWeight);
    glm::vec3 mix2 = glm::mix(pixel3, pixel4, xWeight);

    //Return the interpolation
    return glm::mix(mix1, mix2, yWeight);
}