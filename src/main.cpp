#include "vectorMath.hpp"

#include <iostream>
#include <string>

#include <png++/png.hpp>

template <typename T, std::size_t D>
std::string toString(const Vector<T, D>& vector) {

    std::string output = "[ ";
    for (unsigned int i_dimension = 0; i_dimension < D; ++i_dimension)
        output += std::to_string(vector[i_dimension]) + " ";
    output += "]";

    return output;
}

int main(int argc, char** argv) {

    // Used to define the position, orientation, and size of the viewport.
    Vector<double, 3> viewportCenter = { 0, 0, 0 };         // The center of the viewport relative to the origin.
    // Both Vectors below must be a right angle to one another.
    Vector<double, 3> viewportHorizontalAxis = normalize(Vector<double, 3> { 1, 0, 0 }); // The direction the viewport spans horizontally, relative to viewportCenter. Points "right" in the picture.
    Vector<double, 3> viewportVerticalAxis = normalize(Vector<double, 3> { 0, 1, 0 });   // The direction the viewport spans vertically, relative to viewportCenter. Points "up" in the picture.

    unsigned int viewportWidth = 2000;
    unsigned int viewportHeight = 2000;

    unsigned int viewerDistance = 3000;
    Vector<double, 3> viewerPosition = viewportCenter + viewerDistance * normalize(viewportHorizontalAxis % viewportVerticalAxis);

    Vector<double, 3> sphereCenter = { 0, 0, -1000 };
    unsigned int sphereRadius = 1000;

    Vector<double, 3> lightPosition = { 2000, 2000, 2000 };

    png::image<png::rgba_pixel_16> image(viewportWidth, viewportHeight);

    for (unsigned int i_viewportWidth = 0; i_viewportWidth < viewportWidth; ++i_viewportWidth) {

        for (unsigned int i_viewportHeight = 0; i_viewportHeight < viewportHeight; ++i_viewportHeight) {

            Vector<double, 3> pixelPosition = viewportCenter + viewportHorizontalAxis * (i_viewportWidth - (viewportWidth - 1) / 2.0) + viewportVerticalAxis * ((viewportHeight - 1) / 2.0 - i_viewportHeight);

            double a = (pixelPosition - viewerPosition) * (pixelPosition - viewerPosition);
            double b = 2 * (viewerPosition - sphereCenter) * (pixelPosition - viewerPosition);
            double c = (viewerPosition - sphereCenter) * (viewerPosition - sphereCenter) - pow(sphereRadius, 2);

            double raystepFar = (-b + sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);
            double raystepNear = (-b - sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);

            Vector<double, 3> rayImpact;
            if (std::isnan(raystepNear) || raystepNear <= 1)
                rayImpact = viewerPosition + (pixelPosition - viewerPosition) * raystepFar;
            else
                rayImpact = viewerPosition + (pixelPosition - viewerPosition) * raystepNear;

            Vector<double, 3> N = normalize(rayImpact - sphereCenter);      // Impact surface normal.
            Vector<double, 3> L = normalize(lightPosition - rayImpact);     // Vector to light source.
            Vector<double, 3> V = normalize(viewerPosition - rayImpact);    // Vector to viewer.
            Vector<double, 3> R = normalize(2 * L * N * N - L);             // Reflection of L along N.

            double illuminationAmbient = 0.4;                           // Portion of light without direct illumination.
            double illuminationDiffuse = fmax(N * L, 0) * 0.6;          // Portion of light which is scattered off the object.
            double illuminationSpecular = pow(fmax(V * R, 0), 9) * 0.0; // Portion of light which is reflected directly to the viewer.

            if ((!std::isnan(raystepFar) && raystepFar >= 1) || (!std::isnan(raystepNear) && raystepNear >= 1)) {

                uint16_t R = (uint16_t) round(fmin(65535 * illuminationAmbient + 65535 * illuminationDiffuse + 65535 * illuminationSpecular, 65535));
                uint16_t G = (uint16_t) round(fmin(65535 * illuminationAmbient + 65535 * illuminationDiffuse + 65535 * illuminationSpecular, 65535));
                uint16_t B = (uint16_t) round(fmin(65535 * illuminationAmbient + 65535 * illuminationDiffuse + 65535 * illuminationSpecular, 65535));
                image[i_viewportHeight][i_viewportWidth] = png::rgba_pixel_16(R, G, B, pow(2, 16) - 1);
            }
            else {

                image[i_viewportHeight][i_viewportWidth] = png::rgba_pixel_16(0, 0, 0, 0);
            }
        }
    }

    image.write("output.png");

    return 0;
}