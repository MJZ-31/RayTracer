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
    Vector<double, 3> viewportHorizontalAxis = { 1, 0, 0 }; // The direction the viewport spans horizontally, relative to viewportCenter. Points "right" in the picture.
    Vector<double, 3> viewportVerticalAxis = { 0, 1, 0 };   // The direction the viewport spans vertically, relative to viewportCenter. Points "up" in the picture.

    unsigned int viewportWidth = 10;
    unsigned int viewportHeight = 10;

    unsigned int viewerDistance = 100;
    Vector<double, 3> viewerPosition = viewportCenter + viewerDistance * normalize(viewportHorizontalAxis % viewportVerticalAxis);
    std::cout << "Viewer: " << toString(viewerPosition) << std::endl;

    Vector<double, 3> sphereCenter = { 0, 0, -20 };
    unsigned int sphereRadius = 5;

    for (unsigned int i_viewportWidth = 0; i_viewportWidth < viewportWidth; ++i_viewportWidth) {

        for (unsigned int i_viewportHeight = 0; i_viewportHeight < viewportHeight; ++i_viewportHeight) {

            Vector<double, 3> pixel = viewportCenter + Vector<double, 3> { i_viewportWidth - (viewportWidth - 1) / 2.0, (viewportHeight - 1) / 2.0 - i_viewportHeight, 0.0 };
            // std::cout << "Pixel " << std::to_string(i_viewportWidth) << ", " << std::to_string(i_viewportHeight) << ": " << toString(pixel) << std::endl;
            // std::cout << "Ray: " << toString(pixel - viewerPosition) << std::endl;

            double a = 0;
            for (unsigned int i_dimension = 0; i_dimension < 3; ++i_dimension)
                a += pow(pixel[i_dimension], 2) - 2 * viewerPosition[i_dimension] * pixel[i_dimension] + pow(viewerPosition[i_dimension], 2);

            double b = 0;
            for (unsigned int i_dimension = 0; i_dimension < 3; ++i_dimension)
                b += 2 * (viewerPosition[i_dimension] * pixel[i_dimension] - pow(viewerPosition[i_dimension], 2) - pixel[i_dimension] * sphereCenter[i_dimension] + viewerPosition[i_dimension] * sphereCenter[i_dimension]);

            double c = 0;
            for (unsigned int i_dimension = 0; i_dimension < 3; ++i_dimension)
                c += pow(viewerPosition[i_dimension], 2) - 2 * viewerPosition[i_dimension] * sphereCenter[i_dimension] + pow(sphereCenter[i_dimension], 2);
            c -= pow(sphereRadius, 2);

            // std::cout << std::to_string(a) << " " << std::to_string(b) << " " << std::to_string(c) << std::endl;

            double raystep = (-b - sqrt(pow(b, 2) - 4 * a * c)) / 2 * a;
            // if (!std::isnan(raystep))
                std::cout << raystep << " ";
        }
        std::cout << std::endl;
    }

    png::image<png::rgba_pixel> image(viewportWidth, viewportHeight);
    image.write("output.png");

    return 0;
}