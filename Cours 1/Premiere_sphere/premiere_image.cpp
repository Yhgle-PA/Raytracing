#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define M_PI 3.14159265359

class Vector{
    public:
        explicit Vector(double x=0, double y=0, double z=0) {
            coords[0] = x;
            coords[1] = y;
            coords[2] = z;
        };
        double operator[](int i) const {return coords[i]; };
        double &operator[](int i) {return coords[i]; };
        double sqrNorm() {
            return coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2];
        }
        Vector get_normalized(){
            double n = sqrt(sqrNorm());
            return Vector(coords[0]/n, coords[1]/n, coords[2]/n);
        }


    private:
        double coords[3];
    };

Vector operator+(const Vector& a, const Vector& b){
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b){
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(double a, const Vector& b){
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, double b){
    return Vector(b*a[0], b*a[1], b*a[2]);
}
double dot(const Vector& a, const Vector& b){
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}


class Ray
{
    public:
        explicit Ray(const Vector& C, const Vector& u): C(C), u(u) {
        };
        Vector C, u;
};

class Sphere
{
    public:
        explicit Sphere(const Vector& O, double R): O(O), R(R) {
        };
        bool intersect(const Ray& r){
            double a = 1;
            double b = 2*dot(r.u, r.C-O);
            double c = (r.C-O).sqrNorm() - R*R;
            double delta = b*b - 4*a*c;
            if (delta >= 0) return true;
            return false;
        };

    private:
        Vector O;
        double R;
};

int main() {
    int W = 512;
    int H = 512;
    Vector C(0, 0, 55);
    Vector O(0, 0, 0);
    double R = 10;
    Sphere S(O, R);
    double fov = 60 * M_PI /180;
    std::vector<unsigned char> image(W*H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            Vector u(j - W/2, i - H/2, - W / (2.*tan(fov/2)));
            u = u.get_normalized();
            Ray r(C, u);
            bool inter = S.intersect(r);
            Vector color(0, 0, 0);

            if (inter) {
                color = Vector(255, 255, 255);
            };
            image[(i*W + j) * 3 + 0] = color[0];
            image[(i*W + j) * 3 + 1] = color[1];
            image[(i*W + j) * 3 + 2] = color[2];
            }
        }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);

    return 0;
}