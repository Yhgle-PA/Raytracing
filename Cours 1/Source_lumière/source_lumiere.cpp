#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define M_PI 3.14159265359

#include <algorithm>


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
Vector operator/(const Vector& a, double b){
    return Vector(a[0]/b, a[1]/b, a[2]/b);
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
        bool intersect(const Ray& r, Vector& P, Vector& N){
            double a = 1;
            double b = 2*dot(r.u, r.C-O);
            double c = (r.C-O).sqrNorm() - R*R;
            double delta = b*b - 4*a*c;
            if (delta < 0) return false;

            double sqDelta = sqrt(delta);
            double t2 = (-b + sqDelta)/ (2*a);
            if (t2 < 0) return false;

            double t;
            double t1 = (-b - sqDelta) / (2*a);
            if (t1 > 0)
                t = t1;
            else
                t = t2;

            P = r.C + t*r.u;
            N = (P - O).get_normalized();

            return true;
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
    double I = 1E7;
    Vector rho(1, 0, 0);
    Vector L(-10, 20, 40);

    std::vector<unsigned char> image(W*H * 3, 0);

    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector u(j - W/2, i - H/2, - W / (2.*tan(fov/2)));
            u = u.get_normalized();
            Ray r(C, u);
            Vector P, N;
            bool inter = S.intersect(r, P, N);
            Vector color(0, 0, 0);

            if (inter) {
                Vector PL = L - P;
                double d = sqrt(PL.sqrNorm());
                color = I/(4*M_PI*d*d)*rho/M_PI*std::max(0., dot(N, PL/d));
            };
            image[((H - i - 1)* W + j) * 3 + 0] = std::min(255., color[0]);
            image[((H - i - 1)* W + j) * 3 + 1] = std::min(255., color[1]);
            image[((H - i - 1)* W + j) * 3 + 2] = std::min(255., color[2]);
            }
        }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);

    return 0;
}