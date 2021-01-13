#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define M_PI 3.14159265359

#include <algorithm>
#include <complex>


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
Vector operator-(const Vector& a){
    return Vector(-a[0], -a[1], -a[2]);
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
double sqr(const float a){
    return a*a;
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
        explicit Sphere(const Vector& O, double R, const Vector& albedo, bool isMirror = false, bool isTransp = false, double n = 1.4, bool isCreux = false): O(O), R(R), albedo(albedo), isMirror(isMirror), isTransp(isTransp), n(n), isCreux(isCreux){
        };
        bool intersect(const Ray& r, Vector& P, Vector& N, double& t){
            double a = 1;
            double b = 2*dot(r.u, r.C-O);
            double c = (r.C-O).sqrNorm() - R*R;
            double delta = b*b - 4*a*c;
            if (delta < 0) return false;

            double sqDelta = sqrt(delta);
            double t2 = (-b + sqDelta)/ (2*a);
            if (t2 < 0) return false;

            double t1 = (-b - sqDelta) / (2*a);
            if (t1 > 0)
                t = t1;
            else
                t = t2;

            P = r.C + t*r.u;
            N = (P - O).get_normalized();
            return true;
        };
        bool isMirror, isTransp, isCreux;
        Vector albedo;
        double n;
    private:
        Vector O;
        double R;
};

class Scene{
public:
    Scene() {};
    bool intersect(const Ray& r, Vector& P, Vector& N, Vector &albedo, bool &mirror, bool &transp, double &n, double &t){
        t = 1E10;
        bool has_inter = false;
        for(int i = 0; i < objects.size(); i++){
            Vector localP, localN;
            double localt;
            if (objects[i].intersect(r, localP, localN, localt) && localt < t){
                has_inter = true;
                t = localt;
                P = localP;
                N = localN;
                if (objects[i].isCreux) {
                    N = -N;
                }
                albedo = objects[i].albedo;
                mirror = objects[i].isMirror;
                transp = objects[i].isTransp;
                n = objects[i].n;
            }
        }
        return has_inter;
     }

     Vector getColor(const Ray& r, int rebond) {
            Vector P, N, albedo;
            double t, n2;
            bool mirror, transp;
            bool inter=intersect(r, P, N, albedo, mirror, transp, n2, t);
            Vector color(0, 0, 0);
            if (rebond > 10) return color;

            if (inter) {
                if (mirror) {
                    Vector reflectedDir = r.u - 2*dot(r.u, N)*N;
                    Ray reflectedRay(P + 1E-5*N, reflectedDir);
                    return getColor(reflectedRay, rebond + 1);
                } else {
                    if (transp) {

                        double n1 = 1;
                        Vector N2 = N;
                        if (dot(r.u, N) > 0) {
                            std::swap(n1, n2);
                            N2 = -N;
                        }
                        Vector Tt = n1 / n2 * (r.u - dot(r.u, N2)*N2);
                        double rad = 1 - sqr(n1 / n2) * (1 - sqr(dot(r.u, N2)));
                        if (rad < 0) {
                            Vector reflectedDir = r.u - 2*dot(r.u, N)*N;
                            Ray reflectedRay(P + 1E-5*N, reflectedDir);
                            return getColor(reflectedRay, rebond + 1);

                        }
                        Vector Tn = -sqrt(rad)*N2;
                        Vector refractedDir = Tt + Tn;
                        Ray refractedRay(P - 1E-5*N2, refractedDir);
                        return getColor(refractedRay, rebond + 1);
                    }
                    else {
                        Vector PL = L - P;
                        double d = sqrt(PL.sqrNorm());
                        Vector shadowP, shadowN, shadowAlbedo;
                        double shadowt, shadown2;
                        Ray shadowRay(P+1E-5*N, PL/d);
                        bool shadowMirror, shadowTransp;
                        bool shadowInter = intersect(shadowRay, shadowP, shadowN, shadowAlbedo, shadowMirror, shadowTransp, shadown2, shadowt);
                        if (shadowInter && shadowt < d) {
                        } else {
                            color = I/(4*M_PI*d*d)*albedo/M_PI*std::max(0., dot(N, PL/d));
                        }
                        return color;
                    }
                }
            };
            return color;
        }

    std::vector<Sphere> objects;
    Vector L;
    double I;
};

int main() {
    int W = 1000;
    int H = 1000;

    Scene scene;
    Vector C(0, 0, 55);
    Sphere S1(Vector(-12, 0, 0), 10, Vector(1., 1., 1.), true);
    Sphere S2(Vector(8, 4, 0), 10, Vector(1., 0., 0.), false, true, 1.4);
    Sphere S3(Vector(8, 4, 0), 9.8, Vector(1., 0., 0.), false, true, 1.4, true);
    Sphere S4(Vector(0, -5, 15), 5, Vector(1., 1., 0.), false, true);
    Sphere Ssol(Vector(0, -1000, 0), 990, Vector(1., 1., 1.));
    Sphere Smur1(Vector(-1000, 0, 0), 940, Vector(1., 0., 0.));
    Sphere Smur2(Vector(1000, 0, 0), 940, Vector(1., 0., 0.));
    Sphere Smur3(Vector(0, 0, 1000), 940, Vector(0., 0., 0.));
    Sphere Smur4(Vector(0, 0, -1000), 940, Vector(0., 1., 0.));
    Sphere Splafond(Vector(0, 1000, 0), 940, Vector(0., 0., 1.));
    scene.objects.push_back(S1);
    scene.objects.push_back(S2);
    scene.objects.push_back(S3);
    scene.objects.push_back(S4);
    scene.objects.push_back(Ssol);
    scene.objects.push_back(Splafond);
    scene.objects.push_back(Smur1);
    scene.objects.push_back(Smur2);
    scene.objects.push_back(Smur3);
    scene.objects.push_back(Smur4);
    double fov = 60 * M_PI / 180;
    scene.I = 5E9;
    scene.L = Vector(-10, 20, 40);

    std::vector<unsigned char> image(W*H * 3, 0);

    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector u(j - W/2, i - H/2, - W / (2.*tan(fov/2)));
            u = u.get_normalized();
            Ray r(C, u);
            Vector color = scene.getColor(r, 0);
            image[((H - i - 1)* W + j) * 3 + 0] = std::min(255., std::pow(color[0], 0.45));
            image[((H - i - 1)* W + j) * 3 + 1] = std::min(255., std::pow(color[1], 0.45));
            image[((H - i - 1)* W + j) * 3 + 2] = std::min(255., std::pow(color[2], 0.45));
            }
        }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);

    return 0;
}
