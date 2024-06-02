#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<algorithm>
#include<vector>
#include<random>
#include<limits>

const int image_height = 600;
const int image_width = 600;
const float fov = 40;

#define pi 3.1415926535897932384626

inline float random_float() {
    //Generate a float in [0, 1).
    static std::uniform_real_distribution<float> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
}

class vec3 {
public:
    vec3() :x(0), y(0), z(0) {}
    vec3(float a, float b, float c) : x(a), y(b), z(c) {}
    vec3 operator+(const vec3 &v2) const { return vec3(x + v2.x, y + v2.y, z + v2.z); }
    vec3 &operator+=(const vec3 &v2) { x += v2.x; y += v2.y; z += v2.z; return *this; }
    vec3 operator-(const vec3 &v2) const { return vec3(x - v2.x, y - v2.y, z - v2.z); }
    vec3 operator*(const vec3 &v2) const { return vec3(x * v2.x, y * v2.y, z * v2.z); }
    vec3 operator*(float t) const { return vec3(x * t, y * t, z * t); }
    vec3 operator/(float t) const {
        float inv = 1.0 / t;
        return (*this) * inv;
    }
    vec3 operator-() const { return vec3(-x, -y, -z); }

    float LengthSquared() const { return x * x + y * y + z * z; }
    float length() const { return std::sqrt(LengthSquared()); }

    float x, y, z;
};

using rgbcolor = vec3;
using point3 = vec3;
class scene;

inline float dot(const vec3 &v1, const vec3 &v2) { return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z; }
inline vec3 cross(const vec3 &v1, const vec3 &v2) { return vec3(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x); }
inline vec3 normalize(const vec3 &v) { return v / v.length(); }
inline vec3 operator*(float t, vec3 v) { return v * t; }

class ray {
public:
    ray() = default;
    ray(const vec3 &_o, const vec3 &_d) :o(_o), d(_d) {}

    vec3 p(float t) const { return o + d * t; }

    vec3 o, d;
    mutable float tmax = std::numeric_limits<float>::infinity();
};

rgbcolor Li(const ray &r, const scene &ssene, int depth);

class camera {
public:
    camera(int _spp) :spp(_spp) { initialize(); }

    void write_file(std::ofstream &ofs, const scene &myscene) const {
        ofs << "P3" << '\n'
            << image_width << " " << image_height << '\n'
            << 255 << '\n';

        for (int j = 0; j < image_height; ++j) {
            std::clog << "\rScanlines Remaining: " << image_height - j << std::flush;
            for (int i = 0; i < image_width; ++i) {
                //ofs << i % 256 << " " << j % 256 << " " << 0 << '\n';
                write_color_in_pixel_ij(ofs, myscene, i, j);
            }
        }
        }

    void write_color_in_pixel_ij(std::ofstream &ofs, const scene &myscene, int i, int j) const {
        float r, g, b;
        rgbcolor color(0, 0, 0);

        for (int kk = 0; kk < spp; ++kk) {
            color += Li(random_ray_in_pixel_ij(i, j), myscene, 200);
            //color = color + Li(random_ray_in_pixel_ij(i, j), myscene, 5);
        }
        color = color / static_cast<float>(spp);

        //clamping
        r = std::clamp(color.x, 0.0f, 0.9999f);
        g = std::clamp(color.y, 0.0f, 0.9999f);
        b = std::clamp(color.z, 0.0f, 0.9999f);

        //approximate gamma correction
        r = std::sqrt(r);
        g = std::sqrt(g);
        b = std::sqrt(b);

        //map [0, 1) to [0, 255]
        int ri = static_cast<int>(256 * r);
        int gi = static_cast<int>(256 * g);
        int bi = static_cast<int>(256 * b);

        ofs << ri << " " << gi << " " << bi << "\n";
    }

    ray random_ray_in_pixel_ij(int i, int j) const {
        point3 pixelij_corner = pixel00_corner + i * pixel_delta_u + j * pixel_delta_v;
        float a = random_float();
        float b = random_float();
        point3 random_in_pixelij = pixelij_corner + pixel_delta_u * a + pixel_delta_v * b;
        return ray(p, random_in_pixelij - p);
    }

private:
    void initialize() {
        screen_height = 2 * std::tan(pi * fov / 360) * vec3(p2s - p).length();
        neg_screen_width = -screen_height * (static_cast<double>(image_width) / image_height);

        pixel_delta_u = u * neg_screen_width / image_width;
        pixel_delta_v = -v * screen_height / image_height;
        pixel00_loc = p2s - 0.5 * u * neg_screen_width + 0.5 * v * screen_height + 0.5 * pixel_delta_u + 0.5 * pixel_delta_v;
        pixel00_corner = p2s - 0.5 * u * neg_screen_width + 0.5 * v * screen_height;
    }

    point3 p{ 278, 278, -800 }, p2s{ 278, 278, 0 };
    vec3 u{ 1, 0, 0 }, v{ 0, 1, 0 }, w{ 0, 0, 1 };

    vec3 pixel_delta_u, pixel_delta_v;
    point3 pixel00_loc, pixel00_corner;
    float neg_screen_width, screen_height;

    public:
    int spp;

};

struct interaction {
    vec3 p;
    vec3 normal;
    rgbcolor matcolor;
    bool emissive;
};

class quad {
public:
    quad(const vec3 &_p, const vec3 &_u, const vec3 &_v, const rgbcolor &_c, bool emit = false) :p(_p), u(_u), v(_v), c(_c), emissive(emit) {
        n = cross(u, v);
        Q = dot(p, n);
        w = n / n.LengthSquared();
    }
    quad() : quad(vec3(0, 0, 0), vec3(0, 0, 0), vec3(0, 0, 0), rgbcolor(0, 0, 0), false) {}

    bool intersect(const ray &r, interaction &intr) const {
        float denom = dot(r.d, n);
        
        if (std::abs(denom) < 2e-10) {
            return false;
        }
        
        float t = (Q - dot(r.o, n)) / denom;

        //std::cout << t << "     r.d: " << r.d.x << " " << r.d.y << " " << r.d.z << "   n: " << n.x << " " << n.y << ' ' << n.z << std::endl;

        if (t >= r.tmax)return false;
        if (t < 1e-10)return false;

        vec3 point = r.p(t);
        vec3 vec2decomp = point - p;
        float u_weight = dot(w, cross(vec2decomp, v));
        float v_weight = dot(w, cross(u, vec2decomp));

        if (u_weight < 0 || u_weight>1 || v_weight < 0 || v_weight>1) {
            return false;
        }

        r.tmax = t;
        intr.p = point;
        intr.normal = normalize(n);
        intr.matcolor = c;
        intr.emissive = emissive;

        return true;
    }

private:
    vec3 p;
    vec3 u, v;
    rgbcolor c;
    bool emissive;

    vec3 n;
    float Q;
    vec3 w;
};

class scene {
public:
    scene() = default;
    scene(std::shared_ptr<quad> q) {
        quad_list.push_back(q);
    }
    scene(const std::vector<std::shared_ptr<quad>> &list) {
        for (const auto &ptr : list) {
            quad_list.push_back(std::make_shared<quad>(*ptr));
        }
    }

    bool intersect(const ray &r, interaction &intr) const {
        int n = quad_list.size();
        if (n == 0)return false;
        bool flag = false;
        for (const auto &ptr : quad_list) {
            if (ptr->intersect(r, intr)) {
                flag = true;
            }
        }
        return flag;
    }

private:
    std::vector<std::shared_ptr<quad>> quad_list;
};

vec3 cosine_sample(const vec3 &n) {
    float r1 = random_float(), r2 = random_float();
    float phi = 2 * pi * r1;
    float a = std::cos(phi) * std::sqrt(r2);
    float b = std::sin(phi) * std::sqrt(r2);
    float c = std::sqrt(1 - r2);

    vec3 u = normalize(cross(std::abs(n.y) < 0.5 ? vec3(0, 1, 0) : vec3(0, 0, 1), n));
    vec3 v = cross(n, u);
    return normalize(vec3(a * u + b * v + c * n));
}

rgbcolor Li(const ray &r, const scene &ssene, int depth) {
    if (depth <= 0)return rgbcolor(0, 0, 0);
    interaction intr;
    if (!ssene.intersect(r, intr)) {
        return rgbcolor(0.0, 0.0, 0.0);
    }

    float cos_flag = dot(intr.normal, r.d);
    if (intr.emissive) {
        if (cos_flag < 0)
            return intr.matcolor;
        else return rgbcolor(0, 0, 0);
    }

    if (random_float() < 0.05) return rgbcolor(0, 0, 0);
    
    vec3 wo = cosine_sample(intr.normal);
    point3 po = intr.p + intr.normal * (cos_flag < 0 ? 0.01 : -0.01);
    return intr.matcolor * Li(ray(po, wo), ssene, depth - 1) * 0.95;
}

int main(int argc, char *argv[]) {
    std::ofstream ofs("CornellBox.ppm");

    rgbcolor white(.73, .73, .73);
    rgbcolor red(.65, .05, .05);
    rgbcolor green(.12, .45, .15);

    std::vector<std::shared_ptr<quad>> quads;
    quads.push_back(std::make_shared<quad>(point3(213, 554, 227), vec3(130, 0, 0), vec3(0, 0, 105), rgbcolor(15, 15, 15), true));

    quads.push_back(std::make_shared<quad>(point3(555, 0, 0), vec3(0, 0, 555), vec3(0, 555, 0), green));
    quads.push_back(std::make_shared<quad>(point3(0, 0, 555), vec3(0, 0, -555), vec3(0, 555, 0), red));
    quads.push_back(std::make_shared<quad>(point3(0, 555, 0), vec3(555, 0, 0), vec3(0, 0, 555), white));
    quads.push_back(std::make_shared<quad>(point3(0, 0, 555), vec3(555, 0, 0), vec3(0, 0, -555), white));
    quads.push_back(std::make_shared<quad>(point3(555, 0, 555), vec3(-555, 0, 0), vec3(0, 555, 0), white));

    quads.push_back(std::make_shared<quad>(point3(130, 165, 65), vec3(-48, 0, 160), vec3(160, 0, 49), white));
    quads.push_back(std::make_shared<quad>(point3(290, 0, 114), vec3(0, 165, 0), vec3(-50, 0, 158), white));
    quads.push_back(std::make_shared<quad>(point3(130, 0, 65), vec3(0, 165, 0), vec3(160, 0, 49), white));
    quads.push_back(std::make_shared<quad>(point3(82, 0, 225), vec3(0, 165, 0), vec3(48, 0, -160), white));
    quads.push_back(std::make_shared<quad>(point3(240, 0, 272), vec3(0, 165, 0), vec3(-158, 0, -47), white));

    quads.push_back(std::make_shared<quad>(point3(265, 330, 296), vec3(49, 0, 160), vec3(158, 0, -49), white));
    quads.push_back(std::make_shared<quad>(point3(423, 0, 247), vec3(0, 330, 0), vec3(49, 0, 159), white));
    quads.push_back(std::make_shared<quad>(point3(472, 0, 406), vec3(0, 330, 0), vec3(-158, 0, 50), white));
    quads.push_back(std::make_shared<quad>(point3(314, 0, 456), vec3(0, 330, 0), vec3(-49, 0, -160), white));
    quads.push_back(std::make_shared<quad>(point3(265, 0, 296), vec3(0, 330, 0), vec3(158, 0, -49), white));

    scene myscene(quads);

    int spp = 1024;
    if (argc >= 2)spp = atoi(argv[1]);
    std::cout << "spp = " << spp << std::endl;
    camera cm1(spp);
    cm1.write_file(ofs, myscene);
    std::cout << std::endl << "Done. " << std::endl;
    return 0;
}