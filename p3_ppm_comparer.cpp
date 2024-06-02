#include<iostream>
#include<fstream>

int main(int argc, char *argv[]) {

    if (argc < 3) {
        std::cout << "argc < 3" << std::endl;
        return 0;
    }

    std::ifstream f1(argv[1]), f2(argv[2]);
    std::ofstream f3;
    if (argc >= 4) {
        f3.open(argv[3]);
    }
    else {
        f3.open("difference.ppm");
    }

    std::string s1, s2;
    f1 >> s1, f2 >> s2;
    if (s1 != "P3" || s2 != "P3") {
        f3 << "They are not both P3" << std::endl;
        return 0;
    }

    int w1, h1, i1, w2, h2, i2;
    f1 >> w1 >> h1 >> i1, f2 >> w2 >> h2 >> i2;
    if (w1 != w2 || h1 != h2 || i1 != i2) {
        f3 << "Incompatible in image resolution or color Maxval!" << std::endl;
        return 0;
    }

    f3 << "P3\n" << w1 << ' ' << h1 << "\n" << i1 << "\n";

    long long pixel_count = w1 * h1;

    for (long long i = 0; i < pixel_count; ++i) {
        int r1, g1, b1, r2, g2, b2;
        f1 >> r1 >> g1 >> b1, f2 >> r2 >> g2 >> b2;
        //int r3 = r1 - r2, g3 = g1 - g2, b3 = b1 - b2;
        int r3 = std::abs(r1 - r2), g3 = std::abs(g1 - g2), b3 = std::abs(b1 - b2);
        f3 << r3 << " " << g3 << " " << b3 << "\n";
    }

    return 0;
}