#include <random>
#include <functional>
#include <cmath>
#include <numbers>
#include <complex>
#include <iostream>
#include <algorithm>
#include <initializer_list>
#include <fstream>
#include <sstream>
#include <stdlib.h>


template <typename T, int sz>
struct ComplexArray {
    T real[sz], imag[sz];
};



template <int n> using ComplexArrayFloat = ComplexArray<float, n>;
using ComplexT = ComplexArrayFloat<1>;


constexpr float pi = std::numbers::pi_v<float>;
constexpr int arrsz(int k, int d) { return 1 << (k * d); }
inline int pow2(int n) { return 1 << n; };



template <typename T>
inline T sqre(T x) { return x * x; };


template<int N>
struct OmegaTabel {
    inline static float omega[N + N / 4];
    inline static OmegaTabel<N / 2> next;

    static void make_omega() {
        for (int i = 0;i < N; i++)
            omega[i] = cos(2 * pi * i / N);

        for (int i = 0;i < N / 4; i++)
            omega[N + i] = omega[i];

        if constexpr (N > 1)
            next.make_omega();
    }
    
    static float* cos() {
        return omega;
    }

    static float* sin() {
				    return omega + N/4;
    }
};


template<>
struct OmegaTabel<1> {
    inline static float omega[1];

    static void make_omega() {
        omega[0] = 1.0f;
    }
};


template<int n>
struct FFTPack : ComplexArrayFloat<n>
{
    static constexpr int nby2 = n / 2;
    static FFTPack<n / 2> even, odd;

    template<int N>
    inline void butterfly(int s, FFTPack<N>& input) {
        even.butterfly(s * 2, input);
        odd.butterfly(s * 2 + 1, input);
        float* cos = OmegaTabel<n>::cos();
        float* sin = OmegaTabel<n>::sin();

        for (int i = 0; i < nby2;i++) {

            float even_real = sin[i] * even.real[i] - cos[i] * even.imag[i];
            float even_imag = sin[i] * even.imag[i] + cos[i] * even.real[i];

            this->real[i] = odd.real[i] + even_real;
            this->imag[i] = odd.imag[i] + even_imag;

            this->real[i + nby2] = odd.real[i] - even_real;
            this->imag[i + nby2] = odd.imag[i] - even_imag;
        }

    }
};


template<int n>
FFTPack<n / 2> FFTPack<n>::even;

template<int n>
FFTPack<n / 2> FFTPack<n>::odd;


template<>
struct FFTPack<1> : ComplexArrayFloat<1>
{
    template<int N>
    inline void butterfly(int s, FFTPack<N>& input) {
        real[0] = input.real[s];
        imag[0] = input.imag[s];
    }
};


struct Indx { int i, j; };


template <int sz, int d> 
struct inverse2DFFT {
    static constexpr int sz_sq = sz * sz;
    FFTPack<sz> input[sz];
    inverse2DFFT<sz, d - 1> next_axis;
    inverse2DFFT() {}

    FFTPack<sz>* output() {
        if constexpr (d == 0) return &input[0];
        else return next_axis.output();
    }


    void fft() {
        if constexpr (d > 0) {
            for (int i = 0; i <  sz;i++)
                next_axis.input[i].template butterfly<sz>(0, input[i]);
            next_axis.fft();
        }
    }
};


template<int sz>
struct inverse2DFFT<sz, -1> {};// d dimentional fft requires d+1 arrays of FFTpack<sz

constexpr float mean = 0.f;
constexpr float standard = 1;



template <int k>
struct ComplexNoise {
    static constexpr int sz = 1 << k;
    static constexpr int sz_sq = sz*sz;

    ComplexArrayFloat<sz>* noise;
    float spectral_bias[sz_sq];
    inverse2DFFT<sz, 2> fft;
    std::random_device seed_gen;
    std::normal_distribution<float> normal{ mean, standard };
    int seed;

    void set_seed(int val) { seed = val; }
    void gen_seed() {
        seed = seed_gen();
    }

    void inverse_fft() {
        fft.fft();
    }

    inline int index(int i, int j) { return i + sz * j; }


    template<typename T>
    void set_vals(T* arr, T val, std::initializer_list <Indx> lst) {
        for (auto& x : lst) arr[index(x.i, x.j)] = val;
    }

    void init_arrays() {

        auto normal_gen = std::bind(normal, std::default_random_engine(seed));
        int z = sz / 2;
        for (int i = 0; i < z; i++) {
            for (int j = 0; j < z; j++) {

                int is = sz - i - 1;
                int js = sz - j - 1;

                float val = std::exp(-(sqre(z - i - 1) + sqre(z - j - 1)) * 0.005);
                set_vals(spectral_bias, val, { {i,j} , {is,j},{i,js},{is,js} });

                noise[i].real[j] = noise[is].real[js] = normal_gen();
                noise[is].real[j] = noise[i].real[js] = normal_gen();

                float conj1 = normal_gen();

                noise[i].imag[j] = conj1;
                noise[is].imag[js] = -conj1;

                float conj2 = normal_gen();

                noise[is].imag[j] = conj2;
                noise[i].imag[js] = -conj2;

            }
        }


        noise[0].imag[0] = noise[z].imag[z] = noise[z].imag[sz - 1] = noise[sz - 1].imag[z] = 0;
    }

    void apply_spectral_bias() {
        for (int i = 0; i < sz;i++) {
            for (int j = 0; j < sz; j++) {
                noise[i].real[j] *= 1.1*spectral_bias[index(i,j)];
                noise[i].imag[j] *= 1.1*spectral_bias[index(i,j)];
            }
        }
    }
    
    ComplexNoise() {
        OmegaTabel<sz>::make_omega();
        noise = fft.input;
        //memset(noise, 0.f,sizeof(float)*sz*sz); // temporay for degugging;
        gen_seed();
        init_arrays();
    }
};


template<typename T>
void write_var_to(T& var,std::stringstream& stream, const char* lst_name) {
    stream<< lst_name << " = ";
    stream << var<<"\n";
}

void write_lst_to(float* arr, int sz, int d, int stride, std::stringstream& lst_string, const char* lst_name) {

    lst_string << lst_name;
    lst_string << " = [\n";
    for (int i = 0; i < sz; i++)
    {
        lst_string << "[";
        for (int j = 0; j < sz; j++)
        {
            float lst_val = arr[j * stride + d * i];
            lst_string << lst_val;
            if (j != sz - 1) lst_string << ",";
        }
        lst_string << "]\n";
        if (i != sz - 1) lst_string << ",";
    }
    lst_string << "]\n";

}

template<typename T>
T clip (T val ,T min, T max ) {
    if (min < val < max) return val;
    else if (val < min) return min;
    return max;

}
int main() {

    static ComplexNoise<7> cn;
    int  z = cn.sz / 2;
    int sz = cn.sz;
    cn.apply_spectral_bias();
    cn.inverse_fft();
    std::stringstream lst_string;
    std::ofstream list_file("bitmaplst.py");

    float mag[cn.sz * cn.sz];
    for (int i = 0; i < sz;i++) {
        for (int j = 0; j < sz; j++) {
            mag[i * sz + j] = std::sqrt( sqre(cn.noise[i].real[j]) + sqre(cn.noise[i].imag[j]))/(standard*2);
        }
    }

    float hue[cn.sz * cn.sz];
    for (int i = 0; i < sz;i++) {
        for (int j = 0; j < sz; j++) {
            hue[i * sz + j] = ( std::atan(cn.noise[i].imag[j] / cn.noise[i].real[j])+ pi ) /(2*pi);
        }
    }


    float out_mag[cn.sz * cn.sz];
    FFTPack<cn.sz>* out = cn.fft.output();

    for (int i = 0; i < sz;i++) {
        for (int j = 0; j < sz; j++) {
            mag[i * sz + j] = std::sqrt(sqre(out[i].real[j]) + sqre(out[i].imag[j])) / (standard * 5);
        }
    }


    write_lst_to(mag, cn.sz, cn.sz, 1, lst_string, "inten");
    write_lst_to(cn.spectral_bias,sz,sz,1,lst_string,"spec");
    write_lst_to(hue, cn.sz, cn.sz, 1, lst_string, "hue");

    write_var_to(mean, lst_string, "mean");
    write_var_to(sz, lst_string, "sz");

    write_var_to(standard, lst_string, "standard");
    
    list_file << lst_string.rdbuf();

    list_file.close();
    system("py bitmap.py");
    
}
