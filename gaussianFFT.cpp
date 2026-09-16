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

    static float* sin() { return (omega + N / 4); }
    static float* cos() { return omega; }

    static void make_omega() {
        for (int i = 0;i < N; i++)
            omega[i] = std::cos(2 * pi * i / N);

        for (int i = 0;i < N / 4; i++)
            omega[N + i] = omega[i];

        if constexpr (N > 8)
            next.make_omega();
    }
};


template<>
struct OmegaTabel<1> {};


template<int n>
struct FFTPack : ComplexArrayFloat<n>
{
    static constexpr int nby2 = n / 2;
    static FFTPack<n / 2> even, odd;



    template<int N>
    inline void butterfly(int s, FFTPack<N>& input) {
        even.butterfly(s, input);
        odd.butterfly(s + 1, input);

        float* cos = OmegaTabel<n>::sin();
        float* sin = OmegaTabel<n>::cos();


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


template <int sz> 
struct inverse2DFFT {
    static constexpr int sz_sq = sz * sz;

    FFTPack<sz> input[sz*sz], fx[sz*sz],fy[sz*sz];

    inverse2DFFT() {}

    FFTPack<sz>* output() { return &fy[0]; }
    FFTPack<sz>* x_arr() { return &fx[0]; }
    FFTPack<sz>* y_arr() { return &fy[0]; }

    void fft_of(FFTPack<sz>(& butrfly_inp)[sz * sz], FFTPack<sz>(&butrfly_out)[sz * sz]) {

        for (int i = 0; i < sz;i++)
            butrfly_out[i].template butterfly<sz>(0, butrfly_inp[i]);

        //transpose
       /* ComplexArrayFloat<sz>* arr = &butrfly_out[0];
        for (int i = 0; i < sz; i++) {
            for (int j = i; j < sz; j++) {

                float temp_real = arr[i].real[j];
                float temp_imag = arr[i].imag[j];

                arr[i].real[j] = arr[j].real[i];
                arr[i].imag[j] = arr[j].imag[i];

                arr[j].real[i] = temp_real;
                arr[j].imag[i] = temp_real;

            }
        }*/ 
    }
    
    void fft() {
        fft_of(input,fx);
        fft_of(fx,fy);
    }
};





template <int k>
struct ComplexNoise {
    static constexpr int sz = 1 << k;
    static constexpr int sz_sq = sz*sz;
    static constexpr float standard =  1/ (float)(sz * sz);
    ComplexArrayFloat<sz>* noise;
    float spectral_bias[sz_sq];
    inverse2DFFT<sz> fft;
    std::random_device seed_gen;
    std::normal_distribution<float> normal{ 0.f, standard};
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

                float val = std::exp(-(sqre(z - i - 1) + sqre(z - j - 1)) * 0.0001);
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
                noise[i].real[j] = 0*std::sin(2*pi * j / sz) + 0*std::sin(2*pi * i / sz);
                noise[i].imag[j] = 1*std::sin(2*pi * j / sz);
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


void write_plain_arr(float* arr, int sz,int stride, std::stringstream& lst_string) {

    lst_string << "[";
    for (int j = 0; j < sz; j++)
    {
        float lst_val = arr[j * stride];
        lst_string << lst_val;
        if (j != sz - 1) lst_string << ",";
    }
    lst_string << "]\n";
}

void write_plain_lst(float* arr, int sz, int stride, std::stringstream& lst_string, const char* lst_name) {
    lst_string << lst_name;
    lst_string << " = [\n";
    write_plain_arr(arr , sz, stride, lst_string);
}

void write_lst_to(float* arr, int sz, int d, int stride, std::stringstream& lst_string, const char* lst_name) {

    lst_string << lst_name;
    lst_string << " = [\n";
    for (int i = 0; i < sz; i++)
    {
        write_plain_arr(arr+d*i,sz,stride,lst_string);
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

template <int sz>
void make_arr(ComplexArrayFloat<sz>* input, float* arr, float (*func)(float,float)) {
    for (int i = 0; i < sz;i++) {
        for (int j = 0; j < sz; j++) {
            arr[i * sz + j] = func(input[i].real[j],input[i].imag[j]);
        }
    }
}

float mag_f(float a , float b) {
    return std::sqrt(sqre(a)+sqre(b));
}

float real(float a, float b) {
    return a;
}

float imag(float a, float b) {
    return a;
}

float hue_func(float a, float b) {
    return (std::atan(b / a) + pi) / (2 * pi);
}

 
#define OmegaTabelWrite(tabelno,iterno) \
lst_string << "omega" << tabelno << " = ";\
write_plain_arr(OmegaTabel<tabelno>::cos(), tabelno, 1, lst_string);\
\


int main() {
    constexpr int grid_len = 4;
    constexpr int grid_len_sq0 = 1<<2*grid_len;

    static ComplexNoise<grid_len> cn;
    int sz = cn.sz;
    cn.apply_spectral_bias();
    cn.inverse_fft();

    std::stringstream lst_string;
    std::ofstream list_file("bitmaplst.py");

    float scaling = cn.standard;

    float mag[cn.sz * cn.sz];
    make_arr(cn.noise, mag, real);

  
    float hue[cn.sz * cn.sz];
    make_arr(cn.noise, hue, hue_func);

   
    float out_mag[cn.sz * cn.sz];

    ComplexArrayFloat<cn.sz>* out = cn.noise;
    make_arr(cn.fft.x_arr(), out_mag, mag_f);
    //write_lst_to(mag, cn.sz, cn.sz, 1, lst_string, "input");
    //write_lst_to(cn.spectral_bias,sz,sz,1,lst_string,"spec");
    //write_lst_to(hue, cn.sz, cn.sz, 1, lst_string, "hue");

    //write_lst_to(out_mag, cn.sz, cn.sz, 1, lst_string, "out");
    //write_var_to(sz, lst_string, "sz");
    //write_var_to(scaling, lst_string, "scaling");

    int iterno = 0;
    OmegaTabelWrite(grid_len_sq0,iterno)
    list_file << lst_string.rdbuf();

    list_file.close();
    //system("py bitmap.py");

}
