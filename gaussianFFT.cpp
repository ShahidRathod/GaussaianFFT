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

inline int sqre(int x) { return x * x; };

using ComplexT = std::complex<float>;
constexpr float pi = std::numbers::pi_v<float>;


constexpr int arrsz(int k, int d) { return 1 << (k * d); }

inline int pow2(int n) { return 1 << n; };

struct Indx { int i, j; };

template <int k>
struct inverseFFT {
    static constexpr int sz = 1 << k;
    ComplexT omega[sz];
    ComplexT* output = nullptr;
    ComplexT* input = nullptr;

    inverseFFT() = default;

    inverseFFT(ComplexT* inp, ComplexT* out) {
        input = inp;
        output = out;
        set_omega();
    }

    void set_omega() {
        for (int i = 0; i < sz; i++) {
            omega[i] = std::polar(1.f, 2 * (pi * i) / sz);
        }
    }

    void eval_fft() {
        for (int i = 0; i < sz;i++) fft(sz, 0, output + i * sz);
    }

    void fft(int n, int s, ComplexT* write) {
        if (n == 1) {
            write[0] = input[s];
            return;
        }

        int n2 = n / 2;
        ComplexT* odd = write + n;
        ComplexT* even = odd + n2;
       
        fft(n2, s, odd);
        fft(n2, s + 1, even);

        for (int i = 0; i < n2; i++) {
            ComplexT ei = even[i] * omega[(i * sz / n2) % sz];
            ComplexT oi = odd[i];
            write[i] = oi + ei;
            write[n2 + i] = oi - ei;
        }
    }

};


template <int k>
struct ComplexNoise {
    static constexpr int sz = 1 << k;
    static constexpr int sz_sq = 1 << (2 * k);

    ComplexT noise[sz_sq];
    float spectral_bias[sz_sq];

    ComplexT fx[sz_sq], fy[sz_sq], fft_buffer[3 * sz_sq];

    inverseFFT<k> fftx, ffty;

    std::random_device seed_gen;
    std::normal_distribution<float> normal{ 0, 1 };
    int seed;

    void set_seed(int val) { seed = val; }
    void gen_seed() {
        seed = seed_gen();
    }

    void inverse_fft() {

        fftx.input;
        ffty.input;
        fftx.eval_fft();
        ffty.eval_fft();
    }

    inline int index(int i, int j) { return i + sz * j; }

    template<typename T>
    void set_vals(T* arr, T val, std::initializer_list <Indx> lst) {
        for (auto& x : lst) arr[index(x.i,x.j)] = val;
    }

    void init_arrays() {

        auto normal_gen = std::bind(normal, std::default_random_engine(seed));
        int z = sz / 2;
        for (int i = 0; i < z; i++) {
            for (int j = 0; j < z; j++) {
                int indx = i * sz + j;
                int is = sz - i-1;
                int js = sz - j-1;
                //float val = std::exp(-((z-i)<<1 + (z-j)<<1)*0.1)*20;
                float val  = std::exp(- ( sqre(z - i - 1) + sqre(z - j - 1))*0.1 );
                set_vals(spectral_bias, val, { {i,j} , {is,j},{i,js},{is,js} });
                ComplexT c(normal_gen(), normal_gen());
                noise[indx] = noise[(i + z) * sz + (j + z)] = c;
                noise[indx + z] = noise[(i + z) * sz + j] = std::conj(c);
            }
        }

        set_vals(noise, ComplexT{}, { {z, z}, {z ,sz - 1}, {0,0} });
    }

    ComplexNoise() {
        fftx = inverseFFT<k>{ noise, fx };
        ffty = inverseFFT<k>{ fx,fy };
        gen_seed();
        init_arrays();
    }
};

void write_lst_to(float* arr, int sz , int d,int stride, std::stringstream & lst_string,const char * lst_name){
    
    lst_string << lst_name;
    lst_string << " = [\n";
    for (int i = 0; i < sz; i++)
    {
        lst_string << "[";
        for (int j = 0; j < sz; j++)
        {
            float lst_val = arr[j*stride+ d * i];
            lst_string << lst_val;
            if (j != sz - 1) lst_string << ",";
        }
        lst_string << "]\n";
        if (i != sz - 1) lst_string << ",";
    }
    lst_string << "]\n";

}

int main() {

    ComplexNoise<5> cn;
    int  z = cn.sz / 2;

    cn.inverse_fft();

  

    std::stringstream lst_string;
    std::ofstream list_file("bitmaplst.py");
    

    write_lst_to(cn.spectral_bias, cn.sz, cn.sz, 1, lst_string, "lst");
    write_lst_to((float*)cn.fy, cn.sz, cn.sz, 2, lst_string, "lst2");

    list_file << lst_string.rdbuf();

    list_file.close();
    system("py bitmap.py");


}