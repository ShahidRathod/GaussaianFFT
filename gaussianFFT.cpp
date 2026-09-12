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


template <typename T>
struct Complex {
    T real;
    T imag;

    Complex operator+(Complex& c) {
        return { real + c.real,imag + c.imag };
    }


    Complex operator-(Complex& c) {
        return { real + c.real,imag + c.imag };
    }
};

template <typename T>
struct ComplexRef {
    T* real;
    T* imag;

    void operator=(ComplexRef<T> val) {
        *real = *val.real;
        *imag = *val.imag;
    }
    void operator=(Complex val) {
        *real = val.real;
        *imag = val.imag;
    }


    Complex operator+(Complex& c) {
        return { real + c.real,imag + c.imag };
    }

    Complex operator-(Complex& c) {
        return { real + c.real,imag + c.imag };
    }
};


template <typename T, int sz>
struct ComplexArray {
    T real[sz], imag[sz];

    ComplexRef<T> operator[](int i) {
        return ComplexRef{ real + i,imag + i };
    }


};



template <int n> using ComplexArrayFloat = ComplexArray<float, n>;
using ComplexT = ComplexArrayFloat<1>;


constexpr float pi = std::numbers::pi_v<float>;
constexpr int arrsz(int k, int d) { return 1 << (k * d); }
inline int pow2(int n) { return 1 << n; };



inline int sqre(int x) { return x * x; };

template<int N, int n>
struct FFTPack
{
    constexpr static int nby2 = N / 2, nby4 = n / 4, Nby4 = N / 4;
    inline static float omega[N +Nby4];

    ComplexArrayFloat<n> write;
    FFTPack<N,n / 2> even, odd;
    constexpr static void make_omega() {

        for (int i = 0;i < N; i++)
            omega[i] = cos(2 * pi * i / N);

        for (int i = 0;i < Nby4; i++)
            omega[N + i] = omega[i];

    }

    inline void evaluate() {

        for (int i = 0; i < nby4;i++) {
            float cos = omega[i];
            float sin = omega[i + nby4];
            float res_sin = sin * even.imag[i];
            float res_cos = cos * even.real[i];

            float even_real = res_sin - res_cos;
            float even_imag = res_sin + res_cos;

            even.real[i] = even_real;
            even.imag[i] = even_imag;

            write.real[i] = odd.real[i] + even.real[i];
            write.imag[i] = odd.imag[i] + even.imag[i];


            write.real[i + nby4] = odd.real[i] - even.real[i];
            write.imag[i + nby4] = odd.imag[i] - even.real[i];
        }

    }
};


struct Indx { int i, j; };

template <int sz>
struct inverse2DFFT {
    static constexpr int sz_sq = sz*sz;
    using ComplexArrayTsz = ComplexArrayFloat<sz>;

    ComplexArrayTsz input[sz], fx[sz], fy[sz], fft_buffer_last;
    ComplexT* input = nullptr;
    inverse2DFFT() = default;

    ComplexArrayTsz* output() { return fy; }

    ComplexArrayTsz* input() {
        return input
    }

    template <int n>
    void fft_hlpr(int s, FFTPack<sz, n>* pack) {

        if constexpr (n == 1) {
            *pack = input[s];
            return;
        }
        constexpr int nby2 = n / 2;
        fft_hlpr<nby2>(s + 1, pack->even);
        fft_hlpr<nby2>(s, pack->odd);
        pack->evaluate();
    }

    template <int n>
    void fft(int s, ComplexArrayTsz* inp, ComplexArrayTsz* out) {
        fft_hlpr<nby2>(s + 1, (FFTPack<sz, n>*)inp);
    }


    void eval_fft(ComplexArrayTsz* inp = input,ComplexArrayTsz* out = fx) 
        for (int i = 0; i < 2 * sz;i++) {
            fft_hlpr<sz>(0, out[i]);
        }
    }

};

template <int k>
struct ComplexNoise {
    static constexpr int sz = 1 << k;
    static constexpr int sz_sq = 1 << (2 * k);

    ComplexArrayFloat<sz>* noise;
    float spectral_bias[sz_sq];
    inverse2DFFT<sz> fft;
    std::random_device seed_gen;
    std::normal_distribution<float> normal{ 0, 1 };
    int seed;

    void set_seed(int val) { seed = val; }
    void gen_seed() {
        seed = seed_gen();
    }

    void inverse_fft() {
        fft.input;
        fft.eval_fft();
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
                int indx = i * sz + j;
                int is = sz - i - 1;
                int js = sz - j - 1;
                float val = std::exp(-(sqre(z - i - 1) + sqre(z - j - 1)) * 0.1);
                set_vals(spectral_bias, val, { {i,j} , {is,j},{i,js},{is,js} });
                ComplexT c(normal_gen(), normal_gen());
                noise[indx] = noise[(i + z) * sz + (j + z)] = c;
                noise[indx + z] = noise[(i + z) * sz + j] = ComplexT{ *c.real,-*c.imag };
            }
        }

        set_vals(noise, ComplexT{}, { {z, z}, {z ,sz - 1}, {0,0} });
    }

    void apply_spectral_bias() {
        for (int i = 0; i < sz_sq; i++) {
            *noise[i].real /= spectral_bias[i];
            *noise[i].imag /= spectral_bias[i];
        }
    }

    ComplexNoise() {
        noise = fft.input();
        gen_seed();
        init_arrays();
    }
};



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

int main() {

    ComplexNoise<5> cn;
    int  z = cn.sz / 2;

    cn.inverse_fft();


    std::stringstream lst_string;
    std::ofstream list_file("bitmaplst.py");


    write_lst_to(cn.spectral_bias, cn.sz, cn.sz, 1, lst_string, "lst");
    write_lst_to((float*)cn.fft.output(), cn.sz, cn.sz, 2, lst_string, "lst2");

    list_file << lst_string.rdbuf();

    list_file.close();
    system("py bitmap.py");


}
