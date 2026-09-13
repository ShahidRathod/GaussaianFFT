#+includ---++e <random>
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



inline int sqre(int x) { return x * x; };

template<int N> 
struct OmegaTabel {
    inline static float omega[N + Nby4];
    bool omega_made = false;
    static void make_omega() {
        if (omega_made) return;
        for (int i = 0;i < N; i++)
            omega[i] = cos(2 * pi * i / N);

        for (int i = 0;i < Nby4; i++)
            omega[N + i] = omega[i];
        omega_made = true;
    }
};


template<int n>
struct FFTPack
{
    ComplexArrayFloat<n> write;
    static FFTPack<n / 2> even, odd;
   
    inline void butterfly() {
        for (int i = 0; i < nby2;i++) {

            float cos = OmegaTabel<N>::omega[Nbyn * i];
            float sin = OmegaTabel<N>::omega[i * Nbyn + Nby4];

            even.real[i] = sin * even.real[i] - cos * even.imag[i];
            even.imag[i] = sin * even.imag[i] + cos * even.real[i];

            write.real[i] = odd.real[i] + even.real[i];
            write.imag[i] = odd.imag[i] + even.imag[i];

            write.real[i + nby2] = odd.real[i] - even.real[i];
            write.imag[i + nby2] = odd.imag[i] - even.real[i];
        }

    }
};


struct FFTPack<1> 
{
    ComplexArrayFloat<1> write;
    inline void butterfly() {}
};

template <typename LstT, int I> struct IndexedLst {
    static LstT item;
    IndexedLst<LstT,I - 1>;
    inline void butterfly() { 
        pack.butterfly(); 
    }
};

template<typename T> struct IndexedLst<T, 0> {};

template <int L , int n>
using IndexedFFTPack = IndexLst<FFTPack<n>,L>;

struct Indx { int i, j; };

template <int sz, int d>
struct inverse2DFFT {
    static constexpr int sz_sq = sz * sz;
    FFTPack<sz> input[sz];
    inverse2DFFT<sz, d - 1> output;
    inverse2DFFT() {}

    template <int n>
    void fft_hlpr(int s, FFTPack<sz, n>& pack) {
        if constexpr (n == 1) {
            pack->write.real[0] = input.write.real[s];
            pack->write.imag[0] = input.write.imag[s];
        }
        else {
            constexpr int nby2 = n / 2;
            fft_hlpr<nby2>(s + 1, pack.even);
            fft_hlpr<nby2>(s, pack.odd);
            pack->evaluate();
        }
    }

    void fft() {
        for (int i = 0; i < 2 * sz;i++) fft_hlpr(0, output[i]);
        output.fft();
    }

};

template<int sz> 
struct inverse2DFFT<sz, 0> {
    void fft() {}
};

template <int k>
struct ComplexNoise {
    static constexpr int sz = 1 << k;
    static constexpr int sz_sq = 1 << (2 * k);

    ComplexArrayFloat<sz>* noise;
    float spectral_bias[sz_sq];
    inverse2DFFT<sz,2> fft;
    std::random_device seed_gen;
    std::normal_distribution<float> normal{ 0, 1 };
    int seed;

    void set_seed(int val) { seed = val; }
    void gen_seed() {
        seed = seed_gen();
    }

    void inverse_fft() {
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

                int is = sz - i - 1;
                int js = sz - j - 1;
                float val = std::exp(-(sqre(z - i - 1) + sqre(z - j - 1)) * 0.1);
                set_vals(spectral_bias, val, { {i,j} , {is,j},{i,js},{is,js} });

                noise[i].real[j] = noise[z - i].real[z - j] =
                    noise[i + z].real[j + z] = noise[i + z].real[z] = normal_gen();

                float val = noise[i + z].imag[i] = noise[i].imag[i + z] = normal_gen();
                noise[i].imag[i] = noise[i + z].imag[i + z] = -val;

            }

        }
        noise[0].real[0] = noise[z].real[z] = noise[z].real[sz - 1] = noise[sz - 1].real[z] = 0;

    }

    void apply_spectral_bias() {
        for (int i = 0; i < sz;i++) {
            for (int j = 0; j < sz; j++) {
                noise[i].real[j] /= spectral_bias[i];
                noise[i].imag[j] /= spectral_bias[i];
            }
        }
    }

    ComplexNoise() {
        FFTCompute<sz,sz>::make_omega();
        noise = fft.input;
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
