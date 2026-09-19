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
#include <numbers>


template <typename T, int sz>
struct ComplexArray {
    T real[sz], imag[sz];
};

std::stringstream content;

template <int n> using ComplexArrayFloat = ComplexArray<float, n>;
using ComplexT = ComplexArrayFloat<1>;

constexpr float pi = std::numbers::pi_v<float>;
constexpr float root2f = std::numbers::sqrt2_v<float>;

constexpr int arrsz(int k, int d) { return 1 << (k * d); }
inline int pow2(int n) { return 1 << n; };

template <typename T>
inline T sqre(T x) { return x * x; };


float mag_f(float a, float b) { return std::sqrt(sqre(a) + sqre(b)); }
float real(float a, float b) { return 0 * a + 1 * std::abs(a); }
float imag(float a, float b) { return std::abs(a); }
float hue_func(float a, float b) { return (std::atan(b / a) + pi) / (2 * pi); }


template <int sz>
void make_arr(ComplexArrayFloat<sz>* input, float* arr, float (*func)(float, float)) {
    for (int i = 0; i < sz; i++) {
        for (int j = 0; j < sz; j++) {
            arr[i * sz + j] = func(input[i].real[j], input[i].imag[j]);
        }
    }
}

template<typename T>
void write_var_to(T& var, std::stringstream& stream, const char* lst_name) {
    stream << lst_name << " = ";
    stream << var << "\n";
}


void write_plain_arr(float* arr, int sz, int stride, std::stringstream& lst_string) {
    lst_string << "[";
    for (int j = 0; j < sz; j++) {
        float lst_val = arr[j * stride];
        lst_string << lst_val;
        if (j != sz - 1) lst_string << ",";
    }
    lst_string << "]\n";
}

void write_plain_lst(float* arr, int sz, int stride, std::stringstream& lst_string, const char* lst_name, int k = -1) {
    lst_string << lst_name;
    if (k != -1) lst_string << k;
    lst_string << " = [\n";
    write_plain_arr(arr, sz, stride, lst_string);
}

void write_lst_to(float* arr, int sz, int d, int stride, std::stringstream& lst_string, const char* lst_name) {
    lst_string << lst_name;
    lst_string << " = [\n";
    for (int i = 0; i < sz; i++) {
        write_plain_arr(arr + d * i, sz, stride, lst_string);
        if (i != sz - 1) lst_string << ",";
    }
    lst_string << "]\n";
}


template<int N, int n>
struct OmegaTabel {
    inline static float omega[n + n / 4];
    inline static OmegaTabel<N, n / 2> next;

    static float* get_minus_sin() { return (omega + n / 4); }
    static float* get_cos() { return omega; }

    static void make_omega() {
        for (int i = 0; i < n; i++)
            omega[i] = std::cos(2 * pi * i / n);

        for (int i = 0; i < n / 4; i++)
            omega[n + i] = omega[i];

        if constexpr (n > 2)
            next.make_omega();
    }
};

template<int N>
struct OmegaTabel<N, 2> {
    inline static float minus_sin[2];
    inline static float cos[2];

    static float* get_minus_sin() { return minus_sin; }
    static float* get_cos() { return cos; }

    static void make_omega() {
        cos[0] = 1;
        cos[1] = -1;
        minus_sin[0] = minus_sin[1] = 0;
    }
};



template<int n>
struct FFTPack : ComplexArrayFloat<n>
{
    static constexpr int nby2 = n / 2;
    static FFTPack<n / 2> even, odd;

    template<int N>
    inline void butterfly(int s, FFTPack<N>& input, int sign) {
        even.butterfly(s, input, sign);
        odd.butterfly(s + (N / n), input, sign);

        float* cos = OmegaTabel<N, n>::get_cos();
        float* minus_sin = OmegaTabel<N, n>::get_minus_sin();

        for (int i = 0; i < nby2; i++) {
            float even_real = cos[i] * even.real[i] + sign * (-minus_sin[i]) * even.imag[i];
            float even_imag = cos[i] * even.imag[i] - sign * (-minus_sin[i]) * even.real[i];

            this->real[i] = odd.real[i] + even_real;
            this->imag[i] = odd.imag[i] + even_imag;

            this->real[i + nby2] = odd.real[i] - even_real;
            this->imag[i + nby2] = odd.imag[i] - even_imag;
        }
    }

    template<int N>
    inline void invrs_fft(FFTPack<N>& input) {
        butterfly(0, input, -1); // -minus_sin = sin
    }

    template<int N>
    inline void fft(FFTPack<N>& input) {
        butterfly(0, input, 1); // minus_sin = -sin
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
    inline void butterfly(int s, FFTPack<N>& input, int sign) {
        real[0] = input.real[s];
        imag[0] = input.imag[s];
    }
};

template <typename T, int sz>
void FFT(ComplexArray<T,sz>& input , ComplexArray<T,sz>& output) {

}
struct Indx { int i, j; };

template <int sz>
struct FFT2D {
    static constexpr int sz_sq = sz * sz;

    FFTPack<sz> input[sz * sz], fx[sz * sz], fy[sz * sz];

    FFT2D() {}

    FFTPack<sz>* output() { return &fy[0]; }
    FFTPack<sz>* x_arr() { return &fx[0]; }
    FFTPack<sz>* y_arr() { return &fy[0]; }

    void transpose(ComplexArrayFloat<sz>* arr) {
        for (int i = 0; i < sz; i++) {
            for (int j = i; j < sz; j++) {
                float temp_real = arr[i].real[j];
                arr[i].real[j] = arr[j].real[i];
                arr[j].real[i] = temp_real;

                float temp_imag = arr[i].imag[j];
                arr[i].imag[j] = arr[j].imag[i];
                arr[j].imag[i] = temp_imag;
            }
        }
    }


    void fft_of(FFTPack<sz>(&butrfly_inp)[sz * sz],
        FFTPack<sz>(&butrfly_out)[sz * sz],
        int sign) {
        for (int i = 0; i < sz; i++)
            butrfly_out[i].template butterfly<sz>(0, butrfly_inp[i], sign);
        transpose(&butrfly_out[0]);
    }

    void make_out_in() {
        memcpy(input, fy, sizeof(input));
    }

    // Forward FFT: sign = -1
    void fft() {
        fft_of(input, fx, 1);
        fft_of(fx, fy, 1);
    }

    // Inverse FFT: sign = +1
    void inverse_fft() {
        fft_of(input, fx, -1);
        fft_of(fx, fy, -1);
    }
  
};

constexpr double spec_radii = 0.0001;

template <int k>
struct ComplexNoise {
    static constexpr int sz = 1 << k;
    static constexpr int sz_sq = sz * sz;
    static constexpr float standard = 1 / (root2f * sz);

    ComplexArrayFloat<sz>* noise;
    float ref_landscp[sz_sq];
    float output[sz_sq];
    float spectral_bias[sz_sq];
    FFT2D<sz> fft;
    std::random_device seed_gen;
    std::normal_distribution<float> normal{ 0.f, standard };
    std::default_random_engine engine;
    int seed;
    void set_seed(int val) { seed = val; engine.seed(seed); }
    void gen_seed() { set_seed(seed_gen()); }

    void inverse_fft() { fft.inverse_fft(); }

    inline int index(int i, int j) { return i + sz * j; }

    template<typename T>
    void set_vals(T* arr, T val, std::initializer_list<Indx> lst) {
        for (auto& x : lst) arr[index(x.i, x.j)] = val;
    }
    
    void init_noise() {
        gen_seed();
        int z = sz / 2;
        for (int i = 0; i < z; i++) {
            for (int j = 0; j < z; j++) {
                int is = sz - i - 1;
                int js = sz - j - 1;

                noise[i].real[j] = noise[is].real[js] = normal(engine);
                noise[is].real[j] = noise[i].real[js] = normal(engine);

                float conj1 = normal(engine);
                noise[i].imag[j] = conj1;
                noise[is].imag[js] = -conj1;

                float conj2 = normal(engine);
                noise[is].imag[j] = conj2;
                noise[i].imag[js] = -conj2;
            }
        }

        noise[0].imag[0] = 0;
        noise[z].imag[z] = 0;
        noise[z].imag[sz - 1] = 0;
        noise[sz - 1].imag[z] = 0;
    }

    void init_spectral_bias() {
        int z = sz / 2;
        for (int i = 0; i < z; i++) {
            for (int j = 0; j < z; j++) {
                int is = sz - i - 1;
                int js = sz - j - 1;

                double dis = sqre(z - i - 1) + sqre(z - j - 1);
                double scaled_dis = dis / (double)(z * z);

                double val = std::pow(1 + std::pow(scaled_dis / spec_radii, 1), -1);
                //double val = (i == j == z - 1) ? 1 : 0.1;

                set_vals(spectral_bias, (float)val,
                    { {i,j}, {is,j}, {i,js}, {is,js} });
            }
        }
    }

    float spectral_power() {
        double power = 0;
        for (int i = 0; i < sz_sq; i++) power += sqre((double)spectral_bias[i]);

        double res = pi * std::sqrt(power) / sz;
        //std::cout << "power: " << res << "\n";
        return res;
    }

    void apply_scaling(float* arr) {
        float scaling = spectral_power();
        for (int i = 0; i < sz_sq; i++) arr[i] /= scaling;
    }

    void apply_spectral_bias() {
        for (int i = 0; i < sz; i++) {
            for (int j = 0; j < sz; j++) {
                noise[i].real[j] *= spectral_bias[index(i, j)];
                noise[i].imag[j] *= spectral_bias[index(i, j)];
            }
        }
    }

    void colored_noise(float* hue, float* inten) {
        make_arr(noise, hue, hue_func);
        make_arr(noise, inten, mag_f);
        apply_scaling(inten);
    }

    void grayscale_noise(float* hue, float* inten) {
        make_arr(noise, hue, hue_func);
        make_arr(noise, inten, mag_f);
        apply_scaling(inten);
    }
    void output_grayscale(float* mag_arr) {
        make_arr(&(fft.output())[0], mag_arr, mag_f);
        apply_scaling(mag_arr);
    }
    
    void output_colored(float* hue,float* inten) {
        make_arr(&(fft.output())[0],hue, hue_func);
        make_arr(&(fft.output())[0], inten, mag_f);
        apply_scaling(inten);
    }
    
    void new_landscp() {
        for (int i = 0; i < sz; i++) {
            memcpy(&fft.fx[i].real, ref_landscp + i, sizeof(float) * sz);
        }

        fft.fft();
        make_arr(&(fft.output())[0],spectral_bias,mag_f);
        init_noise();
        apply_spectral_bias();
        fft.inverse_fft();
    }

    ComplexNoise() {
        OmegaTabel<sz, sz>::make_omega();
        noise = fft.input;
        init_noise();
        init_spectral_bias();
    }
};


constexpr long long int grid_pow = 8;
constexpr long int grid_len = 1 << 2 * grid_pow;


static ComplexNoise<grid_pow> cn;
constexpr int sz = cn.sz;



int main() {
    float out_mag[sz * sz];
    std::stringstream lst_string;
    std::ofstream list_file("lists.py");
        write_var_to(sz, lst_string, "sz");

        cn.apply_spectral_bias();
        cn.fft.inverse_fft();
        cn.output_grayscale(out_mag);
        write_lst_to(out_mag, sz, sz, 1, lst_string, "out");
        write_lst_to(cn.spectral_bias, sz, sz, 1, lst_string, "spec");

        cn.fft.make_out_in();
        cn.fft.fft();
        cn.output_grayscale(out_mag);
        write_lst_to(out_mag, sz, sz, 1, lst_string, "spec2");

        memcpy(cn.spectral_bias, out_mag, sizeof(out_mag));

        cn.init_noise();
        cn.apply_spectral_bias();
        cn.fft.inverse_fft();
        cn.output_grayscale(out_mag);
        write_lst_to(out_mag, sz, sz, 1, lst_string, "out2");

        list_file << lst_string.rdbuf();
        list_file.close();
        system("py bitmap.py");

}
