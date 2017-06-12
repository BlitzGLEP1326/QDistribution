#include <QObject>
#include <iostream>
//
#include "qdistribution.h"

using namespace std;

static double stairWidth[257], stairHeight[256];
// константы для нормального распределения
const double x1_Normal = 3.6541528853610088;
const double A_Normal = 4.92867323399e-3; /// area under rectangle
// константы для экспоненциального распределения
const double x1_Exp = 7.69711747013104972;
const double A_Exp = 3.9496598225815571993e-3;

#define M_SQRT3 1.7320508075689
//-----------------------------------------------------------------
//
QDistribution::QDistribution(QObject *parent) : QObject(parent)
{

}
//-----------------------------------------------------------------
//
QDistribution::~QDistribution()
{

}
double QDistribution::UniformDistibRNDGenerator(double a, double b)
{
    try{
        return double(a + qrand() * (b - a) / RAND_MAX);
    }catch(...){}
}
//-----------------------------------------------------------------
//
void QDistribution::setupNormalTables()
{
    try{
        // coordinates of the implicit rectangle in base layer
        stairHeight[0] = exp(-.5 * x1_Normal * x1_Normal);
        stairWidth[0] = A_Normal / stairHeight[0];
        // implicit value for the top layer
        stairWidth[256] = 0;
        for (unsigned i = 1; i <= 255; ++i)
        {
            // such x_i that f(x_i) = y_{i-1}
            stairWidth[i] = sqrt(-2 * log(stairHeight[i - 1]));
            stairHeight[i] = stairHeight[i - 1] + A_Normal / stairWidth[i];
        }
    }catch(...){}
}
//-----------------------------------------------------------------
//
double QDistribution::NormalZiggurat() {
    try{
        setupNormalTables();
        int iter = 0;
        do {
            unsigned long long B = -1 + qrand() * 2;
            int stairId = B & 255;
            double x = UniformDistibRNDGenerator(0, stairWidth[stairId]); // get horizontal coordinate
            if (x < stairWidth[stairId + 1])
                return ((signed)B > 0) ? x : -x;
            if (stairId == 0) // handle the base layer
            {
                static double z = -1;
                double y;
                if (z > 0) // we don't have to generate another exponential variable as we already have one
                {
                    x = ExponentialDistibRNDGenerator(x1_Normal);
                    z -= 0.5 * x * x;
                }
                if (z <= 0) // if previous generation wasn't successful
                {
                    do {
                        x = ExponentialDistibRNDGenerator(x1_Normal);
                        y = ExponentialDistibRNDGenerator(1);
                        z = y - 0.5 * x * x; // we storage this value as after acceptance it becomes exponentially distributed
                    } while (z <= 0);
                }
                x += x1_Normal;
                return ((signed)B > 0) ? x : -x;
            }
            // handle the wedges of other stairs
            if (UniformDistibRNDGenerator(stairHeight[stairId - 1], stairHeight[stairId]) < exp(-.5 * x * x))
                return ((signed)B > 0) ? x : -x;
        } while (++iter <= 1e9); /// one billion should be enough
        return NAN; /// fail due to some error
    }catch(...){}
}
//-----------------------------------------------------------------
//
double QDistribution::NormalDistibRNDGenerator(double mu, double sigma) {
    try{
        return double(mu + NormalZiggurat() * sigma);
    }catch(...){}
}
//-----------------------------------------------------------------
//
void RepresentationSetCreatorDialog::setupExpTables() {
    // coordinates of the implicit rectangle in base layer
    stairHeight[0] = exp(-x1_Exp);
    stairWidth[0] = A_Exp / stairHeight[0];
    // implicit value for the top layer
    stairWidth[256] = 0;
    for (unsigned i = 1; i <= 255; ++i)
    {
        // such x_i that f(x_i) = y_{i-1}
        stairWidth[i] = -log(stairHeight[i - 1]);
        stairHeight[i] = stairHeight[i - 1] + A_Exp / stairWidth[i];
    }
}
//-----------------------------------------------------------------
//
double QDistribution::ExpZiggurat() {
    try{
        setupExpTables();
        int iter = 0;
        do {
            int stairId = qrand() & 255;
            double x = UniformDistibRNDGenerator(0, stairWidth[stairId]); // get horizontal coordinate
            if (x < stairWidth[stairId + 1]) /// if we are under the upper stair - accept
                return x;
            if (stairId == 0) // if we catch the tail
                return x1_Exp + ExpZiggurat();
            if (UniformDistibRNDGenerator(stairHeight[stairId - 1], stairHeight[stairId]) < exp(-x)) // if we are under the curve - accept
                return x;
            // rejection - go back
        } while (++iter <= 1e9); // one billion should be enough to be sure there is a bug
        return NAN; // fail due to some error
    }catch(...){}
}
//-----------------------------------------------------------------
//
double QDistribution::ExponentialDistibRNDGenerator(double rate) {
    try{
        return double(ExpZiggurat() / rate);
    }catch(...){}
}
//-----------------------------------------------------------------
/*Если сложить две случайные величины с гамма-распределением с параметрами k1 и k2,
 * то получится случайная величина с гамма-распределением и с параметром k1+k2.
 * Еще одно свойство — если theta = k = 1, то легко проверить, что распределение будет экспоненциальным.
 * Поэтому, если k целое — то можно просто просуммировать k случайных величин со стандартным экспоненциальным распределением.
 */
double QDistribution::GA1(int k) {
    double x = 0;
    for (int i = 0; i < k; ++i)
        x += ExponentialDistibRNDGenerator(1);
    return x;
}
//-----------------------------------------------------------------
/*Если k не целое, но 2k — целое, то можно вместо одной из экспоненциальных
 * случайных величин в сумме использовать половину квадрата нормальной величины.
 * Почему так возможно, станет ясно позднее.
 */
double QDistribution::GA2(double k) {
    double x = NormalDistibRNDGenerator(0, 1);
    x *= 0.5 * x;
    for (int i = 1; i < k; ++i)
        x += ExponentialDistibRNDGenerator(1);
    return x;
}
//-----------------------------------------------------------------
/* Генерируем стандартную экспоненциально распределенную величину Е и величину U,
 * равномерно распределенную от 0 до 1 + k / e. Если U <= 1, то переходим к шагу 2. Иначе, к шагу 3.
 * Задаем x = U1/k. Если x <= E, то возвращаем x, иначе возвращаемся назад на шаг 1.
 * Задаем x = -ln((1 — U) / k + 1 / e). Если (1 — k) * ln(x) <= E, то возвращаем x, иначе возвращаемся на шаг 1.
 */
double QDistribution::GammaDistibRNDGenerator(double k) {
    // Assume that k < 1
    double x = 0;
    int iter = 0;
    do {
        // M_E is base of natural logarithm
        double U = UniformDistibRNDGenerator(0, 1 + k / M_E);
        double W = ExponentialDistibRNDGenerator(1);
        if (U <= 1)
        {
            x = pow(U, 1.0 / k);
            if (x <= W)return x;
        }
        else
        {
            x = -log((1 - U) / k + 1.0 / M_E);
            if ((1 - k) * log(x) <= W)return x;
        }
    } while (++iter < 1e9); // excessive maximum number of rejections
    return NAN; // shouldn't end up here
}
//-----------------------------------------------------------------
/*Насколько мне известно, автор этого алгоритма, профессор Университета Северной Каролины Дж. С. Фишман,
 * не стал опубликовывать свое достижение. Сам алгоритм работает для k > 1, однако среднее время его выполнения
 * возрастает пропорционально sqrt(k), поэтому он эффективен только для k < 3. Алгоритм:
 * Генерируем две независимых случайных величины E1 и E2 со стандартным экспоненциальным распределением.
 * Если E2 < (k — 1) * (E1 — ln(E1) — 1), то возвращаемся назад на шаг 1.
 * Возвращаем x = k * E1
 */
double QDistribution::GF(double k) {
    // Assume that 1 < k < 3
    double E1, E2;
    do {
        E1 = ExponentialDistibRNDGenerator(1);
        E2 = ExponentialDistibRNDGenerator(1);
    } while (E2 < (k - 1) * (E1 - log(E1) - 1));
    return double(k * E1);
}
//-----------------------------------------------------------------
//
void QDistribution::setupConstants(double k) {
    try
    {

    }catch(...){}
}
//-----------------------------------------------------------------
/*
 * Алгоритм Дитера и Аренса основан на асимптотической нормальности
 * распределения при увеличении параметра, и поэтому более быстр для больших k.
 * Он работает для k > 2.533, однако не столь эффективен как алгоритм Фишера для k < 3.
 * Для начала нужно задать некоторые константы (зависящие только от k).
 */
double QDistribution::GO(double k) {
    try{
        // Assume that k > 3
        double x = 0;
        int iter = 0;
        double m = k - 1;
        double s_2 = sqrt(8.0 * k / 3) + k;
        double s = sqrt(s_2);
        double d = M_SQRT2 * M_SQRT3 * s_2;
        double b = d + m;
        double w = s_2 / (m - 1);
        double v = (s_2 + s_2) / (m * sqrt(k));
        double c = b + log(s * d / b) - m - m - 3.7203285;
        //
        do {
            double U = UniformDistibRNDGenerator(0, 1);
            if (U <= 0.0095722652) {
                double E1 = ExponentialDistibRNDGenerator(1);
                double E2 = ExponentialDistibRNDGenerator(1);
                x = b * (1 + E1 / d);
                if (m * (x / b - log(x / m)) + c <= E2)
                    return x;
            }
            else {
                double N;
                do {
                    N = NormalDistibRNDGenerator(0, 1);
                    x = s * N + m; // ~ Normal(m, s)
                } while (x < 0 || x > b);
                U = UniformDistibRNDGenerator(0, 1);
                double S = 0.5 * N * N;
                if (N > 0) {
                    if (U < 1 - w * S)
                        return x;
                }
                else if (U < 1 + S * (v * N - w))
                    return x;
                if (log(U) < m * log(x / m) + m - x + S)
                    return x;
            }
        } while (++iter < 1e9);
        return NAN; // shouldn't end up here;
    }catch(...){}
}
//-----------------------------------------------------------------
//
double QDistribution::Cauchy(double x0, double gamma)
{
    try
    {
        double x, y;
        do {
            x = UniformDistibRNDGenerator(-1,1);
            y = UniformDistibRNDGenerator(-1,1);
        } while (x * x + y * y > 1.0 || y == 0.0);
        return double(x0 + gamma * x / y);
    }catch(...){}
}
//-----------------------------------------------------------------
//
double QDistribution::Laplace(double mu, double b) {
    try{
        double E = ExponentialDistibRNDGenerator(1.0 / b);
        return double(mu + (((signed)qrand() > 0) ? E : -E));
    }catch(...){}
}
//-----------------------------------------------------------------
//
double QDistribution::Levy(double mu, double c) {
    try{
        double N = NormalDistibRNDGenerator(0, 1);
        return double(mu + c / (N * N));
    }catch(...){}
}
//-----------------------------------------------------------------
//
double QDistribution::ChiSquared(int k) {
    // ~ Gamma(k / 2, 2)
    if (k >= 10) // too big parameter
        return GO(0.5 * k);
    double x = ((k & 1) ? GA2(0.5 * k) : GA1(k >> 1));
    return double(x + x);
}
//-----------------------------------------------------------------
//
double QDistribution::LogNormal(double mu, double sigma) {
    try{
        return exp(NormalDistibRNDGenerator(mu, sigma));
    }catch(...){}
}
//-----------------------------------------------------------------
//
double QDistribution::Logistic(double mu, double s) {
    try{
        return double(mu + s * log(1.0 / UniformDistibRNDGenerator(0, 1) - 1));
    }catch(...){}
}
//-----------------------------------------------------------------
//
double QDistribution::Erlang(int k, double l) {
    try{
        return double(GA1(k) / l);
    }catch(...){}
}
//-----------------------------------------------------------------
//
double QDistribution::Weibull(double l, double k) {
    try{
        return double(l * pow(ExponentialDistibRNDGenerator(1), 1.0 / k));
    }catch(...){}
}
//-----------------------------------------------------------------
//
double QDistribution::Rayleigh(double sigma) {
    try{
        return double(sigma * sqrt(ExponentialDistibRNDGenerator(0.5)));
    }catch(...){}
}
//-----------------------------------------------------------------
//
double QDistribution::Pareto(double xm, double alpha) {
    try{
        return double(xm / pow(UniformDistibRNDGenerator(0, 1), 1.0 / alpha));
    }catch(...){}
}
//-----------------------------------------------------------------
//
double QDistribution::StudentT(int v) {
    try{
        if (v == 1)
            return Cauchy(0, 1);
        return double(NormalDistibRNDGenerator(0, 1) / sqrt(ChiSquared(v) / v));
    }catch(...){}
}
//-----------------------------------------------------------------
//
double QDistribution::FisherSnedecor(int d1, int d2) {
    try{
        double numerator = d2 * ChiSquared(d1);
        double denominator = d1 * ChiSquared(d2);
        return double(numerator / denominator);
    }catch(...){}
}
//-----------------------------------------------------------------
//
double QDistribution::Beta(double a, double b) {
    try{
        double x = gamma(a);
        return double(x / (x + gamma(b)));
    }catch(...){}
}
