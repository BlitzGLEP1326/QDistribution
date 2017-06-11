#ifndef QDISTRIBUTION_H
#define QDISTRIBUTION_H

#include <QObject>

class QDistribution : public QObject
{
    Q_OBJECT
public:
    explicit QDistribution(QObject *parent = 0);
    ~QDistribution();

signals:

public slots:
    double UniformDistibRNDGenerator(double a, double b);
    void setupNormalTables();
    double NormalZiggurat();
    double NormalDistibRNDGenerator(double mu, double sigma);
    void setupExpTables();
    double ExpZiggurat();
    double ExponentialDistibRNDGenerator(double rate);
    double GA1(int k);
    double GA2(double k);
    double GammaDistibRNDGenerator(double k);
    double GF(double k);
    void setupConstants(double k);
    double GO(double k);
    double Cauchy(double x0, double gamma);
    double Laplace(double mu, double b);
    double Levy(double mu, double c);
    double ChiSquared(int k);
    double LogNormal(double mu, double sigma);
    double Logistic(double mu, double s);
    double Erlang(int k, double l);
    double Weibull(double l, double k);
    double Rayleigh(double sigma);
    double Pareto(double xm, double alpha);
    double StudentT(int v);
    double FisherSnedecor(int d1, int d2);
    double Beta(double a, double b);
};

#endif // QDISTRIBUTION_H
