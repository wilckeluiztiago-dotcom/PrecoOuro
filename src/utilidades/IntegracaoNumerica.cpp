/*
 * IntegracaoNumerica.cpp
 * Módulo de integração e derivação numérica
 * 
 * Autor: Luiz Tiago Wilcke
 * Projeto: Análise do Preço do Ouro
 */

#include "IntegracaoNumerica.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <functional>

namespace analise_ouro {

// Constantes para quadratura de Gauss-Legendre (5 pontos)
const std::vector<double> PONTOS_GAUSS = {
    -0.9061798459386640,
    -0.5384693101056831,
     0.0,
     0.5384693101056831,
     0.9061798459386640
};

const std::vector<double> PESOS_GAUSS = {
    0.2369268850561891,
    0.4786286704993665,
    0.5688888888888889,
    0.4786286704993665,
    0.2369268850561891
};

// Construtor padrão
IntegracaoNumerica::IntegracaoNumerica() 
    : tolerancia(1e-8), maxIteracoes(1000) {}

// Construtor com tolerância
IntegracaoNumerica::IntegracaoNumerica(double tol, int maxIter) 
    : tolerancia(tol), maxIteracoes(maxIter) {}

// Destrutor
IntegracaoNumerica::~IntegracaoNumerica() {}

// Regra do trapézio simples
double IntegracaoNumerica::trapezioSimples(FuncaoReal f, double a, double b) const {
    return (b - a) * (f(a) + f(b)) / 2.0;
}

// Regra do trapézio composta
double IntegracaoNumerica::trapezioComposto(FuncaoReal f, double a, double b, int n) const {
    if (n < 1) {
        throw std::invalid_argument("Número de subintervalos deve ser positivo");
    }
    
    double h = (b - a) / n;
    double soma = 0.5 * (f(a) + f(b));
    
    for (int i = 1; i < n; i++) {
        soma += f(a + i * h);
    }
    
    return h * soma;
}

// Regra de Simpson 1/3 simples
double IntegracaoNumerica::simpsonSimples(FuncaoReal f, double a, double b) const {
    double meio = (a + b) / 2.0;
    return (b - a) / 6.0 * (f(a) + 4.0 * f(meio) + f(b));
}

// Regra de Simpson 1/3 composta
double IntegracaoNumerica::simpsonComposto(FuncaoReal f, double a, double b, int n) const {
    if (n < 2 || n % 2 != 0) {
        n = (n / 2) * 2;
        if (n < 2) n = 2;
    }
    
    double h = (b - a) / n;
    double soma = f(a) + f(b);
    
    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        soma += (i % 2 == 0) ? 2.0 * f(x) : 4.0 * f(x);
    }
    
    return h * soma / 3.0;
}

// Regra de Simpson 3/8
double IntegracaoNumerica::simpson38(FuncaoReal f, double a, double b, int n) const {
    if (n < 3 || n % 3 != 0) {
        n = ((n / 3) + 1) * 3;
    }
    
    double h = (b - a) / n;
    double soma = f(a) + f(b);
    
    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        soma += (i % 3 == 0) ? 2.0 * f(x) : 3.0 * f(x);
    }
    
    return 3.0 * h * soma / 8.0;
}

// Quadratura de Gauss-Legendre
double IntegracaoNumerica::gaussLegendre(FuncaoReal f, double a, double b) const {
    // Transforma para intervalo [-1, 1]
    double meio = (a + b) / 2.0;
    double semi = (b - a) / 2.0;
    
    double soma = 0.0;
    for (size_t i = 0; i < PONTOS_GAUSS.size(); i++) {
        double x = meio + semi * PONTOS_GAUSS[i];
        soma += PESOS_GAUSS[i] * f(x);
    }
    
    return semi * soma;
}

// Quadratura de Gauss-Legendre composta
double IntegracaoNumerica::gaussLegendreComposto(FuncaoReal f, double a, double b, int n) const {
    if (n < 1) n = 1;
    
    double h = (b - a) / n;
    double soma = 0.0;
    
    for (int i = 0; i < n; i++) {
        double ai = a + i * h;
        double bi = ai + h;
        soma += gaussLegendre(f, ai, bi);
    }
    
    return soma;
}

// Integração adaptativa de Romberg
double IntegracaoNumerica::romberg(FuncaoReal f, double a, double b, int ordemMax) const {
    if (ordemMax < 1) ordemMax = 1;
    
    std::vector<std::vector<double>> R(ordemMax + 1, std::vector<double>(ordemMax + 1, 0.0));
    
    // R[0][0] é regra do trapézio com 1 intervalo
    R[0][0] = trapezioSimples(f, a, b);
    
    for (int k = 1; k <= ordemMax; k++) {
        // Regra do trapézio com 2^k intervalos
        int n = 1 << k; // 2^k
        R[k][0] = trapezioComposto(f, a, b, n);
        
        // Extrapolação de Richardson
        for (int j = 1; j <= k; j++) {
            double fator = std::pow(4.0, j);
            R[k][j] = R[k][j-1] + (R[k][j-1] - R[k-1][j-1]) / (fator - 1.0);
        }
        
        // Verifica convergência
        if (k > 0 && std::abs(R[k][k] - R[k-1][k-1]) < tolerancia) {
            return R[k][k];
        }
    }
    
    return R[ordemMax][ordemMax];
}

// Integração adaptativa de Simpson
double IntegracaoNumerica::simpsonAdaptativo(FuncaoReal f, double a, double b, 
                                              double tol, int profundidade) const {
    double c = (a + b) / 2.0;
    
    double S1 = simpsonSimples(f, a, b);
    double S2 = simpsonSimples(f, a, c) + simpsonSimples(f, c, b);
    
    if (std::abs(S2 - S1) < 15.0 * tol || profundidade >= 20) {
        return S2 + (S2 - S1) / 15.0;
    }
    
    return simpsonAdaptativo(f, a, c, tol / 2.0, profundidade + 1) +
           simpsonAdaptativo(f, c, b, tol / 2.0, profundidade + 1);
}

// Derivada numérica (diferença central)
double IntegracaoNumerica::derivada(FuncaoReal f, double x, double h) const {
    return (f(x + h) - f(x - h)) / (2.0 * h);
}

// Segunda derivada numérica
double IntegracaoNumerica::segundaDerivada(FuncaoReal f, double x, double h) const {
    return (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h);
}

// Derivada com extrapolação de Richardson
double IntegracaoNumerica::derivadaRichardson(FuncaoReal f, double x, double h, int ordem) const {
    if (ordem < 1) ordem = 1;
    
    std::vector<std::vector<double>> D(ordem + 1, std::vector<double>(ordem + 1, 0.0));
    
    for (int i = 0; i <= ordem; i++) {
        double hi = h / std::pow(2.0, i);
        D[i][0] = (f(x + hi) - f(x - hi)) / (2.0 * hi);
    }
    
    for (int j = 1; j <= ordem; j++) {
        for (int i = j; i <= ordem; i++) {
            double fator = std::pow(4.0, j);
            D[i][j] = D[i][j-1] + (D[i][j-1] - D[i-1][j-1]) / (fator - 1.0);
        }
    }
    
    return D[ordem][ordem];
}

// Gradiente numérico
std::vector<double> IntegracaoNumerica::gradiente(FuncaoMultivariada f, 
                                                   const std::vector<double>& ponto, 
                                                   double h) const {
    size_t n = ponto.size();
    std::vector<double> grad(n);
    
    for (size_t i = 0; i < n; i++) {
        std::vector<double> pontoMais = ponto;
        std::vector<double> pontoMenos = ponto;
        pontoMais[i] += h;
        pontoMenos[i] -= h;
        grad[i] = (f(pontoMais) - f(pontoMenos)) / (2.0 * h);
    }
    
    return grad;
}

// Integração de dados discretos (trapézio)
double IntegracaoNumerica::integrarDadosDiscretos(const std::vector<double>& x, 
                                                   const std::vector<double>& y) const {
    if (x.size() != y.size() || x.size() < 2) {
        throw std::invalid_argument("Vetores devem ter mesmo tamanho e pelo menos 2 elementos");
    }
    
    double integral = 0.0;
    for (size_t i = 1; i < x.size(); i++) {
        integral += (x[i] - x[i-1]) * (y[i] + y[i-1]) / 2.0;
    }
    
    return integral;
}

// Integração de dados discretos (Simpson)
double IntegracaoNumerica::integrarDadosDiscretosSimpson(const std::vector<double>& x,
                                                          const std::vector<double>& y) const {
    if (x.size() != y.size() || x.size() < 3) {
        throw std::invalid_argument("Vetores devem ter mesmo tamanho e pelo menos 3 elementos");
    }
    
    double integral = 0.0;
    size_t n = x.size();
    
    // Usa Simpson para pares de intervalos
    for (size_t i = 0; i + 2 < n; i += 2) {
        double h1 = x[i+1] - x[i];
        double h2 = x[i+2] - x[i+1];
        
        // Simpson para intervalos não-uniformes
        double somah = h1 + h2;
        integral += (somah / 6.0) * (
            y[i] * (2.0 - h2/h1) +
            y[i+1] * (somah * somah) / (h1 * h2) +
            y[i+2] * (2.0 - h1/h2)
        );
    }
    
    // Se número ímpar de intervalos, usa trapézio para o último
    if ((n - 1) % 2 != 0) {
        integral += (x[n-1] - x[n-2]) * (y[n-1] + y[n-2]) / 2.0;
    }
    
    return integral;
}

// Derivada discreta
std::vector<double> IntegracaoNumerica::derivadaDiscreta(const std::vector<double>& x,
                                                          const std::vector<double>& y) const {
    if (x.size() != y.size() || x.size() < 2) {
        throw std::invalid_argument("Vetores devem ter mesmo tamanho e pelo menos 2 elementos");
    }
    
    size_t n = x.size();
    std::vector<double> derivadas(n);
    
    // Primeiro ponto: diferença progressiva
    derivadas[0] = (y[1] - y[0]) / (x[1] - x[0]);
    
    // Pontos internos: diferença central
    for (size_t i = 1; i < n - 1; i++) {
        derivadas[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1]);
    }
    
    // Último ponto: diferença regressiva
    derivadas[n-1] = (y[n-1] - y[n-2]) / (x[n-1] - x[n-2]);
    
    return derivadas;
}

// Integral cumulativa
std::vector<double> IntegracaoNumerica::integralCumulativa(const std::vector<double>& x,
                                                            const std::vector<double>& y) const {
    if (x.size() != y.size() || x.size() < 2) {
        throw std::invalid_argument("Vetores devem ter mesmo tamanho e pelo menos 2 elementos");
    }
    
    std::vector<double> cumsum(x.size(), 0.0);
    
    for (size_t i = 1; i < x.size(); i++) {
        cumsum[i] = cumsum[i-1] + (x[i] - x[i-1]) * (y[i] + y[i-1]) / 2.0;
    }
    
    return cumsum;
}

// Setters
void IntegracaoNumerica::definirTolerancia(double tol) { tolerancia = tol; }
void IntegracaoNumerica::definirMaxIteracoes(int maxIter) { maxIteracoes = maxIter; }

// Getters
double IntegracaoNumerica::obterTolerancia() const { return tolerancia; }
int IntegracaoNumerica::obterMaxIteracoes() const { return maxIteracoes; }

} // namespace analise_ouro
