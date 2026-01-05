/*
 * Interpolacao.cpp
 * Módulo de métodos de interpolação para análise de séries
 * 
 * Autor: Luiz Tiago Wilcke
 * Projeto: Análise do Preço do Ouro
 */

#include "Interpolacao.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <numeric>

namespace analise_ouro {

// Construtor padrão
Interpolacao::Interpolacao() : dadosCarregados(false), numeroPontos(0) {}

// Construtor com dados
Interpolacao::Interpolacao(const std::vector<double>& x, const std::vector<double>& y) {
    carregarDados(x, y);
}

// Destrutor
Interpolacao::~Interpolacao() {}

// Carrega dados para interpolação
void Interpolacao::carregarDados(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("Vetores x e y devem ter o mesmo tamanho");
    }
    
    if (x.size() < 2) {
        throw std::invalid_argument("Necessário pelo menos 2 pontos para interpolação");
    }
    
    pontosX = x;
    pontosY = y;
    numeroPontos = x.size();
    dadosCarregados = true;
    
    // Ordena por x se necessário
    ordenarPontos();
    
    // Calcula coeficientes para spline cúbico
    calcularCoeficientesSpline();
}

// Ordena pontos por coordenada x
void Interpolacao::ordenarPontos() {
    std::vector<std::pair<double, double>> pareado(numeroPontos);
    for (size_t i = 0; i < numeroPontos; i++) {
        pareado[i] = {pontosX[i], pontosY[i]};
    }
    
    std::sort(pareado.begin(), pareado.end());
    
    for (size_t i = 0; i < numeroPontos; i++) {
        pontosX[i] = pareado[i].first;
        pontosY[i] = pareado[i].second;
    }
}

// Interpolação linear básica
double Interpolacao::linear(double x) const {
    verificarDados();
    
    // Encontra intervalo
    size_t indice = encontrarIntervalo(x);
    
    if (indice >= numeroPontos - 1) {
        indice = numeroPontos - 2;
    }
    
    double x0 = pontosX[indice];
    double x1 = pontosX[indice + 1];
    double y0 = pontosY[indice];
    double y1 = pontosY[indice + 1];
    
    if (std::abs(x1 - x0) < 1e-10) {
        return y0;
    }
    
    double t = (x - x0) / (x1 - x0);
    return y0 + t * (y1 - y0);
}

// Interpolação polinomial de Lagrange
double Interpolacao::lagrange(double x) const {
    verificarDados();
    
    double resultado = 0.0;
    
    for (size_t i = 0; i < numeroPontos; i++) {
        double termo = pontosY[i];
        
        for (size_t j = 0; j < numeroPontos; j++) {
            if (i != j) {
                double denominador = pontosX[i] - pontosX[j];
                if (std::abs(denominador) < 1e-10) {
                    throw std::runtime_error("Pontos x muito próximos para Lagrange");
                }
                termo *= (x - pontosX[j]) / denominador;
            }
        }
        
        resultado += termo;
    }
    
    return resultado;
}

// Interpolação de Newton (diferenças divididas)
double Interpolacao::newton(double x) const {
    verificarDados();
    
    // Calcula tabela de diferenças divididas
    std::vector<std::vector<double>> dd(numeroPontos, std::vector<double>(numeroPontos, 0.0));
    
    // Primeira coluna é y
    for (size_t i = 0; i < numeroPontos; i++) {
        dd[i][0] = pontosY[i];
    }
    
    // Calcula diferenças divididas
    for (size_t j = 1; j < numeroPontos; j++) {
        for (size_t i = 0; i < numeroPontos - j; i++) {
            double denominador = pontosX[i + j] - pontosX[i];
            if (std::abs(denominador) < 1e-10) {
                throw std::runtime_error("Pontos x muito próximos para Newton");
            }
            dd[i][j] = (dd[i + 1][j - 1] - dd[i][j - 1]) / denominador;
        }
    }
    
    // Avalia polinômio de Newton
    double resultado = dd[0][0];
    double produto = 1.0;
    
    for (size_t i = 1; i < numeroPontos; i++) {
        produto *= (x - pontosX[i - 1]);
        resultado += dd[0][i] * produto;
    }
    
    return resultado;
}

// Calcula coeficientes para spline cúbico natural
void Interpolacao::calcularCoeficientesSpline() {
    if (numeroPontos < 3) {
        coeficientesA = pontosY;
        coeficientesB.resize(numeroPontos, 0.0);
        coeficientesC.resize(numeroPontos, 0.0);
        coeficientesD.resize(numeroPontos, 0.0);
        return;
    }
    
    size_t n = numeroPontos - 1;
    
    // Vetores h (intervalos)
    std::vector<double> h(n);
    for (size_t i = 0; i < n; i++) {
        h[i] = pontosX[i + 1] - pontosX[i];
    }
    
    // Coeficientes a
    coeficientesA = pontosY;
    
    // Sistema tridiagonal para c
    std::vector<double> alfa(n, 0.0);
    for (size_t i = 1; i < n; i++) {
        alfa[i] = (3.0 / h[i]) * (coeficientesA[i + 1] - coeficientesA[i]) -
                  (3.0 / h[i - 1]) * (coeficientesA[i] - coeficientesA[i - 1]);
    }
    
    // Resolve sistema tridiagonal
    std::vector<double> l(numeroPontos, 1.0);
    std::vector<double> mu(numeroPontos, 0.0);
    std::vector<double> z(numeroPontos, 0.0);
    
    for (size_t i = 1; i < n; i++) {
        l[i] = 2.0 * (pontosX[i + 1] - pontosX[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alfa[i] - h[i - 1] * z[i - 1]) / l[i];
    }
    
    // Back substitution
    coeficientesC.resize(numeroPontos, 0.0);
    coeficientesB.resize(n, 0.0);
    coeficientesD.resize(n, 0.0);
    
    for (int j = static_cast<int>(n) - 1; j >= 0; j--) {
        size_t i = static_cast<size_t>(j);
        coeficientesC[i] = z[i] - mu[i] * coeficientesC[i + 1];
        coeficientesB[i] = (coeficientesA[i + 1] - coeficientesA[i]) / h[i] -
                           h[i] * (coeficientesC[i + 1] + 2.0 * coeficientesC[i]) / 3.0;
        coeficientesD[i] = (coeficientesC[i + 1] - coeficientesC[i]) / (3.0 * h[i]);
    }
}

// Interpolação por spline cúbico
double Interpolacao::splineCubico(double x) const {
    verificarDados();
    
    size_t indice = encontrarIntervalo(x);
    if (indice >= numeroPontos - 1) {
        indice = numeroPontos - 2;
    }
    
    double dx = x - pontosX[indice];
    
    return coeficientesA[indice] + 
           coeficientesB[indice] * dx + 
           coeficientesC[indice] * dx * dx + 
           coeficientesD[indice] * dx * dx * dx;
}

// Derivada do spline cúbico
double Interpolacao::derivadaSpline(double x) const {
    verificarDados();
    
    size_t indice = encontrarIntervalo(x);
    if (indice >= numeroPontos - 1) {
        indice = numeroPontos - 2;
    }
    
    double dx = x - pontosX[indice];
    
    return coeficientesB[indice] + 
           2.0 * coeficientesC[indice] * dx + 
           3.0 * coeficientesD[indice] * dx * dx;
}

// Segunda derivada do spline
double Interpolacao::segundaDerivadaSpline(double x) const {
    verificarDados();
    
    size_t indice = encontrarIntervalo(x);
    if (indice >= numeroPontos - 1) {
        indice = numeroPontos - 2;
    }
    
    double dx = x - pontosX[indice];
    
    return 2.0 * coeficientesC[indice] + 6.0 * coeficientesD[indice] * dx;
}

// Interpolação de Hermite
double Interpolacao::hermite(double x, const std::vector<double>& derivadas) const {
    verificarDados();
    
    if (derivadas.size() != numeroPontos) {
        throw std::invalid_argument("Número de derivadas deve ser igual ao número de pontos");
    }
    
    size_t indice = encontrarIntervalo(x);
    if (indice >= numeroPontos - 1) {
        indice = numeroPontos - 2;
    }
    
    double x0 = pontosX[indice];
    double x1 = pontosX[indice + 1];
    double y0 = pontosY[indice];
    double y1 = pontosY[indice + 1];
    double m0 = derivadas[indice];
    double m1 = derivadas[indice + 1];
    
    double h = x1 - x0;
    double t = (x - x0) / h;
    double t2 = t * t;
    double t3 = t2 * t;
    
    // Funções base de Hermite
    double h00 = 2*t3 - 3*t2 + 1;
    double h10 = t3 - 2*t2 + t;
    double h01 = -2*t3 + 3*t2;
    double h11 = t3 - t2;
    
    return h00 * y0 + h10 * h * m0 + h01 * y1 + h11 * h * m1;
}

// Extrapolação linear
double Interpolacao::extrapolacaoLinear(double x) const {
    verificarDados();
    
    if (x < pontosX.front()) {
        // Extrapola usando os dois primeiros pontos
        double dx = pontosX[1] - pontosX[0];
        double dy = pontosY[1] - pontosY[0];
        double inclinacao = (std::abs(dx) > 1e-10) ? dy / dx : 0.0;
        return pontosY[0] + inclinacao * (x - pontosX[0]);
    } else {
        // Extrapola usando os dois últimos pontos
        size_t n = numeroPontos;
        double dx = pontosX[n-1] - pontosX[n-2];
        double dy = pontosY[n-1] - pontosY[n-2];
        double inclinacao = (std::abs(dx) > 1e-10) ? dy / dx : 0.0;
        return pontosY[n-1] + inclinacao * (x - pontosX[n-1]);
    }
}

// Interpola vetor de pontos
std::vector<double> Interpolacao::interpolarVetor(const std::vector<double>& novosX, 
                                                   MetodoInterpolacao metodo) const {
    verificarDados();
    
    std::vector<double> novosY(novosX.size());
    
    for (size_t i = 0; i < novosX.size(); i++) {
        switch (metodo) {
            case MetodoInterpolacao::LINEAR:
                novosY[i] = linear(novosX[i]);
                break;
            case MetodoInterpolacao::LAGRANGE:
                novosY[i] = lagrange(novosX[i]);
                break;
            case MetodoInterpolacao::NEWTON:
                novosY[i] = newton(novosX[i]);
                break;
            case MetodoInterpolacao::SPLINE_CUBICO:
                novosY[i] = splineCubico(novosX[i]);
                break;
            default:
                novosY[i] = linear(novosX[i]);
        }
    }
    
    return novosY;
}

// Encontra intervalo para x dado
size_t Interpolacao::encontrarIntervalo(double x) const {
    // Busca binária
    if (x <= pontosX.front()) return 0;
    if (x >= pontosX.back()) return numeroPontos - 1;
    
    size_t esquerda = 0;
    size_t direita = numeroPontos - 1;
    
    while (direita - esquerda > 1) {
        size_t meio = (esquerda + direita) / 2;
        if (pontosX[meio] > x) {
            direita = meio;
        } else {
            esquerda = meio;
        }
    }
    
    return esquerda;
}

// Verifica se dados foram carregados
void Interpolacao::verificarDados() const {
    if (!dadosCarregados || numeroPontos == 0) {
        throw std::runtime_error("Dados não carregados para interpolação");
    }
}

// Getters
size_t Interpolacao::obterNumeroPontos() const { return numeroPontos; }
const std::vector<double>& Interpolacao::obterPontosX() const { return pontosX; }
const std::vector<double>& Interpolacao::obterPontosY() const { return pontosY; }

} // namespace analise_ouro
