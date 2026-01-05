/*
 * Interpolacao.hpp
 * Header do módulo de interpolação
 * 
 * Autor: Luiz Tiago Wilcke
 */

#ifndef INTERPOLACAO_HPP
#define INTERPOLACAO_HPP

#include <vector>

namespace analise_ouro {

enum class MetodoInterpolacao {
    LINEAR,
    LAGRANGE,
    NEWTON,
    SPLINE_CUBICO,
    HERMITE
};

class Interpolacao {
private:
    std::vector<double> pontosX;
    std::vector<double> pontosY;
    bool dadosCarregados;
    size_t numeroPontos;
    
    // Coeficientes para spline cúbico
    std::vector<double> coeficientesA;
    std::vector<double> coeficientesB;
    std::vector<double> coeficientesC;
    std::vector<double> coeficientesD;
    
    void ordenarPontos();
    void calcularCoeficientesSpline();
    size_t encontrarIntervalo(double x) const;
    void verificarDados() const;

public:
    // Construtores e destrutor
    Interpolacao();
    Interpolacao(const std::vector<double>& x, const std::vector<double>& y);
    ~Interpolacao();

    // Carregamento de dados
    void carregarDados(const std::vector<double>& x, const std::vector<double>& y);

    // Métodos de interpolação
    double linear(double x) const;
    double lagrange(double x) const;
    double newton(double x) const;
    double splineCubico(double x) const;
    double hermite(double x, const std::vector<double>& derivadas) const;

    // Derivadas
    double derivadaSpline(double x) const;
    double segundaDerivadaSpline(double x) const;

    // Extrapolação
    double extrapolacaoLinear(double x) const;

    // Interpolação vetorizada
    std::vector<double> interpolarVetor(const std::vector<double>& novosX,
                                        MetodoInterpolacao metodo) const;

    // Getters
    size_t obterNumeroPontos() const;
    const std::vector<double>& obterPontosX() const;
    const std::vector<double>& obterPontosY() const;
};

} // namespace analise_ouro

#endif // INTERPOLACAO_HPP
