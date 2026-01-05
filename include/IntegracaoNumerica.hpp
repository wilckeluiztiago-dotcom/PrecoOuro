/*
 * IntegracaoNumerica.hpp
 * Header do módulo de integração numérica
 * 
 * Autor: Luiz Tiago Wilcke
 */

#ifndef INTEGRACAO_NUMERICA_HPP
#define INTEGRACAO_NUMERICA_HPP

#include <vector>
#include <functional>

namespace analise_ouro {

// Tipos de funções
using FuncaoReal = std::function<double(double)>;
using FuncaoMultivariada = std::function<double(const std::vector<double>&)>;

class IntegracaoNumerica {
private:
    double tolerancia;
    int maxIteracoes;

public:
    // Construtores e destrutor
    IntegracaoNumerica();
    IntegracaoNumerica(double tol, int maxIter);
    ~IntegracaoNumerica();

    // Regra do trapézio
    double trapezioSimples(FuncaoReal f, double a, double b) const;
    double trapezioComposto(FuncaoReal f, double a, double b, int n) const;

    // Regra de Simpson
    double simpsonSimples(FuncaoReal f, double a, double b) const;
    double simpsonComposto(FuncaoReal f, double a, double b, int n) const;
    double simpson38(FuncaoReal f, double a, double b, int n) const;
    double simpsonAdaptativo(FuncaoReal f, double a, double b, 
                             double tol, int profundidade = 0) const;

    // Quadratura de Gauss
    double gaussLegendre(FuncaoReal f, double a, double b) const;
    double gaussLegendreComposto(FuncaoReal f, double a, double b, int n) const;

    // Romberg
    double romberg(FuncaoReal f, double a, double b, int ordemMax) const;

    // Derivadas
    double derivada(FuncaoReal f, double x, double h) const;
    double segundaDerivada(FuncaoReal f, double x, double h) const;
    double derivadaRichardson(FuncaoReal f, double x, double h, int ordem) const;
    std::vector<double> gradiente(FuncaoMultivariada f, 
                                  const std::vector<double>& ponto, 
                                  double h) const;

    // Dados discretos
    double integrarDadosDiscretos(const std::vector<double>& x, 
                                  const std::vector<double>& y) const;
    double integrarDadosDiscretosSimpson(const std::vector<double>& x,
                                         const std::vector<double>& y) const;
    std::vector<double> derivadaDiscreta(const std::vector<double>& x,
                                         const std::vector<double>& y) const;
    std::vector<double> integralCumulativa(const std::vector<double>& x,
                                           const std::vector<double>& y) const;

    // Configuração
    void definirTolerancia(double tol);
    void definirMaxIteracoes(int maxIter);
    double obterTolerancia() const;
    int obterMaxIteracoes() const;
};

} // namespace analise_ouro

#endif // INTEGRACAO_NUMERICA_HPP
