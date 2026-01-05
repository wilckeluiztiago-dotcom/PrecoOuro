/*
 * ProcessoWiener.hpp - Header Processo Wiener - C++23
 * Autor: Luiz Tiago Wilcke
 */

#ifndef PROCESSO_WIENER_HPP
#define PROCESSO_WIENER_HPP

#include <vector>
#include <random>
#include <functional>
#include <utility>

namespace analise_ouro {

using FuncaoDerivada = std::function<double(double)>;

class ProcessoWiener {
private:
    double dt;
    mutable std::mt19937_64 gerador;

public:
    ProcessoWiener();
    explicit ProcessoWiener(double deltaT);
    ~ProcessoWiener();

    void definirPassoTempo(double deltaT);
    void definirSemente(unsigned long semente);

    [[nodiscard]] auto simularTrajetoria(int numPassos) const -> std::vector<double>;
    [[nodiscard]] auto simularIncrementos(int numPassos) const -> std::vector<double>;
    [[nodiscard]] auto simularMultiplasTrajetorias(int numTrajetorias, int numPassos) const -> std::vector<std::vector<double>>;
    [[nodiscard]] auto simularWienerCorrelacionado(int numPassos, double correlacao) const -> std::pair<std::vector<double>, std::vector<double>>;

    [[nodiscard]] auto calcularVarianciaQuadratica(const std::vector<double>& W) const -> double;
    [[nodiscard]] auto integralIto(const std::vector<double>& f, const std::vector<double>& W) const -> double;
    [[nodiscard]] auto integralStratonovich(const std::vector<double>& f, const std::vector<double>& W) const -> double;
    [[nodiscard]] auto aplicarLemaIto(const std::vector<double>& S, FuncaoDerivada df, FuncaoDerivada d2f) const -> std::vector<double>;
    [[nodiscard]] auto calcularCovarianciaEmpírica(const std::vector<double>& W1, const std::vector<double>& W2, int passo) const -> double;
    [[nodiscard]] auto propriedadeMartingale(const std::vector<std::vector<double>>& trajetorias, int passo) const -> double;

    [[nodiscard]] auto obterPassoTempo() const noexcept -> double;
};

} // namespace analise_ouro

#endif // PROCESSO_WIENER_HPP
