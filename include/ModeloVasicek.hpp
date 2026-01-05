/*
 * ModeloVasicek.hpp - Header Vasicek - C++23
 * Autor: Luiz Tiago Wilcke
 */

#ifndef MODELO_VASICEK_HPP
#define MODELO_VASICEK_HPP

#include <vector>
#include <random>
#include <tuple>

namespace analise_ouro {

class ModeloVasicek {
private:
    double taxaInicial, velocidadeReversao, mediaLongoPrazo, volatilidade;
    double dt;
    mutable std::mt19937_64 gerador;

public:
    ModeloVasicek();
    ModeloVasicek(double r0, double a, double b, double sigma);
    ~ModeloVasicek();

    void definirParametros(double r0, double a, double b, double sigma);

    [[nodiscard]] auto simularTrajetoria(int numPassos) const -> std::vector<double>;
    [[nodiscard]] auto simularTrajetoriaExata(int numPassos) const -> std::vector<double>;
    [[nodiscard]] auto mediaCondicional(double r0, double t) const -> double;
    [[nodiscard]] auto varianciaCondicional(double t) const -> double;
    [[nodiscard]] auto precoBondZeroCoupon(double taxaAtual, double maturidade) const -> double;
    [[nodiscard]] auto taxaForward(double taxaAtual, double t, double T) const -> double;
    [[nodiscard]] auto estimarParametros(const std::vector<double>& taxasHistoricas) const -> std::tuple<double, double, double>;

    [[nodiscard]] auto obterVelocidadeReversao() const noexcept -> double;
    [[nodiscard]] auto obterMediaLongoPrazo() const noexcept -> double;
    [[nodiscard]] auto obterVolatilidade() const noexcept -> double;
};

} // namespace analise_ouro

#endif // MODELO_VASICEK_HPP
