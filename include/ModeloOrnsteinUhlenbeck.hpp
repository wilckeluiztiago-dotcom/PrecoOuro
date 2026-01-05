/*
 * ModeloOrnsteinUhlenbeck.hpp - Header OU - C++23
 * Autor: Luiz Tiago Wilcke
 */

#ifndef MODELO_ORNSTEIN_UHLENBECK_HPP
#define MODELO_ORNSTEIN_UHLENBECK_HPP

#include <vector>
#include <random>
#include <tuple>

namespace analise_ouro {

class ModeloOrnsteinUhlenbeck {
private:
    double valorInicial, velocidadeReversao, mediaLongoPrazo, volatilidade;
    double dt;
    mutable std::mt19937_64 gerador;

public:
    ModeloOrnsteinUhlenbeck();
    ModeloOrnsteinUhlenbeck(double x0, double theta, double mu, double sigma);
    ~ModeloOrnsteinUhlenbeck();

    void definirParametros(double x0, double theta, double mu, double sigma);

    [[nodiscard]] auto simularTrajetoria(int numPassos) const -> std::vector<double>;
    [[nodiscard]] auto simularTrajetoriaExata(int numPassos) const -> std::vector<double>;
    [[nodiscard]] auto simularMultiplasTrajetorias(int numTrajetorias, int numPassos) const -> std::vector<std::vector<double>>;

    [[nodiscard]] auto mediaCondicional(double x0, double t) const -> double;
    [[nodiscard]] auto varianciaCondicional(double t) const -> double;
    [[nodiscard]] auto varianciaEstacionaria() const -> double;
    [[nodiscard]] auto autoCorrelacao(double lag) const -> double;
    [[nodiscard]] auto meiaVida() const -> double;
    [[nodiscard]] auto estimarParametros(const std::vector<double>& dados) const -> std::tuple<double, double, double>;

    [[nodiscard]] auto obterVelocidadeReversao() const noexcept -> double;
    [[nodiscard]] auto obterMediaLongoPrazo() const noexcept -> double;
    [[nodiscard]] auto obterVolatilidade() const noexcept -> double;
};

} // namespace analise_ouro

#endif // MODELO_ORNSTEIN_UHLENBECK_HPP
