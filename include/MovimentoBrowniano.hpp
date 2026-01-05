/*
 * MovimentoBrowniano.hpp - Header GBM - C++23
 * Autor: Luiz Tiago Wilcke
 */

#ifndef MOVIMENTO_BROWNIANO_HPP
#define MOVIMENTO_BROWNIANO_HPP

#include <vector>
#include <random>
#include <utility>

namespace analise_ouro {

class MovimentoBrowniano {
private:
    double precoInicial;
    double drift;
    double volatilidade;
    double dt;
    mutable std::mt19937_64 gerador;

    [[nodiscard]] auto gerarNormalPadrao() -> double;

public:
    MovimentoBrowniano();
    MovimentoBrowniano(double s0, double mu, double sigma, double deltaT = 1.0/252.0);
    ~MovimentoBrowniano();

    void definirParametros(double s0, double mu, double sigma);
    void definirPassoTempo(double deltaT);
    void definirSemente(unsigned long semente);

    [[nodiscard]] auto simularTrajetoria(int numPassos) const -> std::vector<double>;
    [[nodiscard]] auto simularTrajetoriaExata(int numPassos) const -> std::vector<double>;
    [[nodiscard]] auto simularMultiplasTrajetorias(int numTrajetorias, int numPassos) const -> std::vector<std::vector<double>>;

    [[nodiscard]] auto calcularMediaEstatistica(const std::vector<std::vector<double>>& trajetorias, int passo) const -> double;
    [[nodiscard]] auto calcularVarianciaEstatistica(const std::vector<std::vector<double>>& trajetorias, int passo) const -> double;
    [[nodiscard]] auto valorEsperadoTeorico(double tempo) const -> double;
    [[nodiscard]] auto varianciaTeorica(double tempo) const -> double;
    [[nodiscard]] auto calcularRetornosLogaritmicos(const std::vector<double>& trajetoria) const -> std::vector<double>;
    [[nodiscard]] auto estimarParametros(const std::vector<double>& precos) const -> std::pair<double, double>;

    [[nodiscard]] auto obterDrift() const noexcept -> double;
    [[nodiscard]] auto obterVolatilidade() const noexcept -> double;
    [[nodiscard]] auto obterPrecoInicial() const noexcept -> double;
};

} // namespace analise_ouro

#endif // MOVIMENTO_BROWNIANO_HPP
