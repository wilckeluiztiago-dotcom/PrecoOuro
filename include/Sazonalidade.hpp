/*
 * Sazonalidade.hpp
 * Header do módulo de sazonalidade - C++23
 * Autor: Luiz Tiago Wilcke
 */

#ifndef SAZONALIDADE_HPP
#define SAZONALIDADE_HPP

#include <vector>
#include <utility>

namespace analise_ouro {

struct ResultadoTesteSazonalidade {
    size_t periodo;
    double estatistica;
    double pValor;
    bool temSazonalidade;
};

class Sazonalidade {
private:
    std::vector<double> serieOriginal;
    bool serieCarregada;
    size_t tamanhoSerie;
    
    void verificarDados() const;
    double calcularPValorChiQuadrado(double x, int df) const;

public:
    Sazonalidade();
    explicit Sazonalidade(const std::vector<double>& dados);
    ~Sazonalidade();

    void carregarDados(const std::vector<double>& novosDados);
    [[nodiscard]] auto detectarPeriodo(int maxPeriodo = 365) const -> size_t;
    [[nodiscard]] auto testeSazonalidade(size_t periodo) const -> ResultadoTesteSazonalidade;
    [[nodiscard]] auto calcularIndicesSazonais(size_t periodo) const -> std::vector<double>;
    [[nodiscard]] auto dessazonalizar(size_t periodo, bool multiplicativa = true) const -> std::vector<double>;
    [[nodiscard]] auto analiseEspectral(int maxFrequencias = 50) const -> std::vector<std::pair<double, double>>;
    [[nodiscard]] auto periodograma() const -> std::vector<double>;
    [[nodiscard]] auto regressaoHarmonica(size_t periodo, int numHarmonicos = 3) const -> std::vector<double>;
    [[nodiscard]] auto obterTamanho() const noexcept -> size_t;
};

} // namespace analise_ouro

#endif // SAZONALIDADE_HPP
