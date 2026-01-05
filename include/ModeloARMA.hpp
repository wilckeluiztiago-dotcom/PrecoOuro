/*
 * ModeloARMA.hpp - Header ARMA - C++23
 * Autor: Luiz Tiago Wilcke
 */

#ifndef MODELO_ARMA_HPP
#define MODELO_ARMA_HPP

#include <vector>

namespace analise_ouro {

class ModeloARMA {
private:
    std::vector<double> serieOriginal;
    std::vector<double> coeficientesPhi;
    std::vector<double> coeficientesTheta;
    bool serieCarregada;
    size_t tamanhoSerie;
    int ordemP, ordemQ;
    bool estimado;
    double media, varianciaErro;
    
    [[nodiscard]] auto calcularResiduosInterno() const -> std::vector<double>;
    double somaQuadradosResiduos() const;
    void verificarDados() const;

public:
    ModeloARMA();
    ModeloARMA(const std::vector<double>& dados, int p = 1, int q = 1);
    ~ModeloARMA();

    void carregarDados(const std::vector<double>& novosDados);
    void definirOrdens(int p, int q);
    [[nodiscard]] auto estimar() -> bool;
    [[nodiscard]] auto prever(int horizontePrevisao) const -> std::vector<double>;
    [[nodiscard]] auto calcularResiduos() const -> std::vector<double>;
    [[nodiscard]] auto calcularAIC() const -> double;
    [[nodiscard]] auto calcularBIC() const -> double;
    [[nodiscard]] auto obterCoeficientesPhi() const -> const std::vector<double>&;
    [[nodiscard]] auto obterCoeficientesTheta() const -> const std::vector<double>&;
    [[nodiscard]] auto obterOrdemP() const noexcept -> int;
    [[nodiscard]] auto obterOrdemQ() const noexcept -> int;
};

} // namespace analise_ouro

#endif // MODELO_ARMA_HPP
