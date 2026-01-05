/*
 * ModeloAR.hpp - Header Modelo AR - C++23
 * Autor: Luiz Tiago Wilcke
 */

#ifndef MODELO_AR_HPP
#define MODELO_AR_HPP

#include <vector>

namespace analise_ouro {

class ModeloAR {
private:
    std::vector<double> serieOriginal;
    std::vector<double> coeficientes;
    bool serieCarregada;
    size_t tamanhoSerie;
    int ordemP;
    bool estimado;
    double media;
    double variancia;
    double varianciaErro;
    
    void calcularMomentos();
    [[nodiscard]] auto calcularAutoCorrelacoes(int maxLag) const -> std::vector<double>;
    void verificarDados() const;

public:
    ModeloAR();
    ModeloAR(const std::vector<double>& dados, int ordem = 1);
    ~ModeloAR();

    void carregarDados(const std::vector<double>& novosDados);
    void definirOrdem(int ordem);
    [[nodiscard]] auto estimar() -> bool;
    [[nodiscard]] auto prever(int horizontePrevisao) const -> std::vector<double>;
    [[nodiscard]] auto calcularResiduos() const -> std::vector<double>;
    [[nodiscard]] auto calcularAIC() const -> double;
    [[nodiscard]] auto calcularBIC() const -> double;
    [[nodiscard]] auto obterCoeficientes() const -> const std::vector<double>&;
    [[nodiscard]] auto obterOrdem() const noexcept -> int;
    [[nodiscard]] auto obterVarianciaErro() const noexcept -> double;
};

} // namespace analise_ouro

#endif // MODELO_AR_HPP
