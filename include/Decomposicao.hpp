/*
 * Decomposicao.hpp
 * Header do módulo de decomposição - C++23
 * 
 * Autor: Luiz Tiago Wilcke
 */

#ifndef DECOMPOSICAO_HPP
#define DECOMPOSICAO_HPP

#include <vector>

namespace analise_ouro {

// Estrutura para componentes da decomposição
struct ComponentesDecomposicao {
    std::vector<double> original;
    std::vector<double> tendencia;
    std::vector<double> sazonal;
    std::vector<double> residuo;
};

class Decomposicao {
private:
    std::vector<double> serieOriginal;
    bool serieCarregada;
    size_t tamanhoSerie;
    size_t periodoSazonal;
    
    void verificarDados() const;
    [[nodiscard]] auto extrairSazonalidade(const std::vector<double>& semTendencia) const -> std::vector<double>;
    [[nodiscard]] auto extrairSazonalidadeMultiplicativa(const std::vector<double>& semTendencia) const -> std::vector<double>;

public:
    // Construtores e destrutor
    Decomposicao();
    Decomposicao(const std::vector<double>& dados, size_t periodo = 12);
    ~Decomposicao();

    // Carregamento de dados
    void carregarDados(const std::vector<double>& novosDados);
    void definirPeriodo(size_t novoPeriodo);

    // Decomposições
    [[nodiscard]] auto decomposicaoAditiva() const -> ComponentesDecomposicao;
    [[nodiscard]] auto decomposicaoMultiplicativa() const -> ComponentesDecomposicao;
    [[nodiscard]] auto decomposicaoSTL(int numIteracoes = 2) const -> ComponentesDecomposicao;

    // Componentes individuais
    [[nodiscard]] auto extrairTendencia() const -> std::vector<double>;

    // Holt-Winters
    [[nodiscard]] auto holtWinters(double alfa = 0.2, double beta = 0.1, double gamma = 0.1, 
                                   bool multiplicativo = false) const -> std::vector<double>;

    // Getters
    [[nodiscard]] auto obterTamanho() const noexcept -> size_t;
    [[nodiscard]] auto obterPeriodo() const noexcept -> size_t;
};

} // namespace analise_ouro

#endif // DECOMPOSICAO_HPP
