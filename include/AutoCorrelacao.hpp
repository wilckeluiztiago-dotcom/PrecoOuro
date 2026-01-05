/*
 * AutoCorrelacao.hpp
 * Header do módulo de autocorrelação - C++23
 * 
 * Autor: Luiz Tiago Wilcke
 */

#ifndef AUTO_CORRELACAO_HPP
#define AUTO_CORRELACAO_HPP

#include <vector>
#include <concepts>

namespace analise_ouro {

// Estrutura para resultados de testes estatísticos
struct ResultadoTeste {
    double estatistica;
    double pValor;
    int grausLiberdade;
    bool rejeitaH0;
};

class AutoCorrelacao {
private:
    std::vector<double> serieOriginal;
    bool serieCarregada;
    size_t tamanhoSerie;
    double media;
    double variancia;
    
    void calcularMomento();
    void verificarDados() const;
    double calcularPValorChiQuadrado(double x, int df) const;

public:
    // Construtores e destrutor
    AutoCorrelacao();
    explicit AutoCorrelacao(const std::vector<double>& dados);
    ~AutoCorrelacao();

    // Carregamento de dados
    void carregarDados(const std::vector<double>& novosDados);

    // Funções de autocorrelação
    [[nodiscard]] auto acf(int lag) const -> double;
    [[nodiscard]] auto acfCompleta(int maxLag) const -> std::vector<double>;
    [[nodiscard]] auto pacf(int lag) const -> double;
    [[nodiscard]] auto pacfCompleta(int maxLag) const -> std::vector<double>;

    // Testes estatísticos
    [[nodiscard]] auto testeLjungBox(int maxLag) const -> ResultadoTeste;
    [[nodiscard]] auto testeBoxPierce(int maxLag) const -> ResultadoTeste;
    [[nodiscard]] auto testeADF(int maxLag = 1) const -> ResultadoTeste;

    // Intervalos de confiança
    [[nodiscard]] auto intervaloConfiancaACF(double nivelConfianca = 0.95) const -> double;

    // Getters
    [[nodiscard]] auto obterTamanho() const noexcept -> size_t;
    [[nodiscard]] auto obterMedia() const noexcept -> double;
    [[nodiscard]] auto obterVariancia() const noexcept -> double;
};

} // namespace analise_ouro

#endif // AUTO_CORRELACAO_HPP
