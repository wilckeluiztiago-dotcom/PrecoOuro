/*
 * MediaMovel.hpp
 * Header do módulo de médias móveis - C++23
 * 
 * Autor: Luiz Tiago Wilcke
 */

#ifndef MEDIA_MOVEL_HPP
#define MEDIA_MOVEL_HPP

#include <vector>
#include <concepts>
#include <span>

namespace analise_ouro {

// Concept para tipos numéricos
template<typename T>
concept Numerico = std::is_arithmetic_v<T>;

class MediaMovel {
private:
    std::vector<double> serieOriginal;
    bool serieCarregada;
    size_t tamanhoSerie;
    
    void verificarDados() const;

public:
    // Construtores e destrutor
    MediaMovel();
    explicit MediaMovel(const std::vector<double>& dados);
    ~MediaMovel();

    // Carregamento de dados
    void carregarDados(const std::vector<double>& novosDados);

    // Médias móveis simples
    [[nodiscard]] auto sma(size_t janela) const -> std::vector<double>;
    [[nodiscard]] auto smaComPadding(size_t janela) const -> std::vector<double>;

    // Médias móveis exponenciais
    [[nodiscard]] auto ema(size_t janela) const -> std::vector<double>;
    [[nodiscard]] auto emaCustom(double alfa) const -> std::vector<double>;

    // Médias móveis ponderadas
    [[nodiscard]] auto wma(size_t janela) const -> std::vector<double>;

    // Médias móveis avançadas
    [[nodiscard]] auto tma(size_t janela) const -> std::vector<double>;
    [[nodiscard]] auto hma(size_t janela) const -> std::vector<double>;
    [[nodiscard]] auto kama(size_t janela, size_t janelaRapida = 2, size_t janelaLenta = 30) const -> std::vector<double>;
    [[nodiscard]] auto wilderSmoothing(size_t janela) const -> std::vector<double>;
    [[nodiscard]] auto dema(size_t janela) const -> std::vector<double>;
    [[nodiscard]] auto tema(size_t janela) const -> std::vector<double>;

    // Getters
    [[nodiscard]] auto obterTamanho() const noexcept -> size_t;
    [[nodiscard]] auto obterDados() const noexcept -> const std::vector<double>&;
};

} // namespace analise_ouro

#endif // MEDIA_MOVEL_HPP
