/*
 * GeradorAleatorio.hpp
 * Header do módulo de geração de números aleatórios
 * 
 * Autor: Luiz Tiago Wilcke
 */

#ifndef GERADOR_ALEATORIO_HPP
#define GERADOR_ALEATORIO_HPP

#include <random>
#include <vector>

namespace analise_ouro {

class GeradorAleatorio {
private:
    std::mt19937_64 geradorMersenne;
    unsigned long sementeAtual;
    bool inicializado;
    unsigned long contagemGeracoes;

public:
    // Construtores e destrutor
    GeradorAleatorio();
    explicit GeradorAleatorio(unsigned long semente);
    ~GeradorAleatorio();

    // Reinicialização
    void reinicializar(unsigned long novaSemente);

    // Distribuições uniformes
    double uniformeZeroUm();
    double uniforme(double minimo, double maximo);
    int inteiroUniforme(int minimo, int maximo);

    // Distribuição normal
    double normalBoxMuller(double media, double desvioPadrao);
    double normalPadrao();
    double normal(double media, double desvioPadrao);

    // Outras distribuições
    double exponencial(double taxa);
    int poisson(double lambda);
    double gamma(double forma, double escala);
    double beta(double alfa, double beta_param);
    double logNormal(double mediaLog, double desvioPadraoLog);
    double cauchy(double localizacao, double escala);
    double tStudent(double grausLiberdade);

    // Vetores e correlações
    std::vector<double> normaisCorrelacionados(
        const std::vector<std::vector<double>>& matrizCorrelacao,
        const std::vector<double>& medias,
        const std::vector<double>& desviosPadrao);
    std::vector<double> vetorUniforme(size_t tamanho, double minimo, double maximo);
    std::vector<double> vetorNormal(size_t tamanho, double media, double desvioPadrao);

    // Amostragem
    std::vector<double> amostragemComReposicao(
        const std::vector<double>& populacao, size_t tamanhoAmostra);
    std::vector<int> permutacaoAleatoria(int n);

    // Getters
    unsigned long obterSemente() const;
    unsigned long obterContagemGeracoes() const;
};

} // namespace analise_ouro

#endif // GERADOR_ALEATORIO_HPP
