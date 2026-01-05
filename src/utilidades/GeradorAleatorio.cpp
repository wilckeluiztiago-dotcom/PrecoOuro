/*
 * GeradorAleatorio.cpp
 * Módulo de geração de números aleatórios para simulações estocásticas
 * 
 * Autor: Luiz Tiago Wilcke
 * Projeto: Análise do Preço do Ouro
 */

#include "GeradorAleatorio.hpp"
#include <cmath>
#include <chrono>
#include <stdexcept>
#include <numeric>
#include <algorithm>

namespace analise_ouro {

// Construtor padrão com semente baseada no tempo
GeradorAleatorio::GeradorAleatorio() {
    auto agora = std::chrono::high_resolution_clock::now();
    auto duracao = agora.time_since_epoch();
    sementeAtual = static_cast<unsigned long>(duracao.count());
    geradorMersenne.seed(sementeAtual);
    inicializado = true;
    contagemGeracoes = 0;
}

// Construtor com semente específica
GeradorAleatorio::GeradorAleatorio(unsigned long semente) {
    sementeAtual = semente;
    geradorMersenne.seed(sementeAtual);
    inicializado = true;
    contagemGeracoes = 0;
}

// Destrutor
GeradorAleatorio::~GeradorAleatorio() {
    inicializado = false;
}

// Reinicializa o gerador com nova semente
void GeradorAleatorio::reinicializar(unsigned long novaSemente) {
    sementeAtual = novaSemente;
    geradorMersenne.seed(sementeAtual);
    contagemGeracoes = 0;
}

// Gera número uniforme entre 0 e 1
double GeradorAleatorio::uniformeZeroUm() {
    if (!inicializado) {
        throw std::runtime_error("Gerador não inicializado");
    }
    std::uniform_real_distribution<double> distribuicao(0.0, 1.0);
    contagemGeracoes++;
    return distribuicao(geradorMersenne);
}

// Gera número uniforme entre min e max
double GeradorAleatorio::uniforme(double minimo, double maximo) {
    if (minimo >= maximo) {
        throw std::invalid_argument("Mínimo deve ser menor que máximo");
    }
    std::uniform_real_distribution<double> distribuicao(minimo, maximo);
    contagemGeracoes++;
    return distribuicao(geradorMersenne);
}

// Gera número inteiro uniforme
int GeradorAleatorio::inteiroUniforme(int minimo, int maximo) {
    if (minimo > maximo) {
        throw std::invalid_argument("Mínimo deve ser menor ou igual ao máximo");
    }
    std::uniform_int_distribution<int> distribuicao(minimo, maximo);
    contagemGeracoes++;
    return distribuicao(geradorMersenne);
}

// Método Box-Muller para distribuição normal
double GeradorAleatorio::normalBoxMuller(double media, double desvioPadrao) {
    if (desvioPadrao <= 0) {
        throw std::invalid_argument("Desvio padrão deve ser positivo");
    }
    
    double u1 = uniformeZeroUm();
    double u2 = uniformeZeroUm();
    
    // Evita log(0)
    while (u1 <= 1e-10) {
        u1 = uniformeZeroUm();
    }
    
    // Transformação Box-Muller
    double z0 = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
    
    return media + desvioPadrao * z0;
}

// Gera número com distribuição normal padrão
double GeradorAleatorio::normalPadrao() {
    return normalBoxMuller(0.0, 1.0);
}

// Gera número com distribuição normal
double GeradorAleatorio::normal(double media, double desvioPadrao) {
    std::normal_distribution<double> distribuicao(media, desvioPadrao);
    contagemGeracoes++;
    return distribuicao(geradorMersenne);
}

// Gera número com distribuição exponencial
double GeradorAleatorio::exponencial(double taxa) {
    if (taxa <= 0) {
        throw std::invalid_argument("Taxa deve ser positiva");
    }
    std::exponential_distribution<double> distribuicao(taxa);
    contagemGeracoes++;
    return distribuicao(geradorMersenne);
}

// Gera número com distribuição de Poisson
int GeradorAleatorio::poisson(double lambda) {
    if (lambda <= 0) {
        throw std::invalid_argument("Lambda deve ser positivo");
    }
    std::poisson_distribution<int> distribuicao(lambda);
    contagemGeracoes++;
    return distribuicao(geradorMersenne);
}

// Gera número com distribuição gamma
double GeradorAleatorio::gamma(double forma, double escala) {
    if (forma <= 0 || escala <= 0) {
        throw std::invalid_argument("Forma e escala devem ser positivos");
    }
    std::gamma_distribution<double> distribuicao(forma, escala);
    contagemGeracoes++;
    return distribuicao(geradorMersenne);
}

// Gera número com distribuição beta usando gamma
double GeradorAleatorio::beta(double alfa, double beta_param) {
    if (alfa <= 0 || beta_param <= 0) {
        throw std::invalid_argument("Parâmetros alfa e beta devem ser positivos");
    }
    
    double x = gamma(alfa, 1.0);
    double y = gamma(beta_param, 1.0);
    
    return x / (x + y);
}

// Gera número com distribuição log-normal
double GeradorAleatorio::logNormal(double mediaLog, double desvioPadraoLog) {
    if (desvioPadraoLog <= 0) {
        throw std::invalid_argument("Desvio padrão deve ser positivo");
    }
    std::lognormal_distribution<double> distribuicao(mediaLog, desvioPadraoLog);
    contagemGeracoes++;
    return distribuicao(geradorMersenne);
}

// Gera número com distribuição de Cauchy
double GeradorAleatorio::cauchy(double localizacao, double escala) {
    if (escala <= 0) {
        throw std::invalid_argument("Escala deve ser positiva");
    }
    std::cauchy_distribution<double> distribuicao(localizacao, escala);
    contagemGeracoes++;
    return distribuicao(geradorMersenne);
}

// Gera número com distribuição t de Student
double GeradorAleatorio::tStudent(double grausLiberdade) {
    if (grausLiberdade <= 0) {
        throw std::invalid_argument("Graus de liberdade devem ser positivos");
    }
    std::student_t_distribution<double> distribuicao(grausLiberdade);
    contagemGeracoes++;
    return distribuicao(geradorMersenne);
}

// Gera vetor de números normais correlacionados usando decomposição de Cholesky
std::vector<double> GeradorAleatorio::normaisCorrelacionados(
    const std::vector<std::vector<double>>& matrizCorrelacao,
    const std::vector<double>& medias,
    const std::vector<double>& desviosPadrao) {
    
    size_t n = matrizCorrelacao.size();
    
    if (medias.size() != n || desviosPadrao.size() != n) {
        throw std::invalid_argument("Dimensões incompatíveis");
    }
    
    // Decomposição de Cholesky
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
    
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j <= i; j++) {
            double soma = 0.0;
            for (size_t k = 0; k < j; k++) {
                soma += L[i][k] * L[j][k];
            }
            if (i == j) {
                L[i][j] = std::sqrt(matrizCorrelacao[i][i] - soma);
            } else {
                L[i][j] = (matrizCorrelacao[i][j] - soma) / L[j][j];
            }
        }
    }
    
    // Gera normais independentes
    std::vector<double> z(n);
    for (size_t i = 0; i < n; i++) {
        z[i] = normalPadrao();
    }
    
    // Aplica transformação
    std::vector<double> resultado(n);
    for (size_t i = 0; i < n; i++) {
        resultado[i] = medias[i];
        for (size_t j = 0; j <= i; j++) {
            resultado[i] += desviosPadrao[i] * L[i][j] * z[j];
        }
    }
    
    return resultado;
}

// Gera vetor de números uniformes
std::vector<double> GeradorAleatorio::vetorUniforme(size_t tamanho, double minimo, double maximo) {
    std::vector<double> resultado(tamanho);
    for (size_t i = 0; i < tamanho; i++) {
        resultado[i] = uniforme(minimo, maximo);
    }
    return resultado;
}

// Gera vetor de números normais
std::vector<double> GeradorAleatorio::vetorNormal(size_t tamanho, double media, double desvioPadrao) {
    std::vector<double> resultado(tamanho);
    for (size_t i = 0; i < tamanho; i++) {
        resultado[i] = normal(media, desvioPadrao);
    }
    return resultado;
}

// Amostragem com reposição (bootstrap)
std::vector<double> GeradorAleatorio::amostragemComReposicao(
    const std::vector<double>& populacao, size_t tamanhoAmostra) {
    
    std::vector<double> amostra(tamanhoAmostra);
    for (size_t i = 0; i < tamanhoAmostra; i++) {
        int indice = inteiroUniforme(0, static_cast<int>(populacao.size()) - 1);
        amostra[i] = populacao[indice];
    }
    return amostra;
}

// Permutação aleatória
std::vector<int> GeradorAleatorio::permutacaoAleatoria(int n) {
    std::vector<int> permutacao(n);
    std::iota(permutacao.begin(), permutacao.end(), 0);
    
    // Fisher-Yates shuffle
    for (int i = n - 1; i > 0; i--) {
        int j = inteiroUniforme(0, i);
        std::swap(permutacao[i], permutacao[j]);
    }
    
    return permutacao;
}

// Retorna a semente atual
unsigned long GeradorAleatorio::obterSemente() const {
    return sementeAtual;
}

// Retorna contagem de gerações
unsigned long GeradorAleatorio::obterContagemGeracoes() const {
    return contagemGeracoes;
}

} // namespace analise_ouro
