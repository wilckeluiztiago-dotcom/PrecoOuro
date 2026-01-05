# Análise Estocástica do Preço do Ouro

## Autor
**Luiz Tiago Wilcke**

## Descrição
Projeto abrangente de análise do preço do ouro utilizando **equações diferenciais estocásticas (SDEs)** e **análise de séries temporais**, implementado em **C++23** com recursos avançados da linguagem.

O projeto contém aproximadamente **30 módulos** organizados em categorias, cada um com cerca de 200 linhas de código, cobrindo desde estatísticas básicas até modelos sofisticados de volatilidade estocástica.

## Estrutura do Projeto

```
analiseOuro/
├── include/                    # Headers (.hpp)
├── src/
│   ├── utilidades/            # Módulos utilitários
│   │   ├── GeradorAleatorio.cpp
│   │   ├── MatrizOperacoes.cpp
│   │   ├── EstatisticasBasicas.cpp
│   │   ├── Interpolacao.cpp
│   │   └── IntegracaoNumerica.cpp
│   ├── series_temporais/      # Análise de séries temporais
│   │   ├── MediaMovel.cpp
│   │   ├── AutoCorrelacao.cpp
│   │   ├── Decomposicao.cpp
│   │   ├── Sazonalidade.cpp
│   │   └── Tendencia.cpp
│   ├── modelos_arima/         # Modelos ARIMA e GARCH
│   │   ├── ModeloAR.cpp
│   │   ├── ModeloMA.cpp
│   │   ├── ModeloARMA.cpp
│   │   ├── ModeloARIMA.cpp
│   │   └── ModeloGARCH.cpp
│   ├── equacoes_estocasticas/ # Equações Diferenciais Estocásticas
│   │   ├── MovimentoBrowniano.cpp
│   │   ├── ProcessoWiener.cpp
│   │   ├── ModeloBlackScholes.cpp
│   │   ├── ModeloHeston.cpp
│   │   ├── ModeloVasicek.cpp
│   │   ├── ModeloCIR.cpp
│   │   └── ModeloOrnsteinUhlenbeck.cpp
│   ├── simulacao/             # Métodos de simulação
│   │   ├── SimulacaoMonteCarlo.cpp
│   │   ├── EulerMaruyama.cpp
│   │   ├── Milstein.cpp
│   │   └── BootstrapEstatistico.cpp
│   ├── mercado/               # Análise de mercado
│   │   ├── AnaliseVolatilidade.cpp
│   │   ├── CorrelacaoAtivos.cpp
│   │   ├── IndicadoresTecnicos.cpp
│   │   ├── AnaliseRisco.cpp
│   │   └── PrevisaoPrecos.cpp
│   ├── visualizacao/          # Geração de gráficos e relatórios
│   │   ├── GeradorGraficos.cpp
│   │   ├── RelatorioNumerico.cpp
│   │   └── ExportadorDados.cpp
│   └── Main.cpp               # Programa principal
├── dados/                     # Dados de entrada
├── resultados/                # Resultados e gráficos
├── Makefile
└── README.md
```

## Módulos Principais

### Equações Diferenciais Estocásticas
- **Movimento Browniano Geométrico (GBM)**: dS = μS dt + σS dW
- **Modelo de Heston**: Volatilidade estocástica com esquema QE
- **Vasicek e CIR**: Modelos de taxas de juros com reversão à média
- **Ornstein-Uhlenbeck**: Processo de reversão à média para spreads

### Séries Temporais
- **Modelos ARIMA**: AR(p), MA(q), ARMA(p,q), ARIMA(p,d,q)
- **GARCH(p,q)**: Modelagem de volatilidade condicional
- **Decomposição**: Aditiva, multiplicativa e STL
- **Médias Móveis**: SMA, EMA, WMA, TEMA, KAMA, Hull

### Análise de Risco
- **VaR e CVaR**: Histórico e paramétrico
- **Índices**: Sharpe, Sortino, Calmar, Information Ratio
- **Drawdowns**: Máximo e série temporal

### Indicadores Técnicos
- RSI, MACD, Bandas de Bollinger
- Estocástico, ATR, ADX

## Compilação

### Requisitos
- Compilador C++23 (g++ 12+ ou clang++ 15+)
- Make
- Gnuplot (opcional, para gráficos)

### Comandos
```bash
# Compilar o projeto
make

# Executar
./analise_ouro

# Limpar
make clean
```

## Recursos C++23 Utilizados
- `std::ranges` e `std::views` para processamento de dados
- Concepts (`Numerico<T>`) para programação genérica
- `[[nodiscard]]` para segurança de retorno
- `auto` com dedução de tipo de retorno
- `std::numbers` para constantes matemáticas
- Structured bindings (`auto [a, b] = ...`)

## Exemplos de Uso

```cpp
// Simulação de preços com GBM
MovimentoBrowniano gbm(1800.0, 0.05, 0.15);
auto precos = gbm.simularTrajetoriaExata(252);

// Modelo GARCH para volatilidade
ModeloGARCH garch(precos, 1, 1);
garch.estimar();
auto volatilidades = garch.preverVolatilidade(10);

// Monte Carlo para opções
SimulacaoMonteCarlo mc(10000, 252);
auto call = mc.calcularPrecoOpcaoCall(1800, 1900, 0.05, 0.15, 1.0);
```

## Licença
Este projeto é de uso educacional e de pesquisa.

## Contato
**Autor**: Luiz Tiago Wilcke
