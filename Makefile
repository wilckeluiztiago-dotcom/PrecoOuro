# Makefile para Análise do Preço do Ouro
# Autor: Luiz Tiago Wilcke

CXX = g++
CXXFLAGS = -std=c++23 -Wall -Wextra -O2 -I./include -fconcepts -fcoroutines
LDFLAGS = -lm

# Diretórios
SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build

# Módulos fonte
UTILIDADES_SRC = $(wildcard $(SRC_DIR)/utilidades/*.cpp)
SERIES_SRC = $(wildcard $(SRC_DIR)/series_temporais/*.cpp)
ARIMA_SRC = $(wildcard $(SRC_DIR)/modelos_arima/*.cpp)
ESTOCASTICAS_SRC = $(wildcard $(SRC_DIR)/equacoes_estocasticas/*.cpp)
SIMULACAO_SRC = $(wildcard $(SRC_DIR)/simulacao/*.cpp)
MERCADO_SRC = $(wildcard $(SRC_DIR)/mercado/*.cpp)
VISUALIZACAO_SRC = $(wildcard $(SRC_DIR)/visualizacao/*.cpp)

# Todos os fontes
ALL_SRC = $(UTILIDADES_SRC) $(SERIES_SRC) $(ARIMA_SRC) $(ESTOCASTICAS_SRC) \
          $(SIMULACAO_SRC) $(MERCADO_SRC) $(VISUALIZACAO_SRC) $(SRC_DIR)/Main.cpp

# Objetos
ALL_OBJ = $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(notdir $(ALL_SRC)))

# Executável
TARGET = analise_ouro

.PHONY: all clean dirs

all: dirs $(TARGET)

dirs:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p resultados

$(TARGET): $(ALL_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Regras para compilar cada diretório
$(BUILD_DIR)/%.o: $(SRC_DIR)/utilidades/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILD_DIR)/%.o: $(SRC_DIR)/series_temporais/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILD_DIR)/%.o: $(SRC_DIR)/modelos_arima/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILD_DIR)/%.o: $(SRC_DIR)/equacoes_estocasticas/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILD_DIR)/%.o: $(SRC_DIR)/simulacao/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILD_DIR)/%.o: $(SRC_DIR)/mercado/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILD_DIR)/%.o: $(SRC_DIR)/visualizacao/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILD_DIR)/Main.o: $(SRC_DIR)/Main.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILD_DIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf $(BUILD_DIR) $(TARGET)
	rm -f resultados/*.png resultados/*.csv resultados/*.txt

run: all
	./$(TARGET)
