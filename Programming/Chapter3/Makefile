# 定义编译器和选项
CXX = g++
CXXFLAGS = -Ilibs/eigen -Isrc

# 定义源文件和可执行文件
SRC_DIR = src
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
EXECUTABLES = $(SOURCES:$(SRC_DIR)/%.cpp=bin/%)

# 创建 bin 目录
BIN_DIR = bin
DOC_DIR = doc
TEX_FILES = $(wildcard $(DOC_DIR)/*.tex)
PDF_FILES = $(addprefix $(DOC_DIR)/, $(notdir $(TEX_FILES:.tex=.pdf)))

# 默认目标
all: $(BIN_DIR) $(EXECUTABLES) $(PDF_FILES)

# 生成 bin 目录
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# 编译每个源文件
$(BIN_DIR)/%: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

# 编译每个 .tex 文件生成 PDF
$(DOC_DIR)/%.pdf: $(DOC_DIR)/%.tex
	cd $(DOC_DIR) && pdflatex -interaction=nonstopmode $(notdir $<)
	cd $(DOC_DIR) && pdflatex -interaction=nonstopmode $(notdir $<)
	cd $(DOC_DIR) && pdflatex -interaction=nonstopmode $(notdir $<)

# 清理生成的文件
clean:
	rm -rf $(BIN_DIR)/* $(DOC_DIR)/*.log $(DOC_DIR)/*.out $(DOC_DIR)/*.toc $(DOC_DIR)/*.aux $(DOC_DIR)/*.snm $(DOC_DIR)/*.nav $(DOC_DIR)/*.synctex.gz

.PHONY: all clean
