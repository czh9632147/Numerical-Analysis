# 定义编译器和选项
CXX = g++
CXXFLAGS = -Wall -g
PDFLATEX = xelatex

# 定义源文件和可执行文件
PROBLEMS = ProblemB ProblemC ProblemD ProblemE ProblemF
REPORT = report.tex

# 目标文件
EXECUTABLES = $(PROBLEMS:%=%.out)

# 默认目标
all: run report

# 编译并运行所有问题
run: $(EXECUTABLES)
	@for exe in $^; do \
		echo "Running $$exe:"; \
		./$$exe; \
	done

# 编译每个问题文件
%.out: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

# 编译报告
report: $(REPORT)
	$(PDFLATEX) $(REPORT)
	$(PDFLATEX) $(REPORT) # 两次调用以确保引用更新

# 清理生成的文件
clean:
	rm -f $(EXECUTABLES) *.aux *.log *.pdf

.PHONY: all run report clean

