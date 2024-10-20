NAME := economics

ifeq ($(NAME),)
$(error Переменная NAME не задана. Установите имя проекта в NAME.)
endif

default: all

# Цель "all" собирает и открывает отчет
all: build post run

# Цель "build" собирает LaTeX-отчет
build:
	latexmk -xelatex -synctex=1 main.tex -jobname=$(NAME)

# Цель "post" удаляет некоторые временные файлы LaTeX
post:
	find . -type f \( \
		-name "*.aux" -o \
		-name "*.fdb_latexmk" -o \
		-name "*.fls" -o \
		-name "*.lof" -o \
		-name "*.lot" -o \
		-name "*.log" \
	\) -print -delete

# Цель "run" открывает сгенерированный PDF-файл
run:
	zathura $(NAME).pdf &

# Цель "clean" удаляет все временные файлы LaTeX
clean:
	find . -type f \( \
		-name "*.aux" -o \
		-name "*.fdb_latexmk" -o \
		-name "*.fls" -o \
		-name "*.lof" -o \
		-name "*.lot" -o \
		-name "*.log" -o \
		-name "*.out" -o \
		-name "*.synctex.gz" -o \
		-name "*.xdv" -o \
		-name "*.toc" \
	\) -print -delete
