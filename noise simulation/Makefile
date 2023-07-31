SHELL:=/bin/bash
.PHONY: all

var_B_T1_list:=$(shell echo 0.{1..9..2}/{1..6}/{1..6})

# 默认用 var_list, B_list 和 T1_list 的第二个数，观察变化规律
var_FIX?=0.3
B_FIX?=2
T1_FIX?=2

data=$(var_B_T1_list:%=result_R1/%.csv)
all: $(var_B_T1_list:%=fig_R1/%.png)

# 导入存好的实验参数，模拟噪声
args=$(shell echo $* | sed 's,/, ,g')
result_R1/%.csv: noise_simulation.py
	mkdir -p $(@D)
	python3 noise_simulation.py $@ --var $(word 1,$(args)) -B $(word 2,$(args)) --T1 $(word 3,$(args))

# 噪声谱的可视化
fig_R1/%.png: result_R1/%.csv plot_result1D.py
	mkdir -p $(@D)
	python3 plot_result1D.py $< $@

# 选特定 T1 和 var, 观察 B 的影响，最后一个参数是输出文件
fig_R1/$(var_FIX)/$(T1_FIX)/B_dependence.png: $(data) plot_B_dependence.py
	mkdir -p $(@D)
	python3 plot_B_dependence.py $(wildcard result_R1/$(var_FIX)/*/$(T1_FIX).csv) $@

# 选特定 var, B 的变化, 观察 T1 的影响，最后一个参数是输出文件
fig_R1/$(var_FIX)/$(B_FIX)/T1_dependence.png: $(data) plot_T1_dependence.py
	mkdir -p $(@D)
	python3 plot_T1_dependence.py $(wildcard result_R1/$(var_FIX)/$(B_FIX)/*.csv) $@

# 选特定 T1, B，画出噪声谱的积分（噪声总功率）随 var 的变化 （var ~ 温度），最后一个参数是输出文件
fig_R1/$(B_FIX)/$(T1_FIX)/var_dependence.png: $(data) plot_var_dependence.py
	mkdir -p $(@D)
	python3 plot_var_dependence.py $(wildcard result_R1/*/$(B_FIX)/$(T1_FIX)) $@

# 生成作业说明
README.tex: README.md
	pandoc --read=markdown --write=latex --output $@ --standalone < $^
	sed '/unicode-math/d' -i $@
README.pdf: README.tex
	latexmk -shell-escape -xelatex $^

clean:
	rm -r param/* result_R1/* fig_R1/*

.DELETE_ON_ERROR:
.SECONDARY:
