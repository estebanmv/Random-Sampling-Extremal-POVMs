
# This takes about 3 minutes in my laptop

data_fig_1b ::
#	for alpha in $$(seq 0 0.24166 3.141592); do \
#		echo alpha " "; for beta in $$(seq 1 1 10); do\
#		./RandomExtremalPOVMs.wl -o Qubit -h $$alpha -s 1000;\
#		done ; \
#	done; \

#	for beta in $$(seq 1 1 10); do \
#	./RandomExtremalPOVMs.wl -o Qubit -h 1.578 -s 1000;\
#done;

# Tiempos-Errores Qubit Algoritmo promedios
	for alpha in 10 30 50 70 100 150 200 500 700; do \
	echo alpha " "; for beta in $$(seq 1 1 10); do\
	./RandomExtremalPOVMs.wl -o Qubit -h 1.5707963267948966 -s $$alpha;\
	done ; \
	done; \

# Tiempos-Errores Qubit Naimark promedios
	for alpha in 10 30 50 100 150 200 300 500 700 1000 1300 1500 2000 2500 3000 4000 5000 ; do \
	echo alpha " "; for beta in $$(seq 1 1 10); do\
	./RandomExtremalPOVMs.wl -o QubitNaimark -h 1.5707963267948966 -s $$alpha;\
	done ; \
	done; \

# Tiempos-Errores Qubit Carlos-Chuang promedios
	for alpha in 10 30 50 100 150 200 300 500 700 1000 1300 1500 2000 2500 3000 4000 5000 ; do \
	echo alpha " "; for beta in $$(seq 1 1 10); do\
	./RandomExtremalPOVMs.wl -o QubitCC -h 1.5707963267948966 -s $$alpha;\
	done ; \
	done; \

#data_fig_2a_pre ::
#	for alpha in $$(seq 0.1 0.17272 2); do \
#		alpha_square=$$(echo " ( $$alpha*$$alpha )" | bc -l) ; \
#  		echo alpha $$alpha " ";\
#		./RandomExtremalPOVMs.wl -o DispTherGaussian -s 200 -hD 7 -od 10 -n $$alpha_square ; \
#	done;

#data_fig_2a_pre ::
#	for alpha in $$(seq 0.1 0.14444 1.4); do \
#	for alpha in $$(seq 0.1 0.17272 2); do \
#		alpha_square=$$(echo " ( $$alpha*$$alpha )" | bc -l) ; \
#  		echo alpha $$alpha " ";\
#		./RandomExtremalPOVMs.wl -o CohPlusTherGaussian -s 150 -hD 7 -od 10 -mix 1 -n $$alpha_square ; \
#	done;
#	./RandomExtremalPOVMs.wl -o QubitNaimark -h 1.578 -s 1000;\

#data_fig_2a_pre ::
#	for alpha in $$(seq 0.1 0.17272 2); do \
#		alpha_square=$$(echo " ( $$alpha*$$alpha )" | bc -l) ; \
#  		echo alpha $$alpha " ";\
#		./RandomExtremalPOVMs.wl -o DispTherGaussian -s 200 -hD 7 -od 10 -n $$alpha_square ; \
#	done;

#		alpha_square=$$(echo " ( $$alpha*$$alpha )" | bc -l) ; \
#$$alpha
#all ::
#	date
#	make data_fig_1b
#	date
	

# 		./RandomExtremalPOVMs.wl -o CohPlusTherGaussian -s 5 -hD 5 -od 7 -mix 1 -n $$alpha ; \

#  		echo $$alpha " ";\


# 	echo $$0
# 	NUMBERS := $(seq 0 0.39 3.12)

# 	num1=1 ; while [[ $$num1 -le 4 ]] ; do \
# 	    num2=1 ; while [[ $$num2 -le 3 ]] ; do \
# 	        echo $$num1 $$num2 ; \
# 	        ((num2 = num2 + 1)) ; \
# 	    done ; \
# 	    ((num1 = num1 + 1)) ; \
# 	done
# De 0 a pi con 12 puntos

# ./RandomExtremalPOVMs.wl -o Qubit --EtaAngle 1.
