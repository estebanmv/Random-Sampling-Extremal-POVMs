
# This takes about 3 minutes in my laptop

#data_fig_1b ::
#	for alpha in $$(seq 0 0.24166 3.141592); do \
#		echo alpha " ";\
#		./RandomExtremalPOVMs.wl -o Qubit -h $$alpha -s 10;\
#	done;
#	./RandomExtremalPOVMs.wl -o Qubit -h 1.578 -s 10;\

#data_fig_2a_pre ::
#	for alpha in $$(seq 0.1 0.17272 2); do \
#		alpha_square=$$(echo " ( $$alpha*$$alpha )" | bc -l) ; \
#  		echo alpha $$alpha " ";\
#		./RandomExtremalPOVMs.wl -o DispTherGaussian -s 200 -hD 7 -od 10 -n $$alpha_square ; \
#	done;

data_fig_2a_pre ::
	for alpha in $$(seq 1.09956 0.17272 2); do \
		alpha_square=$$(echo " ( $$alpha*$$alpha )" | bc -l) ; \
  		echo alpha $$alpha " ";\
		./RandomExtremalPOVMs.wl -o CohPlusTherGaussian -s 150 -hD 9 -od 10 -mix 1 -n $$alpha_square ; \
	done;
#	for alpha in $$(seq 0.1 0.14444 1.4); do \

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
