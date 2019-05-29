
# This takes about 3 minutes in my laptop
data_fig_1b ::
	for alpha in $$(seq 0.1 0.5 5.1); do \
		alpha_square=$$(echo " ( $$alpha*$$alpha )" | bc -l) ; \
		echo alpha $$alpha " ";\
		./RandomExtremalPOVMs.wl -o DispTherGaussian -n $$alpha_square -s 150;\
	done;


all ::
	date
	make data_fig_1b
	date
	

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
