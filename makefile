
# This takes about 3 minutes in my laptop
data_fig_1b ::
	for eta in $$(seq 0 0.39 3.12); do \
		echo -n $$eta " ";\
		./RandomExtremalPOVMs.wl -o Qubit --EtaAngle $$eta ;\
	done;

data_fig_2a_pre ::
	for alpha in $$(seq 0.1 0.2 1.2); do \
		alpha_square=$$(echo " ( $$alpha*$$alpha )" | bc -l) ; \
  		echo alpha $$alpha " ";\
		./RandomExtremalPOVMs.wl -o CohPlusTherGaussian -s 5 -hD 5 -od 7 -mix 1 -n $$alpha_square ; \
	done;

data_fig_2b_pre ::
	for alpha in $$(seq 0.1 0.3 2.2); do \
		alpha_square=$$(echo " ( $$alpha*$$alpha )" | bc -l) ; \
  		echo alpha $$alpha " ";\
		./RandomExtremalPOVMs.wl -o CohPlusTherGamma -s 5 -hD 5 -od 7 -mix 0.5 -n $$alpha_square ; \
	done;

data_fig_2a ::
	for alpha in $$(seq 0.1 0.2 1.2); do \
		alpha_square=$$(echo " ( $$alpha*$$alpha )" | bc -l) ; \
  		echo alpha $$alpha " ";\
		./RandomExtremalPOVMs.wl -o CohPlusTherGaussian -s 150 -mix 1 -n $$alpha_square ; \
	done;

data_fig_2b ::
	for alpha in $$(seq 0.1 0.3 2.2); do \
		alpha_square=$$(echo " ( $$alpha*$$alpha )" | bc -l) ; \
  		echo alpha $$alpha " ";\
		./RandomExtremalPOVMs.wl -o CohPlusTherGamma -s 150 -mix 0.5 -n $$alpha_square ; \
	done;

all ::
	date
	make data_fig_1b
	date
	make data_fig_2a
	date
	make data_fig_2b
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
