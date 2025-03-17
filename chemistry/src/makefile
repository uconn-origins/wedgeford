code_name=disk_chemistry
code_name_flip=disk_chemistry_flip
network_name=rateGR

FC = gfortran
#FC = ifort

formFLAG = -fno-automatic -std=legacy
optFLAG = -O3 -ffixed-line-length-132 -finit-local-zero

#optFLAG = -O1 -unroll -zero -80
#optFLAG = -O3 -finit-local-zero
FLAGS =  $(optFLAG) $(formFLAG)

ODE = dvodpk.o fd05.o id05.o ma28.o ma28_aux.o funs_dvodepk.o
SOLVE_FLIP = Drive_dvode_flip.o Fcn.o Jacobian.o calc_rates_flip.o time_dependent.o heatdiff.o
SOLVE = Drive_dvode.o Fcn.o Jacobian.o calc_rates.o time_dependent.o heatdiff.o
READ = read_$(network_name).o uv_field.o
FILES = $(ODE) $(SOLVE) $(READ)

FILES_FLIP = $(ODE) $(SOLVE_FLIP) $(READ)

all: $(code_name)
	        
$(code_name): $(code_name).f $(FILES)
	$(FC) $(FLAGS) $@.f $(FILES) -o $@
	@echo make complete

flip: $(code_name_flip)

$(code_name_flip): $(code_name_flip).f $(FILES_FLIP)
	$(FC) $(FLAGS) $@.f $(FILES_FLIP) -o $@
	@echo make complete

%.o : %.f
	$(FC) $(FLAGS) -c $<

cleanf:
	-rm $(FILES_FLIP)
    
cleanh:
	-rm $(SOLVE) $(READ) 
	
