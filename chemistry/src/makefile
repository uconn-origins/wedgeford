code_name=disk_chemistry
network_name=rateGR

FC = gfortran
#FC = ifort

#formFLAG = -fno-automatic
#optFLAG = -ffixed-line-length-80 -finit-local-zero -freal-4-real-8
optFLAG = -ffixed-line-length-200 -finit-local-zero  -fno-align-commons -freal-4-real-8
#optFLAG = -O1 -unroll -zero -80
FLAGS =  $(optFLAG) $(formFLAG)

ODE = dvodpk.o fd05.o id05.o ma28.o ma28_aux.o funs_dvodepk.o
# SOLVE = Drive_dvode.o Fcn.o Jacobian.o calc_rates.o time_dependent.o heatdiff.o
SOLVE = Drive_dvode.o Fcn.o Jacobian.o calc_rates.o time_dependent.o
READ = read_$(network_name).o uv_field.o
FILES = $(ODE) $(SOLVE) $(READ)


all: $(code_name)

$(code_name): $(code_name).f $(FILES)
	$(FC) $(FLAGS) $@.f $(FILES) -o $@
	@echo make complete

%.o : %.f
	$(FC) $(FLAGS) -c $<

cleanh:
	-rm $(SOLVE) $(READ)

clean:
	-rm *~ *.o
