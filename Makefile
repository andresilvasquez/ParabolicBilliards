# Rutas de include
EIGEN_INC = /usr/include/eigen3 #ruta biblioteca Eigen
MYINC     = include				#ruta includes personalizados

# Compilador y flags
CXX       = g++										#compilador
CXXFLAGS  = -std=c++17 -Ofast -fopenmp -I$(EIGEN_INC) -I$(MYINC)	#flags c17 + rutas de includes

# Graficar solo ls frontera
BOUNDARY_SRC = plot_boundary/dump_boundary.cpp		#codigo fuente
BOUNDARY_BIN = dump_boundary						#binario generado
BOUNDARY_DAT = boundary.dat							#archivo de datos de salida
BOUNDARY_PY  = plot_boundary/plt_boundary.py		#script de visualizacion

# Graficas de densidad de probabilidad y distribución de fase
COMPUTE_SRC  = plot_dprob_dfase/compute_billiard.cpp	#codigo fuente
COMPUTE_BIN  = compute_billiard							#binario
DENSITY_DAT  = density.dat								#datos de densidad
PHASE_DAT    = phase.dat								#datos de fase
PLOTDP_PY    = plot_dprob_dfase/plot_results.py			#script grafica

# Parámetros por defecto (se pueden sobrescribir en la invocación)
xi0  ?= 3.0
eta0 ?= 2.0
k2   ?= 0.406
angle ?= 0.0
threads ?= 4

.PHONY: all plotb plot clean

# Por defecto solo compila los binarios y vuelca la frontera
all: $(BOUNDARY_BIN) $(BOUNDARY_DAT) $(COMPUTE_BIN)

# 1) Frontera								Compila el codigo
$(BOUNDARY_BIN): $(BOUNDARY_SRC)
	$(CXX) $(CXXFLAGS) $< -o $@

# Volcar la frontera pasando xi0 y eta0		Ejecuta el binario para generar datos
$(BOUNDARY_DAT): $(BOUNDARY_BIN)
	./$(BOUNDARY_BIN) $(xi0) $(eta0) > $(BOUNDARY_DAT)

# 2) Compute (densidad/fase)				Compila codigo
$(COMPUTE_BIN): $(COMPUTE_SRC)
	$(CXX) $(CXXFLAGS) $< -o $@

# Graficar solo frontera. Ej: make plotb xi0=3.0 eta0=2.0
plotb: $(BOUNDARY_BIN)
	@echo "Regenerando boundary.dat con xi0=$(xi0), eta0=$(eta0)"
	./$(BOUNDARY_BIN) $(xi0) $(eta0) > $(BOUNDARY_DAT)
	python3 $(BOUNDARY_PY) 

# Graficar densidad y fase. Ej: make plot xi0=3.0 eta0=2.0 k2=0.406
plot: all
	./$(COMPUTE_BIN) $(xi0) $(eta0) $(k2) $(angle) $(threads)
	@echo "Ejecución de $(threads) threads completada"
	python3 $(PLOTDP_PY) $(xi0) $(eta0) $(k2) $(angle) $(threads)

# Graficar espectro de T(k)
SPEC_SRC       = plot_spectrum/scan_spectrum.cpp	#codigo fuente resonancias
SPEC_BIN       = scan_spectrum						#binario
PLOTSPEC_PY    = plot_spectrum/plot_spectrum.py		#script grafica

$(SPEC_BIN): $(SPEC_SRC)
	$(CXX) $(CXXFLAGS) $< -o $@

plotspec: $(SPEC_BIN) #Ej: make plotspec xi0=3.0 eta0=2.0
	./$(SPEC_BIN) $(xi0) $(eta0)
	python3 $(PLOTSPEC_PY)

#Grafica comparativa variando el tamano de malla con optimizers: make run_optimized
run_optimized: 
	@echo "Ejecutando pruebas optimizadas (CPU+GPU)..."
	cd ComputationalTimes && chmod +x tiempos_Opt.sh && ./tiempos_Opt.sh
	cd ../

#Grafica comparativa variando el tamano de malla sin optimizers: make run_not_optimized
run_not_optimized: 
	@echo "Ejecutando pruebas no optimizadas (CPU+GPU)..."
	cd ComputationalTimes && chmod +x tiempos_noOpt.sh && ./tiempos_noOpt.sh
	cd ../

omp_times:
	@echo "Evaluacion de threads de OMP"
	cd ComputationalTimes && chmod +x ./tiempos_threads.sh
	cd ../

animation_k:
	@echo "Ejecutando video variando numero de onda ..."
	cd Videos && chmod +x video_k.sh && ./video_k.sh
	cd ../



clean:
	rm -f $(BOUNDARY_BIN) $(BOUNDARY_DAT) \
	      $(COMPUTE_BIN) $(DENSITY_DAT) $(PHASE_DAT) \
	      $(SPEC_BIN) spectrum.dat resonances.dat \
		  ComputationalTimes/compute_billiard_cpp.o \
          ComputationalTimes/compute_billiard_cuda.o \
          ComputationalTimes/tiempos_cpp.txt \
          ComputationalTimes/tiempos_cuda.txt \
          ComputationalTimes/boundary.dat

