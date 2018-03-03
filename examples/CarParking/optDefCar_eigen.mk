# Here you can specify additional compiler flags
USER_FLAGS:= -O2 -ggdb -DEIGEN_NO_MALLOC -DEIGEN_NO_DEBUG -DEIGEN_DONT_PARALLELIZE -ffast-math -static-libstdc++ -static-libgcc -static -march=native -mfpu=neon -mfloat-abi=hard 
#-fopenmp -pthread