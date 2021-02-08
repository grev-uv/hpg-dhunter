/*
*  haar_vX is the well synchronized version of haar transform into de GPU.
*  This version work with all the samples as a matrix into de GPU
*  with dimension SAMPLES x (sample_num + data_adjust) (rows x cols)
*  Copyright (C) 2018 Lisardo Fernández Cordeiro <lisardo.fernandez@uv.es>
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2, or (at your option)
*  any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
*
*/
 
/** \file
*  \brief Archivo para procesamiento de diferentes muestras metiladas de ADN.
*
*  Este archivo contiene la definición de las funciones para:
*         ..carga de datos en GPU
*         ..lanzamiento de proceso de transformación en GPU
*         ..kernel en GPU para control de transformación en niveles definidos
*         ..kernel de transformación wavelet del vector seleccionado
*         ..kernel para copiar coeficientes desde vector auxiliar de sincronización
*/

#include <stdio.h>
#include <GL/gl.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include "data_pack.h"

#define BLOCK_SIZE  1024		// número de hilos por bloque de GPU
#define AJUSTE_PLOT 1//.70        // ajusta eje Y de gráfica a AJUSTE_PLOT
#define DESPLAZAMIENTO_DIBUJO 0.97
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); } // para gestión de errores en GPU


#define DESPLAZAMIENTO_DIBUJO_X 1           //0.97
#define DESPLAZAMIENTO_DIBUJO_Y 0//0.97        //1
#define WAVELET_COEF            4// coeficientes para wavelet bior-3.1


/** ***********************************************************************************************
  * \fn void gpuAssert(cudaError_t, char*, int, bool)
  *  \brief Función responsable de recoger error en GPU y mostrarlo
  *  \param code	código de error de la GPU
  *  \param *file	fichero donde se produce el error
  *  \param line	línea de código donde se produce el error
  *  \param abort	indica si se sale del programa
  * ***********************************************************************************************
  */
extern "C"
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
	if (code != cudaSuccess)
	{
        fprintf(stderr, "GPUassert: %s %s %d\n\n", cudaGetErrorString(code), file, line);
		if (abort)
			exit(code);
	}
}

/** ***********************************************************************************************
  * \fn void copyValuesTotal(float*, float *, int, int)
  *  \brief función "hija" en GPU responsable de la copia de los datos del segmento a transformar
  *         proporcionando sincronización a nivel GRID
  *  \param *haar	puntero a vector de datos original
  *  \param *aux	puntero a vector de datos a transformar
  *  \param num		numero datos totales a transformar
  *  \param pi  	posicion inicial de copia - offset
  * ***********************************************************************************************
  */
extern "C"
__global__
void copyValuesTotal(float *haar, float *aux, int num, int posicion_inicial)
{
	// variables ------------------------------------------------------------------------------
	int index = threadIdx.x + blockIdx.x * blockDim.x;	// índice sobre todo el vector

	// copiar todos los valores de haar en aux
	if (index < num)
        	aux[index] = haar[index + posicion_inicial];
}


/** ***********************************************************************************************
  * \fn void copyValues(float*, float *, int)
  *  \brief función "hija" en GPU responsable de la copia de los valores escalados
  *         proporcionando sincronización a nivel GRID
  *  \param *aux	puntero a vector de datos a transformar
  *  \param *temp	puntero a vector de datos temporales a copiar
  *  \param num     numero datos a copiar
  * ***********************************************************************************************
  */
extern "C"
__global__
void copyValues(float *aux, float *temp, int num)
{
    // variables ----------------------------------------------------------------------------------
    int index = threadIdx.x + blockIdx.x * blockDim.x;	// índice sobre todo el vector

    // copiar todos los valores de haar en aux
    if (index < num)
        aux[index] = temp[index];
}


/** ***********************************************************************************************
  * \fn void transform(float*, int, int, int)
  *  \brief función "hija" en GPU responsable de la transformación wavelet de un vector
  *         proporcionando sincronización a nivel GRID.
  *  \param *aux	puntero a vector de datos a transformar
  *  \param *temp	puntero a vector de resultados intermedios
  *  \param num		número de posiciones del vector
  * ***********************************************************************************************
  */
extern "C"
__global__
void transform(float *aux, float *temp, int num)
{
    // variables ----------------------------------------------------------------------------------
    int index = threadIdx.x + blockIdx.x * blockDim.x;	// índice sobre todo el vector
/*    float f   = 0.7071067811865476;                     // coeficiente haar wavelet
    float aux1;                                         // variables auxiliares de sincronización
    int idx;                                            // indice auxiliar para guardar dato

    // transformada haar en paralelo sobre el vector recibido -------------------------------------
    if (index < num)
    {
        if ((index & 0x01) == 0)	// solo los hilos con índice par (0, 2, 4, ...)
        {
            idx = index * 0.5;

            aux1 = (aux[index] + aux[index + 1]) * f;	// escalado (filtro paso-bajo)

            temp[idx]  = aux1;
        }
    }
*/


    float f[WAVELET_COEF];
    float bior31[4]  = {-0.3535533905932738, 1.0606601717798214, 1.0606601717798214, -0.3535533905932738};
    float spline[5]  = {-0.1767766953, 0.3535533906, 0.7071067812, 0.3535533906, -0.1767766953};
    float bior33[8]  = {0.06629126073623884, -0.19887378220871652, -0.15467960838455727, 0.9943689110435825,
                        0.9943689110435825, -0.15467960838455727, -0.19887378220871652, 0.06629126073623884};
    float bior35[12] = {-0.013810679320049757, 0.04143203796014927, 0.052480581416189075, -0.26792717880896527,
                       -0.07181553246425874, 0.966747552403483, 0.966747552403483, -0.07181553246425874,
                       -0.26792717880896527, 0.052480581416189075, 0.04143203796014927, -0.01381067932004975};

    switch (WAVELET_COEF)
    {
    case 4:
        for (int i = 0; i < WAVELET_COEF; i++)
            f [i] = bior31 [i];
        break;
    case 5:
        for (int i = 0; i < WAVELET_COEF; i++)
            f [i] = spline [i];
        break;
    case 8:
        for (int i = 0; i < WAVELET_COEF; i++)
            f [i] = bior33 [i];
        break;
    case 12:
        for (int i = 0; i < WAVELET_COEF; i++)
            f [i] = bior35 [i];
        break;
    default:
        ;
    }

    float aux1 = 0.0;                                      // variable auxiliar de almacenamiento intermedio de resultado
    int idx;                                               // indice auxiliar para guardar dato

    if (index < num + (WAVELET_COEF - 2))
    {
        if ((index & 0x01) == 0)
        {
            idx = floorf(index * 0.5);

            for (int i = -(WAVELET_COEF - 2); i < 2; i++)
            {
                if (index + i < 0)
                    aux1 += 0.0;
                else if (index + i > num)
                    aux1 += 0.0;
                else
                    aux1 += aux[index + i] * f[i + (WAVELET_COEF - 2)];
            }

            temp[idx] = aux1;
        }
    }

}

/** ***********************************************************************************************
  * \fn void array2Plot(float*, float*, int)
  *  \brief función "hija" en GPU responsable de crear el array para ploteado.
  *  \param *aux	puntero a vector de datos tranformados
  *  \param *glPtr	puntero a vector de datos para plotear
  *  \param num		número de datos a plotear
  *  \param max     valor máximo de cálculo para escalar todos los valores entre 0 y 2
  *  \param hilo    hilo que gestiona el ploteado
  * ***********************************************************************************************
  */
extern "C"
__global__
void array2Plot(float *aux, float *glPtr, int num, float *max, int hilo)
{
    // variables ----------------------------------------------------------------------------------
    int index = threadIdx.x + blockIdx.x * blockDim.x;	// índice sobre todo el vector
    int idx;

    // copia todos los valores de aux en glPtr duplicando las posiciones anexas para crear una
    // gráfica en escalón propia de la transformada wavelet haar
    if (index < num)
    {
        idx = hilo * num + index;

        glPtr[idx * 4]     = (index * 2.0 / num) - DESPLAZAMIENTO_DIBUJO_X;     // eje x de -1 a 1
        glPtr[idx * 4 + 1] = aux[index] / max[0] * AJUSTE_PLOT - DESPLAZAMIENTO_DIBUJO_Y;              // eje y de 0 a 2
        glPtr[idx * 4 + 2] = ((index+1) * 2.0 / num) - DESPLAZAMIENTO_DIBUJO_X; // eje x de -1 a 1
        glPtr[idx * 4 + 3] = aux[index] / max[0] * AJUSTE_PLOT - DESPLAZAMIENTO_DIBUJO_Y;              // eje y de 0 a 2

    }
}


/** ***********************************************************************************************
  * \fn void maxVal(float*, float*, int, int)
  *  \brief función "hija" en GPU responsable de encontrar el máximo valor de cálculo
  *                 por método de reducción con memoria local a nivel de hilo, encontrando
  *                 el máximo valor por bloque.
  *  \param *aux_c	puntero a vector de datos
  *  \param *max	puntero a vector de maximos encontrados por bloque
  *  \param num		número de datos
  * ***********************************************************************************************
  */
extern "C"
__global__
void maxVal(float *aux_c, float *max, int num)
{
    extern __shared__ float sdata[];

    // para aprovechar el total de los hilos en la primera operación de búsqueda
    // se direcciona al doble de la capacidad de un bloque (blockDim.x * 2)
    unsigned int index = blockIdx.x * (blockDim.x * 2) + threadIdx.x;
    unsigned int tid = threadIdx.x;

    if (index < num)
    {
        // la primera carga de datos a la memoria local, se realiza buscando el máximo
        // de cada parte de los datos direccionados
        if (aux_c[index] >= aux_c[index + blockDim.x])
            sdata[tid] = aux_c[index];
        else
            sdata[tid] = aux_c[index + blockDim.x];

        __syncthreads();

        // a partir de la carga condicionada de los datos en memoria locas,
        // se busca el máximo del bloque por reducción
        for (unsigned int i = blockDim.x / 2; i > 0; i >>= 1)
        {
            if (tid < i) // && (tid + i + (blockDim.x * blockIdx.x)) < index) // / 2)
                if (sdata[tid] < sdata[tid + i])
                    sdata[tid] = sdata[tid + i];

            __syncthreads();
        }

        // guarda el resultado de cada bloque para el siguiente paso
        if (tid == 0)
            max[blockIdx.x] = sdata[0];
    }
}

/** ***********************************************************************************************
  * \fn void maxGlobal(float*, int, int)
  *  \brief función "hija" en GPU responsable de encontrar el máximo valor de cálculo
  *                 por método de reducción con memoria local a nivel de hilo, encontrando
  *                 el máximo valor entre los máximos encontrados de cada muestra.
  *  \param *temp	puntero a vector de máximos
  *  \param pitch	número de bytes por muestra reservados para el vector temporal
  *  \param samples número de muestras
  * ***********************************************************************************************
  */
extern "C"
__global__
void maxGlobal(float *temp, size_t pitch, int samples)
{
    // variables ----------------------------------------------------------------------------------
    int index = threadIdx.x;

    // solo un hilo se encarga de buscar el máximo entre los máximos encontrados para cada muestra
    if (index == 0)
    {
        for (int i = 1; i < samples; i++)
        {
            if (temp[0] < temp[i * pitch / sizeof(float)])
                temp[0] = temp[i * pitch / sizeof(float)];
        }
    }
}


/** ***********************************************************************************************
  * \fn void wavedec(float*, float**, int, int, int, int, int)
  *  \brief Función principal en GPU responsable de calcular y ordenar las partes del vector
  *         para su transformación wavelet multinivel.
  *  \param *haar	puntero a matriz de datos a transformar
  *  \param *aux	puntero a matriz de coeficiente auxiliares para ayuda a la sincronización
  *  \param *temp	puntero a matriz de cálculos temporales de ayuda a la sincronización
  *  \param pitch	desplazamiento óptimo en memoria GPU para alojar cada muestra	
  *  \param pitch_2	desplazamiento óptimo en memoria GPU para alojar cálculo auxiliar
  *  \param pitch_3	desplazamiento óptimo en memoria GPU para alojar cálculo temporal
  *  \param n		número total de posiciones del vector
  *  \param l		número de niveles a computar
  *  \param samples número de muestras a analizar
  *  \param pi      posición inicial del segmento de datos a analizar
  *  \param *glPtr  puntero a array de datos para dibujar con openGL
  * ***********************************************************************************************
  */
extern "C"
__global__ 
void wavedec(float *haar, float *aux, float *temp,
             size_t pitch, size_t pitch_2, size_t pitch_3,
             int n, int l, int samples, int pi, float *glPtr)
{
    // variables ----------------------------------------------------------------------------------
	int index_X = threadIdx.x + blockIdx.x * blockDim.x;	// indice de hilos sobre todo el vector
    int level   = 0;                                        // número de nivel
    int num     = n;                                        // número de posiciones en vector
    int hilo;                                               // guarda el hilo asignado para que se resposabilice de todo el proceso
    int nume    = 0;                                        // numero de posiciones antes de ajuste para cálculo de nuevo nivel
                                                            //     evita que los datos de ploting hagan rayas por descuadre en posición


    // limita el número de hilos al de muestras ---------------------------------------------------
	if (index_X < samples)
	{
		hilo = index_X;		// cada hilo se responsabiliza de una misma muestra

		if (hilo == index_X)		
		{
            // separar los datos por muestras - - - - - - - - - - - - - - - - - - - - - - - - - - -
			float *haar_c = (float *)((char *)haar + index_X * pitch);
            float *aux_c  = (float *)((char *)aux  + index_X * pitch_2);
            float *temp_c = (float *)((char *)temp + index_X * pitch_3);

			__syncthreads(); 

            // llamada a función hija para copiar segmento de vector a transformar
            copyValuesTotal<<<(num + BLOCK_SIZE-1) / BLOCK_SIZE, BLOCK_SIZE>>>(haar_c,
                                                                               aux_c,
                                                                               num,
                                                                               pi);


            // procesamiento multinivel del vector de datos ---------------------------------------
			// repite la transformación tantas veces como niveles se han solicitado
            while (level < l && num >= 2)
            {
				// llamada a función hija para transformación del nivel correspondiente
				// con esta división en padre-hijo, se consigue sincronizar cada nivel 
				// \param	<<<((datos_x_muestra + num_hilos_bloque-1) / num_hilos_bloque),
				// 		numero hilos por bloque>>>
                transform<<<(num + BLOCK_SIZE-1) / BLOCK_SIZE, BLOCK_SIZE>>>(aux_c,
                                                                             temp_c,
                                                                             num);


                // actualizar variables de nivel  - - - - - - - - - - - - - - - - - - - - - - - - -
                level += 1;
                num    = ceilf(num * 0.5);


                // llamada a función hija para copiar resultados en vector auxiliar
                copyValues<<<(num + BLOCK_SIZE-1) / BLOCK_SIZE, BLOCK_SIZE>>>(aux_c,
                                                                              temp_c,
                                                                              num);


                // actualiza el número de datos para el siguiente nivel - - - - - - - - - - - - - -
                nume = num;
                if ((num & 01) == 1)
                {
                    num++;
                    aux_c[num] = 0;
                }
            }



            // hallar el valor máximo de todas las muestras para ajustar ploteado, para ello:
            // busca el máximo donde están los datos calculados (aux_c) y aprovecha para pasarlos a temp_c
            maxVal<<<(nume + BLOCK_SIZE-1) / BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE * sizeof(float)>>>(aux_c,
                                                                                                   temp_c,
                                                                                                   nume);

            // una vez los datos en temp_c, si el número de datos es mayor que la capcidad de un bloque,
            // repite la operación anterior pero dentro de temp_c, dividiendo el tramo en dos
            float maximos = ceilf(nume * 1.0 / BLOCK_SIZE);
            while (maximos > 1.0)
            {
                // halla el máximo por parejas almacenando el resultado en la mitad superior
                maxVal<<<(maximos + BLOCK_SIZE-1) / BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE * sizeof(float)>>>(temp_c,
                                                                                                          &temp_c[pitch_3 / (sizeof(float) * 2)],
                                                                                                          ceilf(maximos));
                __syncthreads();

                maximos = ceilf(maximos / BLOCK_SIZE);

                // copia los datos en la mitad inferior para volver a reducirlos
                copyValues<<<(maximos + BLOCK_SIZE-1) / BLOCK_SIZE, BLOCK_SIZE>>>(temp_c,
                                                                                  &temp_c[pitch_3 / (sizeof(float) * 2)],
                                                                                  ceilf(maximos));

            }

            __syncthreads();

            // una vez obtenido el máximo por cada muestra, se busca el máximo de todos los valores
            if (hilo == 0)
                maxGlobal<<< 1, 1 >>>(temp,
                                      pitch_3,
                                      samples);

            __syncthreads();

            // rellenar el array con datos para dibujar desde openGL
            // llamada a función hija para rellenar datos de gráfica
            array2Plot<<<(num + BLOCK_SIZE-1) / BLOCK_SIZE, BLOCK_SIZE>>>(aux_c,
                                                                          glPtr,
                                                                          nume,
                                                                          temp,
                                                                          hilo);
        }
	}
}

/** ***********************************************************************************************
  * \fn void cuda_send_data(datos_cuda &)
  *  \brief Función para enviar los datos a la GPU
  *  \param &cuda_data  estructura con variables de control de datos
  * ***********************************************************************************************
  */
void cuda_send_data(datos_cuda &cuda_data)
{
    // reserva espacio en GPU para el vector a transformar y copia matriz de datos ----------------
    // devuelve valor de desplazamiento (pitch) óptimo para gestión de memoria adecuada
    // en función de la cantdad de datos a alojar
    // \param 	puntero a posición memoria GPU,
    //          desplazamiento óptimo devuelto por CUDA,
    //          cantidad de bytes a reservar por fila,
    //          número de muestras (filas)
    gpuErrchk(cudaMallocPitch(&cuda_data.d_haar,
                              &cuda_data.pitch,
                              (cuda_data.sample_num + cuda_data.data_adjust) * sizeof(float),
                              cuda_data.samples));


    gpuErrchk(cudaMallocPitch(&cuda_data.d_aux,
                              &cuda_data.pitch_2,
                              (cuda_data.sample_num + cuda_data.data_adjust) * sizeof(float),
                              cuda_data.samples));


    // envío de datos a GPU -----------------------------------------------------------------------
    // \param	puntero a posición de memoria GPU,
    //          desplazamiento óptimo,
    //          puntero a posición de datos en CPU a enviar a GPU,
    //          cantidad de bytes a enviar por muestra,
    //          cantidad de bytes a alojar por muestra,
    //          número de filas (muestras)
    gpuErrchk(cudaMemcpy2D( cuda_data.d_haar,
                            cuda_data.pitch,
                            cuda_data.mc_full[0],
                            cuda_data.sample_num * sizeof(float),
                            cuda_data.sample_num * sizeof(float),
                            cuda_data.samples,
                            cudaMemcpyHostToDevice));

}


/** ***********************************************************************************************
  * \fn void cuda_main(datos_cuda &)
  *  \brief Función para procesar los datos en la GPU
  *  \param &cuda_data  estructura con variables de control de datos
  * ***********************************************************************************************
  */
void cuda_main(datos_cuda &cuda_data)
{
    // reserva TODA la memoria CONTIGUA para la matriz de muestras tranformadas -------------------
    // para trasvase de datos entre GPU y CPU con CUDA, la matriz debe ser contigua completa
    cuda_data.h_haar_C = new float*[cuda_data.samples];                             // reservar punteros a filas
    cuda_data.h_haar_C[0] = new float[cuda_data.samples * (cuda_data.h_haar_L[0])];	// reservar toodos los datos (rows * cols)
    for (int i = 1; i < cuda_data.samples; i++)                                     // asignar valor a punteros de fila
        cuda_data.h_haar_C[i] = cuda_data.h_haar_C[i - 1] + cuda_data.h_haar_L[0];


    // reserva memoria para cálculos temporales en GPU --------------------------------------------
    float *d_temp;
    size_t pitch;
    gpuErrchk(cudaMallocPitch(&d_temp,
                              &pitch,
                              (cuda_data.sample_num + 1) * sizeof(float) * 0.7,
                              cuda_data.samples));


    // transforma el número de muestras elegida ---------------------------------------------------
	// realiza la transformación en la GPU del conjunto de muestras cargado
	// \param	<<< número de bloques a utilizar,
    //          número de hilos por bloque >>> (máximo 1024 para PASCAL GTX 1080)
	// \param	puntero a datos a transformar alojados en GPU,
    //          desplazamiento óptimo de datos por fila,
    //          número de datos por muestra (fila) a transformar,
    //          ajuste de longitud de muestra por número impar al dividir la muestra
    wavedec<<<1, cuda_data.samples>>>(cuda_data.d_haar,
                                      cuda_data.d_aux,
                                      d_temp,
                                      cuda_data.pitch,
                                      cuda_data.pitch_2,
                                      pitch,
                                      cuda_data.sample_num,
                                      cuda_data.levels,
                                      cuda_data.samples,
                                      cuda_data.rango_inferior,
                                      (float *)cuda_data.d_glPtr);

    // espera a que la GPU termine el trabajo - - - - - - - - - - - - - - - - - - - - - - - - - - -
    gpuErrchk(cudaDeviceSynchronize());


    // recupera el resultado de la transformación en memoria GPU a memoria CPU- - - - - - - - - - -
	// \param	puntero a matriz de datos a guardar en CPU,
    //          cantidad de bytes a guardar por muestra,
    //          puntero a datos para copiar de GPU,
    //          desplazamiento óptimo de datos por fila en GPU,
    //          cantidad de bytes en GPU a copiar por muestra,
    //          número de muestras (filas)
    gpuErrchk(cudaMemcpy2D(	cuda_data.h_haar_C[0],
                            cuda_data.h_haar_L[0] * sizeof(float),
                            cuda_data.d_aux,
                            cuda_data.pitch,
                            cuda_data.h_haar_L[0] * sizeof(float),
                            cuda_data.samples,
                            cudaMemcpyDeviceToHost));

    gpuErrchk(cudaMemcpy( cuda_data.d_max,
                          d_temp,
                          cuda_data.samples * sizeof(float),
                          cudaMemcpyDeviceToHost));


    //libera la memoria temporal utilizada para cálculos intemedios
    cudaFree(d_temp);
}

/** ***********************************************************************************************
  * \fn void *cuda_init()
  *  \brief Función para inicializar la gpu
  * ***********************************************************************************************
  */
void cuda_init()
{
    int deviceCount = 0;
    int cudaDevice  = 0;
    char cudaDeviceName [100];
    cudaDeviceProp prop;
    cuInit(0);
    cuDeviceGetCount(&deviceCount);
    cuDeviceGet(&cudaDevice, 0);
    cuDeviceGetName(cudaDeviceName, 100, cudaDevice);
    cudaGetDeviceProperties(&prop, cudaDevice);

    if (cudaChooseDevice(&cudaDevice, &prop) != cudaSuccess)
        puts("failed to choose device");
    if (cudaGLSetGLDevice(cudaDevice) != cudaSuccess)
        puts("failed to set gl device");

    printf("Number of devices: %u \t cuda device: %d\n", deviceCount, cudaDevice);
    printf("Device name: %s\n", cudaDeviceName);
    printf("Warp size: %u\n", prop.warpSize);
}

/** ***********************************************************************************************
  * \fn void cuda_end(data buf)
  *  \brief Función para liberar memoria de la GPU
  *  \param &cuda_data  estructura con variables de control de datos
  * ***********************************************************************************************
  */
void cuda_end(datos_cuda &cuda_data)
{
    //libera la memoria de la gpu utilizada para cálculos intemedios
    cudaFree(cuda_data.d_haar);
    cudaFree(cuda_data.d_aux);
}

/** ***********************************************************************************************
  * \fn void *cuda_registerBuffer(GLuint buf)
  *  \brief Función para registrar el vínculo de cuda con opengl
  *  \param buf buffer donde se alojan los datos para procesar
  * ***********************************************************************************************
  */
void *cuda_registerBuffer(GLuint buf)
{
    cudaGraphicsResource *res = 0;

    if (cudaGraphicsGLRegisterBuffer(&res, buf, cudaGraphicsRegisterFlagsNone) != cudaSuccess)
        printf("Fallo en el registro del buffer %u\n", buf);

    return res;
}

/** ***********************************************************************************************
  * \fn void cuda_unregisterBuffer(void *res)
  *  \brief Función para desvincular opengl de cuda
  *  \param *res    referencia al vínculo
  * ***********************************************************************************************
  */
void cuda_unregisterBuffer(void *res)
{
    if (cudaGraphicsUnregisterResource((cudaGraphicsResource *) res) != cudaSuccess)
        puts("Fsllo eliminando el registro de recursos para el buffer");
}

/** ***********************************************************************************************
  * \fn void *cuda_map(void *res)
  *  \brief Función para mapear los datos cuda sobre opengl
  *  \param *res    variable de vínculo con los datos cuda
  * ***********************************************************************************************
  */
void *cuda_map(void *res)
{
    if (cudaGraphicsMapResources(1, (cudaGraphicsResource **) &res) != cudaSuccess, 0)
    {
        puts("Fallo en el mapeado de recursos");
        return 0;
    }

    void *devPtr = NULL;

    size_t size;

    if (cudaGraphicsResourceGetMappedPointer(&devPtr, &size, (cudaGraphicsResource *) res) != cudaSuccess)
    {
        puts("Fallo en la adquisición del puntero del dispositivo ");
        return 0;
    }

    return devPtr;
}

/** ***********************************************************************************************
  * \fn void cuda_unmap(void *res)
  *  \brief Función para liberar el mapeo de datos
  *  \param *res    variable de vínculo con los datos cuda
  * ***********************************************************************************************
  */
void cuda_unmap(void *res)
{
    if (cudaGraphicsUnmapResources(1,(cudaGraphicsResource **) &res, 0) != cudaSuccess)
        puts("Fallo en el desmapeado de recursos");
}

/** ***********************************************************************************************
  * \fn void calculo_haar_L(datos_cuda &cuda_data)
  *  \brief Función para calcular el número de datos en el nivel dado y el ajuste por impares
  *  \param &cuda_data  estructura con variables de control de datos
  * ***********************************************************************************************
  */
void cuda_calculo_haar_L(datos_cuda &cuda_data)
{
    // cálculo de número de coeficientes por nivel y del ajuste de paso entre escala y coeficiente
    cuda_data.h_haar_L.push_front(cuda_data.sample_num);	// última posición guarda el total de posiciones por muestra

    // para cada nivel se divide por dos la cantidad de posiciones del nivel anterior -------------
    // redondeando al alza y actualizando el ajuste cuando sea impar
    for (int fila = cuda_data.levels; fila > 0; fila--)
    {
        if (ceil(cuda_data.h_haar_L.front() * 0.5 >= 2))
        {
            cuda_data.h_haar_L.push_front(ceil(cuda_data.h_haar_L.front() * 0.5));
            if (fila > 0 && cuda_data.h_haar_L[1] != cuda_data.sample_num)
                cuda_data.data_adjust += size_t(2 * cuda_data.h_haar_L.front() - cuda_data.h_haar_L[1]);
        }
        else
            break;
    }
    cuda_data.h_haar_L.push_front(cuda_data.h_haar_L.front());	// primera posición coincide con el número de datos de escala
}
