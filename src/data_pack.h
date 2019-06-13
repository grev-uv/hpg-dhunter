#ifndef DATA_PACK_H
#define DATA_PACK_H

#include <deque>
#include <string>

using namespace std;

// estructura de datos para transferir entre visualizador, controlador y GPU.

struct datos_cuda
{
    float      **mc_full;       // matriz de datos completos de todas las muestras
    float      **h_haar_C;      // matriz de datos procesados en GPU
    deque<int> h_haar_L;        // vector con número de datos por nivel
    float      *d_haar;         // vector de datos en GPU
    float      *d_aux;          // vector de ayuda a mantener la sincronizació
    size_t     pitch;           // ajuste óptimo de memoria GPU para datos de cada muestra
    size_t     pitch_2;         // ajuste óptimo de memoria GPU para auxiliar
    size_t     sample_num;      // número de datos por muestra
    int        samples;         // número de muestras a trasnformar
    int        levels;          // número de niveles a transformar
    int        data_adjust;     // ajuste desfase en división por nivel para número impar de datos
    size_t     rango_inferior;  // límite inferior ventana de datos a transformar
    size_t     rango_superior;  // límite superior ventana de datos a transformar
    void       *d_glPtr;        // puntero a array de datos para visualización directa opengl desde gpu
    string     **refGen;        // matriz de referencias cromosómicas del cromosoma analizado
    float      *d_max;          // vector con el valor máximo a dibujar en la posición [0]
};

#endif // DATA_PACK_H
