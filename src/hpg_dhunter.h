/*
*  HPG_Dhunter is the main class to control and manipulate data files
*  Copyright (C) 2018 Lisardo Fernández Cordeiro <lisardo.fernandez@uv.es>
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 3, or (at your option)
*  any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
*  or see <https://www.gnu.org/licenses/>.
*
*/

/** \file
*  \brief Programa para procesamiento y visualización de diferentes
*         muestras metiladas de ADN.
*
*  Este archivo contiene la definición de las funciones para:
*         ..declaración de funciones externas de procesamiento en GPU
*         ..control de rango de datos a analizar y visualizar
*         ..selección de fichero a analizaar
*/

#ifndef HPG_DHUNTER_H
#define HPG_DHUNTER_H

#include <QMainWindow>
#include <QTextCursor>
#include <QThread>
#include <QMutex>
#include "data_pack.h"
#include "ogl_graphic.h"
#include "files_worker.h"
#include "refgen.h"
#include <cuda_runtime.h>
#include <cuda.h>
#include <chrono>

// temporizadores para control de tiempos de los procesos de análisis
#define TIMING

#ifdef TIMING
#define INIT_TIMER        auto start = std::chrono::high_resolution_clock::now();
#define START_TIMER       start = std::chrono::high_resolution_clock::now();
#define STOP_TIMER(name)  qDebug() << "DURACION de " << name << ": " << \
                          std::chrono::duration_cast<std::chrono::milliseconds>( \
                          std::chrono::high_resolution_clock::now()-start \
                          ).count() << " ms ";
#else
#define INIT_TIMER
#define START_TIMER
#define STOP_TIMER(name)
#endif

#define DMR_THRESHOLD   0.3 // valor inicial para umbral de cálculo de DMRs



using namespace std;


/** ***********************************************************************************************
  *  \brief declaración de funciones externas (en haar_vX.cu) para compilación con nvcc
  *         funciones en GPU de cálculo de transformada y preparaadción de resultado para dibujado
  *  \fn    void cuda_send_data(datos_cuda &)
  *  \fn    void cuda_main(datos_cuda &)
  *  \fn    void cuda_calculo_haar_L(datos_cuda &)
  *  \fn    void cuda_init()
  *  \fn    void cuda_end(datos_cuda &)
  * ***********************************************************************************************
  */
//extern
void cuda_send_data(datos_cuda &);
//extern
void cuda_main(datos_cuda &);
//extern
void cuda_calculo_haar_L(datos_cuda &);
//extern
void cuda_init();
//extern
void cuda_end(datos_cuda &);

namespace Ui {
class HPG_Dhunter;
}

class HPG_Dhunter : public QMainWindow
{
    Q_OBJECT

public:
    explicit HPG_Dhunter(QWidget *parent = 0);
    ~HPG_Dhunter();

    /** ***********************************************************************************************
      *  \brief variable string responsable de capturar/guardar el nombre de/los fichero a analizar
      *  \param fichero     string con nombre fichero seleccionado por usuario
      *  \param ficheros    stringlist con nombre de todos los ficheros seleccionados
      * ***********************************************************************************************
      */
    QString fichero;
    QStringList ficheros_case;
    QStringList ficheros_control;
    QStringList parametros;

    /** ***********************************************************************************************
      *  \brief variable de tipo estructura con variables para control de datos a analizar
      *  \param cuda_data   estrutura con variables de control y datos
      * ***********************************************************************************************
      */
    datos_cuda cuda_data;

private slots:
    /** ***********************************************************************************************
      * \fn void on_open_file_clicked()
      *  \brief Función responsable de abrir explorador de archivos y capturar fichero a analizar
      * ***********************************************************************************************
      */
    void on_open_file_clicked();

    /** ***********************************************************************************************
      * \fn void on_scroll_adn_valueChanged(int value)
      *  \brief Función responsable capturar movimiento de la barra de scroll
      *         actualizando rango de posiciones a analizar y visualizar
      *  \param value	posición de rango inferior de datos a transformar y analizar
      * ***********************************************************************************************
      */
    void on_scroll_adn_valueChanged(int value);

    /** ***********************************************************************************************
      * \fn void on_page_valueChanged(int position)
      *  \brief Función responsable de capturar movimiendo de deslizador de ancho de ventana de datos
      *  \param position	ancho de ventana de datos a transformar y mostrar
      * ***********************************************************************************************
      */
    void on_page_valueChanged(int position);

    /** ***********************************************************************************************
      * \fn void on_analiza_clicked()
      *  \brief Función responsable de lanzar una solicitud de transformación
      * ***********************************************************************************************
      */
    void on_analiza_clicked();

    /** ***********************************************************************************************
      * \fn void on_rango_inf/superior_editingFinished()
      *  \brief Función responsable de actualizar los rangos de la ventana de datos a transformar
      *         una vez se ha introducido alguno de los límites, inf o sup, de la ventana
      * ***********************************************************************************************
      */
    void on_rango_inferior_editingFinished();
    void on_rango_superior_editingFinished();

    /** ***********************************************************************************************
      * \fn void on_slider_nivel_sliderReleased()
      *  \brief Función responsable de actualizar el número de niveles en la estructura de datos
      *         y de solicitar una nueva tranformación
      * ***********************************************************************************************
      */
    void on_slider_nivel_sliderReleased();
    void on_slider_nivel_valueChanged(int value);

    /** ***********************************************************************************************
      * \fn void on_scroll_adn_sliderReleased()
      *  \brief Función responsable de solicitar una nueva tranformación cuando se suelta el slider
      * ***********************************************************************************************
      */
    void on_scroll_adn_sliderReleased();
    void on_page_sliderReleased();

    /** ***********************************************************************************************
      * \fn void on_threshold_valueChanged(int value)
      *  \brief Función responsable de ajustar el umbral de búsqueda de DMRs a la diferencia seleccionada
      * ***********************************************************************************************
      */
    void on_threshold_valueChanged(int value);
    void on_threshold_sliderReleased();

    /** ***********************************************************************************************
      * \fn void on_dmrs_clicked()
      *  \brief Función responsable de realizar la primera búsqueda de DMRs en el segmento completo
      * ***********************************************************************************************
      */
    void on_dmrs_clicked();

    /** ***********************************************************************************************
      * \fn void on_dmr_position_cursorPositionChanged()
      *  \brief Función responsable de mostrar en gráfica la DMRs seleccionada en ventana
      * ***********************************************************************************************
      */
    void on_dmr_position_cursorPositionChanged();

    /** ***********************************************************************************************
      * \fn void on_fine_tunning_XXX(int value)
      *  \brief Funciones responsables del ajuste fino del segmento a analizar y visualizar
      * ***********************************************************************************************
      */
    void on_fine_tunning_valueChanged(int value);
    void on_fine_tunning_sliderPressed();
    void on_fine_tunning_sliderReleased();

    /** ***********************************************************************************************
      * \fn void on_load_files_clicked() and nine more
      *  \brief Funciones responsables de cargar y ordenar los ficheros a analizar
      * ***********************************************************************************************
      */
    void on_load_files_clicked();
    void on_wavelet_file_cursorPositionChanged();
    void on_delete_file_clicked();
    void on_up_file_clicked();
    void on_down_file_clicked();

    void on_open_control_clicked();
    void on_control_file_cursorPositionChanged();
    void on_delete_control_clicked();
    void on_up_control_clicked();
    void on_down_control_clicked();

    /** ***********************************************************************************************
      * \fn void mouse_coordinates_ogl(int value)
      *  \brief Funciones responsables de recibir la señal de la ventana gráfica
      *         cuando el ratón mueve o pulsa para localizar página web con datos del ADN o hacer zoom
      * ***********************************************************************************************
      */
    void mouse_coordinates_ogl(int);
    void mouse_coordinates_ogl(int, int, int, int);

    /** ***********************************************************************************************
      * \fn void on_match_clicked()
      *  \brief Función responsable de filtrar los DMRs por intersección con genes conocidos
      * ***********************************************************************************************
      */
    void on_match_clicked();

    /** ***********************************************************************************************
      * \fn void on_dmr_dwt_level_sliderReleased()
      *  \brief Función responsable de modificar el nivel de transformación para identificación
      *         de DMRs y lanzar el nuevo análisis
      * ***********************************************************************************************
      */
    void on_dmr_dwt_level_sliderReleased();
    void on_dmr_dwt_level_valueChanged(int value);

    /** ***********************************************************************************************
      * \fn void on_save_dmr_list_clicked()
      *  \brief Función responsable de guardar los DMRs identificados con estadísticas de cada uno
      *         de los ficheros involucrados
      * ***********************************************************************************************
      */
    void on_save_dmr_list_clicked();

    /** ***********************************************************************************************
      * \fn void on_cobertura_sliderReleased()
      *  \brief Función responsable de modificar la cobertura mínima necesaria para formar
      *         la matriz de datos por muestra para realizar los cálculos
      * ***********************************************************************************************
      */
    void on_cobertura_sliderReleased();
    void on_cobertura_sliderMoved(int position);
    void on_min_coverage_textEdited(const QString &arg1);

    /** ***********************************************************************************************
      * \fn void on_mC_clicked() and three more
      *  \brief Funciones responsables actualizar el tipo de análisis a realizar
      * ***********************************************************************************************
      */
    void on_mC_clicked();
    void on_hmC_clicked();
    void on_forward_clicked();
    void on_reverse_clicked();

    /** ***********************************************************************************************
      * \fn void on_dmr_por_lote_clicked()
      *  \brief Función responsable abrir la ventana de identificación por lotes de DMRs
      * ***********************************************************************************************
      */
    void on_dmr_por_lote_clicked();

    /** ***********************************************************************************************
      * \fn void fichero_leido(int, int, int, int)
      *  \brief Función responsable de capturar los datos de los hilos de lectura
      *  \param sample   muestra que se ha leído
      *  \param chrom    cromosoma que se ha leído
      *  \param inicio   posición inicial de la muestra leída
      *  \param final    posición final de la muestra leída
      * ***********************************************************************************************
      */
    void fichero_leido(int, int, int, int);

    /** ***********************************************************************************************
      * \fn void cromosoma_leido(int)
      *  \brief Función responsable de controlar la elctura e identificación de DMRs por cromosoma
      *  \param chrom    cromosoma que se ha leído
      * ***********************************************************************************************
      */
    void cromosoma_leido(int);

    /** ***********************************************************************************************
      * \fn void refGen_worker_acabado(ulong)
      *  \brief Función responsable de recibir el número de genes definidos para un cromosoma dato
      *  \param ref_genes   número de genes
      * ***********************************************************************************************
      */
    void refGen_worker_acabado(ulong);

    /** ***********************************************************************************************
      * \fn void on_chrXX_clicked()
      *  \brief Función responsable de seleccionar el cromosoma a visualizar
      * ***********************************************************************************************
      */
    void on_chr01_clicked();
    void on_chr02_clicked();
    void on_chr03_clicked();
    void on_chr04_clicked();
    void on_chr05_clicked();
    void on_chr06_clicked();
    void on_chr07_clicked();
    void on_chr08_clicked();
    void on_chr09_clicked();
    void on_chr10_clicked();
    void on_chr11_clicked();
    void on_chr12_clicked();
    void on_chr13_clicked();
    void on_chr14_clicked();
    void on_chr15_clicked();
    void on_chr16_clicked();
    void on_chr17_clicked();
    void on_chr18_clicked();
    void on_chr19_clicked();
    void on_chr20_clicked();
    void on_chr21_clicked();
    void on_chr22_clicked();
    void on_chr23_clicked();
    void on_chr24_clicked();
    void on_chr25_clicked();

    /** ***********************************************************************************************
      * \fn void on_num_CpG_x_region_sliderReleased()
      *  \brief Función responsable de controlar el número mínimo de reads con cobertura en región
      *                 para valizar el DMR
      * ***********************************************************************************************
      */
    void on_num_CpG_x_region_sliderReleased();


private:
    Ui::HPG_Dhunter *ui;

    /** ***********************************************************************************************
      * \fn void on_input_file_textChanged(const QString &arg1)
      *  \brief Función responsable de abrir el explorador de ficheros, capturar el fichero de datos
      *         y actualizar la estrutura de variables de control de datos
      * ***********************************************************************************************
      */
//    void on_input_file_textChanged();


    /** ***********************************************************************************************
      *  \brief variables para control de posiciones extremas de fichero a analizar
      *  \param limite_inferior posición menor metilada en fichero de datos
      *  \param limite_superior posición mayor metilada en fichero de datos
      * ***********************************************************************************************
      */
    uint limite_inferior;
    uint limite_superior;


    /** ***********************************************************************************************
      *  \brief variables para búsqueda y muestra de DMRs
      *  \param threshold       limite inferior en cálculo de diferencias para DMRs
      *  \param **dmr_data      datos de transformación wavelet de todo el cromosoma en nivel 9
      *  \param **dmr_diff      datos de diferencias
      *  \param dmr_diff_rows   número de vectores de valores con los que busscar DMRs
      *  \param dmr_diff_cols   número de valores por vector con los que buscar DMRs por columna
      *  \param dmr_listo       señal para habilitar el tratamiento de coloreado en ventana de DMRs
      *  \param hallar_dmrs()   función de cálculo de diferencias
      *  \param *cursor         puntero a la línea en la ventana de DMRs para colorear y capturar su info
      *  \param num_genes       número de DMRs detectados
      *  \param dmrs            lista de todas las posiciones DMRs encontradas
      * ***********************************************************************************************
      */
    float       threshold;
    float       **dmr_data;
    float       *dmr_diff;
    uint        dmr_diff_cols;
    bool        dmr_listo;
    void        hallar_dmrs();
    QTextCursor *cursor;
    ulong       num_genes;
    QStringList dmrs;

    /** ***********************************************************************************************
      *  \brief variable para control de ancho de segmento a analizar
      *  \param paso_visualizacion  ancho ventana posiciones por nivel de visualización
      *                             útil para comenzar los cálculos de transformada en posición correcta
      *                             y se visualice la señal siempre igual para los mismos tramos
      *  \param ancho_ventana       valor del ancho de ventana en el momento se aactiva el ajuste fino
      * ***********************************************************************************************
      */
    uint paso_visualizacion;
    int  ancho_ventana;

    /** ***********************************************************************************************
      *  \brief variables para control de datos de cromosoma y hardware
      *  \param cromosoma           número de cromosoma a analizar
      *  \param memory_available    cantidad de memoria GPU disponible en el PC para controlar capacidad
      * ***********************************************************************************************
      */
    int cromosoma;
    int cromosoma_grid;
    int memory_available;

    /** ***********************************************************************************************
      *  \brief variables para control ventanas de visualización de ficheros a analizar
      *  \param wavelet_file                señal para habilitar el tratamiento de coloreado en ventana de DMRs
      *  \param color                       objeto para asignar color a fondo de línea
      *  \param color_char                  objeto para asignar color a texto
      *  \param *cursor_files               puntero a la línea en la ventana de ficheros para colorear y capturar su info
      *  \param primera_seleccion_casos     controla los ficheros a visualizar
      *  \param to_load                     controla si los ficheros seleccionados se han cargado en memoria
      *  \param directorio                  controla si se ha seleccionado el primer fichero para guardar path
      *  \param path                        guarda el último path del que se ha cargado un fichero
      *  \param wavelet_control_file        señal para habilitar el tratamiento de coloreado en ventana de DMRs
      *  \param *cursor_control_files       puntero a la línea en la ventana de ficheros para colorear y capturar su info
      *  \param primera_seleccion_control   controla los ficheros a visualizar
      *  \param visualiza_casos             guarda posiciones de muestras de casos a visualizar
      *  \param visualiza_control           guarda posiciones de muestras de control a visualizar
      *  \param h_haar_C_distribucion       guarda muestras calculadas en GPU
      * ***********************************************************************************************
      */
    bool             wavelet_file;
    QTextBlockFormat color;
    QTextCharFormat  color_char;
    QTextCursor      *cursor_files;
    bool             primera_seleccion_casos;
    bool             to_load;
    bool             directorio;
    QString          path;
    bool             wavelet_control_file;
    QTextCursor      *cursor_control_files;
    bool             primera_seleccion_control;
    vector<int>      visualiza_casos;
    vector<int>      visualiza_control;
    vector<int>      h_haar_C_distribucion;

    /** ***********************************************************************************************
      *  \brief variable para control de ancho de ventana visualizada
      *  \param page_pressed    bandera para actualizar ancho de ventana sobre la que actúa fine_tunning
      * ***********************************************************************************************
      */
    bool page_pressed;
    bool fine_tunning_pressed;

    /** ***********************************************************************************************
      *  \brief variable con los datos de metilación, cobertura y conteo de las muestras a analizar
      *  \param mc    datos de las muestras en matriz 3D (muestras->posicion->datos)
      * ***********************************************************************************************
      */
    vector<vector<vector<double>>> mc;                  // matriz con posiciones entre límites de todas las muestras
    vector<vector<vector<double>>> mc_aux;              // matriz ordenada de las muestras
    vector<vector<float>>          h_haar_C;            // matriz con los resultados wavelet de las muestras
    vector<vector<uint>>           posicion_metilada;   // acumula posiciones metiladas para cálculo DMR

    /** ***********************************************************************************************
      * \fn void dibuja()
      *  \brief función responsable de actualizar la visualización de los datos en GPU
      * ***********************************************************************************************
      */
    void dibuja();

    /** ***********************************************************************************************
      *  \brief variables para control de procesos en hilos
      *  \param hilo_files_worker   vector de hilos que albergan la función de lectura y procesamiento previo
      *  \param files_worker        vector de funciones de lectura y procesamiento previo de ficheros
      *  \param *hilo_refGen        hilo que alberga la función de lectura de genes por cromosoma
      *  \param *refgen_worker      función de lectura de genes por cromosoma
      * ***********************************************************************************************
      */
    QVector<QThread*>      hilo_files_worker;
    QVector<Files_worker*> files_worker;
    QThread               *hilo_refGen;
    RefGen                *refGen_worker;

    /** ***********************************************************************************************
      *  \brief variable de control de acceso a memoria compartida
      *  \param mutex   control de acceso a memoria compartida por los hilos
      * ***********************************************************************************************
      */
    QMutex mutex;

    /** ***********************************************************************************************
      *  \brief variable para control de evolución de la lectura de ficheros
      *  \param contador    aumenta su valor conforme avanza la lectura
      * ***********************************************************************************************
      */
    int contador;

    /** ***********************************************************************************************
      *  \brief variables para selección de parámetros con los que realizar el análisis
      *  \param _mc                 selecciona análisis por metilación
      *  \param _hmc                selecciona análisis por hidroximetilación
      *  \param _forward            selecciona análisis de ficheros forward
      *  \param _reverse            selecciona análisis de ficheros reverse
      * ***********************************************************************************************
      */
    bool    _mc;
    bool    _hmc;
    bool    _forward;
    bool    _reverse;
};

#endif // HPG_Dhunter_H
