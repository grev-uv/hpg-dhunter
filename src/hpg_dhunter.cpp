#include "hpg_dhunter.h"
#include "ui_hpg_dhunter.h"
#include <QFileDialog>
#include <QDebug>
#include <QFile>
#include <QTimer>
#include <QProcess>
#include <QRegularExpression>
#include <QMessageBox>
#include <QDesktopServices>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

HPG_Dhunter::HPG_Dhunter(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::HPG_Dhunter)
{
    ui->setupUi(this);

    // inicialización de variables ------------------------------------------------------------
    fichero                  = "";
    limite_inferior          = 0;
    limite_superior          = 1000;
    cuda_data.rango_inferior = size_t(ui->scroll_adn->value());
    cuda_data.rango_superior = size_t(ui->scroll_adn->value() + ui->scroll_adn->pageStep());
    cuda_data.mc_full        = nullptr;
    cuda_data.h_haar_C       = nullptr;
    cuda_data.d_aux          = nullptr;
    cuda_data.d_haar         = nullptr;
    cuda_data.refGen         = nullptr;
    cuda_data.d_glPtr        = nullptr;
    cuda_data.d_max          = new float[2];
    dmr_diff                 = nullptr;
    page_pressed             = true;
    fine_tunning_pressed     = true;
    _mc                      = true;
    _hmc                     = false;
    _forward                 = true;
    _reverse                 = false;
    cromosoma                = 0;
    cromosoma_grid           = 0;
    contador                 = 0;

    // inicialización de valores y posiciones de controles en ventana de aplicación -
    ui->limite_inferior->setText("0 - ");
    ui->limite_superior->setText(" - 1000");
    ui->ancho_ventana->setText((QString::number(ui->scroll_adn->pageStep())));
    ui->rango_inferior->setText(QString::number(ui->scroll_adn->value()));
    ui->rango_superior->setText(QString::number(ui->scroll_adn->value() + ui->scroll_adn->pageStep()));

    // inicialización de variables para cálculo de DMRs
    threshold = DMR_THRESHOLD;
    dmr_listo = false;
    cursor    = new QTextCursor();

    // cursor para ventanas con lista de ficheros a visualizar
    cursor_files              = new QTextCursor();
    wavelet_file              = false;
    primera_seleccion_casos   = true;
    cursor_control_files      = new QTextCursor();
    wavelet_control_file      = false;
    primera_seleccion_control = true;


    // comprueba la memoria disponible en la tarjeta gráfica para controlar los ficheros a cargar
    // ..captura la información suministrada por el comando "nvidia-smi"
    QProcess p;
    p.start("nvidia-smi");
    p.waitForFinished();
    QString data = p.readAllStandardOutput();
    p.close();

    // ..busca los datos que concuerdan con "[[:digit:]]+MiB" y se queda con el segundo dato
    //   que informa de la capacidad total de memoria de la tarjeta
    QRegularExpression re("(\\d+)MiB");
    QRegularExpressionMatchIterator i = re.globalMatch(data);
    i.next();
    QRegularExpressionMatch match = i.next();
    qDebug() << match.captured(0) << match.captured(1);

    // ..se asigna a la variable el valor en MiB
    memory_available = match.captured(1).toInt();

    ui->statusBar->showMessage("System available GPU RAM: " + QString::number(memory_available));


    // mantiene el último path en el explorador de ficheros
    directorio = false;
    to_load    = true;

    // conexiones con objetos
    connect(ui->ventana_opengl, SIGNAL(ogl_coordinates(int)), this, SLOT(mouse_coordinates_ogl(int)));
    connect(ui->ventana_opengl, SIGNAL(ogl_coordinates(int, int, int, int)), this, SLOT(mouse_coordinates_ogl(int, int, int, int)));
}

// ************************************************************************************************
HPG_Dhunter::~HPG_Dhunter()
{
    ui->ventana_opengl->unregisterBuffer();
    cuda_end(cuda_data);
    delete [] cuda_data.mc_full[0];
    delete [] cuda_data.mc_full;
    delete ui;
}

// ************************************************************************************************
// **************************************ZONA CARGA************************************************
// ************************************************************************************************
// -------------VENTANA CASOS-----------------
void HPG_Dhunter::on_open_file_clicked()
{
    // abre ventana de explorador de ficheros para seleccionar el archivo de datos
    fichero = QFileDialog::getExistingDirectory(this,
                                                tr("Select a case sample folder"),
                                                (directorio) ? path : QDir::homePath(),
                                                QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks
                                               );

    if(fichero.isEmpty() || fichero.isNull())
        fichero = "";
    else
    {
        wavelet_file = true;
        ui->wavelet_file->appendPlainText(fichero);
        QStringList lista = ui->wavelet_file->toPlainText().split("\n");
        if (lista.size() + ui->control_file->toPlainText().split("\n")[0].size() > 1)
            ui->load_files->setEnabled(true);
        ui->analiza->setEnabled(false);
        ui->dmrs->setEnabled(false);
        ui->num_CpG_x_region->setEnabled(false);
        ui->threshold->setEnabled(false);
        ui->dmr_dwt_level->setEnabled(false);
        ui->save_dmr_list->setEnabled(false);
        ui->delete_file->setEnabled(true);
        ui->up_file->setEnabled(true);
        ui->down_file->setEnabled(true);

        directorio = true;
        to_load    = true;
        primera_seleccion_casos = true;
        fichero.chop(fichero.split("/").back().size());
        path = fichero;

        // colorear las lineas
        cursor_files->movePosition(QTextCursor::Start);
        for (int i = 0; i < lista.size(); i++)
        {
            color.setBackground(Qt::white);
            color_char.setForeground(QColor(255, (40 * i <= 255)? 40 * i : 255, 0));
            cursor_files->movePosition(QTextCursor::StartOfBlock);
            cursor_files->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
            cursor_files->select(QTextCursor::BlockUnderCursor);
            cursor_files->setCharFormat(color_char);
            cursor_files->setBlockFormat(color);
            cursor_files->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);
        }
    }
}

//*************************************************************************************************
void HPG_Dhunter::on_wavelet_file_cursorPositionChanged()
{
    if (wavelet_file)
    {
        // para marcar fila seleccionada en color gris de fondo
        if (to_load)
        {
            // desmarca la línea previa quitando color de fondo
            color.setBackground(Qt::white);
            cursor_files->select(QTextCursor::LineUnderCursor);
            cursor_files->setBlockFormat(color);

            // adquiere el cursor de la línea seleccionada con el ratón
            *cursor_files = ui->wavelet_file->textCursor();
            cursor_files->movePosition(QTextCursor::StartOfBlock);
            cursor_files->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);

            // marca con color de fondo la línea seleccionada
            color.setBackground(Qt::gray);
            cursor_files->select(QTextCursor::LineUnderCursor);
            cursor_files->setBlockFormat(color);
        }
        // para rellenar la lista de ficheros con colores diferentes por fila
        else
        {
            if (primera_seleccion_casos)
            {
                visualiza_casos.assign(uint(ficheros_case.size()), 0);
                primera_seleccion_casos = false;
            }

            qDebug() << visualiza_casos;

            *cursor_files = ui->wavelet_file->textCursor();
            uint muestra = uint(cursor_files->blockNumber());
            visualiza_casos[muestra] = !visualiza_casos[muestra];

            qDebug() << visualiza_casos;

            int fondo = 0;
            cursor_files->movePosition(QTextCursor::Start);
            color_char.setForeground(Qt::black);
            for (uint i = 0; i < visualiza_casos.size(); i++)
            {
                if (visualiza_casos[i])
                {
                    color.setBackground(QColor(255, (40 * fondo <= 255)? 40 * fondo : 255, 0));
                    fondo == 0 ? fondo = 2 : fondo++;
                }
                else
                    color.setBackground(Qt::white);

                cursor_files->movePosition(QTextCursor::StartOfBlock);
                cursor_files->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
                cursor_files->select(QTextCursor::LineUnderCursor);
                cursor_files->setBlockFormat(color);
                cursor_files->select(QTextCursor::BlockUnderCursor);
                cursor_files->setCharFormat(color_char);
                cursor_files->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);
            }

            ui->analiza->setEnabled(true);
        }
    }
}

//*************************************************************************************************
void HPG_Dhunter::on_delete_file_clicked()
{
    // variables internas
    int line_number  = cursor_files->blockNumber();                 // número de línea actual del cursor
    QStringList list = ui->wavelet_file->toPlainText().split("\n"); // lista de strings con ficheros

    // deshabilita el cambio de color por acción sobre el cursor
    wavelet_file = false;
    // borra el fichero seleccionado
    list.removeAt(line_number);
    // limpia la lista de visualización
    ui->wavelet_file->clear();
    // si quedan ficheros en la lista, la copia en la lista de visualización
    if (!list.isEmpty())
    {
        wavelet_file = true;
        ui->wavelet_file->appendPlainText(list.join('\n'));
    }
    else
    {
        // si no hay ficheros deshabilita el botón de ejecutar alineamiento
        ui->load_files->setEnabled(false);
        ui->analiza->setEnabled(false);
        ui->dmrs->setEnabled(false);
        ui->save_dmr_list->setEnabled(false);
        ui->delete_file->setEnabled(false);
        ui->up_file->setEnabled(false);
        ui->down_file->setEnabled(false);
    }

    // aplica fondo blanco a línea bajo el cursor
    color.setBackground(Qt::white);
    cursor_files->select((QTextCursor::LineUnderCursor));
    cursor_files->setBlockFormat(color);

    // colorear las lineas
    QStringList lista = ui->wavelet_file->toPlainText().split("\n");
    cursor_files->movePosition(QTextCursor::Start);
    for (int i = 0; i < lista.size(); i++)
    {
        color.setBackground(Qt::white);
        color_char.setForeground(QColor(255, (40 * i <= 255)? 40 * i : 255, 0));
        cursor_files->movePosition(QTextCursor::StartOfBlock);
        cursor_files->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
        cursor_files->select(QTextCursor::BlockUnderCursor);
        cursor_files->setCharFormat(color_char);
        cursor_files->setBlockFormat(color);
        cursor_files->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);
    }

    // mueve el cursor al principo de la lista de visualización
    cursor_files->movePosition(QTextCursor::Start);
    // traslada el cursor a la posición que tenía antes de borrar
    for (int i = 0; i < line_number; i++)
    {
        cursor_files->movePosition(QTextCursor::StartOfBlock);
        cursor_files->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
        cursor_files->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);

        if (i == list.size() - 2)
            break;
    }

    // actualiza la lista con la posición del cursor adecuada para que se resalte
    ui->wavelet_file->setTextCursor(*cursor_files);
}

//*************************************************************************************************
void HPG_Dhunter::on_up_file_clicked()
{
    if (cursor_files->blockNumber() > 0)
    {
        // variables internas
        int line_number  = cursor_files->blockNumber();                 // número de línea del cursor
        QStringList list = ui->wavelet_file->toPlainText().split("\n"); // lista de ficheros
        QString line     = list.at(line_number);                        // fichero en línea seleccionada

        // deshabilita el cambio de color por acción sobre el cursor
        wavelet_file = false;
        // borra el fichero seleccionado
        list.removeAt(line_number);
        // inserta el fichero selecionado en una posición más arriba
        list.insert(line_number - 1, line);
        // limpia la lista de visualización
        ui->wavelet_file->clear();
        // si quedan ficheros en la lista, la copia en la lista de visualización
        if (!list.isEmpty())
        {
            wavelet_file = true;
            ui->wavelet_file->appendPlainText(list.join('\n'));
        }

        // aplica fondo blanco a línea bajo el cursor
        color.setBackground(Qt::white);
        cursor_files->select((QTextCursor::LineUnderCursor));
        cursor_files->setBlockFormat(color);

        // colorear las lineas
        QStringList lista = ui->wavelet_file->toPlainText().split("\n");
        cursor_files->movePosition(QTextCursor::Start);
        for (int i = 0; i < lista.size(); i++)
        {
            qDebug() << "pintando " << i << " " << int(i/2);
            color.setBackground(Qt::white);
            color_char.setForeground(QColor(255, (40 * i <= 255)? 40 * i : 255, 0));
            cursor_files->movePosition(QTextCursor::StartOfBlock);
            cursor_files->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
            cursor_files->select(QTextCursor::BlockUnderCursor);
            cursor_files->setCharFormat(color_char);
            cursor_files->setBlockFormat(color);
            cursor_files->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);
        }

        // mueve el cursor al principio de la lista de visualización
        cursor_files->movePosition(QTextCursor::Start);
        // traslada el cursor a la posición que ocupa el fichero trasladado
        for (int i = 0; i < line_number - 1; i++)
        {
            cursor_files->movePosition(QTextCursor::StartOfBlock);
            cursor_files->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
            cursor_files->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);
        }

        // actualiza la lista con la posición del cursor adecuada para que se resalte
        ui->wavelet_file->setTextCursor(*cursor_files);
    }
}

//*************************************************************************************************
void HPG_Dhunter::on_down_file_clicked()
{
    // variables internas
    QStringList list = ui->wavelet_file->toPlainText().split("\n");   // lista de ficheros

    // proceder con el desplazamiento abajo si no es el último elemento de la lista
    if (cursor_files->blockNumber() < list.size())
    {
        // variables internas
        int line_number = cursor_files->blockNumber();  // número de línea del cursor
        QString line    = list.at(line_number);         // fichero en línea seleccionada

        // deshabilita el cambio de color por acción sobre el cursor
        wavelet_file = false;
        // borra el fichero seleccionado
        list.removeAt(line_number);
        // inserta el fichero seleccionado en una posición más abajo
        list.insert(line_number + 1, line);
        // limpia la lista de visualización
        ui->wavelet_file->clear();
        // se quedan ficheros en la lista, la copia en la lista de visualización
        if (!list.isEmpty())
        {
            wavelet_file = true;
            ui->wavelet_file->appendPlainText(list.join('\n'));
        }

        // aplica fondo blanco a línea bajo el cursor
        color.setBackground(Qt::white);
        cursor_files->select((QTextCursor::LineUnderCursor));
        cursor_files->setBlockFormat(color);

        // colorear las lineas
        QStringList lista = ui->wavelet_file->toPlainText().split("\n");
        cursor_files->movePosition(QTextCursor::Start);
        for (int i = 0; i < lista.size(); i++)
        {
            qDebug() << "pintando " << i << " " << int(i/2);
            color.setBackground(Qt::white);
            color_char.setForeground(QColor(255, (40 * i <= 255)? 40 * i : 255, 0));
            cursor_files->movePosition(QTextCursor::StartOfBlock);
            cursor_files->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
            cursor_files->select(QTextCursor::BlockUnderCursor);
            cursor_files->setCharFormat(color_char);
            cursor_files->setBlockFormat(color);
            cursor_files->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);
        }

        // mueve el cursor al principio de la lista de visualización
        cursor_files->movePosition(QTextCursor::Start);
        // traslada el cursor a la posición que ocupa el fichero trasladado
        for (int i = 0; i <= line_number; i++)
        {
            cursor_files->movePosition(QTextCursor::StartOfBlock);
            cursor_files->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
            cursor_files->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);

            // comprueba que no ha alcanzado el final de la lista
            if (i == list.size() - 2)
                break;
        }

        // actualiza la lista con la posición del cursor adecuada para que se resalte
        ui->wavelet_file->setTextCursor(*cursor_files);
    }
}


// -------------VENTANA CONTROL-----------------
// ************************************************************************************************
void HPG_Dhunter::on_open_control_clicked()
{
    // abre ventana de explorador de ficheros para seleccionar el archivo de datos
    fichero = QFileDialog::getExistingDirectory(this,
                                                tr("Select a control sample folder"),
                                                (directorio) ? path : QDir::homePath(),
                                                QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks
                                               );

    if(fichero.isEmpty() || fichero.isNull())
        fichero = "";
    else
    {
        wavelet_control_file = true;
        ui->control_file->appendPlainText(fichero);
        QStringList lista = ui->control_file->toPlainText().split("\n");
        if (lista.size() + ui->wavelet_file->toPlainText().split("\n")[0].size() > 1)
            ui->load_files->setEnabled(true);
        ui->analiza->setEnabled(false);
        ui->dmrs->setEnabled(false);
        ui->num_CpG_x_region->setEnabled(false);
        ui->threshold->setEnabled(false);
        ui->dmr_dwt_level->setEnabled(false);
        ui->save_dmr_list->setEnabled(false);
        ui->delete_control->setEnabled(true);
        ui->up_control->setEnabled(true);
        ui->down_control->setEnabled(true);

        directorio = true;
        to_load    = true;
        primera_seleccion_control = true;
        fichero.chop(fichero.split("/").back().size());
        path = fichero;

        // colorear las lineas en dos grupos de color alternativamente
        cursor_control_files->movePosition(QTextCursor::Start);
        for (int i = 0; i < lista.size(); i++)
        {
            color.setBackground(Qt::white);
            color_char.setForeground(QColor(0, (40 * i <= 255)? 40 * i : 255, 255));
            cursor_control_files->movePosition(QTextCursor::StartOfBlock);
            cursor_control_files->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
            cursor_control_files->select(QTextCursor::BlockUnderCursor);
            cursor_control_files->setCharFormat(color_char);
            cursor_control_files->setBlockFormat(color);
            cursor_control_files->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);
        }
    }
}

//*************************************************************************************************
void HPG_Dhunter::on_control_file_cursorPositionChanged()
{
    if (wavelet_control_file)
    {
        if (to_load)
        {
        // desmarca la línea previa quitando color de fondo
        color.setBackground(Qt::white);
        cursor_control_files->select(QTextCursor::LineUnderCursor);
        cursor_control_files->setBlockFormat(color);

        // adquiere el cursor de la línea seleccionada con el ratón
        *cursor_control_files = ui->control_file->textCursor();
        cursor_control_files->movePosition(QTextCursor::StartOfBlock);
        cursor_control_files->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);

        // marca con color de fondo la línea seleccionada
        color.setBackground(Qt::gray);
        cursor_control_files->select(QTextCursor::LineUnderCursor);
        cursor_control_files->setBlockFormat(color);
        }
        else
        {
            if (primera_seleccion_control)
            {
                visualiza_control.assign(uint(ficheros_control.size()), 0);
                primera_seleccion_control = false;
            }

            qDebug() << visualiza_control;

            *cursor_control_files = ui->control_file->textCursor();
            uint muestra = uint(cursor_control_files->blockNumber());
            visualiza_control[muestra] = !visualiza_control[muestra];

            qDebug() << visualiza_control;

            int fondo = 0;
            cursor_control_files->movePosition(QTextCursor::Start);
            color_char.setForeground(Qt::black);
            for (uint i = 0; i < visualiza_control.size(); i++)
            {
                if (visualiza_control[i])
                {
                    color.setBackground(QColor(0, (40 * fondo <= 255)? 40 * fondo : 255, 255));
                    fondo == 0 ? fondo = 2 : fondo++;
                }
                else
                    color.setBackground(Qt::white);

                cursor_control_files->movePosition(QTextCursor::StartOfBlock);
                cursor_control_files->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
                cursor_control_files->select(QTextCursor::LineUnderCursor);
                cursor_control_files->setBlockFormat(color);
                cursor_control_files->select(QTextCursor::BlockUnderCursor);
                cursor_control_files->setCharFormat(color_char);
                cursor_control_files->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);
            }

            ui->analiza->setEnabled(true);
        }
    }
}

//*************************************************************************************************
void HPG_Dhunter::on_delete_control_clicked()
{
    // variables internas
    int line_number  = cursor_control_files->blockNumber();         // número de línea actual del cursor
    QStringList list = ui->control_file->toPlainText().split("\n"); // lista de strings con ficheros

    // deshabilita el cambio de color por acción sobre el cursor
    wavelet_control_file = false;
    // borra el fichero seleccionado
    list.removeAt(line_number);
    // limpia la lista de visualización
    ui->control_file->clear();
    // si quedan ficheros en la lista, la copia en la lista de visualización
    if (!list.isEmpty())
    {
        wavelet_control_file = true;
        ui->control_file->appendPlainText(list.join('\n'));
    }
    else
    {
        // si no hay ficheros deshabilita el botón de ejecutar alineamiento
        ui->load_files->setEnabled(false);
        ui->analiza->setEnabled(false);
        ui->dmrs->setEnabled(false);
        ui->save_dmr_list->setEnabled(false);
        ui->delete_control->setEnabled(false);
        ui->up_control->setEnabled(false);
        ui->down_control->setEnabled(false);
    }

    // aplica fondo blanco a línea bajo el cursor
    color.setBackground(Qt::white);
    cursor_control_files->select((QTextCursor::LineUnderCursor));
    cursor_control_files->setBlockFormat(color);

    // colorear las lineas
    QStringList lista = ui->control_file->toPlainText().split("\n");
    cursor_control_files->movePosition(QTextCursor::Start);
    for (int i = 0; i < lista.size(); i++)
    {
        color.setBackground(Qt::white);
        color_char.setForeground(QColor(0, (40 * i <= 255)? 40 * i : 255, 255));
        cursor_control_files->movePosition(QTextCursor::StartOfBlock);
        cursor_control_files->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
        cursor_control_files->select(QTextCursor::BlockUnderCursor);
        cursor_control_files->setCharFormat(color_char);
        cursor_control_files->setBlockFormat(color);
        cursor_control_files->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);
    }

    // mueve el cursor al principo de la lista de visualización
    cursor_control_files->movePosition(QTextCursor::Start);
    // traslada el cursor a la posición que tenía antes de borrar
    for (int i = 0; i < line_number; i++)
    {
        cursor_control_files->movePosition(QTextCursor::StartOfBlock);
        cursor_control_files->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
        cursor_control_files->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);

        if (i == list.size() - 2)
            break;
    }

    // actualiza la lista con la posición del cursor adecuada para que se resalte
    ui->control_file->setTextCursor(*cursor_control_files);
}

//*************************************************************************************************
void HPG_Dhunter::on_up_control_clicked()
{
    if (cursor_control_files->blockNumber() > 0)
    {
        // variables internas
        int line_number  = cursor_control_files->blockNumber();         // número de línea del cursor
        QStringList list = ui->control_file->toPlainText().split("\n"); // lista de ficheros
        QString line     = list.at(line_number);                        // fichero en línea seleccionada

        // deshabilita el cambio de color por acción sobre el cursor
        wavelet_control_file = false;
        // borra el fichero seleccionado
        list.removeAt(line_number);
        // inserta el fichero selecionado en una posición más arriba
        list.insert(line_number - 1, line);
        // limpia la lista de visualización
        ui->control_file->clear();
        // si quedan ficheros en la lista, la copia en la lista de visualización
        if (!list.isEmpty())
        {
            wavelet_control_file = true;
            ui->control_file->appendPlainText(list.join('\n'));
        }

        // aplica fondo blanco a línea bajo el cursor
        color.setBackground(Qt::white);
        cursor_control_files->select((QTextCursor::LineUnderCursor));
        cursor_control_files->setBlockFormat(color);

        // colorear las lineas
        QStringList lista = ui->control_file->toPlainText().split("\n");
        cursor_control_files->movePosition(QTextCursor::Start);
        for (int i = 0; i < lista.size(); i++)
        {
            qDebug() << "pintando " << i << " " << int(i/2);
            color.setBackground(Qt::white);
            color_char.setForeground(QColor(0, (40 * i <= 255)? 40 * i : 255, 255));
            cursor_control_files->movePosition(QTextCursor::StartOfBlock);
            cursor_control_files->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
            cursor_control_files->select(QTextCursor::BlockUnderCursor);
            cursor_control_files->setCharFormat(color_char);
            cursor_control_files->setBlockFormat(color);
            cursor_control_files->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);
        }

        // mueve el cursor al principio de la lista de visualización
        cursor_control_files->movePosition(QTextCursor::Start);
        // traslada el cursor a la posición que ocupa el fichero trasladado
        for (int i = 0; i < line_number - 1; i++)
        {
            cursor_control_files->movePosition(QTextCursor::StartOfBlock);
            cursor_control_files->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
            cursor_control_files->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);
        }

        // actualiza la lista con la posición del cursor adecuada para que se resalte
        ui->control_file->setTextCursor(*cursor_control_files);
    }
}

//*************************************************************************************************
void HPG_Dhunter::on_down_control_clicked()
{
    // variables internas
    QStringList list = ui->control_file->toPlainText().split("\n");   // lista de ficheros

    // proceder con el desplazamiento abajo si no es el último elemento de la lista
    if (cursor_control_files->blockNumber() < list.size())
    {
        // variables internas
        int line_number = cursor_control_files->blockNumber();  // número de línea del cursor
        QString line    = list.at(line_number);         // fichero en línea seleccionada

        // deshabilita el cambio de color por acción sobre el cursor
        wavelet_control_file = false;
        // borra el fichero seleccionado
        list.removeAt(line_number);
        // inserta el fichero seleccionado en una posición más abajo
        list.insert(line_number + 1, line);
        // limpia la lista de visualización
        ui->control_file->clear();
        // se quedan ficheros en la lista, la copia en la lista de visualización
        if (!list.isEmpty())
        {
            wavelet_control_file = true;
            ui->control_file->appendPlainText(list.join('\n'));
        }

        // aplica fondo blanco a línea bajo el cursor
        color.setBackground(Qt::white);
        cursor_control_files->select((QTextCursor::LineUnderCursor));
        cursor_control_files->setBlockFormat(color);

        // colorear las lineas
        QStringList lista = ui->control_file->toPlainText().split("\n");
        cursor_control_files->movePosition(QTextCursor::Start);
        for (int i = 0; i < lista.size(); i++)
        {
            qDebug() << "pintando " << i << " " << int(i/2);
            color.setBackground(Qt::white);
            color_char.setForeground(QColor(0, (40 * i <= 255)? 40 * i : 255, 255));
            cursor_control_files->movePosition(QTextCursor::StartOfBlock);
            cursor_control_files->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
            cursor_control_files->select(QTextCursor::BlockUnderCursor);
            cursor_control_files->setCharFormat(color_char);
            cursor_control_files->setBlockFormat(color);
            cursor_control_files->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);
        }

        // mueve el cursor al principio de la lista de visualización
        cursor_control_files->movePosition(QTextCursor::Start);
        // traslada el cursor a la posición que ocupa el fichero trasladado
        for (int i = 0; i <= line_number; i++)
        {
            cursor_control_files->movePosition(QTextCursor::StartOfBlock);
            cursor_control_files->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
            cursor_control_files->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);

            // comprueba que no ha alcanzado el final de la lista
            if (i == list.size() - 2)
                break;
        }

        // actualiza la lista con la posición del cursor adecuada para que se resalte
        ui->control_file->setTextCursor(*cursor_control_files);
    }
}

// ************************************************************************************************
// **************************************ZONA CONTROLES INTERFAZ***********************************
// ************************************************************************************************
void HPG_Dhunter::on_scroll_adn_valueChanged(int value)
{
    // actualizar los valores de los límites de la ventana de datos a analizar
    ui->rango_inferior->setText(QString::number(value));
    ui->rango_superior->setText(QString::number(value + ui->scroll_adn->pageStep()));

    // actualizar los valores de los límites del rango de datos a analizar en GPU
    cuda_data.rango_inferior = size_t(value) - size_t(limite_inferior);
    cuda_data.rango_superior = size_t(value) - size_t(limite_inferior) + size_t(ui->scroll_adn->pageStep());
}

// ************************************************************************************************
void HPG_Dhunter::on_page_valueChanged(int position)
{
    // actualiza los valores del control de scroll_adn
    // tanto el ancho de ventana como su valor inicial
    uint nuevo_maximo = limite_superior - uint(position);
    uint nuevo_valor = uint(ui->scroll_adn->value() + ((ui->scroll_adn->pageStep() - position) * 0.5));

    ui->scroll_adn->setPageStep(position);
    ui->scroll_adn->setMaximum(int(nuevo_maximo));
    ui->scroll_adn->setValue(int(nuevo_valor));

    // actualiza los valores en ventana de visualización
    ui->ancho_ventana->setText((QString::number(ui->scroll_adn->pageStep())));
    ui->rango_inferior->setText(QString::number(ui->scroll_adn->value()));
    ui->rango_superior->setText(QString::number(ui->scroll_adn->value() + ui->scroll_adn->pageStep()));

    // actualiza en número máximo de niveles que puede transformar
    ui->slider_nivel->setMaximum(int(log2(ui->scroll_adn->pageStep())));
    if (ui->slider_nivel->maximum() <= ui->slider_nivel->value())
        ui->slider_nivel->setValue(ui->slider_nivel->maximum() - 1);;

    // actualizar los valores de los límites del rango de datos a analizar en GPU
    // dentro de la estructura de variables de control cuda
    if (nuevo_valor > limite_inferior)
        cuda_data.rango_inferior = nuevo_valor - limite_inferior;
    else
        cuda_data.rango_inferior = 0;
    cuda_data.rango_superior = nuevo_valor - limite_inferior + uint(ui->scroll_adn->pageStep());
}

// ************************************************************************************************
void HPG_Dhunter::on_rango_superior_editingFinished()
{
    uint arg1 = ui->rango_superior->text().toUInt();
    uint inicio = uint(ui->scroll_adn->value());

    if (arg1 < limite_superior && arg1 > limite_inferior)
    {
        cuda_data.rango_superior = arg1 - limite_inferior;
        if (arg1 > inicio)
        {
            ui->page->setValue(int(cuda_data.rango_superior - cuda_data.rango_inferior));
            ui->scroll_adn->setValue(int(inicio));
        }
        else
        {
            if (int(arg1 - uint(ui->page->value())) < int(limite_inferior))
            {
                cuda_data.rango_inferior = 0;
                ui->page->setValue(int(cuda_data.rango_superior));
                ui->scroll_adn->setValue(int(limite_inferior));
            }
            else
            {
                ui->scroll_adn->setValue(int(arg1) - ui->page->value());
                cuda_data.rango_inferior = uint(ui->scroll_adn->value()) - limite_inferior;
            }
        }
    }

    // lanza la transformación
    if (ui->wavelet_file->blockCount() != 1)
        dibuja();
}

// ************************************************************************************************
void HPG_Dhunter::on_rango_inferior_editingFinished()
{
    uint arg1 = ui->rango_inferior->text().toUInt();

    if (arg1 < limite_superior &&  arg1 > limite_inferior)
    {
        cuda_data.rango_inferior = arg1 - limite_inferior;

        if (arg1 < uint(ui->scroll_adn->value() + ui->scroll_adn->pageStep()))
        {
            ui->page->setValue(int(cuda_data.rango_superior - cuda_data.rango_inferior));
            ui->scroll_adn->setValue(int(arg1));
        }
        else
        {
            if (arg1 + uint(ui->scroll_adn->pageStep()) > limite_superior)
            {
                cuda_data.rango_superior = limite_superior - limite_inferior;
                ui->page->setValue(int(cuda_data.rango_superior - cuda_data.rango_inferior));
                ui->scroll_adn->setValue(int(arg1));
            }
            else
            {
                ui->scroll_adn->setValue(int(arg1));
                cuda_data.rango_superior = cuda_data.rango_inferior + uint(ui->scroll_adn->pageStep());
            }
        }
    }

    // lanza la transformación
    if (ui->wavelet_file->blockCount() != 1)
        dibuja();
}

// ************************************************************************************************
void HPG_Dhunter::on_slider_nivel_sliderReleased()
{
    // actualiza valor en la estructura de variables
    cuda_data.levels = ui->slider_nivel->value();

    // lanza la transformación
    if (ui->wavelet_file->blockCount() != 1)
        dibuja();
}

// ************************************************************************************************
void HPG_Dhunter::on_scroll_adn_sliderReleased()
{
    // lanza la transformación
    if (ui->wavelet_file->blockCount() != 1)
        dibuja();
}

// ************************************************************************************************
void HPG_Dhunter::on_page_sliderReleased()
{
    // inicializa lider de ajuste fino
    ui->fine_tunning->setValue(0);
    if (ui->ancho_ventana->text().toInt() < 100000)
        ui->fine_tunning->setMinimum(-1 * ui->ancho_ventana->text().toInt());
    else
        ui->fine_tunning->setMinimum(-100000);
    ui->fine_tunning->setValue(0);
    page_pressed = true;

    // lanza la transformación
    if (ui->wavelet_file->blockCount() != 1)
        dibuja();
}


// ************************************************************************************************
void HPG_Dhunter::on_fine_tunning_valueChanged(int value)
{
    if (fine_tunning_pressed)
    {
        if (ancho_ventana + value < 100)
            ui->ancho_ventana->setText("100");
        else
            ui->ancho_ventana->setText(QString::number(ancho_ventana + value));

        // actualiza en número máximo de niveles que puede transformar
        ui->slider_nivel->setMaximum(int(log2(ui->scroll_adn->pageStep())));
        if (ui->slider_nivel->maximum() <= ui->slider_nivel->value())
            ui->slider_nivel->setValue(ui->slider_nivel->maximum() - 1);
    }
}

// ************************************************************************************************
void HPG_Dhunter::on_fine_tunning_sliderPressed()
{
    fine_tunning_pressed = true;

    if (page_pressed)
    {
        ancho_ventana = ui->ancho_ventana->text().toInt();
        page_pressed = false;
    }
}

// ************************************************************************************************
void HPG_Dhunter::on_fine_tunning_sliderReleased()
{
    ui->page->setValue(ui->ancho_ventana->text().toInt());
    fine_tunning_pressed = false;

    // lanza la transformación si hay ficheros para mostrar
    if (ui->wavelet_file->blockCount() != 1)
        dibuja();
}

// ************************************************************************************************
void HPG_Dhunter::on_threshold_valueChanged(int value)
{
    threshold = float(value * 0.01);
    ui->threshold_label->setText(QString::number(double(threshold), 'f', 2));
}

// ************************************************************************************************
void HPG_Dhunter::on_threshold_sliderReleased()
{
    //desactiva la lectura desde la lista de dmrs
    dmr_listo = false;

    // lanza el cálculo
    if (ui->dmr_position->blockCount() != 1)
        hallar_dmrs();
}

// ************************************************************************************************
void HPG_Dhunter::on_num_CpG_x_region_sliderReleased()
{
    //desactiva la lectura desde la lista de dmrs
    dmr_listo = false;

    // lanza el cálculo
    if (ui->dmr_position->blockCount() != 1)
        on_dmrs_clicked();
}

// ************************************************************************************************
void HPG_Dhunter::on_slider_nivel_valueChanged(int value)
{
    ui->dwt_level->setText(QString::number(uint(pow(2,value))));
}


// ************************************************************************************************
void HPG_Dhunter::on_dmr_dwt_level_valueChanged(int value)
{
    ui->label_8->setText(QString::number(uint(pow(2, value))));
}

// ************************************************************************************************
void HPG_Dhunter::on_dmr_dwt_level_sliderReleased()
{
    ui->dmrs->setEnabled(true);
    //ui->save_dmr_list->setEnabled(true);
}

// ************************************************************************************************
void HPG_Dhunter::on_cobertura_sliderReleased()
{
    ui->analiza->setEnabled(true);
}

// ************************************************************************************************
void HPG_Dhunter::on_mC_clicked()
{
    if (ui->wavelet_file->blockCount() > 1)
        ui->analiza->setEnabled(true);

    _mc  = ui->mC->isChecked();
    _hmc = ui->hmC->isChecked();
}

// ************************************************************************************************
void HPG_Dhunter::on_hmC_clicked()
{
    if (ui->wavelet_file->blockCount() > 1)
        ui->analiza->setEnabled(true);

    _mc  = ui->mC->isChecked();
    _hmc = ui->hmC->isChecked();
}

// ************************************************************************************************
void HPG_Dhunter::on_dmr_por_lote_clicked()
{
    this->hide();

    // lanza la detección de DMRs
    QProcess p;
    int dato = p.execute("/home/lifercor/programas/Bio-Informatica/qt/HPG_Dhunter/HPG_Dhunter_v2-6_FDD/HPG_Dhunter/hpg_dhunter");
    p.waitForFinished(-1); // evitamos problemas de terminación más allá de los 30 segundos
    p.close();

    this->show();
    qDebug() << "info sobre la ejecución de la visualización: " << dato;
}

// ************************************************************************************************
void HPG_Dhunter::on_forward_clicked()
{
    _forward = ui->forward->isChecked();

    if (!_forward && !_reverse)
    {
        ui->reverse->setChecked(true);
        _reverse = true;
    }

    ui->load_files->setEnabled(true);
}

// ************************************************************************************************
void HPG_Dhunter::on_reverse_clicked()
{
    _reverse = ui->reverse->isChecked();

    if (!_forward && !_reverse)
    {
        ui->forward->setChecked(true);
        _forward = true;
    }

    ui->load_files->setEnabled(true);
}

// ************************************************************************************************
void HPG_Dhunter::mouse_coordinates_ogl(int x)
{
    QString linea = "";
    uint factor_escala = uint(ui->ancho_ventana->text().toFloat() / ui->ventana_opengl->width());
    int cobertura_maxima = 0;

    // posición del puntero en la ventana de visualización
    uint posicion = ui->rango_inferior->text().toUInt() +
                    uint(x * (ui->ancho_ventana->text().toFloat()
                            / ui->ventana_opengl->width()));
    linea.append(QString::number(posicion) + " coverage: ");

    // cobertura por muestra en cada posición

    for (uint i = 0; i < mc.size(); i++)
    {
        // ..búsqueda binaria sobre el fichero de cobertura de la posición
        uint inicio = 0;
        uint fin    = uint(mc[i].size());
        uint mitad  = 0;
        bool match  = false;

        while (inicio < fin && !match)
        {
            mitad  = uint((inicio + fin) * 0.5);

            if (uint(mc[i][mitad][0]) == posicion)
                break;
            else
            {
                if (uint(mc[i][mitad][0]) < posicion)
                    inicio = mitad + 1;
                else
                    fin = mitad - 1;
            }
        }

        // una vez posicionado, se busca el valor máximo de cobertura en la franja de posiciones
        // del genoma visualizado que comprende el puntero en cada posición de la ventana
        posicion += factor_escala;
        while (posicion > mc[i][mitad][0])
        {
            if (cobertura_maxima < mc[i][mitad][2])
                cobertura_maxima = int(mc[i][mitad][2]);
            mitad++;
        }

        linea.append(" s" + QString::number(i+1) + " " + QString::number(cobertura_maxima));
    }

    if (ui->wavelet_file->blockCount() != 1 && !to_load && !ui->analiza->isEnabled())
        ui->mouse_xpos->setText(linea);
}

// ************************************************************************************************
void HPG_Dhunter::mouse_coordinates_ogl(int x, int y, int xR, int yR)
{
    qDebug() << "coordenadas press x, y: " << x << ", " << y << "  -  coordenadas release xR, yR: " << xR << ", " << yR;

    int inicio = ui->rango_inferior->text().toInt() +
                 int(x * (ui->ancho_ventana->text().toFloat()
                       / ui->ventana_opengl->width()));
    int fin    = ui->rango_inferior->text().toInt() +
                 int(xR * (ui->ancho_ventana->text().toFloat()
                       / ui->ventana_opengl->width()));

    // si la distancia en el eje x no es grande se estima pulsación para abrir explrador genético en web
    if (xR - x < 5)
    {
        // componer mensaje con número de cromosoma y una ventana alrededor de la coordenada X seleccionada.
        // posición del genoma a visualizar en el navegador

        qDebug() << "coordenada pasada a url: " << QString::number(inicio);
        QString url = "https://grch37.ensembl.org/Homo_sapiens/Location/View?r=" +
                       QString::number(cromosoma) + ":" +
                       QString::number(inicio - 500) + "-" + QString::number(inicio + 500) +
                       ";db=core";

        if (ui->wavelet_file->blockCount() != 1 && !to_load && !ui->analiza->isEnabled())
            QDesktopServices::openUrl(QUrl(url, QUrl::TolerantMode));
    }
    // en caso contrario se estima que se trata de una acción de zoom
    else
    {
        page_pressed = true;

        // calculo del ancho de datos visualizado
        int pos_inferior = ui->rango_inferior->text().toInt();
        ancho_ventana    = ui->ancho_ventana->text().toInt();

        qDebug() << "posicion inferior inicial: " << pos_inferior <<
                    " -  posición inicial ratón: " << inicio <<
                    " -  posicion final ratón: " << fin <<
                    " -  ancho adn visualizado: " << ancho_ventana <<
                    " -  ancho ventana visualización: " << ui->ventana_opengl->width() ;

        // asignación de nuevos límites
        // ..en interface
        ui->rango_inferior->setText(QString::number(inicio));
        if (fin - inicio > 100)
        {
            ui->rango_superior->setText(QString::number(fin));
            ui->ancho_ventana->setText(QString::number(fin - inicio));
        }
        else
        {
            ui->rango_superior->setText(QString::number(inicio + 100));
            ui->ancho_ventana->setText(QString::number(100));
        }

        // ..actualización de valores en slider y scroll para realizar el análisis de la nueva región
        ui->page->setValue(fin - inicio);
        ui->scroll_adn->setValue(inicio);

        // actualización de señal visualizada
        if (ui->wavelet_file->blockCount() != 1 && !to_load && !ui->analiza->isEnabled())
            dibuja();

        qDebug() << "posicion inferior final: " << ui->rango_inferior->text() <<
                    " -  ancho adn visualizado: " << ui->ancho_ventana->text() <<
                    " -  ancho ventana visualización: " << ui->ventana_opengl->width() ;
    }
}

// ************************************************************************************************
// **************************************ZONA EJECUCION ACCIONES***********************************
// ************************************************************************************************
void HPG_Dhunter::on_load_files_clicked()
{
    // libera memoria GPU
    cuda_end(cuda_data);

    // inhabilita botón de análisis mientras carga muestras
    ui->analiza->setEnabled(false);
    ui->dmrs->setEnabled(false);
    ui->num_CpG_x_region->setEnabled(false);
    ui->threshold->setEnabled(false);
    ui->dmr_dwt_level->setEnabled(false);
    ui->save_dmr_list->setEnabled(false);


    // inicializa los parámetros a enviar a los hilos
    parametros = (QStringList() << QString::number(_forward) <<    // se informa forward reads 0/1
                                   QString::number(_reverse) <<    // se informa reverse reads 0/1
                                   "0" <<                          // se informa del número de cromosoma
                                   "0"                             // se informa del número de hilo asignado
                 );

    ficheros_case    = ui->wavelet_file->toPlainText().split("\n");
    ficheros_control = ui->control_file->toPlainText().split("\n");
    contador         = 0;

    // obtiene el cromosoma a analizar
    // -------------------------------------------------------------------------------------------
    if (!ui->chr_other->text().isEmpty())
    {
        QString lista = ui->chr_other->text();
        lista.replace(QRegularExpression("[a-zA-Z,.-]+"), " ");
        QStringList llista =  lista.split(QRegularExpression("\\s+"));
        int n = 0;
        if (llista.size() > 0)
            while (llista.at(n).toInt() <= 0 && n < llista.size())
                n++;
        cromosoma = llista.at(n).toInt();
    }
    else
        cromosoma = cromosoma_grid;

    // verifica que existe cromosoma
    if (cromosoma == 0)
    {
        QMessageBox::warning(this,
                                     "ERROR: no chromosome selected",
                                     "Please, select a chromosome from grid or other"
                                     );
    }
    else
    {
        // se comprueba que el cromosoma a cargar existe
        // en todos los directorios seleccionados
        bool chr_existe = false;
        foreach(auto n, ficheros_case)
        {
            QFile archivo;
            QString file = n + "/methylation_map_" +
                           (_forward? "forward_" : "reverse_") +
                           QString::number(cromosoma) + ".csv";

            archivo.setFileName(file);
            if (!archivo.open(QIODevice::ReadOnly))
            {
                chr_existe = false;
                break;
            }
            else
            {
                chr_existe = true;
                archivo.close();
            }
        }

        if (chr_existe)
        {
            foreach(auto n, ficheros_control)
            {
                QFile archivo;
                QString file = n + "/methylation_map_" +
                               (_forward? "forward_" : "reverse_") +
                               QString::number(cromosoma) + ".csv";

                archivo.setFileName(file);
                if (!archivo.open(QIODevice::ReadOnly))
                {
                    chr_existe = false;
                    break;
                }
                else
                {
                    chr_existe = true;
                    archivo.close();
                }
            }
        }

        if (!chr_existe)
        {
            QMessageBox::warning(this,
                                         "ERROR: some files does not exist",
                                         "Please, verify the chromosome number\n"
                                         "or the files into the folders"
                                         );
            return;
        }


        // solicita verificación de cromosoma a analizar
        // y procede a la carga de ficheros
        QMessageBox::StandardButton reply;
        reply = QMessageBox::question(this,
                                      "HPG-Dhunter - visualizer",
                                      "The chromosome to analyze is:\n" +
                                      QString::number(cromosoma) +
                                      "\nIs it right?",
                                      QMessageBox::Yes|QMessageBox::No);

        // confirmación de inicio de proceso de mapeado
        if (reply == QMessageBox::Yes)
        {
            // inhabilita los botones de modificar listado
            ui->load_files->setEnabled(false);
            ui->delete_file->setEnabled(false);
            ui->up_file->setEnabled(false);
            ui->down_file->setEnabled(false);
            ui->delete_control->setEnabled(false);
            ui->up_control->setEnabled(false);
            ui->down_control->setEnabled(false);


            // limpia la matrix de datos del cromosoma anterior
            ui->statusBar->showMessage("freeing memory");
            for (vector<vector<double>> n : mc)
            {
                for (vector<double> m : n)
                {
                        m.clear();
                        m.shrink_to_fit();
                        vector<double>().swap(m);
                }
                n.clear();
                n.shrink_to_fit();
                vector<vector<double>>().swap(n);
            }

            mc.clear();
            mc.shrink_to_fit();
            vector<vector<vector<double>>>().swap(mc);

            limite_inferior = 100000000;        // control de límite inferior
            limite_superior = 0;                // control de límite superior

            // borra todos los posibles hilos creados anteriormente
            foreach(QThread *i, hilo_files_worker)
                delete i;

            hilo_files_worker.clear();
            hilo_files_worker.shrink_to_fit();
            files_worker.clear();
            files_worker.shrink_to_fit();

            ui->statusBar->showMessage("loading files...");

            // crea nuevos hilos para lectura de ficheros
            for (int i = 0; i < ficheros_case.size() + ficheros_control.size(); i++)
            {
                hilo_files_worker.append(new QThread());
                files_worker.append(new Files_worker());
            }

            // moviendo workers a los hilos y conexiones para el proceso de lectura de ficheros
            for (int i = 0; i < hilo_files_worker.size(); i++)
            {
                files_worker[i]->moveToThread(hilo_files_worker[i]);
                connect(files_worker[i], SIGNAL(fichero_leido(int, int, int, int)), SLOT(fichero_leido(int, int, int, int)));
                connect(hilo_files_worker[i], &QThread::finished, files_worker[i], &QObject::deleteLater);
                files_worker[i]->connect(hilo_files_worker[i], SIGNAL(started()), SLOT(lectura()));
                hilo_files_worker[i]->connect(files_worker[i],SIGNAL(lectura_solicitada()), SLOT(start()));
                hilo_files_worker[i]->connect(files_worker[i], SIGNAL(finished()), SLOT(quit()), Qt::DirectConnection);

                // arranque del hilo de lectura
                if (hilo_files_worker[i]->isRunning())
                    hilo_files_worker[i]->wait();

                // se le asigna el número de cromosoma
                parametros[2] = QString::number(cromosoma);
                parametros[3] = QString::number(i);
                qDebug() << "cromosoma a leer:" << parametros[2] << parametros[3];

                // se lanza el hilo de lectura de ficheros para el cromosoma seleccionado
                files_worker[i]->solicitud_lectura(ficheros_case, ficheros_control, parametros, mc_aux, mutex);
            }

            // se lanza el hilo de carga de referencias genéticas, si se dispone de ellas
            switch (ui->genome_reference->currentIndex())
            {
            case 0:
                break;
            case 1:
                hilo_refGen   = new QThread();
                refGen_worker = new RefGen();
                refGen_worker->moveToThread(hilo_refGen);
                connect(refGen_worker, SIGNAL(terminado(ulong)), SLOT(refGen_worker_acabado(ulong)));
                connect(hilo_refGen, &QThread::finished, refGen_worker, &QObject::deleteLater);
                refGen_worker->connect(hilo_refGen, SIGNAL(started()), SLOT(lectura()));
                hilo_refGen->connect(refGen_worker,SIGNAL(lectura_solicitada()), SLOT(start()));
                hilo_refGen->connect(refGen_worker, SIGNAL(finished()), SLOT(quit()), Qt::DirectConnection);

                refGen_worker->solicitud_lectura(cuda_data, parametros[2].toInt());
                break;
            default:
                ;
            }
        }
        else
            return;
    }
}


// ************************************************************************************************
void HPG_Dhunter::on_analiza_clicked()
{
    // deshabilita el pulsador de analiza, carga de ficheros y umbral de cobertura mínima
    ui->analiza->setEnabled(false);
    ui->load_files->setEnabled(false);

    // comprueba la dimensión de la matriz mc para ver si caben todos los datos en GPU
    // en caso contrario, avisa al usuario para que seleccione un conjunto de muestras adecuado.
    uint dimension = limite_superior - limite_inferior + 1;
    if ((dimension & 0x01) == 1)
        dimension++;

    int casos = int(count(visualiza_casos.begin(), visualiza_casos.end(), true));
    int control = int(count(visualiza_control.begin(), visualiza_control.end(), true));
    int samples2visualize = casos + control;
    uint tamanyo = (dimension * sizeof(float) * uint(samples2visualize)) / (1024 * 1024); // tamaño en MiB

    if (tamanyo > 0.25 * memory_available)
    {
        QMessageBox::warning(this,
                             "ERROR: too much samples to visualize",
                             "The GPU memory available will overloaded\n"
                             "Please, select a minor number of samples\n"
                             "to visualize"
                            );
        return;
    }

    qDebug() << "numero de samples a visualizar:" << samples2visualize;


    // actualiza estructura de datos
    cuda_data.samples        = samples2visualize;                       // número de ficheros a analizar
    cuda_data.sample_num     = dimension;                               // cantidad de datos por fichero
    cuda_data.rango_inferior = 0;                                       // primer valor cromosoma
    cuda_data.rango_superior = limite_superior - limite_inferior;       // último valor
    cuda_data.levels         = ui->dmr_dwt_level->value();              // número de niveles a transformar
    cuda_data.data_adjust    = 0;                                       // ajuste desfase en división por nivel para número impar de datos
    cuda_data.pitch          = 0;
    cuda_data.pitch_2        = 0;
    cuda_data.h_haar_L.clear();

    // informa de proceso de lectura de datos en marcha
    ui->statusBar->showMessage("forming the full ADN segment...");

    // borra la memoria previa utilizada
    if (cuda_data.mc_full != nullptr)
    {
        delete [] cuda_data.mc_full[0];
        delete [] cuda_data.mc_full;
    }

    if (cuda_data.h_haar_C != nullptr)
    {
        delete [] cuda_data.h_haar_C[0];
        delete [] cuda_data.h_haar_C;
    }

    // crea matriz ampliada -------------------------------------------------------------------
    //      -> vectores con todas las posiciones contiguas
    //      -> con ceros en las posiciones sin metilación
    // reserva TODA la memoria CONTIGUA con todos los datos de todas las muestras
    // para trasvase de datos entre GPU y CPU con CUDA, la matriz debe ser contigua completa
    // reserva la memoria para la matriz de datos extendida
    cuda_data.mc_full = new float*[cuda_data.samples];
    cuda_data.mc_full[0] = new float[uint(cuda_data.samples) * cuda_data.sample_num];
    for (int i = 1; i < cuda_data.samples; i++)
            cuda_data.mc_full[i] = cuda_data.mc_full[i - 1] + cuda_data.sample_num;

    posicion_metilada.clear();
    posicion_metilada.assign(uint(cuda_data.samples), vector<uint> (cuda_data.sample_num, 0));

    // copia de todos los datos a la matriz ampliada
    // --------------------------------------------------------------------------------------------
    uint muestra_seleccionada = 0;
    uint cuenta_posiciones_metiladas;
    uint pos_met;
    for (uint m = 0; m < mc.size(); m++)
    {
        // comprueba si la posición m corresponde a caso o control
        if (uint(mc[m][0][11]) == 0)
        {
            cuenta_posiciones_metiladas = 0;
            pos_met                     = 1;

            // comprueba si la posición está seleccionada para visualizar
            if (visualiza_casos[uint(mc[m][0][10])])
            {
                // rellena con ceros la matriz de datos
                for (uint n = 0; n < cuda_data.sample_num + 1; ++n)
                    cuda_data.mc_full [muestra_seleccionada][n] = 0.0;

                // rellena con datos las posiciones metiladas si supera la cobertura
                for (uint k = 0; k < mc[m].size(); k++)
                    if(mc[m][k][ui->mC->isChecked() ? 2 : 8] >= ui->cobertura->value())
                    {
                        cuda_data.mc_full [muestra_seleccionada][uint(mc[m][k][0]) - limite_inferior] = float(mc[m][k][ui->mC->isChecked() ? 1 : 7]);

                        cuenta_posiciones_metiladas++;

                        while (pos_met < uint(mc[m][k][0]) - limite_inferior)
                        {
                            posicion_metilada[muestra_seleccionada][pos_met] = cuenta_posiciones_metiladas - 1;
                            pos_met++;
                        }
                        posicion_metilada[muestra_seleccionada][uint(mc[m][k][0]) - limite_inferior] = cuenta_posiciones_metiladas;
                    }


                while (pos_met < posicion_metilada[muestra_seleccionada].size())
                {
                    posicion_metilada[muestra_seleccionada][pos_met] = posicion_metilada[muestra_seleccionada][pos_met - 1];
                    pos_met++;
                }

                h_haar_C_distribucion.push_back(0);
                muestra_seleccionada++;
            }
        }
        else
        {
            cuenta_posiciones_metiladas = 0;
            pos_met                     = 1;

            // comprueba si la posición está seleccionada para visualizar
            if (visualiza_control[uint(mc[m][0][10])])
            {
                // rellena con ceros la matriz de datos
                for (uint n = 0; n < cuda_data.sample_num + 1; ++n)
                    cuda_data.mc_full [muestra_seleccionada][n] = 0.0;

                // rellena con datos las posiciones metiladas si supera la cobertura
                for (uint k = 0; k < mc[m].size(); k++)
                    if(mc[m][k][ui->mC->isChecked() ? 2 : 8] >= ui->cobertura->value())
                    {
                        cuda_data.mc_full [muestra_seleccionada][uint(mc[m][k][0]) - limite_inferior] = float(mc[m][k][ui->mC->isChecked() ? 1 : 7]);

                        cuenta_posiciones_metiladas++;

                        while (pos_met < uint(mc[m][k][0]) - limite_inferior)
                        {
                            posicion_metilada[muestra_seleccionada][pos_met] = cuenta_posiciones_metiladas - 1;
                            pos_met++;
                        }
                        posicion_metilada[muestra_seleccionada][uint(mc[m][k][0]) - limite_inferior] = cuenta_posiciones_metiladas;
                    }


                while (pos_met < posicion_metilada[muestra_seleccionada].size())
                {
                    posicion_metilada[muestra_seleccionada][pos_met] = posicion_metilada[muestra_seleccionada][pos_met - 1];
                    pos_met++;
                }

                h_haar_C_distribucion.push_back(1);
                muestra_seleccionada++;
            }
        }
    }

    qDebug() << "matriz de datos totalmente llena" << muestra_seleccionada;

    // actualización de límites, rangos de sliders y demás datos de interfaz
    // --------------------------------------------------------------------------------------------

    // ancho inicial ventana muestras al 10% longitud muestra, limites modificados  - - - - - -
    ui->scroll_adn->setPageStep(int((limite_superior - limite_inferior) * 0.1));
    ui->scroll_adn->setMaximum(int(limite_superior) - ui->scroll_adn->pageStep());
    ui->scroll_adn->setMinimum(int(limite_inferior));
    ui->scroll_adn->setValue(int(limite_inferior));

    // rango ancho de ventana en slider selector de ancho - - - - - - - - - - - - - - - - - - -
    ui->page->setMaximum(int(limite_superior - limite_inferior));
    ui->page->setSingleStep(5000);
    ui->page->setPageStep(100000);
    ui->page->setValue(ui->scroll_adn->pageStep());

    // etiquetas información límites y ancho ventana  - - - - - - - - - - - - - - - - - - - - -
    ui->limite_inferior->setText(QString::number(limite_inferior) + " - ");
    ui->limite_superior->setText(" - " + QString::number(limite_superior));
    ui->ancho_ventana->setText((QString::number(ui->scroll_adn->pageStep())));
    ui->rango_inferior->setText(QString::number(ui->scroll_adn->value()));
    ui->rango_superior->setText(QString::number(ui->scroll_adn->value() + ui->scroll_adn->pageStep()));

    // actualiza rango de slider selector de nivel de transformación  - - - - - - - - - - - - -
    ui->slider_nivel->setMaximum(int(log2(ui->scroll_adn->pageStep())));
    ui->slider_nivel->setValue(ui->slider_nivel->maximum() - 3);

    // envía los datos a la memoria global de la GPU
    // --------------------------------------------------------------------------------------------
    // informa de proceso de copia de datos a memoria GPU
    ui->statusBar->showMessage("loading full segments to GPU memory...");

    // libera la memoria de la GPU
    cuda_end(cuda_data);

    qDebug() << "memoria GPU liberada";

    // envía el total de los datos a la GPU
    cuda_send_data(cuda_data);

    qDebug() << "datos cargados en GPU";

    // informa de proceso de lectura de datos finalizado
    ui->statusBar->showMessage("samples loaded, ready to transform -> data in GPU global memory");

    // habilita detección de DMRs
    ui->dmrs->setEnabled(true);
    ui->num_CpG_x_region->setEnabled(true);
    ui->threshold->setEnabled(true);
    ui->dmr_dwt_level->setEnabled(true);
    ui->save_dmr_list->setEnabled(false);
    to_load = false;

    // deshabilita el movimiento de los ficheros para preservar la correspondiencia
    // de los colores de las líneas con los colores de la gráfica.
    ui->up_file->setEnabled(false);
    ui->down_file->setEnabled(false);
    ui->delete_file->setEnabled(false);
    ui->up_control->setEnabled(false);
    ui->down_control->setEnabled(false);
    ui->delete_control->setEnabled(false);

    // limpia la ventana de DMRs
    dmr_listo = false;
    ui->dmr_position->clear();
    ui->dmr_detail->clear();
    switch (ui->genome_reference->currentIndex())
    {
    case 0:
        ui->label_6->setText("DMRs  -  unknown genome reference");
        break;
    case 1:
        ui->label_6->setText("DMRs  -  GENE-names -  distance (for " + ui->genome_reference->currentText() + ")");
        break;
    default:
        ;
    }
    dmr_listo = true;

    dibuja();

}

// ************************************************************************************************
void HPG_Dhunter::dibuja()
{
//    INIT_TIMER

    // inicializa la estructura de variables para el segmento a analizar ----------------------
    cuda_data.h_haar_L.clear();                                     // vector con número de datos por nivel
    cuda_data.sample_num  = ui->ancho_ventana->text().toUInt();     // número de datos por muestra
    cuda_data.levels      = ui->slider_nivel->value();              // número de niveles a transformar
    cuda_data.data_adjust = 0;                                      // ajuste desfase en división por nivel para número impar de datos

//    START_TIMER

    // informa de proceso en barra inferior ---------------------------------------------------
    ui->statusBar->showMessage("working ...");

    cuda_data.d_max[0] = 1.0;

    // tranforma la ventana de datos correspondiente ------------------------------------------
    cuda_calculo_haar_L(cuda_data);

    ui->ventana_opengl->setNum_L(cuda_data.h_haar_L[0]);
    ui->ventana_opengl->setNum_samples(cuda_data.samples, int(count(visualiza_casos.begin(), visualiza_casos.end(), true)));

    ui->ventana_opengl->registerBuffer();
    ui->ventana_opengl->mapResource(cuda_data);

    cuda_main(cuda_data);

    ui->ventana_opengl->unmapResource();

    // plotea la señal transformada -----------------------------------------------------------
    ui->ventana_opengl->update();   // con QOpenGLWidget

    qDebug() << "datos último nivel: " << cuda_data.h_haar_L[0] << " ------ valor máximo escala dibujado:" << cuda_data.d_max[0];

    // actualiza valores de escala
    ui->escala_1->setText(QString::number(double(cuda_data.d_max[0]), 'f', 2));
    ui->escala_075->setText(QString::number(double(cuda_data.d_max[0]) * 0.75, 'f', 2));
    ui->escala_05->setText(QString::number(double(cuda_data.d_max[0]) * 0.5, 'f', 2));
    ui->escala_025->setText(QString::number(double(cuda_data.d_max[0]) * 0.25, 'f', 2));



//    STOP_TIMER("tiempo dibujado")

    // informa de proceso en barra inferior ---------------------------------------------------
    ui->statusBar->showMessage("transform finished -> data in GPU and CPU memory. Plotted");

    // deshabilita el pulsador de analiza, carga de ficheros y umbral de cobertura mínima
    ui->analiza->setEnabled(false);
    ui->load_files->setEnabled(false);
}

// ************************************************************************************************
void HPG_Dhunter::on_dmrs_clicked()
{
    INIT_TIMER

    QString linea = "";

    // informa de proceso en barra inferior ---------------------------------------------------
    ui->statusBar->showMessage("Working on matrix of chromosome samples wavelet values differences");

    ui->dmrs->setEnabled(false);
    ui->save_dmr_list->setEnabled(true);
    ui->open_file->setFocus();

    // inicializa la estructura de variables para el segmento a analizar ----------------------
    START_TIMER

    // limpia matriz de resultados de procesamiento en GPU
    foreach (vector<float> n, h_haar_C)
    {
        n.clear();
        n.shrink_to_fit();
        vector<float>().swap(n);
    }
    h_haar_C.clear();
    h_haar_C.shrink_to_fit();
    vector<vector<float>>().swap(h_haar_C);

    // actualiza datos
    cuda_data.h_haar_L.clear();                                         // vector con número de datos por nivel
    cuda_data.rango_inferior = 0;                                       // primer valor cromosoma
    cuda_data.rango_superior = limite_superior - limite_inferior;       // último valor
    cuda_data.sample_num     = limite_superior - limite_inferior + 1;   // número de datos por muestra
    cuda_data.levels         = ui->dmr_dwt_level->value();              // número de niveles a transformar
    cuda_data.data_adjust    = 0;                                       // ajuste desfase en división por nivel para número impar de datos

    // tranforma la ventana de datos correspondiente ------------------------------------------
    cuda_calculo_haar_L(cuda_data);

    ui->ventana_opengl->setNum_L(cuda_data.h_haar_L[0]);
    ui->ventana_opengl->registerBuffer();
    ui->ventana_opengl->mapResource(cuda_data);

    cuda_main(cuda_data);

    ui->ventana_opengl->unmapResource();

    // recoge los resultados en una matriz, acumulando todos los resultados
    vector<float> aux(ulong(cuda_data.h_haar_L[0]), 0.0);
    for (int i = 0; i < cuda_data.samples; i++)
    {
        for (size_t j = 0; j < size_t(cuda_data.h_haar_L[0]); j++)
            aux[j] = cuda_data.h_haar_C[i][j];

        h_haar_C.push_back(aux);
    }

    // restituye la estructura de variables para el segmento analizado ------------------------
    cuda_data.rango_inferior = uint(ui->scroll_adn->value()) - limite_inferior;
    cuda_data.rango_superior = uint(ui->scroll_adn->value()) - limite_inferior + uint(ui->scroll_adn->pageStep());
    cuda_data.levels         = ui->slider_nivel->value();

    // reserva la matriz de diferencias de medias por grupos de control y casos
    if (dmr_diff != nullptr)
        delete[] dmr_diff;

    dmr_diff = new float[cuda_data.h_haar_L[0]];
    for (int i = 0; i < cuda_data.h_haar_L[0]; i++)
        dmr_diff[i] = 0.0;

    uint paso = uint(pow(2, ui->dmr_dwt_level->value()));

    // realiza el cálculo de medias de las muestras de control de los casos
    for (uint m = 0; m < uint(cuda_data.h_haar_L[0]); m++) // ...en cada posición
    {
        float media_casos   = 0.0;
        float media_control = 0.0;
        uint numero_casos   = 0;
        uint numero_control = 0;
        for (uint i = 0; i < uint(cuda_data.samples); i++)
        {
            if (h_haar_C_distribucion[i])
            {
                if (((m + 1) * paso < posicion_metilada[i].size() ? posicion_metilada[i][(m + 1) * paso] : posicion_metilada[i].back()) >
                        posicion_metilada[i][m * paso] + (paso * uint(ui->num_CpG_x_region->value()) / 100))
                {
                    media_control += h_haar_C[i][m];
                    numero_control++;
                }
            }
            else
            {
                if (((m + 1) * paso < posicion_metilada[i].size() ? posicion_metilada[i][(m + 1) * paso] : posicion_metilada[i].back()) >
                        posicion_metilada[i][m * paso] + (paso * uint(ui->num_CpG_x_region->value()) / 100))
                {
                    media_casos += h_haar_C[i][m];
                    numero_casos++;
                }
            }
        }

        if (numero_casos > 0 && numero_control > 0)
            dmr_diff[m] = (media_casos / numero_casos) - (media_control / numero_control);
    }


    dmr_diff_cols = uint(cuda_data.h_haar_L[0]);

    // llama a función de cálculo de diferencias entre muestras para ponerlas en la tabla
    hallar_dmrs();

    STOP_TIMER("análisis DMR")
}

// ************************************************************************************************
void HPG_Dhunter::hallar_dmrs()
{
    QString linea      = "";
    uint paso          = uint(pow(2, ui->dmr_dwt_level->value()));
    uint *posicion_dmr = new uint[dmr_diff_cols];    // crear array de posición

    // informa de proceso en barra inferior ---------------------------------------------------
    ui->statusBar->showMessage("Looking for DMRs");

    // encontrar DMRs en función del threshold establecido ------------------------------------
    // llenar de 0s el vector de posicion de DMRs
    for (uint p = 0; p < dmr_diff_cols; p++)
        posicion_dmr[p] = 0;

    // rellenar las posiciones con diferencias válidas
    for (uint j = 0; j < dmr_diff_cols; j++)
        if (dmr_diff[j] > threshold || dmr_diff[j] < -threshold)
            posicion_dmr[j] = j * paso + limite_inferior;

    // rellenar ventana de datos - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    dmr_listo = false;
    ui->dmr_position->clear();
    dmrs.clear();

    // busca y rellena la lista de DMRs
    int inicio = 0;
    int fin    = int(num_genes);
    for (uint p = 0; p < dmr_diff_cols; p++)
    {
        if (posicion_dmr[p] >= limite_inferior)
        {
            linea.clear();
            uint q = p;     // para ayuda en la zona de detección de referencia de genoma

            // busca las posiciones inicial y final de la DMR
            //-----------------------------------------------
            linea.append(QString::number(posicion_dmr[p]));

            while (p + 1 < dmr_diff_cols && posicion_dmr[p + 1] >= limite_inferior)
               p++;

            linea.append("-" + QString::number(uint(posicion_dmr[p]) + paso));


            // búsqueda del nombre del GEN implicado o más cercano a los DMRs encontrados
            //---------------------------------------------------------------------------
            // se realiza sobre datos de la genome.ucsc.edu data base sobre genes conocidos
            // ..previamente se han cargado los nombres y posiciones de los genes correspondientes
            // al cromosoma que se está analizando
            // ..por búsqueda binaria sobre este fichero se determina el nombre del gen.
            switch (ui->genome_reference->currentIndex())
            {
                case 0:
                    break;

                case 1:
                    int mitad        = inicio;
                    bool match       = false;
                    uint gen_ini     = 0;
                    uint gen_ant_fin = uint(stoul(cuda_data.refGen[0][4]));

                    while (uint(stoul(cuda_data.refGen[mitad][3])) < posicion_dmr[q] && mitad < fin - 1)
                        mitad++;

                    gen_ini = uint(stoul(cuda_data.refGen[mitad][3]));
                    if (mitad > 0)
                        gen_ant_fin = uint(stoul(cuda_data.refGen[mitad - 1][4]));

                    // el inicio dmr es igual que inicio del gen
                    if (gen_ini == posicion_dmr[q])
                    {
                        match = true;
                        linea.append(" " + QString::fromStdString(cuda_data.refGen[mitad][0]) +
                                     " " + QString::fromStdString(cuda_data.refGen[mitad][1]) +
                                     " 0");
                    }
                    // el inicio dmr es menor que inicio del gen pero el final dmr es mayor que el inicio del gen
                    else if (gen_ini <= posicion_dmr[p] + paso)
                    {
                        match = true;
                        linea.append(" " + QString::fromStdString(cuda_data.refGen[mitad][0]) +
                                     " " + QString::fromStdString(cuda_data.refGen[mitad][1]) +
                                     " -" + QString::number(stoul(cuda_data.refGen[mitad][3]) - posicion_dmr[q]));
                    }
                    // el inicio dmr es mayor que inicio del gen anterior pero es menor que el fin del gen anterior
                    else if (gen_ant_fin > posicion_dmr[q])
                    {
                        match = true;
                        linea.append(" " + QString::fromStdString(cuda_data.refGen[mitad > 0? mitad - 1 : 0][0]) +
                                     " " + QString::fromStdString(cuda_data.refGen[mitad > 0? mitad - 1 : 0][1]) +
                                     " +" + QString::number(posicion_dmr[q] - stoul(cuda_data.refGen[mitad > 0? mitad - 1 : 0][3])));
                    }


                    // si se encuentra entre genes, ver de qué gen está más cerca
                    // se elige la distancia más pequeña entre:
                    // ..distancia inicio dmr y fin gen anterior
                    // ..distancia fin dmr e inicio gen posterior
                    if (!match)
                    {
                        ulong dif1 = posicion_dmr[q] - gen_ant_fin;
                        ulong dif2 = gen_ini - posicion_dmr[p] + paso;

                        if (dif1 >= dif2)
                        {
                            linea.append(" " + QString::fromStdString(cuda_data.refGen[mitad][0]) +
                                         " " + QString::fromStdString(cuda_data.refGen[mitad][1]) +
                                         " --" + QString::number(dif2));
                        }
                        else
                        {
                            linea.append(" " + QString::fromStdString(cuda_data.refGen[mitad > 0? mitad - 1 : 0][0]) +
                                         " " + QString::fromStdString(cuda_data.refGen[mitad > 0? mitad - 1 : 0][1]) +
                                         " ++" + QString::number(dif1));
                        }
                    }

                    if (mitad > 0)
                        inicio = mitad - 1;
                    break;
            }


            // define si está hipermetilado o hipometilado el control frente al caso
            //-----------------------------------------------------------------------
            linea.append((dmr_diff[p] > 0)? " hiper" : " hipo");

            // añade resultado de análisis DWT
            //-----------------------------------------------------------------------
            linea.append(" " + QString::number(double(dmr_diff[q])));
            // añade la información a la lista de DMRs
            //-----------------------------------------------------------------------
            ui->dmr_position->appendPlainText(linea);

            // añade la información de posicion al vector de DMRs
            //-----------------------------------------------------------------------
            linea.append("//" + QString::number(q) + " " + QString::number(p));
            dmrs.append(linea);
        }
    }

    delete[] posicion_dmr;

    // selección de DMR desde listado
    qDebug() << "tamaño del fichero que guarda los dmrs localizados: " << dmrs.size() << "x" << dmrs[0].size();
    dmr_listo = true;

    // informa de proceso en barra inferior ---------------------------------------------------
    ui->statusBar->showMessage(QString::number(dmrs.size()) + " DMRs found");
    switch (ui->genome_reference->currentIndex())
    {
    case 0:
        ui->label_6->setText(QString::number(dmrs.size()) + " DMRs found | range - methylation - dwt-diff (for unknown)");
        break;
    case 1:
        ui->label_6->setText(QString::number(dmrs.size()) + " DMRs found | range - GENE-names - distance - methylation - dwt-diff (for " + ui->genome_reference->currentText() + ")");
        break;
    default:
        ;
    }
    ui->match->setEnabled(true);
    ui->match->setChecked(false);
    ui->save_dmr_list->setEnabled(true);

    ui->slider_nivel->setValue(ui->dmr_dwt_level->value());
}

// ************************************************************************************************
void HPG_Dhunter::on_dmr_position_cursorPositionChanged()
{
    if (dmr_listo && ui->dmr_position->blockCount() >= 1)
    {
        QString linea, linea_detail;
        uint pos_inf;
        uint pos_sup;
        uint ancho_dmr;
        QTextBlockFormat color;

        uint entorno = 2;            // zona lateral de dmr detectada para mostrar centrada
        uint paso    = uint(pow(2, ui->dmr_dwt_level->value()));

        // desmarcar la línea previa quitando color de fondo
        color.setBackground(Qt::white);
        cursor->select((QTextCursor::LineUnderCursor));
        cursor->setBlockFormat(color);

        // adquirir el cursor de la línea seleccionada
        *cursor = ui->dmr_position->textCursor();
        cursor->movePosition(QTextCursor::StartOfBlock);
        cursor->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);

        // marcar con color de fondo la línea seleccionada
        color.setBackground(Qt::yellow);
        cursor->select((QTextCursor::LineUnderCursor));
        cursor->setBlockFormat(color);

        linea     = cursor->selectedText();
        pos_inf   = linea.split("-")[0].toUInt() - entorno * paso;
        pos_sup   = linea.split("-")[1].split(" ")[0].toUInt() + (entorno+1) * paso;


        // actualizar información de interfaz
        //------------------------------------
        ui->page->setValue(int(pos_sup - pos_inf));
        ui->scroll_adn->setValue(int(pos_inf));
        if (ui->slider_nivel->value() >= 9)
            ui->slider_nivel->setValue(9);
        ui->rango_inferior->setText(QString::number(pos_inf));
        ui->rango_superior->setText(QString::number(pos_sup));
        ui->ancho_ventana->setText(QString::number(pos_sup - pos_inf));
        ui->fine_tunning->setMinimum(-1 * ui->ancho_ventana->text().toInt());
        ui->fine_tunning->setValue(0);

        // posición inicial y final de zona dwt en h_haar_C correspondiente al DMR identificado
        int num_linea = cursor->blockNumber();
        qDebug() << "numero de linea dmr seleccionada en la ventana:" << num_linea;
        uint pos_dwt_ini = dmrs.at(num_linea).split("//")[1].split(" ")[0].toUInt();
        uint pos_dwt_fin = dmrs.at(num_linea).split("//")[1].split(" ")[1].toUInt();


        // rellena la ventana de información del DMR seleccionado "dmr_detail"
        //--------------------------------------------------------------------
        pos_inf   = linea.split("-")[0].toUInt();
        pos_sup   = linea.split("-")[1].split(" ")[0].toUInt();
        ancho_dmr = linea.split("-")[1].split(" ")[0].toUInt() - linea.split("-")[0].toUInt();

        ui->dmr_detail->clear();
        // una línea de información por cada muestra
        for (uint j = 0; j < mc.size(); j++)
        {
            if (int(mc[j][0][11]) == 0)
            {
                if (visualiza_casos.at(uint(mc[j][0][10])))
                {
                    int cobertura_minima = 500000;
                    int cobertura_maxima = 0;
                    int cobertura_media  = 0;
                    int distancia_minima = 500000;
                    int distancia_maxima = 0;
                    int distancia_media  = 0;
                    int sites_C          = 0;
                    int sites_nC         = 0;
                    int sites_mC         = 0;
                    int sites_hmC        = 0;
                    int posiciones       = 0;
                    float dwt_valor      = 0.0;
                    float ratio_medio    = 0.0;

                    linea_detail.clear();

                    // busca la posición inical
                    uint posicion = 0;
                    while (pos_inf > mc[j][posicion][0])
                        posicion++;

                    // búsqueda de valores a lo largo del DMR
                    while (pos_sup > mc[j][posicion][0] && posicion < mc[j].size())
                    {
                        // cobertura
                        if (cobertura_minima >= mc[j][posicion][ui->mC->isChecked() ? 2 : 8])
                            cobertura_minima = int(mc[j][posicion][ui->mC->isChecked() ? 2 : 8]);
                        if (cobertura_maxima < mc[j][posicion][ui->mC->isChecked() ? 2 : 8])
                            cobertura_maxima = int(mc[j][posicion][ui->mC->isChecked() ? 2 : 8]);
                        cobertura_media += mc[j][posicion][ui->mC->isChecked() ? 2 : 8];
                        if (mc[j][posicion][ui->mC->isChecked() ? 2 : 8] > 0)
                            ratio_medio += float(mc[j][posicion][ui->mC->isChecked() ? 1 : 7]);


                        // distancia
                        if (posicion + 2 < mc[j].size() && ancho_dmr > mc[j][posicion + 1][0] - mc[j][posicion][0])
                        {
                            if (distancia_minima >= mc[j][posicion + 1][0] - mc[j][posicion][0])
                                distancia_minima = int(mc[j][posicion + 1][0] - mc[j][posicion][0]);
                            if (distancia_maxima < mc[j][posicion + 1][0] - mc[j][posicion][0])
                                distancia_maxima = int(mc[j][posicion + 1][0] - mc[j][posicion][0]);
                            distancia_media += mc[j][posicion + 1][0] - mc[j][posicion][0];
                        }

                        // número de posiciones detectadas por tipo de mononucleótico
                        sites_C   += (mc[j][posicion][3] > 0) ? 1 : 0;
                        sites_nC  += (mc[j][posicion][4] > 0) ? 1 : 0;
                        sites_mC  += (mc[j][posicion][5] > 0) ? 1 : 0;
                        sites_hmC += (mc[j][posicion][6] > 0) ? 1 : 0;

                        // número de posiciones detectadas con algún tipo de nucleótido sensible
                        if (mc[j][posicion][ui->mC->isChecked() ? 2 : 8] > 0)
                            posiciones++;

                        posicion++;
                    //    qDebug() << j << " -> " << mc[j].size() << " - " << posicion;
                    }

                    // valor medio dwt en la región identificada
                    for (uint i = pos_dwt_ini; i <= pos_dwt_fin; i++)
                        dwt_valor += h_haar_C[j][i];                    // (h_haar_C[f][i] / (pos_dwt_fin - pos_dwt_ini + 1));
                    dwt_valor /= (pos_dwt_fin - pos_dwt_ini + 1);

                    // carga de resultado en línea de texto para mostrar
                    linea_detail.append(ficheros_case.at(int(mc[j][0][10])).split("/").back() +
                                        " " +   QString("%1").arg(double(dwt_valor)) +
                                        " " +   QString("%1").arg(posiciones > 1 ? double(ratio_medio / posiciones) : double(ratio_medio)) +
                                        " | " + QString::number(posiciones) +
                                        " | " + QString::number(cobertura_minima) +
                                        " - " + QString::number((posiciones > 1) ? cobertura_media / posiciones : cobertura_media) +
                                        " - " + QString::number(cobertura_maxima) +
                                        " | " + QString::number(sites_C) +
                                        " - " + QString::number(sites_nC) +
                                        " - " + QString::number(sites_mC) +
                                        " - " + QString::number(sites_hmC) +
                                        " | " + QString::number(distancia_minima) +
                                        " - " + QString::number((posiciones > 1) ? distancia_media / posiciones : distancia_media) +
                                        " - " + QString::number(distancia_maxima));

                    ui->dmr_detail->appendPlainText(linea_detail);
                }
            }
        }

        for (uint j = 0; j < mc.size(); j++)
        {
            if (int(mc[j][0][11]) == 1)
            {
                if (visualiza_control.at(uint(mc[j][0][10])))
                {
                    int cobertura_minima = 500000;
                    int cobertura_maxima = 0;
                    int cobertura_media  = 0;
                    int distancia_minima = 500000;
                    int distancia_maxima = 0;
                    int distancia_media  = 0;
                    int sites_C          = 0;
                    int sites_nC         = 0;
                    int sites_mC         = 0;
                    int sites_hmC        = 0;
                    int posiciones       = 0;
                    float dwt_valor      = 0.0;
                    float ratio_medio    = 0.0;

                    linea_detail.clear();

                    // busca la posición inical
                    uint posicion = 0;
                    while (pos_inf > mc[j][posicion][0])
                        posicion++;

                    // búsqueda de valores a lo largo del DMR
                    while (pos_sup > mc[j][posicion][0] && posicion < mc[j].size())
                    {
                        // cobertura
                        if (cobertura_minima >= mc[j][posicion][ui->mC->isChecked() ? 2 : 8])
                            cobertura_minima = int(mc[j][posicion][ui->mC->isChecked() ? 2 : 8]);
                        if (cobertura_maxima < mc[j][posicion][ui->mC->isChecked() ? 2 : 8])
                            cobertura_maxima = int(mc[j][posicion][ui->mC->isChecked() ? 2 : 8]);
                        cobertura_media += mc[j][posicion][ui->mC->isChecked() ? 2 : 8];
                        if (mc[j][posicion][ui->mC->isChecked() ? 2 : 8] > 0)
                            ratio_medio += float(mc[j][posicion][ui->mC->isChecked() ? 1 : 7]);


                        // distancia
                        if (posicion + 2 < mc[j].size() && ancho_dmr > mc[j][posicion + 1][0] - mc[j][posicion][0])
                        {
                            if (distancia_minima >= mc[j][posicion + 1][0] - mc[j][posicion][0])
                                distancia_minima = int(mc[j][posicion + 1][0] - mc[j][posicion][0]);
                            if (distancia_maxima < mc[j][posicion + 1][0] - mc[j][posicion][0])
                                distancia_maxima = int(mc[j][posicion + 1][0] - mc[j][posicion][0]);
                            distancia_media += mc[j][posicion + 1][0] - mc[j][posicion][0];
                        }

                        // número de posiciones detectadas por tipo de mononucleótico
                        sites_C   += (mc[j][posicion][3] > 0) ? 1 : 0;
                        sites_nC  += (mc[j][posicion][4] > 0) ? 1 : 0;
                        sites_mC  += (mc[j][posicion][5] > 0) ? 1 : 0;
                        sites_hmC += (mc[j][posicion][6] > 0) ? 1 : 0;

                        // número de posiciones detectadas con algún tipo de nucleótido sensible
                        if (mc[j][posicion][ui->mC->isChecked() ? 2 : 8] > 0)
                            posiciones++;

                        posicion++;
                    }

                    // valor medio dwt en la región identificada
                    for (uint i = pos_dwt_ini; i <= pos_dwt_fin; i++)
                        dwt_valor += h_haar_C[j][i];                    // (h_haar_C[f][i] / (pos_dwt_fin - pos_dwt_ini + 1));
                    dwt_valor /= (pos_dwt_fin - pos_dwt_ini + 1);

                    // carga de resultado en línea de texto para mostrar
                    linea_detail.append(ficheros_control.at(int(mc[j][0][10])).split("/").back() +
                                        " " +     QString("%1").arg(double(dwt_valor)) +
                                        " " +     QString("%1").arg(posiciones > 1 ? double(ratio_medio / posiciones) : double(ratio_medio)) +
                                        " | " + QString::number(posiciones) +
                                        " | " + QString::number(cobertura_minima) +
                                        " - " + QString::number((posiciones > 1) ? cobertura_media / posiciones : cobertura_media) +
                                        " - " + QString::number(cobertura_maxima) +
                                        " | " + QString::number(sites_C) +
                                        " - " + QString::number(sites_nC) +
                                        " - " + QString::number(sites_mC) +
                                        " - " + QString::number(sites_hmC) +
                                        " | " + QString::number(distancia_minima) +
                                        " - " + QString::number((posiciones > 1) ? distancia_media / posiciones : distancia_media) +
                                        " - " + QString::number(distancia_maxima));

                    ui->dmr_detail->appendPlainText(linea_detail);
                }
            }
        }

        // colorear las lineas con el color de la muestra a que pertenecen
        *cursor_files = ui->dmr_detail->textCursor();
        cursor_files->movePosition(QTextCursor::Start);
        for (uint i = 0; i < mc.size(); i++)
        {
            uint cuenta = uint(count(visualiza_casos.begin(), visualiza_casos.end(), 1));
            color.setBackground(Qt::white);
            if (i < cuenta)
                color_char.setForeground(QColor(255, (40 * i <= 255)? int(40 * i) : 255, 0));
            else
                color_char.setForeground(QColor(0, (40 * (i - cuenta) <= 255)? int(40 * (i - cuenta)) : 255, 255));
            cursor_files->movePosition(QTextCursor::StartOfBlock);
            cursor_files->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
            cursor_files->select(QTextCursor::BlockUnderCursor);
            cursor_files->setCharFormat(color_char);
            cursor_files->setBlockFormat(color);
            cursor_files->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);
        }

        // mostrar gráfica de DMR
        // lanza la transformación
        if (ui->wavelet_file->blockCount() != 1)
            dibuja();
    }
}

// ************************************************************************************************
void HPG_Dhunter::on_match_clicked()
{
    dmr_listo = false;
//        ui->dmr_position->clear();    NO FUNCIONA. ERROR OUT OF RANGE ¿?
    switch (ui->genome_reference->currentIndex())
    {
    case 0:
        break;

    case 1:
        ui->dmr_position->setPlainText("");
        if (ui->match->isChecked())
        {
            QStringList dmrs_reducida;
            for (QString linea : dmrs)
            {
                if (!(linea.split("  ").last().contains("--") || linea.split("  ").last().contains("++")))
                    dmrs_reducida.append(linea);
            }
            ui->dmr_position->appendPlainText(dmrs_reducida.join('\n'));
            ui->statusBar->showMessage(QString::number(dmrs_reducida.size()) + " DMRs matching with known genes found");
            ui->label_6->setText(QString::number(dmrs_reducida.size()) + " DMRs found  -  GENE-names -  distance        (for human hg19)");
        }
        else
        {
            ui->dmr_position->appendPlainText(dmrs.join('\n'));
            ui->statusBar->showMessage(QString::number(dmrs.size()) + " DMRs found yeah!");
            ui->label_6->setText(QString::number(dmrs.size()) + " DMRs found  -  GENE-names -  distance        (for human hg19)");
        }
        break;

    default:
        ;
    }
    dmr_listo = true;
}

// ************************************************************************************************
void HPG_Dhunter::on_save_dmr_list_clicked()
{
    //ui->save_dmr_list->setEnabled(false);

    // solicita nombre de fichero y directorio para guardar la lista de dmrs
    fichero = QFileDialog::getSaveFileName( this,
                                            tr("Select path and name to save the DMRs list"),
                                            (directorio) ? path : QDir::homePath() ,
                                            "CSV files (*.csv);; All files (*.*)"
                                            );

    if(fichero.isEmpty() || fichero.isNull())
    {
        qDebug() << "no ha recogido el nombre del fichero. " << fichero;
        fichero = "";
    }
    else
    {
        if (fichero.contains('.'))
        {
            if (fichero.split('.').last() != "csv" && fichero.split('.').last() != "txt")
                fichero.append(".csv");
        }
        else
            fichero.append(".csv");


        QFile data;
        data.setFileName(fichero);

        // comprueba que el fichero se ha abierto correctamente
        if (!data.open(QIODevice::WriteOnly))// | QIODevice::Text))
        {
            QMessageBox::warning(this,
                                 "ERROR Opening files",
                                 "An error occurred opening the file: " + fichero +
                                 "\nPlease, check the file for corrupted"
                                );
            qDebug() << "ERROR opening file: " << fichero;
            return;
        }
        else
        {
            dmr_listo = false;
            QTextStream s(&data);
            if (dmrs[0].size() > 0)
            {
                QString linea, linea_detail;
                uint pos_inf;
                uint pos_sup;
                uint ancho_dmr;

                // encabezado de la información del dmr
                switch (ui->genome_reference->currentIndex())
                {
                case 0:
                    s << "pos_init-pos_end methylation dwt_diff\n";
                    break;
                case 1:
                    s << "pos_init-pos_end name_1 name_2 distance methylation dwt_diff\n";
                    break;
                }

                // añade una línea por dmr detectaado y línea de características por muestra en cada dmr
                for (int i = 0; i < dmrs.size(); i++)
                {
                    // escribe zona dmr detectada
                    s << dmrs.at(i).split("//")[0] << '\n';

                    // posición inicial y final de zona dwt en h_haar_C correspondiente al DMR identificado
                    uint pos_dwt_ini = dmrs.at(i).split("//")[1].split(" ")[0].toUInt();
                    uint pos_dwt_fin = dmrs.at(i).split("//")[1].split(" ")[1].toUInt();

                    // encabezado de las características por fichero dentro de la zona dmr
                    s << " sample dwt_value ratio C_positions cov_min cov_mid cov_max sites_C sites_noC sites_mC sites_hmC dist_min dist_mid dist_max\n";

                    // información del dmr para obtener las características de cada muestra
                    linea     = dmrs.at(i);
                    pos_inf   = linea.split("-")[0].toUInt();
                    pos_sup   = linea.split("-")[1].split(" ")[0].toUInt();
                    ancho_dmr = linea.split("-")[1].split(" ")[0].toUInt() - linea.split("-")[0].toUInt();

                    // guarda información de cada fichero de la zona dmr detectada
                    for (uint j = 0; j < uint(mc.size()); j++)
                    {
                        if (int(mc[j][0][11]) == 0)
                        {
                            if (visualiza_casos.at(uint(mc[j][0][10])))
                            {
                                s << " " << ficheros_case.at(int(mc[j][0][10])).split("/").back() << " ";

                                int cobertura_minima = 500000;
                                int cobertura_maxima = 0;
                                int cobertura_media  = 0;
                                int distancia_minima = 500000;
                                int distancia_maxima = 0;
                                int distancia_media  = 0;
                                int sites_C          = 0;
                                int sites_nC         = 0;
                                int sites_mC         = 0;
                                int sites_hmC        = 0;
                                int posiciones       = 0;
                                float dwt_valor      = 0.0;
                                float ratio_medio    = 0.0;

                                linea_detail.clear();

                                // busca la posición inical
                                uint posicion = 0;
                                while (pos_inf > mc[j][posicion][0])
                                    posicion++;

                                // búsqueda de valores a lo largo del DMR
                                while (pos_sup > mc[j][posicion][0] && posicion < mc[j].size())
                                {
                                    // cobertura
                                    if (cobertura_minima >= mc[j][posicion][ui->mC->isChecked() ? 2 : 8])
                                        cobertura_minima = int(mc[j][posicion][ui->mC->isChecked() ? 2 : 8]);
                                    if (cobertura_maxima < mc[j][posicion][ui->mC->isChecked() ? 2 : 8])
                                        cobertura_maxima = int(mc[j][posicion][ui->mC->isChecked() ? 2 : 8]);
                                    cobertura_media += mc[j][posicion][ui->mC->isChecked() ? 2 : 8];
                                    if (mc[j][posicion][ui->mC->isChecked() ? 2 : 8] > 0)
                                        ratio_medio += float(mc[j][posicion][ui->mC->isChecked() ? 1 : 7]);


                                    // distancia
                                    if (posicion + 2 < mc[j].size() && ancho_dmr > mc[j][posicion + 1][0] - mc[j][posicion][0])
                                    {
                                        if (distancia_minima >= mc[j][posicion + 1][0] - mc[j][posicion][0])
                                            distancia_minima = int(mc[j][posicion + 1][0] - mc[j][posicion][0]);
                                        if (distancia_maxima < mc[j][posicion + 1][0] - mc[j][posicion][0])
                                            distancia_maxima = int(mc[j][posicion + 1][0] - mc[j][posicion][0]);
                                        distancia_media += mc[j][posicion + 1][0] - mc[j][posicion][0];
                                    }

                                    // número de posiciones detectadas por tipo de mononucleótico
                                    sites_C   += (mc[j][posicion][3] > 0) ? 1 : 0;
                                    sites_nC  += (mc[j][posicion][4] > 0) ? 1 : 0;
                                    sites_mC  += (mc[j][posicion][5] > 0) ? 1 : 0;
                                    sites_hmC += (mc[j][posicion][6] > 0) ? 1 : 0;

                                    // número de posiciones detectadas con algún tipo de nucleótido sensible
                                    if (mc[j][posicion][ui->mC->isChecked() ? 2 : 8] > 0)
                                        posiciones++;

                                    posicion++;
                                }

                                // valor medio dwt en la región identificada
                                for (uint i = pos_dwt_ini; i <= pos_dwt_fin; i++)
                                    dwt_valor += h_haar_C[j][i];                    // (h_haar_C[f][i] / (pos_dwt_fin - pos_dwt_ini + 1));
                                dwt_valor /= (pos_dwt_fin - pos_dwt_ini + 1);

                                // carga de resultado en línea de texto para mostrar
                                if (cobertura_maxima >= ui->cobertura->value())
                                {
                                    s << QString("%1").arg(double(dwt_valor)) << " " << //::number(double(dwt_valor), 'f', 3) << " " <<
                                         QString("%1").arg(posiciones > 1 ? double(ratio_medio / posiciones) : double(ratio_medio)) << " " <<
                                         QString::number(posiciones) << " " <<
                                         QString::number(cobertura_minima >= 500000 ? 0 : cobertura_minima) << " " <<
                                         QString::number((posiciones > 1) ? cobertura_media / posiciones : cobertura_media) << " " <<
                                         QString::number(cobertura_maxima) << " " <<
                                         QString::number(sites_C) << " " <<
                                         QString::number(sites_nC) << " " <<
                                         QString::number(sites_mC) << " " <<
                                         QString::number(sites_hmC) << " " <<
                                         QString::number(distancia_minima >= 500000 ? 0 : distancia_minima) << " " <<
                                         QString::number((posiciones > 1) ? distancia_media / posiciones : distancia_media) << " " <<
                                         QString::number(distancia_maxima) << "\n";
                                }
                                else
                                {
                                    s << QString("%1").arg(double(dwt_valor)) << " " << //::number(double(dwt_valor), 'f', 3) << " " <<
                                         QString("%1").arg(posiciones > 1 ? double(ratio_medio / posiciones) : double(ratio_medio)) << " " <<
                                         "0 0 0 0 0 0 0 0 0 0 0\n";
                                }
                            }
                        }
                    }

                    for (uint j = 0; j < uint(mc.size()); j++)
                    {
                        if (int(mc[j][0][11]) == 1)
                        {
                            if (visualiza_control.at(uint(mc[j][0][10])))
                            {
                                s << " " << ficheros_control.at(int(mc[j][0][10])).split("/").back() << " ";

                                int cobertura_minima = 500000;
                                int cobertura_maxima = 0;
                                int cobertura_media  = 0;
                                int distancia_minima = 500000;
                                int distancia_maxima = 0;
                                int distancia_media  = 0;
                                int sites_C          = 0;
                                int sites_nC         = 0;
                                int sites_mC         = 0;
                                int sites_hmC        = 0;
                                int posiciones       = 0;
                                float dwt_valor      = 0.0;
                                float ratio_medio    = 0.0;

                                linea_detail.clear();

                                // busca la posición inical
                                uint posicion = 0;
                                while (pos_inf > mc[j][posicion][0])
                                    posicion++;

                                // búsqueda de valores a lo largo del DMR
                                while (pos_sup > mc[j][posicion][0] && posicion < mc[j].size())
                                {
                                    // cobertura
                                    if (cobertura_minima >= mc[j][posicion][ui->mC->isChecked() ? 2 : 8])
                                        cobertura_minima = int(mc[j][posicion][ui->mC->isChecked() ? 2 : 8]);
                                    if (cobertura_maxima < mc[j][posicion][ui->mC->isChecked() ? 2 : 8])
                                        cobertura_maxima = int(mc[j][posicion][ui->mC->isChecked() ? 2 : 8]);
                                    cobertura_media += mc[j][posicion][ui->mC->isChecked() ? 2 : 8];
                                    if (mc[j][posicion][ui->mC->isChecked() ? 2 : 8] > 0)
                                        ratio_medio += float(mc[j][posicion][ui->mC->isChecked() ? 1 : 7]);


                                    // distancia
                                    if (posicion + 2 < mc[j].size() && ancho_dmr > mc[j][posicion + 1][0] - mc[j][posicion][0])
                                    {
                                        if (distancia_minima >= mc[j][posicion + 1][0] - mc[j][posicion][0])
                                            distancia_minima = int(mc[j][posicion + 1][0] - mc[j][posicion][0]);
                                        if (distancia_maxima < mc[j][posicion + 1][0] - mc[j][posicion][0])
                                            distancia_maxima = int(mc[j][posicion + 1][0] - mc[j][posicion][0]);
                                        distancia_media += mc[j][posicion + 1][0] - mc[j][posicion][0];
                                    }

                                    // número de posiciones detectadas por tipo de mononucleótico
                                    sites_C   += (mc[j][posicion][3] > 0) ? 1 : 0;
                                    sites_nC  += (mc[j][posicion][4] > 0) ? 1 : 0;
                                    sites_mC  += (mc[j][posicion][5] > 0) ? 1 : 0;
                                    sites_hmC += (mc[j][posicion][6] > 0) ? 1 : 0;

                                    // número de posiciones detectadas con algún tipo de nucleótido sensible
                                    if (mc[j][posicion][ui->mC->isChecked() ? 2 : 8] > 0)
                                        posiciones++;

                                    posicion++;
                                }

                                // valor medio dwt en la región identificada
                                for (uint i = pos_dwt_ini; i <= pos_dwt_fin; i++)
                                    dwt_valor += h_haar_C[j][i];                    // (h_haar_C[f][i] / (pos_dwt_fin - pos_dwt_ini + 1));
                                dwt_valor /= (pos_dwt_fin - pos_dwt_ini + 1);

                                // carga de resultado en línea de texto para mostrar
                                if (cobertura_maxima >= ui->cobertura->value())
                                {
                                    s << QString("%1").arg(double(dwt_valor)) << " " << //::number(double(dwt_valor), 'f', 3) << " " <<
                                         QString("%1").arg(posiciones > 1 ? double(ratio_medio / posiciones) : double(ratio_medio)) << " " <<
                                         QString::number(posiciones) << " " <<
                                         QString::number(cobertura_minima >= 500000 ? 0 : cobertura_minima) << " " <<
                                         QString::number((posiciones > 1) ? cobertura_media / posiciones : cobertura_media) << " " <<
                                         QString::number(cobertura_maxima) << " " <<
                                         QString::number(sites_C) << " " <<
                                         QString::number(sites_nC) << " " <<
                                         QString::number(sites_mC) << " " <<
                                         QString::number(sites_hmC) << " " <<
                                         QString::number(distancia_minima >= 500000 ? 0 : distancia_minima) << " " <<
                                         QString::number((posiciones > 1) ? distancia_media / posiciones : distancia_media) << " " <<
                                         QString::number(distancia_maxima) << "\n";
                                }
                                else
                                {
                                    s << QString("%1").arg(double(dwt_valor)) << " " << //::number(double(dwt_valor), 'f', 3) << " " <<
                                         QString("%1").arg(posiciones > 1 ? double(ratio_medio / posiciones) : double(ratio_medio)) << " " <<
                                         "0 0 0 0 0 0 0 0 0 0 0\n";
                                }
                            }
                        }
                    }

                    s << "\n";
                }
            }
            else
                s << "no DMRs were found\n";

            data.close();
            qDebug() << "cerrando fichero";
            dmr_listo = true;
        }
    }

    qDebug() << fichero;
}

// HILO LECTURA DE DATOS DE FICHEROS
// ************************************************************************************************
void HPG_Dhunter::fichero_leido(int sample, int chrom, int inicio, int final)
{

    // paso de los datos leídos al hilo de procesamiento de datos
    qDebug() << "fichero - datos disponibles: " << contador << sample << " " << chrom;

    mutex.lock();
    if (limite_inferior > uint(inicio))
        limite_inferior = uint(inicio);
//    mutex.unlock();
//    mutex.lock();
    if (limite_superior < uint(final))
        limite_superior = uint(final);
    mutex.unlock();

    // contador de evolución de lectura y análisis
    mutex.lock();
    contador ++;
    mutex.unlock();

    if (contador % ((ficheros_case.size() + ficheros_control.size())) == 0) // * (_forward + _reverse)) == 0)
        cromosoma_leido(chrom);

}

// ************************************************************************************************
void HPG_Dhunter::cromosoma_leido(int chrom)
{
    QEventLoop loop;

    QTimer::singleShot(500, &loop, SLOT(quit()));
    loop.exec();


    // rellena los vectores de visualización con '1' en todas las posiciones
    visualiza_casos.assign(uint(ficheros_case.size()), 1);
    visualiza_control.assign(uint(ficheros_control.size()), 1);

    // ordena mc para que se gestione mejor en el programa
    // ..primero los casos ordenados según la lista
    for (uint m = 0; m < uint(ficheros_case.size()); m++)
        for (uint i = 0; i < mc_aux.size(); i++)
            if (uint(mc_aux[i][0][10]) == m && uint(mc_aux[i][0][11]) == 0)
                    mc.push_back(mc_aux[i]);
    // ..después los controles ordenados según la lista
    for (uint m = 0; m < uint(ficheros_control.size()); m++)
        for (uint i = 0; i < mc_aux.size(); i++)
            if (uint(mc_aux[i][0][10]) == m && uint(mc_aux[i][0][11]) == 1)
                    mc.push_back(mc_aux[i]);

/*
    // AJENO A APLICACION -------------------------------------------------------
    // graba un fichero en formato adecuado para pruebas con BSmooth y DSS-single
    // con formato: chr pos cov1 cov2 ... cov6 met1 met2 ... met6
    {
        qDebug() << "------------------generando fichero BSmooth------------";

        vector<uint>last_pos(6,0);
        vector<bool>size_ok(6,true);
        vector<bool>pos_coincidente(6,0);
        uint suma = 0;

        QFile data;
        data.setFileName("/home/lifercor/Descargas/toBSmoothTest_" + QString::number(cromosoma) + ".csv");
        data.open(QIODevice::WriteOnly);
        QTextStream s(&data);

        for (uint i = 0; i < 6; i++)
            suma += size_ok[i];

        while (suma)
        {
            uint pos_chr = 100000000;
            pos_coincidente.assign(6,0);

            for (uint m = 0; m < mc.size(); m++)
                if (size_ok[m])
                    if (pos_chr > mc[m][last_pos[m]][0])
                        pos_chr = uint(mc[m][last_pos[m]][0]);

            for (uint m = 0; m < mc.size(); m++)
                if (size_ok[m])
                    if (pos_chr == uint(mc[m][last_pos[m]][0]))
                    {
                        pos_coincidente[m] = 1;
                        if (last_pos[m] + 1 < mc[m].size())
                            last_pos[m]++;
                        else
                            size_ok[m] = false;
                    }

            // escribe cromosoma y posicion
            s << cromosoma << "," << pos_chr;

            // escribe cobertura
            for (uint m = 0; m < mc.size(); m++)
                if (size_ok[m])
                    s << "," << (pos_coincidente[m] ? mc[m][last_pos[m]][2] : 0);
                else
                    s << "," << 0;

            // escribe metilación
            for (uint m = 0; m < mc.size(); m++)
                if (size_ok[m])
                    s << "," << (pos_coincidente[m] ? mc[m][last_pos[m]][5] : 0);
                else
                    s << "," << 0;

            s << "\n";

            suma = 0;
            for (uint i = 0; i < 6; i++)
                suma += size_ok[i];
        }

        data.close();
        qDebug() << "------------------fichero BSmooth generado------------";
    }
    // --------------------------------------------------------------------------
    // --------------------------------------------------------------------------
*/

    qDebug() << "tamaño de las matrices de datos" << mc_aux.size() << mc.size();

    ui->statusBar->showMessage("chromosoma " + QString::number(chrom) + " read. Freeing temporary memory... please, wait a bit");

    // elimina la matriz auxiliar
    for (vector<vector<double>> n : mc_aux)
    {
        for (vector<double> m : n)
        {
                m.clear();
                m.shrink_to_fit();
                vector<double>().swap(m);
        }
        n.clear();
        n.shrink_to_fit();
        vector<vector<double>>().swap(n);
    }
    mc_aux.clear();
    mc_aux.shrink_to_fit();
    vector<vector<vector<double>>>().swap(mc_aux);

    // informa de carga terminada
    ui->statusBar->showMessage("files loaded at RAM memory. Please, select MIN COVERAGE and press ANALIZE SAMPLES");

    // habilita analizar samples
    ui->load_files->setEnabled(false);
    ui->analiza->setEnabled(true);
    ui->cobertura->setEnabled(true);
    ui->min_coverage->setEnabled(true);

    // incialización de variables y liberación de memoria
    // --------------------------------------------------------------------------------------------
    cuda_data.h_haar_L.clear();     // vector con número de datos por nivel
    cuda_data.pitch          = 0;   // ajuste óptimo de memoria GPU para datos de cada muestra
    cuda_data.pitch_2        = 0;   // ajuste óptimo de memoria GPU para auxiliar
    cuda_data.sample_num     = 0;   // número de datos por muestra
    cuda_data.samples        = 0;   // número de muestras a trasnformar
    cuda_data.levels         = 0;   // número de niveles a transformar
    cuda_data.data_adjust    = 0;   // ajuste desfase en división por nivel para número impar de datos
    cuda_data.rango_inferior = 0;   // límite inferior ventana de datos a transformar
    cuda_data.rango_superior = 1;

    // libera la memoria de la GPU
    cuda_end(cuda_data);
}

// ************************************************************************************************
void HPG_Dhunter::refGen_worker_acabado(ulong num_genex)
{
    num_genes = num_genex;

    ui->statusBar->showMessage("chromosome: " + parametros[2] +
                               ", with " + QString::number(num_genes) + " known genes");
}


// ************************************************************************************************
// ****************CONTROL DE SELECCION DE CROMOSOMA***********************************************
// ************************************************************************************************
void HPG_Dhunter::on_chr01_clicked()
{
    if (ui->chr01->isChecked())
        cromosoma_grid = 1;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr02_clicked()
{
    if (ui->chr02->isChecked())
        cromosoma_grid = 2;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr03_clicked()
{
    if (ui->chr03->isChecked())
        cromosoma_grid = 3;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr04_clicked()
{
    if (ui->chr04->isChecked())
        cromosoma_grid = 4;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr05_clicked()
{
    if (ui->chr05->isChecked())
        cromosoma_grid = 5;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr06_clicked()
{
    if (ui->chr06->isChecked())
        cromosoma_grid = 6;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr07_clicked()
{
    if (ui->chr07->isChecked())
        cromosoma_grid = 7;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr08_clicked()
{
    if (ui->chr08->isChecked())
        cromosoma_grid = 8;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr09_clicked()
{
    if (ui->chr09->isChecked())
        cromosoma_grid = 9;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr10_clicked()
{
    if (ui->chr10->isChecked())
        cromosoma_grid = 10;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr11_clicked()
{
    if (ui->chr11->isChecked())
        cromosoma_grid = 11;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr12_clicked()
{
    if (ui->chr12->isChecked())
        cromosoma_grid = 12;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr13_clicked()
{
    if (ui->chr13->isChecked())
        cromosoma_grid = 13;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr14_clicked()
{
    if (ui->chr14->isChecked())
        cromosoma_grid = 14;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr15_clicked()
{
    if (ui->chr15->isChecked())
        cromosoma_grid = 15;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr16_clicked()
{
    if (ui->chr16->isChecked())
        cromosoma_grid = 16;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr17_clicked()
{
    if (ui->chr17->isChecked())
        cromosoma_grid = 17;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr18_clicked()
{
    if (ui->chr18->isChecked())
        cromosoma_grid = 18;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr19_clicked()
{
    if (ui->chr19->isChecked())
        cromosoma_grid = 19;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr20_clicked()
{
    if (ui->chr20->isChecked())
        cromosoma_grid = 20;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr21_clicked()
{
    if (ui->chr21->isChecked())
        cromosoma_grid = 21;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr22_clicked()
{
    if (ui->chr22->isChecked())
        cromosoma_grid = 22;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr23_clicked()
{
    if (ui->chr23->isChecked())
        cromosoma_grid = 23;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr24_clicked()
{
    if (ui->chr24->isChecked())
        cromosoma_grid = 24;

    ui->load_files->setEnabled(true);
}
void HPG_Dhunter::on_chr25_clicked()
{
    if (ui->chr25->isChecked())
        cromosoma_grid = 25;

    ui->load_files->setEnabled(true);
}


// ************************************************************************************************
void HPG_Dhunter::on_cobertura_sliderMoved(int position)
{
    ui->min_coverage->setText(QString::number(position));
}

// ************************************************************************************************
void HPG_Dhunter::on_min_coverage_textEdited(const QString &arg1)
{
    ui->cobertura->setValue(arg1.toInt());
}
