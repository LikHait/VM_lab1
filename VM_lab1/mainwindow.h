#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "Methods.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    void preparations();
    ~MainWindow();
private slots:
    void changedSize();
    void genMatrix();
    void run();
private:
    SLE* s;
    int size;
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
