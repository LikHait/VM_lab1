#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "QStandardItemModel"
#include "QStandardItem"
#include <QMessageBox>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    size = 3;
    s = new SLE(3, 0.00001);
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::preparations()
{
    ui->sizeBox->setText(QString::number(size));
    ui->pushButton_3->clicked();
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size + 1; ++j)
            ui->tableWidget->setItem(i, j, new QTableWidgetItem(QString::number(0)));
}

void MainWindow::changedSize()
{
    bool ok;
    size = (ui->sizeBox->text()).toDouble(&ok);
    if (ok)
    {
        QStandardItemModel *model = new QStandardItemModel;
        QStringList horizontalHeader;
        QString str = "x";
        for (int i = 0; i < size; ++i)
            horizontalHeader.append(str + QString::number(i));
        horizontalHeader.append("val");

        model->setHorizontalHeaderLabels(horizontalHeader);


        ui->tableWidget->setRowCount(size);
        ui->tableWidget->setColumnCount(size + 1);

        ui->tableWidget->resizeRowsToContents();
        ui->tableWidget->resizeColumnsToContents();

        ui->tableWidget->setHorizontalHeaderLabels(horizontalHeader);

    }
}

void MainWindow::genMatrix()
{
    s->resize(size);
    if (ui->prevalence->isChecked())
        s->prevalenceToTrue();
    else
        s->prevalenceToFalse();

    if (ui->symmetric->isChecked())
        s->symmetricToTrue();
    else
        s->symmetricToFalse();

    s->GenerateSuitableMatrix();

    vector<vector<double>> matrix = s->GetMatrix();

    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size + 1; ++j)
            ui->tableWidget->setItem(i, j, new QTableWidgetItem(QString::number(matrix[i][j])));
}

void MainWindow::run()
{
    QStandardItemModel *model = new QStandardItemModel;
    QStandardItem *item;
    QStringList horizontalHeader;
    QStringList verticalHeader;

    QMessageBox msgBox;
    msgBox.setIcon(QMessageBox::Warning);


    int time[5];
    vector<vector<double>> matrix(size);

    for (int i = 0; i < size; ++i)
    {
        matrix[i].resize(size + 1);
        for (int j = 0; j < size + 1; ++j)
           matrix[i][j] = (ui->tableWidget->item(i, j)->text()).toDouble();
    }
    s->resize(size);
    s->ManualInput(matrix);
    s->resetResult();

    if (ui->method->isChecked())
    {
        time[0] = s->Gauss();
        verticalHeader.append("Метод Гаусса");
    }
    if (ui->method_2->isChecked())
    {
        time[1] = s->Cramer();
        verticalHeader.append("Метод Крамера");
    }
    if (ui->method_3->isChecked())
    {
        try {
        time[2] = s->FixedPointIteration();
        }

        catch(const char * str)
        {
            msgBox.setText(str);

            msgBox.exec();
        }

        verticalHeader.append("Метод простой итераций");
    }
    if (ui->method_4->isChecked())
    {
        try {
        time[3] = s->Seidel();
        }

        catch(const char * str)
        {
            msgBox.setText(str);

            msgBox.exec();
        }

        verticalHeader.append("Метод Зейделя");
    }
    if (ui->method_5->isChecked())
    {
        try {
        time[4] = s->Relaxation();
        }

        catch(const char * str)
        {
            msgBox.setText(str);

            msgBox.exec();
        }
        verticalHeader.append("Метод верхней релаксации");
    }


    QString str = "x";

    for (int i = 0; i < size; ++i)
        horizontalHeader.append(str + QString::number(i));

    //horizontalHeader.append("time");

    model->setHorizontalHeaderLabels(horizontalHeader);
    model->setVerticalHeaderLabels(verticalHeader);

    int st = 0;
    vector<vector<double>> result = s->GetResult();
    if (ui->method->isChecked())
    {
        for (int j = 0; j < size; ++j)
        {
            item = new QStandardItem(QString::number(result[GAUSS][j]));
            model->setItem(st, j, item);
        }
        //item = new QStandardItem(QString::number(time[GAUSS]));
        //model->setItem(st, size, item);
        st++;
    }
    if (ui->method_2->isChecked())
    {
        for (int j = 0; j < size; ++j)
        {
            item = new QStandardItem(QString::number(result[CRAMER][j]));
            model->setItem(st, j, item);
        }
        //item = new QStandardItem(QString::number(time[CRAMER]));
        //model->setItem(st, size, item);
        st++;
    }
    if (ui->method_3->isChecked())
    {
        for (int j = 0; j < size; ++j)
        {
            item = new QStandardItem(QString::number(result[FIXEDPOINT][j]));
            model->setItem(st, j, item);
        }
        //item = new QStandardItem(QString::number(time[FIXEDPOINT]));
        //model->setItem(st, size, item);
        st++;
    }
    if (ui->method_4->isChecked())
    {
        for (int j = 0; j < size; ++j)
        {
            item = new QStandardItem(QString::number(result[SEIDEL][j]));
            model->setItem(st, j, item);
        }
        //item = new QStandardItem(QString::number(time[SEIDEL]));
        //model->setItem(st, size, item);
        st++;
    }
    if (ui->method_5->isChecked())
    {
        for (int j = 0; j < size; ++j)
        {
            item = new QStandardItem(QString::number(result[RELAXATION][j]));
            model->setItem(st, j, item);
        }
        //item = new QStandardItem(QString::number(time[RELAXATION]));
        //model->setItem(st, size, item);
        st++;
    }

    if (st == 0)
    {
    QMessageBox msgBox;
    msgBox.setIcon(QMessageBox::Information);
    msgBox.setText("Вы не выбрали ни одного метода. Результат отсутствует.");

    msgBox.exec();
    }

    ui->tableView->setModel(model);

    ui->tableView->resizeRowsToContents();
    ui->tableView->resizeColumnsToContents();

}
