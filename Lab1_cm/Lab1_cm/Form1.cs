using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.IO;

namespace Lab1_cm
{
    public partial class Form1 : Form
    {
        public double[,] m;
        public Form1()
        {
            InitializeComponent();
        }
        //R = B - AX

        public double[] Relaxation(double[,] matr) // Метод Релаксации решения СЛАУ лабараторная работа 1 задние 4,4
        {
            // Разделим матрицу на две матрицы, матрицу коэф. и матрицу свободных членов
            int n = matr.GetLength(0);
            int iter = 0;
            double norma;
            double[,] A = new double[matr.GetLength(0), matr.GetLength(0)];
            double[] B = new double[matr.GetLength(0)];
            double[] previousVariableValues = new double[matr.GetLength(0)]; //previousVariableValues -  предидущее приближение R(k -1) 
            double[] currentVariableValues = new double[matr.GetLength(0)]; //currentVariableValues - текущие приближение R(k)
            
            string s = "Метод Релаксации\nПриближение:\n";
            string s2 = "";

            for (int i = 0; i < n; i++)
            {
                s += "x[" + (i + 1) + "]        ";
                for (int j = 0; j < n; j++)
                {
                    A[i, j] = matr[i, j];
                }
            }
            for (int i = 0; i < n; i++)
                B[i] = matr[i, n];
            s += "\n";

            double eps = Convert.ToDouble(textBox1.Text);
            double w = Convert.ToDouble(textBox2.Text);
            
            //Задаем нулевое начальное приближение
            for (int i = 0; i < n; i++)
            {
                currentVariableValues[i] = Convert.ToDouble(textBox3.Text);
                previousVariableValues[i] = currentVariableValues[i];
            }
            //Реализация алгоритма 
            do
            {
                iter++;
                norma = 0;
                double[] R = new double[n];
                int count = 0;
                for (int i = 0; i < n; i++) // Вывод невязок
                {
                    R[i] = B[i];
                    for (int j = 0; j < n; j++)
                    {
                        R[i] -= A[i, j] * previousVariableValues[j];
                        count++;
                    }
                    s2 += Math.Round(R[i], 7).ToString() + " || ";
                }

                for (int i = 0; i < n; i++)
                {
                    s += Math.Round(previousVariableValues[i], 5) + "   ";
                    previousVariableValues[i] = B[i];
                    for (int j = 0; j < n; j++) // Построение приближеия Ri(k - 1)
                    {
                        if (i != j)
                            previousVariableValues[i] = previousVariableValues[i] - A[i, j] * previousVariableValues[j]; // Записаь предидущего приближения для i != j
                    }

                    previousVariableValues[i] /= A[i, i]; // Приведем матрицу к удобному виду для решения
                    previousVariableValues[i] = w * previousVariableValues[i] + (1 - w) * currentVariableValues[i];
                    

                    if (Math.Abs(previousVariableValues[i] - currentVariableValues[i]) > norma) // Условие окончания алгоритма, достижение задаанной точности 
                        norma = Math.Abs(previousVariableValues[i] - currentVariableValues[i]);


                    currentVariableValues[i] = previousVariableValues[i]; // Переход к следующей Rk, переход к следующей итерации

                }
                s += "\n";
                s2 += "\n";
            } while (norma > eps);

            s += "\nИтераций " + iter.ToString() + "\n";
            for (int i = 0; i < n; i++)
            {
                i++;
                s += "x[" + i + "] = ";
                i--;
                s += previousVariableValues[i].ToString();
                s += "\n";
            }
            richTextBox2.Clear();
            richTextBox2.AppendText(s);
            richTextBox2.AppendText("\nНевязки\n" + s2);

            return previousVariableValues;
        }
        private double[,] GetMatrixFromFile(string fileName)
        {
            double[,] matrix;
            if (File.Exists(fileName) == false)
            {

                MessageBox.Show("ERROR \nНеправельное имя файла");
                return null;
            }
            using (StreamReader file = new StreamReader(fileName))
            {
                string line;
                string str = file.ReadLine();
                string[] NM = str.Split(new char[] { ' ', 'x' }, StringSplitOptions.RemoveEmptyEntries);
                int n = int.Parse(NM[0]);
                int m = int.Parse(NM[1]);
                matrix = new double[n, m];
                for (int i = 0; (i < n) && ((line = file.ReadLine()) != null); i++)
                {
                    char[] delimiterChars = { ' ', '\t' };
                    string[] numbers = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);
                    int j = 0;
                    foreach (string numString in numbers)
                    {
                        double x;
                        bool CanConvert = double.TryParse(numString, out x);
                        if (CanConvert)
                        {
                            matrix[i, j] = x;
                            j++;
                        }
                    }
                }
            }
            int k = 0;
            for(int i = 0; i < matrix.GetLength(0); i++)
                for(int j = 0; j <matrix.GetLength(1); j++ )
                {
                    if (matrix[i, j] == 0) k++;
                    if(k == matrix.GetLength(0) * matrix.GetLength(1))
                    {
                        MessageBox.Show("Расширенная матрица нулевая\nВыбирите другой файл");
                        break;
                    }
                }

            MessageBox.Show("Ввод матрицы завершен");
            m = matrix;
            return matrix;
        }
        private double GaussMethod(double[,] matr)
        {
            int n = matr.GetLength(0);
            double[,] A = new double[matr.GetLength(0), matr.GetLength(0)];
            double[] B = new double[matr.GetLength(0)];
            double s;
            double[] x = new double[n];

            for (int i = 0; i < n; i++) x[i] = 0; //Массив ответов

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    A[i, j] = matr[i, j];
                }
            } //Перезапись массивов 
            for (int i = 0; i < n; i++)
                B[i] = matr[i, n];

            double max;
            int k, index;
            double eps = 0.00001;  // точность
            k = 0;
            string l = "";
            while (k < n)
            {
                // Поиск строки с максимальным a[i][k]
                max = Math.Abs(A[k,k]);
                index = k;
                for (int i = k + 1; i < n; i++)
                {
                    if (Math.Abs(A[i, k]) > max)
                    {
                        max = Math.Abs(A[i, k]);
                        index = i;
                    }
                }

                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        l += Math.Round(A[i, j], 5).ToString() + " || ";
                    }
                    l += Math.Round(B[i], 5).ToString() + "\n";
                }
                //Перестановка строк
                if (max < eps)
                {
                    // нет ненулевых диагональных элементов
                    richTextBox2.Clear();
                    richTextBox2.AppendText("Решение получить невозможно \n");
                }

                double temp = 0;
                for (int j = 0; j < n; j++)
                {
                    temp = A[k,j];
                    A[k,j] = A[index,j];
                    A[index,j] = temp;
                }
                temp = B[k];
                B[k] = B[index];
                B[index] = temp;
                // Нормализация уравнений
                for (int i = k; i < n; i++)
                {
                    temp = A[i,k];
                    if (Math.Abs(temp) < eps) continue; // для нулевого коэффициента пропустить
                    for (int j = 0; j < n; j++)
                        A[i,j] = A[i,j] / temp;
                    B[i] = B[i] / temp;
                    if (i == k) continue; // уравнение не вычитать само из себя
                    for (int j = 0; j < n; j++)
                        A[i,j] = A[i,j] - A[k,j];
                    B[i] = B[i] - B[k];
                }
                k++;
                l += "================================\n";
            }
            //if(B[n] !=0)
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    l += Math.Round(A[i, j], 5).ToString() + " || ";
                }
                l += Math.Round(B[i], 5).ToString() + "\n";
            }
            // обратная подстановка
            for (k = n - 1; k >= 0; k--)
            {
                x[k] = B[k];
                for (int i = 0; i < k; i++)
                    B[i] = B[i] - A[i,k] * x[k];
            }
            bool flag = false;
            for (int i = 0; i < n; i++)
            {
                flag = false;
                for (int j = 0; j < n; j++)
                {
                    if (A[i, j] == 0)
                        flag = false;
                    else
                    {
                        flag = true;
                        break;
                    } 
                }
                if (!flag)
                {
                    if (B[i] == 0) break;
                    else flag = true;
                }
            }
            string m = "Метод Гаусса\n";
            if (!flag) m += "Ситсема имеет множество решений\n Одно из них\n";
            for (int i = 0; i < n; i++)
            {
                i++;
                m += "x[" + i + "] = ";
                i--;
                m += x[i].ToString();
                m += "\n";
            }//Вывод на экран 
            richTextBox2.Clear();
            richTextBox2.AppendText(m);
            richTextBox2.AppendText(l);
            return 1;
        }
        

        private void button1_Click(object sender, EventArgs e)
        {
            OpenFileDialog OPF = new OpenFileDialog();
            OPF.Filter = "Файлы txt|*.txt";
            double[,] matrix;
            string s = "";
            richTextBox1.Clear();
            if (OPF.ShowDialog() == DialogResult.OK)
            {
                matrix = GetMatrixFromFile(OPF.FileName);
                for (int i = 0; i < matrix.GetLength(0); i++, Console.WriteLine())
                {

                    for (int j = 0; j < matrix.GetLength(1); j++)
                        s += "  " + matrix[i, j] + " ";
                    s += "\n";
                }
                richTextBox1.AppendText(s);
            }
        }

        private void button2_Click(object sender, EventArgs e)
        {
            var test = Relaxation(m);
        }

        private void button3_Click(object sender, EventArgs e)
        {
            var test = GaussMethod(m);
        }
    }
}
// выводить вектор невязки 
// исправить выраждегые матрицы 
//
