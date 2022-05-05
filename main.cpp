#include <iostream>
#define iteracje 50
#define TOLF 1e-10
#define TOLX 1e-10
#define n 4

double max(double *vector)
{
    double maksimum = fabs(vector[0]);
    for (int i = 1; i < n; i++) {
        if (fabs(vector[i]) > maksimum)
            maksimum = fabs(vector[i]);
    }
    return maksimum;
}

double **nowaMacierz()
{
    double **matrix;
    matrix = new double *[n];
    for (int i = 0; i < n; i++)
        matrix[i] = new double[n];
    return matrix;
}
void wypelnijMarcierz(double **matrix)
{
    double tmp[4][4] = { { 100.0, -1.0, 2.0, -3.0 },
                         { 1.0, 200.0, -4.0, 5.0 },
                         { -2.0, 4.0, 300.0, -6.0 },
                         { 3.0, -5.0, 6.0, 400.0 } };
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            matrix[i][j] = tmp[i][j];
}
void usunMacierz(double **matrix)
{
    for (int i = n - 1; i >= 0; i--)
        delete[]matrix[i];
    delete[]matrix;
}
void wypiszMacierz(double **matrix)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            std::cout << matrix[i][j] << "\t";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
double *nowyWektor()
{
    double *vector = new double[n];
    return vector;
}
void wypelnijWektor(double *v, double *vtmp)
{
    for (int i = 0; i < n; i++)
        v[i] = vtmp[i];
}
void wypiszWektor(double *vector)
{
    for (int i = 0; i < n; i++)
        std::cout << vector[i] << "\t";
    std::cout << std::endl;
}

void metodaJacobiego(double** A, double* b, double* x){
    std::cout << "\n\n";
    std::cout << "Metoda Jacobiego" << std::endl;
    double estymator[n], residuum[n];
    double est, res;
    double nowyx[n], suma;

    for(int licznik = 0; licznik < iteracje; licznik++){
        for(int i = 0; i < n; i++){
            suma = 0.0;
            for(int j = 0; j < n; j++){
                if(i != j){
                    suma += x[j] * A[i][j];
                }
            }
            nowyx[i] = (b[i] - suma) / A[i][i];
        }

        for(int i = 0; i < n; i++){
            residuum[i] = 0.0;
            estymator[i] = fabs(nowyx[i] - x[i]);

            for(int j = 0; j < n; j++){
                residuum[i] += A[i][j] * x[j];
            }
            residuum[i] = fabs(residuum[i] - b[i]);

            x[i] = nowyx[i]; //zamiana dopiero po obliczeniu estymatora i residuum
        }

        est = max(estymator);
        res = max(residuum);

        std::cout << "Licznik: " << licznik << "\t" << "Residuum: " << res << "\t" << "Estymator: " << est << std::endl;
        std::cout << "Wektor: ";
        wypiszWektor(x);
        std::cout << "\n\n";
        if(est < TOLX && res < TOLF)
            break;
    }
}

void gaussaSeidela(double** A, double* b, double* x){
    std::cout << "\n\n";
    std::cout << "Metoda Gaussa-Seidela" << std::endl;
    double estymator[n], residuum[n];
    double est, res;
    double nowyx[n], suma;

    for(int licznik = 0; licznik < iteracje; licznik++) {
        for(int i = 0; i < n; i++) {
            suma = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    suma += A[i][j] * x[j];
                }
            }
            nowyx[i] = x[i];
            x[i] = (b[i] - suma) / A[i][i];
            estymator[i] = fabs(nowyx[i] - x[i]);
        }

        for(int i = 0; i < n; i++){
            residuum[i] = 0.0;
            for(int j = 0; j < n; j++){
                residuum[i] += A[i][j] * x[j];
            }
            residuum[i] = fabs(residuum[i] - b[i]);
        }
        est = max(estymator);
        res = max(estymator);

        std::cout << "Licznik: " << licznik << "\t" << "Residuum: " << res << "\t" << "Estymator: " << est << std::endl;
        std::cout << "Wektor: ";
        wypiszWektor(x);
        std::cout << "\n\n";
        if(est < TOLX && res < TOLF)
            break;
    }

}

void sor(double** A, double* b, double* x, double omega){
    std::cout << "\n\n";
    std::cout << "Metoda SOR" << std::endl;
    static const int size = n;
    double estymator[size], residuum[size];
    double est, res;
    double nowyx[size], suma;

    for(int licznik = 0; licznik < iteracje; licznik++){
        for(int i = 0; i < n; i++) {
            suma = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    suma += A[i][j] * x[j];
                }
            }
            nowyx[i] = x[i];
            x[i] = (1.0 - omega) * x[i] + (omega * (b[i] - suma) / A[i][i]);
            estymator[i] = fabs(nowyx[i] - x[i]);
        }

        for(int i = 0; i < n; i++){
            residuum[i] = 0.0;
            for(int j =0; j < n; j++){
                residuum[i] += A[i][j] * x[j];
            }
            residuum[i] = fabs(residuum[i] - b[i]);
        }
        est = max(estymator);
        res = max(residuum);

        std::cout << "Licznik: " << licznik << "\t" << "Residuum: " << res << "\t" << "Estymator: " << est << std::endl;
        std::cout << "Wektor: ";
        wypiszWektor(x);
        std::cout << "\n\n";
        if(est < TOLX && res < TOLF)
            break;
    }
}

int main() {
    double** macierz = nowaMacierz();
    double pb[n] = {116.0, -226.0, 912.0, -1174.0};
    double px[n] = {2.0, 2.0, 2.0, 2.0};
    double* b = nowyWektor();
    double* x = nowyWektor();
    double omega = 0.5;
    wypelnijWektor(b, pb);
    wypelnijWektor(x, px);
    wypelnijMarcierz(macierz);

    std::cout << "Macierz A:" << std::endl;
    wypiszMacierz(macierz);

    std::cout << "Wektor x:" << std::endl;
    wypiszWektor(x);

    std::cout << "Wektor b:" << std::endl;
    wypiszWektor(b);

    metodaJacobiego(macierz, b, x);


    wypelnijWektor(b, pb);
    wypelnijWektor(x, px);
    wypelnijMarcierz(macierz);

    gaussaSeidela(macierz, b, x);


    wypelnijWektor(b, pb);
    wypelnijWektor(x, px);
    wypelnijMarcierz(macierz);

    sor(macierz, b, x , omega);

    usunMacierz(macierz);
}
