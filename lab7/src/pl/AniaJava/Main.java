package pl.AniaJava;

import java.io.PrintWriter;

import static java.lang.Math.pow;

public class Main {

    final private static double delta = 0.01;
    final private static double p = 1;
    final private static int u = 1;
    final private static int nx = 200;
    final private static int ny = 90;
    final private static int i1 = 50;
    final private static int j1 = 55;
    final private static int IT_MAX = 20000;

    public static void WB_psi(double Qwe, double Qwy, double[] y, double[][] psi) {

        //A
        for(int j = j1; j <= ny; j++){
            psi[0][j] = Qwe / (2*u) * (y[j] * y[j] * y[j] / 3 - y[j] * y[j] / 2
                    * (y[j1] + y[ny]) + y[j] * y[j1] * y[ny]);
        }

        //C
        for(int j = 0; j <= ny; j++){
            psi[nx][j] = Qwy / (2*u) * (y[j] * y[j] * y[j] / 3 - y[j] * y[j] / 2
                    * y[ny]) +  Qwe * y[j1] * y[j1] * (-y[j1] + 3 * y[ny]) / (12 * u);
        }

        //B
        for(int i = 1; i < nx; i++){
            psi[i][ny] = psi[0][ny];
        }

        //D
        for(int i = i1; i < nx; i++){
            psi[i][0] = psi[0][j1];
        }

        //E
        for(int j = 1; j <= j1; j++){
            psi[i1][j] = psi[0][j1];
        }

        //F
        for(int i = 1; i <= i1; i++){
            psi[i][j1] = psi[0][j1];
        }
    }

    public static void WB_zeta(double Qwe, double Qwy, double[] y, double[][] psi, double[][] zeta) {

        //A
        for(int j = j1; j <= ny; j++){
            zeta[0][j] = Qwe / (2*u) * (y[j] * 2 - y[j1] - y[ny]);
        }

        //C
        for(int j = 0; j <= ny; j++){
            zeta[nx][j] = Qwy / (2*u) * (2 * y[j] - y[ny]);
        }

        //B
        for(int i = 1; i < nx; i++){
            zeta[i][ny] = 2 / (pow(delta, 2)) * (psi[i][ny - 1] - psi[i][ny]);
        }

        //D
        for(int i = i1 + 1; i < nx; i++){
            zeta[i][0] = 2 / (pow(delta, 2)) * (psi[i][1] - psi[i][0]);
        }

        //E
        for(int j = 1; j < j1; j++){
            zeta[i1][j] = 2 / (pow(delta, 2)) * (psi[i1 + 1][j] - psi[i1][j]);
        }

        //F
        for(int i = 1; i <= i1; i++){
            zeta[i][j1] = 2 / (pow(delta, 2)) * (psi[i][j1 + 1] - psi[i][j1]);
        }

        //E/F
        zeta[i1][j1] = 0.5 * (zeta[i1 - 1][j1] + zeta[i1][j1 - 1]);
    }

    public static void setY(double[] y) {

        for(int j=0; j<=ny; j++){
            y[j] = delta * j;
        }
    }

    public static double calculateQwy(double Qwe, double[] y) {

        return Qwe * (pow(y[ny], 3) - pow(y[j1], 3) - 3 * pow(y[ny], 2)
                * y[j1] + 3 * pow(y[j1], 2) * y[ny]) / (pow(y[ny], 3));
    }

    public static boolean checkBox(int i, int j) {

        return i > i1 || j >= j1;
    }

    public static boolean checkBorder(int i, int j){

        if(i==0 || i == nx) return true;
        if(j==0 || j == ny) return true;
        if(i == i1 && j <=j1) return true;
        return i <= i1 && j == j1;
    }

    public static void relaksacjaNS(double Qwe, double Qwy, double[] y, double[][] psi, double[][] zeta) {

        WB_psi(Qwe, Qwy, y, psi);

        int omega;

        for (int it = 1; it <= IT_MAX; it++) {

            if (it < 2000) {
                omega = 0;
            } else {
                omega = 1;
            }

            for (int i = 1; i < nx; i++) {
                for (int j = 1; j < ny; j++) {
                    if (!checkBorder(i, j) && checkBox(i, j)) {
                        psi[i][j] = 0.25 * (psi[i + 1][j] + psi[i - 1][j] + psi[i][j + 1] + psi[i][j - 1]
                                - pow(delta, 2) * zeta[i][j]);

                        zeta[i][j] = 0.25 * (zeta[i + 1][j] + zeta[i - 1][j] + zeta[i][j + 1] + zeta[i][j - 1]) -
                                omega * p / (16 * u) * ((psi[i][j + 1] - psi[i][j - 1]) * (zeta[i + 1][j] - zeta[i - 1][j])
                                        - (psi[i + 1][j] - psi[i - 1][j]) * (zeta[i][j + 1] - zeta[i][j - 1]));
                    }
                }

                WB_zeta(Qwe, Qwy, y, psi, zeta);
            }
        }
    }

    public static void getResults(double Qwe, double Qwy, double[] y) {

        double[][] psi = new double[nx+1][ny+1];
        double[][] zeta = new double[nx+1][ny+1];

        WB_psi(Qwe, Qwy, y, psi);
        relaksacjaNS(Qwe, Qwy, y, psi, zeta);

        double u = 0;
        double v = 0;

        try {
            PrintWriter data = new PrintWriter("Qwe" + Qwe + ".dat");
            for(int i = 1; i < nx; i++) {

                for (int j = 1; j < ny; j++) {

                    if (!checkBorder(i, j) && checkBox(i, j)) {

                        u = (psi[i][j + 1] - psi[i][j - 1]) / (2 * delta);
                        v = -(psi[i + 1][j] - psi[i - 1][j]) / (2 * delta);
                    } else {

                        u = 0;
                        v = 0;
                    }

                    data.println(delta * i + " " + delta * j + " " + psi[i][j] + " " + zeta[i][j] + " " + u + " " + v);
                }
                data.println();
            }
            data.close();
        } catch (Exception e) {
            System.out.println(e);
        }
    }

    public static void main(String[] args) {

        double[] y = new double[ny+1];
        setY(y);

        double Qwe = -1000;
        double Qwy = calculateQwy(Qwe, y);
        getResults(Qwe, Qwy, y);

        Qwe = -4000;
        Qwy = calculateQwy(Qwe, y);
        getResults(Qwe, Qwy, y);

        Qwe = 4000;
        Qwy = calculateQwy(Qwe, y);
        getResults(Qwe, Qwy, y);
    }
}
