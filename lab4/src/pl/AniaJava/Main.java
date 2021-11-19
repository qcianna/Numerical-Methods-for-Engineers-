package pl.AniaJava;

import java.io.FileNotFoundException;
import java.io.PrintWriter;

import static java.lang.Math.*;

public class Main {

    static private final double eps = 1;
    static private final double delta = 0.1;
    static private final int nx = 150;
    static private final int ny = 100;
    static private final double V1 = 10;
    static private final double V2 = 0;
    static private final double xMax = delta * nx;
    static private final double yMax = delta * ny;
    static private final double dX = 0.1 * xMax;
    static private final double dY = 0.1 * yMax;
    static private final double TOL = pow(10, -8);

    public static double p(double i, double j){
        double x = i*delta;
        double y = j*delta;
        double p1 = exp(-(pow((x - 0.35*xMax), 2)/pow(dX, 2)) - (pow((y - 0.5*yMax), 2)/pow(dY, 2)));
        double p2 = -exp(-(pow((x - 0.65*xMax), 2)/pow(dX, 2)) - (pow((y - 0.5*yMax), 2)/pow(dY, 2)));
        return p1 + p2;
    }

    public static void warunekPoczatkowy(double[][] V){
        for (int i = 0; i <= nx; i++) {
            V[i][0] = V1;
            V[i][ny] = V2;
        }
    }

    public static double warunekStopu(double[][] V){
        double S = 0;
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                S += pow(delta, 2) * (0.5 * pow(((V[i+1][j] - V[i][j])/delta), 2) + 0.5 * pow(((V[i][j+1] - V[i][j])/delta), 2) - p(i, j) * V[i][j]);
            }
        }
        return S;
    }

    public static void relaksacjaGlobalna(double wG){

        double[][] Vn = new double[nx+1][ny+1];
        double[][] Vs = new double[nx+1][ny+1];
        double S = 0;
        double Stemp;
        int iter = 0;

        warunekPoczatkowy(Vs);
        warunekPoczatkowy(Vn);

        try {
            PrintWriter dataIt = new PrintWriter("glob_iter.txt");
            PrintWriter dataS = new PrintWriter("glob_S.txt");
            PrintWriter dataVXY = new PrintWriter("glob_vxy.txt");
            PrintWriter dataErr = new PrintWriter("glob_err.txt");
            do {
                for (int i = 1; i < nx; i++) {
                    for (int j = 1; j < ny; j++) {
                        Vn[i][j] = 0.25 * (Vs[i + 1][j] + Vs[i - 1][j] + Vs[i][j + 1] + Vs[i][j - 1] + pow(delta, 2) / eps * p(i, j));
                    }
                }

                for (int j = 0; j < ny; j++) {
                    Vn[0][j] = Vn[1][j];
                    Vn[nx][j] = Vn[nx - 1][j];
                }

                for (int i = 0; i <= nx; i++) {
                    for (int j = 1; j < ny; j++) {
                        Vs[i][j] = (1.0 - wG) * Vs[i][j] + wG * Vn[i][j];
                    }
                }
                Stemp = S;
                S = warunekStopu(Vs);
                iter++;
                dataIt.println(iter);
                dataS.println(S);
            }while(abs((S - Stemp)/Stemp) > TOL);
            for(int i=0; i<=nx; i++){
                for(int j=0; j<=ny; j++){
                    dataVXY.println(i*delta + " " + j*delta + " " + Vs[i][j]);
                }
            }
            double err;
            for(int i = 1; i < nx; i++){
                for(int j = 1; j < ny; j++) {
                    err = ((Vs[i+1][j]-2*Vs[i][j]+Vs[i-1][j])/(pow(delta,2))+(Vs[i][j+1]-2*Vs[i][j]+Vs[i][j-1])/(pow(delta,2)))+(p(i,j)/eps);
                    dataErr.println(i*delta + " " + j*delta + " " + err);
                }
            }
            dataIt.close();
            dataS.close();
            dataVXY.close();
            dataErr.close();
        } catch (FileNotFoundException e) {
            System.out.println("Error: File Not Found");
        }
    }

    public static void relaksacjaLokalna(double wL) {

        double[][] V = new double[nx + 1][ny + 1];
        double S = 0;
        double Stemp;
        int iter = 0;

        warunekPoczatkowy(V);

        try {
            PrintWriter dataIt = new PrintWriter("lok_iter.txt");
            PrintWriter dataS = new PrintWriter("lok_S.txt");
            do {
                for (int i = 1; i < nx; i++) {
                    for (int j = 1; j < ny; j++) {
                        V[i][j] = (1 - wL) * V[i][j] + (wL / 4.0) * (V[i + 1][j] + V[i - 1][j] + V[i][j + 1] + V[i][j - 1] + pow(delta, 2.0) / eps * p(i, j));
                    }
                }

                for (int j = 1; j < ny; j++) {
                    V[0][j] = V[1][j];
                    V[nx][j] = V[nx - 1][j];
                }

                Stemp = S;
                S = warunekStopu(V);
                iter++;
                dataIt.println(iter);
                dataS.println(S);
            } while (abs((S - Stemp) / Stemp) >= TOL);
            dataIt.close();
            dataS.close();
        } catch (FileNotFoundException e) {
            System.out.println("Error: File Not Found");
        }
    }

    public static void main(String[] args) throws FileNotFoundException {

        double[] wG = {0.6, 1.0};
        relaksacjaGlobalna(wG[0]);

        //  double[] wL = {1.0, 1.4, 1.8, 1.9};
        //  relaksacjaLokalna(wL[0]);
    }
}
