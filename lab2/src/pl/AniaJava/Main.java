package pl.AniaJava;

import java.io.PrintWriter;
import java.io.FileNotFoundException;

public class Main {

    final static private double beta = 0.001;
    final static private int N = 500;
    final static private double gamma = 0.1;
    final static private double tMax = 100;
    final static private double dT = 0.1;
    final static private double u0 = 1; //poczatkowa ilosc osob zarazonych
    final static private double TOL = 10e-6;
    final static private int maxIter = 20;
    final static private double alfa = beta*N - gamma;

    public static double f(double u){
        return alfa * u - beta * Math.pow(u, 2);
    }

    public static void metodaTrapezowPicarda(){
        double uTemp = u0;
        double u = uTemp;
        int iter;
        double um;

        try {
            PrintWriter data = new PrintWriter("picard.txt");
            for (double t = 0; t < tMax; t += dT) {
                uTemp = u;
                iter = 0;
                um = 0;
                data.println(N-u);
                while (iter <= maxIter && Math.abs(u - um) >= TOL) {
                    um = u;
                    u = uTemp + dT/2 * (f(uTemp) + f(um));
                    iter++;
                }
            }
            data.close();
        } catch (FileNotFoundException e) {
            System.out.println("Error: File Not Found");
        }
    }

    public static void metodaTrapezowNewtona(){
        double uTemp = u0;
        double u = uTemp;
        int iter;
        double um;
        try {
            PrintWriter data = new PrintWriter("newton.txt");
            for (double t = 0; t < tMax; t += dT) {
                uTemp = u;
                iter = 0;
                um = 0;
                data.println(N-u);
                while (iter <= maxIter && Math.abs(u - um) >= TOL) {
                    um = u;
                    u = um - (um - uTemp - dT/2 * ((alfa*uTemp - beta*Math.pow(uTemp, 2)) + (alfa*um - beta*Math.pow(um, 2)))
                            / (1 - dT/2 * (alfa - 2*beta*um)));
                    iter++;
                }
            }
            data.close();
        } catch (FileNotFoundException e) {
            System.out.println("Error: File Not Found");
        }
    }

    public static void RK2(){

        double a11 = 0.25;
        double a12 = 0.25 - Math.sqrt(3)/6;
        double a21 = 0.25 + Math.sqrt(3)/6;
        double a22 = 0.25;

        double b1 = 0.5;
        double b2 = 0.5;

        double U1, U2;
        double Um1, Um2;

        double uTemp = u0;
        double u = uTemp;
        int iter;

        double F1, F2;
        double dU1 = 0;
        double dU2 = 0;

        double m11, m12, m21, m22;

        try {
            PrintWriter data = new PrintWriter("rk2.txt");
            for (double t = 0; t < tMax; t += dT) {
                uTemp = u;
                iter = 0;
                U1 = u;
                U2 = u;
                data.println(N-u);
                while (iter <= maxIter && Math.min(dU1, dU2)>=TOL) { //tutaj dodac warunek stopu
                    Um1 = U1;
                    Um2 = U2;

                    F1 = U1 - uTemp - dT*(a11*f(U1) + a12*f(U2));
                    F2 = U2 - uTemp - dT*(a21*f(U1) + a22*f(U2));

                    m11 = 1 - dT*a11*(alfa - 2*beta*U1);
                    m12 = -dT*a12*(alfa - 2*beta*U2);
                    m21 = -dT*a21*(alfa - 2*beta*U1);
                    m22 = 1 - dT*a22*(alfa - 2*beta*U2);

                    dU1 = (F2*m12 - F1*m22) / (m11*m22 - m12*m21);
                    dU2 = (F1*m21 - F2*m11) / (m11*m22 - m12*m21);

                    U1 = Um1 + dU1;
                    U2 = Um2 + dU2;

                    iter++;
                }
                u = uTemp + dT*(b1*f(U1) + b2*f(U2));
            }
            data.close();
        } catch (FileNotFoundException e) {
            System.out.println("Error: File Not Found");
        }
    }

    public static void main(String[] args) throws FileNotFoundException {
        metodaTrapezowPicarda();
        metodaTrapezowNewtona();
        RK2();
    }
}
