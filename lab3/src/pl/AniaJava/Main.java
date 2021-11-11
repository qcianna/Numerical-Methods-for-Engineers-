package pl.AniaJava;

import java.io.FileNotFoundException;
import java.io.PrintWriter;

import static java.lang.Math.*;

public class Main {

    final private static double x0 = 0.01;
    final private static double v0 = 0.0;
    final private static double dT0 = 1.0;
    final private static double S = 0.75;
    final private static int p = 2; //rząd dokładności metod
    final private static double tMax = 40;
    final private static double alfa = 5.0;

    public static double g(double x, double v){
        return alfa * (1 - pow(x, 2.0))*v - x;
    }

    public static double[] metodaTrapezow(double xn, double vn, double dT, double alfa){

        double xn1 = xn;
        double vn1 = vn;
        double delta = pow(10,-10);
        double dx, dv;
        double a11, a12, a21, a22;
        double F, G;
        System.out.println("\n" + xn);
        do{

            a11 = 1.0;
            a12 = -dT/2.0;
            a21 = -(dT/2.0) * (-2 * alfa * xn1 * vn1 - 1);
            a22 = 1 - (dT/2.0) * alfa * (1 - pow(xn1, 2.0));

            F = xn1 - xn - (dT/2.0) * (vn + vn1);
            G = vn1 - vn - (dT/2.0) * (g(xn, vn) + g(xn1, vn1));

            dx = ((-F)*a22 - (-G)*a12) / (a11*a22 - a12*a21);
            dv = (a11*(-G) - a21*(-F)) / (a11*a22 - a12*a21);

            xn1 += dx;
            vn1 += dv;

        }while(max(abs(dx), abs(dv)) < delta);

        double [] res = {xn1, vn1};
        return res;
    }

    public static double[] rk2(double xn, double vn, double dT, double alfa){

        double k1x = vn;
        double k1v = alfa*(1-pow(xn,2))*vn - xn;

        double k2x = vn + dT*k1v;
        double k2v = alfa*(1 - pow((xn+dT*k1x),2))*(vn + dT*k1v) - (xn + dT*k1x);

        double xn1 = xn + dT/2*(k1x + k2x);
        double vn1 = vn + dT/2*(k1v + k2v);

        double [] res = {xn1, vn1};
        return res;
    }

    public static void timeStepControl(double TOL, int met){

        double t = 0.0;
        double dT = dT0;
        double xn = x0;
        double vn = v0;
        double xn12, vn12;
        double xn22, vn22;
        double xn21, vn21;
        double Ex, Ev;

        if(met == 0) {
            try {
                PrintWriter data = new PrintWriter("rk2.txt");
                do {

                    xn12 = rk2(xn, vn, dT, alfa)[0];
                    vn12 = rk2(xn, vn, dT, alfa)[1];

                    xn22 = rk2(xn12, vn12, dT, alfa)[0];
                    vn22 = rk2(xn12, vn12, dT, alfa)[1];

                    xn21 = rk2(xn, vn, 2.0 * dT, alfa)[0];
                    vn21 = rk2(xn, vn, 2.0 * dT, alfa)[1];

                    Ex = (xn22 - xn21) / (pow(2.0, p) - 1);
                    Ev = (vn22 - vn21) / (pow(2.0, p) - 1);

                    double maxE = max(abs(Ex), abs(Ev));

                    if (maxE < TOL) {
                        t = t + 2.0 * dT;
                        xn = xn22;
                        vn = vn22;
                        data.println(vn);
                    }

                    dT = pow(((S * TOL) / maxE), (1.0 / (p + 1))) * dT;

                } while (t < tMax);
                data.close();
            } catch (FileNotFoundException e) {
                System.out.println("Error: File Not Found");
            }
        }

        else {
            try {
                PrintWriter data2 = new PrintWriter("trapez.txt");
                do {

                    xn12 = metodaTrapezow(xn, vn, dT, alfa)[0];
                    vn12 = metodaTrapezow(xn, vn, dT, alfa)[1];

                    xn22 = metodaTrapezow(xn12, vn12, dT, alfa)[0];
                    vn22 = metodaTrapezow(xn12, vn12, dT, alfa)[1];

                    xn21 = metodaTrapezow(xn, vn, 2 * dT, alfa)[0];
                    vn21 = metodaTrapezow(xn, vn, 2 * dT, alfa)[1];

                    Ex = (xn22 - xn21) / (pow(2, p) - 1);
                    Ev = (vn22 - vn21) / (pow(2, p) - 1);

                    double maxE = max(abs(Ex), abs(Ev));

                    if (maxE < TOL) {
                        t += 2 * dT;
                        xn = xn22;
                        vn = vn22;
                        data2.println(vn);
                    }

                    dT = pow(((S * TOL) / maxE), (1.0 / (p + 1))) * dT;

                } while (t < tMax);
                data2.close();
            } catch (FileNotFoundException e) {
                System.out.println("Error: File Not Found");
            }
        }
    }

    public static void main(String[] args) throws FileNotFoundException {

        timeStepControl(1e-2, 0);
        timeStepControl(1e-2, 1);

//        timeStepControl(1e-5, 0);
//        timeStepControl(1e-5, 1);
    }
}
