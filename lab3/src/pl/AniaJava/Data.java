package pl.AniaJava;

import java.io.FileNotFoundException;
import java.io.PrintWriter;

import static java.lang.Math.*;

public class Data {

    final protected static double x0 = 0.01;
    final protected static double v0 = 0.0;
    final protected static double dT0 = 1.0;
    final protected static double S = 0.75;
    final protected static int p = 2;
    final protected static double tMax = 40;
    final protected static double alfa = 5.0;

    public static double g(double x, double v) {
        return alfa * (1 - pow(x, 2.0)) * v - x;
    }

    public static void timeStepControl(double TOL, Function f) {

        double t = 0.0;
        double dT = dT0;
        double xn = x0;
        double vn = v0;
        double xn12, vn12;
        double xn22, vn22;
        double xn21, vn21;
        double Ex, Ev;

        try {
            PrintWriter data = new PrintWriter("rk2.txt");
            do {

                xn12 = f.getFunction(xn, vn, dT, alfa)[0];
                vn12 = f.getFunction(xn, vn, dT, alfa)[1];

                xn22 = f.getFunction(xn12, vn12, dT, alfa)[0];
                vn22 = f.getFunction(xn12, vn12, dT, alfa)[1];

                xn21 = f.getFunction(xn, vn, 2.0 * dT, alfa)[0];
                vn21 = f.getFunction(xn, vn, 2.0 * dT, alfa)[1];

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
}
