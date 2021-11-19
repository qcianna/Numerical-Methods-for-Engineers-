package pl.AniaJava;

import static java.lang.Math.*;

public class MetodaTrapezow extends Data implements Function{

    public double[] getFunction(double xn, double vn, double dT, double alfa){
        double xn1 = xn;
        double vn1 = vn;
        double delta = pow(10,-10);
        double dx, dv;
        double a11, a12, a21, a22;
        double F, G;
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
}
