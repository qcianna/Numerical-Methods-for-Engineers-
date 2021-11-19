package pl.AniaJava;

import static java.lang.Math.pow;

public class RK2 extends Data implements Function{

    public double[] getFunction(double xn, double vn, double dT, double alfa){
        double k1x = vn;
        double k1v = alfa*(1-pow(xn,2))*vn - xn;

        double k2x = vn + dT*k1v;
        double k2v = alfa*(1 - pow((xn+dT*k1x),2))*(vn + dT*k1v) - (xn + dT*k1x);

        double xn1 = xn + dT/2*(k1x + k2x);
        double vn1 = vn + dT/2*(k1v + k2v);

        double [] res = {xn1, vn1};
        return res;
    }
}
