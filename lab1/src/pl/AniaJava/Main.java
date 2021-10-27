import java.io.File;
import java.io.PrintWriter;
import java.io.FileNotFoundException;
import java.lang.Math;

public class Main {

    static double y0 = 1;
    static double lambda = -1;
    static double tMin = 0;
    static double tMax = 5;
    static double[] deltaT = {0.01, 0.1, 1};

    static double L = 0.1;
    static double C = 0.001;
    static double omega = 1/(Math.sqrt(L*C));
    static double[] omegaTab = { 0.5*omega, 0.8*omega, 1.0*omega, 1.2*omega };

    public static void check(){

        new File("check.txt");
        try{
            PrintWriter data = new PrintWriter("check.txt");
            for(int i=0; i<deltaT.length; i++){
                for(double t=tMin; t<tMax; t+=deltaT[i]){
                    data.println(t + "\t" + Math.exp(-t));
                }
                data.println();
            }
            data.close();
        } catch (FileNotFoundException e) {
            System.out.println("Error: File Not Found");
        }
    }

    public static void euler(){

        double y = y0;
        double yTemp = 0;

        new File("euler.txt");
        try{
            PrintWriter data = new PrintWriter("euler.txt");
            for(int i=0; i<deltaT.length; i++){
                for(double t=tMin; t<tMax; t+=deltaT[i]){
                    data.println(t + "\t" + y);
                    yTemp = y;
                    y += deltaT[i]*lambda*yTemp;
                }
                data.println();
                y=y0;
            }
            data.close();
        } catch (FileNotFoundException e) {
            System.out.println("Error: File Not Found");
        }
    }

    public static void rk2(){

        double k1;
        double k2;
        double y = y0;

        new File("rk2.txt");
        try{
            PrintWriter data = new PrintWriter("rk2.txt");
            for(int i=0; i<deltaT.length; i++){
                for(double t=tMin; t<tMax; t+=deltaT[i]){
                    data.println(t + "\t" + y);
                    k1 = lambda * y;
                    k2 = lambda * (y+deltaT[i]*k1);
                    y += (deltaT[i]/2)*(k1+k2);
                }
                data.println();
                y = y0;
            }
            data.close();
        } catch (FileNotFoundException e) {
            System.out.println("Error: File Not Found");
        }
    }

    public static void rk4(){

        double k1,k2,k3,k4;
        double y = y0;

        new File("rk4.txt");
        try{
            PrintWriter data = new PrintWriter("rk4.txt");
            for(int i=0; i<deltaT.length; i++){
                for(double t=tMin; t<tMax; t+=deltaT[i]){
                    data.println(t + "\t" + y);
                    k1 = lambda * y;
                    k2 = lambda * (y+(deltaT[i]/2)*k1);
                    k3 = lambda * (y+(deltaT[i]/2)*k2);
                    k4 = lambda * (y+deltaT[i]*k3);
                    y += (deltaT[i]/6) * (k1+2*k2+2*k3+k4);
                }
                data.println();
                y = y0;
            }
            data.close();
        } catch (FileNotFoundException e) {
            System.out.println("Error: File Not Found");
        }
    }

    public static double V(double n, double t, int i){
        return 10 * Math.sin(omegaTab[i] * (t+n));
    }

    public static void rrz2(){

        double dT = 1E-4;
        double R = 100;
        double T0 = 2 * Math.PI/omega;
        double Q = 0;
        double I = 0;

        double k1q, k2q, k3q, k4q;
        double k1i, k2i, k3i, k4i;

        double min = 0;
        double max = 4*T0;

        new File("rrz2.txt");
        try{
            PrintWriter data = new PrintWriter("rrz2.txt");
            for(int i=0; i<omegaTab.length; i++){
                for(double t=min; t<max; t+=dT){
                    data.println(t + "\t\t" + I + "\t" + Q);
                    k1q = I;
                    k1i = V(0, t, i)/L - Q/(L*C) - (R*I)/L;
                    k2q = I + (dT*k1i)/2;
                    k2i = V(0.5, t, i)/L - (Q+dT*k1q/2)/(L*C) - R/L*(I+dT*k1i/2);
                    k3q = I + dT*k2i/2;
                    k3i = V(0.5, t, i)/L - (Q+dT*k2q/2)/(L*C) - R/L*(I+dT*k2i/2);
                    k4q = I + dT*k3i;
                    k4i = V(1, t, i)/L - (Q+dT*k3q)/(L*C) - R/L*(I+dT*k3i);
                    Q += dT/6 * (k1q + 2*k2q + 2*k3q + k4q);
                    I += dT/6 * (k1i + 2*k2i + 2*k3i + k4i);
                }
                data.println();
                Q = 0;
                I = 0;
            }
            data.close();
        } catch (FileNotFoundException e) {
            System.out.println("Error: File Not Found");
        }
    }

    public static void main(String[] args) {
        check();
        euler();
        rk2();
        rk4();
        rrz2();
    }
}