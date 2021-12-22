import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Scanner;

import static java.lang.Math.*;


public class Main {

    final public static int nx = 400;
    final public static int ny = 90;
    final public static int i1 = 200;
    final public static int i2 = 210;
    final public static int j1 = 50;
    final public static double delta = 0.01;
    final public static double sigma = 10 * delta;
    final public static double xa = 0.45;
    final public static double ya = 0.45;
    final public static int ITMAX = 10000; //ile to ma wynosic
    final public static int step = ITMAX/5;

    public static void setV(double[][] vx, double[][] vy, double[][] psi) {
        for(int i=1; i<nx; i++){
            for(int j=1; j<ny; j++){
                vx[i][j] = (psi[i][j+1] - psi[i][j-1]) / (2.0*delta);
                vy[i][j] = - (psi[i+1][j] - psi[i-1][j]) / (2.0*delta);
            }
        }

        for(int i=i1; i<=i2; i++){
            for(int j=0; j<=j1; j++){
                vx[i][j] = 0;
                vy[i][j] = 0;
            }
        }

        for(int i=1; i<nx; i++){
            vx[i][0] = 0.0;
            vy[i][ny] = 0.0;
        }

        for(int j=0; j<=ny; j++){
            vx[0][j] = vx[1][j];
            vx[nx][j] = vx[nx-1][j];
        }
    }

    public static double findMax(double[][] vx, double[][] vy) {

        double vMax = 0;
        for(int i=0; i<=nx; i++){
            for(int j=0; j<=ny; j++){
                if(sqrt(pow(vx[i][j], 2) + pow(vy[i][j], 2)) > vMax) {
                    vMax = sqrt(pow(vx[i][j], 2) + pow(vy[i][j], 2));
                }
            }
        }

        return vMax;
    }

    public static void setU0(double[][] u0) {

        for(int i=0; i<=nx; i++){
            for(int j=0; j<=ny; j++){
                u0[i][j] = 1.0 / (2 * PI * pow(sigma, 2)) * exp(-(pow((delta*i) - xa, 2)
                        + pow((delta*j) - ya, 2)) / (2 * pow(sigma, 2)));
            }
        }
    }

    public static void copyU(double[][] u, double[][] uCopy) {

        for(int i=0; i<=nx; i++){
            System.arraycopy(u[i], 0, uCopy[i], 0, ny + 1);
        }
    }

    public static void algorithm(double[][] u0, double[][] u1, double D, double deltaT, double[][] vx, double[][] vy) {

        double tmax = ITMAX * deltaT;

        try {
//            PrintWriter date = new PrintWriter("CXsr.dat");
            for(int it=1; it<=ITMAX; it++){

                copyU(u0, u1);

                for(int k=1; k<=20; k++){
                    for(int i=0; i<=nx; i++){
                        for(int j=1; j<ny; j++){
                            if (i >= i1 && i <= i2 && j <= j1) {
                                continue;
                            } else if(i==0 || i==nx ) {
                                int iPrev;
                                int iNext;
                                if(i == 0) {
                                    iPrev = nx;
                                    iNext = i+1;
                                } else {
                                    iPrev = i-1;
                                    iNext = 0;
                                }
                                u1[i][j] = (1.0 / (1 + ((2*D*deltaT) / pow(delta, 2)))) * (u0[i][j] - (deltaT/2.0) * vx[i][j] *
                                        (((u0[iNext][j] - u0[iPrev][j]) / (2*delta)) + (u1[iNext][j] - u1[iPrev][j]) / (2*delta)) - (deltaT / 2.0) * vy[i][j] *
                                        ((u0[i][j+1] - u0[i][j-1]) / (2*delta) + (u1[i][j+1] - u1[i][j-1])/(2*delta)) + (deltaT/2.0) *
                                        D * ((u0[iNext][j] + u0[iPrev][j] + u0[i][j+1] + u0[i][j-1] - 4*u0[i][j]) / pow(delta,2)
                                        + (u1[iNext][j] + u1[iPrev][j] + u1[i][j+1] + u1[i][j-1] ) / pow(delta,2)));
                            } else {
                                u1[i][j] = (1.0 / (1 + ((2*D*deltaT) / pow(delta, 2)))) * (u0[i][j] - (deltaT/2.0) * vx[i][j] *
                                        (((u0[i+1][j] - u0[i-1][j]) / (2.0*delta)) + (u1[i+1][j] - u1[i-1][j]) / (2.0*delta)) - (deltaT/2.0) * vy[i][j] *
                                        ((u0[i][j+1] - u0[i][j-1] ) / (2.0*delta) + (u1[i][j+1] - u1[i][j-1]) / (2.0*delta)) + (deltaT/2.0)
                                        * D * ((u0[i+1][j] + u0[i-1][j] + u0[i][j+1] + u0[i][j-1] - 4*u0[i][j]) / pow(delta,2)
                                        + (u1[i+1][j] + u1[i-1][j] + u1[i][j+1] + u1[i][j-1]) / pow(delta,2)));
                            }
                        }
                    }
            }
            copyU(u1, u0);
            double c = calculateCXsr(u0, it, deltaT)[0];
            double xSr = calculateCXsr(u0, it, deltaT)[1];
//            date.println(it*deltaT + " " + c + " " + xSr);
            if(it % step == 0){
                PrintWriter dateU = new PrintWriter("U" + it + "0.1.dat");
                for(int l=0; l<=nx; l++){
                    for(int m=0; m<=ny; m++){
                        dateU.println(l + " " + m + " " + u0[l][m]);
                    }
                    dateU.println();
                }
                dateU.close();
            }
        }
//            date.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    public static double[] calculateCXsr(double[][] u0, int it, double deltaT) {

        double c = 0;
        double xSr = 0;

        for(int i=0; i<=nx; i++){
            for(int j=0; j<=ny; j++){
                c += u0[i][j];
                xSr += (i*delta) * u0[i][j];
            }
        }
        c *= pow(delta, 2);
        xSr *= pow(delta, 2);

        return new double[]{c, xSr};
    }

    public static void main(String[] args) {

        double[][] psi = new double[nx+1][ny+1];
        double[][] vx = new double[nx+1][ny+1];
        double[][] vy = new double[nx+1][ny+1];
        double[][] u0 = new double[nx+1][ny+1];
        double[][] u1 = new double[nx+1][ny+1];

        try {
            Scanner scanner = new Scanner(new File("psi.dat"));
            while(scanner.hasNext()) {

                String line = scanner.nextLine();
                String[] tokens = line.split("\s+");
                int i = Integer.parseInt(tokens[1]);
                int j = Integer.parseInt(tokens[2]);

                psi[i][j] = Double.parseDouble(tokens[3]);
            }

            setV(vx, vy, psi);

//            PrintWriter dateVx = new PrintWriter("Vx.dat");
//            PrintWriter dateVy = new PrintWriter("Vy.dat");
//            for(int i=0; i<=nx; i++){
//                for(int j=0; j<=ny; j++){
//                    dateVx.println(i + " " + j + " " + vx[i][j]);
//                    dateVy.println(i + " " + j + " " + vy[i][j]);
//                }
//            }
//            dateVx.close();
//            dateVy.close();

            double vMax = findMax(vx, vy);
            double deltaT = delta / (4*vMax);

            setU0(u0);

//            algorithm(u0, u1, 0, deltaT, vx, vy);
            algorithm(u0, u1, 0.1, deltaT, vx, vy);



        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
}
