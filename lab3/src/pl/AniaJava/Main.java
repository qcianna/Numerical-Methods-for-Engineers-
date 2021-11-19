package pl.AniaJava;

import java.io.FileNotFoundException;

public class Main {

    public static void main(String[] args) throws FileNotFoundException {

        MetodaTrapezow trapez = new MetodaTrapezow();
        RK2 rk2 = new RK2();

        Data.timeStepControl(1e-2, trapez);
        Data.timeStepControl(1e-2, rk2);

//        Data.timeStepControl(1e-5, trapez);
//        Data.timeStepControl(1e-5, rk2);
    }
}
