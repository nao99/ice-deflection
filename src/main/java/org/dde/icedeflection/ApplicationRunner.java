package org.dde.icedeflection;

import org.dde.icedeflection.deflection.DeflectionsCalculator;
import org.dde.icedeflection.deflection.parameter.*;

import java.io.InputStream;
import java.util.Arrays;
import java.util.Scanner;

/**
 * ApplicationRunner class
 *
 * @author  Nikolai Osipov <nao99.dev@gmail.com>
 * @version 1.0.0
 * @since   2021-04-04
 */
public class ApplicationRunner {
    private static final String KSI_FILENAME    = "ksi.txt";
    private static final String P1_FILENAME     = "p1.txt";
    private static final String P2_FILENAME     = "p2.txt";
    private static final String PHYSIC_FILENAME = "physic.txt";
    private static final String PSI_FILENAME    = "psi.txt";
    private static final String X_FILENAME      = "x.txt";

    private static final double[] LAMBDA_EVEN = new double[] {
        2.365020372431352,
        5.497803919000836,
        8.639379828699740,
        11.78097245102023,
        14.92256510455163,
        18.06415775814131,
    };

    private static final double[] LAMBDA_ODD = new double[] {
        3.927378719118806,
        7.068584195523235,
        10.21017612552063,
        13.35176877775915,
        16.49336143134642,
        19.63495408493621,
    };

    /**
     * Runs a deflections calculation process: <br>
     *  - Loads all required parameters
     *  - Starts process
     *
     * @param args application arguments
     */
    public static void run(String[] args) {
        ClassLoader classloader = getClassLoader();

        Psi psi = loadPsi(classloader);
        if (psi.getNumber() > LAMBDA_EVEN.length) {
            throw new IllegalArgumentException("Unable to run application. Maximum for %d vogues calculation possible");
        }

        double[] x = loadX(classloader);

        double[] lambdaEven = Arrays.copyOf(LAMBDA_EVEN, psi.getNumber());
        double[] lambdaOdd = Arrays.copyOf(LAMBDA_EVEN, psi.getNumber());

        Ksi ksi = loadKsi(classloader);
        P P1 = loadP1(classloader);
        P P2 = loadP2(classloader);
        Physic physic = loadPhysic(classloader);

        DeflectionsCalculator calculator = new DeflectionsCalculator();
        calculator.calculate(lambdaEven, lambdaOdd, psi, ksi, physic, P1, P2, x);
    }

    /**
     * Gets {@link ClassLoader} of a current {@link Thread}
     *
     * @return a class loader
     */
    private static ClassLoader getClassLoader() {
        return Thread.currentThread().getContextClassLoader();
    }

    /**
     * Reads psi data from psi file
     *
     * @param classloader a class loader
     * @return psi data
     */
    private static Psi loadPsi(ClassLoader classloader) {
        InputStream is = classloader.getResourceAsStream(PSI_FILENAME);

        PsiReader psiReader = new PsiReader(is);
        Psi psi = psiReader.read();

        psiReader.close();

        return psi;
    }

    /**
     * Reads ksi data from ksi file
     *
     * @param classloader a class loader
     * @return ksi data
     */
    private static Ksi loadKsi(ClassLoader classloader) {
        InputStream is = classloader.getResourceAsStream(KSI_FILENAME);

        KsiReader ksiReader = new KsiReader(is);
        Ksi ksi = ksiReader.read();

        ksiReader.close();

        return ksi;
    }

    /**
     * Reads P1 data from p1 file
     *
     * @param classloader a class loader
     * @return p1 data
     */
    private static P loadP1(ClassLoader classloader) {
        InputStream is = classloader.getResourceAsStream(P1_FILENAME);

        PReader pReader = new PReader(is);
        P p1 = pReader.read();

        pReader.close();

        return p1;
    }

    /**
     * Reads P2 data from p2 file
     *
     * @param classloader a class loader
     * @return p2 data
     */
    private static P loadP2(ClassLoader classloader) {
        InputStream is = classloader.getResourceAsStream(P2_FILENAME);

        PReader pReader = new PReader(is);
        P p2 = pReader.read();

        pReader.close();

        return p2;
    }

    /**
     * Reads physic data from physic file
     *
     * @param classloader a class loader
     * @return physic data
     */
    private static Physic loadPhysic(ClassLoader classloader) {
        InputStream is = classloader.getResourceAsStream(PHYSIC_FILENAME);

        PhysicReader physicReader = new PhysicReader(is);
        Physic physic = physicReader.read();

        physicReader.close();

        return physic;
    }

    /**
     * Reads x data from physic file
     *
     * @param classloader a class loader
     * @return x data
     */
    private static double[] loadX(ClassLoader classloader) {
        InputStream is = classloader.getResourceAsStream(X_FILENAME);
        assert is != null;

        Scanner scanner = new Scanner(is);

        double[] x = new double[1000];

        int i = 0;
        while (scanner.hasNext()) {
            double value = Double.parseDouble(scanner.next());

            x[i] = value;
            i++;
        }

        scanner.close();

        return Arrays.copyOf(x, i);
    }
}
